#ifndef RAD_VET_RAD_VET_HPP_
#define RAD_VET_RAD_VET_HPP_
// VET_DEBUG (cmake -D Athena_VET_DEBUG=ON): the debug hooks (restart fingerprints, poison
// gates) are compiled in; 0 = a production build, where naming them is a FATAL
#ifndef VET_DEBUG
#define VET_DEBUG 0
#endif
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet.hpp
//! \brief definitions for the grey photon two-moment (M1) radiation module, <rad_vet>.
//! See docs/dev/rad_vet_design.md.  MILESTONE 1a: explicit transport only -- the closure,
//! the HLL flux with thick_flux = none, the PD-ARS time integration with a null source,
//! sub-cycling and the cell-centred plumbing.  The implicit matter coupling (design
//! sect. 4) is milestone 1b and enters at the marked Coupling task in both stages.
//!
//! MILESTONE 1b adds: analytic (constant / power-law / user-hook) opacities stored as
//! rho*kappa per cell, the two thick-limit fluxes (ap_hll, scaled), and the implicit
//! local matter coupling to <hydro>, called inside BOTH PD-ARS stages.
//!
//! MILESTONE 1c adds: the advective enthalpy-flux split of the E equation for moving
//! media (<rad_vet>/advect_split, design sect. 3 "Moving fluid"), the unified form of the
//! ap_hll E-flux (<rad_vet>/ap_form, kept as an option -- the 1b pair measures better,
//! design sect. 11), and the O(v/c) source truncation used as the failing control of
//! T4/T4b (<rad_vet>/source_form).  Still absent: implicit transport (stage 3) and any
//! geometry but Cartesian (stage 4).

#include "rad_vet/rad_vet_vgdlg.hpp"
#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
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
#include "utils/time_table.hpp"

// forward declarations
class Driver;

namespace radvet {
struct CleanArgs;   // rad_vet_clean.hpp
struct CleanLevel;  // rad_vet_clean.hpp

//! \struct CleanW
//! \brief the clean work vectors cln_w (CLW_*) in two arrays: slots n < ns in s, the
//! layout's own slots in o (n - ns; the staggered layout's lie in mem_arena between the
//! builds, MemArenaSeal)

struct CleanW {
  DvceArray5D<Real> s, o;
  int ns = 1 << 30;
  CleanW() = default;
  CleanW(const DvceArray5D<Real> &a) : s(a) {}   // NOLINT (one array)
  KOKKOS_FORCEINLINE_FUNCTION
  Real &operator()(const int m, const int n, const int k, const int j,
                   const int i) const {
    return (n < ns) ? s(m,n,k,j,i) : o(m,n-ns,k,j,i);
  }
};

// cubed-sphere geometry (array cln_gn), Cartesian components
constexpr int CLG_N = 0;     // unit normal of the lower face along dir: CLG_N + 3 dir
constexpr int CLG_X = 9;     // cell centroid
constexpr int CLG_R = 12;    // r_hat of the cell centre (the hydro's, hydro_am.hpp)
constexpr int CLG_E = 15;    // the hydro's unit panel tangents e_xi, e_eta (+3)
constexpr int CLG_M = 21;    // first moment int x dA of the lower x2, x3 face (+3)
constexpr int CLG_F = 27;    // staggered, exchanged (a seam ghost's are its owner's): the
                             // normal of the lower x2 face, the volume (+3), the normal
                             // of the lower x3 face (+4) and the MeshBlock gid (+7)
constexpr int CLG_A = 35;    // staggered, exchanged: the mean area of the cell's two
                             // faces of family x1, x2, x3 (+2), CleanStagArea
constexpr int CLG_Q = 38;    // staggered: the face's family and orientation in the
                             // owner's axes of its two cells (CleanCsSlotQ), 4 bits per
                             // x2 (low) and x3 face: q - 1 and s > 0 of the low cell (+2)
constexpr int CLG_NV = 39;
constexpr int CLG_NA = 25;   // components of CleanGeo::a

//! \struct CleanGeo
//! \brief cln_gn (CLG_*) in two arrays: the components that do not change along x1 (the
//! directions: CLG_N, CLG_R, CLG_E and the exchanged normals CLG_F, CLG_F + 4..6) once
//! per (m,k,j) and x1 range (the low ghosts, the MeshBlock's cells, the high ghosts: a
//! corner ghost's exchanged ones are another cell's) in a, the others (CLG_X, CLG_M, the
//! volume CLG_F + 3, the gid CLG_F + 7 and CLG_A) in r, CleanCsGeometry

struct CleanGeo {
  DvceArray5D<Real> r, a;
  int is = 0, ie = 0;
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr int Slot(const int n) {   // a: n >= 0, r: -1 - n
    // (a table, not a chain of branches on n: with n a loop index, the face kernels keep
    // their registers on the GPU)
    constexpr signed char t[CLG_NV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, 9, 10, 11,
      12, 13, 14, 15, 16, 17, -4, -5, -6, -7, -8, -9, 18, 19, 20, -10, 21, 22, 23, -11,
      -12, -13, -14, CLG_NA - 1};
    return t[n];
  }
  // component B + o, B known at compile time, o < the run of components from B that lie
  // together in one of the arrays (CLG_N 9, CLG_X 3, CLG_R 9, CLG_E 6, CLG_M 6, CLG_F 3,
  // CLG_F + 4 3, CLG_A 3): the array and the slot without a lookup
  template <int B>
  KOKKOS_FORCEINLINE_FUNCTION
  Real At(const int m, const int o, const int k, const int j, const int i) const {
    constexpr int s = Slot(B);
    if constexpr (s >= 0) {
      return a(m,s+o,k,j,(i < is) ? 0 : ((i > ie) ? 2 : 1));
    } else {
      return r(m,-1-s+o,k,j,i);
    }
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real operator()(const int m, const int n, const int k, const int j, const int i) const {
    const int s = Slot(n);
    return (s >= 0) ? a(m,s,k,j,(i < is) ? 0 : ((i > ie) ? 2 : 1)) : r(m,-1-s,k,j,i);
  }
};


// indices of the evolved moments in u0(m,n,k,j,i)
constexpr int M1_E  = 0;
constexpr int M1_F1 = 1;
constexpr int M1_F2 = 2;
constexpr int M1_F3 = 3;
constexpr int M1_NVAR = 4;

// number of explicit stages of the PD-ARS (Chu et al. 2019) integrator used per substep.
// With S == 0 it reduces to Heun / SSP-RK2.
constexpr int M1_NSTAGE = 2;

//! ONE canonical key name per concept for both formal solutions
//! (vet_sc, vet_col / vet_gd), the older closure-specific spellings accepted as aliases.
//! Returns the name to read in <rad_vet>: the canonical key if the input has it, else the
//! first alias it has, else the canonical key (a GetOrAdd* on it then records the
//! resolved default).  An alias is read under its own name, so an input (or a restart
//! file) that names only the old spelling behaves exactly as before.  The canonical key
//! wins over an alias given with a different value (printed).
inline std::string M1KeyName(ParameterInput *pin, const std::string &canon,
                             const std::vector<std::string> &alias) {
  const bool hc = pin->DoesParameterExist("rad_vet", canon);
  for (const auto &a : alias) {
    if (!pin->DoesParameterExist("rad_vet", a)) {continue;}
    if (!hc) {return a;}
    if (pin->GetString("rad_vet", a).compare(pin->GetString("rad_vet", canon)) != 0 &&
        global_variable::my_rank == 0) {
      std::cout << "<rad_vet> NOTE: " << canon << " = "
                << pin->GetString("rad_vet", canon)
                << " overrides its old spelling " << a << " = "
                << pin->GetString("rad_vet", a) << std::endl;
    }
  }
  return canon;
}

//! clean_lz_ledger, clean_lz_ledger_k1: the clean L_z ledger and K1's sum (restart), one
//! name on every mesh; a file with the old one (cubed sphere clean_cs_lz_ledger /
//! clean_cs_lz_k1, spherical polar clean_sp_lz_bnd / _k1) keeps it (M1KeyName)
inline std::string M1LzKey(ParameterInput *pin, const bool cs, const bool k1) {
  return M1KeyName(pin, k1 ? "clean_lz_ledger_k1" : "clean_lz_ledger",
                   {cs ? (k1 ? "clean_cs_lz_k1" : "clean_cs_lz_ledger")
                       : (k1 ? "clean_sp_lz_bnd_k1" : "clean_sp_lz_bnd")});
}

// components of the per-cell opacity array, all stored as rho*kappa (inverse length)
// components of RadiationVet::vet_cell (closure = vet_sc)
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
constexpr int M1_VET_HP1 = 23;  // sum w I max(mu_a, 0), a = 1 2 3 (half-range)
constexpr int M1_VET_NL  = 26;  // sum w |I - S| (C1 non-locality), fs only
constexpr int M1_VET_EPS = 27;  // thermal fraction min(kappa_P/chi, 1) (fs)
constexpr int M1_VET_NC  = 28;
// scheme = clean: K+^d_bd = sum w I Omega_b mu_d over mu_d > 0 (slot KP + 3 d + b), after
// the NC slots above (allocated only for the clean scheme)
constexpr int M1_VET_KP  = 28;
constexpr int M1_VET_NKP = 9;
// clean_consistency, after the KP slots (VetShortCharSweep): the rays' flux through the
// low x2 face of the cell, sum w g Phi and sum w g Phi Omega_b (b = 1 2 3), g = |mu_1|
// dx2/dx1 and Phi the intensity of the cell's plane the foot interpolation shifts across
// the face; their flux sum w mu_1 I and sum w mu_1 I Omega_b through the x1 boundary face
// of the cell (entering and leaving) and through its low x1 face, I the intensity of the
// segment where it crosses the face; the flux through the low x3 face as that of x2
constexpr int M1_VET_F2  = 37;
constexpr int M1_VET_F1  = 41;
constexpr int M1_VET_FX  = 45;
constexpr int M1_VET_F3  = 49;
constexpr int M1_VET_NFS = 16;  // the face sums F2, F1, FX, F3
// tau_ten slots of the LATERAL correction (vet_col_lat; the cs CS2 interface): dD_ab in
// the mesh basis, order rr, r-a, r-b, aa, ab, bb (sp: a = theta, b = phi)
constexpr int M1_TT_LAT0 = 4;
constexpr int M1_TT_NLAT = 6;
// the extra vgd_hr components of the thin switch on vet_gd (VetGdHalfRange)
constexpr int M1_VGD_HRC1 = 6;     // C1 = sum w |I - S| / J (twin-corrected as r+-)
constexpr int M1_VGD_HREP = 7;     // eps = min(kappa_P/chi, 1) (1 where chi = 0)
constexpr int M1_VGD_HRJ  = 8;     // J (E units; twin-corrected by J_tw/<J_tw>)
constexpr int M1_VGD_HRNX = 9;
#ifdef M1_VGD_NMAX_BUILD
constexpr int M1_VGD_NMAX = M1_VGD_NMAX_BUILD;   // cmake -D Athena_M1_VGD_NMAX
#else
constexpr int M1_VGD_NMAX = 432;   // vet_gd: at most 432 directions.  kept at 432
#endif
                                    // (sizes the per-thread Real buf[M1_VGD_NMAX] in
                                    // VetGdWall; 2048 would spill on A100).  The 4x-GL
                                    // angular study (sptests P2) needs a local build with
                                    // 2048.
constexpr int M1_OP_P = 0;  // rho kappa_P, Planck (emission) mean
constexpr int M1_OP_E = 1;   // rho kappa_E, energy (absorption) mean
constexpr int M1_OP_T = 2;   // rho (kappa_F + kappa_s), the TRANSPORT opacity: what the
                             // flux source relaxes F with, and what tau_face is built of
constexpr int M1_NOPAC = 3;
// vet_scatter: a 4th component, allocated ONLY when the key is on
// (opac then has M1_NOPAC + 1 components): rho kappa_e, the electron-scattering part of
// M1_OP_T.  M1_OP_T stays the TOTAL extinction rho (kappa_a + kappa_e) = rho kappa_R.
constexpr int M1_OP_S = 3;

// <rad_vet>/opacity
constexpr int M1_OPAC_CONST    = 0;
constexpr int M1_OPAC_POWERLAW = 1;
constexpr int M1_OPAC_USER     = 2;
constexpr int M1_OPAC_TABLE    = 3;   // stellar Rosseland + Planck tables (milestone 2a)

// <rad_vet>/force_reference: what the momentum coupling SUBTRACTS from the radiative
// force before handing it to the gas.  See RadiationVet::arad_ref and the bookkeeping
// note at the head of rad_vet_fluid.cpp.
constexpr int M1_FREF_NONE    = 0;
constexpr int M1_FREF_WB_ARAD = 1;
constexpr int M1_FREF_CLEAN   = 2;   // scheme = clean: its own force of the first solve
constexpr int M1_FREF_WB_MD   = 3;   // wb_arad_md: staggered clean, the steady start's
                                     // force in every direction (cln_fref, hydro walks)

// <rad_vet>/clean_yres: the y of a cell (cell layout) whose balance the gas momentum
// takes

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
  // <rad_vet>/opac_abs_rho_max (default 0 = off): code densities below amask_rho get no
  // absorption (kappa_P = kappa_E = 0) and the transport opacity amask_kes (code units
  // per mass, = opac_abs_kappa_s [cm^2/g] x kunit); see rad_vet.cpp
  Real amask_rho = 0.0;
  Real amask_kes = 0.0;
  // vet_scatter: the split of the tabulated TOTAL kappa_R into true
  // absorption kappa_a = kappa_R - kappa_e and electron scattering kappa_e (Ma+2026
  // sect. 2.2.4).  scat 0 = no split (os = 0, bitwise the old lookup), 1 = constant
  // kappa_e = kes, 2 = kappa_e from the EOS electron fraction, kxe x_e/mu (M1ScatterEos,
  // rad_vet_opacity.hpp).  kes, kfl, kxe in CODE opacity per unit mass (times kunit).
  int scat = 0;
  Real kes = 0.0;    // constant kappa_e, vet_kappa_es (0.34 cm^2/g = 0.2(1+X), X = 0.7)
  Real kfl = 0.0;    // floor: kappa_e <= kappa_R - kfl, vet_kappa_floor (1e-5)
  Real kxe = 0.0;    // sigma_T/m_u = 0.40062 cm^2/g: kappa_e = kxe n_e/(n_tot mu)
};

//----------------------------------------------------------------------------------------
//! \fn VetScatterMix
//! \brief <rad_vet>/vet_scatter (read only when named, default off =
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

// entries of the implicit-solve diagnostic counter (see RadiationVet::cnt)
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
//! \struct RadiationVetTaskIDs
//! \brief container to hold TaskIDs of all rad_vet tasks

struct VetMBState;   // rad_vet_vet.cpp: the multi-block short-characteristics sweep

struct RadiationVetTaskIDs {
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
//! \class RadiationVet

class RadiationVet {
 public:
  RadiationVet(MeshBlockPack *ppack, ParameterInput *pin);
  ~RadiationVet();

  // fin4-dtretry (driver/dt_retry.cpp): <time>/dt_retry.  dtr_on: the driver checkpoints
  // the cycle and takes a failed one again with a smaller dt; then a solve that would
  // stop the run (DIVERGED) or end NON-CONVERGED outside a hesdirk2 stage (whose own redo
  // is the backward-Euler step) raises dtr_fail with the reason instead.
  bool dtr_on = false;
  bool dtr_fail = false;
  // 1 NON-CONVERGED, 2 DIVERGED, 3 DIVERGED (masked), 4 eint<=0, 5 time2_reject
  int dtr_why = 0;
  Real dtr_resid = 0.0;
  // armed by DtRetry (driver/dt_retry.cpp) once <time>/dt_retry_t2_be cuts (default 2)
  // did not cure a time2_reject (the cycle is then retaken at its own dt), or when the
  // previous cycle took that redo: a non-admissible hesdirk2 stage is then redone with
  // backward Euler (as time2_reject = be) instead of raising another retry
  bool dtr_t2be = false;
  bool dtr_t2hit = false;    // this attempt took that backward-Euler redo
  template <class V> void DtRetryVisit(V &v);   // rad_vet_dtretry_visit.hpp (generated)
  void DtRetryQuiesce(bool restoring);          // no host thread writes while restoring
  void DtRetryAfterRestore();                   // async table prefetch state made safe

  // parameters read from <rad_vet>
  Real c_light{};      // speed of light in code units (required)
  Real chat_over_c{};  // reduced speed of light factor (default 1)
  Real chat{};         // = chat_over_c * c_light
  Real e_floor{};      // floor on E
  int chi_kind{};      // M1_CHI_*: levermore (closure = m1) | minerbo | kershaw
  bool eddington{};    // closure = eddington (chi = 1/3); default is the M1 closure
  bool advect_split{};          // the advective enthalpy-flux split of the E equation

  // opacities, per unit mass, in code units (design sect. 4).  kappa = kappa0 for
  // opacity = const, kappa0 (rho/rho_ref)^opac_a (T/t_ref)^opac_b for powerlaw.
  int opacity_type{};           // M1_OPAC_*
  Real kappa_p{}, kappa_e{}, kappa_f{}, kappa_s{};
  Real opac_rho_ref{}, opac_t_ref{}, opac_a{}, opac_b{};
  bool opac_zero{};             // all four opacities are identically zero
  M1OpacTab otab;               // opacity = table: filled by SetOpacityTables

  // <rad_vet>/force_reference (milestone 2a).  With M1_FREF_WB_ARAD the momentum handed
  // to the gas is the RESIDUAL between the instantaneous radiative force and the
  // reference acceleration arad_ref(z) that an external well-balanced scheme already
  // carries (box_convection's Phi_eff).  arad_ref is a per-cell device array filled by
  // the problem generator through SetForceReference.
  int force_ref{};
  DvceArray4D<Real> arad_ref;   // (m,k,j,i), a reference x1 acceleration, code units
  // OPTIONAL explicit energy source of the radiation (erg/cm^3/s, code units), (m,k,j,i).
  // Unallocated and esrc_on = false by default (every path bitwise).  A problem generator
  // that allocates and fills it sets esrc_on = true (transport = implicit only);
  // ImplicitSolve then adds dt_solve*(chat/c)*esrc to the OLD vector of E and the same
  // rate to the hesdirk2 slope, so a constant esrc raises E by exactly dt*esrc per step
  // under be and hesdirk2 alike (stage retries redo the step from U^n with be).
  // he_star_vet: the frozen MLT flux.
  DvceArray4D<Real> esrc;
  bool esrc_on = false;
  // <rad_vet>/vet_source_noesrc (default false = bitwise): the vet_col /
  // vet_col_lat / vet_gd source reads E^n WITHOUT the explicit deposit dt (chat/c) esrc
  // that the implicit solve's old vector E^n carries (the frozen MLT scaffold): the
  // deposit is transported away inside the step by the solve, but a LOCAL relaxed
  // source keeps it in the cell (AG Car A: S up to 10 E at 0.985-0.994 R_ph, L_fs 3 L)
  // F1: every reader (incl. the vet_col inner-boundary intensity) takes
  // the physical E^n = u0 directly (be only; dt-free, right for nsub > 1)
  bool vsrc_noes = false;
  // under time_scheme = hesdirk2 the BACKWARD-EULER solves (the BESTORE
  // start step, the fallbacks, the ic_rad passes) read the polluted E^n + dt esrc too
  // (F1 was be-only).  vet_source_noesrc_h2be (default true; read only when named):
  // noes also for those solves; the stage solves keep iw EN (predict puts P there)
  bool vsrc_noes_h2be = true;
  bool VsrcNoes() const;

  // m1-mhd (docs/dev/m1_mhd_0927.md): the fluid is <hydro> or <mhd> (FluidRef,
  // vet_fluid.hpp).  Under MHD the conserved energy carries |B|^2/2, which every kernel
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
  void FluidRefresh(Driver *pdrive, bool lat = true);
  // <rad_vet>/floorcap_rho (read only when named; 0 = off = bitwise):
  // FLOOR-GAS CAP after the implicit solve, active cells with rho < floorcap_rho only.
  // |v| <= floorcap_v (0 = no velocity cap) and T_gas <= T_rad.  The clipped kinetic
  // energy and the excess internal energy go to E_rad (gas + radiation energy conserved
  // per cell).  Counters printed every floorcap_report cycles (0 = never).  See
  // FloorGasCap.
  Real fcap_rho = 0.0, fcap_v = 0.0;
  int fcap_report = 0;
  bool fcap_mvol = false;   // floorcap_mesh_vol: counters by CleanVol and CleanLzArm
  DualArray1D<Real> fcap_cnt;   // 0 n_v, 1 n_t, 2 dKE, 3 dEint (energy units x volume),
                                // 4 |dp| V, 5 dL_z V (spherical polar;
                                // floorcap_mesh_vol: every mesh) removed
  void FloorGasCap();

  // RSLA start-up check (design sect. 2): v_max*tau_max/chat, evaluated once, on the
  // first filled opacity array.  rsla_vmax <= 0 means "measure max |v| over the domain".
  Real rsla_warn{}, rsla_vmax{};
  bool rsla_force{}, rsla_done{};

  // matter coupling (design sect. 4)
  bool coupling{};     // run the implicit local solve at all
  bool gas_feedback{}; // write the solve's reaction back into hydro u0.  With this off
                       // the medium is a PRESCRIBED static background: the radiation
                       // still feels the opacity, the gas does not move or heat.  That
                       // is what the static-medium tests T3/T3b/T3c mean by "no hydro".
  Real arad{};         // radiation constant in CODE units: the equilibrium energy
                       // density is arad*T^4 with T the EOS's code temperature

  // ---- DEBUG SWITCH.  Defaults to the production value, so that leaving it out of an
  // input file is bitwise inert; it takes the coupling apart when a configuration is
  // unstable and it has to be decided WHICH piece drives it.
  bool dbg_gas_force{}; // <rad_vet>/dbg_gas_force = false: the radiation field still
                       // relaxes its flux, but the gas receives NO momentum and no work

  // ---- VET by SHORT CHARACTERISTICS (rad_vet_vet.cpp), <rad_vet>/closure = vet_sc: the
  // Eddington tensor of the multi-D implicit solve is K/J of a grey formal solution for
  // the specific intensity, computed ONCE per hydro step from the start-of-step state
  bool vet_sc{};             // closure = vet_sc (default false)
  int vet_nmu = 0, vet_nphi = 0, vet_nray = 0;  // mu nodes per hemisphere, azimuths, rays
  bool vsc_gd = false;       // vet_sc_set = gd: the vet_gd direction set, pole on x1
  int vsc_gdmir = 0;         // ... ray r + nray/2 the x1 mirror of ray r (1) or not (-1)
  bool vsc_rayk = false;     // vet_sc_kernel = ray: the layer launch ray-parallel
  DvceArray4D<Real> vet_rayg;   // ... its per-ray face-sum scalars (block, 5 nray, k, j)
  bool vet_milne = false;    // DIAGNOSTIC: source = exact grey Milne S(tau), gate 3
  bool vet_bc_bath = false;  // a MARSHAK x1 end gives the SC rays its bath (m1-h2div)
  // the x1 boundary intensity of the formal solutions, one rule
  // for vet_sc and vet_gd (VetBcInit, rad_vet.cpp).  vet_bc_inner (alias vet_gd_inner on
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
                             // vet_gd, vet_col_lat, vet_col core)
  bool vet_seg_split = false;  // vet_seg = split: the SC segment u -> 0
  bool vet_seg_emit = false;   // vet_seg = emit: split, limiter only where S_u > S_0
                               // as two sub-segments split at the face
  bool vet_seg_jump = false;   // vet_seg = jump (coreL-1012): emit + the contrast rule
  bool vet_seg_plateau = false; // vet_seg_plateau: jump's contrast only across a plateau
  bool vet_xcr = false;        // vet_seg_cross: the cell crossed before the x1 face
  // the segment mode VetSegSplit takes: 0 split, 1 emit, 2 jump with the plateau weights,
  // 3 jump without them (rad_vet_vetseg.hpp)
  int VetSegMode() const {
    return vet_seg_emit ? (vet_seg_jump ? (vet_seg_plateau ? 2 : 3) : 1) : 0;
  }
  bool vet_x1per = false;    // vet_x1_periodic (default false): periodic x1 sweep
  int vet_x1npass = 0;       // ...and its number of passes (runs_3r_radwave)
  // adaptive periodic sweep: stop when the estimated error of the entry plane
  // d_p/(1 - q) < vet_x1_ptol (0 = exactly vet_x1npass passes); counters for VetReport
  Real vet_x1ptol = 0.0;
  DvceArray4D<Real> vet_x1ent;   // (nmb, nray, k, j): the previous pass's exit plane
  Real vet_x1_nbuild = 0.0, vet_x1_npsum = 0.0, vet_x1_ncap = 0.0;
  int vet_x1_npmax = 0;
  int vet_path = 0;          // 1 = the x1 step subdivided at the row nodes
  int vet_interp = 0;        // upwind stencil 0 lin, 1 log, 2 corner
  int vet_mirror_npass = 4;  // passes of a cavity mirrored at both x1
                             // ends (vet_bc_reflect; vet_x1_npass when named)
  bool vet_x1src = false;    // vet_x1_periodic_inflow = source
  // a NON-periodic x2 (x3) mesh boundary (outflow, reflect, user, ...) is
  // a VACUUM face for the formal solution: a ray whose foot leaves the box through it
  // enters with I = 0 (no incoming radiation; single MeshBlock sweep only).  Periodic
  // x2 (x3) keeps the wrap, bit for bit.
  bool vet_x2open = false, vet_x3open = false;
  // the vet_sc angle set and its ROTATION on the vet_gd clock (the same
  // keys vet_gd_rotate_dt | _t0 | _align | _stretch | _slope and the same fresh default
  // auto = 20 x the first dt; VetGdBlock / VetGdAngleBlk / VetGdAlignDt / VetGdRotStep).
  // Double-Gauss: the azimuths about x1 turn by (2 pi/nphi) q; ring: the in-plane angles
  // turn by 0.8 (q - 1/2) (2 pi/N) (never past +-x2, where an x1 sweep has no ray);
  // q = alpha/(pi/2) in [0, 1) the low-discrepancy fraction of the block.
  int vet_ring_n = 0, vet_ring_g = 1;   // vet_ring, vet_ring_groups (0: double-Gauss)
  bool vsc_rot = false;                 // vet_sc rotates (physical clock)
  Real vsc_beta = 0.0;                  // the current rotation of the vet_sc set
  void VetBuildAngles(const Real beta);
  void VetGdAngles(const Real beta);
  void VgdBuildSet(std::vector<double> &d) const;
  void VetRotInit(ParameterInput *pin);   // the shared rotation keys (rad_vet_vetgd.cpp)
  void VetRotUpdate();
  Real VetRotBeta(const int64_t blk) const;   // the set of block blk
  void VetRotSet(const Real beta);
  void VetScRotAvgSetup();     // vet_sc's averaged fields (vet_cell)
  void VetShortCharSweep();    // the build with the current set
  bool vet_axis_flux = false;  // uniaxial axis: the SC flux (true) or the principal axis
  bool vet_full = false;     // vet_tensor = full: the solve reads all six D_ab = K_ab/J
  static constexpr Real vet_eig_min = 0.0;   // vet_tensor = full: eigenvalue floor of D
  Real vet_nguard{}, vet_ncell{}; // full: cell-calls the guard changed D / all cell-calls
  Real vet_guard_max{};      // full: largest |D_guarded - D_raw| (max norm) of the run
  Real vet_d11_min{};        // full: smallest guarded D_11 of the run (x1 line diagonal)
  Real vet_time{}, vet_itime{}, vet_ncall{}; // SC seconds, implicit-solve seconds, calls
  DvceArray2D<Real> vet_ang;   // (nray, 4): mu1, mu2, mu3, weight (sum of weights = 1)
  DvceArray5D<Real> vet_cell;  // (nmb, M1_VET_NC, k, j, i): see M1_VET_* below
  DvceArray5D<Real> vet_ipl;   // (nmb, 2, nray, k, j): intensity of the last two layers
  // several MeshBlocks / MPI ranks (rad_vet_vet.cpp): the halo-banded plane buffers, the
  // neighbour tables and the exchange state of the exact layer-by-layer sweep; null on
  // one MeshBlock (the single-block sweep above is then run unchanged)
  VetMBState *vet_mbs;
                       // of the transverse face flux; 1 = backward Euler (bitwise)

  // true when this module's dtnew must be folded into Mesh::NewTimeStep, i.e. when
  // sub-cycling is off OR when no other module sets the mesh timestep
  bool sets_mesh_dt{};

  // implicit transport (rad_vet_implicit.hpp)
  int nstage{};                 // 1: one stage solve per task-list pass
  Real impl_cfl{};              // <rad_vet>/implicit_cfl: dtnew = impl_cfl*dx/c (<=0 off)
  DvceArray1D<Real> opn_nskip_d;
  // vet-positivity (default true on a fresh run, off on a restart whose file lacks them):
  //  implicit_pos_gas: energy-conserving positivity limiter of the written-back gas eint
  //  implicit_pos_floor: the E floor takes the added energy from the gas (conserving)
  // pos_cnt_d: device counters (M1_POS_*), pos_cnt the all-rank totals at the report
  bool impl_pos_gas = false;
  bool impl_pos_floor = false;
  Real impl_pos_gas_frac = 1.0e-3;
  DvceArray1D<Real> pos_cnt_d;
  Real pos_cnt[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // the E floor's created energy of discarded attempts (dt_retry, backward-Euler redo,
  // helper solves; rank-local, CounterKeep::Restore), not in pos_cnt; all-rank total at
  // the report
  Real pos_un_drop = 0.0;
  int dtr_ncut = 0, dtr_maxcut = 8;   // DtRetry bookkeeping (cuts of this cycle, limit)
  // this cycle's step failed and is redone (clean: solved strictly)
  bool dtr_redo = false;
  DvceArray4D<Real> ktd;        // d(rho kappa_T)/dT at the iterate (implicit_opac_newton)
  Real marshak_q{};             // free-surface condition F_f = c*marshak_q*E
  int ibc_x1min{}, ibc_x1max{}; // M1_IBC_*
  // the transverse (x2, x3) physical faces of the implicit solve,
  // <rad_vet>/implicit_bc_x2min|x2max|x3min|x3max = auto (default: periodic on a
  // periodic mesh, else REFLECT = zero flux, as always) | reflect | marshak (the vacuum
  // free surface F_f = c marshak_q E of the boundary cell, outward; the formal solution
  // enters there with I = 0).  Read only when named.
  int ibc_x2min = 2, ibc_x2max = 2, ibc_x3min = 2, ibc_x3max = 2;
  Real iflux_x1min{}, iflux_x1max{}; // the imposed face flux of M1_IBC_FLUX
  Real fmod_per = 0.0, fmod_amp = 0.0, fmod_f0 = 0.0;   // (test only)
  Real fmod_rt0 = 0.0, fmod_rw = 0.0;   // smooth switch-on (test only)
  // <rad_vet>/implicit_flux_x1min_table = <file> (read only when named;
  // absent = bitwise): the inner face flux F(t) from a (t, value) table, piecewise linear
  // in t - implicit_flux_x1min_table_t0, held at its end values; value = the flux (units
  // = flux, default) or the luminosity (units = luminosity: F = L/(4 pi x1min^2)).  It
  // replaces the constant F0 of the hooks above (the sine / tanh factors, if named, still
  // multiply it) at the same stage times.  InnerFluxAt(t) = the flux these hooks give.
  bool ftab_on = false;
  timetable::Table ftab;
  Real ftab_t0 = 0.0, ftab_scale = 1.0;
  Real InnerFluxAt(const Real t) const;
  Real InnerFluxIntegral(const Real t0, const Real t1) const;   // table only (exact)
  Real iebath_x1min{}, iebath_x1max{}; // M1_IBC_MARSHAK: the INCIDENT bath
  // PER-COLUMN INFLOW columns of the outer x1 Marshak face (S5, the
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
  // <rad_vet>/implicit_bc_x2min|max = flux: the x2 physical face carries an IMPOSED
  // radiative flux F_theta(m,k,i) (lab frame, + = toward larger theta), set by the
  // problem generator through SetX2FaceFlux (an equilibrium lateral flux, e.g. the von
  // Zeipel flux of a tilted stellar photosphere at a wedge wall); explicit in the solve
  bool ifx2_on = false;
  DvceArray3D<Real> ifx2_lo, ifx2_hi;
  DvceArray3D<Real> ibadv_x1max;
  // (implicit_beam_n > 0, vet_sc only): up to 2 collimated beams entering
  // through the x1min face on an x2 patch, carried by one SC ray (2-D: the +-mu3 mirror
  // pair).  ibm_i(m,q,k,j): the entering intensity of beam ray q at column (k,j) (local,
  // active index); ibm_f(m,k,j): its x1 face flux c sum w mu1 I, which replaces the
  // Marshak row of every column it is > 0 on (an imposed-flux face).
  int ibm_n = 0, ibm_nr = 0;
  int ibm_ray[4] = {-1, -1, -1, -1};
  DvceArray4D<Real> ibm_i;
  DvceArray3D<Real> ibm_f;
  DvceArray4D<Real> ibm_s;   // the beams' stress column c sum w Omega_b mu1 I (clean)
  void ImplBeamInit(ParameterInput *pin);
  // implicit_marshak_face = linear (m1-sp-order2, spherical-polar wedge only): the
  // Marshak face flux takes the face E extrapolated from r^2 E of the two end cells
  // (M1SphMarshakCoef in the row, limiter remainder deferred) instead of the end cell's
  // E (the default `cell`, first order).  Read only when named.  The outer face keeps
  // the cell form under vet_col_surface_q, whose q = H(face)/J(top cell) is built for
  // the cell E.
  bool impl_mface_lin = false;
                                // energy density, F_f = +-c q (E - E_bath);
                                // 0 is the plain free surface
  DvceArray4D<Real> f0x1;       // face-normal comoving flux on x1 faces, (m,k,j,i);
                                // PERSISTENT state, carried by the restart file
  DvceArray5D<Real> iw;         // per-cell work array of the solve, M1_NIW components
  Real impl_nstep{}, impl_itsum{}, impl_itmax{};  // Picard counters
  // implicit_thin_freeze (m1-picard-aa; default 1e-2 on fresh runs, 0 = off): cells with
  // c dt rho kappa_P < this at the start of the solve keep their start-of-step opacities
  // for every pass and take the frozen-opacity gas-T find
  DvceArray4D<Real> thin_frz;   // 1 = frozen cell of this solve (implicit_thin_freeze)
  // implicit_resid_fatal (read only when named; 0 = off): stop the run when a Picard
  // solve ends non-converged with resid above this value or non-finite (a DIVERGED
  // solve, not the usual 1e-7 stall), before its state is written back and spreads
  Real impl_res_fatal = 0.0;
  // implicit_timers = N (m1-fast4; default 0 = off): fenced host timers of the M1 stage
  // tasks and of the parts of ImplicitSolve, accumulated over cycles >= N and printed
  // by ImplicitReport.  The fences change the timing a little; only a diagnostic.
  int tmr_c0 = 0;
  bool tmr_on = false;
  double tmr_acc[12] = {}, tmr_cnt[8] = {}, tmr_last = 0.0;
  Kokkos::Timer tmr_t;
  void TmrMark(int c);
  // implicit_mr_every = k (m1-fast4, rad_vet_mr.cpp; default 1 = off): MULTI-RATE
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
  Real mr_thmax = 0.0, mr_ksum = 0.0;   // largest theta seen, sum of window lengths
  // vet_sc_every = N (m1-fast4, rad_vet_mr.cpp; default 1): the vet_sc formal solution
  // on every N-th stage-1 (or multi-rate) tensor build only; in between the tensor is
  // extrapolated linearly in time from the last two formal solutions (O(dt^2))
  int vsc_every = 1;
  int vsc_cnt = 0, vsc_nb = 0;  // builds requested, formal solutions stored (0-2)
  Real vsc_t0 = 0.0, vsc_t1 = 0.0;   // their times, this one
  Real vsc_nskip = 0.0, vsc_nclip = 0.0;          // extrapolated builds, clipped cells
  DvceArray5D<Real> vsc_d0, vsc_d1;               // (m,M1_T2_NVET,k,j,i)
  bool VscSkip();
  void VscStore();
  // coupling_split (ke-dt-0926, rad_vet_mr.cpp; read only when named, default = off):
  // how the implicit radiation + coupling R is placed relative to the Heun hydro H.
  //   0 default  hesdirk2 in the Heun stages (or implicit_mr_every)
  //   1 strang   H(dt) then R(dt) by the multi-rate SDIRK2 every step (k = 1)
  //   2 alternate  H R on even steps, R H on odd steps
  //   3 mix      coupling_mix_h2 hesdirk2 steps, then coupling_mix_s strang steps
  // The step counter travels in the multi-rate restart header (mr_cnt).
  int csplit = 0;
  bool dbg_hydro_off = false;   // DIAGNOSTIC (ke-dt-0926): the Driver skips hydro stages
  // force_reference_work (ke-dt-0926; read only when named, force_reference = wb_arad
  // only): `full` (default) = the solve hands the gas the work of the FULL force and the
  // WB source none (rad_vet_fluid.cpp header; its internal-energy cancellation is
  // O(dt^2) per step); `split` = the WB source gives the gas the work of rho arad_ref at
  // its own stage (the pgen reads fref_wsplit), the solve gives only the residual's
  // work, and the radiation pays both
  bool fref_wsplit = false;
  bool fref_wsplit_ok = false;  // set by a pgen whose WB source gives the reference work
                                // (BEFORE SetForceReference, which resolves `auto`)
  bool fref_wsplit_auto = false;  // force_reference_work = auto, not yet resolved
  // force_reference_work = split with time_scheme = be: the WB
  // source gives the gas rho a_ref v at each hydro stage, and the Heun combination keeps
  // sum_s c_s beta_s dt (rho a_ref v)_s of it; the be solve, which runs once after the
  // whole step, paid v' dt rho a_ref at the post-kick velocity instead, an O(dt^2) per
  // step energy leak (-1.8e-3 L_top at dt 88 s in the BSG).  The pgen now also adds its
  // increment to fref_wacc (acc = gam0_s acc + increment, the same linear combination
  // the gas energy undergoes), and the be write-back makes E pay exactly fref_wacc
  // (times dt_solve/dt_mesh).  hesdirk2 and multi-rate keep theirs.
  DvceArray4D<Real> fref_wacc;
  bool fref_wacc_t2 = false;   // <rad_vet>/time2_fref_wacc (see FrefWaccOn)
  int fref_hnst = 0;
  Real fref_hgam0[4] = {0.0, 0.0, 0.0, 0.0}, fref_hbeta[4] = {0.0, 0.0, 0.0, 0.0};
  bool FrefWaccOn() const;
  // clean staggered, a hesdirk2 stage step (time2_fref_wacc): each stage solve pays the
  // reference work of the hydro stage before it (fref_wacc with gam0 = 0; the Heun
  // average of the radiation's stage start combines the payments as the gas's works),
  // and no slope carries a payment (CleanWriteBack), so the stage increments re-apply
  // none: E pays exactly the work the gas got, as a backward-Euler step does
  bool FrefWaccStage() const;
  // gam0 of the hydro stage whose beta dt is bdt (FATAL if none matches)
  Real FrefWaccGam0(Real bdt, Real dt_mesh) const;
  // the accumulator, (re)allocated to the shape of arad_ref
  DvceArray4D<Real> FrefWacc();
  bool F4RstHdr() const {return (mr_every > 1) || (vsc_every > 1) || (csplit > 0);}
  bool MRActive() const {return mr_every > 1;}
  void MRInit(ParameterInput *pin);
  void MRReport();

  // ---- MILESTONE 3a2 (the limits of 3a).  See rad_vet_implicit.hpp for what each
  // constant means and docs/dev/rad_vet_implicit_design.md sect. 7 for the measurements.
  int impl_flux{};              // M1_IFLUX_*: the spatial form of the implicit E-flux
  int impl_recon{};             // M1_IRECON_*: dc, or plm by deferred correction
  int impl_enth{};              // M1_IENTH_*: face value of the enthalpy flux a E
  // implicit_hr_wscale_x1 / _lat (diagnostic, read only when named): the half-range face
  // coefficients AL, HCL, HCR, DG of the x1 / transverse faces times s (0 = central faces
  // there); off = bitwise implicit_hr_fcimp (read only when named, default false): the
  // half-range face flux in DISSIPATION form, F = F_c(E') + w (c/2) (s_L E_L - s_R E_R),
  // s = h+ - h- (the |mu| half-range moment ratio): the advective part is the evolved
  // central face flux (implicit, as implicit_aphll_fcimp) instead of the lagged
  // formal-solution ratios h+ + h-; AL = 0, HCL = w c s_L/2 >= 0, HCR = -w c s_R/2 <= 0
  // (M-matrix kept); the thin correction DG unchanged (hc = the formal solution's own
  // dissipation at J, with Hf = the mean of H).  Thin-switch faces only.
  // implicit_hr_recon_sfgate (read only when named, default false): plm on h
  // (implicit_hr_recon_q = h) only OFF the beams.  A cell whose own J does not carry the
  // SC face flux (Jh > J, the sf state of M1HrSf: a beam passing through) and its two
  // neighbours along d keep sig = 0 (donor cell, sf ratios, hc, js unchanged); everywhere
  // else (Jh = J: every smooth field) the plm slopes act.  Allows q = h with the sf face
  // ratios and with implicit_thin_recon = js / jse.  Angular data only.
  // implicit_hr_recon_lim = smooth0: as smooth, but where the smoothness test was
  // evaluable (cells c-2..c+2 in the block) and FAILED the slope is 0 (donor cell)
  // instead of van Leer: jumps (beam / shadow edges, fronts) stay first order
  // implicit_hr_recon_hjump = X > 0 (read only when named; 0 off): plm on h off where the
  // half-range ratios jump between a cell and a neighbour along d, |h+_c - h+_n| + |h-_c
  // - h-_n| > X (an opaque-body surface: isotropic inside, outward peaked outside; a
  // shadow edge): those cells keep donor cell.  Angular data only. (opt-in, read only
  // when named): implicit_hr_recon_hjump_w = W (default 1): the hjump gate looks W cells
  // each way along d (a dilation of the gate by W - 1); implicit_hr_recon_kjump = K > 1:
  // plm on h also off where the extinction chi of a cell within W along d differs from
  // the cell's by more than a factor K (opaque surfaces) implicit_hr_recon_hsmooth = C >
  // 0 (opt-in, read only when named): plm on h off where the RELATIVE smoothness r = |dh+
  // - dh-| / (|dh+| + |dh-| + eps) of h+ or h- exceeds C at the cell or within W cells
  // along d (W = implicit_hr_recon_hjump_w; lateral reach <= ng - 1); r = O(dx) in a
  // smooth field at any resolution, -> 1 at a jump or kink.  eps =
  // implicit_hr_recon_hsmooth_eps (default 1e-3).  Per axis d. implicit_hr_recon_hlim =
  // true (opt-in, read only when named): the plm-on-h slope weight sig (van Leer of E)
  // also limited by the van Leer weights of H+ = h+ E and H- = h- E themselves (sig =
  // min(sig_E, sig_H+, sig_H-)): a ratio-preserving limiter, no H+- overshoot where h
  // jumps under a smooth E (shadow / opaque-body edges) implicit_hr_recon_qe_tau = T > 0
  // (opt-in, read only when named; plm on h, not led): a slope centred on a cell OPAQUE
  // along d (tau_d = chi dx_d > T) is the donor ratio times the E slope, h_c (E_c+1 -
  // E_c-1) sig/4, instead of the slope of H = h E (flag in the M1_IM_B slots, read by the
  // row and every M1HqFace) implicit_thin_gface = max (default) | mean | jump (read only
  // when named): the face thinness g_f of the implicit_thin_switch weight S(g_f).
  // max(g_L, g_R) moves the switch by dx/2 against the gradient of g (an O(dx) shift of
  // the blend region: a first-order space error wherever 0 < w < 1); mean = (g_L + g_R)/2
  // (second order); jump = mean + (g_R - g_L)^2/implicit_thin_gface_d (O(dx^2) in a
  // smooth field, -> 1 at a beam edge where g jumps by O(1), as max).  Angular data only
  // (g = |H|/J). (opt-in, read only when named, default false = bitwise):
  // implicit_thin_gface_shell: the x1-face thinness g of the thin switch from the SHELL
  // MEANS of the cell g (repro::ShellSums, layout invariant) of the two shells instead of
  // the two cells (then the gface mean / jump rule on the two shell values);
  // implicit_hr_fsc: inside the blend band the central branch takes the formal-solution
  // net flux c (h_L E_L + h_R E_R)/2, h = r+ + r- (M1HrCell), with weight w: F = w F_hr +
  // (1 - w) [w F_fs + (1 - w) F_central].  x1 faces only; gface_shell: vet_gd (src 2),
  // fsc: vet_gd or vet_sc (src 1, the Cartesian / sp vet_sc rays).
  // implicit_thin_gface_filter = L > 0 (opt-in, read only when named): the x1-face
  // thinness g from the cell g FILTERED along x1 with a hat kernel of fixed physical
  // half-width L (max(0, 1 - |r' - r|/L), renormalised over the valid cells of the
  // column: a resolution-free g profile, O(dx^2) for a smooth g), then the gface rule on
  // the two filtered cell values.  x1 faces only; the column is one MeshBlock along x1
  // (implicit_x1 needs it). diagnostic implicit_thin_gw_dump: the x1-face g_L, g_R (as
  // fed to M1ThinFace) and the weight AL of the last thin-switch face build (output
  // m1_gw)
  DvceArray5D<Real> gwd;
  // time2_face_fct: the restart file carries K1's face energy flux t2pf1 (one
  // cell-sized slab per direction, kM1FaxRstMagic, behind the pgen slabs) once a solve
  // stored it, so a restart limits as the run would.  (The slabs of the removed
  // implicit_hr_cmem_rst in an old file: read and IGNORED, pgen.cpp)
  int FaxRstNslab() const {
    const int nd = trans_on ? (trans_x3 ? 3 : 2) : 1;
    return (t2_pfst && t2pf1_ok && t2pf1.extent_int(0) > 0) ? nd : 0;
  }
  void FaxRstPack(DvceArray5D<Real> out, int nmb);
  bool FaxRstUnpack(DvceArray5D<Real> in, int nmb, int nsl);
  // stage 2b: the transverse (x2, x3; cs seams) faces.  fmm_t (m,8,k,j,i)
  // at the face index: x2 [0] Llo, [1] Lhi (the row's E'-linear coefficients of the face
  // flux, chat units, as M1TransRowCore / m1_impl_tcell build them), [2] G = cr F(E^k) -
  // Llo E_lo - Lhi E_hi (its explicit part, the lagged off-diagonal / lateral term
  // included), [3] lam; x3 [4..7]; [8..11] the row's F_L coefficients of a limited
  // face (x2, x3), F_L = chat ((1 + s) E_lo - (1 - s) E_hi)/2 with s the face's reduced
  // flux at E^k.  iw slots iw_fmm+3..+6: dt A/V of the cell's x2+, x2-, x3+, x3- face.
  // [12], [13] the lagged off-diagonal / lateral part of the x2 / x3 face
  // flux (lab units), kept out of G and unscaled (the operator form re-adds it), and G
  // is built with TRHS (the same iterate) in ImplicitTransverseTerms.
  DvceArray5D<Real> fmm_t;
  // the Zalesak rounds of time2_face_fct / implicit_face_mmat exchange R only
  // when a cell within 2 ng of a MeshBlock face changed (else the ghosts already hold
  // what the exchange would write: bitwise).  <rad_vet>/zalesak_halo_skip (default true;
  // false = exchange every round as before).  Counters: rounds, exchanges done.
  bool zl_hskip = true;
  // <rad_vet>/zalesak_active (default false): a Zalesak round evaluates only
  // the cells whose 7-point stencil changed in the previous round or in the R exchange
  // after it (time2 face FCT: those next to an R that fell, its sign bit set).  A cell
  // whose inputs did not change cannot be newly limited (it was evaluated with the same
  // inputs before): BITWISE.
  bool zl_active = false;
  Real zl_rounds = 0.0, zl_halos = 0.0;
  // <rad_vet>/implicit_thin_switch (M1_THSW_*, rad_vet_implicit.hpp): a
  // dt-free thin/thick weight frozen per solve (formal-solution |H|/J, its non-locality
  // C1, the E^n Knudsen number); implicit_thin_corr adds the lagged formal-solution flux
  int impl_thsw = 0;
  bool thsw_c1 = false;           // the sweep accumulates M1_VET_NL
  Real thsw_h0 = 0.3, thsw_h1 = 0.6;
  // implicit_thin_h0_x1 / _h1_x1 (read only when named; default = h0 / h1):
  // the thin-switch band of the x1 faces only; implicit_thin_hx1_curv = true applies them
  // only on curvilinear meshes (sp / cs radial faces), Cartesian x1 faces keep h0 / h1
  // implicit_thin_recon = js (M1ThinJs: half-range ratios from the face
  // mean of the formal-solution H+-), implicit_thin_bq = lin (outer Marshak q from the
  // face-extrapolated H+ of the top cell), the js factor limit
  // implicit_thin_hcsf = jh: hc/hcs take the sf donor scale Jh as the cell's
  // formal-solution scale on faces with sf ratios (vet_sc); default j = the old J
  // (jh opt-in: +16 % in vacuum crossing beams, NOTE-unify)
  DvceArray5D<Real> ifw;        // per-x1-FACE work array, M1_NIFW components
  // <rad_vet>/implicit_flux_faces = all.  berthon | blend on EVERY face of
  // the multi-D solve: Cartesian x1/x2/x3 and the sp LATERAL x2/x3 faces (the sp radial
  // faces keep the sp-blend-1008 path).  ifw2/ifw3 hold M1_IFW_AL/HCL/HCR of the x2/x3
  // faces (DG unused); the central part of a transverse face keeps 1 - AL through
  // thx2/thx3.  Off (the default): nothing allocated, every kernel bitwise unchanged.
  bool impl_beam_fs = false;    // implicit_flux_beam = fs: |f| = |H|/J of vet_sc
  // implicit_flux_beam = halfrange.  Every face flux of the all-faces blend
  // is (1 - w) F_central + w c (r+_L E_L + r-_R E_R), r+- = H+-/J the HALF-RANGE ratios
  // along the face normal (vet_sc / vet_gd rays; otherwise the isotropic + beam model of
  // the lagged closure), w = exp(-(tau_f/tau0)^2) (implicit_blend = tau) or 1 (berthon)
  bool impl_beam_hr = false;
  // implicit_hr_face = sf (default with vet_sc): the half-range face ratios
  // from the face mean of the formal-solution half-range fluxes over the DONOR SCALE
  // Jh = max(J, the outgoing face SC fluxes) wherever the cell's own J does not carry
  // the SC flux through its faces (M1HrSf); cell = the hrup-1009 cell ratios (opt-out)
  bool impl_hr_sf = false;
  // P2 (opt-in, vet_gd only): implicit_hr_excess = iso gives the part of the
  // lagged E^n that exceeds the formal-solution J of the cell the ISOTROPIC half-range
  // shape (+-1/4) instead of the formal solution's (M1HrExcess); none = off (default)
  double hsf_ncell = 0.0, hsf_nsf = 0.0;   // cell-solves, of which Jh > J (rank 0)
  bool hr_q = false;            // the per-column outer q from h+_1 (vcol_q)
  DvceArray5D<Real> vgd_hr, vgd_hrt;   // (m, 6, k, j, i): r+_r r-_r r+_t r-_t r+_p r-_p
  // with implicit_thin_switch = fs* or an implicit_thin_corr (vgd_thx) vgd_hr
  // carries M1_VGD_HRNX components: + C1 = sum w |I - S|/J, eps, J (E units)
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_pin_cnt;
  DvceArray5D<Real> ifw2, ifw3;
  // <rad_vet>/implicit_gas_src = matched (default: the term-by-term deposit
  // above, bitwise) | conservative: the gas takes minus the radiation's own source of the
  // solve, dm = -(F'_c - F^old_c)/(c chat) - (1/2) sum_faces Z_f, with F_c the STORED
  // cell flux and Z_f the face's transport (and non-momentum-equation) part, written by
  // the face kernels of the solve (rad_vet_implicit.cpp, NOTE-consrc.md)
  bool impl_gas_src_cons = false;
  DvceArray5D<Real> csk1, csk2, csinc;  // hesdirk2: slopes / increment of the cell F (3)
  bool csk1_ok = false;                 // csk1 holds the slope that pairs with t2k1
  bool csinc_ok = false;                // csinc was formed for the coming stage solve
  bool gcs_first = true;                // the first conservative solve of this run
  // csk1 (+ csk1_ok, gcs_first) in the restart file (kM1CskRstMagic):
  // the header whenever a conservative solve has run, the 3 slabs when csk1 holds a slope
  bool CskRstHdr() const {return impl_gas_src_cons && !gcs_first;}
  int CskRstNch() const {
    return (CskRstHdr() && csk1_ok && csk1.extent_int(0) > 0) ? 3 : 0;
  }
  int cs2_vnp = 0;                      // implicit_cs2_vimp_npass (fresh default 3)
  // time2_enth_be: this (non-stage) solve takes the stage a_f rule
  bool T2EnthBe() const;

  // ---- MILESTONE 3b phase B: transport = implicit, the TRANSVERSE (x2/x3) couplings.
  // All inert with transport = explicit | implicit_x1.
  int impl_solver{};            // M1_ISOLV_*: line_jacobi | bicgstab
  bool trans_on{};              // transport == implicit AND the mesh is multi-D
  bool trans_x3{};              // ...and three-dimensional
  DvceArray4D<Real> f0x2, f0x3;   // face-normal comoving fluxes on the x2 / x3 faces,
                                // (m,k,j,i).  PERSISTENT state, carried by the restart
                                // file next to f0x1; allocated only when trans_on.

  // ---- MILESTONE 3b phase C: implicit_solver = bicgstab.  All null under line_jacobi.
  bool bicg_on{};               // impl_solver == M1_ISOLV_BICGSTAB and trans_on
  int impl_prec_fb{};           // <rad_vet>/implicit_precond_fallback: -1 = none (the
                                // default), else the impl_prec code (0 line, 1 rbgs)
  bool bcg_capped{};            // the last BiCGStab stopped at lin_maxit, test unmet
  // ---- GPU cost of the Krylov iteration (bench/m1_fast_0923).  Every key defaults OFF;
  // off, nothing below is allocated or referenced and the path is bitwise HEAD.
  DvceArray5D<Real> odc;        // (m,3,k,j,i) that cache
  DvceArray5D<Real> ost;        // (m,19,k,j,i) that stencil
  bool st_edges{};              // ...and whether any edge coefficient is non-zero
  int impl_onep{};              // <rad_vet>/implicit_one_pass: check period N (0 = off;
                                // default 8 for be, see ImplicitInit)
  Real onep_qa[3]{}, onep_qb[3]{};  // last two measured contractions (be, stage 1, 2)
  Real onep_cnt[3]{};           // solves since the last measurement
  bool impl_onep_auto{};        // <rad_vet>/implicit_one_pass_auto (default true): switch
                                // one_pass off per kind of solve when it does not pay
  Real onep_off[3]{};           // 1: one_pass switched off for this kind (be, s1, s2)
  Real onep_actr[3]{};          // on: eligible solves since the last one-pass acceptance;
                                // off: switched-off solves since it was switched off
  Real onep_prb[3]{};           // 1: on as a re-probe, not yet confirmed by an acceptance
  int impl_pord{};              // <rad_vet>/implicit_predictor_order: 1 or 2 (default 2
                                // for be, see ImplicitInit)
  int impl_prec{};              // <rad_vet>/implicit_precond: 0 = x1 line Jacobi (the
                                // original), 1 = symmetric red-black transverse line
                                // Gauss-Seidel, block-local (no communication), 2 = its
                                // forward half (red, black), 3 = mg (rbgs_fwd + block-
                                // local semicoarsened levels, rad_vet_precond.cpp), 4 =
                                // mg_gc (mg + the global band coarse space)
  // implicit_precond = mg_gc (rad_vet_precond.cpp, tests_m1/runs_5p_coarse2): a GLOBAL
  // coarse space of per-x1-layer (x2,x3)-band means in front of mg; used only under it
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gc_part;  // (m,lb,q,c,i) band sums
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gc_x;     // (gi*nb + b) x_g
  // implicit_precond = mg_gf (rad_vet_precond.cpp, tests_m1/runs_5t_fast5box): a GLOBAL
  // coarse space of the lowest (x2,x3) Fourier modes of each x1 layer in front of mg
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gf_gh;   // (gi*nm + q) multi-rank
  DvceArray1D<int> gf_okd;          // (q) the mode-q factorisation is usable
  Real bcg_r0rel{};             // max|b - A x0|/max|b| of the last BiCGStab call
  DvceArray1D<Real> lay_part;   // layout_invariant_mode = block: per-MeshBlock
                                // partials of M1LayReduce (sized on first use)
  Kokkos::View<Real*, Kokkos::HostSpace> lay_part_h;   // ... their host copy
  Real ew_fprev{}, ew_etaprev{}; // max|r0| and eta of the previous pass of this step
  bool impl_pred{};             // <rad_vet>/implicit_predictor = step (default step*):
                                // start the Picard loop from E^n + (dt/dt_prev) dE_prev
                                // (and T likewise); none = cold start
  bool pred_ok{};               // a previous-step increment is stored
  Real pred_dt{};               // the dt of that step
  DvceArray5D<Real> ipred;      // (m,3,k,j,i): dE_prev, dT_prev, T^n of this step

  // ---- <rad_vet>/time_scheme = hesdirk2 (docs/dev/rad_vet_time2_design.md sect. 5,
  // tests_m1/runs_3x_hesdirk2): two implicit stage solves inside the Heun stages.  With
  // time_scheme = be (default) nothing below is allocated or read.
  int time_scheme{};            // M1_TIME_BE | M1_TIME_HESDIRK2
  bool t2_ok{};                 // a valid FSAL slope K1 is stored
  bool t2_enth_be{};            // time2_enth_be: the non-stage (backward-
                                // Euler) solves of a hesdirk2 run take the stage solves'
                                // a_f rule
  int t2_solve{};               // what the next ImplicitSolve does (M1_T2S_*)
  bool t2_fail{};               // the last stage solve was not admissible
  int t2_dbg_fail{};            // DEBUG time2_dbg_fail: fail stage 1 at this cycle
  // <rad_vet>/time2_reject = be (default) | retry (read only when named): a
  // non-admissible hesdirk2 stage raises the dt_retry flag (reason time2_reject) instead
  // of the backward-Euler redo, so the cycle is taken again at a smaller dt with hesdirk2
  int t2_reject = 0;            // 0 be, 1 retry
  // time2_k1_stiff = s: move the exchange part of K1 (gas -> E) in cells
  // with relaxation k dt > s (Time2FormStage); 0 = off
  Real t2_k1stiff = 0.0;
  Real t2_k1stiff_n = 0.0;      // cell-steps changed (this rank, whole run)
  void Time2Stash();
  // Ebnd (an energy-budget diagnostic, off unless a problem generator sets
  // ebnd_on; sph_shell: <problem>/am_hist_ebnd): the radiation energy the scheme put in
  // through the physical x1 faces (inner minus outer), summed over the accepted steps
  // since the run start or restart.  A hesdirk2 step applies the face flux
  // L = sum A F1 of its three slopes with the tableau weights of U^{n+1} (K1 = the step
  // start, K2 = the accepted stage-1 solve, K3 = the stage-2 solve; a time2_posfix step
  // keeps these weights, its dropped K1 part re-enters stage 2); a backward-Euler step
  // applies L of its solve.  Read-only: no number of the run changes.
  bool ebnd_on = false;
  Real ebnd_acc = 0.0;            // accumulated energy (code units)
  // ebnd_split (set by the he_star_vet drive ledger only): the inner and
  // the outer faces accumulated apart with the same weights (ebnd_acc_in = energy put in
  // through the physical inner x1 faces, ebnd_acc_out = energy LEAVING through the outer
  // ones, both incl. the ebnd_adv term); ebnd_acc is then their difference
  bool ebnd_split = false;
  Real ebnd_acc_in = 0.0, ebnd_acc_out = 0.0;
  // with ebnd_split also the explicit source esrc as the solves apply it: each solve adds
  // its dt (chat/c) esrc to the old vector and to its slope, so the step deposits the
  // slopes' (chat/c) sum esrc V with the same tableau weights (K1 carries the esrc of
  // the previous step's last solve, FSAL)
  Real ebnd_acc_src = 0.0;
  Real ebnd_s1 = 0.0, ebnd_s2 = 0.0;
  // the esrc sum of the last accepted solve = the one inside the next step's K1 (a
  // restart sets it from the drive trailer: the esrc array rebuilt at a restart need not
  // be the one the stored K1 carries)
  Real ebnd_s3 = 0.0;
  bool ebnd_s3_ok = false;
  Real EbndEsrc();
  void EbndStart();
  void EbndStage1();
  void EbndAccept(Real dt, bool t2step);
  // <rad_vet>/time2_accept_fallback (DEFAULT true under hesdirk2 since
  // a run that never takes the fallback is bitwise): a stage whose
  // solve took an operator positivity fallback (implicit_vimp, offdiag operator) for the
  // rest of the step is admissible when it converged and is
  // positive -- the backward-Euler redo it would otherwise take uses the same fallback
  // operator (backward Euler accepts the fallback), so rejecting it only loses order
  bool t2_accfb = false;
  // <rad_vet>/time2_adm_poslim (read only when named; default false = bitwise):
  // the hesdirk2 admissibility test of the written-back gas eint counts a cell as a fault
  // only when the write-back's own conserving limiter (implicit_pos_gas: raise eint to
  // e_min with energy from the cell's radiation above e_floor) cannot cover it.  The
  // stage slope K is taken after that limiter (rad_vet_impl_wb_kernel.hpp), and the
  // backward-Euler redo applies the same limiter, so such a stage is positive and
  // conservative; rejecting it only replaced the whole step by backward Euler (cs
  // accretor: lit floor cells in front of the nozzle, 1-20 cells per rejected stage).
  // E or T <= 0 stays a fault.  Needs implicit_pos_gas = true.
  bool t2_adm_pl = false;
  Real t2_adm_pl_n = 0.0;       // stages admitted by it (would have been rejected)
  // <rad_vet>/implicit_hold (DEFAULT true; a pgen that
  // marks no cell is bitwise either way): cells a problem generator holds (hold_m = 1,
  // set by the pgen where it writes the radiation) are Dirichlet rows of every implicit
  // solve, as an
  // M1_IBC_EFIX end cell: E' = the held E, no stage increment, no vimp kick, no work row;
  // the write-back leaves their E, F and gas to the pgen (NOTE pgrst_1021 implicit_hold)
  bool impl_hold = false;
  DvceArray4D<Real> hold_m;
  Real t2_accfb_n = 0.0;        // stages accepted with a fallback (counter)
  Real t2_nstep{}, t2_nbe{}, t2_nfall{}; // stage steps, BE steps, fallbacks (counters)
  // <rad_vet>/time2_auto_be = N (default 16 USER 10-29; 0 = off,
  // bitwise): after a step whose stages fell back to the backward-Euler redo, the next
  // steps skip the doomed stage attempt and take that same redo path at once (Heun hydro
  // from U^n + the BESTORE backward-Euler solve, t2_ok = false, cold predictor, as after
  // Time2Restore); the stages are tried again after 1, 2, 4, ... N skipped steps (the run
  // doubles while the retries keep failing; a stage step that succeeds resets it).  The
  // state of a skip run in progress is in the restart file.
  int t2_abe = 16;              // N (0 = off; DEFAULT 16, USER 10-29)
  int t2_abe_len = 0;           // current skip-run length (0 = stage mode)
  int t2_abe_left = 0;          // skipped steps left in the current run
  bool t2_abe_now = false;      // this cycle skips the stages
  // plain <rad_vet>/time_scheme = be runs as hesdirk2 with EVERY step a time2_auto_be
  // skip step (Time2AutoBeBegin), i.e. the hesdirk2 backward-Euler path (BESTORE solve,
  // cold predictor, time2_* keys, work fix, enth_be, ...); also selectable under hesdirk2
  // by <rad_vet>/time2_be_only = true (the reference); <rad_vet>/be_legacy = true keeps
  // the old Lie-split be step
  bool t2_beonly = false;
  // <rad_vet>/time2_cfl_post (default true, false = old): the hydro / MHD CFL dt
  // of the next step is taken again after the radiation update, on the gas it wrote.
  // Without it the dt is the stage-2 hydro CFL of the gas BEFORE the coupling: a
  // backward-Euler step (first step, stage fallback redo, time2_auto_be) kicks the gas
  // by a whole dt of the radiation force after it, and the next step ran past its CFL
  // (flz qFB cfl 0.31: x1 CFL dt 1.06 of the kicked gas, step 2.16; a ring emptied to
  // gas E < 0 and the hydro EOS floor created dE/E 1.8e-6, EBUG)
  bool t2_cfl_post = true;
  Real t2_abe_n = 0.0;          // steps taken by time2_auto_be (counter)
  Real t2_abe_ntry = 0.0;       // stage retries after a skip run
  Real t2_abe_nok = 0.0;        // of which succeeded
  // hesdirk2 wall-time split (rank 0, fenced, s; diagnostic only), accepted
  // cycles by kind: clean stage steps, failed stage attempts, their backward-Euler redo,
  // auto_be steps, other BE steps; discarded dt_retry attempts
  Real t2w_clean = 0.0, t2w_fail = 0.0, t2w_redo = 0.0;
  Real t2w_abe = 0.0, t2w_be = 0.0, t2w_retry = 0.0;
  bool Time2AutoBeBegin();
  void Time2AutoBeEnd(int kind, Real tw_a, Real tw_b);
  Real t2_dtprev{};             // the dt of the previous step (vet_sc extrapolation)
  bool t2_vprev{};              // vet_prev holds the tensor of the previous step
  bool t2_vext{};               // time2_vet_extrap (default false: D^n)
  int t2_vsmode = 0;            // time2_vet_sc: 0 lag, 1 predict, 2 rebuild
  // time2_tableau (ke-dt-0926; read only when named): 0 = hesdirk2 (Heun + H-ESDIRK2),
  // 1 = trbdf2 (TR-BDF2 implicit, stage order 2, with the RK2 hydro of c2 = 2 - sqrt 2;
  // the Driver sets the matching gam0/gam1/beta).  Stage coefficients of Time2FormStage:
  // stage 1 t2inc = t2_a21 dt K1; stage 2 radiation start = t2_w0 Y2 + t2_w1 U^n,
  // t2inc = dt (t2_c31 K1 + t2_c32 K2); the diagonal is g = 1 - 1/sqrt 2 in both.
  int t2_tab = 0;
  Real t2_a21 = 0.0, t2_w0 = 0.5, t2_w1 = 0.5, t2_c31 = 0.0, t2_c32 = 0.0;
  Real t2_onep_s{};             // time2_one_pass_safety: implicit_one_pass_safety of the
                                // stage solves (0: implicit_one_pass_safety)
  Real t2_lin_save{};           // the implicit_lin_tol put back after a stage solve
  Real t2_nclip{};              // vet_sc cells whose extrapolated tensor was clipped
  DvceArray5D<Real> t2k1;       // (m,M1_T2_NK,k,j,i) the FSAL slope K1 (restart state)
  DvceArray5D<Real> t2k2;       // (m,M1_T2_NK,k,j,i) the stage-2 slope K2
  DvceArray5D<Real> t2inc;      // (m,M1_T2_NK,k,j,i) old vector - stage start state
  DvceArray4D<Real> t2f1, t2f2, t2f3;  // the face fluxes of U^n (Heun average)
  bool pred2_ok{};
  Real pred2_dt{};
  DvceArray5D<Real> vet_prev;   // vet_sc: the start-of-step tensor of the previous step
  DvceArray5D<Real> vet_now;    // vet_sc: the start-of-step tensor of this step
  DvceArray5D<Real> vet_opac;   // vet_sc: the stage-start opacities, saved over the
                                // formal solution at U^n
  void Time2Init(ParameterInput *pin);
  bool Time2Active();           // the next step runs the stages (else BE)
  void Time2FormStage(int stage, Real dt);
  // sdirk2 with the Heun hydro: the stage-1 solve on the gas at t^n + g dt, the hydro's
  // stage-2 flux at t^{n+1} (Time2GasTime; t2hs the gas U* of Heun stage 1)
  void Time2GasTime(Driver *pdrive, int which);
  DvceArray5D<Real> t2hs;
  int t2_gt = 0;   // the last Time2GasTime step of this stage loop
  bool t2_gtm = true;   // time2_gtime_mom: the flux state's momentum f = 1/g
  void Time2Restore(Driver *pdrive);
  // implicit_vimp_bc (ImplicitVimpBuild; default true on fresh runs with implicit_vimp,
  // false on a restart whose file lacks it): the implicit_bc_advect enthalpy flux of a
  // Marshak x1 end face takes the stage's radiative kick dv of its cell, as the interior
  // faces do under vimp
  bool t2_inprobe = false;      // inside a helper solve (Time2RotCorrect)
  void Time2VetStart();        // vet_sc: formal solution at U^n, then extrapolate
  void Time2VetExtrapolate();
  // closure = vet_col under hesdirk2 (m1-sp-order2b): time2_vet_col = lag (D^n in both
  // stages, first order in time) | predict (ONE build per step at the predicted state
  // U^n + dt K1, at t^{n+1} to O(dt^2)) | rebuild (predict for stage 1, a second build
  // at the stage-1 solution Y1 for stage 2) | extrap (DIAGNOSTIC: D^n + r (D^n -
  // D^{n-1}), history in memory only, not in the restart)
  int t2_vcmode = 0;
  // <rad_vet>/cs_vertex_mono (read only when named; default false = bitwise):
  // on the cubed sphere, a MONOTONE (two-point) lateral transport at the 8 cube vertices.
  // The lateral face flux of the implicit M1 there carries three non-monotone pieces that
  // only the three-panel corner combines: the lagged skew cross term (|cos| -> 1/2 at the
  // vertex, one-sided M1CsTanDer along both axes), the seam-midpoint quadratic of
  // (M1CsSeamMid, one-sided at the vertex: weights 0.657 / 0.437 / -0.095,
  // applied to the gradient, to v g0 and MULTIPLICATIVELY to theta = 1/(1 + chat dt
  // kt), whose along-seam contrast at the photosphere is orders of magnitude, so the
  // extrapolant can change sign) and D_r,lat (central radial difference inside a lateral
  // flux).  RSG E prodG 1681-1691 d: vertex cells driven to the E floor (T_rad 339 K)
  // next to hot vertex cells across the vertex seam faces, then 5.4e42 erg created by the
  // E floor in two days -> 48 kK bubble at the vertex -> the 1704 d plume.  With the key
  // on, every lateral face with a cell within cs_vertex_mono_n cells of a cube vertex
  // takes the two-point gradient (no cross term), its seam faces the foot values (no
  // midpoint quadratic), and D_r,lat is zero in the cells within cs_vertex_mono_n + 1
  // (so every face of those cells sees none).  Conservative (per face), local.
  Real impl_pred_tcap = 2.0;   // seamE-1106 implicit_predictor_tcap
  bool cs_vtx_mono = false;    // DEFAULT = cs geometry (USER 10-29)
  int cs_vtx_n = 1;
  // <rad_vet>/implicit_pos_floor_cons (read only when named; default false =
  // bitwise; needs implicit_pos_floor): the E floor of the write-back creates no energy.
  // The raise to e_floor is paid first by the cell's gas above its floor (as before); the
  // part the gas cannot cover (M1_POS_FLR_UN, which used to be CREATED: RSG E 1691-1693 d
  // 5.4e42 erg at a cube vertex) is taken from the radiation energy above e_floor of the
  // same radial column, nearest cells first (i -+ 1, -+ 2, ..., the pair split by their
  // E above floor), with the donors' stage slopes and |F| <= c E kept consistent.  The
  // column lies in one MeshBlock (nx1 of the block), so the move is local, exact and
  // decomposition-invariant; M1_POS_FLR_MV counts it, M1_POS_FLR_UN only a column with
  // no radiation left above the floor.
  bool impl_pos_flr_cons = true;   // DEFAULT true (USER 10-29)
  // s3c <rad_vet>/implicit_pos_floor_gas_keep th (read only when named; default -1 = off,
  // bitwise): the clean write-back's E floor raise takes from the gas at most (1 - th) of
  // its eint (the rest from the column's radiation, implicit_pos_floor_cons); th = 1:
  // never from the gas
  Real impl_pos_flr_gk = -1.0;
  // s3c <rad_vet>/implicit_pos_gas_cons (read only when named; default false, bitwise):
  // the part of implicit_pos_gas's raise that the cell's E cannot cover (an explicit
  // stage increment's eint deficit, a wind cell near tfloor) is taken from the radiation
  // of the cell's radial column, as implicit_pos_floor_cons does for E, instead of being
  // left to the hydro EOS floor (which creates it); needs implicit_pos_floor_cons
  bool impl_pos_gas_cons = false;
  DvceArray4D<Real> pgc_def;      // (m,k,j,i) the uncovered gas raise x volume
  DvceArray4D<Real> pfc_def;      // (m,k,j,i) the uncovered floor energy of the cell
  bool t2_vcprev = false;       // extrap: vcol_prev holds the previous step's tensor
  Real t2_vcnfb = 0.0;          // predict/rebuild builds that fell back to U^n
  // time2_vet_col_source (read only when named, default relaxed = bitwise):
  // state = the predict/rebuild build takes its VET source from the predicted state
  // itself (vet_col_source = gas there), not relaxed by a further dt (rad_vet_time2.cpp)
  bool t2_vcsrc = false;
  bool t2_vcsbl = false;       // time2_vet_col_source = blend
  void Time2VetSourceInit(ParameterInput *pin, const bool act, const bool col);
  // time2_pred_cap: the predicted P = U^n + theta dt K1 of a predict /
  // rebuild build takes per cell theta = min(1, cap/r), r the relative change of E and
  // of the gas internal energy (0 = off, the unlimited predictor)
  Real t2_pcap = 0.0;
  Real t2_pnlim = 0.0, t2_pncell = 0.0;   // limited / all cell-builds (Time2Report)
  DvceArray4D<Real> t2pth;      // (m,k,j,i) theta of the current build
  DvceArray4D<Real> vcol_prev;  // extrap: (m,k,j,i) f_K of the previous step
  DvceArray3D<Real> vcol_qprev; // extrap: (m,k,j) q of the previous step
  void Time2VetColAt(int which);    // vet_col build at 1 = U^n + dt K1, 2 = Y1
  void Time2VetColSaveY1();         // end of the stage-1 solve: E, T of Y1 (rebuild)
  void Time2PredTheta(const Real dt, const Real gdt, const bool dbg);
  void Time2PredDebugSc();
  void Time2VetColExtrap();
  // the old solver's time2_vet_col = extrap2 restart channels (version < 4, read past):
  // its tau_ten (10) and vgd_hr (10) builds, q, time, alpha, dts
  static constexpr int M1_X2_NCH = 24;
  void Time2Report();
  void RstFingerprint(const char *tag);
  int dbg_rst_fp = 0;
  int dbg_rst_fp_cycle = -1;    // dbg_rst_fp_cycle: fingerprints inside that cycle
  bool dbg_rst_fp_dump = false; // dbg_rst_fp_dump: vgd_cs written at every fingerprint
  bool dbg_rst_fp_done = false; // (the first clean solve of that cycle has printed)
  int Time2RstNchWant();        // restart channels this run keeps (0 unless hesdirk2)
  int Time2RstNch();            // ...and writes (0 unless a slope is stored)
  bool Time2RstAccept(const int nfile, const int ver);   // a file's channels fit (sets
                                // t2_nfv, t2_rstd3)
  int t2_nfv = 0;               // conservative-ray intensity channels of the file
  int t2_rcl = -1;              // the closure channels of the file read (clean_closure =
                                // iterate, CleanRstNcl; -1: this run's own layout)
  int t2_rskip = 0;             // a file's clean_rot_ring channels read with the ring off
  bool t2_rstd3 = false;        // a version-1 file: the fkin block is d only (no K1_d)
  int t2_rstv = 5;              // the clean file's version: <= 3 zero predictor / extrap2
                                // channels and the ring in CLN_N ghosted slabs per
                                // record, 4 the cubed sphere's sc rays' cln_fvi ghosted
  void Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0 = 0, int w0 = 0);
  // <rad_vet>/restart_pack_per_mb (read only when named; default true since
  // was false): the restart writer stages the hesdirk2 channels one
  // MeshBlock at a time (a device buffer of one block instead of all; same bytes)
  bool rst_pack_mb = true;
  void Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb);

  // ---- MILESTONE 3b phase D: the OFF-DIAGONAL Eddington terms and the closure lag.
  // All inert with transport = explicit | implicit_x1 and on a 1-D mesh.
  int impl_offdiag{};           // M1_OD_*: lagged | operator | none
  int od_now{};                 // the mode this STEP is running: impl_offdiag, or
  // <rad_vet>/implicit_closure_thin_relax = C (runs_5c_thinstab; 0 = off, nothing is
  // allocated or read): under implicit_closure_lag = step, the (chi, n) a cell starts a
  // step with is relaxed from the one it used in the previous step with a weight
  // w = min(1, 2/(1 + G^2)), G = C max(chi', b/f)/tau_c + f chi'/min(chi, 1-chi), the
  // bound of the lagged-closure gain.  ctr_mem holds (chi, n1, n2, n3) of the last step.
  Real impl_ctrelax{};
  bool ctr_init{};
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
  int RssRstMode() const;   // RssMode + 16 (the extended restart state)
  int RssRstNch(int mode) const;  // channels written: opac(M1_OP_T) + the tensor's
  //! the array and channel holding restart channel n of a given mode
  void RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c);
  void RssRstPack(DvceArray5D<Real> &dst, int nmb, int mode, int m0 = 0);
  void RssPrime(bool restart);
  DvceArray5D<Real> rss_stage;   // (nmb, nch, k, j, i) read from the restart file
  int rss_stage_mode = -1;       // mode of rss_stage, -1 = nothing staged

  // ---- the TRANSVERSE realizability limiter, <rad_vet>/implicit_trans_limit.
  // Default `none` = not allocated and not referenced, so every earlier configuration
  // is bitwise unchanged.  See ImplicitTransTheta.
  int impl_tlim{};              // M1_TLIM_*: none | lp
  DvceArray4D<Real> thx2, thx3;  // theta_f = 1/(1 + c dt (kt_f + klim_f)) on the x2/x3

  // <rad_vet>/implicit_vimp (tests_m1/runs_3v_vimplicit): the gas velocity of the
  // enthalpy flux implicit through the radiative force of the solve (Newton form)
  bool impl_vimp{};             // the switch (default false: nothing below is touched)
  // <rad_vet>/doppler_face = donor | central (read only when named; default donor =
  // bitwise).  central: the x1 enthalpy (Doppler) flux a E, a = v (1 + chi), of the
  // implicit E row is the face MEAN 0.5 (a_L E_L + a_R E_R) of the two cells' fluxes,
  // with no donor switch and no plm correction (implicit_enthalpy), in the Cartesian and
  // sp rows and at the fin4 ghost-bath wall face.  That is the Athena++ implicit DO
  // near-central face (tau_scheme 1 HLL, f -> 1/tau in thick cells) applied to the
  // cell-centred lab flux (4/3) v E.  Forces implicit_vimp = false.  Athena++
  // reproduction only.
  int dop_face = 0;
  // <rad_vet>/implicit_vimp_ewf (read only when named; default false =
  // bitwise): the lagged enthalpy-work term of a face, E_f (1 + D_f) dv_f with dv_f the
  // face mean of the two cells' predicted kicks, is taken as
  // (1 + D_f) (E_L dv_L + E_R dv_R)/2: each cell's kick carries its OWN E (same order; a
  // strong kick of 1e-8 floor gas beside a
  // pgen-held 73 kK nozzle cell no longer charges it E_nozzle dv_floor: S2 -16.7 there)
  bool vimp_now{};              // on for this step (the positivity fallback drops it)
  // (cubed sphere only): a 3-component exchange (scalar, face-normal pair)
  // of iw slots c0..c0+2, for the dA / s vectors of the conservative sources
  DvceArray5D<Real> d3w;
  MeshBoundaryValuesCC *pbval_d3 = nullptr;
  // the implicit_hr_recon limiter ghosts (hr quad map across a seam)
  MeshBoundaryValuesCC *pbval_sg = nullptr;
  // step 1 <rad_vet>/time2_face_fct = true (read only when named; default
  // false = bitwise; no threshold): Zalesak positivity limiter on the FACE fluxes of the
  // explicit stage increment t2inc.  Each stage solve also stores the E transport flux of
  // its solution per face, t2pf (m,d,k,j,i) = (chat/c) F_applied (the cmem applied flux
  // fax* where it exists), next to its slope (t2pf2 with K2, t2pf1 with K1, rotated).
  // The increment's face part is G_f = A_f (c1 t2pf1 + c2 t2pf2); the rest of the cell's
  // E increment (exchange, enthalpy, metric) is the cell-local remainder, so the split is
  // exact.  A cell whose old vector E + dE is not > 0 gets R = (1 - 1e-12) Q/P (Q = V E +
  // V dE + P - I + I(R): what it holds without face transport plus the inflow its donors
  // still let through, P/I = its out/inflow), a face takes the R of its donor cell,
  // and the E flux and the face F slope of the increment are scaled by it together on
  // both sides: energy exact, gas momentum/energy and the conservative2 cell F untouched
  // (momentum exact).  Rounds until no cell is newly limited (a receiver loses inflow),
  // R halo-exchanged (MeshBlocks, ranks, cs seams).  Unlimited faces/cells bitwise.  What
  // is left (Q <= 0, gas eint) goes to time2_stage_local_be / posfix as before.
  bool t2_ffct = false;
  DvceArray5D<Real> t2pf1, t2pf2;   // (m,3,k,j,i) (chat/c) F_applied of K1's / K2's solve
  bool t2pf1_ok = false, t2pf2_ok = false;
  // t2pf is stored with the limiter or with clean staggered (the history's boundary
  // energy, CleanBndAcc): cln_lacc (m,4,k,j) = the energy the step put through the
  // physical x1 faces so far, inner (0, + inward) and outer (1, + outward), h A Phi of
  // its solves and the old vector's face part (2, 3: that of the stage's old vector);
  // cln_ldt = the step's dt (0: none yet; the history then reads the faces' F)
  bool t2_pfst = false;
  DvceArray4D<Real> cln_lacc;
  Real cln_ldt = 0.0;
  void CleanBndAcc(int stage, Real c1, Real c2, Real w, Real dt);
  Real t2_ffct_calls = 0.0, t2_ffct_cells = 0.0, t2_ffct_left = 0.0;
  Real t2_ffct_rmax = 0.0, t2_ffct_de = 0.0, t2_ffct_nop = 0.0, t2_ffct_fill = 0.0;
  // t2_ffct_emin = th (0.1 for clean staggered, else 0 = as before): the limiter keeps
  // th V E^n in a limited cell instead of draining it to zero (a stage started from an
  // emptied cell refills a vacuum: the gas takes its radiation pressure gradient as a
  // kick); counted: limited cells left
  // below th E^n (their own inflow-free content is less) and min E/E^n of limited cells
  Real t2_ffct_emin = 0.0, t2_ffct_nemp = 0.0, t2_ffct_qmin = 1.0e300;
  void Time2FaceFct(int stage, Real c1, Real c2);
  bool am_vec = false;             // + <hydro>/am_vector: the vector form
  DvceArray4D<Real> amk1, amk2, aminc;   // spin slopes (FSAL / stage 1) and increment

  // ---- MILESTONE 3e: ANDERSON ACCELERATION of the Picard map, <rad_vet>/implicit_accel.
  // Default `none` = nothing below is allocated and nothing is called, so every earlier
  // configuration is bitwise unchanged.  See ImplicitAccelApply.
  DvceArray5D<Real> aa_xc;      // (m,nc,k,j,i) x_k, the SCALED state entering the pass
  DvceArray5D<Real> aa_fc;      // (m,nc,k,j,i) g_k = G(x_k) - x_k, scaled
  DvceArray5D<Real> aa_xp, aa_fp;  // the previous pass' pair, for the differences
  int aa_nh{};                  // history columns in use
  int aa_head{};                // ring head (the column the next difference overwrites)
  bool aa_hasp{};               // aa_xp/aa_fp hold a usable previous pass
  Real aa_fnp{};                // ||g|| of the previous pass (< 0: none)

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
  Real dtnew{};   // implicit_cfl*min(dx)/c (implicit_cfl > 0), else FLT_MAX
  Real dt_sub{};  // the dt of the stage solve being taken
  int nsub{};     // always 1 (no sub-cycling)

  RadiationVetTaskIDs id;

  // functions
  void AssembleRadVetTasks(std::map<std::string, std::shared_ptr<TaskList>> tl);
  //! set dt_sub = dt_mesh (nsub = 1); returns nsub
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
  // ---- milestone 3a (rad_vet_implicit.cpp)
  //! read the implicit-solver parameters and allocate its arrays; a no-op in
  //! transport = explicit
  void ImplicitInit(ParameterInput *pin);
  // (rad_vet_icrad.cpp): opt-in IC radiation start, fresh starts only
  int icr_mode = 0;             // <rad_vet>/ic_rad: 0 none, 1 fs, 2 relax, 3 fs_relax
                                // (default 3 when eligible, IcRadInit)
  // (set by a pgen, default false = nothing): at the end of IcRadRun, ONE more
  // pass at the run's dt with the radiative force on the gas, whose net momentum
  // deposit dm - dt rho a_ref per cell (x1 component) is kept in icr_dm (icr_dm_dt = the
  // dt it took); the gas momentum and kinetic energy are then restored (v as before)
  bool icr_fmeas = false;
  DvceArray4D<Real> icr_dm;
  Real icr_dm_dt = 0.0;
  //! STAGE S1 (rad_vet_sph.cpp): the configuration the spherical-polar wedge supports;
  //! fatal on anything else.  Called after ImplicitInit when sph_geom is set.
  void SphericalS1Check(ParameterInput *pin);
  // true on a spherical-polar mesh (a wedge clear of the poles, S1): the implicit
  // kernels then take the face areas, cell volumes and centre-to-centre distances of
  // Coordinates instead of the uniform mb_size.dx1..3 (appended as overwrites, so
  // the Cartesian arithmetic is untouched)
  bool sph_geom = false;
  // STAGE S2 (rad_vet_sph.cpp): sph_geom with a closure that is not Eddington (m1,
  // minerbo, kershaw): the radial faces take the integrating factor (M1SphDrr) and every
  // face equation the lagged curvature M1SphCurv; false keeps the S1 arithmetic
  bool sph_q = false;
  // STAGE CS0 (rad_vet_sph.cpp, m1-cs-implicit): the cubed sphere.  sph_geom is set there
  // too (x1 = r and the Coordinates areas, volumes and face distances exist on both
  // meshes), and cs_geom restricts it to transport = implicit, closure = eddington,
  // time_scheme = be, bicgstab + rbgs_fwd; the transverse faces are the plain two-point
  // form (no skew 1/sin, no cross term: stage CS1)
  bool cs_geom = false;
  // STAGE CS1 (rad_vet_sph.cpp CubedS1Init): the transverse operator on the cubed sphere.
  //   m1bcs     mb_bcs with BoundaryFlag::panel read as block, so that a panel seam is an
  //             OPEN face of the implicit operator (CS0 closed it); the kernels take it
  //             in place of mb_bcs only on cs (elsewhere the View they read is mb_bcs).
  //   cs_seam   per MeshBlock: its inner_x2, outer_x2, inner_x3, outer_x3 face is a seam.
  //   csg2/3    the canonical seam geometry (M1CsSeamGeom) of the x2 / x3 seam faces,
  //             (m, {dth, angm, sfoot, smid}, along-seam cell incl. ghosts, side lo/hi).
  DualArray2D<BoundaryFlag> m1bcs;
  DualArray2D<int> cs_seam;
  DvceArray4D<Real> csg2, csg3;
  // STAGE CS3: the effective transverse geometry of the cs operator, for kernels that
  // read Coordinates through the sp branches (ImplicitVimpBuild): dxface with the x2/x3
  // two-point distances (centre arc x sin of the face angle; r angm at a panel seam)
  // and area with the canonical seam areas; x1 parts are copies
  DvceFaceFld4D<Real> cs_dxf_eff{"m1csdxf", 1, 1, 1, 1};
  DvceFaceFld4D<Real> cs_area_eff{"m1csarea", 1, 1, 1, 1};
  // the planar unit normals (Cartesian) of the x2 faces (m, j = 0..n2, c)
  // and x3 faces (m, k = 0..n3, c) of each block, oriented along +x2 / +x3 (M1CsLat0)
  DvceArray3D<Real> cs_fn2, cs_fn3;
  // the same normals at the CELL CENTRES (m, j or k incl. ghosts, c): the
  // lateral axes of the cs half-range ratios (VetGdHalfRange on cs)
  DvceArray3D<Real> cs_cn2, cs_cn3;
  void CubedS1Init();
  //! the seam average of the stored face state f0x2/f0x3 (stage CS1, C5)
  void CubedSeamFaceAverage();
  bool cs_seam_avg_on = true;
  ParameterInput *pin_cs_ = nullptr;
  //! print the timers and positivity totals of the implicit transport (destructor)
  void ImplicitReport();
  // <rad_vet>/dbg_energy_tally (default false, output only): rank 0 prints
  // at the start and end of every be solve the global sums of the gas total energy, E,
  // the reference-work accumulator and esrc (x volume): line ETALLY
  bool dbg_etally = false;
  // <rad_vet>/dbg_energy_audit (default false, output only): per implicit
  // solve, the change of U = gas total energy + (c/chat) E (x volume, active cells) in
  // the segments between the marks EAuditMark(0 closure-limit entry, 1 exit, 2 solve
  // entry, 3 solve exit, 4 after the pgen m1post hook, 5 after the fluid refresh), per
  // cell class (0 rest, 1 floor gas rho < 10 dfloor, 2 reset: E overwritten outside the
  // M1 chain since its last end, i.e. pgen hold/reset cells and closure-limit raises),
  // the boundary energy fluxes of the solve and its hesdirk2 old-vector increment; rank
  // 0 prints line EAUDIT at mark 5 (see EAuditMark)
  // mark 6 = the state was put back outside the M1 chain (Time2Restore):
  // the change since the last mark (the discarded attempt) is summed here
  // instead of into a segment; printed as D= at mark 5
  // F4: ApplyClosureLimits clips over ACTIVE cells (0 = E raised to
  // e_floor, 1 = |F| scaled to c E), every stage, and the event-log totals (collective)
  DvceArray1D<Real> acl_cnt_d;
  void EventTotals(Real *v);
  //! VET (rad_vet_vet.cpp): read <rad_vet>/vet_*, check the mesh, allocate
  void VetInit(ParameterInput *pin);
  //! VET: the short-characteristics formal solution and the uniaxial tensor (chi, n)
  void VetShortChar();
  //! VET, vet_tensor = full: the guarded D = K/J (and its ghosts) for the solve
  void VetFullTensor();
  //! VET: write the column of j = js (m = 0) of the last formal solution to `fname`
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
  //! the whole backward-Euler step, in place of the explicit stage chain
  // <rad_vet>/scheme = clean: implicit stage solve that replaces ImplicitSolve
  // (rad_vet_clean*.cpp, docs/dev/rad_clean_scheme_design.md)
  TaskStatus CleanSolve(Driver *d, int stage);
  void CleanInit();
  void CleanPreflight(ParameterInput *pin);
  void CleanCsGeometry();
  void CleanSteadyStart(Driver *pdrive);
  void CleanInitFaceFluxes();
  void CleanRstFaces();
  void CleanOldVector();
  void CleanClosure();
  void CleanPredictState(const bool use_k1);
  void CleanRayMoments(DvceArray5D<Real> dst, const bool solve = true);
  void CleanClosureExtrap();
  void CleanClosureAfterSolve();
  void CleanClosureIterate();
  void CleanStageOpacity();
  void CleanRays();
  void CleanClosureFromRays();
  void CleanClosureFromGd();
  void CleanRayCut();
  void CleanGdFaceSums();
  const CleanArgs *cln_gdx = nullptr;   // their CleanArgs in device memory (slot 3)
  void CleanFvSets(const Real tc, const int cyc);
  void CleanRingMoments();
  void CleanRingShift(const int64_t b);
  int64_t CleanRingBlk(const Real dt) const;
  bool CleanRingTwo() const;
  bool CleanFhBuilt() const;
  void CleanRingPreShift();   // dt_retry: the shift before the checkpoint
  void CleanRingJump();
  void CleanRingSolveSet();
  void CleanClosureFill();
  void CleanClosureEddington();
  void CleanFaceCoefficients();
  void CleanAdvectiveCoefficients();
  void CleanFirstIterate();
  // (beta: CLW_BT, CLW_RE only; lat false: CLW_L holds E_k's lateral sums)
  void CleanRightHandSide(const bool beta = false, const bool lat = true);
  void CleanPreconditionerSetup();
  void CleanPrecondition(const int src, const int dst);
  CleanLevel CleanStagLevel(const int l);
  void CleanStagLineSolve(const CleanLevel &lv, const int src, const int dst,
                          const int col, const int lim = -1);
  void CleanStagCoarseSetup();
  void CleanStagVCycle(const int l, const int src, const int dst);
  void CleanApply(const int src, const int dst);
  int CleanBiCGStab();
  Real CleanNewtonUpdate();
  Real CleanVelocityUpdate();
  void CleanCopyIterate();
  int CleanCountNegative();
  int CleanStagRealizable();
  void CleanStagBracket();
  void CleanStagFSolve(const int e, const bool tie = true, const bool lat = true);
  void CleanDMem(const bool save);
  // the rays' history a failed stage's builds change (intensities, face half ranges, the
  // mirrored top, the ring, the swept sets' tables, the ray cut, a
  // restart's first-fill flags): kept at the stage-loop start, put back with U^n
  // (Time2Stash, Time2Restore): the backward-Euler redo builds from the step's.  Also
  // where the sets are not built yet (a restart's first step): the redo builds them anew
  struct CleanRayKeep {
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
  void CleanRayHistory(const bool save);
  void CleanDInc(const bool two, const Real c1, const Real c2, const Real w1,
                 const Real w0);
  void CleanDtrD(const bool save);   // dt_retry: d to / from cln_rstd
  int CleanRingFields(int *fm, const bool file = false) const;   // clean_rot_ring
  // time2_face_fct's two scratch slots of iw (scheme = clean: the line-solve slots, its
  // iw ends before the Krylov vectors)
  int Time2FctSlot(const int q) const;
  void Time2FctHalo(const int c);
  int CleanRingRstN() const;            // ... in the restart files (version 4)
  // (vet_rotate_avg_ema, production's EMA: a fill takes 2 sets; cln_rr = the EMA offsets
  // e, the aux J and D11 and the set count n, cln_rpd = the Delta, a_ch and the next
  // set's J and D11 (CleanRingEma))
  int CleanRingK() const {return vra_ema ? 2 : vra_k;}
  int CleanRingNr(const int nf) const {return vra_ema ? nf + 3 : (vra_k - 1)*nf;}
  int CleanRingNd(const int nf) const {return vra_ema ? nf + 3 : nf;}
  void CleanRingRstCopy(DvceArray5D<Real> &a, const int c, const int nmb, const int m0,
                        const bool pk, const bool have);
  // closure components in the restart files of the extrap history and clean_rot_ring:
  // push 16's order (CleanRstSlot); the cell layout keeps its first CLN_N - 1, CLN_JR is
  // staggered only
  int CleanRstNq() const;
  // the rays' intensity channels of the restart files: the cubed sphere's sc rays'
  // (cln_fvi), and their last build when a reflecting top mirrors it (cln_scmem: its top
  // shell)
  int CleanRstNfv() const;
  int CleanRstNcl() const;   // clean_closure = iterate: the closure in place (cln_s)
  void CleanScUnion2();
  void CleanFaceHalf(const DvceArray2D<Real> &dv, const int o, const bool zero,
                     const int fo);
  DvceArray5D<Real> &CleanRstFv() {return (cln_scmem && !cln_cssc) ? cln_top : cln_fvi;}
  void CleanStagDrag(const bool first, const bool vel = false, const bool lat = true);
  void CleanStagFCsShared();
  void CleanStagFCone(const int e, const int kd);
  void CleanStagFwbProject(const int d, const int slot, const DvceArray4D<Real> &rs,
                           const DvceArray4D<Real> &rw, const bool cnt);
  void CleanStagFaceTie(const int e, const int kd, const bool pt = false);
  void CleanStagCsTie(const int nc, const int c0, const int c1 = 0, const int c2 = 0);
  Real CleanStagFSweepChange();
  void CleanStagFSchur();
  void CleanStagFSchurTile(const int gr[3], const int gs[3], const bool per[3],
                           const int ngd);
  void CleanStagCsLateral(const int e, const bool st = false, const int part = 0);
  void CleanWriteBack();
  void CleanWbFct(const CleanArgs &x, const CleanArgs *px, const DvceArray5D<Real> &fp,
                  const bool pfo, const DvceArray5D<Real> &pf);
  CleanArgs CleanMakeArgs();
  Real MeshBlock0Dx(const int d) const;
  const void *CleanDevCopy(const void *p, const size_t n, const int s);
  template <class T> const T *CleanDev(const T &x, const int s) {
    return static_cast<const T *>(CleanDevCopy(&x, sizeof(T), s));
  }
  const CleanArgs *CleanDevArgs(const CleanArgs &x, const int s);
  template <typename T>
  void CleanGhostFillT(DvceArray5D<T> a, const int c0, const int nc, const bool phys,
                       const int part, const int dep = 0,
                       const std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u});
  void CleanGhostFill(DvceArray5D<Real> a, const int c0, const int nc,
                      const bool phys = true, const int part = 0, const int dep = 0,
                      const std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u});
  // (dep: the ghost layers its readers need, 0 = all ng; clean_halo_direct sends those)
  // (fm[f], clean_halo_direct: the components the ghosts of family f take, bit n for
  // component c0 + n; f the axis of their face in the axes of the cells they copy, 3 at
  // edges and corners (CleanHaloInit); the others keep their values, POISON: NaN)
  template <typename T> struct CleanGSeg {
    DvceArray5D<T> a; int c0, nc; bool phys; int dep = 0;
    std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u};
  };
  template <typename T>
  void CleanGhostExchange(const std::vector<CleanGSeg<T>> &sg, const int part);
  template <typename T>
  void CleanGhostPhys(DvceArray5D<T> a, const int c0, const int nc);
  void CleanGhostDefer(DvceArray5D<Real> a, const int c0, const int nc,
                       const int dep = 0);
  void CleanGhostFill(const CleanW &a, const int c0, const int nc,   // (in one array)
                      const int part = 0, const int dep = 0,
                      const std::array<uint32_t, 4> fm = {~0u, ~0u, ~0u, ~0u});
  void CleanGhostDefer(const CleanW &a, const int c0, const int nc, const int dep = 0);
  void CleanGhostFlush();
  std::array<Real, 2> CleanDots(CleanW w, const int nc, const int n,
                                const std::array<int, 4> ab);
  Real CleanReduce(Real v, const bool max);
  Real CleanSumCells(const DvceArray4D<Real> &a);
  void CleanClosureSphere();

  Real cln_c0 = 0.1;            // clean_w_c0: C1 below which w = 0
  Real cln_c1 = 0.5;            // clean_w_c1: C1 above which w = 1
  Real cln_tol = 1.0e-10;       // clean_tol: relative residual of the linear solve
  // clean_init_cnorm: the steady start's per-cell test (0 off)
  Real cln_cnorm = 0.0;
  bool cln_insteady = false;    // in CleanSteadyStart's passes
  Real cln_lcnorm = 0.0;        // clean_lin_cnorm: the steps' strict solves (0 off)
  int cln_maxit = 800;          // clean_maxit: BiCGStab iterations
  bool cln_pcr = true;          // clean_stag_pcr: PC line solves by PCR on a GPU
  int cln_optn = 3;             // clean_stag_opac_passes: refreshed in passes < N
  bool cln_optdr = false;       // clean_stag_opac_drag: the refresh rebuilds the drag
  bool cln_hface = false;       // clean_face_opac = harmonic: the x1 face chi
  bool cln_fspcr = false;      // clean_stag_fs_line = pcr: F solve x1 lines on a GPU
  // the operator's rows (CLN_ST per active cell, CleanPreconditionerSetup)
  DvceArray5D<Real> cln_st;
  //   (cubed sphere: CLN_CSY coefficients of the lateral sums)
  DvceArray5D<Real> cln_csy;
  bool cln_csyok = false;       //   cln_csy is of the current closure (no rebuild)
  bool cln_ysum = false;        // clean_stag_lat_sums: the sums Y off the cubed sphere
  DvceArray5D<Real> cln_yl;     //   (in their own array: CLW_L is the F solve's scratch)
  int cln_stall_n = 0;          // clean_bicg_stall_its: stall window (0 = off)
  int cln_fsweeps = 1;          // clean_stag_fsweeps: line-solve sweeps of the F solve
  bool cln_fschur = false;      // clean_stag_fschur: the F lines exact across MeshBlocks
  bool cln_tgid = false;        // clean_stag_tie_gid: sp/C owner the lower gid (cs tie)
  Real cln_stall_f = 0.99;      // clean_bicg_stall_fac: least residual drop per window
  int cln_maxnl = 50;           // clean_maxnl: Newton iterations
  Real cln_tol_nl = 1.0e-11;    // clean_tol_nl: Newton tolerance (relative dE, dT, dv)
  bool cln_gtf = false;         // clean_gas_tfloor: gas T >= tfloor, energy to/from E
  bool cln_wbfct = false;       // clean_wb_fct
  bool cln_wbff = false;        // clean_wb_floor_fail
  bool cln_wbxl = false;        // clean_wb_exch_local
  Real cln_flrun = 0.0;         // the last write-back's implicit_pos_floor creation
  int cln_flrid = -1;           // its cycle and t2s, 4 ncycle + t2_solve
  Real cln_wbffn = 0.0;         // its steps accepted at the last cut
  int cln_wbfdg = 0;            // clean_wb_floor_diag (diagnostic)
  Real cln_raymin = 1.0e-3;     // clean_ray_min: rays describe a cell when J >= this E
  // clean_ray_cut (CleanRayCut; default on where it applies): no rays below the shell
  // where every deeper cell is in the diffusion regime, R = |grad S|/(chi S) <
  // clean_ray_cut_r and no cell thinner than vet_col_lat_cut_thin, less a margin of
  // clean_ray_cut_margin optical depths and at tau_top >= clean_ray_cut_taumin; it rises
  // only after clean_ray_cut_wait builds
  bool cln_rcut_on = false;
  Real cln_rcut_r = 0.01, cln_rcut_mg = 10.0, cln_rcut_tau = 0.0;
  int cln_rcut_wait = 10, cln_rcut = -1, cln_rcut_n = 0;
  ParameterInput *pin_rcut_ = nullptr;
  int cln_fillpass = 8;         // clean_fill_passes: closure fill layers
  bool cln_fillcs = false;      // clean_fill_cs_remap (CleanClosureFill)
  bool cln_jsec = false;        // clean_face_jse_cell (CleanFaceCoefficients)
  bool cln_pred = true;         // clean_predict: closure of the predicted state
  bool cln_pst = false;         // clean_pred_stage (CleanFirstIterate)
  int cln_gsec = 0;             // clean_gas_secant (CleanRightHandSide)
  bool cln_pxc = false;         // clean_pred_xch (CleanFirstIterate)
  bool cln_kpred = true;        // implicit_predictor: the solve starts from E + h K1_E
  bool cln_k1ok = false;        // K1 = f(U^n) valid (an accepted solve wrote it): the
                                // predictor and the budget's time error also in BE steps
  int cln_bdry_form = 0;        // clean_bc: 0 rays' half range, 1 Marshak
  int cln_maxpic = 12;          // clean_maxpic: passes that update the lagged slopes
  bool cln_gflux = false;       // clean_stag_gflux: rays' face flux as flux factors
  bool cln_csvr = false;        // clean_stag_csvec_radial (CleanStagCsVector)
  bool cln_rayratio = false;    // clean_stag_rayratio: the bracket of ray ratios x E
  bool cln_ringjmin = true;     // clean_ring_jmin: ring J ratios only of ray-held cells
  bool cln_rfrst = false;       // clean_rot_ring_fill_rst: full window at an old restart
  bool cln_rfdone = false;      // its fill taken
  bool cln_rmsh = false;        // clean_rst_mshare: the M1 share's memory at an old restart
  bool cln_rmdone = false;      // its start taken
  bool cln_linpass = false;   // clean_lin_pass: a pass goes on past a failed linear solve
  void CleanStagXDonor(const bool first, const bool lat = true, const bool hold = false);
  void CleanWorstCell(const Real dnl, const Real dnl0);
  bool cln_srccvx = false;      // clean_source_convex: vet_scatter S between J and B
  bool cln_freal = false;       // clean_stag_freal: realisable sb, restart faces
  // clean_stag_xdonor (staggered, default off): the cross flux of the iterate as a
  // donor-cell coefficient of its face (CleanStagXDonor); its face slot, set in a solve
  bool cln_xd = false, cln_xd_on = false, cln_xlag = false, cln_xsrc = false;
  int cln_xqs = -1;   // (5 slots: qL >= 0, qR <= 0, the lagged explicit flux, C(E_k),
                    //  C's last change)
  bool cln_eicry = false;   // clean_stag_eicarry (CleanOldVector)
  bool cln_colf = false;    // (clean_stag_colflux: column moves in the stored Phi)
  int cln_lza = 0;          // clean_cs_lz_abs
  bool cln_lzbook = false;  // (the L_z ledger: cubed sphere, or 3-D spherical polar)
  bool cln_lzt = false;     // clean_cs_lz_tel
  bool cln_lzbf = false;    // clean_cs_lz_bnd_free
  bool cln_thinf = false;   // clean_stag_thinf
  DvceArray1D<Real> cln_mclc;    // (the Newton passes' momentum cap: the terms it moved)
  bool cln_fwb = false;     // clean_stag_fwb
  DvceArray3D<Real> cln_fsjt[3];   // (layout tile: the F join per tile, CleanFsjt)
  DvceArray1D<Real> cln_fwbc;    // (fwb: the faces, their momentum sums |dp|, x1 dp)
  bool cln_fcone = false;   // clean_stag_fcone
  bool cln_fcgas = false;   // clean_stag_fcone_gas
  bool cln_fwbst = false;   // clean_stag_fwb_stage
  bool cln_fwbpj[3] = {false, false, false};   // clean_stag_fwb_line, -, _ring
  int cln_fwbs = -1;        // (fwb_stage: the face slot of s; cln_fwbok: s is the round's)
  bool cln_fwbok = false;
  bool cln_fwbt2 = false;   // clean_stag_fwb_t2
  DvceArray1D<Real> cln_fconc;   // (fcone: the faces scaled, their momentum sum)
  int cln_xdfz = 0;         // clean_stag_xdonor_freeze
  bool cln_xdh = false;     // (its hold in this pass)
  // (clean_cs_lz_abs = 2: the L_z V R (dp_torque - dp) the gas does not take, booked as
  // the radiation's, cln_lzl (<rad_vet>/clean_lz_ledger); its sums over the cells as
  // the gas momentum's slopes K1 (cln_lzs1, <rad_vet>/clean_lz_ledger_k1), K2 (cln_lzs2)
  // and stage start increment (cln_lzc), and the stage-1 part pending until stage 2
  // (cln_lzp). 3-D spherical polar: the same ledger, whose cells' sum there is the L_z
  // that crossed the boundaries, the open or imposed radial ones (M1LzKey))
  Real cln_lzl = 0.0, cln_lzs1 = 0.0, cln_lzs2 = 0.0, cln_lzc = 0.0, cln_lzp = 0.0;
  ParameterInput *cln_lzpin = nullptr;
  void CleanLzBook(const Real s, const bool start);
  Real CleanLzTotal();   // (the history's L_z: gas + radiation + ledger)
  bool cln_xtnorm = false, cln_xtcons = false;   // (CleanFirstIterate, CleanSolve)
  int cln_xte = 0;   // (clean_xt_exact: 0 off, 1 cons, 2 x2)
  int cln_gnsw = 0, cln_gnswmin = 20;   // clean_gas_newton_switch, _min
  Real cln_stfac = 0.5;   // clean_stall_fac
  bool cln_gpic = false;   // (a stalled solve's passes without the gas Newton slope)
  Real cln_ngnsw = 0.0;   // (solves switched so far)
  int cln_fdg = -1;             // staggered bracket: face slot of its lagged flux DG
  bool cln_opt = false;         // clean_stag_opac_update: opacities at the iterate's T
  bool cln_consist = false;     // clean_consistency: bracket C = F_ray - f(Ur), cell
  bool cln_closit = false;      // clean_closure: predict | iterate
  int cln_extrap = 0;           // clean_closure = extrap (1) | extrap2 (2)
                                // | extrap1 (3, cell default)
  bool cln_edd = false;         // clean_closure = eddington (D = I/3, no rays)
  bool cln_sphere = false;      // clean_closure = sphere (spherical polar)
  Real cln_sph_r = 0.0;         // clean_sphere_radius
  bool cln_nl_cell = false;     // clean_nl_norm = cell: per-cell relative Newton changes
  Real cln_trtol = 0.0;        // clean_troot_rtol: residual stop of the exchange root
                               // (staggered default 1e-13)
  Real cln_nl_escale = 0.0;     // that E scale, clean_nl_efloor max E of the iterate
  int cln_wo_nc = 0;            // the slots of cln_w.o (CleanWoBind)
  static constexpr int cln_hmax = 16;   // most variables of one ghost exchange
  bool cln_budget = false;      // clean_conv = budget (staggered default): CleanSolve
  bool cln_bon = false;         // this solve takes it (Newton passes with a valid K1)
  Real cln_tolb = 0.0;          // its tolerance of the pass (kappa x time error)
  Real cln_bkap = 1.0;          // clean_conv_kappa
  Real cln_bmax = 1.0e-5;       // clean_conv_max
  DvceArray5D<Real> cln_wci;    // E before a backward-Euler solve's last pass
  Real cln_ewr = -1.0, cln_eta = 0.1;   // Eisenstat-Walker: last |r0|_S, forcing
  int cln_verbose = 0;          // clean_verbose: 1 = also the per-solve diagnostics
  int cln_stagmg = 0;           // the old production mg (staggered): 1 its line colours
                                //   (j + k) parity, 2 mirror (clean_stag_mg_color)
  int cln_pcsweep = 0;          // clean_pc_sweep: +1 a 2nd red pass, +2 a post-sweep
  bool cln_steady = false;      // clean_init = steady: radiation relaxed before step 1
  int cln_steady_n = 50;        // clean_init_n: most pre-solve passes
  Real cln_steady_dtfac = 1.0e6, cln_steady_tol = 1.0e-8;  // clean_init_dtfac, _tol
  std::vector<DvceArray5D<Real>> cln_mg;          // coarse levels of the multigrid
  std::vector<std::array<int, 3>> cln_mgn;        //   and their cells along x1..x3
  std::array<int, 3> cln_mgt = {0, 0, 0};       //   its tiles' cells (layout_tile)
  bool cln_fsolve = false;      // staggered: the face F of a separate F solve
  bool cln_fkin = false;        // staggered fsolve: streaming share of the rays' kinetic
                                //   momentum flux
  bool cln_dok = true;          // staggered fkin: the last solve converged (memory d)
  bool cln_fwork = false;       // staggered fsolve: also in the Newton loop (the work W)
  bool cln_fsk = false;         //   CLF_N holds the F solve's F of the iterate
  bool cln_vimp = false;        // staggered: implicit gas velocity (drag), per face
  bool cln_upw = false;         // staggered: donor-cell advective flux, lagged rest
  bool cln_ucp = false;         //   clean_stag_central
  bool cln_ucen = false;        //   clean_stag_upwind_central: the rest central
  bool cln_split = false;       // staggered: energy flux apart from the face momentum
  int cln_nr1 = 0, cln_nr2 = 0;  // staggered: cells E < 0, faces |F| > c E_donor
  int cln_ncwin = 0;   // solves accepted by the per-cell measure's stall window (cnorm)
  Real cln_fe = 0.0, cln_fv = 0.0;  // round-off levels of the Newton changes dE, dv
  bool cln_feedback = false;    // gas coupling active
  bool cln_init = false;        // arrays allocated
  bool cln_faces_set = false;   // starting face fluxes set
  bool cln_rstold = false;      // restart from a file without the clean state
  bool cln_havec = false;       // a closure has been built
  bool cln_cip = false;         // a closure is in place (clean_closure = iterate: the
                                // last iterate's, restart state)
  // opac holds predicted opacities during the formal solution
  bool cln_opsaved = false;
  Real cln_tc = -1.0e300;       // time of the closure build
  Real cln_dtb = 0.0;           // dt of the closure build
  DvceArray5D<Real> cln_s;      // closure (CLN_*)
  bool cln_cssc = false;        // the cubed sphere, clean_rays = sc (vet_gd sweep)
  bool cln_wfse = false;        // clean_w_fse
  bool cln_wfst = false;        // clean_w_fse_tau
  bool cln_hctau = false;       // clean_face_hctau
  bool cln_jse = false;         // clean_face_jse
  bool cln_farea = false;       // clean_face_area
  bool cln_bcarea = false;      // clean_bc_area
  DvceArray5D<Real> cln_fvi;    // the sc rays' intensities (m, d, k, j, i)
  DvceArray2D<Real> cln_fvd;    // the swept directions (x, y, z, weight; CleanFvSets)
  int cln_fvn = 0;              // their number
  std::vector<Real> cln_fvkey;  // the set (block key) in each slot of a rotating union
  // clean_rot_ring: the rotating union swept as the current and the next set only, the
  // K - 1 older sets carried as offsets of their ray moments (CleanRingMoments)
  int cln_ring = 0;             // clean_rot_ring = live (1: + the next set) | jump (2)
  bool cln_rfuse = false;       // clean_rot_fuse: live ring's two sets in one sweep
  bool cln_rfuse_said = false;
  bool cln_ring_named = false;  // clean_rot_ring named in the input
  bool cln_rlast = false;       // jump: the build closes its block (the next set swept)
  bool cln_rjump = false;       // jump: cln_rnew holds the next block's moments
  DvceArray5D<Real> cln_rnew;   // jump: the moments of block b+1's window at the boundary
  bool cln_rfix = false;        // no build of this cycle shifts (CleanRingPreShift)
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
  DvceArray5D<Real> cln_fh;     // the sc rays' face half-range sums (CleanFaceHalf)
  DvceArray5D<Real> cln_gdgl;   // their ghost layers (CleanGdFaceSums, scratch)
  // the slots' moments g = M/J (CLN_J: J); a build's scratch
  DvceArray5D<Real> cln_rl;
  DvceArray5D<Real> cln_rr;     // offsets g_{b-j} - g_b, J_{b-j}/J_b (< 0: none), K-1
  DvceArray5D<Real> cln_rpd;    // the pending Delta g_{b+1} - g_b, J_{b+1}/J_b
  Real cln_rguard_sum = 0.0, cln_roff_sum = 0.0;   // guard falls, off-boundary shifts
  DvceArray5D<Real> cln_v;      // closure fill, bracket, E < 0 flags (m,3,k,j,i)
  DvceArray5D<Real> cln_opb;    // opacities of the stage, kept during the closure build
  // clean_closure = extrap(2): ray moments (CLN_* layout, not normalised by J) of the
  // converged states at t^n (cln_mn), t^{n-1} (cln_mp), the end of the last step solved
  // (cln_mc) and the stage-1 solution (cln_m1), with their times and validity
  DvceArray5D<Real> cln_mn, cln_mp, cln_mc, cln_m1;
  Real cln_tn = -1.0e300, cln_tp = -1.0e300, cln_tm = -1.0e300, cln_t1 = -1.0e300;
  bool cln_have_n = false, cln_have_p = false, cln_have_c = false, cln_have_1 = false;
  bool cln_cand_be = false;     // the step that made cln_mc ended with backward Euler
  int cln_neint = 0;            // cells of the write-back left with gas eint < 0
  int cln_nguard = 0;           // cells whose extrapolated moments fell back to t^n
  Real cln_nguard_sum = 0.0, cln_next_sum = 0.0;   // guard falls, extrapolated closures
  CleanW cln_w;                 // work vectors (CLW_*)
  DvceArray5D<Real> cln_fc1, cln_fc2, cln_fc3;   // face variables (CLF_*)
  DvceArray5D<Real> cln_fref;   // force_reference = clean: reference force (m,3,k,j,i)
  // staggered: the cell flux of the faces' energy flux Phi (m,3,k,j,i; the write-back,
  // output m1 as m1_phi), u0 F holding the face momenta's (the radiation momentum)
  DvceArray5D<Real> cln_phic;
  DvceArray5D<Real> cln_rstd;   // staggered fkin: the restart's d and its K1 (m,6,k,j,i)
  bool cln_fref_set = false;    // cln_fref holds the force of the first coupled solve
  bool cln_mdcap = false;       // wb_arad_md: the write-back captures cln_fref (steady)
  bool cln_vzero = false;       // the iterates keep the gas momentum (steady start)
  // cubed sphere: face areas, centroid distances across the faces, normals and centroids
  DvceArray4D<Real> cln_ga[3], cln_gd[3];
  CleanGeo cln_gn;
  DvceArray3D<Real> cln_tg;   // spherical polar: CleanArgs::tg
  struct CleanHalo {            // ghost exchange of nc components: scratch array and
    DvceArray5D<Real> a, ca;    // its coarse dummy
  };
  std::map<int, CleanHalo> cln_halo;   // by nc
  DvceArray1D<char> cln_xdev;   // device copies of kernel arguments (CleanDevCopy)
  // the exchange object of every clean ghost fill (clean_halo_vars variables; a fill
  // sends only its own), the fills deferred to the next exchange (CleanGhostDefer)
  std::unique_ptr<MeshBoundaryValuesCC> cln_hpb;
  std::vector<CleanGSeg<Real>> cln_hpend;
  DvceArray5D<Real> cln_hbuf;   // the storage of every scratch array (as used, <= 32)
  // clean_halo_direct (only true): the exchange in old production's halo_mpi form,
  // cell lists read once off cln_hpb's own exchange (CleanHaloInit): on-rank ghosts by
  // one copy, one message per neighbour rank, receives posted before the pack
  struct CleanHaloList {
    int state = 0;              // 0 lists to make at first use, 1 on, -1 while made
    int n1 = 0, n2 = 0, n3 = 0, nl = 0, ns = 0, nr = 0, par = 0;
    DvceArray1D<int> ld, ls, sc, rc;   // on-rank ghost <- its source; sent; received
    std::vector<int> rank, soff, slen, roff, rlen;
    // (the lists in order of depth, the ghost layer, 1..ng; per message
    // and depth d the entries of depth <= d, slend/rlend[s ng + d - 1], nld[d - 1];
    // rd, ldp the depth of each received / on-rank ghost; dcur the exchange's depth)
    int ng = 0, dcur = 0;
    bool poison = false;
    Real pval = std::numeric_limits<Real>::quiet_NaN();
    std::vector<int> slend, rlend, nld;
    DvceArray1D<int> rd, ldp, sdp;   // (sdp: the depth of each sent entry)
    // (device buffers on HIP too: the messages are MB-size, ROCm-aware MPI moves them
    // GPU to GPU; host-pinned staging is CPU-copy bound on MI300A)
    DvceArray1D<Real> sbuf, rbuf;      // (rbuf: two halves, alternate exchanges)
    // (CleanGSeg::fm: the family of each on-rank, sent and received ghost (hsf, hrf on
    // the host); per mask set (its fm over all components, nc) each sent and received
    // entry's buffer offset and the messages' lengths by depth)
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
  CleanHaloList cln_hdl;
  void CleanHaloInit(const int n3, const int n2, const int n1);
  template <typename T> void CleanHaloPost(const std::vector<CleanGSeg<T>> &sg);
  template <typename T> void CleanHaloFinish(const std::vector<CleanGSeg<T>> &sg);
  template <typename T>
  std::array<uint32_t, 5> CleanHaloKey(const std::vector<CleanGSeg<T>> &sg);
  CleanHaloList::Mask &CleanHaloMask(const std::array<uint32_t, 5> &key);
  Real cln_nsolve = 0.0, cln_nit = 0.0, cln_nfail = 0.0, cln_nneg = 0.0;
  //! exchange the x1 ghost layers the multi-block solve needs; `eponly` picks the
  //! one-quantity set (E of the new iterate) instead of the six lagged ones
  //! gather the assembled rows onto the roots, run Thomas there, scatter back
  //! milestone 3b phase B: the lagged transverse (x2/x3) divergence of one Picard pass,
  //! split into its diagonal part (M1_IW_TDIA) and its right-hand side (M1_IW_TRHS)
  //! the TRANSVERSE realizability limiter: fill thx2/thx3 (and, when newk, the lagged
  //! limiter opacity klx2/klx3) for the current pass.  Inert under
  //! implicit_trans_limit = none, where thx2/thx3 are not even allocated.
  //! milestone 3e: snapshot the SCALED state entering a Picard pass (x_k) into aa_xc.
  //! Called at the top of every pass when implicit_accel = anderson.
  //! milestone 3e: at the end of Picard pass `it` the iw state is G(x_k).  Form the
  //! fixed-point residual g_k, push (dX, dG) onto the history, solve the small least-
  //! squares problem and overwrite the iw state with the accelerated iterate x_{k+1}
  //! (realizability re-applied).  The statistics go into aa_nacc / aa_nrst.
  //! milestone 3b phase C: the x1 line solve of the assembled rows (M1_IW_TA..TR ->
  //! M1_IW_S2), Thomas or cyclic Thomas, gathered over the stack when part_nblk > 1.
  //! It IS the preconditioner of the BiCGStab wrapper and the line-Jacobi pass itself.
  //! implicit_line_solver = pcr: the same line solve by parallel cyclic reduction,
  //! one Kokkos team per x1 column (GPU)
  //! implicit_krylov_fuse: the pcr line solve with the right-hand side read from `rc`
  //! (upd = 0), or made on the fly as the BiCGStab p (upd = 1) / s (upd = 2) update,
  //! and the answer written to `zc`.  `col` >= 0 solves only the columns of that
  //! red-black colour (implicit_precond = rbgs); `sub` >= 0 then subtracts the
  //! transverse 5-point coupling to component `sub` from the right-hand side.
  //! mgrows = true reads the rows from mgr (implicit_mg_full_rows) and
  //! subtracts the other-colour edge couplings of the stencil with `sub`
  //! implicit_od_cache: fill odc with sum_{e!=d} d_e P_de(x) of component xc
  //! the preconditioner of the fast path: z = M^{-1} r (see implicit_precond)
  //! implicit_od_cache: y += sgn L_off(x) (with7 = false) or y = A x (with7 = true) from
  //! the cache; red > 0 also returns reductions in out[4] (see the definition)
  //! implicit_krylov_fuse = 3: BiCGStab with TWO blocking reductions per iteration
  //! implicit_op_stencil: build the 19-point operator of the pass / apply it
  // red_only (m1-fast3): y is already in yc (the overlap parts): the reductions only
  //! the fast-path operator product: stencil or od cache
  //! the Thomas / cyclic-Thomas sweep (implicit_line_solver = thomas)
  //! milestone 3b phase C: y = A x with the frozen 7-point operator (one halo of x)
  //! milestone 3b phase D: accumulate sgn * L_off(x) into the component yc, with
  //! L_off the cell-row contribution of the OFF-DIAGONAL Eddington terms of every face
  //! equation, evaluated at the component xc with the closure of this Picard pass
  //! FROZEN (so it is linear in x).  No communication: the caller has already put x in
  //! the ghost zones (ImplicitApplyOp) or x IS the iterate, whose halo is current.
  //! milestone 3b phase C: z = M^{-1} r, the x1 line solve applied to an arbitrary
  //! right-hand side (it overwrites M1_IW_TR and M1_IW_S1..S3).  rc < 0: the caller
  //! has already staged r in M1_IW_TR.
  //! milestone 3b phase C: solve the frozen 7-point system by BiCGStab; the answer is
  //! left in M1_IW_S2, exactly where the line-Jacobi pass leaves it.  Returns the number
  //! of inner iterations taken.
  // returns 0: go on (implicit_final_resid_stop failed, this pass continues), 1: accept,
  // 2: failed with no pass left (non-converged); 1 and 2 put the saved state back
  //! the same BiCGStab with fused reductions and vector updates (implicit_bcg_sync > 0)

  // ---- <rad_vet>/closure = tau (rad_vet_tau.cpp): the Eddington tensor of the multi-D
  // implicit solve from the COLUMN OPTICAL DEPTH measured from the top of the domain
  // (x1max): chi = exact grey plane-parallel K/J(tau) (Hopf fit), axis along -grad tau.
  // Never reads the local flux.  Rebuilt once per hydro step from the start-of-step
  // transport opacity; x1 stacks of MeshBlocks across MPI ranks are supported.
  bool tau_closure{};          // = vet_col: the column tensor tau_ten is in use
  bool tau_ready{};            // VetColInit has run
  DvceArray5D<Real> tau_ten;   // (m,4,k,j,i): chi, n1, n2, n3 read by ImplicitSolve (b)

  // ---- <rad_vet>/closure = vet_col (rad_vet_vetcol.cpp, design stage S5, option D): the
  // Eddington factor f_K = K/J of a 1-D formal solution per RADIAL COLUMN (spherical
  // impact-parameter rays on the sp wedge, Gauss rays in plane-parallel on a Cartesian
  // mesh) with the column's own extinction and source, handed to the implicit solve as a
  // FIXED uniaxial tensor about r_hat (or the M1 flux axis) through the tau closure's
  // tau_ten (tau_closure is set as well; only the tensor build differs).
  bool vet_col{};              // closure = vet_col (default false)
  bool vcol_sph{};             // spherical rays (sp) or plane-parallel (Cartesian)
  bool vcol_axis_flux{};       // vet_col_axis = flux: n = the M1 cell flux direction
  int vcol_nc{}, vcol_np{}, vcol_nmu{}, vcol_every{}, vcol_nray{};
  // vet_col_fk_min (m1-sp-order2b): the lower bound of f_K = K/J.  1/3 (the old clamp,
  // inherited from vet_sc's eigen/flux projection) or 0 (sp default: f_K < 1/3 is
  // physical behind a radiating front and D = diag(f, (1-f)/2, (1-f)/2) stays
  // realizable for 0 <= f <= 1)
  Real vcol_fkmin = 1.0/3.0;
  // vet_source = relaxed | relaxed_esc (vcol_sesc) | gas, vet_col_source its alias: the
  // relaxed ray source was removed (the clean rays take the state's own source), the
  // key is only checked and echoed
  bool vcol_srelax = true;
  bool vcol_sesc = false;
  Real vcol_sx0 = 1.0;         // time2_vet_col_source_x0 (blend)
  // vet_scatter (read only when named; default off = bitwise): ONE
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
  // sweep
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
  bool vcol_built{};           // a tensor exists (vet_col_every > 1 skips builds)
  Real vcol_time{}, vcol_ncall{}, vcol_nskip{};
  // vet_col_surface_q (runs_5e): the outer-x1 Marshak q per column from the formal
  // solution, q = H(top face)/J(top cell), lagged like the tensor
  bool vcol_sq{};
  Real vcol_qmin{}, vcol_qmax{};
  int vcol_lc{};               // shells per chunk of the team build (scratch budget)
  DvceArray3D<Real> vcol_q;    // (m, k, j): the column's Marshak q at the top face
  void VetColInit();           // checks, tau_ten, vcol_q (first call)
  void VetColReport();

  // ---- <rad_vet>/vet_col_lat (m1-vetcol-lat, rad_vet_vetlat.cpp; design
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
  // the Athena++ lateral ratio keys <mesh>/x2rat, x3rat named != 1 in
  // the input (AthenaK ignores them: its theta/phi grid is uniform); VetLatInit refuses
  // them under vet_col_lat / vet_gd, so a ported deck cannot think it stretched theta
  std::string vlat_rat_keys;
  // vet_col_lat_offdiag: what D_r,lat does in the solve.  0 none, 1 lagged (a Picard
  // term, C1), 2 operator (default, C2: in the Krylov operator, M1SphLat of the
  // Krylov vector, through the stored stencil and VetLatOp)
  int vlat_odm = 2;
  bool vlat_now = false;       // the term is on in this step (positivity fallback)
  Real vlat_odmax = 0.0;       // max |D_r,lat| of the last sweep (all ranks)
  // the operator form is active: on in this step and a non-zero D_r,lat (a zero one
  // leaves the operator, its stencil and the right-hand side untouched, bitwise)
  bool VlatOp() const {
    return vlat_now && !vlat_lagfb && (vlat_odm == 2) && (vlat_odmax > 0.0);
  }
  // vet_col_lat_posfix = lag (read only when named; default drop =
  // bitwise): at a positivity trigger the lateral D_r,lat and Q terms leave the OPERATOR
  // (which is then the M-matrix form) but stay in full as the lagged Picard term of the
  // face equations, so the fixed point keeps the whole tensor.  Only a second trigger in
  // the same step drops them (vlat_nfall, as before)
  bool vlat_posfix_lag = false;
  bool vlat_lagfb = false;     // this step: the lateral terms lagged after a trigger
  // (rad_vet_meshfb.cpp): implicit_fb_diag, implicit_fb_local
  bool fbl_any = false;
  DvceArray5D<Real> vlat_tt_loc;
  DvceArray5D<Real> VlatTT() const {return fbl_any ? vlat_tt_loc : tau_ten;}
  // theta of the data each x2 row holds, for the trigonometric factors of
  // M1SphTan (= x2v; the periodic source cell's theta in the ghost rows at a
  // theta-periodic seam).  Built on first use, rebuilt when the pack size changes.
  const Real *vtan_x2t_src = nullptr;
  static constexpr int vlat_nmu = 4, vlat_npsi = 4;   // vet_col_lat directions
  Real vlat_taucut = 30.0;     // vet_col_lat_taucut: the sweep starts where every column
                               // has tau_top >= this (global shell index vlat_icut)
  int vlat_icut = 0;
  Real vlat_cut_thin = 0.0;    // vet_col_lat_cut_thin: no cut above a
                               // cell thinner than this (0 off)
  bool vlat_ready = false;     // VetLatInit done
  int vlat_nbuild = 0;         // builds done
  Real vlat_time = 0.0, vlat_ncall = 0.0, vlat_nclamp = 0.0;
  DvceArray5D<Real> vlat_i;    // (m, 2 nmu npsi, k, j, i): 3-D intensities
  DvceArray5D<Real> vlat_t;    // (m, 2 nmu npsi, k, j, i): twin intensities
  DvceArray5D<Real> vlat_t_c;
  DvceArray5D<Real> vlat_d;    // (m, 2 nmu npsi, k, j, i): 3-D minus twin, ghosts
                               // lagged (the lateral part of the inflow)
  DvceArray5D<Real> vlat_d_c;
  DvceArray5D<Real> vlat_cs;   // (m, 2, k, j, i): ln chi, ln S (ghosts exchanged)
  DvceArray5D<Real> vlat_cs_c;
  DvceArray3D<Real> vlat_geo;  // (shell l, 2 a + dir, 0..3): per-shell ray geometry
  DvceArray1D<Real> vlat_mu;   // (a): Gauss nodes on [0, 1], descending
  DvceArray1D<Real> vlat_w;    // (a): hemisphere weights, sum 1/2 (J = sum w I)
  DvceArray1D<Real> vlat_cnt;  // (1): lateral reads clamped to the ghost band
  std::vector<double> vlat_nodes;   // the mu nodes (host)
  int vlat_geo_icut = -1;       // the first shell vlat_geo was built for
  MeshBoundaryValuesCC *pbval_vl = nullptr;   // vlat_i
  MeshBoundaryValuesCC *pbval_ea = nullptr;   // EAUDITF diagnostic
  MeshBoundaryValuesCC *pbval_vs = nullptr;   // vlat_cs
  void VetLatInit();
  void VetLatBuild();
  void VetLatSweep(const int stage);   // 0 twin, 1 3-D, 2 moments
  void VetLatOdMax();
  void VetLatExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                      MeshBoundaryValuesCC *pb);
  void VetLatDump();
  // vet_col_lat_offdiag = operator: y += sgn L_lat(x) (the direct form, legacy
  // operator path and the right-hand side), and the same terms added to the stored
  // 19-point stencil of the pass (ImplicitStencilBuild)
  // the cubed-sphere forms (M1CsLat0 / M1CsLatF, Cartesian D_r,lat)

  // ---- <rad_vet>/vet_gd (m1-vet-gd, rad_vet_vetgd.cpp; design
  // rt_design_1003/SC_PROPER_SP.md, candidate b): the Eddington tensor from an SC sweep
  // along GLOBALLY FIXED (Cartesian) HEALPix directions, handed over through the
  // vet_col_lat interface (sets vlat_on; the vet_col_lat keys apply).  Read only when
  // named; sp wedge only.
  bool vgd_on = false;         // vet_gd
  // HEALPix nside of the set that names none: 12 nside^2 = 108 directions (only on old
  // restarts / recipe none)
  static constexpr int vgd_nside = 3;
  // vet_gd_tangential (read with vet_gd): spherical polar the lagged M1SphTan term; the
  // cubed sphere Q in LAT4, LAT5 about the per-cell tangent frame cs_qf (m, 6, k, j incl.
  // ghosts; M1CsQab), M1CsTanFaceLag / M1CsTanFaceOp
  bool vgd_tan = true;
  DvceArray4D<Real> cs_qf;
  void CsTanFrameInit(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                      MeshBoundaryValuesCC *pb);
  bool vgd_replace = false;    // vet_gd_replace: no vet_col sweep (D = I/3 below the
                               // first gd shell, q from the gd top cell)
  int vgd_n = 0;
  int vgd_ls = 0;              // vet_gd_ls = 6/8/10/12: level-symmetric LQ_N set instead
  double vgd_zper = 0.5*M_PI;  // the set's z-rotation period (VetGdInit)
  int vgd_gl_nmu = 0;          // vet_gd_gl_nmu > 0: product set GL(nmu in n_z) x nphi
  int vgd_gl_nphi = 12;        // vet_gd_gl_nphi (C_nphi about z)
  bool vgd_gl_stag = true;     // vet_gd_gl_stagger: alternate rings half a step
  int vgd_rot = 0;             // vet_gd_rotate_every: z-rotation every N cycles
  // vet_gd_rotate_dt > 0 (read only when named): the rotation in PHYSICAL
  // time instead (block b of the step = floor((t + dt/2 - t0)/dt_rot), the same angle
  // at the same t on every dt), vet_gd_rotate_t0 its origin; vet_gd_rotate_align (default
  // true): Mesh::NewTimeStep lands a step on every block boundary (VetGdAlignDt)
  Real vgd_rot_dt = 0.0, vgd_rot_t0 = 0.0;
  // vet_gd_rotate_stretch: VetGdAlignDt may lengthen dt by <= this
  Real vgd_rot_stretch = 0.0;
  // vet_gd_rotate_dt = auto (the fresh vet_gd default): resolved by the
  // first VetGdAlignDt to vgd_rot_auto_fac x that dt and recorded in pin_rot_
  bool vgd_rot_auto = false;
  // the unnamed fresh default resolves from the first CFL dt before the
  // dt_fixed / dt_retry cap (Mesh::dt_uncapped), so a fixed-dt ladder shares one clock
  bool vgd_rot_unc = false;
  Real vgd_rot_auto_fac = 20.0;
  ParameterInput *pin_rot_ = nullptr;
  bool vgd_rot_align = false;
  double vgd_rot_nalign = 0.0;   // steps shortened by VetGdAlignDt (diagnostic)
  // vet_gd_rotate_slope (with vet_gd_rotate_dt; be | keep | correct): the hesdirk2 FSAL
  // slope K1 of the first step of a new rotation block belongs to the OLD set (an O(dt)
  // stage-1 error per jump = first order); be: that step is backward Euler, keep: as
  // before (be measured WORSE than keep, vgdrot2 NOTE: be is Lie split); probe (a probe
  // solve of the new set) is removed
  int vgd_rot_slope = 0;
  int64_t vgd_rot_lastblk = INT64_MIN;
  void VetGdRotStep();
  // vet_gd_rotate_slope = correct (3): on the first step of a new block a
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
  // (USER decision 1): <rad_vet>/vet_gd_rotate_axis = z | 3d (default z on
  // the sp wedge: bitwise; 3d on the cubed sphere).  3d: the set of rotation block b is
  // R_b n with R_b a 3-D rotation about a quasi-random axis (VgdRot3); vgd_alpha is then
  // the block KEY (zig-zag code of b, or of b mod K with vet_rotate_avg = K), not an
  // angle
  bool vgd_rot3 = false;
  void VgdRot3(const Real key, double r[9]) const;
  std::vector<double> VgdRotate(const Real alpha) const;
  // <rad_vet>/vet_rotate_avg = K >= 2 (read only when named; default 0 = off,
  // bitwise): the closure handed to the moment step is the mean over the LAST K rotated
  // direction sets, not the current set's alone (rad_vet_rotavg.cpp).  Closure-agnostic:
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
    int arr[M1_RA_NF]{}, comp[M1_RA_NF]{}, kind[M1_RA_NF]{};
    Real lo[M1_RA_NF]{}, hi[M1_RA_NF]{};
    int farr = -1, fcomp = 0;  // validity flag field: invalid where < 0 (-1: none)
    int jf = -1;               // the J field (index into the table) of kind-3 fields
    int jw = -1;               // the ratio field weighting the sets (-1: equal)
  };
  int vra_k = 0;               // vet_rotate_avg (K; 0/1 off)
  // <rad_vet>/vet_rotate_avg_ema = N >= 1 (default 128 on a fresh run whose
  // set rotates, else read only when named; 0 = off, the K window): the K-window mean
  // replaced by an exponential moving average over the rotated sets, weight a = 1/N per
  // rotation block, kept in the same offset form (ONE
  // slot: e = avg - g_b(now), r = Jbar/J_b, + an aux slot for the change test; vra_k = 3
  // internally): effective window
  // ~2N - 1 sets at the memory of K = 2.  vet_rotate_avg_ema_tau > 0 (physical time) and
  // vet_rotate_avg_ema_lc > 0 (x the local light-crossing time r/c) tie a to a time
  // window instead: a = max(1/N, 1 - exp(-dt_blk/tau)).  The rotation sequence does not
  // repeat (3d: the R_3 sequence of b, z: golden phases)
  bool vra_ema = false;
  Real vra_ema_n = 0.0, vra_ema_tau = 0.0, vra_ema_lc = 0.0;
  void RotAvgEmaInit(const int ilo);
  void RotAvgEmaShift(const int ds);
  void RotAvgEmaChange(const int ilo);
  Real vra_ema_eps = 0.05;     // vet_rotate_avg_ema_change (0: no change test)
  bool vra_ema_chj = true;     // vet_rotate_avg_ema_change_j: the J ratio in the test
  bool vra_ema_chstep = false; // vet_rotate_avg_ema_change_form = step (default) | linear
  // vet_rotate_avg_ema_wb (read only when named, default false = bitwise; clean
  // ring only): a = a_ch alone, no 1/N floor and no 1/(n + 1) count, so the pattern
  // average moves only where the field changed and a static state keeps its closure
  bool vra_ema_wb = false;
  int vra_ema_diag = 0;        // vet_rotate_avg_ema_diag = M: print the change weights
  // vet_rotate_avg_ema_dema (read only when named, default false): the lag-
  // corrected (double) EMA: slot 2 holds the EMA of the EMA (e2, same weights, same frame
  // shifts), the closure uses 2 e1 - e2 (ratio fields 2 r1 - r2, kept in [r1/4, 4 r1]);
  // vra_k = 4 (slot 3 the Delta)
  bool vra_ema_dema = false;
  int64_t vra_ema_nsh = 0;     // the boundaries seen (diag cadence only)
  void RotAvgEmaApply(const int ilo);
  // <rad_vet>/vet_gd_rotate_xfade = true (read only when named; default false; needs
  // vet_gd_rotate_dt > 0): the clean ring live (clean_rot_ring = live), the closure
  // cross-faded from the current to the next set over each rotation block
  // (CleanFvSets); the two sets' tables stay on the device (VgdXfTables)
  bool vgd_xf = false;
  // vgd_xff (with vet_gd_rotate_xfade): the two sets of a cross-fade build swept in
  // ONE sweep: per shell one kernel launch over both sets' (cell, direction) lists and
  // (sp, MPI exact lists) one message per partner carrying both sets' values.  Every
  // value is computed by the same code from the same inputs as by the two sweeps, and the
  // moments / blend run in the same order afterwards (bitwise).  The second set keeps its
  // own intensities (vgd_i_xb, the size of vgd_i: the extra memory), shell lists and
  // exact halo lists; they travel with the set's tables (VgdXfTables).  No band ring
  // (vgd_ring off: the ragged array).  Falls back to the two sweeps where
  // a combination is not covered (VetGdXffOk)
  bool vgd_xff = false;
  bool vgd_xff_now = false;      // this sweep carries the next set too (VetGdSweep)
  bool vgd_xff_prep = false;     // VetGdSweep: the set's lists only, no passes
  DvceArray1D<int> vgd_shl_xb[2];
  int vgd_shn_xb[2] = {0, 0};
  Real vgd_shl_alpha_xb = -3.0;
  DvceArray1D<char> vgd_bbufx;     // the fused shell closure (device copy)
  DvceArray1D<int> vgd_cslst_xb;   // cs: the other set's pass lists
  DvceArray1D<Real> vgd_xff_sb, vgd_xff_rb;   // the merged messages (sp)
  double vgd_xff_n = 0.0;          // fused builds (diagnostic)
  double vgd_xff_nmsg = 0.0;       // merged exchanges (diagnostic)
  double vgd_xff_nfb = 0.0;        // builds that fell back to two sweeps (diagnostic)
  bool VetGdXffOk();
  void VgdXffSwap();
  void VgdEs2Make();   // the cubed sphere's early send's second stream
  void VgdGbSwap() {   // the compact global band's kept lists go with their set
    std::swap(vgd_gbcf, vgd_gbcf_xb);
    std::swap(vgd_gbcb, vgd_gbcb_xb);
    std::swap(vgd_gbrec, vgd_gbrec_xb);
  }
  bool VetGdHlXff(const int i, const int ws);
  bool VetGdGbHaloX(const int i, const int ws);   // (vgd_raw: both sets, one message)
  // early send of the cubed sphere (VetGdCsSweep): this shell's SEND cells first, on a
  // second high-priority stream, with their pack; the host exchanges while the shell
  // kernel (which skips the stamped cells) runs
  DevExeSpace vgd_es2;
  bool vgd_es2_made = false;
  DvceArray4D<int> vgd_est[2];     // per set (m, kb, jb, d): the stamp of the early cells
  int vgd_esid = 0;                // the stamp of the current early shell
  Real vgd_xf_alt = -9.0;        // the angle of the tables in the second device set
  void VgdXfTables(const Real alpha);
  bool vra_eq = false;        // vet_rotate_avg_weight = equal (default j: J-weighted)
  bool vra_on = false;
  RotAvgTab vra_tab;
  DvceArray5D<Real> vra_a0, vra_a1;   // the closure's field arrays (set per call)
  DvceArray3D<Real> vra_c;     // a column field (m, k, j) at the top cell (arr 2)
  DvceArray5D<Real> vra_d;     // (m, (K-1)(nf+1), k, j, i) interior: offsets + valid
  int64_t vra_gen = 0;        // bumped by every writer of vra_d
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
  // the history in the restart file (kM1RaRstMagic, rad_vet_rotavg.cpp)
  static constexpr int M1_RA_RST_NH = 11;   // int64 header entries (10)
  int RotAvgNf() const;
  int RotAvgRstNslab(const int ncell) const;
  int RotAvgRstI0() const;     // the first written radial cell
  void RotAvgRstHdr(int64_t *h) const;
  void RotAvgRstPack(DvceArray5D<Real> &dst, const int nmb);
  void RotAvgRstPackHost(HostArray5D<Real> &out, const int nmb, const int nsc);
  bool RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &src, const int nmb);
  bool RotAvgRstCheck();
  int vra_rst_nf = -1;         // nf of a history loaded from the file (-1: none/checked)
  int64_t vra_blk = INT64_MIN; // the block of the last build
  bool vra_last = false;       // this build is on the last step of its block
  DvceArray5D<Real> vra_bk_vc;              // vet_sc: vet_cell kept over the next-set
                                            // sweep of a last build
  void VetGdSweepPost();       // twin, sweep(s), fused moments, VetGdPost
  // vet_gd_inner = emit (read only when named; default diff = bitwise): with the first gd
  // shell at the physical inner face (icut = 0) the outward rays start from the inner
  // boundary's own emission (flux: black sphere I = 4 F_in/c; marshak: the bath) instead
  // of the diffusion intensity
  bool vgd_tbath = false;   // the top (x1max Marshak) bath enters the sweep
  std::vector<double> vgd_base;   // the unrotated set (x, y, z, w)
  DvceArray2D<Real> vgd_dir;   // (d, 0..3): n_x, n_y, n_z, weight (sum 1)
  // (m, d, k, j, i): intensities, ghosts lagged.  LayoutLeft (i slowest), so
  // that one shell of the per-shell sweep and halo is one contiguous ~nd*nk*nj slab (the
  // LayoutRight array spread every shell over all of its GBs: strided, TLB-bound); the
  // values and every arithmetic are unchanged (bitwise)
  // (memory): RAGGED per-shell band.  Shell i stores only the lateral band
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
  VgdIView vgd_i_xb;   // the other cross-fade set's intensities
  void VgdRagAlloc(VgdRag &a);
  // the memory arena (rad_vet_memarena.cpp; with vet_gd): the intensity arrays and solve
  // scratch that is dead across a VetGdBuild in one allocation
  bool marena_on = false, marena_sealed = false;
  DvceArray1D<Real> marena;   // the intensities, then (first checkpoint) solve scratch
  int64_t marena_sw = 0;      // Reals of the intensity arrays bound so far
  void MemArenaInit(ParameterInput *pin, int nmb, int niw, int n3, int n2, int n1);
  bool MemArenaVgd(VgdRag &a, int64_t n);
  void MemArenaGrow(int64_t need, bool keep = true);
  void MemArenaSeal();
  void CleanWoBind();   // cln_w.o, .s at the first solve (mem_arena)
  void MemcutPoison(int where);
  // vet_gd_band_ring (where it applies: MPI band with the exact halo lists, not the sp
  // raw band; on the cubed sphere with the cs compact lists, vet_gd_halo_compact > 0):
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
  int vgd_rg_k = VGD_RING_K;          // the slots of this allocation
  struct VgdRingOp {
    int kind = 0;      // 1 save (inward band -> store), 2 release, 3 acquire
    int i = -1, slot = -1;
    bool copyin = false, restore = false;
  };
  bool vgd_ring = false;              // the ring layout is allocated (both arrays)
  bool vgd_ring_named = false, vgd_ring_par = true;   // vet_gd_band_ring (named)
  bool vgd_ring_act = false;          // the ring runs in this sweep
  bool vgd_ring_poison = false;       // env VGD_RING_POISON=1: NaN in lent slabs (test)
  bool vgd_rg_two = false;            // this sweep carries a second array (vgd_i_xb)
  bool vgd_rg_xb = false;             // ... the other xfade set's array
  std::vector<int64_t> vgd_rg_sn2;    // (i) the other set's saved count
  bool vgd_rg_lastx = false;          // the last halo used an exact list
  int vgd_rg_lastx_e = -1;
  std::vector<int64_t> vgd_rg_offc;   // (m c1 + i): compact slab offsets (host)
  std::vector<int64_t> vgd_rg_soff;   // (s nmb + m): slot offsets (host)
  std::vector<int> vgd_rg_wfull;      // (i): a shell's band depth when lent a slot
  std::vector<int> vgd_rg_state;      // (i): -1 compact, else slot
  std::vector<int> vgd_rg_sown;       // (s): the shell holding slot s, -1 free
  std::vector<int64_t> vgd_rg_so, vgd_rg_sn;   // (i): save store offset, count/array
  std::vector<int> vgd_rg_sk;         // (i): save kind: 1 exact list, 2 dense band, 3 raw
  std::vector<int> vgd_rg_se;         // (i): the exact list record (3: the raw key)
  // vgd_raw: the band positions (encoded as the exact lists) each raw exchange writes,
  // keyed as vgd_gbcf (rad_vet_vetgd_sph2.inc VetGdRawRec)
  std::map<int, DvceArray1D<int>> vgd_gbrec, vgd_gbrec_xb;
  int VetGdRawKey(const int i, const bool inw, const int ws) const;
  DvceArray1D<int> VetGdRawRec(const int key);
  int64_t vgd_rg_stop = 0;            // store fill this sweep
  DvceArray1D<Real> vgd_rg_store;     // saved inward band values (main, then twin)
  Kokkos::View<int64_t**, DevMemSpace> vgd_rg_doffc;   // (m, i)
  Kokkos::View<int64_t**, DevMemSpace> vgd_rg_dsoff;   // (s, m)
  std::vector<VgdRingOp> vgd_rg_pend;
  Real vgd_rg_nk = 0.0, vgd_rg_ndense = 0.0;
  void VetGdRingFlush();
  void VetGdRingSweepBegin(const bool two);
  void VetGdRingSweepEnd();
  int VetGdRingFreeSlot(const int avoid) const;
  Kokkos::View<int ***, LayoutWrapper, DevMemSpace> vgd_wall;    // (m, k, j)
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_map;    // (m, k, j, d)
  DvceArray4D<Real> vgd_mr;    // (m, k, j, d): n . r_hat of the value stored there
  // vet_gd_wall_interp (read only when named, default false): wall ghosts from the 3
  // nearest same-branch directions with weights (m, k, j, d, 3)
  bool vgd_wint = false;
  // vet_gd_wall_sym (read only when named, default false; needs vet_gd_wall_interp):
  // the set turned so that a mirror plane is the wedge mid-plane, phi-wall ghosts linear
  // in azimuth on their ring (VgdMirPick): phi-symmetric on any wedge width
  bool vgd_wsym = false;
  bool vgd_wsx = false;       // and its phi wall map inexact: both-frame seam faces
  Real vgd_wsw = 0.0;         // the phi period W
  DvceArray1D<int> vgd_mir;   // (n): the mirror image of each direction
  // vet_gd_theta_seam_sym (read only when named, default false; needs
  // vet_gd_wall_interp): the clean scheme's theta seam faces in both frames (paired)
  bool vgd_wst = false;
  DvceArray1D<int> vgd_mirz;  // (n): the equatorial mirror image of each direction
  VgdLG vgd_lg;   // the vet_gd sweep geometry from the global index
  Kokkos::View<int *****, LayoutWrapper, DevMemSpace> vgd_m3;
  DvceArray5D<Real> vgd_w3;
  // vet_gd_theta_wall (read only when named; periodic | mirror, default
  // periodic = the old local-frame theta periodicity): mirror makes the theta walls
  // specular mirrors of the formal solution.  A theta ghost (band row beyond a wall at
  // depth g) reads the interior row at depth g on the other side of the wall, direction
  // (n.r, n.theta, n.phi) -> (n.r, -n.theta, n.phi) in the local frames (the exact
  // equatorial plane: z -> -z), from up to 3 same-branch set directions (the ghost's
  // own wall tables, vgd_wall = 2, vgd_map / vgd_m3 / vgd_w3; nearest only without
  // vet_gd_wall_interp); ln chi, ln S of the ghost rows are the mirror rows'
  // (VetGdThMirCs, from VetGdHaloCs)
  bool vgd_thmir = false;
  bool cln_scmem = false;      // clean sc rays, reflecting top: its mirror reads cln_top
  DvceArray5D<Real> cln_top;   // its last top shell: slot s, direction d at (s*nt + d/n1,
                               // i = d%n1), nt = ceil(n/n1) (VetGdTopCopy)
  void VetGdTopCopy();
  int cln_rq = 0;               // clean_rot_ring, sc rays: the slot swept (CleanFvSets)
  int cln_rqm = 0;              // and its sp top-mirror slot (blk mod 2)
  void VetGdThMirCs(DvceArray5D<Real> &a, const int i0, const int i1);
  // vet_gd_theta_wall = vacuum (1) | extrap (2), 0 otherwise: the theta
  // ghosts beyond a PHYSICAL theta wall (wall flag 3 / 4) are not the far wall's.
  // vacuum: zero incoming intensity at the wall (a ray whose upwind point lies beyond
  // the wall starts at its wall crossing with I = 0, as the vacuum top face; ghost reads
  // of the bilinear stencil are 0); extrap: a ghost reads the edge row's intensity of
  // the same direction (zero gradient).  ln chi, ln S of the ghost rows = the edge row
  // (VetGdThVacCs).  Swept by the wall kernel (bodym) with the mirror
  int vgd_thvac = 0;
  void VetGdThVacCs(DvceArray5D<Real> &a, const int i0, const int i1);
  MeshBoundaryValuesCC *pbval_gd = nullptr;
  // EXACT per-shell lateral halo (decomposition-invariant sweep): vgd_i and vgd_cs carry
  // a lateral ghost band of vgd_w cells (the deepest reach of an upwind point), filled
  // shell by shell from the 8 lateral neighbours (on-rank copies + MPI)
  int vgd_w = 0;
  // the clean closure on a wedge: the band of any depth keyed on the global cell
  // (VetGdGlobalBandInit, theta and phi periodic), the intensities unmapped; the sweep
  // reads an on-rank owner (vgd_own) directly, a remote one from the band, and applies
  // the wall map at the read, as one MeshBlock reads itself
  bool vgd_raw = false;
  DvceArray4D<int> vgd_own;   // (m, band k, band j, 3): owner's local index (-1 remote),
                              // its band k, j
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
  bool vlat_src_ep = false;
  Real vgd_nrb = 0.0, vgd_trb = 0.0;
  void VetGdIterRebuild();
  // the cs half-range data at the start of a solve (|H|/J, seam ghosts)
  void VetGdCsHrPrep();
  void VetGdSpHrHalo();   // sp ghosts of vgd_hr
  MeshBoundaryValuesCC *pbval_hr = nullptr;
  DvceArray5D<Real> vgd_hr_c;
  void VetGdRealDiag();
  void VetGdHaloInit();
  // (CS_VETGD.md sect. 4A, stage S3; rad_vet_vetgd_cs.cpp): vet_gd on the
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
  // message buffers.  in PINNED HOST memory (the pack/unpack kernels write/
  // read them directly), so MPI never sees a GPU pointer here: no GPUDirect RDMA and
  // no UCX cuda staging for the vet_gd band (pure copies: bitwise the device buffers)
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_cssb, vgd_csrb;
  int vgd_csnloc = 0, vgd_csnsnd = 0, vgd_csnrcv = 0;
  void VetGdCsInit();
  void VetGdCsSource();
  void VetGdCsSweep();
  void VetGdCsMoments(const int part);   // 0: D_rr -> LAT0, 1: surface q
  void VetGdCsTab(const std::vector<double> &d);   // into vgd_tabh (walls off, mr)
  template <class V>
  void VetGdCsHalo(V &a, const int nv, const int i0, const int i1, const bool loc,
                   const int ws = -1);
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_csown;   // (m, k, j, 3) owner
  // (USER decision 2): vet_gd_halo_compact > 0 on cs -- the chart-aware
  // exact lists: per shell and pass only the (remote band cell, direction) values the
  // sweep reads (VetGdCsCompact: a read-mask kernel that replays the sweep's upwind
  // stencils for the direction set, the receivers' masks sent to the owners once per set,
  // both sides' lists in the same order).  Bitwise the dense band copies (the same values
  // are read; every remote read is checked against the mask in the sweep: FATAL on a
  // miss)
  struct VgdCsCl {
    Real alpha = -9.0;
    int ng = 1;                      // direction groups of the lists
    int64_t use = 0;
    std::vector<int64_t> so, ro;     // ((pass*n1 + s)*(np + 1) + p): list offsets
    Kokkos::View<int *, LayoutWrapper, DevMemSpace> si, ri;   // e*n + d per value
    // (e_rcv, d, s): 1 in, 2 out, 0 not read; 2 bits per value, 16 values
    // per word (value t = (e_rcv*n + d)*n1 + s in word t/16, bits 2 (t%16))
    Kokkos::View<uint32_t *, LayoutWrapper, DevMemSpace> msk;
  };
  bool vgd_cscomp = false;
  // vet_gd_cs_list_pad; the read-mask stride per receive entry (n*nx1, or
  // rounded up to 16 with the pad); sub-timers of the list build (rank 0, fenced):
  // mask kernel, mask download + check, mark exchange, host lists, upload
  int64_t vgd_csper = 0;
  double vgd_cstb_k = 0.0, vgd_cstb_d = 0.0, vgd_cstb_m = 0.0, vgd_cstb_h = 0.0,
         vgd_cstb_u = 0.0;
  // persistent host buffers of the read-mask exchange (VetGdCsCompact):
  // own marks (sent) and the receivers' marks (received); grown, never freed
  std::vector<uint32_t> vgd_csmw_h, vgd_csmw_r;
  // the cache holds max(2, vet_rotate_avg) sets (vgd_csncl): with K-set
  // rotation (keys b mod K) every set's lists are made once, not at every rotation
  static constexpr int M1_VGD_CSCL_MAX = 16;
  VgdCsCl vgd_cscl[M1_VGD_CSCL_MAX];
  // hgmem: the read masks of two entries (the set in use and the one made or re-masked
  // next); the others keep their lists alone (VetGdCsCompact)
  Kokkos::View<uint32_t *, LayoutWrapper, DevMemSpace> vgd_csmb[2];
  int vgd_csncl = 2;
  int vgd_cscur = 0;
  int64_t vgd_csuse = 0;
  Kokkos::View<int ***, LayoutWrapper, DevMemSpace> vgd_csrix;   // (m, k, j): e_rcv / -1
  Kokkos::View<int *, LayoutWrapper, DevMemSpace> vgd_csmiss;   // reads not in the mask
  double vgd_csnbuild = 0.0, vgd_cstbuild = 0.0, vgd_csvsent = 0.0, vgd_csvdense = 0.0;
  void VetGdCsCompact();
  void VetGdCsHaloC(const int pass, const int i);
  bool VetGdCsHaloCX(const int pass, const int i, const int eb);
  // vet_gd_cs_xgroups = G > 1: the sweep's directions in G contiguous groups,
  // each (shell, group) its own kernel and compact exchange, software-pipelined (the
  // exchange of one group runs while the next groups compute; every value computed by
  // the same code from the same inputs: bitwise G = 1).  Lists ordered (pass, shell,
  // group, partner, entry, d); pinned-host slot buffers
  int vgd_csxg = 1;
  int64_t vgd_csxmax = 0;              // largest (pass, shell, group) list, values
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_csxsb, vgd_csxrb;
  DvceArray1D<Real> vgd_csdsb, vgd_csdrb;   // NCCL device buffers
  DvceArray1D<int> vgd_cslst;               // per-pass branch lists
  // vet_gd_halo_compact = 0 (off), 1 (pass branch), 2 (+ ray reach): compact messages
  int vgd_hcomp = 0;
  bool vgd_capped = false, vgd_allow_clamp = false;
  Real vgd_nclamp_seen = 0.0;
  Real vgd_npath = 0.0;   // vet_path subdivided steps (ALL ranks, all sweeps)
  bool vgd_hinw = true;
  int vgd_nb2 = 1, vgd_nb3 = 1;
  // several MeshBlocks along x1: the levels swept one after the other
  int vgd_nb1 = 1;
  std::vector<int> vgd_x1nb, vgd_x1lv;   // (2 m + 0/1): gid below / above, -1 none; lx1
  DvceArray1D<int> vgd_x1f;   // (m): 1 innermost, 2 outermost level, 4 reflecting top,
                              // + 8 lx1
  DvceArray1D<int> vgd_x1mo;  // the blocks by level (the shell list order)
  std::vector<int> vgd_shlv[2];   // (level): first code of the level in vgd_shl[pass]
  DvceArray2D<Real> vgd_x1b, vgd_x1r;    // (m, v): one shell to / from the x1 neighbour
  void VetGdX1Hand(const bool up);
  void VetGdX1Pass(const bool up);
  DvceArray4D<int> vgd_rmap;   // (m, k, j, d): mirror direction at a reflecting top
  void VetGdTopMirror();
  DvceArray3D<Real> vgd_isc;   // (m, k, j): flux ratio of a reflecting inner x1
  void VetGdInnerMirror();
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_lxy;   // (m, 0/1): lx2, lx3
  DvceArray1D<Real> vgd_csb, vgd_crb;
  // the compact mask of the NEXT shell is built (bitwise)
  // while this shell's messages travel; 2 slots of flags/positions (vgd_hfs/hfr) and of
  // the partner boundaries in pinned host memory (vgd_hpb), tagged (i, inw, ws, icut,
  // sweep id)
  DvceArray1D<int> vgd_hfs[2], vgd_hfr[2];
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hpb[2];
  int vgd_htag[2][5] = {{-1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1} };
  int vgd_hlast = 1, vgd_hsweep = 0;
  // vet_gd_halo_cache_mb (default 0 = off): byte budget (MB per rank) of a
  int vgd_hnext[3] = {-1, 0, 0};
  void VetGdHcPrep(const int slot, const int i0, const bool inw0, const int ws);
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_pbd;   // (ws, 2 (np+1))
  std::vector<double> vgd_r1v, vgd_r1f;
  template <class V>
  bool VetGdHaloCompact(V &a, const int nv, const int i0, const int ws);
  // vet_gd_halo_list_mb (read only when named; default 4096, 0 = off):
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
    std::vector<int> dsto;        // (np): start of my values in each
                                  // partner's receive layout (batch-made exact lists)
    Real alpha = -2.0;
    int scut = -1, ws = -1, ns = 0, nr = 0, nw = 0;
    int64_t gen = -1;             // build id (> 0) while valid
    bool nofit = false;           // over the budget for (alpha, scut, ws)
    size_t bytes = 0;
    // vet_gd_halo_exact: the lists hold only the values the sweep READS
    // (sr / sw = -1: a read ghost whose source is not in the mask: written 0)
    bool exact = false;
  };
  static constexpr int vgd_hl_mb = 4096;   // MB per rank: budget of the halo lists
  // the shell kernel runs over the (m, k, j, d) of the pass's branch only (VetGdSweep)
  DvceArray1D<int> vgd_shl[2];
  int vgd_shn[2] = {0, 0};
  Real vgd_shl_alpha = -3.0;
  DvceArray5D<Real> vcol_hr;     // (m, 3, k, j, i): vet_col h+_r, h-_r, f_K (unclamped)
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
  // (all bitwise unless stated)
  DvceArray4D<Real> vgd_trig;   // (m, k, j, 4): sin, cos of theta, phi (shell kernel)
  // (rad_vet_vetgd_sph2.inc): vet_gd on a STRETCHED theta and/or the FULL
  // sphere (polar x2 faces).  vgd_sth: the (theta) positions of the band rows come from
  // the node table vgd_thn (m, band row) = the face midpoint of the row in the EXTENDED
  // chart (beyond a pole theta continues to -theta / 2 pi - theta, i.e. the point
  // (theta, phi + pi)); the upwind fractional row is linear between the two nodes that
  // bracket it.  false = the uniform index formula (bitwise as before).  vgd_pol: the
  // full sphere: the halo slots beyond a pole take the pole partner (lx2, lx3 + nb3/2)
  // with the rows REVERSED (vgd_hpol(8 m + o) = 1 north, 2 south; on-rank partners are
  // read in-kernel with the reversed row), an upwind point behind the axis is taken in
  // the extended chart, and reads beyond the shell's stored band are clamped (counted)
  bool vgd_sth = false, vgd_pol = false;
  static constexpr Real vgd_bcolat = 30.0;   // (deg) the full sphere's band depth: sin
  DvceArray2D<Real> vgd_thn;              // (m, c2w) band-row theta (extended chart)
  std::vector<double> vgd_thn_h;          // host copy, m-major
  DvceArray1D<int> vgd_hpol;              // (8 nmb) pole slot flag
  DvceArray2D<int> vgd_wsd;               // (2, c1) per-shell halo depth: inward, outward
  // <rad_vet>/vet_gd_band_global (read only when named; default false =
  // bitwise; full sphere only, rad_vet_vetgd_sph2.inc): a decomposition-INDEPENDENT band.
  // (1) every band position of every block is filled by an exact copy of its owner's
  // interior value keyed on the GLOBAL cell (the cs band lists vgd_cs* and VetGdCsHalo,
  // shared with the cubed sphere; any depth, several blocks away: the band is no longer
  // capped at the MeshBlock size), per shell only the positions within the shell's halo
  // depth (the lists sorted by depth, prefix counts per depth); (2) the upwind stencil is
  // clamped to the CELL's own window (k +- wu, j +- wu, wu the read shell's halo depth)
  // instead of to the block's band edge, so a clamped polar read is the same wherever the
  // MeshBlock boundaries lie
  bool vgd_gband = false;
  DvceArray2D<Real> vgd_phc;              // (m, band k): band-column phi, global index
  std::vector<double> vgd_phc_h;          // host copy, m-major
  DvceArray1D<int> vgd_gdps, vgd_gdpr;    // depth of every send / receive entry
  std::vector<int> vgd_gdnl;              // (ws): on-rank entries with depth <= ws
  std::vector<std::vector<int>> vgd_gdsn, vgd_gdrn;   // (ws)(partner): send / recv counts
  void VetGdGlobalBandInit();
  // G1: vet_gd_halo_compact > 0 with vet_gd_band_global (the full sphere,
  // stretched theta): per shell and pass only the (band cell, direction) values that pass
  // VgdHaloKeep, evaluated at the SOURCE cell's global (row J, column K) on both sides
  // (node-table theta, global phi, the row's own theta spacing); rows within the polar
  // zone (node colatitude < vet_gd_band_colat) + the band depth keep their whole branch
  bool vgd_gbc = false;
  DvceArray2D<int> vgd_gbjs, vgd_gbjr;    // (entry, 2): source cell global J, K
  DvceArray1D<Real> vgd_gbth, vgd_gbdt;   // (J): node theta, theta cell width
  DvceArray1D<int> vgd_gbcap;             // (J): 1 = polar zone (branch only)
  DvceArray1D<int> vgd_gbfs, vgd_gbfr;    // flags -> positions (entry*n + v; +1)
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_gbpb;   // segment bounds
  // vgd_raw: the scanned flags and bounds per (shell, pass, depth) (VetGdGbHaloC)
  std::map<int, std::array<DvceArray1D<int>, 2>> vgd_gbcf, vgd_gbcf_xb;   // (_xb: the
  std::map<int, std::vector<int>> vgd_gbcb, vgd_gbcb_xb;   // other cross-fade set's)
  Real vgd_gbdmin = 0.0;
  template <class V>
  bool VetGdGbHaloC(V &a, const int i, const int ws);
  void VetGdThetaInit();                  // the node table and the checks
  double VetGdThetaDmin() const;          // the narrowest lateral spacing for the band
  int VgdPoleSlot(const int o, const int lx2) const;  // 0, 1 north, 2 south
  Kokkos::View<int64_t*, DevMemSpace> vgd_scan_sums;   // VgdScan2 chunk sums
  // vet_gd_halo_exact (read only when named; default false; bitwise): the compact-halo
  // lists of each (pass, shell) cut down to the band values the sweep actually reads
  // (marked from the sweep's own ray geometry, VetGdSweep), the senders told by one
  // request message per list build; no zeros of the other pass's entries needed.
  bool vgd_hx_on = false;
  DvceArray1D<int> vgd_hx_mk;              // (m, v, k, j) dense band: read marks
  int vgd_hx_stamp = 0;                    // the marks of the current build
  Kokkos::View<int, DevMemSpace> vgd_hx_oc;   // marked ghosts outside the region
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hx_pin, vgd_hx_pin2;
  DvceArray1D<int> vgd_hx_need, vgd_hx_nc, vgd_hx_ns;   // build scratch
  DvceArray1D<int> vgd_hx_t0, vgd_hx_t1;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hx_hb;   // bounds, counts
  Real vgd_hx_nmade = 0.0, vgd_hx_nout = 0.0, vgd_hx_nfall = 0.0, vgd_hx_nun = 0.0;
  Real vgd_hx_vmask = 0.0, vgd_hx_vexact = 0.0;   // received values: mask lists, exact
  std::vector<VgdHlEntry> vgd_hxe;         // (2 n1): the exact lists per (pass, shell)
  std::vector<VgdHlEntry> vgd_hxe_alt;     // the other xfade set's lists
  bool VetGdHxKey(const int e, const int ws) const;   // record e matches the masks' key
  bool VetGdHxCut(const int e, const VgdHlEntry &X) const;   // ray-cut move test
  void VetGdHxBuild2(const int e, const int i, const int ws);
  // the exact lists of many shells in a few collectives (VetGdHxBatch);
  // per-slot scratch pools, the per-depth partner / piece tables, pinned bounds
  DvceArray1D<char> vgd_hxb_tmp;
  int vgd_hxb_g = 0, vgd_hxb_npc = 0;
  size_t vgd_hxb_rcap = 0, vgd_hxb_scap = 0;
  DvceArray1D<int> vgd_hxb_nd, vgd_hxb_f0, vgd_hxb_f1, vgd_hxb_rq, vgd_hxb_so, vgd_hxb_oc;
  Kokkos::View<int**, LayoutWrapper, DevMemSpace> vgd_hxb_pst, vgd_hxb_pcs, vgd_hxb_pci,
                                                  vgd_hxb_pds;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hxb_hb, vgd_hxb_sp;
  Real vgd_hxb_t = 0.0, vgd_hxb_ngroup = 0.0, vgd_hxb_nrealloc = 0.0;
  void VetGdHxBatch(const std::vector<std::array<int, 4>> &jobs);
  // vet_gd_halo_ipc (read only when named; default false; bitwise): the
  // per-shell exchange of the batch-made exact lists over node-local CUDA IPC, driven by
  // the GPU (spin / flag kernels on the sweep's stream; the host never waits in the
  // sweep).
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
  void VetGdIpcInit(const std::vector<int> &prk, size_t cs, size_t half);
  void VetGdIpcCheck();
  // the pending list-path exchange of the overlap (VetGdHlBegin / VetGdHlEnd)
  struct VgdHlPend {
    bool on = false, xe = false;                // xe: X is an exact (vgd_hxe) entry
    bool ipc = false;                           // the IPC path
    bool nccl = false;                          // NCCL (stream-ordered) path
    // (CUDA IPC, lean): VetGdHalo's begin + end as ONE kernel (pack, flags, wait,
    // scatter).  now: the end follows
    // the begin (VetGdHalo); done: the scatter ran; cnt: the kernel's teams done packing
    bool now = false, done = false;
    Kokkos::View<int, DevMemSpace> cnt;
    long long seq = 0;
    int ex = -1, px = 0, i = -1;
    std::array<std::vector<int64_t>, 2> *gwa = nullptr;
    VgdRag a;
#if MPI_PARALLEL_ENABLED
    std::vector<MPI_Request> req;
#endif
    Kokkos::Timer tq;
  };
  VgdHlPend vgd_hlp;
  // (compile -DM1_VGD_NCCL; not used, VgdNcclOn = false): the per-shell list exchange of
  // the MPI path as stream-ordered NCCL send/recv on the sweep's stream
  int acc_nccl = -1;            // -1 not decided, 0 off, 1 on
  void *acc_ncomm = nullptr;    // ncclComm_t
  bool VgdNcclOn();
  template <class V>
  bool VetGdHlBegin(V &a, const int nv, const int i0, const int ws);
  void VetGdHlEnd();
  std::vector<int> vgd_wsh;    // (i): band depth the shell's data needs (<= vgd_w)
  std::vector<int> vgd_wsi, vgd_wso;   // (i): the same per pass (inward, outward)
  template <class V>
  void VetGdHalo(V &a, const int nv, const int i0, const int i1, const int ws,
                 const bool mapd);
  // the multi-shell (ln chi, ln S) halo in chunks of shells that fit the
  // per-shell message buffers (the same values; the buffers sized for one shell of all
  // directions instead of all shells of the 2 source fields)
  void VetGdHaloCs(DvceArray5D<Real> &a, const int i0, const int i1);
  Real vgd_tsrc = 0.0, vgd_tswp = 0.0, vgd_tmom = 0.0;
  void VetGdInit();
  Real VetGdAngle(const int cyc) const;
  bool VgdRotOn() const {return vgd_rot > 0 || vgd_rot_dt > 0.0 || vgd_rot_auto;}
  int64_t VetGdBlock(const int cyc) const;       // the rotation block
  Real VetGdAngleBlk(const int64_t b) const;
  double RotPhase(const int64_t b) const;   // block -> fraction of period
  Real VetGdAlignDt(const Real time, const Real dt);
  void VetGdTables(const Real alpha);
  // the tables' host part in persistent pinned buffers (no per-rotation
  // mirrors), and the next rotation's tables prefetched (host thread + upload into a
  // second device set on vgd_tex, swapped in at the rotation)
  struct VgdTabHost {
    Kokkos::View<Real **, LayoutWrapper, Kokkos::SharedHostPinnedSpace> dir;
    Kokkos::View<int ***, LayoutWrapper, Kokkos::SharedHostPinnedSpace> w;
    Kokkos::View<int ****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> mp;
    Kokkos::View<Real ****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> mr;
    Kokkos::View<int *****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> m3;
    Kokkos::View<Real *****, LayoutWrapper, Kokkos::SharedHostPinnedSpace> w3;
    std::vector<double> mt;
    Real alpha = -9.0;
    int nwall = 0, nexact = 0, nmap = 0, nphx = 0;   // nphx: inexact phi-only wall maps
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
  int vgd_tab_state = 0;       // 0 idle, 1 host thread running, 3 uploaded (vgd_tabd)
  Real vgd_tab_nswap = 0.0;
  DevExeSpace vgd_tex;
  void VgdTabCompute(const Real alpha, const int nthr);
  void VgdTabAlloc(const bool dev);
  void VetGdTabPrefetch(const int cyc);
  void VetGdTabPoll();
  void VetGdBuild();
  void VetGdSweep();
  void VetGdWall(const int i0, const int i1);
  void VetGdMoments();
  // the shell body closure in device memory, called through a pointer
  // (VgdClosureToDevice)
  DvceArray1D<char> vgd_bbuf0, vgd_bbuf1;
#if defined(KOKKOS_ENABLE_CUDA)
  DvceArray1D<char> vgd_bbuf2, vgd_bbuf3;   // CUDA gx-false closures
#endif
  int vgd_scut = 0;            // vlat_icut the running sweep uses (sync: = vlat_icut)
  DevExeSpace vgd_cur;         // the instance the sweep/halo use
#if MPI_PARALLEL_ENABLED
  MPI_Comm vgd_comm = MPI_COMM_WORLD;
#endif
  void VetGdSource();
  void VetGdPost();
  void VetLatCsFinish();

  // ...in "vet_before_stagen"
  TaskStatus InitRecv(Driver *d, int stage);
  // ...in "vet_stagen"
  TaskStatus CopyCons(Driver *d, int stage);
  TaskStatus ApplyClosureLimits(Driver *d, int stage);
  TaskStatus Opacity(Driver *d, int stage);
  TaskStatus CalculateFluxes(Driver *d, int stage);
  TaskStatus SendFlux(Driver *d, int stage);
  TaskStatus RecvFlux(Driver *d, int stage);
  TaskStatus HydroConToPrim(Driver *d, int stage);
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
