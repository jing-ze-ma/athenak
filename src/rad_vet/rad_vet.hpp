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
//! \brief the clean work vectors cln_w (CLW_*): slots n < ns in s, the layout's own
//! slots in o (n - ns; staggered: in mem_arena, MemArenaSeal)

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
constexpr int CLG_F = 27;    // staggered, exchanged: normal of the lower x2 face, volume
                             // (+3), normal of the lower x3 face (+4), MeshBlock gid (+7)
constexpr int CLG_A = 35;    // staggered, exchanged: mean area of the two faces of family
                             // x1, x2, x3 (+2), CleanStagArea
constexpr int CLG_Q = 38;    // staggered: face family and orientation in the owner's axes
                             // (CleanCsSlotQ), 4 bits per x2 (low) and x3 face
constexpr int CLG_NV = 39;
constexpr int CLG_NA = 25;   // components of CleanGeo::a

//! \struct CleanGeo
//! \brief cln_gn (CLG_*) in two arrays: components constant along x1 (CLG_N, CLG_R,
//! CLG_E, CLG_F, CLG_F + 4..6) once per (m,k,j) and x1 range (low ghosts, active, high
//! ghosts) in a; the others (CLG_X, CLG_M, CLG_F + 3, CLG_F + 7, CLG_A) in r

struct CleanGeo {
  DvceArray5D<Real> r, a;
  int is = 0, ie = 0;
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr int Slot(const int n) {   // a: n >= 0, r: -1 - n
    // a table, not a branch chain: keeps the face kernels' registers on the GPU
    constexpr signed char t[CLG_NV] = {0, 1, 2, 3, 4, 5, 6, 7, 8, -1, -2, -3, 9, 10, 11,
      12, 13, 14, 15, 16, 17, -4, -5, -6, -7, -8, -9, 18, 19, 20, -10, 21, 22, 23, -11,
      -12, -13, -14, CLG_NA - 1};
    return t[n];
  }
  // component B + o with B a compile-time start of a contiguous run: no table lookup
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

// explicit stages of the PD-ARS (Chu et al. 2019) integrator per substep
constexpr int M1_NSTAGE = 2;

//! \fn M1KeyName
//! \brief the <rad_vet> name to read for a key with aliases: canon if present, else the
//! first alias present, else canon. The canonical key wins over a differing alias.
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

//! \fn M1LzKey
//! \brief key of the clean L_z ledger (and K1's sum) on every mesh; old files' cs / sp
//! spellings are read as aliases
inline std::string M1LzKey(ParameterInput *pin, const bool cs, const bool k1) {
  return M1KeyName(pin, k1 ? "clean_lz_ledger_k1" : "clean_lz_ledger",
                   {cs ? (k1 ? "clean_cs_lz_k1" : "clean_cs_lz_ledger")
                       : (k1 ? "clean_sp_lz_bnd_k1" : "clean_sp_lz_bnd")});
}

// components of RadiationVet::vet_cell (closure = vet_sc)
constexpr int M1_VET_CHX = 0;   // extinction rho (kappa_F + kappa_s) [1/length]
constexpr int M1_VET_SRC = 1;   // source function in E units (4 pi S / c)
constexpr int M1_VET_J   = 2;   // E of the formal solution (4 pi J / c)
constexpr int M1_VET_K11 = 3;   // K_ab in E units: 11 22 33 12 13 23
constexpr int M1_VET_H1  = 9;   // F_a / c of the formal solution: 1 2 3
constexpr int M1_VET_CHI = 12;  // the uniaxial projection: chi, n1, n2, n3
constexpr int M1_VET_N1  = 13;
constexpr int M1_VET_D11 = 16;  // vet_tensor = full: guarded D = K/J, 11 22 33 12 13 23
constexpr int M1_VET_GD  = 22;  // full: |D_guarded - D_raw| (max norm) of the cell
constexpr int M1_VET_HP1 = 23;  // sum w I max(mu_a, 0), a = 1 2 3 (half-range)
constexpr int M1_VET_NL  = 26;  // sum w |I - S| (C1 non-locality), fs only
constexpr int M1_VET_EPS = 27;  // thermal fraction min(kappa_P/chi, 1) (fs)
constexpr int M1_VET_NC  = 28;
// scheme = clean: K+^d_bd = sum w I Omega_b mu_d over mu_d > 0 (slot KP + 3 d + b)
constexpr int M1_VET_KP  = 28;
constexpr int M1_VET_NKP = 9;
// clean_consistency (VetShortCharSweep): ray fluxes sum w g Phi, sum w g Phi Omega_b
// (g = |mu_1| dx2/dx1) through the low x2 face; sum w mu_1 I, sum w mu_1 I Omega_b
// through the x1 boundary face and the low x1 face; low x3 face as x2
constexpr int M1_VET_F2  = 37;
constexpr int M1_VET_F1  = 41;
constexpr int M1_VET_FX  = 45;
constexpr int M1_VET_F3  = 49;
constexpr int M1_VET_NFS = 16;  // the face sums F2, F1, FX, F3
// tau_ten slots of the lateral correction (vet_col_lat): dD_ab in the mesh basis,
// rr, r-a, r-b, aa, ab, bb (sp: a = theta, b = phi)
constexpr int M1_TT_LAT0 = 4;
constexpr int M1_TT_NLAT = 6;
// extra vgd_hr components of the vet_gd thin switch (VetGdHalfRange)
constexpr int M1_VGD_HRC1 = 6;     // C1 = sum w |I - S| / J
constexpr int M1_VGD_HREP = 7;     // eps = min(kappa_P/chi, 1) (1 where chi = 0)
constexpr int M1_VGD_HRJ  = 8;     // J (E units)
constexpr int M1_VGD_HRNX = 9;
#ifdef M1_VGD_NMAX_BUILD
constexpr int M1_VGD_NMAX = M1_VGD_NMAX_BUILD;   // cmake -D Athena_M1_VGD_NMAX
#else
constexpr int M1_VGD_NMAX = 432;   // vet_gd: max directions (sizes buf[] in VetGdWall)
#endif
// components of opac, all stored as rho*kappa (inverse length)
constexpr int M1_OP_P = 0;  // rho kappa_P, Planck (emission) mean
constexpr int M1_OP_E = 1;   // rho kappa_E, energy (absorption) mean
constexpr int M1_OP_T = 2;   // rho (kappa_F + kappa_s): transport opacity (flux source,
                             // tau_face)
constexpr int M1_NOPAC = 3;
// vet_scatter only: rho kappa_e, the electron-scattering part of M1_OP_T (which stays
// the total rho kappa_R)
constexpr int M1_OP_S = 3;

// <rad_vet>/opacity
constexpr int M1_OPAC_CONST    = 0;
constexpr int M1_OPAC_POWERLAW = 1;
constexpr int M1_OPAC_USER     = 2;
constexpr int M1_OPAC_TABLE    = 3;   // stellar Rosseland + Planck tables

// <rad_vet>/force_reference: what the momentum coupling subtracts from the radiative
// force (arad_ref; see the head of rad_vet_fluid.cpp)
constexpr int M1_FREF_NONE    = 0;
constexpr int M1_FREF_WB_ARAD = 1;
constexpr int M1_FREF_CLEAN   = 2;   // scheme = clean: its own force of the first solve
constexpr int M1_FREF_WB_MD   = 3;   // wb_arad_md: staggered clean, steady-start force in
                                     // every direction (cln_fref)


//----------------------------------------------------------------------------------------
//! \struct M1OpacTab
//! \brief the stellar opacity tables as device Views held by value (capturable in a
//! device lambda); one (log10 T, log10 rho) grid, read with RosselandTable.
//! Lookup in K and g/cm^3: T[K] = T_code*tunit, rho = rho_code*dunit; kunit converts
//! cm^2/g to code units.

struct M1OpacTab {
  DvceArray2D<Real> kr;        // log10 kappa_Rosseland [cm^2/g], (iT, iD)
  DvceArray2D<Real> kp;        // log10 kappa_Planck    [cm^2/g], same grid
  DvceArray1D<Real> lT, lD;    // log10 T[K], log10 rho[g/cm^3]
  int nT = 0, nD = 0;
  Real tunit = 1.0;
  Real dunit = 1.0;
  Real kunit = 1.0;
  // opac_abs_rho_max: below amask_rho no absorption and transport opacity amask_kes
  Real amask_rho = 0.0;
  Real amask_kes = 0.0;
  // vet_scatter: split of kappa_R into kappa_a = kappa_R - kappa_e and kappa_e.
  // scat 0 none, 1 constant kes, 2 kxe x_e/mu (M1ScatterEos); code units (times kunit)
  int scat = 0;
  Real kes = 0.0;    // constant kappa_e (vet_kappa_es)
  Real kfl = 0.0;    // floor: kappa_e <= kappa_R - kfl (vet_kappa_floor)
  Real kxe = 0.0;    // sigma_T/m_u: kappa_e = kxe n_e/(n_tot mu)
};

//----------------------------------------------------------------------------------------
//! \fn VetScatterMix
//! \brief vet_scatter: grey source S (E units) of the VET formal solutions for total
//! extinction chi = rho kappa_R, rkp = rho kappa_P, rkes = rho kappa_e, thermal source b
//! and mean intensity e.
//!   form 1 (ma, Jiang 2021 eq. 6): S = max(J + (rho k_P/chi)(B - J), 0)
//!   form 2 (absorption): S = [(chi - rkes) B + rkes J]/chi
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
//! \brief guarded kappa_e (per mass) for total kt: kappa_a = kt - kappa_e >= kfl
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

// bracketed-Newton controls
constexpr int  M1_MAXIT = 100;
constexpr Real M1_RTOL  = 1.0e-12;

//----------------------------------------------------------------------------------------
//! \struct RadiationVetTaskIDs
//! \brief container to hold TaskIDs of all rad_vet tasks

//! \struct CleanRingMates
//! \brief the MeshBlocks of each local block's line along a direction (CleanRingMatesInit)
struct CleanRingMates {
  int nmb = -1, nt = -1, nslot = 0;   // blocks here and in the mesh; received blocks
  DvceArray2D<int> src;               // (m, position): local block or -1 - received slot
  std::vector<std::array<int, 3>> snd, rcv;   // (block or slot, rank, source gid)
};

struct VetMBState;   // multi-block short-characteristics sweep (rad_vet_vet.cpp)

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
//! \brief the grey M1 / VET radiation module

class RadiationVet {
 public:
  RadiationVet(MeshBlockPack *ppack, ParameterInput *pin);
  ~RadiationVet();

  // <time>/dt_retry: dtr_on: a failed cycle is retaken with a smaller dt; a solve that
  // would stop the run raises dtr_fail with reason dtr_why instead
  bool dtr_on = false;
  bool dtr_fail = false;
  // 1 NON-CONVERGED, 2 DIVERGED, 3 DIVERGED (masked), 4 eint<=0, 5 time2_reject
  int dtr_why = 0;
  Real dtr_resid = 0.0;
  // dt_retry_t2_be cuts did not cure a time2_reject: redo the hesdirk2 stage with
  // backward Euler instead of another retry
  bool dtr_t2be = false;
  bool dtr_t2hit = false;    // this attempt took that backward-Euler redo
  template <class V> void DtRetryVisit(V &v);   // rad_vet_dtretry_visit.hpp (generated)
  void DtRetryQuiesce(bool restoring);          // no host thread writes while restoring
  void DtRetryAfterRestore();                   // make the async prefetch state safe

  // parameters read from <rad_vet>
  Real c_light{};      // speed of light in code units (required)
  Real chat_over_c{};  // reduced speed of light factor (default 1)
  Real chat{};         // = chat_over_c * c_light
  Real e_floor{};      // floor on E
  int chi_kind{};      // M1_CHI_*: levermore | minerbo | kershaw
  bool eddington{};    // closure = eddington (chi = 1/3); default is the M1 closure
  bool advect_split{};          // the advective enthalpy-flux split of the E equation

  // opacities per unit mass (code units): kappa0, or kappa0 (rho/rho_ref)^a (T/t_ref)^b
  int opacity_type{};           // M1_OPAC_*
  Real kappa_p{}, kappa_e{}, kappa_f{}, kappa_s{};
  Real opac_rho_ref{}, opac_t_ref{}, opac_a{}, opac_b{};
  bool opac_zero{};             // all four opacities are identically zero
  M1OpacTab otab;               // opacity = table: filled by SetOpacityTables

  // <rad_vet>/force_reference: with M1_FREF_WB_ARAD the gas gets the residual between the
  // radiative force and the reference acceleration arad_ref of an external well-balanced
  // scheme (filled by the pgen through SetForceReference)
  int force_ref{};
  DvceArray4D<Real> arad_ref;   // (m,k,j,i), a reference x1 acceleration, code units
  // optional explicit radiation energy source (code units), set by a pgen with
  // esrc_on = true (implicit only): adds dt*(chat/c)*esrc to E's old vector and slope
  DvceArray4D<Real> esrc;
  bool esrc_on = false;
  // vet_source_noesrc: the vet_col / vet_col_lat / vet_gd source reads E^n without the
  // explicit esrc deposit (be only)
  bool vsrc_noes = false;
  // vet_source_noesrc_h2be: the same for the backward-Euler solves under hesdirk2
  bool vsrc_noes_h2be = true;
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
  bool fcap_mvol = false;   // floorcap_mesh_vol: counters by CleanVol and CleanLzArm
  DualArray1D<Real> fcap_cnt;   // removed: 0 n_v, 1 n_t, 2 dKE, 3 dEint, 4 |dp|V, 5 dL_zV
  void FloorGasCap();

  // RSLA check: v_max*tau_max/chat once, on the first opacity (rsla_vmax <= 0: measure)
  Real rsla_warn{}, rsla_vmax{};
  bool rsla_force{}, rsla_done{};

  // matter coupling
  bool coupling{};     // run the implicit local solve at all
  bool gas_feedback{}; // write the solve's reaction back to the gas (off: static medium)
  Real arad{};         // radiation constant in CODE units: the equilibrium energy
                       // density is arad*T^4 with T the EOS's code temperature

  // debug: dbg_gas_force = false gives the gas no momentum and no work
  bool dbg_gas_force{};

  // ---- vet_sc short characteristics (rad_vet_vet.cpp): Eddington tensor K/J of a grey
  // formal solution, once per hydro step from the start-of-step state
  bool vet_sc{};             // closure = vet_sc (default false)
  int vet_nmu = 0, vet_nphi = 0, vet_nray = 0;  // mu nodes per hemisphere, azimuths, rays
  bool vsc_gd = false;       // vet_sc_set = gd: the vet_gd direction set, pole on x1
  int vsc_gdmir = 0;         // ray r + nray/2 is the x1 mirror of ray r (1) or not (-1)
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
  bool vet_seg_split = false;  // vet_seg = split: SC segment u -> 0
  bool vet_seg_emit = false;   // vet_seg = emit: split, limiter only where S_u > S_0
  bool vet_seg_jump = false;   // vet_seg = jump: emit + the contrast rule
  bool vet_seg_plateau = false; // vet_seg_plateau: jump's contrast only across a plateau
  // VetSegSplit mode: 0 split, 1 emit, 2 jump with plateau weights, 3 jump without
  int VetSegMode() const {
    return vet_seg_emit ? (vet_seg_jump ? (vet_seg_plateau ? 2 : 3) : 1) : 0;
  }
  bool vet_x1per = false;    // vet_x1_periodic: periodic x1 sweep
  int vet_x1npass = 0;       // its number of passes
  // vet_x1_ptol: stop the periodic sweep when d_p/(1 - q) < vet_x1_ptol (0: npass)
  Real vet_x1ptol = 0.0;
  DvceArray4D<Real> vet_x1ent;   // (nmb, nray, k, j): the previous pass's exit plane
  // a build has run, so the next adaptive sweep starts warm from vet_x1ent; both are in
  // the restart file (kM1X1eRstMagic; a file without them restarts the sweep cold)
  bool vet_x1_warm = false;
  int X1eRstNslab(const int n1) const {
    return vet_x1_warm ? (vet_x1ent.extent_int(1) + n1 - 1)/n1 : 0;
  }
  void X1eRstPack(DvceArray5D<Real> &dst, const int nmb) const;
  bool X1eRstUnpack(const DvceArray5D<Real> &src, const int nmb, const int nsl);
  Real vet_x1_nbuild = 0.0, vet_x1_npsum = 0.0, vet_x1_ncap = 0.0;
  int vet_x1_npmax = 0;
  int vet_interp = 0;        // upwind stencil 0 lin, 1 log, 2 corner
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
  void VetBuildAngles(const Real beta);
  void VetGdAngles(const Real beta);
  void VgdBuildSet(std::vector<double> &d) const;
  void VetRotInit(ParameterInput *pin);   // shared rotation keys (rad_vet_vetgd.cpp)
  void VetRotUpdate();
  Real VetRotBeta(const int64_t blk) const;   // rotation of block blk
  void VetRotSet(const Real beta);
  void VetScRotAvgSetup();     // vet_sc's averaged fields (vet_cell)
  void VetShortCharSweep();    // build with the current set
  bool vet_axis_flux = false;  // uniaxial axis: the SC flux (true) or the principal axis
  bool vet_full = false;     // vet_tensor = full: the solve reads all six D_ab = K_ab/J
  static constexpr Real vet_eig_min = 0.0;   // full: eigenvalue floor of D
  Real vet_nguard{}, vet_ncell{}; // full: cell-calls the guard changed D / all cell-calls
  Real vet_guard_max{};      // full: largest |D_guarded - D_raw| (max norm)
  Real vet_d11_min{};        // full: smallest guarded D_11
  Real vet_time{}, vet_itime{}, vet_ncall{}; // SC seconds, implicit-solve seconds, calls
  DvceArray2D<Real> vet_ang;   // (nray, 4): mu1, mu2, mu3, weight (sum of weights = 1)
  DvceArray5D<Real> vet_cell;  // (nmb, M1_VET_NC, k, j, i), M1_VET_*
  DvceArray5D<Real> vet_ipl;   // (nmb, 2, nray, k, j): intensity of the last two layers
  // multi-block sweep state (plane buffers, neighbour tables, exchange); null on one
  // MeshBlock
  VetMBState *vet_mbs;

  // fold dtnew into Mesh::NewTimeStep (no sub-cycling, or no other module sets dt)
  bool sets_mesh_dt{};

  // implicit transport (rad_vet_implicit.hpp)
  int nstage{};                 // 1: one stage solve per task-list pass
  Real impl_cfl{};              // <rad_vet>/implicit_cfl: dtnew = impl_cfl*dx/c (<=0 off)
  DvceArray1D<Real> opn_nskip_d;
  // implicit_pos_gas: energy-conserving positivity limiter of the written-back gas eint;
  // implicit_pos_floor: the E floor takes its energy from the gas. pos_cnt_d: device
  // counters (M1_POS_*), pos_cnt the all-rank totals
  bool impl_pos_gas = false;
  bool impl_pos_floor = false;
  Real impl_pos_gas_frac = 1.0e-3;
  DvceArray1D<Real> pos_cnt_d;
  Real pos_cnt[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // E-floor energy created by discarded attempts (rank-local), not in pos_cnt
  Real pos_un_drop = 0.0;
  int dtr_ncut = 0, dtr_maxcut = 8;   // DtRetry cuts of this cycle, limit
  bool dtr_redo = false;   // this cycle's failed step is redone (clean: solved strictly)
  DvceArray4D<Real> ktd;        // d(rho kappa_T)/dT at the iterate (implicit_opac_newton)
  Real marshak_q{};             // free-surface condition F_f = c*marshak_q*E
  int ibc_x1min{}, ibc_x1max{}; // M1_IBC_*
  // implicit_bc_x2min|x2max|x3min|x3max: auto (periodic on a periodic mesh, else
  // reflect) | reflect | marshak (F_f = c marshak_q E outward; formal solution I = 0)
  int ibc_x2min = 2, ibc_x2max = 2, ibc_x3min = 2, ibc_x3max = 2;
  Real iflux_x1min{}, iflux_x1max{}; // imposed face flux of M1_IBC_FLUX
  Real iebath_x1min{}, iebath_x1max{}; // M1_IBC_MARSHAK: incident bath
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
  // implicit_marshak_face = linear (spherical-polar wedge): the Marshak face E from
  // r^2 E of the two end cells (M1SphMarshakCoef); cell (default): the end cell's E.
  // The outer face keeps cell under vet_col_surface_q.
  bool impl_mface_lin = false;
  DvceArray4D<Real> f0x1;       // face-normal comoving x1 flux (m,k,j,i), restart state
  DvceArray5D<Real> iw;         // per-cell work array of the solve, M1_NIW components
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
  bool fref_wacc_t2 = false;   // <rad_vet>/time2_fref_wacc (see FrefWaccOn)
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
  int impl_flux{};              // M1_IFLUX_*: the spatial form of the implicit E-flux
  int impl_recon{};             // M1_IRECON_*: dc, or plm by deferred correction
  int impl_enth{};              // M1_IENTH_*: face value of the enthalpy flux a E
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
  //  implicit_thin_gw_dump: g_L, g_R and weight AL of the last build (output m1_gw)
  DvceArray5D<Real> gwd;
  // time2_face_fct: the restart file carries K1's face energy flux t2pf1 (one slab per
  // direction, kM1FaxRstMagic). Old files' implicit_hr_cmem_rst slabs: read, ignored.
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
  // zalesak_halo_skip: Zalesak rounds exchange R only when a cell within 2 ng of a face
  // changed. Counters: rounds, exchanges done.
  bool zl_hskip = true;
  // zalesak_active: a Zalesak round evaluates only cells whose stencil changed in the
  // previous round or its R exchange (others cannot be newly limited: bitwise)
  bool zl_active = false;
  Real zl_rounds = 0.0, zl_halos = 0.0;
  // implicit_thin_switch (M1_THSW_*): dt-free thin/thick weight frozen per solve
  // (formal-solution |H|/J, C1, E^n Knudsen number); implicit_thin_corr adds the lagged
  // formal-solution flux
  int impl_thsw = 0;
  bool thsw_c1 = false;           // the sweep accumulates M1_VET_NL
  Real thsw_h0 = 0.3, thsw_h1 = 0.6;
  // implicit_thin_h0_x1 / _h1_x1: thin-switch band of the x1 faces (_hx1_curv: only on
  // curvilinear meshes); implicit_thin_recon = js (M1ThinJs: ratios from the face mean
  // of H+-); implicit_thin_bq = lin (outer Marshak q from the face-extrapolated H+);
  // implicit_thin_hcsf = jh: hc/hcs take the sf donor scale Jh (vet_sc), j: the cell J
  DvceArray5D<Real> ifw;        // per-x1-face work array, M1_NIFW components
  // implicit_flux_faces = all: berthon | blend on every face of the multi-D solve
  // (Cartesian and sp lateral; sp radial keeps its own). ifw2/ifw3: AL/HCL/HCR of the
  // x2/x3 faces; a transverse face's central part keeps 1 - AL through thx2/thx3.
  bool impl_beam_fs = false;    // implicit_flux_beam = fs: |f| = |H|/J of vet_sc
  // implicit_flux_beam = halfrange: F = (1 - w) F_central + w c (r+_L E_L + r-_R E_R),
  // r+- = H+-/J the half-range ratios along the normal, w = exp(-(tau_f/tau0)^2)
  // (implicit_blend = tau) or 1 (berthon)
  bool impl_beam_hr = false;
  // implicit_hr_face = sf (default with vet_sc): face ratios from the face mean of the
  // half-range fluxes over the donor scale Jh = max(J, outgoing SC face fluxes) where J
  // does not carry the SC flux (M1HrSf); cell: the cell ratios
  bool impl_hr_sf = false;
  // implicit_hr_excess = iso (vet_gd): the part of E^n above the formal-solution J takes
  // the isotropic half-range shape (+-1/4) (M1HrExcess)
  double hsf_ncell = 0.0, hsf_nsf = 0.0;   // cell-solves, of which Jh > J (rank 0)
  bool hr_q = false;            // per-column outer q from h+_1 (vcol_q)
  DvceArray5D<Real> vgd_hr, vgd_hrt;   // (m, 6, k, j, i): r+_r r-_r r+_t r-_t r+_p r-_p
  // with implicit_thin_switch = fs* or implicit_thin_corr vgd_hr has M1_VGD_HRNX
  // components: + C1, eps, J
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> vgd_pin_cnt;
  DvceArray5D<Real> ifw2, ifw3;
  // implicit_gas_src = matched (term-by-term deposit) | conservative: the gas takes minus
  // the radiation's own source, dm = -(F'_c - F^old_c)/(c chat) - (1/2) sum_f Z_f, with
  // F_c the stored cell flux and Z_f the face's transport part
  bool impl_gas_src_cons = false;
  DvceArray5D<Real> csk1;               // hesdirk2: slope of the cell F (3)
  bool csk1_ok = false;                 // csk1 holds the slope that pairs with t2k1
  bool gcs_first = true;                // the first conservative solve of this run
  // restart (kM1CskRstMagic): the header once a conservative solve ran, 3 slabs when
  // csk1 holds a slope
  bool CskRstHdr() const {return impl_gas_src_cons && !gcs_first;}
  int CskRstNch() const {
    return (CskRstHdr() && csk1_ok && csk1.extent_int(0) > 0) ? 3 : 0;
  }
  int cs2_vnp = 0;                      // implicit_cs2_vimp_npass
  // time2_enth_be: this (non-stage) solve takes the stage a_f rule
  bool T2EnthBe() const;

  // ---- implicit transverse (x2/x3) coupling
  int impl_solver{};            // M1_ISOLV_*: line_jacobi | bicgstab
  bool trans_on{};              // transport == implicit and the mesh is multi-D
  bool trans_x3{};              // ...and three-dimensional
  DvceArray4D<Real> f0x2, f0x3;   // face-normal comoving x2 / x3 fluxes (m,k,j,i),
                                  // restart state (trans_on only)

  // ---- implicit_solver = bicgstab (null under line_jacobi)
  bool bicg_on{};               // impl_solver == M1_ISOLV_BICGSTAB and trans_on
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
  int impl_pord{};              // implicit_predictor_order: 1 or 2
  int impl_prec{};              // implicit_precond: 0 x1 line Jacobi, 1 symmetric
                                // red-black transverse line Gauss-Seidel (block-local),
                                // 2 its forward half, 3 mg (rbgs_fwd + block-local
                                // semicoarsening), 4 mg_gc (mg + global coarse space)
  // implicit_precond = mg_gc: global coarse space of per-x1-layer (x2,x3)-band means
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gc_part;  // (m,lb,q,c,i) band sums
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gc_x;     // (gi*nb + b) x_g
  // implicit_precond = mg_gf: global coarse space of the lowest (x2,x3) Fourier modes
  Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> gf_gh;   // (gi*nm + q) multi-rank
  DvceArray1D<int> gf_okd;          // (q) the mode-q factorisation is usable
  Real bcg_r0rel{};             // max|b - A x0|/max|b| of the last BiCGStab call
  DvceArray1D<Real> lay_part;   // layout_invariant_mode = block: per-MeshBlock
                                // partials of M1LayReduce (sized on first use)
  Kokkos::View<Real*, Kokkos::HostSpace> lay_part_h;   // their host copy
  Real ew_fprev{}, ew_etaprev{}; // max|r0| and eta of the previous pass of this step
  bool impl_pred{};             // implicit_predictor = step: start the Picard loop from
                                // E^n + (dt/dt_prev) dE_prev (and T); none: cold start
  bool pred_ok{};               // a previous-step increment is stored
  Real pred_dt{};               // dt of that step
  DvceArray5D<Real> ipred;      // (m,3,k,j,i): dE_prev, dT_prev, T^n of this step

  // ---- time_scheme = hesdirk2: two implicit stage solves inside the Heun stages
  int time_scheme{};            // M1_TIME_BE | M1_TIME_HESDIRK2
  bool t2_ok{};                 // a valid FSAL slope K1 is stored
  bool t2_enth_be{};            // time2_enth_be: the backward-Euler solves of a hesdirk2
                                // run take the stage a_f rule
  int t2_solve{};               // what the next ImplicitSolve does (M1_T2S_*)
  bool t2_fail{};               // the last stage solve was not admissible
  int t2_dbg_fail{};            // debug time2_dbg_fail: fail stage 1 at this cycle
  // time2_reject = be | retry: a non-admissible stage raises the dt_retry flag instead
  // of the backward-Euler redo
  int t2_reject = 0;            // 0 be, 1 retry
  // time2_k1_stiff = s: move the exchange part of K1 in cells with k dt > s (0 = off)
  Real t2_k1stiff = 0.0;
  Real t2_k1stiff_n = 0.0;      // cell-steps changed (this rank)
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
  // time2_accept_fallback: a stage whose solve took an operator positivity fallback is
  // admissible when converged and positive (the backward-Euler redo uses the same
  // operator)
  bool t2_accfb = false;
  // time2_adm_poslim: the admissibility test counts a gas eint fault only when the
  // write-back's limiter (implicit_pos_gas) cannot cover it; E or T <= 0 stays a fault.
  // Needs implicit_pos_gas.
  bool t2_adm_pl = false;
  Real t2_adm_pl_n = 0.0;       // stages admitted by it (would have been rejected)
  // implicit_hold: cells a pgen holds (hold_m = 1) are Dirichlet rows of every implicit
  // solve (E' = held E, no increment, kick or work); the write-back leaves them to the
  // pgen
  bool impl_hold = false;
  DvceArray4D<Real> hold_m;
  Real t2_accfb_n = 0.0;        // stages accepted with a fallback (counter)
  Real t2_nstep{}, t2_nbe{}, t2_nfall{}; // stage steps, BE steps, fallbacks (counters)
  // time2_auto_be = N (0 off): after a stage fallback the next steps take the
  // backward-Euler redo path at once; stages are retried after 1, 2, 4, ... N skipped
  // steps (reset on success). The skip-run state is in the restart file.
  int t2_abe = 16;              // N
  int t2_abe_len = 0;           // current skip-run length (0 = stage mode)
  int t2_abe_left = 0;          // skipped steps left in the current run
  bool t2_abe_now = false;      // this cycle skips the stages
  // time_scheme = be runs as hesdirk2 with every step a time2_auto_be skip step
  // (Time2AutoBeBegin); also time2_be_only = true under hesdirk2
  bool t2_beonly = false;
  // time2_cfl_post: the next step's hydro CFL dt is taken again after the radiation
  // update, on the gas it wrote (a backward-Euler kick can exceed the pre-coupling CFL)
  bool t2_cfl_post = true;
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
  bool t2_vext{};               // time2_vet_extrap (false: D^n)
  int t2_vsmode = 0;            // time2_vet_sc: 0 lag, 1 predict, 2 rebuild
  // time2_tableau: 0 hesdirk2 (Heun + H-ESDIRK2), 1 trbdf2 (TR-BDF2, RK2 hydro with
  // c2 = 2 - sqrt 2). Time2FormStage: stage 1 t2inc = t2_a21 dt K1; stage 2 start =
  // t2_w0 Y2 + t2_w1 U^n, t2inc = dt (t2_c31 K1 + t2_c32 K2); diagonal g = 1 - 1/sqrt 2.
  int t2_tab = 0;
  Real t2_a21 = 0.0, t2_w0 = 0.5, t2_w1 = 0.5, t2_c31 = 0.0, t2_c32 = 0.0;
  Real t2_onep_s{};             // time2_one_pass_safety of the stage solves (0: the
                                // implicit_one_pass_safety)
  Real t2_lin_save{};           // implicit_lin_tol put back after a stage solve
  Real t2_nclip{};              // vet_sc cells whose extrapolated tensor was clipped
  DvceArray5D<Real> t2k1;       // (m,M1_T2_NK,k,j,i) the FSAL slope K1 (restart state)
  DvceArray5D<Real> t2k2;       // (m,M1_T2_NK,k,j,i) the stage-2 slope K2
  DvceArray5D<Real> t2inc;      // (m,M1_T2_NK,k,j,i) old vector - stage start state
  DvceArray4D<Real> t2f1, t2f2, t2f3;  // face fluxes of U^n (Heun average)
  bool pred2_ok{};
  Real pred2_dt{};
  DvceArray5D<Real> vet_prev;   // vet_sc: start-of-step tensor of the previous step
  DvceArray5D<Real> vet_now;    // vet_sc: start-of-step tensor of this step
  DvceArray5D<Real> vet_opac;   // vet_sc: stage-start opacities, saved over the formal
                                // solution at U^n
  void Time2Init(ParameterInput *pin);
  bool Time2Active();           // the next step runs the stages (else BE)
  void Time2FormStage(int stage, Real dt);
  // sdirk2 with Heun hydro: the stage-1 solve on the gas at t^n + g dt, the hydro's
  // stage-2 flux at t^{n+1} (t2hs: the gas U* of Heun stage 1)
  void Time2GasTime(Driver *pdrive, int which);
  DvceArray5D<Real> t2hs;
  int t2_gt = 0;   // last Time2GasTime step of this stage loop
  bool t2_gtm = true;   // time2_gtime_mom: the flux state's momentum factor 1/g
  bool t2_gts = false;  // time2_gtime_slope: its eint by the limited exchange slope
  void Time2Restore(Driver *pdrive);

  bool t2_inprobe = false;      // inside a helper solve (Time2RotCorrect)
  void Time2VetStart();        // vet_sc: formal solution at U^n, then extrapolate
  void Time2VetExtrapolate();
  // closure = vet_col under hesdirk2, time2_vet_col: lag (D^n) | predict (one build at
  // U^n + dt K1) | rebuild (predict for stage 1, a build at Y1 for stage 2) | extrap
  // (diagnostic: D^n + r (D^n - D^{n-1}), not in the restart)
  int t2_vcmode = 0;
  // cs_vertex_mono (cubed sphere): two-point lateral transport near the 8 cube vertices.
  // Lateral faces with a cell within cs_vertex_mono_n of a vertex take the two-point
  // gradient (no skew cross term), seam faces the foot values (no midpoint quadratic),
  // and D_r,lat = 0 within cs_vertex_mono_n + 1. Conservative per face, local.
  Real impl_pred_tcap = 2.0;   // implicit_predictor_tcap
  bool cs_vtx_mono = false;    // default on for cs geometry
  int cs_vtx_n = 1;
  // implicit_pos_floor_cons (needs implicit_pos_floor): the E floor creates no energy.
  // The raise is paid by the cell's gas above its floor, then by the radiation above
  // e_floor of the same radial column, nearest cells first (donor slopes and
  // |F| <= c E kept). The column is in one MeshBlock: local, exact, layout invariant.
  // M1_POS_FLR_MV counts moves, M1_POS_FLR_UN columns with nothing left.
  bool impl_pos_flr_cons = true;
  // implicit_pos_floor_gas_keep th (-1 off): the clean E floor raise takes at most
  // (1 - th) of the gas eint, the rest from the column's radiation
  Real impl_pos_flr_gk = -1.0;
  // implicit_pos_gas_cons (needs implicit_pos_floor_cons): the part of the
  // implicit_pos_gas raise the cell's E cannot cover comes from the radial column's
  // radiation instead of the hydro EOS floor
  bool impl_pos_gas_cons = false;
  DvceArray4D<Real> pgc_def;      // (m,k,j,i) uncovered gas raise x volume
  DvceArray4D<Real> pfc_def;      // (m,k,j,i) uncovered floor energy of the cell
  bool t2_vcprev = false;       // extrap: vcol_prev holds the previous step's tensor
  Real t2_vcnfb = 0.0;          // predict/rebuild builds that fell back to U^n
  // time2_vet_col_source = state: the predict/rebuild build takes its VET source from the
  // predicted state itself (vet_col_source = gas); relaxed (default)
  bool t2_vcsrc = false;
  bool t2_vcsbl = false;       // time2_vet_col_source = blend
  void Time2VetSourceInit(ParameterInput *pin, const bool act, const bool col);
  // time2_pred_cap: predicted P = U^n + theta dt K1 with theta = min(1, cap/r), r the
  // relative change of E and eint (0 = off)
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
  // old files (version < 4) carry time2_vet_col = extrap2 channels: read past
  static constexpr int M1_X2_NCH = 24;
  void Time2Report();
  int Time2RstNchWant();        // restart channels this run keeps (0 unless hesdirk2)
  int Time2RstNch();            // ...and writes (0 unless a slope is stored)
  bool Time2RstAccept(const int nfile, const int ver);   // a file's channels fit (sets
                                                         // t2_nfv, t2_rstd3)
  int t2_nfv = 0;               // conservative-ray intensity channels of the file
  int t2_rcl = -1;              // closure channels of the file read (clean_closure =
                                // iterate, CleanRstNcl; -1: this run's own layout)
  int t2_rskip = 0;             // a file's clean_rot_ring channels read with the ring off
  bool t2_rstd3 = false;        // version-1 file: the fkin block is d only (no K1_d)
  int t2_rstv = 5;              // clean file version: <= 3 zero predictor / extrap2
                                // channels and the ring in CLN_N ghosted slabs per
                                // record, 4 the cubed sphere's sc rays' cln_fvi ghosted
  void Time2RstPack(DvceArray5D<Real> &a, int nmb, int m0 = 0, int w0 = 0);
  // restart_pack_per_mb: the restart writer stages the hesdirk2 channels one MeshBlock
  // at a time
  bool rst_pack_mb = true;
  void Time2RstSet(int ch, const HostArray4D<Real> &w, int nmb);

  // ---- off-diagonal Eddington terms and the closure lag (multi-D implicit)
  int impl_offdiag{};           // M1_OD_*: lagged | operator | none
  int od_now{};                 // mode this step runs (impl_offdiag, or none/lagged)
  // implicit_closure_thin_relax = C (0 off; closure_lag = step): relax the step's
  // (chi, n) from the previous step's with w = min(1, 2/(1 + G^2)),
  // G = C max(chi', b/f)/tau_c + f chi'/min(chi, 1-chi); ctr_mem: (chi, n1, n2, n3)
  Real impl_ctrelax{};
  bool ctr_init{};
  DvceArray5D<Real> ctr_mem;
  // <hydro>/rad_signal_speed reads opac(M1_OP_T) and the closure tensor of the last M1
  // step: they travel in the restart (kM1RssRstMagic), else are filled from the initial
  // state before the first dt (RssPrime).
  //! closure tensor the signal speed reads: 0 E/3, 1 M1Chi(F/cE), 2 tau_ten,
  //! 3 vet_cell (chi, n), 4 vet_cell (D_dd), 5 cln_s (D, clean)
  int RssMode() const;
  int RssRstMode() const;   // RssMode + 16 (extended restart state)
  int RssRstNch(int mode) const;  // channels written: opac(M1_OP_T) + the tensor's
  //! array and channel of restart channel n of a mode
  void RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c);
  void RssRstPack(DvceArray5D<Real> &dst, int nmb, int mode, int m0 = 0);
  void RssPrime(bool restart);
  DvceArray5D<Real> rss_stage;   // (nmb, nch, k, j, i) read from the restart file
  int rss_stage_mode = -1;       // mode of rss_stage, -1 = nothing staged

  // ---- transverse realizability limiter, implicit_trans_limit (see ImplicitTransTheta)
  int impl_tlim{};              // M1_TLIM_*: none | lp
  DvceArray4D<Real> thx2, thx3;  // theta_f = 1/(1 + c dt (kt_f + klim_f)) on the x2/x3

  // implicit_vimp: the gas velocity of the enthalpy flux implicit through the solve's
  // radiative force (Newton form)
  bool impl_vimp{};
  // doppler_face = donor | central: central takes the x1 enthalpy flux a E,
  // a = v (1 + chi), as the face mean 0.5 (a_L E_L + a_R E_R) (Athena++ reproduction;
  // forces implicit_vimp = false)
  int dop_face = 0;
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
  bool t2_ffct = false;
  DvceArray5D<Real> t2pf1, t2pf2;   // (m,3,k,j,i) (chat/c) F_applied of K1 / K2's solve
  bool t2pf1_ok = false, t2pf2_ok = false;
  // t2pf is stored with the limiter or with clean staggered (CleanBndAcc): cln_lacc
  // (m,4,k,j) = energy through the physical x1 faces so far, inner (0, + inward), outer
  // (1, + outward), and the old vector's face part (2, 3); cln_ldt = the step's dt
  // (0: none yet, the history reads the faces' F)
  bool t2_pfst = false;
  DvceArray4D<Real> cln_lacc;
  Real cln_ldt = 0.0;
  void CleanBndAcc(int stage, Real c1, Real c2, Real w, Real dt);
  Real t2_ffct_calls = 0.0, t2_ffct_cells = 0.0, t2_ffct_left = 0.0;
  Real t2_ffct_rmax = 0.0, t2_ffct_de = 0.0, t2_ffct_nop = 0.0, t2_ffct_fill = 0.0;
  // t2_ffct_emin = th: the limiter keeps th V E^n in a limited cell; counted: cells
  // left below th E^n, min E/E^n of limited cells
  Real t2_ffct_emin = 0.0, t2_ffct_nemp = 0.0, t2_ffct_qmin = 1.0e300;
  void Time2FaceFct(int stage, Real c1, Real c2);
  bool am_vec = false;             // <hydro>/am_vector: the vector form
  DvceArray4D<Real> amk1, amk2, aminc;   // spin slopes (FSAL / stage 1) and increment

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
  DvceArray5D<Real> opac;       // rho*kappa per cell, M1_NOPAC components
  DvceArray5D<Real> ugas1;      // hydro (IEN, IM1, IM2, IM3) at the start of the substep

  // implicit-solve counters (M1_CNT_*), Reals so that long-run sums cannot wrap
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
  void RSLACheck();
  // ---- implicit transport (rad_vet_implicit.cpp)
  //! read the implicit-solver parameters, allocate its arrays
  void ImplicitInit(ParameterInput *pin);
  // IC radiation start, fresh starts only (rad_vet_icrad.cpp)
  int icr_mode = 0;             // ic_rad: 0 none, 1 fs, 2 relax, 3 fs_relax (IcRadInit)
  // icr_fmeas (set by a pgen): one more IcRadRun pass with the force on the gas; its net
  // x1 momentum deposit dm - dt rho a_ref goes to icr_dm (dt icr_dm_dt), the gas is
  // then restored
  bool icr_fmeas = false;
  DvceArray4D<Real> icr_dm;
  Real icr_dm_dt = 0.0;
  //! spherical-polar wedge: fatal on an unsupported configuration (after ImplicitInit)
  void SphericalS1Check(ParameterInput *pin);
  // spherical-polar or cubed-sphere mesh: the implicit kernels take the Coordinates face
  // areas, volumes and centre distances
  bool sph_geom = false;
  // sph_geom with a non-Eddington closure: radial faces take the integrating factor
  // (M1SphDrr), every face the lagged curvature M1SphCurv
  bool sph_q = false;
  // cubed sphere (sets sph_geom too)
  bool cs_geom = false;
  // cubed-sphere transverse operator (CubedS1Init):
  //   m1bcs     mb_bcs with BoundaryFlag::panel read as block (a seam is an open face)
  //   cs_seam   per MeshBlock: its inner_x2, outer_x2, inner_x3, outer_x3 face is a seam
  //   csg2/3    seam geometry (M1CsSeamGeom) of the x2 / x3 seam faces,
  //             (m, {dth, angm, sfoot, smid}, along-seam cell incl. ghosts, side lo/hi)
  DualArray2D<BoundaryFlag> m1bcs;
  DualArray2D<int> cs_seam;
  DvceArray4D<Real> csg2, csg3;
  // effective transverse geometry of the cs operator for kernels that read Coordinates
  // (ImplicitVimpBuild): two-point x2/x3 distances and the seam areas; x1 copied
  DvceFaceFld4D<Real> cs_dxf_eff{"m1csdxf", 1, 1, 1, 1};
  DvceFaceFld4D<Real> cs_area_eff{"m1csarea", 1, 1, 1, 1};
  // planar unit normals of the x2 faces (m, j = 0..n2, c) and x3 faces (m, k = 0..n3, c),
  // along +x2 / +x3 (M1CsLat0)
  DvceArray3D<Real> cs_fn2, cs_fn3;
  // the same normals at cell centres (m, j or k incl. ghosts, c): lateral axes of the cs
  // half-range ratios
  DvceArray3D<Real> cs_cn2, cs_cn3;
  void CubedS1Init();
  //! seam average of the stored face state f0x2/f0x3
  void CubedSeamFaceAverage();
  bool cs_seam_avg_on = true;
  ParameterInput *pin_cs_ = nullptr;
  //! print the timers and positivity totals of the implicit transport
  void ImplicitReport();
  // dbg_energy_tally (output only): rank 0 prints the global sums of gas energy, E, the
  // reference-work accumulator and esrc at the start and end of every be solve (ETALLY)
  bool dbg_etally = false;
  // dbg_energy_audit (output only): per solve, the change of gas energy + (c/chat) E
  // between the marks of EAuditMark, by cell class (rest, floor gas, reset), with the
  // boundary fluxes; rank 0 prints EAUDIT at mark 5 (mark 6: a discarded attempt)
  // acl_cnt_d: ApplyClosureLimits clips over active cells (0 E raised, 1 |F| scaled)
  DvceArray1D<Real> acl_cnt_d;
  void EventTotals(Real *v);
  //! VET: read <rad_vet>/vet_*, check the mesh, allocate
  void VetInit(ParameterInput *pin);
  //! VET: short-characteristics formal solution and the uniaxial tensor (chi, n)
  void VetShortChar();
  //! VET, vet_tensor = full: guarded D = K/J (and its ghosts) for the solve
  void VetFullTensor();
  //! VET: cost line at the end of the run
  void VetReport();
  //! VET, several MeshBlocks: allocate the banded buffers and neighbour tables
  void VetMBInit(ParameterInput *pin);
  //! VET, several MeshBlocks: exact sweep with per-layer band exchange
  void VetSweepMB();
  //! VET: free the multi-block state
  void VetFree();
  // ---- scheme = clean: the implicit stage solve (rad_vet_clean*.cpp)
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
  const CleanArgs *cln_gdx = nullptr;   // CleanArgs in device memory (slot 3)
  void CleanFvSets(const Real tc, const int cyc);
  void CleanRingMoments();
  void CleanRingShift(const int64_t b);
  int64_t CleanRingBlk(const Real dt) const;
  bool CleanRingTwo() const;
  bool CleanFhBuilt() const;
  void CleanRingPreShift();   // dt_retry: shift before the checkpoint
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
  void CleanStagPcrFactors();
  void CleanRingMatesInit(const int d, CleanRingMates &rm);
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
  // ray history a failed stage's builds change: kept at the stage-loop start, put back
  // with U^n (Time2Stash, Time2Restore) for the backward-Euler redo
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
  // time2_face_fct's two scratch slots of iw (clean: the line-solve slots)
  int Time2FctSlot(const int q) const;
  void Time2FctHalo(const int c);
  int CleanRingRstN() const;            // ring channels in the restart (version 4)
  // vet_rotate_avg_ema: a fill takes 2 sets; cln_rr = EMA offsets e, aux J and D11 and
  // the set count n; cln_rpd = Delta, a_ch and the next set's J and D11 (CleanRingEma)
  int CleanRingK() const {return vra_ema ? 2 : vra_k;}
  int CleanRingNr(const int nf) const {return vra_ema ? nf + 3 : (vra_k - 1)*nf;}
  int CleanRingNd(const int nf) const {return vra_ema ? nf + 3 : nf;}
  void CleanRingRstCopy(DvceArray5D<Real> &a, const int c, const int nmb, const int m0,
                        const bool pk, const bool have);
  // closure components in the restart (extrap history, clean_rot_ring), CleanRstSlot
  // order; the cell layout keeps the first CLN_N - 1, CLN_JR is staggered only
  int CleanRstNq() const;
  // ray intensity channels of the restart: the cs sc rays' (cln_fvi), and their last
  // build when a reflecting top mirrors it (cln_scmem: its top shell)
  int CleanRstNfv() const;
  int CleanRstNcl() const;   // clean_closure = iterate: the closure in place (cln_s)
  void CleanScUnion2();
  void CleanCellHalf(const DvceArray2D<Real> &dv);
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
  // dep: ghost layers its readers need (0 = all ng); fm[f]: components (bit n for c0 + n)
  // the ghosts of face family f take (3 at edges and corners, CleanHaloInit)
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
  bool cln_pcrf = false;        // clean_stag_pcr_factor
  std::vector<DvceArray5D<Real>> cln_pf;   // its PCR factors per mg level (a pass's)
  std::map<std::array<int, 8>, DvceArray1D<int>> cln_pcrl;   // PCR lines of a colour
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
  bool cln_c1b = false;         // clean_clos1_budget
  bool cln_vset = false;        // clean_vel_settled
  bool cln_bmn = false;         // clean_bud_maxnorm
  bool cln_gtf = false;         // clean_gas_tfloor: gas T >= tfloor, energy to/from E
  bool cln_wbfct = false;       // clean_wb_fct
  bool cln_wbff = false;        // clean_wb_floor_fail
  bool cln_wbxl = false;        // clean_wb_exch_local
  Real cln_flrun = 0.0;         // the last write-back's implicit_pos_floor creation
  int cln_flrid = -1;           // its cycle and t2s, 4 ncycle + t2_solve
  bool cln_wbmc = false;        // clean_wb_mcap
  Real cln_wbffn = 0.0;         // its steps accepted at the last cut
  int cln_wbfdg = 0;            // clean_wb_floor_diag (diagnostic)
  int cln_dbginf = -1;          // clean_dbg_inf (test): the cycle of the injected inf
  Real cln_raymin = 1.0e-3;     // clean_ray_min: rays describe a cell when J >= this E
  // clean_ray_cut (CleanRayCut): no rays below the shell where every deeper cell is
  // diffusive (|grad S|/(chi S) < clean_ray_cut_r, none thinner than
  // vet_col_lat_cut_thin), less clean_ray_cut_margin optical depths, at
  // tau_top >= clean_ray_cut_taumin; rises only after clean_ray_cut_wait builds
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
  bool cln_rmsh = false;        // clean_rst_mshare: M1 share memory at an old restart
  bool cln_rmdone = false;      // its start taken
  bool cln_linpass = false;   // clean_lin_pass: a pass goes on past a failed linear solve
  void CleanStagXDonor(const bool first, const bool lat = true, const bool hold = false);
  void CleanWorstCell(const Real dnl, const Real dnl0);
  bool cln_srccvx = false;      // clean_source_convex: vet_scatter S between J and B
  bool cln_freal = false;       // clean_stag_freal: realisable sb, restart faces
  // clean_stag_xdonor: cross flux of the iterate as a donor-cell coefficient of its face
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
  std::array<CleanRingMates, 3> cln_fwr;   // (fwb: the lines' MeshBlocks)
  DvceArray1D<Real> cln_fwbc;    // (fwb: the faces, their momentum sums |dp|, x1 dp)
  bool cln_fcone = false;   // clean_stag_fcone
  bool cln_fcgas = false;   // clean_stag_fcone_gas
  bool cln_fclz = false;    // clean_stag_fcone_lz
  bool cln_fcemx = false;   // clean_stag_fcone_emax
  bool cln_lzsp = false;    // clean_lz_spin
  DvceArray5D<Real> cln_fcx;   // (clean_lz_spin: fcone's cut withheld, by face)
  DvceArray2D<Real> cln_cx[2];   // (clean_lz_spin, Cartesian: the cell centres x, y)
  bool cln_fwbst = false;   // clean_stag_fwb_stage
  bool cln_fwbpj[3] = {false, false, false};   // fwb_line, -, _ring (fwb_cart: all)
  int cln_fwbs = -1;        // (fwb_stage: the face slot of s; cln_fwbok: s is the round's)
  bool cln_fwbok = false;
  bool cln_fwbt2 = false;   // clean_stag_fwb_t2
  DvceArray1D<Real> cln_fconc;   // (fcone: the faces scaled, their momentum sum)
  int cln_xdfz = 0;         // clean_stag_xdonor_freeze
  bool cln_xdh = false;     // (its hold in this pass)
  // clean L_z ledger (clean_cs_lz_abs = 2; M1LzKey): cln_lzl the L_z the gas does not
  // take, booked as the radiation's; its sums as the gas momentum's slopes K1 (cln_lzs1),
  // K2 (cln_lzs2), stage-start increment (cln_lzc) and stage-1 part pending (cln_lzp).
  // 3-D spherical polar: the L_z that crossed the radial boundaries.
  Real cln_lzl = 0.0, cln_lzs1 = 0.0, cln_lzs2 = 0.0, cln_lzc = 0.0, cln_lzp = 0.0;
  ParameterInput *cln_lzpin = nullptr;
  void CleanLzBook(const Real s, const bool start);
  Real CleanLzTotal();   // gas + radiation + ledger L_z
  bool cln_xtnorm = false, cln_xtcons = false;   // (CleanFirstIterate, CleanSolve)
  int cln_xte = 0;   // clean_xt_exact: 0 off, 1 cons, 2 x2
  int cln_gnsw = 0, cln_gnswmin = 20;   // clean_gas_newton_switch, _min
  Real cln_stfac = 0.5;   // clean_stall_fac
  bool cln_gpic = false;   // a stalled solve's passes without the gas Newton slope
  Real cln_ngnsw = 0.0;   // solves switched so far
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
  Real cln_nl_escale = 0.0;     // E scale of clean_nl_efloor: max E of the iterate
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
  int cln_stagmg = 0;           // clean_stag_mg_color: 1 line colours (j + k) parity,
                                //   2 mirror
  int cln_pcsweep = 0;          // clean_pc_sweep: +1 a 2nd red pass, +2 a post-sweep
  bool cln_steady = false;      // clean_init = steady: radiation relaxed before step 1
  int cln_steady_n = 50;        // clean_init_n: most pre-solve passes
  Real cln_steady_dtfac = 1.0e6, cln_steady_tol = 1.0e-8;  // clean_init_dtfac, _tol
  std::vector<DvceArray5D<Real>> cln_mg;          // coarse levels of the multigrid
  std::vector<std::array<int, 3>> cln_mgn;        //   and their cells along x1..x3
  std::array<int, 3> cln_mgt = {0, 0, 0};       //   its tiles' cells (layout_tile)
  int cln_mgnl = 2;             //   most coarse levels (clean_stag_mg_levels)
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
  bool cln_rfh = false;         // cs cell half ranges in the ring (CLN_HO, CleanCellHalf)
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
  DvceArray5D<Real> cln_rnew;   // jump: moments of block b+1's window at the boundary
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
  // staggered: cell flux of the faces' energy flux Phi (m,3,k,j,i; output m1_phi); u0 F
  // holds the face momenta's
  DvceArray5D<Real> cln_phic;
  DvceArray5D<Real> cln_rstd;   // staggered fkin: the restart's d and its K1 (m,6,k,j,i)
  bool cln_fref_set = false;    // cln_fref holds the force of the first coupled solve
  bool cln_mdcap = false;       // wb_arad_md: the write-back captures cln_fref (steady)
  bool cln_vzero = false;       // the iterates keep the gas momentum (steady start)
  // cubed sphere: face areas, centroid distances across the faces, normals and centroids
  DvceArray4D<Real> cln_ga[3], cln_gd[3];
  CleanGeo cln_gn;
  DvceArray3D<Real> cln_tg;   // spherical polar: CleanArgs::tg
  // vet_gd_theta_seam_sym: CleanArgs' geometry with the theta-periodic seam's ghost rows
  // their periodic images' (CleanSeamImage): vol, area x1f / x2f, dxface x2f / x3f,
  // x2v, tg
  DvceArray4D<Real> cln_sgv, cln_sga0, cln_sga1, cln_sgd1, cln_sgd2;
  DvceArray2D<Real> cln_sgx;
  DvceArray3D<Real> cln_sgtg;
  void CleanSeamImage();
  struct CleanHalo {            // ghost exchange of nc components: scratch array and
    DvceArray5D<Real> a, ca;    // its coarse dummy
  };
  std::map<int, CleanHalo> cln_halo;   // by nc
  DvceArray1D<char> cln_xdev;   // device copies of kernel arguments (CleanDevCopy)
  // exchange object of every clean ghost fill (a fill sends only its own variables);
  // fills deferred to the next exchange (CleanGhostDefer)
  std::unique_ptr<MeshBoundaryValuesCC> cln_hpb;
  std::vector<CleanGSeg<Real>> cln_hpend;
  DvceArray5D<Real> cln_hbuf;   // storage of every scratch array (<= 32)
  // clean_halo_direct: cell lists read once off cln_hpb's exchange (CleanHaloInit);
  // on-rank ghosts by one copy, one message per neighbour rank, receives posted first
  struct CleanHaloList {
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
  CleanHaloList cln_hdl;
  void CleanHaloInit(const int n3, const int n2, const int n1);
  template <typename T> void CleanHaloPost(const std::vector<CleanGSeg<T>> &sg);
  template <typename T> void CleanHaloFinish(const std::vector<CleanGSeg<T>> &sg);
  template <typename T>
  std::array<uint32_t, 5> CleanHaloKey(const std::vector<CleanGSeg<T>> &sg);
  CleanHaloList::Mask &CleanHaloMask(const std::array<uint32_t, 5> &key);
  Real cln_nsolve = 0.0, cln_nit = 0.0, cln_nfail = 0.0, cln_nneg = 0.0;

  // ---- closure = tau / vet_col column tensor (rad_vet_tau.cpp)
  bool tau_closure{};          // = vet_col: the column tensor tau_ten is in use
  bool tau_ready{};            // VetColInit has run
  DvceArray5D<Real> tau_ten;   // (m,4,k,j,i): chi, n1, n2, n3 read by ImplicitSolve

  // ---- closure = vet_col (rad_vet_vetcol.cpp): f_K = K/J of a 1-D formal solution per
  // radial column (impact-parameter rays on sp, plane-parallel Gauss rays on Cartesian),
  // a fixed uniaxial tensor about r_hat (or the M1 flux axis) through tau_ten
  bool vet_col{};              // closure = vet_col (default false)
  bool vcol_sph{};             // spherical rays (sp) or plane-parallel (Cartesian)
  bool vcol_axis_flux{};       // vet_col_axis = flux: n = the M1 cell flux direction
  int vcol_nc{}, vcol_np{}, vcol_nmu{}, vcol_every{}, vcol_nray{};
  // vet_col_fk_min: lower bound of f_K = K/J, 1/3 or 0 (D = diag(f, (1-f)/2, (1-f)/2)
  // stays realizable for 0 <= f <= 1)
  Real vcol_fkmin = 1.0/3.0;
  // vet_source = relaxed | relaxed_esc | gas (alias vet_col_source): checked and echoed
  // only (the clean rays take the state's own source)
  bool vcol_srelax = true;
  bool vcol_sesc = false;
  Real vcol_sx0 = 1.0;         // time2_vet_col_source_x0 (blend)
  // vet_scatter: electron scattering, the opacity split (otab.scat, opac M1_OP_S) and
  // VetScatterMix; vscat_form 1 ma (Jiang 2021 eq. 6), 2 absorption
  bool vscat = false;
  int vscat_form = 1;
  // vet_scatter_j: false current (J = E^n / iterate, absorption form only), true
  // relaxed (J = E*, the end-of-step pair of B*)
  bool vscat_jrel = false;
  // vet_col with a reflecting outer x1: a second sweep with I_in = b/(1 - a) per ray
  // (b the vacuum-top outgoing intensity, a the round-trip transmission)
  bool vcol_rtop = false;
  // vet_col_surface_face (sp, with vet_col_surface_q and implicit_marshak_face = linear):
  // outer face flux c H(face) (E/J)_face with E, J extrapolated linearly in r^2 (second
  // order). vcol_jca, vcol_jcb: top-face weights (M1SphMarshakCoef), vcol_jcap the
  // limiter's 2 (r_c/r_f)^2
  bool vcol_sqf = false;
  Real vcol_jca = 0.0, vcol_jcb = 0.0, vcol_jcap = 0.0;
  bool vcol_built{};           // a tensor exists (vet_col_every > 1 skips builds)
  Real vcol_time{}, vcol_ncall{}, vcol_nskip{};
  // vet_col_surface_q: outer-x1 Marshak q = H(top face)/J(top cell) per column
  bool vcol_sq{};
  Real vcol_qmin{}, vcol_qmax{};
  int vcol_lc{};               // shells per chunk of the team build (scratch budget)
  DvceArray3D<Real> vcol_q;    // (m, k, j): the column's Marshak q at the top face
  void VetColInit();           // checks, tau_ten, vcol_q (first call)
  void VetColReport();

  // ---- vet_col_lat (rad_vet_vetcol.cpp): lateral correction of vet_col from a few-angle
  // SC sweep, each direction swept through the 3-D field and through the column's
  // laterally homogeneous twin; their moment difference (0 on a uniform state) is stored
  // in tau_ten M1_TT_LAT0.. (mesh basis). dD_rr enters the diagonal, D_r,lat the face
  // equations (M1SphLat).
  bool vlat_on = false;        // vet_col_lat
  // Athena++ <mesh>/x2rat, x3rat named != 1 (AthenaK grids are uniform): VetLatInit
  // refuses them under vet_col_lat / vet_gd
  std::string vlat_rat_keys;
  // vet_col_lat_offdiag: D_r,lat in the solve, 0 none, 1 lagged Picard term, 2 operator
  // (Krylov operator through the stencil and VetLatOp)
  int vlat_odm = 2;
  bool vlat_now = false;       // the term is on in this step (positivity fallback)
  Real vlat_odmax = 0.0;       // max |D_r,lat| of the last sweep (all ranks)
  // operator form active (a zero D_r,lat leaves the operator untouched)
  bool VlatOp() const {
    return vlat_now && !vlat_lagfb && (vlat_odm == 2) && (vlat_odmax > 0.0);
  }
  // vet_col_lat_posfix = lag: at a positivity trigger the lateral terms leave the
  // operator (M-matrix form) but stay as the lagged Picard term; a second trigger in the
  // step drops them (vlat_nfall)
  bool vlat_posfix_lag = false;
  bool vlat_lagfb = false;     // this step: the lateral terms lagged after a trigger
  // implicit_fb_diag, implicit_fb_local (rad_vet_meshfb.cpp)
  bool fbl_any = false;
  DvceArray5D<Real> vlat_tt_loc;
  DvceArray5D<Real> VlatTT() const {return fbl_any ? vlat_tt_loc : tau_ten;}
  // theta of each x2 row for M1SphTan (x2v; the source cell's theta in periodic ghost
  // rows), built on first use
  const Real *vtan_x2t_src = nullptr;
  int vlat_icut = 0;           // first shell of the sweep (CleanRayCut)
  Real vlat_cut_thin = 0.0;    // vet_col_lat_cut_thin: no cut above a cell thinner than
                               // this (0 off)
  bool vlat_ready = false;     // VetLatInit done
  int vlat_nbuild = 0;         // builds done
  Real vlat_time = 0.0, vlat_ncall = 0.0, vlat_nclamp = 0.0;
  DvceArray5D<Real> vlat_cs;   // (m, 2, k, j, i): ln chi, ln S (ghosts exchanged)
  DvceArray5D<Real> vlat_cs_c;
  DvceArray1D<Real> vlat_cnt;  // (1): lateral reads clamped to the ghost band
  MeshBoundaryValuesCC *pbval_ea = nullptr;   // EAUDITF diagnostic
  MeshBoundaryValuesCC *pbval_vs = nullptr;   // vlat_cs
  void VetLatInit();
  void VetLatSweep();   // ln chi, ln S of a vet_gd build
  void VetLatOdMax();
  void VetLatExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                      MeshBoundaryValuesCC *pb);
  void VetLatDump();


  // ---- vet_gd (rad_vet_vetgd.cpp): Eddington tensor from an SC sweep along globally
  // fixed (Cartesian) HEALPix directions, through the vet_col_lat interface (sets
  // vlat_on); sp wedge / cubed sphere
  bool vgd_on = false;         // vet_gd
  // HEALPix nside of a set that names none (12 nside^2 directions)
  static constexpr int vgd_nside = 3;
  // vet_gd_tangential: sp the lagged M1SphTan term; cs Q in LAT4, LAT5 about the cell
  // tangent frame cs_qf (m, 6, k, j incl. ghosts; M1CsQab)
  bool vgd_tan = true;
  DvceArray4D<Real> cs_qf;
  bool vgd_replace = false;    // vet_gd_replace: no vet_col sweep (D = I/3 below the
                               // first gd shell, q from the gd top cell)
  int vgd_n = 0;
  int vgd_ls = 0;              // vet_gd_ls = 6/8/10/12: level-symmetric LQ_N set instead
  double vgd_zper = 0.5*M_PI;  // z-rotation period of the set (VetGdInit)
  int vgd_gl_nmu = 0;          // vet_gd_gl_nmu > 0: product set GL(nmu in n_z) x nphi
  int vgd_gl_nphi = 12;        // vet_gd_gl_nphi (C_nphi about z)
  bool vgd_gl_stag = true;     // vet_gd_gl_stagger: alternate rings half a step
  int vgd_rot = 0;             // vet_gd_rotate_every: z-rotation every N cycles
  // vet_gd_rotate_dt > 0: rotate in physical time, block
  // b = floor((t + dt/2 - t0)/dt_rot) (t0 = vet_gd_rotate_t0); vet_gd_rotate_align: steps
  // land on block boundaries (VetGdAlignDt)
  Real vgd_rot_dt = 0.0, vgd_rot_t0 = 0.0;
  // vet_gd_rotate_stretch: VetGdAlignDt may lengthen dt by <= this
  Real vgd_rot_stretch = 0.0;
  // vet_gd_rotate_dt = auto: vgd_rot_auto_fac x the first dt, recorded in pin_rot_
  bool vgd_rot_auto = false;
  // the unnamed default resolves from the uncapped first CFL dt (Mesh::dt_uncapped)
  bool vgd_rot_unc = false;
  Real vgd_rot_auto_fac = 20.0;
  ParameterInput *pin_rot_ = nullptr;
  bool vgd_rot_align = false;
  double vgd_rot_nalign = 0.0;   // steps shortened by VetGdAlignDt (diagnostic)
  // vet_gd_rotate_slope = be | keep | correct: what the hesdirk2 FSAL slope K1 (of the
  // old set) does on the first step of a new rotation block; be: a backward-Euler step
  int vgd_rot_slope = 0;
  int64_t vgd_rot_lastblk = INT64_MIN;
  void VetGdRotStep();
  // vet_gd_rotate_slope = correct: a correction solve (Time2RotCorrect) of size
  // dt_c = vet_gd_rotate_slope_frac x dt from U^n with the new set turns K1 into
  // K1 + (I - dt_c A_new)^{-1} (f_new - f_old)(U^n); vet_gd_rotate_slope_reuse: the
  // stage-1 solve keeps that closure build
  bool t2_rotcorr = false;
  bool t2_vcreuse = false;
  bool vgd_rot_creuse = true;
  Real vgd_rot_cfrac = 0.0;
  double vgd_rot_ncorr = 0.0, vgd_rot_ncfail = 0.0, vgd_rot_ncpass = 0.0;
  bool Time2RotCorrect(Driver *pdrive, Real dt);
  Real vgd_alpha = -1.0;       // current z-angle (-1: tables not built)
  // vet_gd_rotate_axis = z | 3d (default z on sp, 3d on the cubed sphere). 3d: the set of
  // block b is R_b n, R_b a 3-D rotation about a quasi-random axis (VgdRot3); vgd_alpha
  // is then the block key, not an angle
  bool vgd_rot3 = false;
  void VgdRot3(const Real key, double r[9]) const;
  std::vector<double> VgdRotate(const Real alpha) const;
  // vet_rotate_avg = K >= 2 (rad_vet_rotavg.cpp): the closure is the mean over the last K
  // rotated direction sets. A closure registers its normalised fields (ratios to J); the
  // mean is kept as offsets d_j = g_{b-j} - g_b, each from two sweeps of the same state:
  // used = g_b(now) + sum_j d_j / K. Needs a rotating closure.
  static constexpr int M1_RA_NF = 16;
  struct RotAvgTab {           // registered fields (POD: captured by the kernels)
    int nf = 0;
    int arr[M1_RA_NF]{}, comp[M1_RA_NF]{}, kind[M1_RA_NF]{};
    Real lo[M1_RA_NF]{}, hi[M1_RA_NF]{};
    int farr = -1, fcomp = 0;  // validity flag field: invalid where < 0 (-1: none)
    int jf = -1;               // J field (index into the table) of kind-3 fields
    int jw = -1;               // the ratio field weighting the sets (-1: equal)
  };
  int vra_k = 0;               // vet_rotate_avg (K; 0/1 off)
  // vet_rotate_avg_ema = N: the K-window mean replaced by an exponential moving average
  // over the rotated sets, weight a = 1/N per block, in the same offset form (one slot
  // e = avg - g_b(now), r = Jbar/J_b, + an aux slot; vra_k = 3). _tau > 0 (physical time)
  // or _lc > 0 (x the light-crossing time r/c): a = max(1/N, 1 - exp(-dt_blk/tau)).
  bool vra_ema = false;
  Real vra_ema_n = 0.0, vra_ema_tau = 0.0, vra_ema_lc = 0.0;
  void RotAvgEmaInit(const int ilo);
  void RotAvgEmaShift(const int ds);
  void RotAvgEmaChange(const int ilo);
  Real vra_ema_eps = 0.05;     // vet_rotate_avg_ema_change (0: no change test)
  bool vra_ema_chj = true;     // vet_rotate_avg_ema_change_j: the J ratio in the test
  bool vra_ema_chstep = false; // vet_rotate_avg_ema_change_form = step | linear
  // vet_rotate_avg_ema_wb (clean ring): a = a_ch alone (no 1/N floor, no 1/(n + 1)
  // count): a static state keeps its closure
  bool vra_ema_wb = false;
  int vra_ema_diag = 0;        // vet_rotate_avg_ema_diag = M: print the change weights
  // vet_rotate_avg_ema_dema: double EMA, slot 2 the EMA of the EMA (e2); the closure uses
  // 2 e1 - e2 (ratios 2 r1 - r2 in [r1/4, 4 r1]); vra_k = 4 (slot 3 the Delta)
  bool vra_ema_dema = false;
  int64_t vra_ema_nsh = 0;     // boundaries seen (diag cadence only)
  void RotAvgEmaApply(const int ilo);
  // vet_gd_rotate_xfade (needs vet_gd_rotate_dt > 0): clean ring live, the closure
  // cross-faded from the current to the next set over each block (CleanFvSets); both
  // sets' tables on the device (VgdXfTables)
  bool vgd_xf = false;
  // vgd_xff: both sets of a cross-fade build in one sweep (one launch per shell, one
  // message per partner), bitwise the two sweeps. The second set keeps its own
  // intensities (vgd_i_xb), shell and halo lists. Falls back where VetGdXffOk fails.
  bool vgd_xff = false;
  bool vgd_xff_now = false;      // this sweep carries the next set too (VetGdSweep)
  bool vgd_xff_prep = false;     // VetGdSweep: the set's lists only, no passes
  DvceArray1D<int> vgd_shl_xb[2];
  int vgd_shn_xb[2] = {0, 0};
  Real vgd_shl_alpha_xb = -3.0;
  DvceArray1D<char> vgd_bbufx;     // fused shell closure (device copy)
  DvceArray1D<int> vgd_cslst_xb;   // cs: the other set's pass lists
  DvceArray1D<Real> vgd_xff_sb, vgd_xff_rb;   // merged messages (sp)
  double vgd_xff_n = 0.0;          // fused builds (diagnostic)
  double vgd_xff_nmsg = 0.0;       // merged exchanges (diagnostic)
  double vgd_xff_nfb = 0.0;        // builds that fell back to two sweeps (diagnostic)
  bool VetGdXffOk();
  void VgdXffSwap();
  void VgdEs2Make();   // second stream of the cubed sphere's early send
  void VgdGbSwap() {   // the compact global band's kept lists go with their set
    std::swap(vgd_gbcf, vgd_gbcf_xb);
    std::swap(vgd_gbcb, vgd_gbcb_xb);
    std::swap(vgd_gbrec, vgd_gbrec_xb);
  }
  bool VetGdHlXff(const int i, const int ws);
  bool VetGdGbHaloX(const int i, const int ws);   // vgd_raw: both sets, one message
  // early send (VetGdCsSweep): the shell's send cells first on a second high-priority
  // stream; the host exchanges while the shell kernel (skipping them) runs
  DevExeSpace vgd_es2;
  bool vgd_es2_made = false;
  DvceArray4D<int> vgd_est[2];     // per set (m, kb, jb, d): stamp of the early cells
  int vgd_esid = 0;                // stamp of the current early shell
  Real vgd_xf_alt = -9.0;        // angle of the tables in the second device set
  void VgdXfTables(const Real alpha);
  bool vra_eq = false;        // vet_rotate_avg_weight = equal (default j: J-weighted)
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
  // history in the restart file (kM1RaRstMagic, rad_vet_rotavg.cpp)
  static constexpr int M1_RA_RST_NH = 11;   // int64 header entries (10)
  int RotAvgNf() const;
  int RotAvgRstNslab(const int ncell) const;
  int RotAvgRstI0() const;     // first written radial cell
  void RotAvgRstHdr(int64_t *h) const;
  void RotAvgRstPack(DvceArray5D<Real> &dst, const int nmb);
  void RotAvgRstPackHost(HostArray5D<Real> &out, const int nmb, const int nsc);
  bool RotAvgRstLoad(const int64_t *h, const DvceArray5D<Real> &src, const int nmb);
  bool RotAvgRstCheck();
  int vra_rst_nf = -1;         // nf of a history loaded from the file (-1: none/checked)
  int64_t vra_blk = INT64_MIN; // block of the last build
  bool vra_last = false;       // this build is on the last step of its block
  DvceArray5D<Real> vra_bk_vc;              // vet_sc: vet_cell kept over the next-set
                                            // sweep of a last build
  void VetGdSweepPost();       // twin, sweep(s), fused moments, VetGdPost
  bool vgd_tbath = false;   // the top (x1max Marshak) bath enters the sweep
  std::vector<double> vgd_base;   // unrotated set (x, y, z, w)
  DvceArray2D<Real> vgd_dir;   // (d, 0..3): n_x, n_y, n_z, weight (sum 1)
  // (m, d, k, j, i) intensities, ghosts lagged, as a ragged per-shell band: shell i keeps
  // only the lateral band it needs (depth wi(i) <= vgd_w) and only on sides whose
  // neighbour is off-rank. Index space (m, d, k, j, i) with k, j in the band space of
  // depth vgd_w; one slab per (m, i), d fastest, then k, then j. An access outside the
  // slab goes to a dummy element and is counted (VgdRag::oos, must stay 0).
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
  // memory arena (with vet_gd): the intensity arrays and the solve scratch dead across a
  // VetGdBuild in one allocation
  bool marena_on = false, marena_sealed = false;
  DvceArray1D<Real> marena;   // intensities, then (first checkpoint) solve scratch
  int64_t marena_sw = 0;      // Reals of the intensity arrays bound so far
  void MemArenaInit(ParameterInput *pin, int nmb, int niw, int n3, int n2, int n1);
  bool MemArenaVgd(VgdRag &a, int64_t n);
  void MemArenaGrow(int64_t need, bool keep = true);
  void MemArenaSeal();
  void CleanWoBind();   // cln_w.o, .s at the first solve (mem_arena)
  // vet_gd_band_ring (MPI band with exact halo lists, or cs compact lists): the ragged
  // array keeps each shell's interior only and lends a shell one of VGD_RING_K slabs of
  // its full band depth while the sweep needs its band (outward pass: interior copied in,
  // the saved inward-exchange band put back). Only off(m, i) / wi(i) switch, so every
  // value read is the dense band's (bitwise).
  static constexpr int VGD_RING_K = 3;
  int vgd_rg_k = VGD_RING_K;          // slots of this allocation
  struct VgdRingOp {
    int kind = 0;      // 1 save (inward band -> store), 2 release, 3 acquire
    int i = -1, slot = -1;
    bool copyin = false, restore = false;
  };
  bool vgd_ring = false;              // ring layout allocated (both arrays)
  bool vgd_ring_named = false, vgd_ring_par = true;   // vet_gd_band_ring (named)
  bool vgd_ring_act = false;          // the ring runs in this sweep
  bool vgd_rg_two = false;            // this sweep carries a second array (vgd_i_xb)
  bool vgd_rg_xb = false;             // the other xfade set's array
  std::vector<int64_t> vgd_rg_sn2;    // (i) the other set's saved count
  bool vgd_rg_lastx = false;          // the last halo used an exact list
  int vgd_rg_lastx_e = -1;
  std::vector<int64_t> vgd_rg_offc;   // (m c1 + i): compact slab offsets (host)
  std::vector<int64_t> vgd_rg_soff;   // (s nmb + m): slot offsets (host)
  std::vector<int> vgd_rg_wfull;      // (i): a shell's band depth when lent a slot
  std::vector<int> vgd_rg_state;      // (i): -1 compact, else slot
  std::vector<int> vgd_rg_sown;       // (s): shell holding slot s, -1 free
  std::vector<int64_t> vgd_rg_so, vgd_rg_sn;   // (i): save store offset, count/array
  std::vector<int> vgd_rg_sk;         // (i): save kind: 1 exact list, 2 dense band, 3 raw
  std::vector<int> vgd_rg_se;         // (i): exact list record (3: the raw key)
  // vgd_raw: band positions (as exact lists) each raw exchange writes, keyed as vgd_gbcf
  // (VetGdRawRec)
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
  // vet_gd_wall_interp: wall ghosts from the 3 nearest same-branch directions
  // (m, k, j, d, 3)
  bool vgd_wint = false;
  // vet_gd_wall_sym (needs _wall_interp): the set turned so a mirror plane is the wedge
  // mid-plane, phi-wall ghosts linear in azimuth on their ring (VgdMirPick)
  bool vgd_wsym = false;
  bool vgd_wsx = false;       // and its phi wall map inexact: both-frame seam faces
  Real vgd_wsw = 0.0;         // phi period W
  DvceArray1D<int> vgd_mir;   // (n): mirror image of each direction
  // vet_gd_theta_seam_sym (needs _wall_interp): clean theta seam faces in both frames
  bool vgd_wst = false;
  DvceArray1D<int> vgd_mirz;  // (n): the equatorial mirror image of each direction
  VgdLG vgd_lg;   // the vet_gd sweep geometry from the global index
  Kokkos::View<int *****, LayoutWrapper, DevMemSpace> vgd_m3;
  DvceArray5D<Real> vgd_w3;
  // vet_gd_theta_wall = periodic | mirror: mirror makes the theta walls specular; a theta
  // ghost at depth g reads the interior row at depth g across the wall, direction
  // (n.r, n.theta, n.phi) -> (n.r, -n.theta, n.phi), from up to 3 same-branch directions
  // (vgd_wall = 2, vgd_map / vgd_m3 / vgd_w3); ln chi, ln S from the mirror rows
  // (VetGdThMirCs)
  bool vgd_thmir = false;
  bool cln_scmem = false;      // clean sc rays, reflecting top: its mirror reads cln_top
  DvceArray5D<Real> cln_top;   // its last top shell: slot s, direction d at (s*nt + d/n1,
                               // i = d%n1), nt = ceil(n/n1) (VetGdTopCopy)
  void VetGdTopCopy();
  int cln_rq = 0;               // clean_rot_ring, sc rays: slot swept (CleanFvSets)
  int cln_rqm = 0;              // its sp top-mirror slot (blk mod 2)
  void VetGdThMirCs(DvceArray5D<Real> &a, const int i0, const int i1);
  // vet_gd_theta_wall = vacuum (1) | extrap (2) at a physical theta wall (flag 3 / 4):
  // vacuum: I = 0 entering at the wall (ghost reads 0); extrap: ghosts read the edge
  // row's intensity (zero gradient); ln chi, ln S from the edge row (VetGdThVacCs)
  int vgd_thvac = 0;
  void VetGdThVacCs(DvceArray5D<Real> &a, const int i0, const int i1);
  MeshBoundaryValuesCC *pbval_gd = nullptr;
  // exact per-shell lateral halo: vgd_i and vgd_cs carry a ghost band of vgd_w cells,
  // filled shell by shell from the 8 lateral neighbours (on-rank copies + MPI)
  int vgd_w = 0;
  // clean closure on a wedge: band of any depth keyed on the global cell
  // (VetGdGlobalBandInit); the sweep reads an on-rank owner (vgd_own) directly, a remote
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
  bool vlat_src_ep = false;
  Real vgd_nrb = 0.0, vgd_trb = 0.0;
  void VetGdIterRebuild();
  // cs half-range data at the start of a solve (|H|/J, seam ghosts)
  void VetGdCsHrPrep();
  void VetGdSpHrHalo();   // sp ghosts of vgd_hr
  MeshBoundaryValuesCC *pbval_hr = nullptr;
  DvceArray5D<Real> vgd_hr_c;
  void VetGdRealDiag();
  void VetGdHaloInit();
  // vet_gd on the cubed sphere (rad_vet_vetgd_cs.cpp): a block's band is its own panel's
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
  void VetGdCsInit();
  void VetGdCsSource();
  void VetGdCsSweep();
  void VetGdCsTab(const std::vector<double> &d);   // into vgd_tabh (walls off, mr)
  template <class V>
  void VetGdCsHalo(V &a, const int nv, const int i0, const int i1, const bool loc,
                   const int ws = -1);
  Kokkos::View<int ****, LayoutWrapper, DevMemSpace> vgd_csown;   // (m, k, j, 3) owner
  // vet_gd_halo_compact > 0 on cs: chart-aware exact lists, per shell and pass only the
  // (remote band cell, direction) values the sweep reads (VetGdCsCompact: a read-mask
  // kernel replays the upwind stencils; masks sent to the owners once per set). Every
  // remote read is checked against the mask: fatal on a miss.
  struct VgdCsCl {
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
  static constexpr int M1_VGD_CSCL_MAX = 16;
  VgdCsCl vgd_cscl[M1_VGD_CSCL_MAX];
  // read masks of two entries (the set in use and the next); the others keep lists only
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
  void VetGdX1Hand(const bool up);
  void VetGdX1Pass(const bool up);
  DvceArray4D<int> vgd_rmap;   // (m, k, j, d): mirror direction at a reflecting top
  void VetGdTopMirror();
  DvceArray3D<Real> vgd_isc;   // (m, k, j): flux ratio of a reflecting inner x1
  void VetGdInnerMirror();
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
  void VetGdHcPrep(const int slot, const int i0, const bool inw0, const int ws);
  Kokkos::View<int **, LayoutWrapper, DevMemSpace> vgd_pbd;   // (ws, 2 (np+1))
  std::vector<double> vgd_r1v, vgd_r1f;
  template <class V>
  bool VetGdHaloCompact(V &a, const int nv, const int i0, const int ws);
  // compact-halo lists per (pass, shell), built from the masks at the first exchange of a
  // mask set (direction set, cut, depth) and reused: gather pack, MPI, scatter (+ zeros
  // on the other pass's entries) instead of flag / scan / pack / expand / dense unpack
  struct VgdHlEntry {
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
  // per band array and pass: build id of the entry whose exchange last wrote each shell's
  // band (0: never written, all zero; -1: unknown)
  std::map<const void *, std::array<std::vector<int64_t>, 2>> vgd_hl_gw;
  bool VetGdHlBuild(const int e, const int i, const bool inw, const int ws,
                    const DvceArray1D<int> &fs, const DvceArray1D<int> &fr,
                    const std::vector<int> &pbh);
  DvceArray4D<Real> vgd_trig;   // (m, k, j, 4): sin, cos of theta, phi (shell kernel)
  // vet_gd on a stretched theta and/or the full sphere (rad_vet_vetgd_sph2.inc).
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
  // vgd_cs*, VetGdCsHalo; any depth), per shell only within its halo depth; the upwind
  // stencil is clamped to the cell's own window (k +- wu, j +- wu), not the block's band
  bool vgd_gband = false;
  DvceArray2D<Real> vgd_phc;              // (m, band k): band-column phi, global index
  std::vector<double> vgd_phc_h;          // host copy, m-major
  DvceArray1D<int> vgd_gdps, vgd_gdpr;    // depth of every send / receive entry
  std::vector<int> vgd_gdnl;              // (ws): on-rank entries with depth <= ws
  std::vector<std::vector<int>> vgd_gdsn, vgd_gdrn;   // (ws)(partner): send / recv counts
  void VetGdGlobalBandInit();
  // vet_gd_halo_compact > 0 with vet_gd_band_global: per shell and pass only the values
  // VgdHaloKeep passes at the source cell's global (J, K); rows in the polar zone
  // (colatitude < vet_gd_band_colat) + the band depth keep their whole branch
  bool vgd_gbc = false;
  DvceArray2D<int> vgd_gbjs, vgd_gbjr;    // (entry, 2): source cell global J, K
  DvceArray1D<Real> vgd_gbth, vgd_gbdt;   // (J): node theta, theta cell width
  DvceArray1D<int> vgd_gbcap;             // (J): 1 = polar zone (branch only)
  DvceArray1D<int> vgd_gbfs, vgd_gbfr;    // flags -> positions (entry*n + v; +1)
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_gbpb;   // segment bounds
  // vgd_raw: scanned flags and bounds per (shell, pass, depth) (VetGdGbHaloC); _xb: the
  // other cross-fade set's
  std::map<int, std::array<DvceArray1D<int>, 2>> vgd_gbcf, vgd_gbcf_xb;
  std::map<int, std::vector<int>> vgd_gbcb, vgd_gbcb_xb;
  Real vgd_gbdmin = 0.0;
  template <class V>
  bool VetGdGbHaloC(V &a, const int i, const int ws);
  void VetGdThetaInit();                  // node table and checks
  double VetGdThetaDmin() const;          // narrowest lateral spacing for the band
  int VgdPoleSlot(const int o, const int lx2) const;  // 0, 1 north, 2 south
  Kokkos::View<int64_t*, DevMemSpace> vgd_scan_sums;   // VgdScan2 chunk sums
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
  std::vector<VgdHlEntry> vgd_hxe;         // (2 n1): exact lists per (pass, shell)
  std::vector<VgdHlEntry> vgd_hxe_alt;     // the other xfade set's lists
  bool VetGdHxKey(const int e, const int ws) const;   // record e matches the masks' key
  bool VetGdHxCut(const int e, const VgdHlEntry &X) const;   // ray-cut move test
  void VetGdHxBuild2(const int e, const int i, const int ws);
  // exact lists of many shells in a few collectives (VetGdHxBatch): scratch pools,
  // per-depth partner / piece tables, pinned bounds
  DvceArray1D<char> vgd_hxb_tmp;
  int vgd_hxb_g = 0, vgd_hxb_npc = 0;
  size_t vgd_hxb_rcap = 0, vgd_hxb_scap = 0;
  DvceArray1D<int> vgd_hxb_nd, vgd_hxb_f0, vgd_hxb_f1, vgd_hxb_rq, vgd_hxb_so, vgd_hxb_oc;
  Kokkos::View<int**, LayoutWrapper, DevMemSpace> vgd_hxb_pst, vgd_hxb_pcs, vgd_hxb_pci,
                                                  vgd_hxb_pds;
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> vgd_hxb_hb, vgd_hxb_sp;
  Real vgd_hxb_t = 0.0, vgd_hxb_ngroup = 0.0, vgd_hxb_nrealloc = 0.0;
  void VetGdHxBatch(const std::vector<std::array<int, 4>> &jobs);
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
  void VetGdIpcInit(const std::vector<int> &prk, size_t cs, size_t half);
  void VetGdIpcCheck();
  // pending list-path exchange of the overlap (VetGdHlBegin / VetGdHlEnd)
  struct VgdHlPend {
    bool on = false, xe = false;                // xe: X is an exact (vgd_hxe) entry
    bool ipc = false;                           // the IPC path
    // IPC: begin + end as one kernel (pack, flags, wait, scatter). now: the end follows
    // the begin; done: the scatter ran; cnt: teams done packing
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
  template <class V>
  bool VetGdHlBegin(V &a, const int nv, const int i0, const int ws);
  void VetGdHlEnd();
  std::vector<int> vgd_wsh;    // (i): band depth the shell's data needs (<= vgd_w)
  std::vector<int> vgd_wsi, vgd_wso;   // (i): same per pass (inward, outward)
  template <class V>
  void VetGdHalo(V &a, const int nv, const int i0, const int i1, const int ws,
                 const bool mapd);
  // multi-shell (ln chi, ln S) halo in chunks of shells that fit the per-shell buffers
  void VetGdHaloCs(DvceArray5D<Real> &a, const int i0, const int i1);
  Real vgd_tsrc = 0.0, vgd_tswp = 0.0, vgd_tmom = 0.0;
  void VetGdInit();
  Real VetGdAngle(const int cyc) const;
  bool VgdRotOn() const {return vgd_rot > 0 || vgd_rot_dt > 0.0 || vgd_rot_auto;}
  int64_t VetGdBlock(const int cyc) const;       // rotation block
  Real VetGdAngleBlk(const int64_t b) const;
  double RotPhase(const int64_t b) const;   // block -> fraction of period
  Real VetGdAlignDt(const Real time, const Real dt);
  void VetGdTables(const Real alpha);
  // host part of the tables in persistent pinned buffers; the next rotation's tables
  // prefetched (host thread + upload into a second device set on vgd_tex)
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
    int nthx = 0;                                     // nthx: theta-only wall maps
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
  // shell body closure in device memory, called through a pointer (VgdClosureToDevice)
  DvceArray1D<char> vgd_bbuf0, vgd_bbuf1;
#if defined(KOKKOS_ENABLE_CUDA)
  DvceArray1D<char> vgd_bbuf2, vgd_bbuf3;   // CUDA gx-false closures
#endif
  int vgd_scut = 0;            // vlat_icut the running sweep uses (sync: = vlat_icut)
  DevExeSpace vgd_cur;         // instance the sweep/halo use
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
