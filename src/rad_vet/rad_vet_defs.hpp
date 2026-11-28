#ifndef RAD_VET_RAD_VET_DEFS_HPP_
#define RAD_VET_RAD_VET_DEFS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_defs.hpp
//! \brief array slots, key values and small structs of <rad_vet>.

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

namespace radvet {
struct SolveArgs;   // rad_vet_args.hpp
struct MgLevel;  // rad_vet_args.hpp

//! \struct WorkVec
//! \brief the clean work vectors cln_w (CLW_*): slots n < ns in s, the layout's own
//! slots in o (n - ns; staggered: in mem_arena, MemArenaSeal)

struct WorkVec {
  DvceArray5D<Real> s, o;
  int ns = 1 << 30;
  WorkVec() = default;
  WorkVec(const DvceArray5D<Real> &a) : s(a) {}   // NOLINT (one array)
  KOKKOS_FORCEINLINE_FUNCTION
  Real &operator()(const int m, const int n, const int k, const int j,
                   const int i) const {
    return (n < ns) ? s(m,n,k,j,i) : o(m,n-ns,k,j,i);
  }
};

//! \struct FaceVars
//! \brief the face variables of one direction (CLF_*): slots n < CLF_NS in s, the others
//! (dead across a closure build) in o (n - CLF_NS; over the intensities, MemArenaBind)

constexpr int CLF_NS = 11;

struct FaceVars {
  DvceArray5D<Real> s, o;
  KOKKOS_FORCEINLINE_FUNCTION
  Real &operator()(const int m, const int n, const int k, const int j,
                   const int i) const {
    return (n < CLF_NS) ? s(m,n,k,j,i) : o(m,n-CLF_NS,k,j,i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  int extent_int(const int r) const {return s.extent_int(r);}   // (not of the slots)
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
                             // x1, x2, x3 (+2), AreaStream
constexpr int CLG_Q = 38;    // staggered: face family and orientation in the owner's axes
                             // (CsSlotQ), 4 bits per x2 (low) and x3 face
constexpr int CLG_NV = 39;
constexpr int CLG_NA = 25;   // components of CsGeoArrays::a

//! \struct CsGeoArrays
//! \brief cln_gn (CLG_*) in two arrays: components constant along x1 (CLG_N, CLG_R,
//! CLG_E, CLG_F, CLG_F + 4..6) once per (m,k,j) and x1 range (low ghosts, active, high
//! ghosts) in a; the others (CLG_X, CLG_M, CLG_F + 3, CLG_F + 7, CLG_A) in r

struct CsGeoArrays {
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
constexpr int VET_E  = 0;
constexpr int VET_F1 = 1;
constexpr int VET_F2 = 2;
constexpr int VET_F3 = 3;
constexpr int VET_NVAR = 4;

// explicit stages of the PD-ARS (Chu et al. 2019) integrator per substep
constexpr int VET_NSTAGE = 2;

//! \fn VetAnalytic
//! \brief the input names an analytic closure_form (eddington | sphere): no rays
inline bool VetAnalytic(ParameterInput *pin) {
  if (!pin->DoesParameterExist("rad_vet", "closure_form")) {return false;}
  const std::string c = pin->GetString("rad_vet", "closure_form");
  return c == "eddington" || c == "sphere";
}

//! \fn VetCurvedRays
//! \brief spherical polar, cubed sphere: the rays of closure = sc.  Off with an analytic
//! closure_form; an input that names the old key vet_gd decides
inline bool VetCurvedRays(ParameterInput *pin) {
  return pin->DoesParameterExist("rad_vet", "vet_gd") ?
         pin->GetBoolean("rad_vet", "vet_gd") : !VetAnalytic(pin);
}

//! \fn VetKeyName
//! \brief the <rad_vet> name to read for a key with aliases: canon if present, else the
//! first alias present, else canon. The canonical key wins over a differing alias.
inline std::string VetKeyName(ParameterInput *pin, const std::string &canon,
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

//! \fn VetLzKey
//! \brief key of the clean L_z ledger (and K1's sum) on every mesh; old files' cs / sp
//! spellings are read as aliases
inline std::string VetLzKey(ParameterInput *pin, const bool cs, const bool k1) {
  return VetKeyName(pin, k1 ? "lz_ledger_k1" : "lz_ledger",
                   {cs ? (k1 ? "lz_cs_k1" : "lz_cs_ledger")
                       : (k1 ? "lz_sp_bnd_k1" : "lz_sp_bnd")});
}

// components of RadiationVet::vet_cell (closure = sc)
constexpr int VETRAY_CHX = 0;   // extinction rho (kappa_F + kappa_s) [1/length]
constexpr int VETRAY_SRC = 1;   // source function in E units (4 pi S / c)
constexpr int VETRAY_J   = 2;   // E of the formal solution (4 pi J / c)
constexpr int VETRAY_K11 = 3;   // K_ab in E units: 11 22 33 12 13 23
constexpr int VETRAY_H1  = 9;   // F_a / c of the formal solution: 1 2 3
constexpr int VETRAY_CHI = 12;  // the uniaxial projection: chi, n1, n2, n3
constexpr int VETRAY_N1  = 13;
constexpr int VETRAY_D11 = 16;  // vet_tensor = full: guarded D = K/J, 11 22 33 12 13 23
constexpr int VETRAY_GD  = 22;  // full: |D_guarded - D_raw| (max norm) of the cell
constexpr int VETRAY_HP1 = 23;  // sum w I max(mu_a, 0), a = 1 2 3 (half-range)
constexpr int VETRAY_NL  = 26;  // sum w |I - S| (C1 non-locality), fs only
constexpr int VETRAY_EPS = 27;  // thermal fraction min(kappa_P/chi, 1) (fs)
constexpr int VETRAY_NC  = 28;
// closure_consistency (ScSlabFaceSum, the slots of VetScFaceSum): ray fluxes sum w g
// Phi, sum w g Phi Omega_b (g = |mu_1| dx2/dx1) through the low x2 face; sum w mu_1 I,
// sum w mu_1 I Omega_b through the x1 boundary face and the low x1 face; x3 as x2
constexpr int VETRAY_F2  = 37;
constexpr int VETRAY_F1  = 41;
constexpr int VETRAY_FX  = 45;
constexpr int VETRAY_F3  = 49;
constexpr int VETRAY_NFS = 16;  // the face sums F2, F1, FX, F3
// extra vgd_hr components of the vet_gd thin switch (ScHalfRange)
constexpr int VET_VGD_HRC1 = 6;     // C1 = sum w |I - S| / J
constexpr int VET_VGD_HREP = 7;     // eps = min(kappa_P/chi, 1) (1 where chi = 0)
constexpr int VET_VGD_HRJ  = 8;     // J (E units)
constexpr int VET_VGD_HRNX = 9;
#ifdef VET_VGD_NMAX_BUILD
constexpr int VET_VGD_NMAX = VET_VGD_NMAX_BUILD;   // cmake -D Athena_M1_VGD_NMAX
#else
constexpr int VET_VGD_NMAX = 432;   // vet_gd: max directions
#endif
// components of opac, all stored as rho*kappa (inverse length)
constexpr int VET_OP_P = 0;  // rho kappa_P, Planck (emission) mean
constexpr int VET_OP_E = 1;   // rho kappa_E, energy (absorption) mean
constexpr int VET_OP_T = 2;   // rho (kappa_F + kappa_s): transport opacity (flux source,
                             // tau_face)
constexpr int VET_NOPAC = 3;
// vet_scatter only: rho kappa_e, the electron-scattering part of VET_OP_T (which stays
// the total rho kappa_R)
constexpr int VET_OP_S = 3;

// <rad_vet>/opacity
constexpr int VET_OPAC_CONST    = 0;
constexpr int VET_OPAC_POWERLAW = 1;
constexpr int VET_OPAC_USER     = 2;
constexpr int VET_OPAC_TABLE    = 3;   // stellar Rosseland + Planck tables

// <rad_vet>/force_reference: what the momentum coupling subtracts from the radiative
// force (arad_ref)
constexpr int VET_FREF_NONE    = 0;
constexpr int VET_FREF_WB_ARAD = 1;
constexpr int VET_FREF_CLEAN   = 2;   // scheme = clean: its own force of the first solve
// wb_arad_md: staggered clean, steady-start force in
constexpr int VET_FREF_WB_MD   = 3;
                                     // every direction (cln_fref)

//----------------------------------------------------------------------------------------
//! \struct VetOpacTab
//! \brief the stellar opacity tables as device Views held by value (capturable in a
//! device lambda); one (log10 T, log10 rho) grid, read with RosselandTable.
//! Lookup in K and g/cm^3: T[K] = T_code*tunit, rho = rho_code*dunit; kunit converts
//! cm^2/g to code units.

struct VetOpacTab {
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
  // scat 0 none, 1 constant kes, 2 kxe x_e/mu (VetScatterEos); code units (times kunit)
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
//!   the ma form (Jiang 2021 eq. 6): S = max(J + (rho k_P/chi)(B - J), 0)
KOKKOS_INLINE_FUNCTION
Real VetScatterMix(const Real rkp, const Real chi, const Real b, const Real e) {
  const Real ep = rkp/chi;
  return fmax(e + ep*(b - e), 0.0);
}

//! \fn VetScatterGuard
//! \brief guarded kappa_e (per mass) for total kt: kappa_a = kt - kappa_e >= kfl
KOKKOS_INLINE_FUNCTION
Real VetScatterGuard(const Real kt, const Real ke, const Real kfl) {
  return (kt < ke + kfl) ? fmax(kt - kfl, 0.0) : ke;
}

// entries of the implicit-solve diagnostic counter (see RadiationVet::cnt)
constexpr int VET_CNT_NSOLVE = 0;   // cells handed to the bracketed Newton
constexpr int VET_CNT_ITSUM  = 1;   // total iterations over those cells
constexpr int VET_CNT_ITMAX  = 2;   // largest iteration count of any one cell
constexpr int VET_CNT_NFAIL  = 3;   // cells that hit maxit without reaching the tolerance
constexpr int VET_CNT_NBRAK  = 4;   // cells for which no bracket [T_lo, T_hi] was found
constexpr int VET_NCNT = 5;

//----------------------------------------------------------------------------------------
//! \struct RadiationVetTaskIDs
//! \brief container to hold TaskIDs of all rad_vet tasks

//! \struct RingMates
//! \brief the MeshBlocks of each local block's line along a direction (RingMatesInit)
struct RingMates {
  int nmb = -1, nt = -1, nslot = 0;   // blocks here and in the mesh; received blocks
  DvceArray2D<int> src;               // (m, position): local block or -1 - received slot
  std::vector<std::array<int, 3>> snd, rcv;   // (block or slot, rank, source gid)
};


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
//! \struct SolvePass
//! \brief state of the passes of one stage solve (StageSolve); host only

struct SolvePass {
  int n = 0, nmax = 1;         // pass, cap of the passes
  bool iterate = false;        // passes (gas coupling, lagged closure), else one solve
  bool clos1 = false;          // the stage-1 solve takes its own closure once (pass 1)
  bool t2stage = false;        // a stage solve (its failure: the backward-Euler redo)
  bool converged = false;
  bool clos = false, opr = false;   // this pass takes a new closure, new opacities
  bool lp = false;             // implicit_lin_maxit: the linear solve ended at its cap
  int nit = 0, nit_sum = 0, nsolve = 0, nneg = 0;
  std::array<int, 2> nfirst = {0, 0};
  Real tol_nl = 0.0, trel = 0.0;    // tolerance of the changes
  Real de = 0.0, dv = 0.0;     // the pass's changes of E and v, the ones before
  Real de0 = 1.0e300, dv0 = 1.0e300;
  Real dnl = 0.0, dnl0 = -1.0; // the last pass's change and the one before
  std::vector<Real> rhist;     // changes of the passes (clean_gas_newton_switch)
  // a change that stopped decreasing within 100 x both its round-off level (cln_fe,
  // cln_fv) and tol_nl is converged too
  bool Settled(const Real d, const Real d0, const Real f) const {
    return d < tol_nl || (d > 0.5*d0 && d < 100.0*fmin(f, tol_nl));
  }
};

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_DEFS_HPP_
