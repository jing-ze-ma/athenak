#ifndef UTILS_TWO_STREAM_RT_HPP_
#define UTILS_TWO_STREAM_RT_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_rt.hpp
//! \brief The two-stream radiative-transfer SOLVER: the picket-fence and correlated-k
//! column sweeps, their configuration, their scratch arrays, and the explicit-source
//! limiter.
//!
//! EXTRACTED from src/pgen/deep_hot_jupiter_rt.cpp, where it grew, so that a second
//! problem generator can use it.  The arithmetic is untouched; the extraction is gated on
//! a bitwise-identical problem/ck_dump_file column.  What changed: the enclosing
//! namespace and `inline` on the namespace-scope definitions.
//!
//! HOW IT IS CONFIGURED.  Through the namespace-scope variables below (rt_ck, rt_split,
//! rt_nchain, rt_ck_pcut, rt_de_max, rt_int_at_cut, the dump controls) and through
//! `pm->pgen->hot_jupiter_param`, which despite its name is a generic carrier for the
//! star and planet numbers this solver needs: Teq, grav, ap, omega, Rgas, met and the
//! gravity flags.  A self-luminous object sets Teq = 0.  Tidying that into an explicit
//! config struct is the obvious next step; it was deliberately NOT done in the same
//! change as the move, so that the move could be proved byte-for-byte.
//!
//! WHAT IT NEEDS.  correlated_k for the band opacities and atm_column for the geometry
//! and gravity; both are included here.  The scratch arrays are the rt_*_ptr pointers,
//! allocated by the caller once the mesh exists.

#include <math.h>
#include <cstdio>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "coordinates/cell_locations.hpp"
#include "utils/eint_from_cons.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"
#include "utils/correlated_k.hpp"
#include "utils/atm_column.hpp"
#include "pgen/roche_common.hpp"
#include "pgen/pgen_eos_utils.hpp"
#include "utils/rad_taper.hpp"
#include "utils/two_stream_warm_rst.hpp"
#include "utils/two_stream_column_implicit.hpp"
#include "utils/two_stream_column_partition.hpp"
#include "utils/two_stream_column_ck.hpp"
#include "utils/two_stream_ck_rst_state.hpp"

namespace two_stream_rt {

//! \fn void par_reduce_clip4 / par_reduce_clip3
//! \brief par_for with an int sum reduction bolted on (clip4: j fastest, see below;
//! clip3: flattened the way athena.hpp's par_for flattens its ranges). These exist only
//! so the source limiter can report how many cells it clipped without a second pass over
//! the grid; that is why they live here rather than in athena.hpp.

template <typename Function>
inline void par_reduce_clip4(const std::string &name, const int ml, const int mu,
                             const int kl, const int ku, const int jl, const int ju,
                             const int il, const int iu, int &nclip,
                             const Function &function, const bool async = false) {
  const int nk = ku-kl+1, nj = ju-jl+1, ni = iu-il+1;
  const int nkji = nk*nj*ni, nji = nj*ni;
  // ck-jlin: j FASTEST, not i.  The only caller (rt_apply) reads ~5 nblk-slot sums per
  // cell from the chain-block arrays laid out (m, blk, i, k, j), against a dozen
  // (m, k, j, i) reads; with i fastest every lane of a warp hit its own cache line.
  // Each cell's arithmetic is untouched, and the int count is order-free.
  auto body = KOKKOS_LAMBDA(const int &idx, int &sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int i = (idx - m*nkji - k*nji)/nj;
    int j = (idx - m*nkji - k*nji - i*nj) + jl;
    function(m+ml, k+kl, j, i+il, sum);
  };
  if (async) {
    // problem/ck_impl_nosync: the count goes to a device scalar that nobody reads, so
    // the launch does not block the host (the caller knows the count is 0)
    static Kokkos::View<int, DevMemSpace> *dcnt = nullptr;
    if (dcnt == nullptr) dcnt = new Kokkos::View<int, DevMemSpace>("rt_nclip_dev");
    Kokkos::parallel_reduce(name,
    Kokkos::RangePolicy<>(DevExeSpace(), 0, (mu-ml+1)*nkji), body, *dcnt);
    return;
  }
  int cnt = 0;
  Kokkos::parallel_reduce(name,
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (mu-ml+1)*nkji), body, cnt);
  nclip += cnt;
}

template <typename Function>
inline void par_reduce_clip3(const std::string &name, const int ml, const int mu,
                             const int kl, const int ku, const int jl, const int ju,
                             int &nclip, const Function &function) {
  const int nk = ku-kl+1, nj = ju-jl+1;
  const int nkj = nk*nj;
  int cnt = 0;
  Kokkos::parallel_reduce(name,
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (mu-ml+1)*nkj),
  KOKKOS_LAMBDA(const int &idx, int &sum) {
    int m = idx/nkj;
    int k = (idx - m*nkj)/nj;
    int j = (idx - m*nkj - k*nj);
    function(m+ml, k+kl, j+jl, sum);
  }, cnt);
  nclip += cnt;
}

using pgen_eos::EintFromP;
using pgen_eos::PresFromEint;
using pgen_eos::TempKelvin;
using pgen_eos::PresTempFromEint;
using pgen_eos::DensFromPT;
using pgen_eos::GradAd;
using pgen_eos::EintFromDensT;

using correlated_k::CK_NB;
using correlated_k::CK_NG;
using correlated_k::CK_DIFFUSIVITY;
using correlated_k::ck_nq;
using correlated_k::ck_nT;
using correlated_k::ck_nP;
using correlated_k::ck_lT_ptr;
using correlated_k::ck_lP_ptr;
using correlated_k::ck_lk_ptr;
using correlated_k::ck_gw_ptr;
using correlated_k::ck_wl_ptr;
using correlated_k::ck_swf_ptr;
using correlated_k::ck_pf_ptr;
using correlated_k::ck_pf_lTmin;
using correlated_k::ck_pf_idlT;
using correlated_k::ce_nT;
using correlated_k::ce_nP;
using correlated_k::ce_lT_ptr;
using correlated_k::ce_lP_ptr;
using correlated_k::ce_ptr;
using correlated_k::cia_nT_ptr;
using correlated_k::cia_T_ptr;
using correlated_k::cia_k_ptr;
using correlated_k::ray_x_ptr;
using correlated_k::ck_tp_index;
using correlated_k::ck_kappa;
using correlated_k::ck_continuum;
using correlated_k::ck_planck_frac;
using correlated_k::ck_planck_bands;
using atm_column::TGuess;
using atm_column::CSCellAngles;
using atm_column::GravAccAt;
using atm_column::GravPotAt;
using atm_column::TideAccR;
using atm_column::TideAccT;
using atm_column::TideAccP;
using atm_column::EffGravAt;

// Number of (band, quadrature) chains stepped together in the IR sweep. This is the
// instruction-level parallelism the sweep gets; it also fixes the private I_ir_down_c
// footprint at NC radial columns, independent of how many chains are requested.
#ifndef RT_NB
#define RT_NB 4
#endif

// Radial extent of the private intensity column in the PICKET-FENCE split kernel
// (rt_chain, rt_ck = false, rt_split = true), a test path that keeps the fixed size,
// guarded.  The correlated-k kernels (rt_chain_ck, rt_chain_ck_lin, rt_chain_ck_lin1) and
// the grey chain (rt_chain_grey) no longer use a compile-time column: their whole-column
// arrays live in the n1-sized global buffer rt_ckscr_ptr (ck-scratch, ck-tiers), so any
// n1 runs with no tier waste and no cap.
#ifndef RT_NNC
#define RT_NNC 72
#endif

// Cache the per-layer two-stream coefficients between the two sweeps of the split
// chain kernel. The down-sweep at layer i-1 and the up-sweep at layer i are the SAME
// layer -- both use dtau = tau[i-1]-tau[i] and the same (T,p) proxy -- so they compute
// the identical e0, and alp(down) == gm(up), bet(down) == bet(up). Caching removes one
// expm1, one divide and (with the table on) one k-table lookup and two logs per cell
// per chain, at the price of three more private [NC][NNC] arrays.
#ifndef RT_CACHE
#define RT_CACHE 0
#endif

// Run the correlated-k chain kernel's recurrence in single precision. Default OFF, but
// the
// reason has changed and is worth stating precisely.
//
// Measured while the kernel was still starved on scattered loads, this bought 2.8 %.
// Re-measured once the array layout was fixed it buys 1.42x on rt_chain, 457 -> 325 ms
// per
// 100 cycles: with the memory fed, the FP64 transcendentals really are the next thing in
// the way. That is a 25 % saving on RT.
//
// It is off anyway because 25 % of RT is only 6 % of the run. RT is 41 % of kernel time
// and kernel time is less than wall, so the whole trade is 8.62 -> 8.10 s of integration
// for a 1.7e-4 relative change in the net longwave flux. Judge it on the total, not on
// the
// kernel.
//
// Accuracy if it is ever wanted: the recurrence is a contraction and single precision
// handles it fine. The exposure is the NET flux, which deep down is the difference of
// I_up and I_down when both are close to B. On a real column the outgoing flux at the top
// agrees to 1e-7 and the flux at the cut to 1e-5; the worst level is 1.7e-4.
#ifndef RT_FP32
#define RT_FP32 0
#endif
#if RT_FP32
using RtF = float;
#define RT_EXPM1(x) expm1f(x)
#define RT_EXP(x)   expf(x)
#else
using RtF = Real;
#define RT_EXPM1(x) expm1(x)
#define RT_EXP(x)   exp(x)
#endif

// THE ck CHAIN KERNEL'S WHOLE-COLUMN ARRAYS LIVE IN GLOBAL MEMORY (ck-scratch).
// rt_chain_ck keeps up to ten [RT_NB][n1] column arrays per thread (I_down, Cmx, Krs,
// Kpc, Cc0/Cci/Cco, Js0..2).  As private arrays they were sized by a compile-time column
// tier and put in the private (scratch) segment: 31 kB per lane for the production
// kernel at the 136 tier, 60 kB at 264 (86 kB with the Jacobian), which the MI300A
// runtime reserves for FULL occupancy (~40 GB at 264) and cannot; every ck run with
// n1 > 136 died with HSA_STATUS_ERROR_OUT_OF_RESOURCES.  They are now slices of one
// device buffer (rt_ckscr_ptr) sized by the actual n1 and thread count, laid out the way
// the hardware lays out scratch: threads in tiles of CKS_W = 64 (one wavefront), each
// tile owning one contiguous region, and inside it group g (one array, all RT_NB chains)
// holding element (cc, i) of lane l at
//   tile + ((g*RT_NB + cc)*n1 + i)*CKS_W + l.
// So a wave's access at one (cc, i) is one contiguous 512-byte line, and the whole of a
// wave's column data is one contiguous block (few pages; a plain (i, thread) layout put
// successive i half a megabyte apart and cost 30 % on the kernel in TLB misses).  The
// thread number is the par_for flattening of (m, blk, k, j), so a tile is a wave.  Only
// the storage moves: every value, and every operation on it, is the one the private
// array held.  A switched-off array (the SPH/BSP/CCH/JAC tag false) stays the
// one-element private array it was.  CkScrCol[cc][i] reads like the array.
constexpr int CKS_W = 64;
template <typename T>
struct CkScrRow {
  T *p;
  KOKKOS_INLINE_FUNCTION T &operator[](const int i) const { return p[i*CKS_W]; }
};
template <typename T>
struct CkScrCol {
  T *p;
  int sc;
  KOKKOS_INLINE_FUNCTION CkScrRow<T> operator[](const int c) const {
    return CkScrRow<T>{p + c*sc};
  }
};
template <typename T>
struct CkScrOne {
  T a[1][1];
  KOKKOS_INLINE_FUNCTION T *operator[](const int c) { return a[c]; }
  KOKKOS_INLINE_FUNCTION const T *operator[](const int c) const { return a[c]; }
};
// group g of thread t; tsz = Reals per tile.  ON false gives the one-element stand-in.
template <typename T, bool ON>
KOKKOS_INLINE_FUNCTION
std::conditional_t<ON, CkScrCol<T>, CkScrOne<T>> CkScrGet(Real *buf, const int g,
                                                         const size_t tsz, const int n1,
                                                         const int t) {
  if constexpr (ON) {
    Real *tile = buf + static_cast<size_t>(t/CKS_W)*tsz;
    T *p = reinterpret_cast<T*>(tile + static_cast<size_t>(g)*RT_NB*n1*CKS_W)
           + (t % CKS_W);
    return CkScrCol<T>{p, n1*CKS_W};
  } else {
    return CkScrOne<T>{};
  }
}
// ck-store: the chain kernel's groups INTERLEAVED by radial index instead: inside a
// tile, element (cc, i) of group g of lane l sits at
//   tile + ((i*S + g)*RT_NB + cc)*CKS_W + l   (in Reals),  S = the kernel's group count,
// so at a given i every (group, chain) of the thread is a compile-time offset from ONE
// address.  The group-major layout above needs a separate pointer per (group, chain),
// runtime multiples of n1, which the chain kernel kept live in registers (about 40 of
// its 164).  A wave's access at one (g, cc, i) is still one contiguous line, and the
// tile is the same size.  An element takes one Real slot whatever T is.
template <typename T, int S>
struct CkScrRowI {
  Real *p;
  KOKKOS_INLINE_FUNCTION T &operator[](const int i) const {
    return *reinterpret_cast<T*>(p + i*(S*RT_NB*CKS_W));
  }
};
template <typename T, int S>
struct CkScrColI {
  Real *p;
  KOKKOS_INLINE_FUNCTION CkScrRowI<T, S> operator[](const int c) const {
    return CkScrRowI<T, S>{p + c*CKS_W};
  }
};
template <typename T, bool ON, int S>
KOKKOS_INLINE_FUNCTION
std::conditional_t<ON, CkScrColI<T, S>, CkScrOne<T>> CkScrGetI(Real *buf, const int g,
                                                              const size_t tsz,
                                                              const int t) {
  if constexpr (ON) {
    Real *tile = buf + static_cast<size_t>(t/CKS_W)*tsz;
    return CkScrColI<T, S>{tile + g*RT_NB*CKS_W + (t % CKS_W)};
  } else {
    return CkScrOne<T>{};
  }
}
// ck-tiers: the same wave-tiled layout for the other column kernels (rt_chain_ck_lin,
// rt_chain_ck_lin1, rt_chain_grey).  Thread t's tile holds rows of n1 values; row r
// (one chain or angle of one array) holds element i of lane l at
//   tile + (r*n1 + i)*CKS_W + l,  tsz = (rows per thread)*n1*CKS_W.
// CkScrRows(...)[c][i] is row r0 + c.
template <typename T>
KOKKOS_INLINE_FUNCTION
CkScrCol<T> CkScrRows(Real *buf, const int r0, const size_t tsz, const int n1,
                      const int t) {
  Real *tile = buf + static_cast<size_t>(t/CKS_W)*tsz;
  T *p = reinterpret_cast<T*>(tile + static_cast<size_t>(r0)*n1*CKS_W) + (t % CKS_W);
  return CkScrCol<T>{p, n1*CKS_W};
}
// the buffer (grown on demand, never shrunk; leaked like the other rt_ Views) and the cap
// on its size: a launch whose groups would exceed rt_ckscr_gb is split over chain blocks
inline DvceArray1D<Real> *rt_ckscr_ptr = nullptr;
inline size_t rt_ckscr_len = 0;
inline Real rt_ckscr_gb = 16.0;
inline void CkScrEnsure(const size_t len) {
  if (rt_ckscr_ptr != nullptr && rt_ckscr_len >= len) return;
  if (rt_ckscr_ptr != nullptr) delete rt_ckscr_ptr;
  rt_ckscr_ptr = new DvceArray1D<Real>(
      Kokkos::ViewAllocateWithoutInitializing("rt_ckscr"), len);
  rt_ckscr_len = len;
}
// ck-determinism: the JAC pass of rt_chain_ck (ck_implicit without ck_impl_jac_lin)
// added every chain block's share of the tridiagonal into ck_jac with atomics, whose
// order -- and so whose rounding -- changed from run to run on a GPU.  Each block now
// adds into its OWN slots (m, 3*blk + s, k, j, i) of this View with plain adds (the one
// thread of (m, blk, k, j) owns them), and ck_jacq_sum adds the blocks into ck_jac in
// block order: the same sums in a fixed order, so the run is reproducible.
inline DvceArray5D<Real> *ck_jacq_ptr = nullptr;
#define SQRT3_VEF 1.7320508075688772
// problem/ck_sph_face = 5 (see ck_sph_face): s = sqrt(3 f) at the frames and the face
// factors gamma (+ the top reflection zeta), per (band, g); the formal solution's work
inline DvceArray5D<Real> *ck_vef_s_ptr = nullptr;
inline DvceArray5D<Real> *ck_vef_g_ptr = nullptr;
inline DvceArray5D<Real> *ck_vef_w_ptr = nullptr;
//! \fn bool CkVefDue
//! \brief ck_sph_face = 5: is this pass the one that refreshes the VEF factors?  With
//! ck_vef_sync only a STORING pass, the first one or one ck_vef_every or more cycles
//! after the last refresh.  The same test decides the restart snapshot (CkVefSnap).
inline bool CkVefDue(Mesh *pm, const bool ckfop, const bool ckfst);
KOKKOS_INLINE_FUNCTION void CkJacQAdd(Real *p, const Real v) { *p += v; }

// Column solves ("chains") the RT kernel steps per cell: 4 for the grey picket fence
// (two IR channels x the two Gauss angles), CK_NB*CK_NG*ck_nquad for correlated-k.
// Derived from the scheme, not an input.
inline int rt_nchain = 4;

// problem/rt_split: run the RT as three kernels with the chain-block index promoted to
// a parallel dimension, instead of one kernel that loops over chains serially. Same
// arithmetic, same block-summation order -- see the use site.
// DEFAULT TRUE in deep_hot_jupiter_rt since 2026-09-22; box_convection and red_giant
// force it true (the grey kernel lives on the split path).  The monolithic path is
// therefore reachable only as the grey PICKET FENCE (rt_ck = false) with the flag set
// false by hand, and is kept for that case: the split grey sweep sizes its column at
// compile time and fatals when n1 > RT_NNC, which the monolithic one does not.
//
// THE "BITWISE IDENTICAL" CLAIM BELOW HOLDS FOR THE CORRELATED-K AND GREY KERNELS, NOT
// FOR THE PICKET FENCE.  MEASURED 2026-09-22 on inputs/mhd/deep_hot_jupiter_rt_ideal_xe
// (rt_ck = false, 20 cycles, same binary, deterministic on rerun): the two paths agree
// on mass and total energy to the printed digits but differ by 7.1 % in the 1-KE
// history column, so a picket-fence run that relied on the old default answers
// differently now.  Unexplained; the difference is confined to that path.
inline bool rt_split = false;
// PER-CYCLE RT DIAGNOSTIC (deep_hot_jupiter_rt's problem/diag_gid).  Off by default and
// costing nothing when off: the array is not even allocated.  (m,slot,k,j,i) with
// slot 0 = the source RATE src [code units, erg/cm^3/s] BEFORE the limiter and before
// multiplying by bdt, slot 1 = the de actually added to u0(IEN) after the semi-implicit
// relaxation and the clip, slot 2 = 1.0 if the limiter changed de on this call,
// slot 3 = Ft, the block-summed net flux on the TOP face i+1 after the taublend factor,
// slot 4 = Fb, the same on the BOTTOM face i -- exactly the two values that make up
// src = -(Ft-Fb)/dx1 -- slot 5 = Qs, the block-summed stellar heating (0 below icut),
// slot 6 = Em, the block-summed emission rate with the semi-implicit branch's taublend
// factor, computed for the diagnostic even when the source is applied explicitly.
inline bool rt_diag = false;
inline DvceArray5D<Real> *rt_diag_ptr = nullptr;
// number of cells the RT source limiter clipped on the LAST call (host side, free)
inline int rt_nclip_last = 0;
// (m,k,j,i) face optical depth from the top
inline DvceArray4D<Real> *rt_tau_ptr = nullptr;
inline DvceArray4D<Real> *rt_B_ptr = nullptr;     // (m,k,j,i) Planck function
inline DvceArray4D<Real> *rt_Qv_ptr = nullptr;    // (m,k,j,i) stellar heating rate
inline DvceArray4D<Real> *rt_cf_ptr = nullptr;    // (m,k,j,{gamir1,gamir2,beta})
// (m,blk,k,j,i) IR flux, one slot per block
inline DvceArray5D<Real> *rt_Fb_ptr = nullptr;

// Correlated-k per-cell fields, filled by rt_pre and read by the chain kernel. Laid out
// (m,band,k,j,i) so that the radial index is contiguous for a thread sweeping one band --
// (m,k,j,i,band) would stride every read of the sweep by CK_NB.
inline DvceArray5D<Real> *rt_kc_ptr = nullptr;    // (m,b,k,j,i) continuum kappa [cm^2/g]
inline DvceArray5D<Real> *rt_Bb_ptr = nullptr;    // (m,b,k,j,i) band Planck intensity
inline DvceArray4D<Real> *rt_T_ptr = nullptr;     // (m,k,j,i) temperature [K]
inline DvceArray4D<Real> *rt_pb_ptr = nullptr;    // (m,k,j,i) pressure [bar]
// (m,k,j,i) continuous k-table index in T
inline DvceArray4D<Real> *rt_xT_ptr = nullptr;
// (m,k,j,i) continuous k-table index in p
inline DvceArray4D<Real> *rt_xP_ptr = nullptr;
// (m,k,j) deepest cell doing correlated-k
inline DvceArray3D<int>  *rt_icut_ptr = nullptr;
// (m,blk,k,j,i) stellar heating, one per block
inline DvceArray5D<Real> *rt_Qb_ptr = nullptr;
// (m,blk,i,k,j) the LOCAL thermal emission rate per unit volume [erg/cm^3/s] that each
// chain block puts into cell i.  Summed over blocks this is 4 sigma kappa_P rho T^4, the
// only part of the radiative source that depends on the cell's OWN temperature.  It is
// what turns the explicit update into a stable exponential one; see the apply kernel.
inline DvceArray5D<Real> *rt_Em_ptr = nullptr;
inline DvceArray5D<Real> *rt_Src_ptr = nullptr;   // per-cell net source, formed DIRECTLY

// --- correlated-k (Lee/Exo-FMS premixed tables, Kataria+2013 11-band grid) --------
// problem/rt_ck turns it on; problem/ck_table is the path to the premixed table and
// problem/ck_pcut_bar the pressure below which correlated-k is used at all (deeper than
// that the atmosphere is optically thick -- grey tau ~ 9e3 at 10 bar for this setup --
// and the interior is far outside any molecular table, up to 12000 K).
//
// Table layout and the traps it carries are documented in data/exo_fms_ck/PROVENANCE.md.
// The two that matter here: the data block runs the band index BACKWARDS, and kappa is
// cgs cm^2/g, which is what this code already works in.
inline bool rt_ck = false;
// problem/ck_star_teff: host effective temperature. If > 0 the stellar spectrum is a
// blackbody at this temperature, which needs no external data since only the SHAPE is
// used. If <= 0, problem/ck_swflux is read instead.
//
// Default 6000 K, the middle of the range the ultra-hot Jupiter GCM literature actually
// uses: Tan et al. (2024, MNRAS 528, 1016) grid over host T_eff of 5500, 6000 and 6500 K,
// and Parmentier et al. (2018, A&A 617, A110) model WASP-121b at 6460 K. That range also
// keeps the flux falling bluer than the grid's 0.26 um edge down to 1-3 %, which is what
// makes the 11-band Kataria structure defensible in the first place -- an A-type host
// would put 18-32 % outside it.
inline Real rt_star_teff = 6000.0;
// problem/albedo: an observed Bond albedo.  When the key is present (>= 0) it replaces
// the Parmentier+2015 fit (get_albedo) everywhere the stellar heating uses the albedo:
// the RT pass below (every (1-albedo) Fstar deposit and the picket-fence Teff) and the
// picket-fence initial condition in the pgen.  Absent (-1, default): the fit, as before.
inline Real rt_albedo_set = -1.0;
// problem/ck_dump_file: write one column's RT solution -- level pressures, temperatures,
// net longwave flux and stellar heating -- straight out of the production kernel, once,
// at
// the first RT call. This exists so the kernel ITSELF can be compared against an external
// code on the same profile, rather than a transcription of it. See
// docs/correlated_k_rt.md.
inline std::string rt_dump_file = "";
inline int rt_dump_m = 0;
inline int rt_dump_j = -1;
inline int rt_dump_k = -1;
inline bool rt_dump_done = false;
// problem/ck_dump_kap (test only, default false): append kappa(b,g) [cm^2/g, incl. the
// continuum] of every cell and the g weights to the ck_dump_file column, so an exact
// ray-traced reference can be built from the kernel's OWN opacities.  Dump path only.
inline bool rt_dump_kap = false;
inline Real rt_ck_pcut = 10.0;                    // bar
// problem/rt_grey: the GREY two-stream, built on the correlated-k machinery rather than
// on the old picket-fence path.  One band, one opacity, taken from the conduction
// module's own kappa table (<hydro|mhd>/rad_kappa_src) so that the two-stream and the
// radiative diffusion it blends with cannot disagree about the opacity, and the full
// Planck function sigma T^4/pi as the source.
//
// WHY IT EXISTS, and why it is not the old grey path.  The old one is the Parmentier
// picket fence with the Freedman fit: two IR channels whose split (gamma_1, gamma_2,
// beta) is a fit to IRRADIATED giant planets, and an analytic opacity that is 1.7x off
// the stellar table on a cool giant.  This one shares every structural fix the
// correlated-k path got -- the layer-integrated two-stream coefficients, the flux
// DIFFERENCE deposition, the tau blend, the icut handover, the semi-implicit source --
// and differs from it in exactly one way: emission and absorption use the SAME opacity,
// so kappa cancels out of the local radiative balance and the equilibrium temperature of
// a thin layer is a fixed point no matter how the opacity varies.  A correlated-k layer
// weights its absorption by the incident spectrum and its emission by its own, so a
// layer that drifts cold can absorb less than it emits and keep drifting; that is the red
// giant top-layer runaway.  Grey cannot do that.  It also cannot represent line
// blanketing, which is the reason the band solver exists -- so this is the control, and
// the diagnosis, not automatically the production choice.
inline bool rt_grey = false;

// problem/rt_plane_parallel: run the GREY sweep on a PLANE-PARALLEL Cartesian mesh.
//
// Everything else in this file assumes a RADIAL mesh: the cell centres and faces come
// from Coordinates::x1v / xx1f and the radial width from Coordinates::dx1, and all three
// are 1x1 placeholder Views on a Cartesian mesh (coordinates.cpp reallocates them only
// for spherical-polar and cubed-sphere grids), so reading them there is an out-of-bounds
// access.  With this flag the three are taken from the MeshBlock's own RegionSize
// instead -- uniform dx1, centres from CellCenterX, faces from LeftEdgeX -- which is
// exactly what a plane-parallel column needs and costs the radial path nothing: the
// substitutions are runtime branches on a `const bool` that is false in every existing
// run, and the spherical expressions are untouched.
//
// The sweep needs no other geometric change, because with this flag on the area and
// volume helpers of the radial path (AFC/ACC/VLS/VLA below, and the SPHERICAL
// DILUTION note that introduces them) collapse to A = 1 and V = dx1, which is exactly
// the plane-parallel operator -(F_top - F_bot)/dx1 this path used to carry
// everywhere.  What DOES have to be right is the
// gravity model -- a plane-parallel box has constant g, so the caller sets
// grav_point_mass = false and stellar_tide = false, and EffGravAt then returns g
// whatever radius it is handed.
//
// RESTRICTED to the grey split path (rt_grey && rt_split): the monolithic picket-fence
// kernel and the correlated-k kernel carry stellar-beam geometry (the substellar angle,
// the slant path) that has no plane-parallel meaning.  The guard is a fatal error.
inline bool rt_plane_parallel = false;

//----------------------------------------------------------------------------------------
// problem/rt_rad_force -- the RADIATIVE MOMENTUM SOURCE that goes with the EOS's
// thin-region radiation taper (<eos>/eos_rad_rho_hi, see utils/rad_taper.hpp).
//
// WHY.  With the taper on, the hydro's pressure carries only w(rho) of the LTE radiation
// pressure, so the momentum equation loses (1-w) of -grad Prad -- and in the thin layer
// that is precisely where the real force lives: MEASURED on the relaxed B-star column,
// grad Prad/(rho g) is 0.10 at tau = 1e-2, 0.15 at tau = 0.1, 0.23 at tau = 2/3 and 0.30
// at tau = 10, and at the very top it agrees with kappa_R F/(c g) to 1.8x and below tau 1
// to a few per cent.  Dropping it would make the top of the box lighter than it is.
//
// WHAT IS ADDED, per cell:
//     f = (1 - w) rho kappa_R F_net/c   +   Prad grad w,     Prad = a T^4/3,
// with F_net the face-averaged net two-stream flux (positive upward, so the force is
// along +x1, against gravity).  The FIRST term is the direct momentum deposition the
// pressure gradient no longer supplies.  The SECOND is a correction, not a force: with
// p = pgas + w(rho) Prad the hydro differences
//     -grad(w Prad) = -w grad Prad - Prad grad w,
// and the last piece is an artefact of the taper: inside the ramp dw/dln rho is of
// order one, so it is of the same order as grad Prad itself, up to ~20 % of g.  Adding
// Prad grad w back cancels it exactly, and the total radiative force is then -grad Prad
// in the deep limit (w = 1), kappa rho F/c in the thin limit (w = 0), and continuous
// through the ramp.  The x1 component of the first term is the only one the column-wise
// two-stream knows; the second is a genuine pressure gradient and is applied in all three
// directions.  The work v.f is added to the total energy alongside.
//
// The step is EXPLICIT and lives in the same place as the tau-blend handover: after the
// semi-implicit relaxation, which must not damp it.  Default OFF; box_convection requires
// the EOS taper to be on with it.
inline bool rt_rad_force = false;
//----------------------------------------------------------------------------------------
// problem/rt_force_tau_gate (default false, BITWISE OFF) -- the OPTICAL-DEPTH-GATED form
// of the momentum source above, with rt_force_tau_lo / rt_force_tau_hi.
//
// WHY.  The form above splits the force by the EOS taper weight w(rho,T): the direct
// deposit is scaled by (1-w) and the taper artefact Prad grad w is added back.  That is
// the right split only where w itself is the statement "the gas no longer carries the
// radiation".  On a star with a density inversion w is a function of the LOCAL density
// (and, under the gate, the local T), so a dense cell sitting at tau < 1 gets w = 1 and
// loses the whole kappa F/c deposit, while a thin cell deep inside the star gets w = 0
// and is given a deposit it must not have: the force then has a dipole across the
// inversion that no optical depth justifies.
//
// WHAT THE GATE DOES.  s(tau) = 1 for tau <= rt_force_tau_lo, 0 for tau >= tau_hi and
// a cubic smoothstep in log10(tau) between, with tau the grey optical depth to the
// TOP of the column.  Where s = 1 (the thin layers) the force is the honest optically
// thin one,
//     f = rho kappa F_net/c + grad(w Prad),
// i.e. the FULL radiative deposit plus the gradient of whatever radiation pressure the
// hydro still carries -- no (1-w) prefactor, because at tau < 1 the gas does not carry
// the radiation whatever its density says.  Where s = 0 (tau above the ramp) the old
// expressions are executed verbatim, so with the gate off, or in any cell with s = 0,
// the arithmetic is bitwise the old one.  Default OFF.
inline bool rt_force_tau_gate = false;
inline Real rt_force_tau_lo = 0.3;
inline Real rt_force_tau_hi = 3.0;
// problem/rt_budget_verbose (box_convection.cpp): when non-null, the radiative force's
// WORK term v.f is accumulated, box-integrated over the step, into slot 10 of this
// array.  Diagnostic only: nothing here changes a source term.
inline DvceArray1D<Real> *rt_bud_ptr = nullptr;
// problem/rt_force_verbose -- print the hydrostatic balance of every cell inside the
// taper ramp for this many RT calls, then stop.  a_p + a_g + a_f normalised by g: inside
// the ramp this is what says whether the Prad grad w correction is doing its job.
inline int rt_force_verbose = 0;
inline Real rt_force_grav = 0.0;    // |g| used only to normalise that print
// problem/rt_top_vacuum: nothing above the top of the domain.  The default top boundary
// is the UNRESOLVED HYDROSTATIC COLUMN, of optical depth kappa p / g_eff, which is the
// right model for a star whose atmosphere continues above x1max.  A local box is cut out
// of a stratification and the cells above it are not part of the problem; more to the
// point, a box with weak gravity would be handed p/g = a huge column and sealed.  With
// this flag the downward intensity entering the top face is exactly zero, which is the
// standard grey-atmosphere boundary and the one solar_convection.cpp's local solver uses.
inline bool rt_top_vacuum = false;

// --- read-only access to the band solver's radial face flux, for diagnostics ---------
// rt_Fb holds the NET longwave flux on the radial faces, (m, block, i, k, j), in code
// flux units (pressure x velocity), for faces i = is..ie+1.  It is a namespace-scope
// pointer and lives for the run once the split band path has allocated it, so a problem
// generator can sum it without re-running the solver.  rt_face_flux_ready() is the guard
// every caller must use: the monolithic and non-band paths never fill it.
inline bool rt_face_flux_ready() {
  return (rt_Fb_ptr != nullptr) && (rt_ck || rt_grey);
}
inline int rt_face_nblk() {
  return (rt_Fb_ptr != nullptr) ? rt_Fb_ptr->extent_int(1) : 0;
}
inline DvceArray5D<Real> rt_face_flux() { return *rt_Fb_ptr; }

// the stellar flux sigma T_irr^4 and the albedo of the LAST RT call, host scalars
// recorded for the problem generator's flux history (problem/flux_hst); nothing reads
// them back into the solver.
inline Real rt_hist_fstar = 0.0;
inline Real rt_hist_albedo = 0.0;

// --- problem/flux_hst_rst: the face-flux cache across a restart (rthst-1101) ----------
// rt_Fb, rt_Qb, rt_icut and rt_cf slot 3 are written by the full RT calls only; under
// the cadence (ck_impl_every) the history rows in between read what the last full call
// left.  A restart starts with no cache, so until the next full call the flux history
// had no radiative part (Lir_top, Lsw_refl, Lsw_abs, Lrad_bot, Lsw_in = 0, Etot_top/bot
// hydro only).  The cache is written to the implicit-ck restart block as three slabs
// (two_stream_ck_rst.hpp: the chain blocks SUMMED in block order, i.e. exactly the sums
// the history forms, and per column icut, mu0, fstar, albedo).  After the restart the
// history reads the staged slabs (RtFhStageViews) until the RT arrays are allocated;
// RtFhRstAfterAlloc then puts the sums into chain block 0 of the new arrays (blocks
// 1.. stay 0, so every block sum is unchanged) and frees the staging.  Diagnostic
// only: no kernel of the solver reads these arrays across calls.
inline bool RtFhLive() {
  return rt_fhst_rst && rt_face_flux_ready() && (rt_Qb_ptr != nullptr) &&
         (rt_cf_ptr != nullptr) && (rt_icut_ptr != nullptr);
}

// one slab (nmb, n3, n2, n1) of the cache, into the host restart array dst
template <typename V>
inline void RtFhCollect(const int id, const V &dst) {
  const int nmb = static_cast<int>(dst.extent(0));
  const int n3 = static_cast<int>(dst.extent(1));
  const int n2 = static_cast<int>(dst.extent(2));
  const int n1 = static_cast<int>(dst.extent(3));
  auto fb = *rt_Fb_ptr;
  auto qb = *rt_Qb_ptr;
  auto cf = *rt_cf_ptr;
  auto icut = *rt_icut_ptr;
  const int nblk = fb.extent_int(1);
  const bool fit = (fb.extent_int(0) >= nmb) && (fb.extent_int(2) == n1) &&
                   (fb.extent_int(3) == n3) && (fb.extent_int(4) == n2);
  if (!fit) {
    Kokkos::deep_copy(dst, 0.0);
    return;
  }
  DvceArray4D<Real> tmp("rt_fh_col", nmb, n3, n2, n1);
  if (id == kCkSlabFhCol) {
    const Real fstar = rt_hist_fstar, alb = rt_hist_albedo;
    par_for("rt_fh_col", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      tmp(m,k,j,0) = static_cast<Real>(icut(m,k,j));
      tmp(m,k,j,1) = cf(m,k,j,3);
      tmp(m,k,j,2) = fstar;
      tmp(m,k,j,3) = alb;
    });
  } else {
    const bool isq = (id == kCkSlabFhQb);
    par_for("rt_fh_sum", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real f = 0.0;
      for (int b=0; b<nblk; ++b) f += isq ? qb(m,b,i,k,j) : fb(m,b,i,k,j);
      tmp(m,k,j,i) = f;
    });
  }
  DeepCopyAcross(dst, tmp);
}

// the staged slab id on the device (nmb, n3, n2, n1); empty View when it does not fit
inline DvceArray4D<Real> RtFhStageDev(const int id) {
  using UnmanagedHost4D = Kokkos::View<Real****, Kokkos::LayoutRight,
                                       Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  const std::size_t nd = static_cast<std::size_t>(ck_rst_nmb)*ck_rst_n3*ck_rst_n2*
                         ck_rst_n1;
  if (!CkRstStaged(id) || ck_rst_stage[id].size() != nd) return DvceArray4D<Real>();
  DvceArray4D<Real> d("rt_fh_stg", ck_rst_nmb, ck_rst_n3, ck_rst_n2, ck_rst_n1);
  UnmanagedHost4D h(ck_rst_stage[id].data(), ck_rst_nmb, ck_rst_n3, ck_rst_n2,
                    ck_rst_n1);
  DeepCopyAcross(d, h);
  return d;
}

inline bool RtFhStaged(const int nmb) {
  return rt_fhst_rst && (ck_rst_nmb == nmb) && (ck_rst_n1 >= 4) &&
         CkRstStaged(kCkSlabFhFb) && CkRstStaged(kCkSlabFhQb) &&
         CkRstStaged(kCkSlabFhCol);
}

// for the flux history before the first RT call after a restart: cache-shaped Views
// (one chain block) built from the staging, which stays staged.  False: nothing staged.
inline bool RtFhStageViews(const int nmb, DvceArray5D<Real> &fb, DvceArray5D<Real> &qb,
                           DvceArray4D<Real> &cf, DvceArray3D<int> &icut, Real &fstar,
                           Real &alb) {
  if (rt_Fb_ptr != nullptr || !RtFhStaged(nmb)) return false;
  auto fs = RtFhStageDev(kCkSlabFhFb);
  auto qs = RtFhStageDev(kCkSlabFhQb);
  auto cs = RtFhStageDev(kCkSlabFhCol);
  if (fs.size() == 0 || qs.size() == 0 || cs.size() == 0) return false;
  const int n3 = ck_rst_n3, n2 = ck_rst_n2, n1 = ck_rst_n1;
  fb = DvceArray5D<Real>("rt_fh_fb", nmb, 1, n1, n3, n2);
  qb = DvceArray5D<Real>("rt_fh_qb", nmb, 1, n1, n3, n2);
  cf = DvceArray4D<Real>("rt_fh_cf", nmb, n3, n2, 5);
  icut = DvceArray3D<int>("rt_fh_icut", nmb, n3, n2);
  auto fb_ = fb, qb_ = qb;
  auto cf_ = cf;
  auto icut_ = icut;
  par_for("rt_fh_stv", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    fb_(m,0,i,k,j) = fs(m,k,j,i);
    qb_(m,0,i,k,j) = qs(m,k,j,i);
    if (i == 0) {
      icut_(m,k,j) = static_cast<int>(cs(m,k,j,0));
      cf_(m,k,j,3) = cs(m,k,j,1);
    }
  });
  fstar = ck_rst_stage[kCkSlabFhCol][2];
  alb = ck_rst_stage[kCkSlabFhCol][3];
  return true;
}

// right after the split path allocated rt_Fb/rt_Qb/rt_icut/rt_cf: the staged cache into
// chain block 0, then free the staging (also when it does not fit)
inline void RtFhRstAfterAlloc(const int nmb) {
  if (RtFhStaged(nmb) && rt_Fb_ptr != nullptr && rt_Qb_ptr != nullptr &&
      rt_cf_ptr != nullptr && rt_icut_ptr != nullptr) {
    auto fb = *rt_Fb_ptr;
    auto qb = *rt_Qb_ptr;
    auto cf = *rt_cf_ptr;
    auto icut = *rt_icut_ptr;
    const int n3 = ck_rst_n3, n2 = ck_rst_n2, n1 = ck_rst_n1;
    const bool fit = (fb.extent_int(2) == n1) && (fb.extent_int(3) == n3) &&
                     (fb.extent_int(4) == n2) && (qb.extent_int(2) == n1);
    auto fs = RtFhStageDev(kCkSlabFhFb);
    auto qs = RtFhStageDev(kCkSlabFhQb);
    auto cs = RtFhStageDev(kCkSlabFhCol);
    if (fit && fs.size() > 0 && qs.size() > 0 && cs.size() > 0) {
      par_for("rt_fh_rst", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        fb(m,0,i,k,j) = fs(m,k,j,i);
        qb(m,0,i,k,j) = qs(m,k,j,i);
        if (i == 0) {
          icut(m,k,j) = static_cast<int>(cs(m,k,j,0));
          cf(m,k,j,3) = cs(m,k,j,1);
        }
      });
      Kokkos::fence();
    }
  }
  for (int s=kCkSlabFhFb; s<=kCkSlabFhCol; ++s) CkRstStageFree(s);
}
// (m,k,j) deepest cell the band solver integrates; the face below it is where the
// interior flux is handed in from the diffusion operator.
inline DvceArray3D<int> rt_cut_index() { return *rt_icut_ptr; }
// problem/rt_de_max: the cap in LimitRTSource, as a fraction of the cell's internal
// energy per RT application. Applies to every EXPLICIT radiative update -- grey and
// correlated-k, split and monolithic. Set <= 0 to disable the limiter entirely.
inline Real rt_de_max = 0.5;
// problem/rt_floor_consistent: THE RADIATIVE SOURCE MAY NOT COOL A CELL BELOW THE STATE
// THE EOS FLOORS WILL PUT IT BACK AT.  Default FALSE = today's arithmetic, bit for bit.
//
// In the tenuous NIGHT-SIDE TOP SLAB of a hot Jupiter -- density on dfloor, pressure
// within a factor of a few of pfloor, every correlated-k band optically thin -- the
// two-stream's equilibrium for the topmost cell lies BELOW the pressure floor.  The two
// operators then fight each other once per stage, forever: the sweep takes energy out,
// ConsToPrim's pressure floor puts it straight back (this is the eos_efloor counter), and
// the only trace left is a spurious energy source in the floor and, when the step is
// large enough, the rt_de_max clip.  MEASURED on the shrunk production grid
// (inputs/tests/dhj_ck_spherical.athinput, 4000 cycles): the centre-to-centre layers give
// 27048 eos_efloor events against 0 for the staggered ones, and lose 1.7e-5 of the total
// energy; from one shared state the same A/B is 16 floor events per cycle against 0.
// The layers are NOT what is wrong -- their emission is the cell's own kappa rho B over
// its own two halves, while the staggered layer averaged it with the GHOST cell above the
// domain and so under-emitted the top cell by 46 % -- and neither is the floor.  What was
// missing is that the apply did not know the floors exist.
//
// With this on, the decrement is bounded by the floor state itself: e^{n+1} >= e_floor =
// max(e(rho, tfloor), e(rho, pfloor)), and a cell already at or below it is handed
// nothing.  A HEATING step is untouched, as is every cell that stays clear of the floors,
// so the change is confined to exactly the cells the floor was catching.  It costs one
// EOS evaluation per COOLING cell and stage.
inline bool rt_floor_consistent = false;
// problem/rt_semi_lin: recover the OLD semi-implicit step, which linearized the emission
// about the current state and relaxed at lambda = 4E/e.  See the long note in rt_apply:
// that form bounds cooling but leaves heating explicit and unbounded, so a cold optically
// thin cell overshoots its equilibrium badly.  Kept only to reproduce runs made before
// the fix.  DEFAULT TRUE, i.e. the old linearization, so that every problem generator
// sharing this header reproduces its pre-fix runs bit for bit; the red-giant inputs set
// it false to select the new equilibrium-relaxation step.
inline bool rt_semi_lin = true;
// problem/rt_use_cons: take the cell's internal energy and density from the CONSERVED
// state u0 instead of from w0.  w0 is the previous stage's ConToPrim output, and by the
// time the user source function runs RKUpdate, the explicit source terms and the
// implicit radial conduction have all moved u0.  In a smooth cell the two agree to
// O(dt); in a runaway cell the hydro step changes e by ~100 % per stage, so the
// equilibrium, the Newton iterate, the positivity guard and LimitRTSource are all
// evaluated on a number that no longer exists while the step is applied to u0 anyway --
// which is how a cell is driven to a negative internal energy, repaired by the floor,
// and handed to the next Riemann solve as a 0.1 K cell beside a 1e4 K one.  The
// namespace value below is only the pre-pgen initialiser: every problem generator that
// uses this header sets it from problem/rt_use_cons, and since 2026-09-22 the default
// there is TRUE in box_convection and deep_hot_jupiter_rt, and the mode-3 flag in
// red_giant.
inline bool rt_use_cons = false;
// problem/rt_test_mu0 and problem/ck_test_kgrey: TEST HOOKS for the well-posed checks
// of the implicit correlated-k solve (tests_ck_implicit/wellposed).  rt_test_mu0 >= -1
// gives every correlated-k column the same stellar cos(zenith) (uniform irradiation);
// ck_test_kgrey > 0 replaces the per-band continuum opacity of every band by this
// constant [cm^2/g], which with a k-table of negligible line opacity makes the sweep
// GREY.
// Defaults (-2, 0) are off, and off they change nothing: both are one untaken branch.
inline Real rt_test_mu0 = -2.0;
inline Real ck_test_kgrey = 0.0;
// problem/rt_bface: the emissivity-weighted far-endpoint Planck source (see BFace).
// It is a red-giant fix -- it exists because the corona/star join put a cell's emission
// on a neighbour's Planck function 1e9 times its own -- and this header is shared with
// solar_convection and the hot-Jupiter problems, which have no such join.  DEFAULT FALSE
// therefore returns the old endpoint b_far bit for bit at every one of the seven call
// sites; the red-giant inputs set it true.
inline bool rt_bface = false;
// how many cells the Newton loop had to be rescued from a non-positive internal energy
// (see rt_apply).  A device counter, read back where the clip count is reported.
inline DvceArray1D<int> *rt_efix_ptr = nullptr;
// ...and the STATE of the first cell it happened to, which is the only way to tell a
// genuine stiff-cooling cell from a poisoned neighbour without re-running.  Filled once
// per run by whichever cell claims efix(0) == 0, printed with the warning.
inline DvceArray1D<Real> *rt_efix_rec = nullptr;
// the column energy budget of the source application, accumulated when rt_cell_report is
// on: slot 0 is the EXPLICIT deposit sum(src*bdt*dx1), which telescopes to F_bot - F_top
// over a column, and slot 1 the deposit actually applied, sum(de*dx1).  The two differ by
// whatever the relaxation, the sub-cycle and LimitRTSource did, so their ratio is the
// conservation statement for problem/rt_relax_sub.
inline DvceArray1D<Real> *rt_desum_ptr = nullptr;
// problem/rt_newton: refine the semi-implicit step with a Newton solve of the exact
// backward-Euler balance, seeded by the closed form.  The closed form already lands on
// the right equilibrium under e ~ T; this drops that assumption and uses the EOS's own
// T(e) and c_v(e), which matters wherever H2 or H is partly dissociated.  Costs one
// Temperature() and one SpecificHeatCv() per iteration, and from that seed it is
// normally one or two.  General EOS only; with an ideal gas e ~ T is exact and the
// closed form is already the answer.  DEFAULT FALSE so shared pgens are unchanged; the
// red-giant inputs set it true.
inline bool rt_newton = false;
// problem/rt_rescue_eq: when the Newton step would leave e <= 0, land the cell on the
// radiative equilibrium it actually sees -- deq = ei((A/Em)^(1/4) - 1), i.e. Em(T_eq) = A
// with A = src + Em the absorption, which in the thin limit is (pi/mu) kappa rho
// (I_up + I_dn) so that sigma T_eq^4 = pi (I_up + I_dn)/2 up to the diffusivity factor --
// instead of the unconditional 99.9 % drop.  Floored at that same 99.9 %, so it can only
// make the rescue less violent.  DEFAULT FALSE = the old 99.9 % policy, so shared pgens
// are unchanged; the red-giant inputs set it true.
inline bool rt_rescue_eq = false;
// problem/rt_relax_sub: SUB-CYCLE the local relaxation of the two-stream source.  The
// single-step form damps each cell toward the equilibrium of the un-relaxed column by a
// per-cell factor (1 - exp(-x))/x; neighbours at different x are damped by different
// amounts, so the non-local exchange between them no longer cancels and what is left is
// a pressure perturbation ~ dt x the compression rate, which pumps the standing acoustic
// modes of a closed box (the He-star box's saturated v_rms was LINEAR in dt).  With this
// set > 1 the local balance is integrated in nsub sub-steps of bdt/nsub instead, the
// absorbed field A frozen for the stage and the cell's own emission re-formed from the
// running energy after each sub-step, so every sub-step is taken at small x.  The value
// is the FLOOR on nsub; the stiffness rule below can raise it.  DEFAULT 1 = the
// single-step form, BIT FOR BIT.
// problem/rt_ali_diag: the ACCELERATED-LAMBDA (ALI) DIAGONAL in the per-cell
// semi-implicit apply.  Without it the step relaxes the cell toward A, the absorption
// the sweep computed, at the rate 4E/e set by the cell's OWN emission -- and in an
// optically thick cell that rate is enormous, so x = src*bdt/deq reaches 1e4-1e5 and the
// factor (1 - exp(-x))/x delivers 1/x of the source to the gas.  MEASURED on the B-star
// relaxed column (bench/bstar_fecz/leakdiag/dump0, arm d, w = 0 everywhere, dtau 5-134):
// the assembled source sums to +0.2931 F_bot = F_bot - F_top exactly, and the APPLIED
// sum is +0.000089 F_bot.  Nothing is clipped, no guard fires, the Newton converges in
// 1-2 steps -- the source is simply damped away, cell by cell, and the energy the sweep
// took out of the radiation field is never given to the gas.  That is the 0.29 F_bot
// deficit: Ftop/F_bot froze at 0.707 because the column can never re-relax.
//
// The mistake is treating A as an EXTERNAL field.  In a thick cell almost all of A is
// the cell's own emission coming back from its immediate neighbours, so the true
// response of the net exchange to a change in this cell's source function is not 4E/e
// but (1 - Lambda*_ii) 4E/e, with Lambda*_ii the diagonal of the Lambda operator.  That
// is textbook accelerated Lambda iteration; the sub-cycle and the Newton refinement are
// both fixes to the wrong rate rather than to the rate itself.
//
// Lambda*_ii here is the local escape-probability diagonal the sweep's own layer
// coefficients already carry, 1 - (1 - e^-x)/x at x = dtau/mu, averaged over the
// quadrature.  It goes to 1 as dtau -> infinity, which sends the relaxation rate to zero
// and the step to the EXACT explicit one, de = src*bdt -- conservative, and stable
// because a thick cell's net imbalance relaxes on the diffusion time, not the thermal
// one -- and to 0 as dtau -> 0, where the old stiff relaxation is what is wanted and is
// recovered bitwise.  rt_ali_diag = false restores the old behaviour exactly.
inline bool rt_ali_diag = true;
inline int rt_relax_sub = 1;
// problem/rt_relax_xcrit: the per-sub-step stiffness the sub-cycle aims for.  nsub is
// ceil(x/rt_relax_xcrit) with x = src_relax*bdt/deq the first step's stiffness, floored
// at rt_relax_sub.  Inert unless rt_relax_sub > 1.
inline Real rt_relax_xcrit = 1.0;
// problem/rt_relax_submax: the cap on nsub, so a single pathological cell cannot cost an
// unbounded number of EOS calls.  Inert unless rt_relax_sub > 1.
inline int rt_relax_submax = 32;
// problem/rt_src_direct: form the per-cell radiative source DIRECTLY as absorption minus
// emission during the sweeps, instead of as the difference of the two face fluxes.  The
// two are the same number algebraically -- the stream update across a layer is
//     I_out = (1 - e0) I_in + (alp B_far + bet B_near),   alp + bet = e0,
// so the flux change over the layer is e0 (I_in - B) exactly, and summing those over a
// column telescopes to F_top - F_bot -- but not numerically.  In a transparent cell the
// two face fluxes agree to every digit and their difference is round-off: measured
// 1.7e-11 relative on 5.67e9, i.e. 5.8e-12 erg/cm^3/s against a cell whose whole
// internal energy is 6.3e-12.  That noise, of either sign, is what heated the ambient
// medium to 3000 K in two steps and then drove it to zero; no limiter or opacity floor
// touches it, because it is not a physical term at all.  Formed directly the source is
// O(dtau) with no large numbers cancelling, and it vanishes as the cell goes transparent
// the way the physics says it must.  Off reproduces the flux-difference form bit for bit,
// and is the DEFAULT for that reason; the red-giant inputs set it true.
inline bool rt_src_direct = false;
// problem/rt_top_clamp: the band solver's top slot i = ie+1 takes its fluid state from
// the top ACTIVE cell ie instead of from the hydro ghost (see rt_pre_tp).  It decouples
// the whole column from whatever the outer boundary condition put in the ghost, which is
// what a poisoned ghost needs, but it changes numbers in any run whose top ghost is a
// legitimate state.  DEFAULT FALSE = read the ghost exactly as before; the red-giant
// inputs set it true.
inline bool rt_top_clamp = false;
// problem/rt_semi_implicit: apply the split two-stream source semi-implicitly, relaxing
// the cell toward radiative equilibrium instead of stepping explicitly (see the apply
// kernel).  Default true, which is the behaviour since 048dff30.  The flag exists
// because the semi-implicit step is NOT answer-preserving in general: on the cubed
// sphere deep_hot_jupiter 09-07 configuration it moves the kinetic energy by 3 % in 20
// cycles (A/B 2026-09-09).  Runs that must reproduce pre-048dff30 results set it false,
// which restores the plain explicit de = src*bdt.
inline bool rt_semi_implicit = true;
// problem/rt_outer_verbose: print the per-cycle clip/rescue census and the mode-3
// convergence lines even without problem/rt_cell_report.
inline bool rt_outer_verbose = false;
// ---- problem/rt_implicit_column: THE IMPLICIT COLUMN SOLVE --------------------------
// 0 = the per-cell semi-implicit relaxation of the apply block (the historical default,
// bitwise).  3 = the EXACT block-tridiagonal column solve below, which is what both
// stellar productions run.  The earlier linearised modes 1 and 2 were removed; the
// problem generators fatal on any other value.
inline int rt_implicit_column = 0;
// convergence tolerance and pass cap of that Newton
inline Real rt_impl_tol = 1.0e-6;
inline int rt_impl_maxit = 5;
// problem/rt_col3_ex_iter: under rt_implicit_column = 3, re-form the tau-blend handover
// src_ex from the COLUMN'S OWN converged face flux instead of the entry sweep's, so that
// the applied source is the divergence of one field, div[(1-w)F_3], and telescopes.  See
// two_stream_column_implicit.hpp, step 4a'.  Off = bitwise the frozen handover.
// Read ONLY at the c3.ex_iter assignment below, inside the rt_implicit_column = 3
// branch, so it is inert in every other mode.  problem/rt_col3_ex_iter defaults TRUE in
// box_convection and red_giant since 2026-09-22; deep_hot_jupiter_rt never reads it
// (it refuses rt_implicit_column != 0 outright).
inline bool rt_col3_ex_iter = false;
// problem/rt_col3_skip_sweep: under rt_implicit_column = 3 with the tau blend weight w
// = 0 on EVERY face (rad_tau_lo deeper than the whole box, i.e. the rt_bottom_flux
// production configuration), SKIP the explicit entry sweep, which the column no longer
// consumes.  With w = 0 the handover src_ex is 0 on every cell and, under
// rt_col3_ex_iter, the column re-forms it from its own flux anyway: mode 3 then reads
// only kc, Bb, icut, Qb (identically 0 on the grey path) and the frozen top-face
// intensity, all of which the PRE-kernels and rt_c3_top build without the sweep.  The
// sweep's own products -- Fb (the face flux every diagnostic, rad_f2s, rt_rad_force and
// the surface/history dumps read), Src and Em -- are then supplied differently: Fb is
// written by the COLUMN SOLVE from its own converged D/U intensities (which is the field
// mode 3 actually applies, not the entry-state one), while Src and Em stay at zero, so
// every consumer of THOSE two must be off.  box_convection.cpp enforces the whole list.
// Default false = bitwise off.
inline bool rt_col3_skip_sweep = false;
// ---- mode 3 CONVERGENCE, the switches of the acceleration pass ----------------------
// problem/rt_impl_exjac: carry the rt_col3_ex_iter handover in the JACOBIAN instead of
// lagging it.  With ex_iter on the applied source is A_i = Src_i + div[w F_3]_i + Q_i,
// and div[w F_3] is LINEAR in the very unknowns the block system already carries, so
// there is no excuse for lagging it: off, the missing term makes the whole Newton a
// PICARD iteration wherever w is neither 0 nor 1, and the measured convergence is linear
// at ~0.4 per pass (13-16 passes to 1e-12, bench/bstar_fecz/m3tol).  On, the block rows
// are the exact derivative.  Off = the old lagged assembly, bitwise.
//
// MEASURED, bench/bstar_fecz/m3acc/sweep, the relaxed 15 Msun B-star column at cycle 0,
// passes taken to reach a given residual (norm 0, no step-size stop, maxit 30):
//     tol     1e-2   1e-4   1e-6   1e-8   1e-10  1e-12
//   exjac 0     2      4      7     10     13     16      LINEAR, rate 0.215 per pass
//   exjac 1     2      3      3      3      4      4      QUADRATIC: 1.6e-4, 4.3e-9,
//                                                         5.4e-14, i.e. round-off at 4
// Both reach the SAME fixed point (max|du/u| 6.350328e-02 either way at tol 1e-14), so
// this is purely the path.  Cost of the whole run, 2000 cycles on 1 MI300A:
//   old default (tol 1e-6, maxit 6, never converged, exits on maxit)   1.53e5 zc/s
//   old converged (tol 1e-12, maxit 20, 13-16 passes)                  8.36e4 zc/s
//   new default (exjac, tol 1e-8, maxit 8, 3 passes, resid 5e-10)      2.15e5 zc/s
// against 3.96e5 zc/s for the standard scheme (rt_implicit_column = 0), so the CONVERGED
// mode 3 now costs 1.84x the standard scheme where the unconverged one cost 2.58x.
inline bool rt_impl_exjac = true;
// problem/rt_impl_norm: the residual norm the stopping rule measures.  0 = the old
// max_i |R_i|/e_i, which the optically thin TOP of a stellar column dominates simply
// because its e is 5e-6 of the base's -- a fixed relative tolerance there demands an
// absolute precision nothing else in the scheme has.  1 = max_i |R_i|/(e_i + eps e_max),
// eps = problem/rt_impl_norm_eps, so a cell holding a negligible share of the column's
// energy is measured against that share and not against itself.
inline int rt_impl_norm = 1;
inline Real rt_impl_norm_eps = 1.0e-3;
// problem/rt_impl_dstop: also stop when the Newton STEP is small, max_i |db_i/b_i| < tol.
// A converged Newton takes a negligible step; insisting on the residual as well costs a
// whole extra pass whose only effect is to confirm it.
inline bool rt_impl_dstop = true;
// problem/rt_impl_rescheck: form the energy-row residual BEFORE the block factorisation
// so that the pass which only confirms convergence never factorises.  Bitwise identical
// to the standard path (same test, same value, same iterate); see the note on RTCol3.
inline bool rt_impl_rescheck = false;
// problem/rt_impl_reuse: REUSE THE BLOCK FACTORISATION across Newton passes.  The
// Jacobian of the column system is constant in the iterate except for the single entry
// de/db (the gas heat capacity), and the right-hand side is one number per cell -- the
// energy-row residual -- so a pass that keeps pass 1's factors costs a residual
// evaluation and two small matrix-vector products instead of the whole assembly,
// factorisation and reduced elimination.  1 = reuse with a contraction check (the pass
// refactorises when the residual failed to fall by a factor 0.3), 2 = reuse always
// (diagnostic).  0 (the default) is the old code, bitwise.  The converged answer is the
// same root of the same residual to the same tolerance; only the iterates between the
// first and the last differ, at O(the frozen de/db).  The
// PARTITIONED path (rt_impl_solver = pcr) only; thomas is refused.
inline int rt_impl_reuse = 0;
// problem/rt_impl_reuse_rho: the contraction a reuse pass must show for the NEXT pass to
// keep the factorisation.  The residual of a full Newton pass falls much faster than of
// a frozen-Jacobian one, so a loose threshold buys cheap passes at the price of more of
// them; tighten it to refactorise sooner.
inline Real rt_impl_reuse_rho = 0.3;
// ---- rt_implicit_column = 3: THE EXACT IMPLICIT COLUMN SOLVE ------------------------
// The intensities become unknowns alongside the gas energy, so the column system is
// exactly block-tridiagonal (5x5) and a block Thomas solves it with no Jacobi/
// Gauss-Seidel iteration at all.  See
// utils/two_stream_column_implicit.hpp for the rows, the M-matrix argument and what is
// frozen.  Mode 3 REPLACES the two-stream's own semi-implicit apply and nothing else:
// the tau-blend handover, the radial conduction tridiagonal (rad_implicit_x1) and the
// transverse operator are untouched, so the conduction wrapper task still runs.
// Grey only; fatal under correlated-k or picket fence.
// problem/rt_impl_solver: HOW the block-tridiagonal column system is solved.
//   thomas (0, the default)  one thread per column, the serial block Thomas of
//                            two_stream_column_implicit.hpp -- the reference bit pattern
//   pcr    (1)               one TEAM per column: the column is cut into
//                            problem/rt_impl_nseg segments, each thread eliminates its
//                            own segment with the incoming unknown carried symbolically,
//                            and the resulting 5x5 block system over the segment
//                            BOUNDARIES is solved by one thread and back-substituted in
//                            parallel.  Same system, same Newton, same clamps; the two
//                            agree to round-off.  See two_stream_column_partition.hpp.
inline int rt_impl_solver = 0;
// problem/rt_impl_mixed: MIXED PRECISION in the partitioned column solve's block algebra.
//   0 (default)  everything double -- bit for bit the code before the switch existed
//   1            the Newton CORRECTION is formed in SINGLE precision: the 5x5 inverses
//                (RTCol3Inv5X), the segment forward eliminations, the p/Q/R recurrence
//                and the back-substitutions.  The residual, the convergence test, the
//                layer/Planck/kappa coefficients, the reduced system and the u0(IEN)
//                update stay double, so the outer Newton is an iterative refinement of
//                the float factor and converges to the same rt_impl_tol.
//   2            ... and the STORED factors (G, H, d, and the reuse pair v, M) live in
//                their own float workspace, which halves their memory traffic.
// A block whose float Gauss-Jordan loses its pivots is redone in double and counted in
// stat(24) (nmixfb under rt_outer_verbose).  Partitioned path only.
inline int rt_impl_mixed = 0;
inline int rt_impl_nseg = 64;
// problem/rt_impl_redpar: on the pcr path, the REDUCED block-tridiagonal system over the
// segment boundaries (nsg 5x5 rows, non-periodic) is by default eliminated serially by
// one lane while the other nsg-1 idle -- ~29 % of a Newton pass at nseg = 32.  With this
// switch it is solved by PARALLEL CYCLIC REDUCTION instead: lane s owns block row s and
// log2(nsg) rounds of neighbour eliminations leave y_s = B_s^-1 r_s on every lane at
// once, with no reduced back-substitution.  Same system, different elimination order, so
// the answer is equal to ROUND-OFF, not bitwise.  Falls back to the serial solve unless
// 2 <= nsg <= 64 and nsg is a power of two.  Costs the larger per-segment workspace
// (RTCOL3_NRDP), so the array is only enlarged when the switch is on.
inline bool rt_impl_redpar = false;
// problem/rt_col3_hybrid_tau: THE HYBRID COLUMN.  Below this COLUMN optical depth the
// cells are solved as DIFFUSION -- one unknown per cell, F = (4 pi/3) dB/dtau across the
// faces, which is the exact deep limit of the same two-stream -- instead of as the full
// five-unknown two-stream, all inside the SAME Newton system and the same block
// elimination (1x1 blocks deep, 5x5 thin, a 1x5/5x1 pair at the interface).  The two
// segments are coupled by flux continuity: the thin segment's lower boundary is the deep
// limit U_q = b_f + mu_q dB/dtau with b_f and dB/dtau carried implicitly by the deep
// unknowns, and the deep segment's top face loses exactly the two-stream's own net flux
// at that face.  0 = off, bitwise the whole-column solve.  See
// utils/two_stream_column_implicit.hpp, RTCol3Hyb.  Supported on BOTH solvers: with
// rt_impl_solver = pcr the deep segments run a scalar partitioned Thomas and the reduced
// system has mixed block sizes.  See two_stream_column_partition.hpp.
inline Real rt_col3_hybrid_tau = 0.0;
// problem/rt_impl_warm: WARM-START the mode-3 Newton from the previous call's converged
// Planck function instead of from the entry state's.  0 = off (bitwise the old code),
// 1 = the previous converged b per cell, 2 = linear extrapolation in time from the last
// two.  Only the ITERATE moves; the residual, the Jacobian, the tolerance and the entry
// state are untouched, so the converged root is the same one and the answer is equal to
// round-off rather than bitwise.  Costs one (warm = 1) or two (warm = 2) extra 4D arrays.
// rt_impl_warm, rt_c3bp_ptr, rt_c3bp2_ptr and rt_c3_bdt_prev are the warm start's own
// state and live in utils/two_stream_warm_rst.hpp, which the restart writer and reader
// also include: the history is carried in the restart file.
inline DvceArray5D<Real> *rt_c3wk_ptr = nullptr;
inline DvceArray5D<float> *rt_c3wkf_ptr = nullptr;
inline DvceArray4D<Real> *rt_c3rd_ptr = nullptr;
inline DvceArray4D<Real> *rt_c3top_ptr = nullptr;
inline DvceArray1D<Real> *rt_c3stat_ptr = nullptr;
inline int rt_c3_lines = 0;
// problem/ck_int_at_cut: deliver the planet's internal flux sigma T_int^4 as an extra
// upward source at the correlated-k cut (the historical behaviour, true). Set false when
// the layers below the cut carry it themselves -- <mhd|hydro>/isotropic_conduction =
// radiative with rad_flux_inner at the bottom wall -- so the cut's upward intensity is
// just the thermalised Planck function and nothing is counted twice.
inline bool rt_int_at_cut = true;
// problem/rt_bottom_flux: THE IMPOSED FLUX ON THE TWO-STREAM'S OWN LOWER BOUNDARY.
//
// With the tau blend switched off (rad_tau_lo/hi beyond the bottom of the box, so the
// blend weight is 0 on every face) the cut falls on the inner x1 wall and the two-stream
// owns the WHOLE column: the radiative-diffusion operator has nothing left to carry, and
// the handover that splits the two is gone.  The energy then has to enter through the
// two-stream's own lower boundary instead of through the conduction wall face.
//
// This is the code-unit flux that boundary carries.  The cut's upward intensity is the
// deep limit I_up(mu) = B + mu dB/dtau with the gradient set BY THE FLUX rather than read
// off the two cells above it,
//     dB/dtau = 3 F_bot/(4 pi),
// which with the two-point Gauss-Legendre quadrature returns sum_q w_q (I_up - I_down) =
// F_bot exactly (rt_nquad = 1 returns 0.904 of it, the hemispheric mean's own error).
// 0 = off, the local deep-limit gradient.  box_convection.cpp sets it
// from <hydro>/rad_flux_inner under <problem>/rt_bottom_flux, and zeroes rad_flux_inner
// at the same time so the flux is not also injected at the conduction wall face.
inline Real rt_bot_flux = 0.0;
// problem/rt_layer_legacy: the STAGGERED layer the grey and correlated-k sweeps used.
//
// Every layer was integrated over a WHOLE cell in optical depth, dtau = kappa rho dz of
// cell i, but with its source running linearly in tau from that cell's CENTRE Planck
// function to the NEIGHBOUR CELL'S CENTRE one.  The interval the optical depth measures
// and the interval the source runs over are then offset by half a cell, and the two
// endpoints of the source are two cells apart while the layer is one cell thick.  The
// error is O(dtau^2) and smooth, but it is not small where it matters: MEASURED against
// an accurate formal solution on the He-star column, the emergent flux came out ~4 % low
// at dtau/cell of 1-3 and wandered between 0.73 and 1.06 of the true flux at dtau/cell
// 3-6, which is exactly the handover ramp of the He and B-star boxes.
//
// THE FIX.  The layers now run between cell CENTRES: layer (i, i+1) has
//     dtau = 0.5 (kappa rho dz)_i + 0.5 (kappa rho dz)_{i+1},
// each half with its OWN opacity, and a source linear in tau between B_i and B_{i+1} --
// which is exact for that layer, since the source really is (piecewise) linear there.
// The layer is integrated as its two halves, split at the face, with the face source
// taken at its own tau within the layer; composing the two halves reproduces the whole
// layer exactly, so this costs nothing in accuracy and hands back the FACE intensity as a
// by-product.  That face intensity is what Fb_g reports, so the fluxes the tau blend and
// the heating read are the formal solution AT the face, not an interpolation of it.  The
// column is closed by two HALF layers: the upper half of cell ie, from the top face to
// its centre, with the source held at B(ie); and the lower half of cell icut, between its
// centre and the cut face, with B continued at the deep-limit gradient dB/dtau.
//
// Each half layer lies entirely inside ONE cell, so the absorbed and emitted amounts that
// build Src_g are attributed unambiguously: a cell's emission is simply split between its
// two halves.  Em_g, the relaxation rate, becomes the cell's own 4 sigma kappa rho T^4
// instead of the centre-to-centre average of two Planck functions.
//
// Default false, i.e. THE FIX IS ON.  Set true to reproduce a pre-fix run bitwise.
inline bool rt_layer_legacy = false;
// problem/ck_spherical: give the CORRELATED-K thermal two-stream the SPHERICAL form of
// 8bca3dfa, the one the grey centre-to-centre sweep and the mode-3 column solve already
// carry.  Since 2026-09-22 problem/ck_spherical DEFAULTS TRUE on a spherical-polar or
// cubed-sphere mesh (and only under rt_ck && !rt_layer_legacy, the two things the guard
// below requires) and is REFUSED on a Cartesian one; false is bitwise the
// plane-parallel arithmetic the ck kernel used to run.
//
// WHY IT IS NEEDED NOW.  8bca3dfa converted the grey path and documented the ck kernel
// as "not converted", justified by r_out/r_in < 1.3 hot-Jupiter domains.  The deep
// hot-Jupiter production grid spans x1 = 9.44e9 .. 2.0556e10 cm, r_out/r_in = 2.18 and
// an AREA RATIO of 4.75, so on that grid the unconverted kernel (a) does not dilute the
// thermal flux, i.e. L = A F grows like r^2 through a transparent shell, and (b) misses
// the -2F/r term of the deposit, so the thin upper atmosphere is irradiated from below
// by an undiluted thermal flux.
//
// THE FORM.  Exactly the one derived on RTCol3::Bt in two_stream_column_implicit.hpp and
// reused by the grey sweep: inside a cell the frame area A_i = V_i/dx_i is constant and
// the layer solve is the plane-parallel one, per unit area, against the source function
// B; ALL the geometry sits in the face conditions, where S = (I+ + I-)/2 and A D =
// A (I+ - I-)/2 are passed continuously,
//     u_above = u_below + c,  d_below = d_above + c,  c = beta (d_above - u_below),
//     beta = (A_above - A_below)/(A_above + A_below),
// which is the BTF/MIXF pair below.  Consequences, all of them shared with the grey
// path: Phi = A F is single-valued at every face, so Src_g stays the plane-parallel
// divergence (nothing in the deposit changes) and telescopes to -(A_t F_t - A_b F_b)/V
// exactly; Fb_g stays a flux PER UNIT AREA at its own face, reported as Phi_f/A(f).
//
// THE DOWN RAY IS LAGGED, as it is in the grey mode-0 sweep and for the same reason: the
// down-sweep's mixing needs the UPWARD intensity at each face, which it does not yet
// have, so a plane-parallel PROBE up-sweep supplies it to O(beta), leaving the down
// ray's mixing exact to O(beta^2).  The probe costs one extra ck_kappa lookup per cell
// and chain -- roughly +50 % on the chain kernel's table traffic -- which is why this is
// a switch and not unconditional.  The probe is stored IN the I_down column (each face
// slot is read one statement before the real down-sweep overwrites it), so it costs no
// extra per-thread memory.
//
// THE STELLAR BEAM IS NOT TOUCHED AT ALL.  This switch converts the THERMAL two-stream
// and nothing else.  The beam is a parallel pencil from a source at infinity, and a
// parallel pencil is NOT confined to a spherically diverging column: the column widens
// with r, the beam does not, so rays entering the column's top face leave through its
// sides.  The power a column intercepts at its top face is therefore NOT the power it
// absorbs, and the only exact local statement is the volumetric absorption rate
// kappa rho F_star exp(-tau_ray) -- which is exactly what the kernel already computes,
// written as a flux difference across the cell so that the column integral comes out
// F_star mu0 (1 - e^-tau_tot) by construction.  Neither an r^-2 dilution nor a
// conserved-power form belongs here.
//
// RETRACTED.  The first revision of this switch carried the beam as a conserved power
// Phi_beam(f) = A(ie+1) F_star mu0 T(f), on the argument that a column must deposit every
// watt it intercepts.  It must not.  That made every dayside column absorb the power
// crossing the DOMAIN TOP, i.e. the planet absorbed pi r_top^2 F_star instead of
// pi r_abs^2 F_star with r_abs the radius where the slant optical depth reaches 1 --
// an over-heating by (r_top/r_abs)^2, a factor ~3 on the production grid.  The claim in
// tests_ck_sph/README.md that the plane-parallel form "under-heats the dayside by a
// factor 2.6" was that error seen from the wrong side, and is retracted there.
//
// NOT CONVERTED, and refused: rt_layer_legacy (the staggered whole-cell ck layers) and
// the mode-1/2 nearest-neighbour Jacobian, which the ck path never fills anyway.
inline bool ck_spherical = false;
// problem/ck_beam_sph: the PSEUDO-SPHERICAL direct stellar beam.  Since 2026-09-22 it
// follows the same rule as ck_spherical: DEFAULT TRUE on a spherical-polar or
// cubed-sphere mesh (under rt_ck && !rt_layer_legacy), REFUSED on a Cartesian one.
// False is bitwise the plane-parallel slant path the correlated-k kernel used to use.
// Independent of ck_spherical (that switch converts the THERMAL two-stream; this one the
// BEAM), and the two compose.
//
// WHAT IS WRONG WITHOUT IT.  The beam's optical depth to a point is tau_vertical/mu0
// through the column's own profile, with mu0 CLAMPED at 0.1 (facsw below), and mu0 <= 0
// is simply dark.  On a 2.18 domain that is wrong three ways: the secant law diverges as
// the terminator is approached and the clamp then silently caps it at 10x the vertical
// depth; the true slant path through curved shells is finite at mu0 = 0 (the Chapman
// function goes to sqrt(pi X/2), ~30 for X = r/H ~ 2000, not to infinity); and the
// TWILIGHT region mu0 < 0, which on a real planet is lit down to the radius where the ray
// grazes the opaque body, receives nothing at all.
//
// THE GEOMETRY (pseudo-spherical: straight rays, the column's own profile assumed on
// every shell -- horizontally homogeneous, which is the standard approximation and is
// what makes this a per-column calculation).  For a target at radius r with stellar
// zenith angle theta0, cos theta0 = mu0, the impact parameter of the ray is
//     b = r sin(theta0) = r sqrt(1 - mu0^2),
// and a shell between faces r_j and r_{j+1} contributes the chord
//     ds_j = sqrt(r_{j+1}^2 - b^2) - sqrt(max(r_j^2 - b^2, 0)),
// the max() handling the shell that contains the tangent point exactly.  Then
//     tau_ray(r) = sum_j w_j ds_j (kappa rho)_j,
// with w_j = 1 for mu0 >= 0 and j from the target's own shell upward (the target's shell
// contributes from r itself, which is what sqrt(r^2 - b^2) = r|mu0| gives), and for
// mu0 < 0 the ray reaches the target through the TANGENT POINT at radius b: a NEAR leg
// from r down to b and a FAR leg from b out to the top, so w_j = 2 for the shells below
// the target (b <= r_j < r) and w_j = 1 above it, the sum starting at the shell that
// contains b.  A mu0 < 0 ray is DARK if b <= the radius at which the column becomes
// opaque, which here is the correlated-k cut face x1f(icut) (below it the grey optical
// depth is ~1e4): the ray would have to pass through the planet.
//
// Above the domain top there is NOTHING: tau_ray(r_top) = 0.  The plane-parallel path
// instead starts the beam with the unresolved hydrostatic ghost column's optical depth
// (RTTopDtau); that term is ~2e-4 of the incident flux on this grid and is dropped here,
// because a pseudo-spherical ray leaves the domain and the ghost column is a
// plane-parallel construction with no radial extent.
//
// THE DEPOSIT IS UNCHANGED IN FORM -- the local absorption rate, no column-power
// bookkeeping (see the retraction on the beam in 2db7095c).  Written as a flux difference
// across the cell so that a thick layer does not under-deposit,
//     Qb_i = (1-albedo) F_star w_b (kappa rho)_i [e^-tau(r_i) - e^-tau(r_{i+1})]
//                                               / [tau(r_{i+1}) - tau(r_i)],
// which is IDENTICALLY the plane-parallel expression when tau(r) = tau_vert(r)/mu0 (the
// bracket is then dtau_i/mu0 and (kappa rho)_i/(dtau_i/mu0) = mu0/dz_i), and which never
// divides by a chord that disagrees with the optical depths it is paired with.
//
// NO MATRIX IS STORED.  The chords depend only on (r, b) and b depends on the column only
// through mu0, so a per-column triangular G_ij = ds_ij would be 128*129/2 entries x 6144
// columns x 4 B = 203 MB on the production cubed-sphere mesh (6 x 32 x 32 x 128).  The
// on-the-fly form costs ONE sqrt per (face, shell) pair, hoisted out of the g-point loop
// because b does not depend on the chain, and needs no memory at all; on a GPU that trade
// is obviously right.  mu0 is a pure function of the grid (rt_pre_geom: the substellar
// point is fixed at phi = 0 in the mesh, this problem generator has no non-synchronous
// option), so nothing here can go stale either.
//
// The only new per-thread storage is the (kappa rho) column, NC x NN, needed because
// tau_ray(r_i) is not a running sum -- every target has its own chord set.  It is sized
// to one element when the switch is off: the kernel is templated on the flag.
inline bool ck_beam_sph = false;
// problem/ck_beam_par (ck-fast2 lever 3): the pseudo-spherical beam of the tm sweep's
// STORING pass in its own kernels, a thread per target face instead of a thread per
// column (the in-kernel block walks O(N^2) chords serially: measured 48 of the storing
// kernel's 94 ms per call at nx1 256).  Same chords, same order, same deposit arithmetic;
// the in-kernel block's "stop once the block's slant depth passes 60" becomes a per-chain
// stop at VERTICAL depth 60 (a lower bound of the slant one), so the two differ only by
// deposits below e^-60 of the incident beam.  Needs ck_impl_lin, lin_thr = 1 (scratch).
inline bool ck_beam_par = false;
inline DvceArray4D<int> *ck_bstop_ptr = nullptr;   // (m, beam chain, k, j)
// ck-next: ck_coef's per-cell inputs xT, xP, rho, dz as (m, q, i, k, j), j fastest
inline DvceArray5D<Real> *ck_cin_ptr = nullptr;
// problem/ck_store_split (ck-store): on the storing pass of the linear re-apply, store the
// opacity-only operator with a per-(band, g) kernel (ck_coef) and take the pass's fluxes
// from the linear kernels, instead of running rt_chain_ck (see ckstsp_ in the pass).  The
// operator is bitwise the chain kernel's; the storing pass's Src and Fb change by
// round-off.  Only under ck_beam_par (and the conditions of ckstsp_); true by default.
inline bool ck_store_split = true;
// The beam's optical depth to the TOP face of a mu0 < 0 (twilight) column: the ray
// enters the domain top on the far side, grazes the tangent radius b = r_top sin(theta0)
// and climbs back to the target, so every shell between b and r_top is crossed TWICE.
// Dark (tau = 1e30, transmission 0) if b <= r_cut: the ray would pass through the body.
// Without this the top face was given tau = 0 (full, unattenuated beam) on every lit
// twilight column, and the top cell absorbed as if the star were overhead-ish.
template <typename XF, typename KR, typename RtFT, int NCT>
KOKKOS_INLINE_FUNCTION
void BeamTopTau(const Real rcut, const Real sinz, const int icut, const int ie,
                const int m, const XF &X1F, const KR &Krs, Real (&tauh)[NCT],
                RtFT (&thi)[NCT]) {
  const Real bt = X1F(m,ie+1)*sinz;
  if (bt <= rcut) {
    for (int cc=0; cc<NCT; ++cc) {
      tauh[cc] = 1.0e30;
      thi[cc] = static_cast<RtFT>(0.0);
    }
    return;
  }
  const Real b2 = bt*bt;
  int jlo = icut;
  for (int jj=ie; jj>=icut; --jj) {
    if (X1F(m,jj) <= bt) {
      jlo = jj;
      break;
    }
  }
  const Real rl = X1F(m,jlo);
  Real prev = (rl*rl > b2) ? sqrt(rl*rl - b2) : 0.0;
  for (int jj=jlo; jj<ie+1; ++jj) {
    const Real ru = X1F(m,jj+1);
    const Real cur = sqrt(ru*ru - b2);
    for (int cc=0; cc<NCT; ++cc) {
      tauh[cc] += 2.0*(cur - prev)*static_cast<Real>(Krs[cc][jj]);
    }
    prev = cur;
  }
  for (int cc=0; cc<NCT; ++cc) thi[cc] = RT_EXP(-static_cast<RtFT>(tauh[cc]));
}
// problem/ck_sweep_cache: CACHE THE PER-(CELL, CHAIN) LAYER OPERATOR IN THE THREAD.
//
// Every half layer this kernel ever steps across is the half of ONE cell taken at ONE
// chain's mu -- the note on ck_impl_frozen_op in the chain below says so, and the
// implicit path already exploits it through the GLOBAL ckc0_g/ckci_g/ckco_g triple.  The
// explicit path did not, and paid for it:
//
//   per (cell, chain)          plane-parallel        ck_spherical
//   ck_kappa look-ups                2                    4
//   expm1 + 3 divides                4                    8
//
// because each of the four column passes (the down probe P1, the up probe P2, the real
// down-sweep and the real up-sweep) rebuilt the whole operator from the k-table, and
// because within one pass the SAME half layer is built twice -- once as the lower half
// of cell i in iteration i and once as the upper half of cell i in iteration i-1.
//
// The switch has three settings, because the two kinds of reuse cost different things:
//
//   0  off -- the code as it was.
//   1  CARRY.  Within one pass the second crossing of a half layer reads the triple the
//      first left in NC registers (cry0/cryi/cryo, exactly like kfar), so each pass
//      builds it once per cell instead of twice: 8 -> 4 expm1 under ck_spherical, 4 -> 2
//      plane-parallel.  Costs NO private memory.
//   2  COLUMN.  The first pass that visits a cell (P1 under ck_spherical, the down-sweep
//      otherwise) also parks kappa and the triple in NC x NN private columns and the
//      other three passes load them: 1 look-up and 1 expm1 per (cell, chain).  Costs
//      4 x NC x NN private storage (17.4 kB/thread at the production tier), which is
//      real -- see tests_ck_sweep_cost/README.md, where it wins 7 % with the spherical
//      form on and loses 1 % with it off.
//
// The stored numbers are the very ones the recomputation would have produced -- same
// operands, same expression, same order -- so the sweep is BITWISE identical at every
// setting, which is the A/B that gates it.
inline int ck_sweep_cache = 2;
// problem/ck_sweep_form: SOLVE THE FACE COUPLING EXACTLY, IN TWO PASSES, NO PROBE.
//
// WHAT THE PROBE WAS FOR.  Under ck_spherical the up and down rays are coupled at every
// face by c = beta (d_above - u_below): the down-sweep needs u_below, which it does not
// have, so the kernel runs two PROBE passes (a plane-parallel down probe, then an up
// probe carrying the mixing) to estimate it, and only then the real down and up sweeps.
// Four column passes per chain instead of two, and a c that is right only to O(beta^2).
//
// WHAT THIS DOES INSTEAD.  The column is a stack of half layers that TRANSMIT with no
// reflection (pure absorption) separated by faces that REFLECT with coefficient beta --
// that is all the mixing is.  So the exact solution is the textbook invariant-imbedding
// (adding) one: carry, from the cut face upward, the linear relation between the two
// rays at each face, then impose the top boundary datum and come back down.  TWO passes,
// EXACT, and the probe disappears.  Both settings do exactly that; they differ only in
// which variable carries the relation, which is a real arithmetic difference at the
// faces and in the layers:
//
//   1  tm -- the TRANSFER-MATRIX / adding form, in (u, d).  The relation is
//         u_below = R d_below + Sc,  R = the reflectance of everything below the face,
//         Sc = what it emits.  A half layer of transmission T maps
//              R <- T^2 R,   Sc <- T (R Q + Sc) + P,
//         P and Q being the layer's emission along the up and down crossings -- NO
//         DIVIDE, and no raw product of exponentials anywhere, which is the whole point
//         of using reflectances rather than 2x2 transfer matrices: T^2 R underflows to
//         zero harmlessly where a transfer matrix would overflow on 1/T.  The face is
//         the Moebius map of the reflection addition,
//              R <- (R + beta)/(1 + R beta),   Sc <- (1 - beta) Sc/(1 + R beta),
//         one divide per face.
//
//   2  sd -- the same problem in the DECOUPLED variables S = (u+d)/2, D = (u-d)/2, where
//         the spherical face conditions are "S passes unchanged, A D passes unchanged"
//         and therefore carry NO cross-coupling at all.  The relation is written
//              D = -G S + H,    G = (1-R)/(1+R),  H = Sc/(1+R),
//         and the face is then just a RESCALING, G <- rho G, H <- rho H with
//         rho = A_below/A_above = (1-beta)/(1+beta) -- two multiplies, no divide, which
//         is exactly the statement that (S, D) diagonalise the face.  The cost moves
//         into the layer, which needs one divide per half layer:
//              den = (1+G) + T^2 (1-G),
//              G <- [(1-T^2) + G(1+T^2)]/den,
//              H <- [T((1-G) Q + 2H) + P(1+G)]/den.
//         Two half layers per cell against one face per cell, so sd trades one divide
//         for two: it is the more elegant form and the more expensive one.
//
// A NOTE ON "SWEEPING IN S AND D".  One cannot sweep S and D as an initial-value problem:
// inside a layer they obey S'' = S - B, whose solutions are e^{+t} as well as e^{-t}, so
// a forward recurrence in (S, D) amplifies like e^{+tau} and dies on the first thick
// cell.  u and d ARE the characteristics; what (S, D) buys is the face, not the layer.
// Hence the Riccati variable above, which is the stable way to carry the same statement.
//
// THE BOUNDARY DATA.  At the cut face the up intensity is given (B_cut + the internal
// flux), i.e. R = 0, Sc = B_cut + Iint -- in (S, D), G = 1 and H = the same datum, the
// relation S + D = u_bc.  At the top face the DOWN intensity is given (the unresolved
// ghost column's emission, in the FACE's own frame), i.e. S - D = Itop above the face;
// pass 2 starts by combining it with the accumulated relation, which is one divide.
//
// WHAT IS CONSERVED, AND HOW.  Exactly as in the four-pass form, A D is made continuous
// by CONSTRUCTION rather than inferred: at each face the flux is reported as
// A_below (u_b - d_b)/A_face and the up ray is handed u_above = d_above + (u_b - d_b)
// A_below/A_above.  The deposit then telescopes to -(A_t F_t - A_b F_b)/V cell by cell
// identically -- and here it does so with NOTHING left over, because the pass-2 up ray is
// closed on the very u_below that the face above it reported (see the one-line correction
// in the sweep).  Fluxes stay per unit area at their own face.
//
// WHAT DIFFERS FROM THE FOUR-PASS FORM, AND IT IS NOT ROUND-OFF.  The probe form's c is
// accurate to O(beta^2) only; these forms are exact.  So the two answers differ by the
// probe's lag, not by rounding -- measured in tests_ck_sweep_form/README.md.
//
// NOT CONVERTED, and refused at startup: ck_spherical = false (with no face coupling
// there is nothing to solve and the plane-parallel sweep is already the exact two-pass
// form), rt_layer_legacy, and ck_implicit with sd (2).  ck_implicit with tm (1) carries
// its own tridiagonal, the JAC instantiation of the kernel (2026-09-22, phase T1 of
// tests_ck_implicit/DESIGN_tm.md); ck_impl_frozen_op is instantiated for it (phase T2),
// and ck_impl_lin re-applies its stored factorisation in a separate kernel (phase T3).
//
// THE DEFAULT.  Since 2026-09-22 deep_hot_jupiter_rt defaults this to 1 (tm) wherever it
// can run, i.e. under ck_spherical, and to 0 in exactly the cases refused above, so that
// no working input is turned into a startup fatal by the flip.  The reason is that the
// four-pass form's probe-lagged mixing constant is only O(beta^2) -- 2.3e-2 of the face
// flux at the production area ratio 3.45 -- while tm and sd solve the coupled column
// exactly (they agree to 3.5e-17), and tm is 1.36x faster on the GPU.  The value below
// remains the flag's own zero-initialisation, i.e. the historical form: it is what any
// other caller gets, and what the pgen falls back to.  See
// tests_ck_sweep_form/README.md and README_default.md.
inline int ck_sweep_form = 0;
// problem/ck_sph_dilute (PROTOTYPE, default false = bitwise the face above): replace the
// spherical face condition of the tm form (S and A D continuous, i.e. the Eddington
// K = J/3 closure, under which J does NOT dilute through a transparent shell and a sharp
// photosphere below a shell of area ratio x = A_top/A_ph emits 2x/(1+x) times sigma T^4)
// by the geometric one for hemispherically isotropic streams:
//     d_below = d_above,   A_a u_above = A_b u_below + (A_a - A_b) d_above,
// i.e. the down ray's intensity is invariant (the rays that reach the sphere below keep
// it) and the down rays that miss the sphere below turn round into the up ray.  A (u - d)
// is still single-valued at every face, an isotropic field (u = d) passes unchanged, and
// through a transparent shell A u is conserved, so the photosphere emits sigma T^4 for
// any x.  In the thick limit it differs from the Eddington face by O(1/(tau r/H)).
// Face maps (tm): R <- 1 - rho (1 - R), Sc <- rho Sc, rho = A_b/A_a = (1-beta)/(1+beta);
// solve: u_b = R d_a + Sc, d_b = d_a.  tm (ck_sweep_form = 1) only, explicit or
// ck_implicit, NOT with ck_impl_lin / ck_impl_jac_lin (their stored-factorisation kernels
// carry the old face).  tests_ck_sph/closure_0930 (cksph_test_0930).
inline bool ck_sph_dilute = false;
// problem/ck_sph_face (default 6 on a fresh spherical ck run whose keys allow it, else 0
// = bitwise the Eddington face above; f6-default 10-02, see the pgen): the face
// condition of the tm spherical ck two-stream, one of a one-parameter family that keeps
// A (u - d) single-valued (so the deposit still telescopes):
//     e (u_a - u_b) + (d_a - d_b) = 0,   A_a (u_a - d_a) = A_b (u_b - d_b),
// e = 1 the Eddington face (S continuous), e = 0 the dilute face of ck_sph_dilute.
//   0  Eddington (e = 1), the historical face;
//   1  dilute (e = 0), = ck_sph_dilute;
//   2  thick/thin blend: e = 1 - exp(-dtau_f/mu) per chain, dtau_f the radial optical
//      depth of the two half cells that meet at the face.  Thin faces dilute (the
//      transparent shell above the photosphere), thick faces keep the Eddington face
//      that matches the diffusion flux.
//   3  dilute face + CONE PATH FACTOR for the up ray;  4  blend face + cone path.
//      Through a transparent shell above a photosphere of area A_ph the real up field
//      fills a cone, sin^2 th_c = A_ph/A, and absorbs 2/(1 + cos th_c) of its flux per
//      unit radial depth, not the 2 of a hemispherically isotropic stream.  The up
//      crossing of a half cell therefore absorbs along dtau/(mu (1 + cos th_c)); its
//      emission and the down crossing are unchanged, so an isotropic source still
//      emits isotropically.  cos^2 th_c = T^4 R, R the relation's reflection just
//      above the face and T the half cell's plain transmission: R = 1 - A_ph/A through
//      a transparent shell over a black photosphere (the dilute face maps 1 - R by
//      A_b/A_a), T^4 R -> 0 in a thick cell, where f -> 1.  Depends on kappa only, so
//      the sweep stays affine in B.  Not meant for ck_dif_dtau handover columns (they
//      start from R = 1).
// General maps (rho = A_b/A_a, Dc = 1 + e R + e rho (1 - R), al = (1 + e)/Dc,
// be = e (1 - rho)/Dc):  R <- 1 - rho (1 - R) al,  Sc <- rho al Sc;  solve
// d_b = al d_a - be Sc_b, u_b = R d_b + Sc_b.  tm (ck_sweep_form = 1) only, explicit or
// ck_implicit; refused with ck_impl_lin / ck_impl_jac_lin.  cksph_test_0930/VARIANTS.md.
inline int ck_sph_face = 0;
// problem/ck_sph_top (default 1 with ck_sph_face = 6, else 0 = bitwise the historical
// datum): the down
// intensity the correlated-k chain kernel hands the top face, (1 - e^-dtau/mu) B over
// the hydrostatic column p/g above the domain.  0 takes kappa and B at the GHOST cell
// (which a pgen IC can leave hot: T_ghost 3415 K on the old dhj grid); 1 sends
// nothing down (d_top = 0); 2 takes kappa and B at the top ACTIVE cell (the column
// depth still from the ghost pressure).  The stellar beam is unchanged.  Chain kernel
// only; refused with ck_impl_lin / ck_impl_jac_lin.
inline int ck_sph_top = 0;
// problem/ck_wall_flux_exact (wallfix-1002; default true on a fresh run, false on a
// restart whose input lacks the key): the internal flux of the ck_sph_face = 5 / 6 wall.
// The Marshak-type wall J + beta H = (B_cut + I_int)(1/2 + beta/4) with the pair at
// mu = mu_eff = sqrt f gives u = R0 d + (B_cut + I_int)(1 + beta/2)/(1 + beta mu), so
// the I_int part of the net flux 2 pi mu (u - d) is 2 pi mu (1 + beta/2)/(1 + beta mu)
// I_int: 4 pi/(2 + sqrt 3) I_int = 1.0718 sigma T_int^4 for the deep field (beta = 2,
// f = 1/3), against pi I_int = sigma T_int^4 at face 0 (f6default_1002 part B: measured
// Etot_bot face 5/6 over face 0 1.07166).  On, I_int in the wall source is scaled by
//   face 5 (chain kernel, rt_chain_ck_lin1p):  c = (1 + beta mu)/(2 mu (1 + beta/2)),
//     exact for the actual beta and mu of the wall when nothing comes back from above
//     (the condition under which face 0's I_int is exactly pi I_int too);
//   face 6 (the ladder, rt_chain_ck_lin6):  c6 = beta (y_lo + y_w)/((2 + beta) y_lo),
//     y_lo the driving-point admittance of the network above the wall node and y_w the
//     wall's conductance: by the Thevenin split the I_int part of the wall current is
//     r^2 c6 I_int (1/2 + beta/4)/beta y_lo/(y_lo + y_w) = r^2 I_int/4 EXACTLY, i.e.
//     the wall flux carries pi I_int whatever the column; for a thick wall cell y_lo =
//     beta mu y_w and c6 = c.
// I_int does not enter any Jacobian row (the rows are dSrc/dB), so only the sources
// change; the formal solution's bottom boundary stays I^+ = B_cut + I_int (isotropic,
// as the exact reference).  false = bitwise the historical wall.
inline bool ck_wall_flux_exact = false;
// problem/ck_sph_face = 5 (PROTOTYPE): the VARIABLE EDDINGTON FACTOR form (Auer 1971;
// Mihalas, Stellar Atmospheres, ch. 7): the moment equations
//     d(A H)/dr = A kappa (B - J),   d(f q r^2 J)/dr = -q r^2 kappa H,
//     ln(q r^2) = int (3f - 1)/(r f) dr   (Phi = q r^2),
// with f = K/J from a FORMAL SOLUTION along impact-parameter (p-z) rays through the
// radial column (straight chords through the shells, S = B, the wall as the bottom
// datum, vacuum at the top), per (band, g), Newton-lagged: recomputed on the first RT
// call of every ck_vef_every-th cycle and frozen in between.  In the tm sweep this is
// the same two-stream with BOTH chains of a (band, g) pair at mu_eff = sqrt(f) inside
// the cell (weights 2 pi w mu_eff, so the pair is exactly the moment system with K =
// f J; s = sqrt(3 f) = mu_eff sqrt 3 below) and the exact jump
// conditions of the moment equations at the faces between frames,
//     Phi f J continuous,   A s (u - d) continuous  (A H),
// i.e. gamma (u_a + d_a) = u_b + d_b, gamma = (Phi f)_a/(Phi f)_b, and
// u_a - d_a = rho' (u_b - d_b), rho' = rho s_b/s_a.  The outer boundary is Auer's
// H = h J: the down ray at the top is reflected, d = zeta u + datum, with
// zeta = (mu_eff - h)/(mu_eff + h).  f = 1/3 (isotropic, diffusion) gives s = 1,
// gamma = 1: the Eddington face (with the pair at 1/sqrt 3 instead of the Gauss
// points, which carries the same diffusion flux).  A (u - d) s is single-valued at
// every face, so the deposit telescopes.  The quadrature in mu is normalised so that
// an isotropic field gives f = 1/3 and a hemispherically isotropic one h = 1/2 to
// round-off.
// Under ck_impl_frozen_op the factors are part of the STORED OPERATOR (the frozen triple
// at mu_eff, the tm factorisation of ck_impl_lin, the ck_impl_jac_lin coefficients), so
// with problem/ck_vef_sync (default true) they are refreshed on a STORING pass only (the
// first pass of a call that does not re-apply an operator stored earlier, see
// ck_impl_xstep): every pass of every call then sees one consistent operator.  false
// is the prototype's timing (the first pass of the cycle, a re-applying one included).
// With ck_impl_lin (production path) face 5 needs ck_impl_lin_thr = 1 and ck_nquad = 2
// (ck_coef then stores the triple at mu_eff for a frozen chain kernel that reads it).
inline int ck_vef_every = 1;
inline int ck_vef_ncore = 8;          // rays through the bottom wall (p < r_cut)
inline bool ck_vef_sync = true;
// problem/ck_vef_report = N > 0 (diagnostic, print only): at every VEF refresh (after
// the first) rank 0 prints, per pressure decade, the range of f = s^2/3 and of the face
// factor gamma over all (band, g, cell), the largest relative change of s and gamma
// against the factors being replaced (the lag of ck_vef_every), and the number of faces
// whose x = gamma rho' (see the kernel) is outside [1/1.1, 1.1].  Nothing reads it.
inline int ck_vef_report = 0;
// problem/ck_f6_report = N > 0 (diagnostic, print only, read only when named): after
// every ck_sph_face = 6 build rank 0 prints the number of ladder elements (g_s, g_p,
// K_eq/B, R_c over all cells and (band, g)), how many are negative or not finite, and
// the smallest ratio row excess / diagonal of the node rows in K (the shunt g_p over
// g_p + the node's two series conductances; > 0: strictly diagonally dominant).
inline int ck_f6_report = 0;
// problem/ck_vef_xlim = X > 1 (default 0 = off, bitwise): limit the VEF face factor
// x = gamma rho' of every interior face to [1/X, X] (gamma = (Phi f)_a/(Phi f)_b,
// rho' = (A_b/A_a) s_b/s_a).  The face map of ck_sph_face = 5 has the reflections
// (1 - x)/(1 + x) and (x - 1)/(1 + x): one of them is NEGATIVE whenever x != 1, and a
// strong f contrast between two cells makes the frozen operator non-monotone (a cell's
// net heating can grow with its own temperature).  On the WASP-121b 10x fresh start the
// factors reach f ~ 0.1 and x = 0.5..1.6 at 1e-5..1e-3 bar within 150 cycles, the
// runaway spreads down and the hydro dt collapses (ckf5_stab_1001/RESULTS.md).
// Limiting x bounds the negative reflection by (X - 1)/(X + 1); it changes gamma only
// where f jumps by a factor ~X^2 between neighbouring cells.  rho' (the flux
// continuity) is untouched, so the deposit still telescopes.
inline Real ck_vef_xlim = 0.0;
// problem/ck_sph_face = 6 (the DEFAULT since f6-default 10-02): the VEF closure of face
// 5 (same formal solution,
// same f, h, beta, same ck_vef_every refresh) discretised as a PASSIVE LADDER NETWORK in
// K = Phi f J (Phi = q r^2), so that it is monotone by construction and needs no limiter
// (ckf6_mono_1001/RESULTS.md).  The moment equations
//     dK/dr = -kappa rho (Phi/r^2) G,   dG/dr = kappa rho r^2 (B - J),   G = r^2 H,
// are solved EXACTLY inside each cell for constant coefficients and the cell's own B
// (the face-5 two-stream at mu_eff = sqrt f is this solution): with
//     a = int kappa rho Phi/r^2 dr,  b = int kappa rho r^2/(Phi f) dr  (over the cell,
//     ln Phi = (3 - 1/f) ln r inside it), Z = sqrt(a/b), Lam = sqrt(a b),
// the cell is the two-port of a uniform lossy line: a series conductance
// g_s = 1/(Z sinh Lam) between its two face nodes and a shunt g_p = tanh(Lam/2)/Z from
// each face node to K_eq = kappa rho V B/b (its emission).  K is continuous at a face
// (the VEF jump condition Phi f J continuous) and so is G.  Face 5 is exactly this
// network, but closes the thick-cell staircase (two thick cells at different constant B
// exchange Delta K/(Z_a + Z_b) instead of the diffusion flux) with a LINEAR source inside
// the cell, B at the faces interpolated from the neighbours: that makes a cell's
// emission depend on its neighbours' temperatures and the operator non-monotone (the
// 10x blow-up, ckf5_stab_1001).  Face 6 keeps the emission cell-local and closes the
// staircase with a CONTACT resistance at every interior face,
//     R_c = w(Lam_b) Z_a phc(Lam_a) + w(Lam_a) Z_b phc(Lam_b),
//     phc(Lam) = Lam/2 - tanh(Lam/2),   w(Lam) = Lam^2/(1 + Lam^2),
// which is the exact value for which a uniform medium with a LINEAR B carries the
// Eddington diffusion flux at any Lam (Lam^3/24 for a thin cell, Z (Lam/2 - 1) per side
// for a thick one), switched off next to an optically thin neighbour so that a thick
// cell under a thin one keeps its surface emission (an unresolved photosphere).
// Every element is a positive conductance and every source is the cell's own B: the
// network is passive, so J >= 0, J is monotone in every B, and the heating of a cell
// does not grow with its own B (the Newton row b = 1 - h dSrc_i/dT_i/cv > 0; only the
// cut cell, whose B also drives the wall, can reach + round-off).  Top: Auer's H = h J
// (h from the formal solution, vacuum above: needs ck_sph_top = 1); wall: face 5's
// Marshak-type J + beta H = (B_cut + I_int)(1/2 + beta/4).  Solved by the positive
// Norton recursion from the top (Y' = g Y/(g + Y), I' = g I/(g + Y): no subtraction) and
// the node potentials upward; Src_i = 4 pi g_w g_p (K_lo + K_hi - 2 K_eq)/V (the cell's
// absorption minus emission), the fluxes accumulated from the top: V Src = (A F)_i -
// (A F)_{i+1} to round-off.  The Newton tridiagonal is the exact (i-1, i, i+1) part of
// dSrc/dB of the network (driving-point and transfer impedances of the ladder).
// Production path only (ck_implicit, frozen operator, ck_impl_lin + ck_impl_jac_lin,
// ck_impl_lin_thr = 1, ck_nquad = 2): ck_coef stores kappa rho, ck_f6_build the ladder,
// rt_chain_ck_lin6 re-applies it, rt_chain_ck_jlin6 forms the tridiagonal.  A storing
// pass that runs the chain kernel (no ck_beam_par) runs it on the face-5 map and then
// takes Src, Fb and Em from rt_chain_ck_lin6.
//! \fn Real CkF6Ex
//! \brief (e^x - 1)/x, 1 at x = 0
KOKKOS_INLINE_FUNCTION
Real CkF6Ex(const Real x) {
  return (fabs(x) < 1.0e-5) ? 1.0 + x*(0.5 + x/6.0) : expm1(x)/x;
}
//! \struct CkF6Cell
//! \brief ck_sph_face = 6: one cell's two-port (see ck_sph_face = 6)
struct CkF6Cell {
  Real lpl;    // ln Phi at the lower face (Phi = 1 at the top face of the column)
  Real z;      // Z
  Real lam;    // Lam
  Real gp;     // shunt conductance tanh(Lam/2)/Z
  Real gs;     // series conductance 1/(Z sinh Lam) (0 beyond Lam = 700)
  Real keq;    // K_eq/B = kappa rho V/b
  Real phc;    // Lam/2 - tanh(Lam/2)
};
//! \fn CkF6Cell CkF6Make
//! \brief the two-port of a cell between radii rl < rh with kappa rho kr, Eddington
//! factor f and ln Phi lph at its upper face
KOKKOS_INLINE_FUNCTION
CkF6Cell CkF6Make(const Real rl, const Real rh, const Real kr, const Real f,
                  const Real lph) {
  CkF6Cell c;
  const Real e = 3.0 - 1.0/f;
  const Real t = log(rh/rl);
  Real lpl = lph - e*t;
  lpl = (lpl > 600.0) ? 600.0 : ((lpl < -600.0) ? -600.0 : lpl);
  c.lpl = lpl;
  const Real phl = exp(lpl);
  const Real a = kr*(t/rl)*CkF6Ex((e - 1.0)*t)*phl;
  const Real b = (kr/f)*rl*rl*rl*t*CkF6Ex(t/f)/phl;
  c.z = sqrt(a/b);
  c.lam = sqrt(a*b);
  const Real lm = (c.lam < 700.0) ? c.lam : 700.0;
  const Real th = tanh(0.5*lm);
  c.gp = th/c.z;
  c.gs = (c.lam < 700.0) ? 1.0/(c.z*sinh(lm)) : 0.0;
  c.keq = kr*(rh*rh*rh - rl*rl*rl)/(3.0*b);
  c.phc = (c.lam < 1.0e-3) ? c.lam*c.lam*c.lam/24.0 : 0.5*c.lam - th;
  return c;
}
//! \fn Real CkF6Rc
//! \brief the contact resistance of the face between cells a (below) and b (above)
KOKKOS_INLINE_FUNCTION
Real CkF6Rc(const CkF6Cell &a, const CkF6Cell &b) {
  const Real la2 = a.lam*a.lam, lb2 = b.lam*b.lam;
  return lb2/(1.0 + lb2)*a.z*a.phc + la2/(1.0 + la2)*b.z*b.phc;
}
template <typename XF>
inline void CkVefReport(Mesh *pm, const DvceArray5D<Real> &vs,
                        const DvceArray5D<Real> &vg,
                        const DvceArray5D<Real> &vso, const DvceArray5D<Real> &vgo,
                        const DvceArray4D<Real> &pb, const DvceArray3D<int> &icc,
                        const XF &X1F,
                        const int nmb1, const int nc, const int ks, const int ke,
                        const int js, const int je, const int ie) {
  constexpr int NBIN = 12;   // decades 1e-9 .. 1e3 bar
  constexpr int NQ = 9;      // fmin fmax gmin gmax dsmax dgmax nout xmin xmax
  DvceArray1D<Real> acc("ck_vef_rep", NBIN*NQ);
  auto acc_ = acc;
  Kokkos::parallel_for("ck_vef_rep_init", Kokkos::RangePolicy<DevExeSpace>(0, NBIN),
  KOKKOS_LAMBDA(const int b) {
    acc_(b*NQ+0) = 1.0e30; acc_(b*NQ+1) = -1.0e30;
    acc_(b*NQ+2) = 1.0e30; acc_(b*NQ+3) = -1.0e30;
    acc_(b*NQ+4) = 0.0; acc_(b*NQ+5) = 0.0; acc_(b*NQ+6) = 0.0;
    acc_(b*NQ+7) = 1.0e30; acc_(b*NQ+8) = -1.0e30;
  });
  par_for("ck_vef_rep", DevExeSpace(), 0, nmb1, 0, nc - 1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int c, const int k, const int j) {
    const int icut = icc(m,k,j);
    for (int i=icut; i<=ie; ++i) {
      const Real p = pb(m,k,j,i);
      int b = static_cast<int>(floor(log10(fmax(p, 1.0e-30)))) + 9;
      b = (b < 0) ? 0 : ((b >= NBIN) ? NBIN - 1 : b);
      const Real s = vs(m,c,i+1,k,j);
      const Real f = s*s/3.0;
      const Real g = vg(m,c,i,k,j);
      Kokkos::atomic_min(&acc_(b*NQ+0), f);
      Kokkos::atomic_max(&acc_(b*NQ+1), f);
      Kokkos::atomic_min(&acc_(b*NQ+2), g);
      Kokkos::atomic_max(&acc_(b*NQ+3), g);
      const Real so = vso(m,c,i+1,k,j), go = vgo(m,c,i,k,j);
      if (so > 0.0) Kokkos::atomic_max(&acc_(b*NQ+4), fabs(s/so - 1.0));
      if (go > 0.0) Kokkos::atomic_max(&acc_(b*NQ+5), fabs(g/go - 1.0));
      // x = gamma rho' of the face below cell i (rho' = (A_b/A_a) s_b/s_a, A = the
      // frame area r_node^2): the face reflections are (1 - x)/(1 + x) and
      // (x - 1)/(1 + x), both >= 0 only at x = 1
      if (i > icut) {
        const Real rm = X1F(m,i-1), r0 = X1F(m,i), rp = X1F(m,i+1);
        const Real ab = rm*rm + rm*r0 + r0*r0, aa = r0*r0 + r0*rp + rp*rp;
        const Real x = g*(ab/aa)*vs(m,c,i,k,j)/s;
        Kokkos::atomic_min(&acc_(b*NQ+7), x);
        Kokkos::atomic_max(&acc_(b*NQ+8), x);
        if (x > 1.1 || x < 1.0/1.1) Kokkos::atomic_add(&acc_(b*NQ+6), 1.0);
      }
    }
  });
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), acc);
  Real v[NBIN*NQ], r[NBIN*NQ];
  for (int q=0; q<NBIN*NQ; ++q) v[q] = h(q);
#if MPI_PARALLEL_ENABLED
  for (int b=0; b<NBIN; ++b) {
    v[b*NQ+0] = -v[b*NQ+0]; v[b*NQ+2] = -v[b*NQ+2]; v[b*NQ+7] = -v[b*NQ+7];
  }
  Real cnt[NBIN];
  for (int b=0; b<NBIN; ++b) {
    cnt[b] = v[b*NQ+6]; v[b*NQ+6] = 0.0;
  }
  MPI_Reduce(v, r, NBIN*NQ, MPI_ATHENA_REAL, MPI_MAX, 0, MPI_COMM_WORLD);
  Real cr[NBIN];
  MPI_Reduce(cnt, cr, NBIN, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  for (int b=0; b<NBIN; ++b) {
    r[b*NQ+0] = -r[b*NQ+0]; r[b*NQ+2] = -r[b*NQ+2]; r[b*NQ+6] = cr[b];
    r[b*NQ+7] = -r[b*NQ+7];
  }
#else
  for (int q=0; q<NBIN*NQ; ++q) r[q] = v[q];
#endif
  if (global_variable::my_rank == 0) {
    for (int b=0; b<NBIN; ++b) {
      if (r[b*NQ+1] < 0.0) continue;
      std::printf("### vefrep ncycle=%d p=1e%+d f=[%.4e,%.4e] gam=[%.4e,%.4e] "
                  "dsmax=%.3e dgmax=%.3e x=[%.4e,%.4e] nx10=%.0f\n",
                  static_cast<int>(pm->ncycle),
                  b - 9, r[b*NQ+0], r[b*NQ+1], r[b*NQ+2], r[b*NQ+3], r[b*NQ+4],
                  r[b*NQ+5], r[b*NQ+7], r[b*NQ+8], r[b*NQ+6]);
    }
  }
}
inline bool CkVefDue(Mesh *pm, const bool ckfop, const bool ckfst) {
  const bool stp = !ck_vef_sync || !ckfop || ckfst;    // a pass that may refresh
  return stp && (ck_vef_last < 0 || ck_vef_s_ptr == nullptr ||
                 pm->ncycle - ck_vef_last >= ck_vef_every);
}
// problem/rt_top_re: what the unresolved column ABOVE the domain sends back down.
// false (historical) makes it radiate at the ghost cell's own temperature. That is safe
// only while the ghost is pinned to something outside the solution: with an open outer
// boundary the ghost IS the top active cell, the column hands back exactly what that cell
// emitted, and the net flux through the top face goes to zero as the column thickens --
// the atmosphere isothermalises and the object stops radiating (measured on the red
// giant: emergent flux 0.3 % of L, top density 22x). true treats the column as a slab in
// RADIATIVE EQUILIBRIUM instead: it absorbs (1-e^-dtau) I_up from below and re-emits half
// up to space and half back down, so its source is I_up/2 whatever its optical depth and
// the face always keeps at least half of its outgoing flux. Grey path only so far.
// NOTE, and the reason for the startup warning below: the back-radiation is delivered
// as (1 - exp(-dtau/mu)) bsrc with dtau = kappa(ghost) p/g.  Above rad_kappa_rmax the
// ghost's opacity is rad_kappa_above, so with rad_kappa_above = 0 -- the exact setting
// every corona run uses -- dtau is 0 and this boundary hands back EXACTLY ZERO however
// carefully bsrc was probed.  Measured in V8_corona: I_dn at the top face of the active
// column is 0.000000e+00 at every cycle.  Behaviour unchanged; it is only documented.
inline bool rt_top_re = false;
// SELF-LUMINOUS objects: > 0 uses this internal temperature directly instead of the
// Thorngren+2019 T_int(T_eq) relation, which is a fit for IRRADIATED giant planets and
// floors at 100 K.  A star or a brown dwarf sets it to its own T_eff (and Teq = 0, which
// zeroes the stellar sweep).  0 keeps the historical behaviour exactly.
inline Real rt_tint_override = 0.0;
// problem/rt_apply_debug: print the energy balance the applied source actually sees, for
// the top rt_apply_debug_n cells of the (rt_dump_m, rt_dump_k, rt_dump_j) column, on the
// first N calls.  The column dump says what the SWEEP produced; this says what the CELL
// then did with it, which is where a thin top layer that will not sit at its equilibrium
// temperature has to be caught.
inline int rt_apply_debug = 0;
inline int rt_apply_debug_n = 8;
// problem/rt_cell_report: the radiative-balance report for the cells that decide the
// top of the active column.  Prints, once per RT call, the FIRST cell that trips the
// Newton positivity rescue -- its location, state, opacity, tau to the top, the up and
// down streams at its two faces, the absorption and emission terms, e_eq, dt and every
// Newton iterate -- and, every rt_report_every cycles, the same for the fixed cell that
// contains rt_report_r on the (ks, js) column.  Off by default: a default run neither
// allocates the two stream arrays nor prints.
inline bool rt_cell_report = false;
inline Real rt_report_r = 3.887e12;               // problem/rt_report_r [cm]
inline int rt_report_every = 100;                 // problem/rt_report_every [cycles]
inline DvceArray4D<Real> *rt_idn_ptr = nullptr;   // downward stream at face i
inline DvceArray4D<Real> *rt_iup_ptr = nullptr;   // upward stream at face i
inline bool rt_srclim_warned = false;             // the one-time warning has been issued
inline int64_t rt_srclim_acc = 0;     // clips on THIS rank not yet reduced (see below)
inline int rt_srclim_ncall = 0;       // outer RT calls since the last reduction
// problem/nan_report: catch the cell whose conserved energy the GREY apply makes
// non-finite or non-positive, IN the kernel, with the inputs that produced it.  Off by
// default, so a default run is bit-identical.  Pointers, not Views, for the same reason
// the pgen's guards use pointers: a file-scope View outlives Kokkos::finalize.
inline bool rt_nan_report = false;
inline DvceArray1D<int> *rt_nanrep_cnt = nullptr;
inline DvceArray1D<Real> *rt_nanrep_rec = nullptr;
inline int rt_nanrep_lines = 0;                   // printed so far on this rank
inline int rt_nanrep_maxlines = 400;              // then stay quiet (the run continues)

// problem/rt_use_cons NaN GUARD.  EintFromCons returns whatever E - KE gives, INCLUDING
// a negative number: in a runaway cell the total energy is kinetic dominated by many
// orders of magnitude and the difference is catastrophic cancellation.  A non-positive
// e makes TempKelvin/PresTempFromEint return NaN, and one NaN Planck function poisons
// the whole column (measured: T(i-1) = nan, B(i-1) = nan, I_dn = nan ahead of the
// red-giant dt collapse).  eiN clamps such a cell to e(rho, tfloor) -- the same floor
// ConsToPrim would apply -- so the sweep never sees a NaN, and counts the clamps here.
// Active only under rt_use_cons; with it off eiN returns w0 and nothing changes.
inline DvceArray1D<int> *rt_eiclamp_cnt = nullptr;
inline bool rt_eiclamp_warned = false;

// The SAME hazard one level up: a (p, T) pair the EOS could not form -- a non-positive
// or non-finite conserved DENSITY in the cell a slot reads, or a solve that did not
// converge -- gives a NaN Planck function and a NaN top-of-column optical depth, and a
// single NaN poisons the WHOLE column below it.  The down-sweep carries the stream as
// I(i) = (1-e0)*I(i+1) + ..., and at large dtau (1-e0) underflows to EXACTLY zero, so
// 0*NaN = NaN survives to arbitrary depth: measured in I8_vceil, 348 rescued cells at
// tau_to_top = 75-90 with I_dn = -nan while their own B, T and kappa were finite.  The
// two guards below are the density and the state analogues of the eiN clamp, counted in
// one counter, and they are no-ops on any state the EOS returned cleanly.
inline DvceArray1D<int> *rt_stclamp_cnt = nullptr;
inline bool rt_stclamp_warned = false;

//----------------------------------------------------------------------------------------
//! \fn bool RTBadState
//! \brief true when a (p, T) pair cannot be used to form a Planck function or a
//! top-of-column optical depth: NaN, infinite, or non-positive.
KOKKOS_INLINE_FUNCTION
bool RTBadState(const Real p, const Real t) {
  return !(t > 0.0) || !(p > 0.0) || !Kokkos::isfinite(t) || !Kokkos::isfinite(p);
}

//----------------------------------------------------------------------------------------
//! \fn Real RTWPrad
//! \brief w(rho,T) a T^4/3 -- the share of the LTE radiation pressure the hydro still
//! carries in a cell, with w the SAME taper weight the EOS uses (rad_taper::WeightGated).
//! Used only by the optical-depth-gated radiative momentum source
//! (problem/rt_force_tau_gate), which differences it across a cell's neighbours.
KOKKOS_INLINE_FUNCTION
Real RTWPrad(const Real dd, const Real tt, const Real xlo, const Real xhi,
             const Real ylo, const Real yhi, const bool tgate, const Real arad) {
  if (!(dd > 0.0) || !(tt > 0.0)) return 0.0;
  Real ww, dwx, dwy;
  rad_taper::WeightGated(log10(dd), xlo, xhi, log10(tt), ylo, yhi, tgate, ww, dwx, dwy);
  return ww*arad*tt*tt*tt*tt/3.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real RTForceTauGate
//! \brief s(tau): 1 at and below tau_lo, 0 at and above tau_hi, a cubic smoothstep in
//! log10(tau) between.  lg_lo/lg_hi are log10 of the two thresholds.
KOKKOS_INLINE_FUNCTION
Real RTForceTauGate(const Real tau, const Real lg_lo, const Real lg_hi) {
  if (!(tau > 0.0)) return 1.0;
  const Real lt = log10(tau);
  if (lt <= lg_lo) return 1.0;
  if (lt >= lg_hi) return 0.0;
  const Real x = (lt - lg_lo)/(lg_hi - lg_lo);
  return 1.0 - x*x*(3.0 - 2.0*x);
}

//----------------------------------------------------------------------------------------
//! \fn Real RTTopDtau
//! \brief the optical depth of the unresolved hydrostatic column above the domain,
//! kappa p / g_eff, with the g_eff = 0 case resolved.
//!
//! g_eff > 0 returns the old expression BIT FOR BIT.  Where the effective gravity is
//! zero (a plane-parallel run that never set one) the old form was kappa*p/0, which is
//! +inf for an opaque ghost -- an opaque lid, harmless -- but 0/0 = NaN as soon as the
//! ghost's opacity is zero, which is exactly what the density gate makes it above the
//! star.  The limit is taken here instead: no opacity above the domain means no column.
KOKKOS_INLINE_FUNCTION
Real RTTopDtau(const Real kap, const Real pcgs, const Real geff) {
  const Real kp = kap*pcgs;
  if (geff > 0.0) return kp/geff;           // the old expression, unchanged
  return (kp > 0.0) ? 1.0e30 : 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RTSourceLimiterWarn
//! \brief COUNT the cells LimitRTSource clipped on this rank; RTSourceLimiterReduce says
//! so, once per run.
//!
//! This is called from inside the pass (picket_fence_two_stream_RT_pass), whose number
//! of calls is RANK-LOCAL under ck_implicit (the Newton loop stops per rank), so it must
//! not hold a collective: it only accumulates.  The count used to be printed from here
//! on rank 0 with rank 0's share alone, so a run whose clipping sat on other ranks
//! showed a clean log.
inline void RTSourceLimiterWarn(const int nclip) {
  if (rt_srclim_warned || nclip <= 0) return;
  rt_srclim_acc += nclip;
}

//----------------------------------------------------------------------------------------
//! \fn void RTSourceLimiterReduce
//! \brief sum the clip count over ALL ranks and say ONCE, from rank 0, that LimitRTSource
//! has clipped cells.
//!
//! Called once per outer call (picket_fence_two_stream_RT), which every rank makes the
//! same number of times, so every rank takes the same branch: the reduction happens on
//! the same calls everywhere and rt_srclim_warned flips on every rank at once (it is set
//! from the reduced sum).  Once, not per call: a run that trips this trips it every
//! cycle, and the message is about the configuration, not the individual step.  After
//! the warning there is no reduction at all; before it one MPI_Allreduce of one integer
//! every kRTSrcLimEvery calls (at most a few tens of microseconds per 16 RT calls).
constexpr int kRTSrcLimEvery = 16;
inline void RTSourceLimiterReduce() {
  if (rt_srclim_warned) return;
  int64_t tot = rt_srclim_acc;
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    if (++rt_srclim_ncall < kRTSrcLimEvery) return;
    rt_srclim_ncall = 0;
    MPI_Allreduce(MPI_IN_PLACE, &tot, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  if (tot <= 0) return;
  rt_srclim_warned = true;
  rt_srclim_acc = 0;
  if (global_variable::my_rank == 0) {
    std::cout << "### WARNING in deep_hot_jupiter_rt: the explicit radiative source was "
              << "clipped in " << tot << " cell update(s) (all ranks";
    if (global_variable::nranks > 1) {
      std::cout << ", summed over the last <= " << kRTSrcLimEvery << " RT calls";
    }
    std::cout << ") by problem/rt_de_max = " << rt_de_max
              << ".\n    The radiative time e/|src| is shorter than the timestep there, "
              << "so the radiation is outside\n    the regime this operator-split scheme "
              << "is valid in. The usual cause is a pressure or\n    density floor "
              << "pinning the top of the atmosphere; see docs/correlated_k_rt.md.\n"
              << "    Reported once per run." << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real LimitRTSource
//! \brief cap one explicit radiative energy update at a fraction of the cell's internal
//! energy.
//!
//! WHY THIS EXISTS. The radiation is operator split and applied EXPLICITLY at the
//! hydrodynamic timestep: u0(IEN) += src*bdt, with src the net flux divergence plus the
//! stellar heating. That is stable only while the local radiative time e/|src| exceeds
//! bdt, and nothing in the code enforces it. Measured on the column that first blew up in
//! the ideal-gas correlated-k run (see docs/correlated_k_rt.md): a healthy column has
//! e/|src| between 86 and 570 timesteps, but once the pressure floor pins the top of the
//! atmosphere at p = pfloor -- which fixes e = pfloor/(gamma-1), 2.1 erg/cm^3 there -- an
//! unremarkable flux divergence of 29 erg/cm^3/s gives e/|src| = 0.014 timesteps. The
//! explicit update then drives e far negative, the floor rescues it, and the next step
//! overshoots harder: adjacent cells at 650, 2081 and 9806 K, and a NaN a few cycles
//! later.
//!
//! WHY A HARD CLAMP, NOT A SMOOTH ONE. A smooth limiter perturbs every cell a little; a
//! clamp is the IDENTITY wherever |de| < de_max*e, so it cannot move an answer in any
//! regime where the explicit step was legitimate in the first place. The correlated-k
//! fluxes are validated against Exo-FMS and that validation has to survive this. At the
//! default de_max = 0.5 the clamp is 40x looser than the worst healthy cell measured
//! above and 140x tighter than the runaway, so it separates the two cleanly.
//!
//! This bounds the damage; it does not make the step accurate. A run that trips it is
//! reporting that its floors, or its timestep, put the radiation outside the regime the
//! scheme is valid in -- which is why tripping it is warned about exactly once.
// The emissivity ratio below which BFace stops trusting the neighbour's Planck function
// as this layer's far source endpoint.  1/10: at V8_corona's join the ratio is exactly 0
// (rad_kappa_above = 0) and in a smooth stellar column it never falls below ~0.5, so this
// separates the two by 5x on the smooth side and infinitely on the pathological one.
#define RT_BFACE_R 0.1

//----------------------------------------------------------------------------------------
//! \fn Real BFace
//! \brief the Planck source at a layer's FAR endpoint, weighted by emitting matter.
//!
//! WHY THIS EXISTS.  Both sweeps take the source function of a layer as B interpolated
//! linearly between the two cell CENTRES -- alp*B_far + bet*B_near -- and Em averages the
//! same two.  Nothing in that asks whether the far cell radiates at all.  Across the join
//! at rad_kappa_rmax it does not: measured in V8_corona at t = 0, the last active cell
//! (i = 439, T = 3363.6 K, B = 2.310e9) sits directly under the first corona cell
//! (T = 6.0e5 K, kappa = 0, B = 2.339e18), so its emission was formed with a Planck
//! function 1.01e9 times its own: 3.60e-8 erg/cm^3/s of emission against an absorption of
//! 5.50e-17 and an internal energy of 7.46e-10, i.e. 1480 times its own energy removed in
//! one 30.6 s step.  All 192 cells of that shell hit the temperature floor on the FIRST
//! RT call and the corona then accreted onto the cold sink.
//!
//! WHAT IT DOES.  It is SURGICAL: while the far cell's emissivity kappa*rho is within a
//! factor RT_BFACE_R of this layer's, the old endpoint B_far is returned unchanged, so
//! every column whose opacity varies smoothly is BIT-IDENTICAL to the pre-fix code
//! (measured on the V7 grey column: all 39 active faces identical).  Below that ratio the
//! endpoint is blended into the emissivity-weighted face source
//!     S = (j_own B_own + j_far B_far)/(j_own + j_far),   j = kappa rho,
//! with a weight that reaches S exactly as the far cell stops radiating -- so an inert
//! neighbour contributes B_own, i.e. the layer emits with its OWN Planck function.  The
//! blend is continuous at the threshold (S is entered with weight 0 there), so a cell
//! cannot jump between the two forms.
KOKKOS_INLINE_FUNCTION
Real BFace(const Real k_own, const Real k_far, const Real b_own, const Real b_far,
           const bool on) {
  if (!on) return b_far;                  // problem/rt_bface off: the pre-fix endpoint
  const Real kt = RT_BFACE_R*k_own;
  if (k_far >= kt) return b_far;          // the old expression, bit for bit
  if (!(k_own > 0.0)) return b_far;       // this layer does not emit at all
  if (!(k_far > 0.0)) return b_own;       // the far cell is radiatively inert
  const Real w = k_far/kt;                // 1 at the threshold, 0 at k_far = 0
  const Real sem = (k_own*b_own + k_far*b_far)/(k_own + k_far);
  return w*b_far + (1.0 - w)*sem;
}

//----------------------------------------------------------------------------------------
//! \fn Real BFaceW
//! \brief d(BFace)/d(b_far): BFace is a CONVEX COMBINATION of b_far and b_own at frozen
//! opacity, so one weight describes it completely and the own-side weight is 1 - this.
//! Used by the correlated-k column solve (problem/ck_implicit) to linearise the sweep.
KOKKOS_INLINE_FUNCTION
Real BFaceW(const Real k_own, const Real k_far, const bool on) {
  if (!on) return 1.0;
  const Real kt = RT_BFACE_R*k_own;
  if (k_far >= kt) return 1.0;
  if (!(k_own > 0.0)) return 1.0;
  if (!(k_far > 0.0)) return 0.0;
  const Real w = k_far/kt;
  return w + (1.0 - w)*k_far/(k_own + k_far);
}

//----------------------------------------------------------------------------------------
//! \fn void CkCinCout
//! \brief ck-next: the half-layer pair (cin, cout) that ck_coef and the chain kernel's
//! cofs store in ck_ci / ck_co, formed again from the stored e0 and kappa rho, the
//! cell's dz and the chain's mu -- the same expressions on the same numbers, so the
//! stored values bitwise.  The linear kernels read two Reals per (cell, chain) less.
KOKKOS_INLINE_FUNCTION
void CkCinCout(const Real kro, const Real dz, const Real mu, const Real e0r, Real &ci,
               Real &co) {
  const Real dtau = 0.5*kro*dz;
  const RtF x = static_cast<RtF>(dtau/mu);
  const RtF e0 = static_cast<RtF>(e0r);
  const RtF one = static_cast<RtF>(1.0);
  const RtF cin = (x > static_cast<RtF>(1.0e-3)) ? (e0 - one + e0/x) : (x/2 - x*x/3);
  const RtF cout = (x > static_cast<RtF>(1.0e-3)) ? (one - e0/x) : (x/2 - x*x/6);
  ci = static_cast<Real>(cin);
  co = static_cast<Real>(cout);
}
//! \fn Real CkIdn
//! \brief ck-next: lP slot 1, 1/(1 + R beta), from slot 0 and beta by ck_lin_build's
//! expression (bitwise the stored value)
KOKKOS_INLINE_FUNCTION
Real CkIdn(const Real rj, const Real bt) {
  const RtF dn = static_cast<RtF>(1.0) + static_cast<RtF>(rj)*static_cast<RtF>(bt);
  return 1.0/static_cast<Real>(dn);
}
//! \fn Real CkEmW
//! \brief ck-next: the emission weight lP slot 5 held, 2 (wfc/mu) kappa rho, from
//! 2 (wfc/mu) (lC row 2, stored by ck_lin_build) as a lone rounded product (CkMulRn):
//! the stored value
KOKKOS_INLINE_FUNCTION
Real CkEmW(const Real emw, const Real kro) {
  return CkMulRn(emw, kro);
}
//! \fn Real CkE0
//! \brief ck-lin2: the half-layer e0 = 1 - e^{-dtau/mu}, dtau = 0.5 kappa rho dz, as
//! ck_coef stores it in ck_c0 (ck_coef forms it by this function): the linear kernels
//! form it from kappa rho and dz instead of reading ck_c0
KOKKOS_INLINE_FUNCTION
Real CkE0(const Real kro, const Real dz, const Real mu) {
  const Real dtau = 0.5*kro*dz;
  const RtF x = static_cast<RtF>(dtau/mu);
  return static_cast<Real>(-RT_EXPM1(-x));
}
//! \fn void CkLayW
//! \brief ck-lin2: lP slots 2-4 of the layer joining cells i and i+1 -- dslv/dB_i,
//! dsuu/dB_{i+1} and the face interpolation dt_l/dtc -- from the two cells' kappa rho and
//! dz by ck_lin_build's expressions (ck_lin_build forms them by this function): the
//! linear kernels form them instead of reading lP
KOKKOS_INLINE_FUNCTION
void CkLayW(const Real kro, const Real kru, const Real dz, const Real dzu, const bool bf,
            Real &wl, Real &wu, Real &ffj) {
  wl = BFaceW(kru, kro, bf);
  wu = BFaceW(kro, kru, bf);
  const Real dt_l = 0.5*kro*dz;
#if defined(__CUDA_ARCH__)
  // ck-lin2-fma: nvcc contracts the sum as fma(0.5 kru, dzu, dt_l) in ck_lin_build,
  // which stored these weights before; in the kernels that now form them it chose
  // otherwise.  The rounding is spelled out, so every kernel forms the numbers the
  // stored ones were (GPU-bitwise).  Host builds do not contract (a*b + c as written).
  const Real dtc = __fma_rn(0.5*kru, dzu, __dmul_rn(0.5*kro, dz));
#else
  const Real dt_u = 0.5*kru*dzu;
  const Real dtc = dt_l + dt_u;
#endif
  ffj = (dtc > 0.0) ? (dt_l/dtc) : 0.5;
}

//----------------------------------------------------------------------------------------
//! \fn void RTLayer
//! \brief one exact short-characteristic step through a layer (see rt_layer_legacy).
//!
//! The layer has optical thickness dtau along the ray (mu the direction cosine) and a
//! source linear in tau from s_in, where the ray ENTERS, to s_out, where it LEAVES.  The
//! exponential coefficients are the ones the sweeps have always used, with their series
//! forms below x = 1e-3; what changes with the new layers is only which interval they are
//! handed.  absorbed and emitted come back separately because the per-cell source needs
//! them apart, and because their difference is exactly the change in the intensity, which
//! is what keeps Src_g and the divergence of the face fluxes the same number.
KOKKOS_INLINE_FUNCTION
void RTLayer(const Real dtau, const Real mu, const Real s_in, const Real s_out,
             Real &intens, Real &absorbed, Real &emitted) {
  const Real x = dtau/mu;
  const Real e0 = -expm1(-x);
  const Real c_in  = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0 - SQR(x)/3.0);
  const Real c_out = (x > 1.0e-3) ? (1.0 - e0/x)      : (x/2.0 - SQR(x)/6.0);
  absorbed = e0*intens;
  emitted  = c_in*s_in + c_out*s_out;
  intens   = (1.0 - e0)*intens + emitted;
}

KOKKOS_INLINE_FUNCTION
Real LimitRTSource(const Real de, const Real eint, const Real de_max) {
  const Real cap = de_max*eint;
  return (de > cap) ? cap : ((de < -cap) ? -cap : de);
}

KOKKOS_INLINE_FUNCTION
void get_albedo(const Real &Teff, const Real &gg, Real &A) {
  // Parmentier+2015
  Real X = log10(Teff);
  Real g = gg*0.01;
  Real a,b;
  if (Teff < 250.0) {
    a = -0.335*pow(g,0.07);
    b = 0.0;
  } else if (Teff < 750.0) {
    a = -0.335*pow(g,0.07) + 2.149*pow(g,0.135);
    b = -0.896*pow(g,0.135);
  } else if (Teff < 1250.0) {
    a = -0.335*pow(g,0.07) - 0.428*pow(g,0.135);
    b = 0.0;
  } else {
    a = 16.947 - 3.174*pow(g,0.07) - 4.051*pow(g,0.135);
    b = -5.472 + 0.917*pow(g,0.07) + 1.170*pow(g,0.135);
  }
  Real log10A = a + b*X;
  A = pow(10.0,log10A);
  return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2,
                            Real &gamv3, Real &beta, Real &gamir1, Real &gamir2) {
  // Parmentier & Giollot 2014; Parmentier+2015; Roth+2024
  Real X = log10(Teff);
  Real a3, a2, a1, b3, b2, b1, ab, bb;

  Real ap = -2.36;
  Real bp = 13.92;
  Real cp = -19.38;
  if (Teff >= 1400.0 && Teq < 1800.0) {
    ap = -12.45;
    bp = 82.25;
    cp = -134.42;
  }

  if (Teff < 2000.0) {
    ab = 0.84;
    bb = 0.0;
  } else {
    ab = 6.21;
    bb = -1.63;
  }
  if (Teff >= 1400.0 && Teq < 1800.0) {
    ab = 3.0;
    bb = -0.69;
  }

  if (Teff < 200.0) {
    a3 = -3.03;
    b3 = -0.2;
    a2 = -7.37;
    b2 = 2.53;
    a1 = -5.51;
    b1 = 2.48;
  } else if (Teff < 300.0) {
    a3 = -13.87;
    b3 = 4.51;
    a2 = 13.99;
    b2 = -6.75;
    a1 = 1.23;
    b1 = -0.45;
  } else if (Teff < 600.0) {
    a3 = -11.95;
    b3 = 3.74;
    a2 = -15.18;
    b2 = 5.02;
    a1 = 8.65;
    b1 = -3.45;
  } else if (Teff < 1400.0) {
    a3 = -6.97;
    b3 = 1.94;
    a2 = -10.41;
    b2 = 3.31;
    a1 = -12.96;
    b1 = 4.33;
  } else if (Teff < 2000.0) {
    a3 = -3.65;
    b3 = 0.89;
    a2 = -19.95;
    b2 = 6.34;
    a1 = -23.75;
    b1 = 7.76;
  } else {
    a3 = -6.02;
    b3 = 1.61;
    a2 = 13.56;
    b2 = -3.81;
    a1 = 12.65;
    b1 = -3.27;
  }
  if (Teff >= 1400.0 && Teq < 1800.0) {
    if (Teff < 2000.0) {
      a3 = 0.02;
      b3 = -0.28;
      a2 = 6.96;
      b2 = -2.21;
      a1 = -1.68;
      b1 = 0.75;
    } else {
      a3 = -16.54;
      b3 = 4.74;
      a2 = -2.4;
      b2 = 0.62;
      a1 = 10.37;
      b1 = -2.91;
    }
  }
  Real log10gamv1 = a1 + b1*X;
  Real log10gamv2 = a2 + b2*X;
  Real log10gamv3 = a3 + b3*X;
  Real log10gamp = ap*SQR(X) + bp*X + cp;
  beta = ab + bb*X;

  gamv1 = pow(10.0,log10gamv1);
  gamv2 = pow(10.0,log10gamv2);
  gamv3 = pow(10.0,log10gamv3);
  Real gamp = pow(10.0,log10gamp);
  Real dum = (gamp-1.0)/(2.0*beta*(1.0-beta));
  Real R = 1.0 + dum + sqrt(SQR(dum)+dum);
  gamir1 = beta + R - beta*R;
  gamir2 = gamir1/R;

  return;
}

KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr) {
  // Freedman+2014, shared with the conduction module (utils/rosseland.hpp)
  kapr = RosselandFreedman2014(T, p, met);
  return;
}

KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint) {
  // Thorngren+2019 +Erratum
  Real boltz_sigma = 5.6704e-5;
  Real F = 4.0*boltz_sigma*SQR(SQR(Teq));
  Tint = 0.39*Teq*exp(-SQR(log10(F)-9.0-0.14)/1.095);
  Tint = (Tint < 100.0) ? 100.0 : Tint;
  return;
}

// defined in utils/two_stream_rt_pass.cpp (tusplit-1102)
void picket_fence_two_stream_RT_pass(Mesh *pm, Real bdt);

//----------------------------------------------------------------------------------------
//! \fn CkConserveFlux
//! \brief problem/ck_impl_conserve (0 = off, bitwise): make the face fluxes a
//! ck_implicit call REPORTS (rt_Fb: Lir_top, Lrad_bot, the budget ledger) consistent with
//! the energy the call actually APPLIED.  The applied increment of a cell is
//! dep = e - e*; the sweep's fluxes book h S (h = bdt).  They differ by the residual
//! r = dep - h S, which is ~0 in a converged cell but NOT in a cell the bound-constrained
//! Newton leaves on a bound (ck_impl_floorbound: a cell on e_floor that the fluxes would
//! cool further, r > 0; ck_impl_kkt_demax: a cell on e* (1 -+ demax), either sign):
//! such a cell is converged FOR THE CONSTRAINED PROBLEM, i.e. by exclusion from the
//! test, and its fluxes still carry the unconstrained emission.  Mode 1 corrects exactly
//! those KKT cells (the same test as the fused kernel's kktcell, on the final e); mode 2
//! corrects every cell of the column (the Newton tolerance and ck_impl_pred included).
//! The correction keeps the deepest face (icut, Lrad_bot) and every cell below the
//! corrected one, whose balance is converged against the downward half of the corrected
//! cell's emission, untouched, and takes the retained power r V/h out of every face
//! ABOVE the cell: A_f F_f -= sum_{i < f} r_i V_i/h.  That is the unique face-flux set
//! with A F(icut) fixed for which A F_i - A F_i+1 + Qb V = dep V/h in EVERY corrected
//! cell and is unchanged elsewhere, so sum over the column of applied change = booked
//! boundary fluxes.  The correction is split over the band blocks in proportion to each
//! block's emission Em in the cell it comes from (block 0 where Em = 0).  Only the
//! REPORTED face fluxes rt_Fb change (Src_g, Em_g and ck_src keep the sweep's values):
//! no pass, no later call and no state reads rt_Fb before the next sweep rewrites it,
//! so the solution is bitwise the old one (gate: dhj-ck-conserve 09-29).  Counted
//! (problem/ck_impl_verbose, "cons=" / "consde=" on the ck_implicit line): the cells
//! corrected in the last call and their sum r V [erg], this rank.

void CkConserveFlux(Mesh *pm, const Real bdt);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

//----------------------------------------------------------------------------------------
//! \fn CkNonconvLoc
//! \brief problem/ck_impl_ncloc: after a NOT-CONVERGED ck_implicit call, print where the
//! columns still active (ck_done = 0) are: gid, m, k, j, the cell i of the column's
//! largest relative residual |e - e* - dt S|/(e + eps e_max) (the fused test's norm;
//! KKT cells skipped, as the test skips them),
//! its pressure [bar] and T from the last sweep, and the column's top pressure.  Print
//! only; nothing is written.

void CkNonconvLoc(Mesh *pm, const Real bdt);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

//----------------------------------------------------------------------------------------
//! \fn CkNcRefSnap / CkNcRefCmp
//! \brief problem/ck_impl_ncref (see its note): latch the running increment of the call
//! after pass N, and at the end of the call compare it with the final one.  Print only.

void CkNcRefSnap(Mesh *pm);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

void CkNcRefCmp(Mesh *pm, const int nsnap, const int npass, const bool conv);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

//----------------------------------------------------------------------------------------
//! \fn picket_fence_two_stream_RT
//! \brief the two-stream source step.

void picket_fence_two_stream_RT_body(Mesh *pm, Real bdt);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

//----------------------------------------------------------------------------------------
//! \fn void picket_fence_two_stream_RT
//! \brief the outer radiative call: the body above, then the rank-uniform once-per-call
//! point where the source-limiter clip count is reduced (see RTSourceLimiterReduce).
void picket_fence_two_stream_RT(Mesh *pm, Real bdt);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)


//----------------------------------------------------------------------------------------
//! \struct TsrtX1F, TsrtEiN, TsrtRhoN, TsrtAFC, TsrtACC, TsrtBTF, TsrtMIXF
//! \brief tusplit-1102: the seven helpers of picket_fence_two_stream_RT_pass that ride
//! in the ctx tuples of TsrtCkChain / TsrtGreyRt (X1F, eiN, rhoN, AFC, ACC, BTF, MIXF).
//! They were KOKKOS_LAMBDAs local to the pass; as named functors (the same bodies, the
//! same by-value captures as members, auto return type as the lambda's) the tuple types
//! can be spelled outside the pass, so the chain and grey templates are explicitly
//! instantiated in their own translation units (two_stream_rt_chain*.cpp).  See the pass
//! for what each one is.
struct TsrtX1F {
  bool pp_;
  int is_pp, nx1_pp;
  DualArray1D<RegionSize> size;
  DvceArray2D<Real> x1f_;
  KOKKOS_INLINE_FUNCTION auto operator()(const int m, const int i) const {
      return pp_ ? LeftEdgeX(i-is_pp, nx1_pp, size.d_view(m).x1min,
                             size.d_view(m).x1max)
                 : x1f_(m,i);
  }
};
struct TsrtEiN {
  bool usecons_;
  DvceArray5D<Real> w0_uc_, u0_uc_;
  bool cs_uc_;
  DvceArray3D<Real> cosc_uc_;
  bool etg_uc_;
  DvceArray4D<Real> phicc_uc_;
  bool mhd_uc_;
  DvceArray5D<Real> bcc_uc_;
  DvceArray1D<int> eicl_g;
  EOS_Data eos_uc_;
  KOKKOS_INLINE_FUNCTION auto operator()(const int m, const int k, const int j,
                                         const int i) const {
      if (!usecons_) return w0_uc_(m,IEN,k,j,i);
      const Real ei_uc = EintFromCons(u0_uc_, m, k, j, i,
                                      cs_uc_ ? cosc_uc_(m,k,j) : 0.0, cs_uc_,
                                      etg_uc_, etg_uc_ ? phicc_uc_(m,k,j,i) : 0.0,
                                      mhd_uc_ ? MagEnergyCC(bcc_uc_,m,k,j,i) : 0.0);
      if (ei_uc > 0.0) return ei_uc;           // false for NaN too, which is the point
      Kokkos::atomic_fetch_add(&eicl_g(0), 1);
      Real ei_fl = eos_uc_.EnergyFromTemperature(u0_uc_(m,IDN,k,j,i), eos_uc_.tfloor);
      if (!(ei_fl > 0.0)) ei_fl = 1.0e-300;
      return ei_fl;
  }
};
struct TsrtRhoN {
  bool usecons_;
  DvceArray5D<Real> u0_uc_, w0_uc_;
  DvceArray1D<int> stcl_g;
  Real dfl_uc_;
  KOKKOS_INLINE_FUNCTION auto operator()(const int m, const int k, const int j,
                                         const int i) const {
      const Real d = usecons_ ? u0_uc_(m,IDN,k,j,i) : w0_uc_(m,IDN,k,j,i);
      if (d > 0.0) return d;                 // false for NaN too, as in eiN
      Kokkos::atomic_fetch_add(&stcl_g(0), 1);
      return dfl_uc_;
  }
};
struct TsrtAFC {
  bool pp_;
  DvceArray4D<Real> area1;
  KOKKOS_INLINE_FUNCTION auto operator()(const int m, const int k, const int j,
                                         const int i) const {
      return pp_ ? 1.0 : area1(m,k,j,i);
  }
};
struct TsrtACC {
  bool pp_;
  DvceArray4D<Real> volume, dx1;
  KOKKOS_INLINE_FUNCTION auto operator()(const int m, const int k, const int j,
                                         const int i) const {
      return pp_ ? 1.0 : volume(m,k,j,i)/dx1(m,k,j,i);
  }
};
struct TsrtBTF {
  bool pp_;
  DvceArray4D<Real> area1, volume, dx1;
  int ie;
  KOKKOS_INLINE_FUNCTION auto operator()(const int m, const int k, const int j,
                                         const int f, const int ic) const {
      if (pp_) return 0.0;
      const Real ab = (f <= ic) ? area1(m,k,j,f) : volume(m,k,j,f-1)/dx1(m,k,j,f-1);
      const Real aa = (f > ie) ? area1(m,k,j,f) : volume(m,k,j,f)/dx1(m,k,j,f);
      return (aa - ab)/(aa + ab);
  }
};
struct TsrtMIXF {
  KOKKOS_INLINE_FUNCTION auto operator()(const Real bt, const Real ub,
                                         const Real da) const {
      return bt*(da - ub);
  }
};

//----------------------------------------------------------------------------------------
//! \fn TsrtCkChain, TsrtGreyRt
//! \brief the correlated-k chain and grey column launches of
//! picket_fence_two_stream_RT_pass (DECLARATIONS).  tusplit-1102: the bodies live in
//! utils/two_stream_rt_chain.hpp and are explicitly instantiated, for the tag sets the
//! pass dispatches, in utils/two_stream_rt_chain_*.cpp; the pass itself is
//! utils/two_stream_rt_pass.cpp.  TsrtCkChainCtx / TsrtGreyRtCtx are the types of the
//! pass's ck_chain_ctx / grey_rt_ctx (std::forward_as_tuple of its locals, in that
//! order; the pass static_asserts the match).
using TsrtCkChainCtx = std::tuple<
    TsrtACC &, TsrtAFC &, Real &, Real &, DvceArray5D<Real> &, const bool &, Real &,
    TsrtBTF &, const bool &, DvceArray4D<Real> &, DvceArray5D<Real> &,
    DvceArray5D<Real> &, DvceArray5D<Real> &, DvceArray5D<Real> &, const bool &,
    DvceArray3D<Real> &, const bool &, bool &, bool &, DvceArray1D<Real> &,
    DvceArray5D<Real> &, bool &, DvceArray5D<Real> &, DvceArray4D<Real> &, const int &,
    DvceArray2D<Real> &, const bool &, DvceArray1D<Real> &, DvceArray4D<Real> &,
    DvceArray4D<Real> &, DvceArray4D<Real> &, DvceArray5D<Real> &, DvceArray5D<Real> &,
    Real &, Real &, const bool &, DvceArray3D<int> &, DvceArray3D<int> &, int &,
    bool &, int &, int &, int &, DvceArray5D<Real> &, int &, int &, const bool &,
    TsrtMIXF &, Real (&)[2], int &, const int &, int &, Real &, DvceArray4D<Real> &,
    const Real &, const Real &, DvceArray5D<Real> &, TsrtRhoN &, DvceArray5D<Real> &,
    DvceArray4D<Real> &, const bool &, Real &, Real &, Real (&)[2], TsrtX1F &,
    DvceArray2D<Real> &, DvceArray4D<Real> &, DvceArray4D<Real> &, const int &,
    const int &, DvceArray5D<Real> &, DvceArray5D<Real> &, const bool &>;
using TsrtGreyRtCtx = std::tuple<
    Real &, Real &, DvceArray4D<Real> &, Real &, Real &, const Real &,
    DvceArray4D<Real> &, TsrtEiN &, EOS_Data &, Real &, Real &, Real &, const bool &,
    Real &, int &, Real &, int &, int &, int &, int &, int &, DualArray1D<int> &,
    Real &, Real (&)[2], const bool &, const int &, int &, int &, DvceArray1D<int> &,
    DvceArray1D<Real> &, Real &, Real &, TsrtRhoN &, DvceArray1D<int> &, Real &,
    const bool &, Real &, Real &, DvceArray5D<Real> &, const bool &, const bool &,
    DvceArray4D<Real> &, Real (&)[2], DvceArray4D<Real> &, DvceArray2D<Real> &,
    DvceArray2D<Real> &, DvceArray2D<Real> &, DvceArray2D<Real> &>;
template <typename Ctx, typename NnTag, typename SphTag, typename BspTag, typename CchTag,
          typename FrmTag, typename FopTag, typename JacTag>
void TsrtCkChain(Ctx &ctx_, NnTag nn_tag, SphTag sph_tag, BspTag bsp_tag,
                 CchTag cch_tag, FrmTag frm_tag, FopTag fop_tag, JacTag jac_tag);
template <typename Ctx, typename NnTag>
void TsrtGreyRt(Ctx &ctx_, NnTag nn_tag);

//----------------------------------------------------------------------------------------
//! \fn void CkCadState
//! \brief problem/ck_impl_every: the pieces every cadence kernel needs to read the
//! thermodynamic state of a cell out of u0 exactly as the RT's eiN/rhoN do.
struct CkCadState {
  DvceArray5D<Real> u0;
  DvceArray4D<Real> wtemp;
  DvceArray4D<Real> phi;
  DvceArray3D<Real> cosc;
  DvceArray5D<Real> bcc;
  EOS_Data eos;
  Real gm1, Rgas;
  bool etg, cs, mhd;
  KOKKOS_INLINE_FUNCTION
  Real Eint(const int m, const int k, const int j, const int i) const {
    return EintFromCons(u0, m, k, j, i, cs ? cosc(m,k,j) : 0.0, cs, etg,
                        etg ? phi(m,k,j,i) : 0.0, mhd ? MagEnergyCC(bcc,m,k,j,i) : 0.0);
  }
  // T [K] and de/dT [energy per volume per K] at (d, e); false for an unusable state
  KOKKOS_INLINE_FUNCTION
  bool TCv(const int m, const int k, const int j, const int i, const Real d,
           const Real e, Real &tk, Real &cv) const {
    if (!(d > 0.0) || !(e > 0.0)) return false;
    Real pp;
    PresTempFromEint(eos, gm1, Rgas, d, e, TGuess(wtemp, m, k, j, i), pp, tk);
    if (!(tk > 0.0)) return false;
    if (eos.IsGeneral()) {
      cv = d*eos.SpecificHeatCv(d, e, tk/eos.temp_cgs)/eos.temp_cgs;
    } else {
      cv = e/tk;
    }
    return (cv > 0.0);
  }
};

inline CkCadState CkCadMakeState(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  CkCadState s;
  const bool hyd = (pmbp->phydro != nullptr);
  s.u0 = hyd ? pmbp->phydro->u0 : pmbp->pmhd->u0;
  s.wtemp = hyd ? pmbp->phydro->wtemp : pmbp->pmhd->wtemp;
  s.phi = hyd ? pmbp->phydro->phicc0 : pmbp->pmhd->phicc0;
  s.etg = hyd ? pmbp->phydro->use_etotgrav : pmbp->pmhd->use_etotgrav;
  s.eos = hyd ? pmbp->phydro->peos->eos_data : pmbp->pmhd->peos->eos_data;
  s.mhd = !hyd;
  s.bcc = hyd ? CkDum<DvceArray5D<Real>>("ck_cad_bcc_d") : pmbp->pmhd->bcc0;
  s.cs = pm->use_cubed_sphere;
  s.cosc = pmbp->pcoord->cos_cell;
  s.gm1 = s.eos.gamma - 1.0;
  s.Rgas = pm->pgen->hot_jupiter_param.Rgas;
  return s;
}

//----------------------------------------------------------------------------------------
//! \fn void CkCadStore
//! \brief after a full implicit call: store Q0, D, T0, rho0 for every cell of every
//! column the call solved (all of them, or the unmasked ones of a guard call).
//! Q0 = (e - e*)/dt is the rate the gas actually received (ck_dep is the call's total
//! applied increment), so the stored rate is exact whatever ended the Newton (tol, the
//! prediction of ck_impl_pred).  D is the Newton diagonal jac1 = dS_i/dT_i (thick
//! cell) or -4 E/T (a thin cell, whose own implicit step CkThinSolve linearises exactly
//! that), clamped to <= 0 so that c_v - D dt >= c_v.  Cells below the cut get Q0 = D = 0.
//! With diag, slots 6/7 of ck_cad_stat get sum dx |Q_new - (Q0 + D (T - T0))_old| and
//! sum dx |Q_new| over the refreshed cells: the linear model's error at the refresh.

void CkCadStore(Mesh *pm, const Real dt, const bool partial, const bool diag);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

//----------------------------------------------------------------------------------------
//! \fn int CkCadLin
//! \brief the linearised step of problem/ck_impl_every.  One team per column.  With
//! ck_impl_every_thr > 0 the column's largest relative change of T or rho since its last
//! full solve is measured first; above thr the column is left untouched and marked for a
//! full solve this step (mask 0), otherwise (mask 2) every cell takes
//!     dT = (Q0 + D (T - T0)) dt / (c_v - D dt),  de = c_v dT,
//! implicit in T, unconditionally stable for D <= 0 (|dT| <= |Q0 dt/c_v| + |T - T0|, and
//! the step moves T toward T0 - Q0/D), with |de| capped at min(dtmax, demax) e.
//! Stat slots: 0 columns flagged, 1 sum dx de, 2 sum dx Q0 dt, 3 capped cells,
//! 4 max deviation, 5 sum dx |Q0| dt.  Returns the columns flagged on this rank.

int CkCadLin(Mesh *pm, const Real dt, const bool diag);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

//----------------------------------------------------------------------------------------
//! \fn void CkCadStep
//! \brief problem/ck_impl_every: the whole correlated-k RT of one hydro step (called
//! once per cycle by the pgen's split hook in place of picket_fence_two_stream_RT).
//! Full call on ncycle % N == 0 (and on the first call of a run or a restart), the
//! linearised step otherwise, with a full call for the columns the guard flags.

void CkCadStep(Mesh *pm, const Real dt);
// (defined in utils/two_stream_rt_driver.cpp, tusplit-1102)

}  // namespace two_stream_rt

#endif  // UTILS_TWO_STREAM_RT_HPP_
