#ifndef RAD_VET_RAD_VET_ARGS_HPP_
#define RAD_VET_RAD_VET_ARGS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_args.hpp
//! \brief array slots, kernel launchers and CleanArgs, what the device functions of the
//! solve read.

#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include "athena.hpp"
#include "globals.hpp"
#include "bvals/bvals.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "reconstruct/plm.hpp"

// CUDA: the face flux is not inlined (inlined, nvcc takes over an hour)
#if defined(__CUDACC__)
#define CLEAN_DEVICE_NOINLINE __noinline__
#else
#define CLEAN_DEVICE_NOINLINE
#endif
// HIP: face flux, stress and source forced inline (as calls they need a large scratch
// frame); on paths a pass does not take they stay calls (CleanFaceFluxCall, ...) to
// keep the register use of the taken path low
#if defined(__HIPCC__)
#define CLEAN_HIP_INLINE __attribute__((always_inline))
#define CLEAN_HIP_NOINLINE __attribute__((noinline))
#else
#define CLEAN_HIP_INLINE
#define CLEAN_HIP_NOINLINE
#endif
// HIP: kernels read CleanArgs from a device copy (CleanDevArgs), since a capture indexed
// at run time is copied to each lane's stack; CUDA keeps the capture
#if defined(KOKKOS_ENABLE_HIP)
#define CLEAN_HIPX(...) __VA_ARGS__
#else
#define CLEAN_HIPX(...)
#endif

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//! \fn void CleanMaxRanks
//! \brief max of v[0..n-1] over all ranks; non-finite -> inf so all ranks agree
inline void CleanMaxRanks(Real *v, const int n) {
  for (int q = 0; q < n; ++q) {
    if (!std::isfinite(v[q])) v[q] = std::numeric_limits<Real>::infinity();
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, v, n, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
}

// CleanFor: a 4-D clean pass.  CUDA: launch bounds of kCudaMb blocks of 256 per SM (a
// register cap), par_for's index arithmetic (bitwise).  HIP: kHipMb blocks if > 0
// (4: 128 VGPRs), else par_for.  Host: par_for
#if defined(ATHENA_CLEAN_FOR_MB)   // build-time override (register cap only, bitwise)
constexpr int kCleanForMb = ATHENA_CLEAN_FOR_MB;
#elif defined(KOKKOS_ARCH_HOPPER90)
constexpr int kCleanForMb = 3;
#else
constexpr int kCleanForMb = 2;
#endif
template <int kHipMb = 0, int kCudaMb = kCleanForMb, typename F>
inline void CleanFor(const char *name, const int nmb1, const int ks, const int ke,
                     const int js, const int je, const int is, const int ie,
                     const F &f) {
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int mb = kHipMb;
#elif defined(KOKKOS_ENABLE_CUDA)
  constexpr int mb = kCudaMb;
#else
  constexpr int mb = 0;
#endif
  if constexpr (mb > 0) {
    const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
    AthFor(name, Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, mb>>(
           DevExeSpace(), 0, (nmb1 + 1)*nkji),
    KOKKOS_LAMBDA(const int idx) {
      const int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/ni;
      f(m, ks + k, js + j, is + idx - m*nkji - k*nji - j*ni);
    });
  } else {
    par_for(name, DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie, f);
  }
}

// CleanFor with k fastest on HIP (the row order of a VgdRag slab), so adjacent lanes
// read adjacent intensity rows; only for passes with independent cells (bitwise)
template <typename F>
inline void CleanForK(const char *name, const int nmb1, const int ks, const int ke,
                      const int js, const int je, const int is, const int ie,
                      const F &f) {
#if defined(KOKKOS_ENABLE_HIP)
  par_for(name, DevExeSpace(), 0, nmb1, is, ie, js, je, ks, ke,
  KOKKOS_LAMBDA(const int m, const int i, const int j, const int k) {f(m, k, j, i);});
#else
  CleanFor(name, nmb1, ks, ke, js, je, is, ie, f);
#endif
}

// a pass of columns (one thread per (m, k, j), marching in i); HIP: blocks of 64 with
// k fastest to spread few columns over more CUs
template <typename F>
inline void CleanForCol(const char *name, const int nmb1, const int ks, const int ke,
                        const int js, const int je, const F &f) {
#if defined(KOKKOS_ENABLE_HIP)
  const int nk = ke - ks + 1, njk = (je - js + 1)*nk;
  Kokkos::parallel_for(name, Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<64,
                       1>>(DevExeSpace(), 0, (nmb1 + 1)*njk),
  KOKKOS_LAMBDA(const int idx) {
    const int m = idx/njk, j = (idx - m*njk)/nk;
    f(m, ks + idx - m*njk - j*nk, js + j);
  });
#else
  par_for(name, DevExeSpace(), 0, nmb1, ks, ke, js, je, f);
#endif
}

// closure variables (arrays cln_c and cln_s); cln_s holds the first CLN_FH (cubed
// sphere: CLN_F1)
constexpr int CLN_D11 = 0;   // Eddington tensor D11 D22 D33 D12 D13 D23
constexpr int CLN_HP = 6;    // h+_a = H+_a/J, a = 1,2,3
constexpr int CLN_HM = 9;    // h-_a = H-_a/J (<= 0)
constexpr int CLN_C1 = 12;   // non-locality C1 = <|I - S|>/J
constexpr int CLN_J = 13;    // J of the rays (energy-density units)
constexpr int CLN_F2 = 14;   // rays' flux through the low x2 face (not normalised)
constexpr int CLN_FX = 18;   // rays' flux through the low x1 face
constexpr int CLN_F3 = 22;   // rays' flux through the low x3 face
constexpr int CLN_G = 26;    // curvature of the rays' stress (y_r, y_theta)
constexpr int CLN_JR = 28;   // J also below clean_ray_min E (CLN_J is 0 there)
constexpr int CLN_FH = 29;   // cubed sphere: flux through the high x2, x3 face of a
                             // MeshBlock's last cell (+4; differs from a seam ghost's)
constexpr int CLN_F1 = 37;   // rays' flux through an x1 boundary face
constexpr int CLN_KP = 41;   // k+^d_b = K+^d_bd/J (stress of mu_d > 0), slot 3 d + b
constexpr int CLN_N = 50;
// vet_rotate_avg_ema_wb on the cubed sphere (sc rays): in place of k+ (not held there)
// each cell's half ranges about its faces' normals, h+ out of its high faces (3), h- of
// its low ones (3) (CleanClosureFromGd), so that the ring carries them as cell moments
constexpr int CLN_HO = CLN_KP;
constexpr int CLN_NHO = CLN_HO + 6;

// the slots of cln_s
inline int CleanNs(const bool cub, const bool ho = false) {
  return cub ? (ho ? CLN_NHO : CLN_F1) : CLN_FH;
}

// cell moments (D, h+-, C1, J, k+); the others are face sums
KOKKOS_INLINE_FUNCTION
bool CleanCellMoment(const int n) {return n <= CLN_J || n >= CLN_KP;}

// slot of restart channel n; restart files order them D h+ h- k+ C1 J F2 F1 FX F3 FH G JR
inline int CleanRstSlot(const int n) {
  constexpr int o[9] = {0, 12, 21, 23, 27, 31, 39, 47, CLN_N};
  constexpr int s[8] = {CLN_D11, CLN_KP, CLN_C1, CLN_F2, CLN_F1, CLN_FX, CLN_FH, CLN_G};
  int g = 0;
  while (n >= o[g+1]) ++g;
  return s[g] + n - o[g];
}

// face variables (cln_fc1, cln_fc2, cln_fc3: CleanF).  Slots below CLF_NS live across a
// closure build; the others are rebuilt after one before they are read
constexpr int CLF_N = 0;     // clean_stag_fsolve: CLF_N .. CLF_N + 3 hold F, the
                             // dissipation term, F^H and (clean_stag_fkin) the energy
                             // flux's memory d
constexpr int CLF_DS = 4;    // clean_stag_fkin: d at the stage-loop start
constexpr int CLF_DB = 5;    // clean_stag_vimp: drag share b = beta/(1 + beta),
constexpr int CLF_DQ = 6;    //   face momentum Q = G - c^2 dm_f of the iterate
constexpr int CLF_KD = 7;    // clean_stag_fkin: slopes K1, K2 of d (2)
constexpr int CLF_XC = 9;    // clean_stag_xdonor: C(E_k) and its last change (2)
// (CLF_NS = 11, rad_vet.hpp)
constexpr int CLF_W = CLF_NS;         // blend weight w
constexpr int CLF_TH = CLF_NS + 1;    // theta = 1/(1 + h c chi_f)
constexpr int CLF_G = CLF_NS + 2;     // flux of the old vector G_f
constexpr int CLF_HP = CLF_NS + 3;    // h+_f
constexpr int CLF_HM = CLF_NS + 4;    // h-_f
constexpr int CLF_QL = CLF_NS + 5;    // al, ar (CleanFaceWeights), cross factor
                                      // -(1-w) theta h c^2 (3)
constexpr int CLF_PT = CLF_QL;        // outside the solve: face momentum terms
                                      // (CleanStagDrag) and the RHS's a E (cln_rhsaf)
constexpr int CLF_SR = CLF_QL + 1;    //   and the RHS's face sources (CleanFaceSource)
constexpr int CLF_NF = CLF_QL + 3;    // (clean_stag_xdonor, the bracket: more slots)
constexpr Real CLN_MEMUP = 1.0/6.0;  // clean_stag_gflux: memory outflow above this share
                                     // of the donor's E_old goes upwind

// work vectors (array cln_w): CLW_NS shared, then the staggered layout's own
constexpr int CLW_EO = 0;    // E of the old vector
constexpr int CLW_EK = 1;    // Newton iterate of E
constexpr int CLW_T = 2;     // gas temperature of the iterate
constexpr int CLW_EI = 3;    // gas internal energy of the old vector
// exchange linearisation r = s/(d cv + s), s = 4 h c kP a T^3
constexpr int CLW_RR = 4;
constexpr int CLW_A1 = 5;    // advective coefficients a_1, a_2, a_3
constexpr int CLW_DM = 8;    // momentum change of the iterate (3)
constexpr int CLW_B = 11;    // right-hand side
constexpr int CLW_X = 12;    // solution
constexpr int CLW_BT = 13;   // beta = v/c of the iterate (3)
constexpr int CLW_RE = 16;   // h S_E = hc kP G a T^4 - hc kE E + RE E + RY.y - W (Ma et
constexpr int CLW_GM = 17;   //   al. 2025, eq. 56, G = CLW_GM): beta terms (RY: 3);
                             //   staggered: RE = hc kP a T^4, RY 1 - r
constexpr int CLW_EP = 18;   // -(RE E + RY.y) of the solution (gas eint)
constexpr int CLW_JV = 19;   // (1 - r) dA/dbeta/(c^2 rho), dS_p/dbeta/(c^2 rho) (3 + 1):
constexpr int CLW_GX = 22;   //   implicit gas velocity; gflux: JV + 1 = hc kE dop
constexpr int CLW_RY = 23;
constexpr int CLW_NS = 24;
constexpr int CLW_AM = 24;   // staggered layout: x1-line coefficient of E_{i-1}
constexpr int CLW_A0 = 25;   // diagonal
constexpr int CLW_AP = 26;   // x1-line coefficient of E_{i+1}
constexpr int CLW_CP = 27;   // tridiagonal solve scratch
constexpr int CLW_R = 28;    // BiCGStab vectors
constexpr int CLW_R0 = 29;
constexpr int CLW_P = 30;
constexpr int CLW_V = 31;
constexpr int CLW_S = 32;
constexpr int CLW_Q = 33;
constexpr int CLW_Y = 34;
constexpr int CLW_Z = 35;
constexpr int CLW_TM = 36;   // mg: coefficients of E_{j-1}, E_{j+1}, E_{k-1}, E_{k+1} (4)
constexpr int CLW_L = 40;    // cubed sphere: lateral stress sums Y of face family q at
                             // CLW_L + 3 q (CleanStagCsLateral)
constexpr int CLW_SE = 49;   // clean_conv = budget: gross terms S of the E and gas rows
constexpr int CLW_SG = 50;   //   (their scales), the weighted B and a preconditioner
constexpr int CLW_BW = 51;   //   input (CleanBiCGStab)
constexpr int CLW_PW = 52;
constexpr int CLW_N = 53;
constexpr int CLW_PD = 55;   // the rows' diagonal of the work (53, 54 unused)
constexpr int CLW_XA = 56;   // B less the linearised exchange,
constexpr int CLW_XW = 57;   //   the row's work W
constexpr int CLW_NX = 58;

//----------------------------------------------------------------------------------------
//! \struct CleanArgs
//! \brief everything the device functions of the clean solve read

struct CleanArgs {
  DvceArray5D<Real> cs;        // closure of the stage
  CleanW w;                    // work vectors
  CleanF fc[3];                // face variables in x1, x2, x3
  DvceArray3D<Real> beam;      // boundary beam flux at x1min (m, k-ks, j-js)
  DvceArray4D<Real> bstress;   // beam stress column c K_b1 at x1min (m, b, k-ks, j-js)
  Real h, c, idx[3];           // stage step, light speed, inverse cell widths
  int bc[6];                   // boundary types x1min, x1max, ..., x3max
  DvceArray2D<BoundaryFlag> mbc;  // boundary flags of the six faces of each MeshBlock
  Real ebath[6];               // boundary bath energy densities
  Real iflux[2];               // imposed fluxes of x1 flux faces, x1min/max
  DvceArray3D<Real> ifx2[2];   // imposed fluxes of x2 flux faces (m, k, i), x2min/max
  bool use_beam;
  bool feedback;               // gas coupling (advective flux, exchange linearisation)
  bool vzero;                  // gas momentum frozen (cln_vzero, steady start)
  bool fd;                     // clean_stag_fkin: energy flux G + F^H - F
  bool gf;                     // clean_stag_gflux: low-order flux factors
  bool fre;                    // clean_stag_freal: realisable comoving sb
  int fdg;                     // bracket: face slot of its lagged flux DG (-1: none)
  bool fs;                     // clean_stag_fsolve: face F of the F solve (CLF_N)
  bool hi;                     // clean_stag_split: face momentum state F^H (w = 0)
  int xq;                      // clean_stag_xdonor: face slots qL, qR of the lagged
                               // donor part, explicit part, C(E_k) (-1: central)
  bool mcl;                    // face momentum term capped at h chi_f E (realisable)
  DvceArray1D<Real> mclc;      // (mcl: count of face terms the cap moved)
  bool fa;                     // store the faces' momentum terms (CLF_PT, CleanStagDrag)
  bool drag;                   // clean_stag_vimp: implicit gas velocity
  bool upw;                    // clean_stag_upwind: donor-cell advective flux
  bool fwc;                    // al, ar and the cross factor from CLF_QL
  bool sar;                    // streaming share's net flux in area form
  DvceArray5D<Real> op;        // opacities (VET_OP_*)
  int nf;                      // number of flux components (dimensions)
  int is, ie, js, je, ks, ke;
  bool multi_d, three_d;
  bool sph;                    // spherical polar (r, theta, phi), orthonormal basis
  DvceArray4D<Real> vol;       // cell volumes
  DvceArray4D<Real> area[3];   // face areas normal to x1, x2, x3
  DvceArray4D<Real> dxf[3];    // distances between the cell centres across those faces
  DvceArray2D<Real> x1v, x1f;  // cell-centre (volume centroid) and face radii
  DvceArray2D<Real> x2v, x2f;  // cell-centre and face polar angles
  DvceArray3D<Real> tg;         // tan, tan, sin, sin of x2v, x2f, x2v, x2f by (m, q, j)
  Real cseam;                  // cot of a theta-periodic seam face: its two copies' mean
  // cubed sphere: y in the cell's frame (r, theta, phi); g the geometry (CLG_*)
  bool cub;
  bool ys = false;    // clean_stag_lat_sums: lateral sums Y off the cubed sphere
  bool lwb = false;   // clean_cs_lat_wb (CleanStagCsLat)
  DvceArray5D<Real> yl;   //   (cln_yl: CLW_L + 3 q on the cubed sphere)
  DvceArray4D<Real> ae0;   // cubed sphere: pcoord's radial face area (CleanDivW)
  CleanGeo g;
  // the gas energy also holds d phi (<hydro>/etotgrav) and |B|^2/2 (<mhd>)
  bool etg, mhd;
  DvceArray4D<Real> phi, emag;
};

// staggered operator row: coefficients of E at the cell and its 18 face/edge neighbours.
// Slot of (dk, dj, di): 0 the cell; 1 + 2 a + (d > 0) faces along axis a (i, j, k);
// 7 + 4 p + 2 (d1 > 0) + (d2 > 0) edges in plane p of axes (j, i), (k, i), (k, j)
constexpr int CLN_ST = 19;
// cubed sphere: coefficients of a cell's lateral sums Y_q,b on E of its neighbour
// across face s of family q (slot 3 (2 q + s) + b) and on its own E (18 + 3 q + b)
constexpr int CLN_CSY = 27;

KOKKOS_INLINE_FUNCTION
constexpr int CleanStSlot(const int dk, const int dj, const int di) {
  return ((dk != 0) + (dj != 0) + (di != 0) == 0) ? 0 :
         ((dk != 0) + (dj != 0) + (di != 0) == 1) ?
         1 + 2*((dk != 0) ? 2 : (dj != 0)) + (dk + dj + di > 0) :
         7 + 4*((dk == 0) ? 0 : ((dj == 0) ? 1 : 2)) + 2*(((dk == 0) ? dj : dk) > 0) +
         (((di == 0) ? dj : di) > 0);
}

KOKKOS_INLINE_FUNCTION
void CleanStOff(const int s, int &dk, int &dj, int &di) {   // (the offset of slot s)
  dk = dj = di = 0;
  if (s > 0 && s < 7) {
    const int d = ((s - 1) % 2) ? 1 : -1;
    ((s < 3) ? di : ((s < 5) ? dj : dk)) = d;
  } else if (s >= 7) {
    const int p = (s - 7)/4, d1 = (((s - 7)/2) % 2) ? 1 : -1, d2 = ((s - 7) % 2) ? 1 : -1;
    ((p == 0) ? dj : dk) = d1;
    ((p == 2) ? dj : di) = d2;
  }
}

//----------------------------------------------------------------------------------------
//! \fn Real CleanGasEint
//! \brief gas internal energy of uh at (k,j,i): total less ekin, d phi (etotgrav) and
//! the magnetic energy, as RadiationVet::Opacity forms it

KOKKOS_INLINE_FUNCTION
Real CleanGasEint(const CleanArgs &x, const DvceArray5D<Real> &uh, const Real ekin,
                  const int m, const int k, const int j, const int i) {
  Real eint = uh(m,IEN,k,j,i) - ekin;
  if (x.etg) eint -= uh(m,IDN,k,j,i)*x.phi(m,k,j,i);
  if (x.mhd) eint -= x.emag(m,k,j,i);
  return eint;
}

//----------------------------------------------------------------------------------------
//! \fn int CleanBc
//! \brief boundary type of side s of MeshBlock m: the mesh's on a physical boundary,
//! VET_IBC_PERIODIC on a face shared with a MeshBlock or panel

KOKKOS_INLINE_FUNCTION
int CleanBc(const CleanArgs &x, const int m, const int s) {
  const BoundaryFlag f = x.mbc(m,s);
  return (f == BoundaryFlag::block || f == BoundaryFlag::periodic ||
          f == BoundaryFlag::panel) ? VET_IBC_PERIODIC : x.bc[s];
}

//----------------------------------------------------------------------------------------
//! \fn bool CleanOpen
//! \brief true if side s of MeshBlock m lets the rays' flux through (not wall or flux)

KOKKOS_INLINE_FUNCTION
bool CleanOpen(const CleanArgs &x, const int m, const int s) {
  const int bc = CleanBc(x, m, s);
  return bc != VET_IBC_REFLECT && bc != VET_IBC_FLUX;
}

//----------------------------------------------------------------------------------------
//! \struct CleanShell
//! \brief cells of an n3 x n2 x n1 array outside the inner box [k0,k1]x[j0,j1]x[i0,i1]
//! (inset by g), numbered plane by plane

struct CleanShell {
  int n1, n2, n3, i0, i1, j0, j1, k0, k1;
  int nio, ncell;
  CleanShell(const int n3_, const int n2_, const int n1_, const int g, const int is,
             const int ie, const int js, const int je, const int ks, const int ke,
             const bool md, const bool td) :
      n1(n1_), n2(n2_), n3(n3_), i0(is + g), i1(ie - g), j0(md ? js + g : 0),
      j1(md ? je - g : n2_ - 1), k0(td ? ks + g : 0), k1(td ? ke - g : n3_ - 1) {
    nio = n1 - (i1 - i0 + 1);
    ncell = n3*n2*n1 - (k1 - k0 + 1)*(j1 - j0 + 1)*(i1 - i0 + 1);
  }
  bool Empty() const {return i1 < i0 || j1 < j0 || k1 < k0;}
};

//----------------------------------------------------------------------------------------
//! \fn void CleanShellCell
//! \brief cell (k,j,i) number q of the shell s

KOKKOS_INLINE_FUNCTION
void CleanShellCell(const CleanShell &s, int q, int &k, int &j, int &i) {
  const int np = s.n2*s.n1, nj = s.j1 - s.j0 + 1;
  const int pl = np - nj*(s.i1 - s.i0 + 1), nb = s.k0*np, nm = (s.k1 - s.k0 + 1)*pl;
  if (q < nb || q >= nb + nm) {
    if (q >= nb) q -= nb + nm - (s.k1 + 1)*np;
    k = q/np;
    j = (q - k*np)/s.n1;
    i = q - k*np - j*s.n1;
    return;
  }
  q -= nb;
  k = s.k0 + q/pl;
  int r = q%pl;
  const int rb = s.j0*s.n1, rm = nj*s.nio;
  if (r < rb || r >= rb + rm) {
    if (r >= rb) r -= rb + rm - (s.j1 + 1)*s.n1;
    j = r/s.n1;
    i = r%s.n1;
    return;
  }
  r -= rb;
  j = s.j0 + r/s.nio;
  const int t = r%s.nio;
  i = (t < s.i0) ? t : s.i1 + 1 + (t - s.i0);
}

//----------------------------------------------------------------------------------------
//! \fn void CleanSlotAdd, CleanSlotSum
//! \brief counters of a kernel over many faces or cells: one atomic per thread in one of
//! CLN_NSL slots (a single counter's atomics serialise all threads on a GPU), added to
//! dst(n) after the kernel

constexpr int CLN_NSL = 4096;

KOKKOS_INLINE_FUNCTION
void CleanSlotAdd(const DvceArray1D<Real> &v, const int n, const int m, const int k,
                  const int j, const int i, const Real a) {
  Kokkos::atomic_add(&v(n*CLN_NSL + ((i + 61*j + 1021*k + 131*m) & (CLN_NSL - 1))), a);
}

template <typename T>
void CleanSlotSum(const DvceArray1D<Real> &v, const T &dst, const int nc) {
  auto w = v;
  auto d = dst;
  Kokkos::parallel_for("cln_slot", Kokkos::RangePolicy<>(DevExeSpace(), 0, nc),
  KOKKOS_LAMBDA(const int n) {
    Real s = 0.0;
    for (int q = 0; q < CLN_NSL; ++q) s += w(n*CLN_NSL + q);
    d(n) += s;
  });
}

//----------------------------------------------------------------------------------------
//! \struct CleanLevel
//! \brief one level of the staggered multigrid: its blocks, vectors and their slots

struct CleanLevel {
  CleanW u;                    // blocks and vectors (m, slot, k, j, i)
  int o[3];                    // index of the first active cell along x1, x2, x3
  int n[3];                    // active cells along x1, x2, x3
  int t[3];                    // cells of a tile along x1, x2, x3 (layout_tile)
  int a0;                      // diagonal block
  int sm[3], sp[3];            // blocks of the lower and upper neighbour along each dir
};

// slots of a coarse mg level (cln_mg): x1-line coefficients, diagonal and Thomas scratch
// as CLW_AM..CLW_CP, then as CLW_TM
constexpr int CLS_A0 = 1;    // diagonal (CLS_A0 -+ 1: E_{i-1}, E_{i+1}; + 2: scratch)
constexpr int CLS_TM = 4;    // E_{j-1}, E_{j+1}, E_{k-1}, E_{k+1}
constexpr int CLS_B = 8;     // right-hand side
constexpr int CLS_X = 9;     // solution
constexpr int CLS_N = 10;

// heavy per-cell reductions at 256 threads per block (HIP's 1024 default caps VGPRs at
// 128 on gfx942 and spills); max reductions and counts are bitwise.
using CleanRp256 = Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, 1>>;

// heavy per-cell passes: 2 blocks of 256 per SM (<= 128 registers) on CUDA, 1 on HIP
#if defined(KOKKOS_ENABLE_CUDA)
static constexpr int CleanLbMb(int mb) {return mb;}
#else
static constexpr int CleanLbMb(int) {return 1;}
#endif

using CleanRpLb2 =
    Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, CleanLbMb(2)>>;

//----------------------------------------------------------------------------------------
//! \fn void CleanFatal
//! \brief stops the run with a message about an unsupported clean-scheme setup

static inline void CleanFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  // _Exit: no static destructors behind live Kokkos views
  std::_Exit(EXIT_FAILURE);
}

}  // namespace radvet

#endif  // RAD_VET_RAD_VET_ARGS_HPP_
