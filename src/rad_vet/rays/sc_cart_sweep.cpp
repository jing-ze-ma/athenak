//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_sweep.cpp
//! \brief closure = sc: variable Eddington tensor for the multi-D implicit VET solve
//! from a grey short-characteristics formal solution on the Cartesian VET grid.
//!
//! Solved per build at the lagged state:
//!     mu . grad I = chi (S - I),   chi = rho (kappa_F + kappa_s)  (opac VET_OP_T),
//!     S = eps_th a T^4 + (1 - eps_th) E,   eps_th = min(kappa_P/chi, 1)
//! (E units; vet_scatter: S = VetScatterMix).  S is a convex blend, never negative.
//! D = K/J is used in full (vet_tensor = full, realizability guarded) or projected on the
//! uniaxial form D ~ (1-chi)/2 delta + (3 chi - 1)/2 n n, chi = n.D.n in [1/3, 1], with
//! n = H/|H| of the formal solution (vet_axis = flux; eigen: the largest eigenvector,
//! which flips between axes in the nearly isotropic interior).
//! Angles: double-Gauss in mu_1 x vet_nphi azimuths about x1 (no horizontal ray), the
//! 2-D ring vet_ring, or the vet_gd set (vet_sc_set = gd); optionally rotated.
//! Sweep: layer by layer in x1, each ray traced back to the plane of cell centres of the
//! upwind layer, I, chi, S (bi)linearly interpolated there (vet_interp, PLM foot), S
//! linear in tau along the segment (vet_seg, rad_vet_vetseg.hpp).  Up and down rays run
//! in the same launches, one par_for per layer, so the moments need no atomics.
//! x1 ends: diffusion intensity E + 3 F.Omega/c (vet_bc_diff), Marshak bath / black
//! sphere (vet_bc_inner, vet_bc_bath), mirror at a reflecting end (vet_bc_reflect), or
//! the periodic sweep (vet_x1_periodic).  x2/x3: periodic, or vacuum on one MeshBlock.
//! Several MeshBlocks / ranks (VetSweepMB): the same sweep over the global layers, exact,
//! with banded planes exchanged after every launch; uniform mesh, periodic x2 (x3).
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_ring.cpp
//!  calls into: rad_vet_imex.cpp, rays/sc_cart_rotavg.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cart_sweep_mb.cpp, rays/sc_directions.cpp

#include <algorithm>
#include <cstdint>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_cart_exchange.hpp"

namespace radvet {

namespace {

// per-cell sums of one vet_sc launch: moments of the up (u = 0) and down (u = 1) rays of
// the cells iu, id and the face sums into them (f2) and into their upwind cells (fn),
// added in ray order (per-cell and ray-parallel kernels give the same bits)
struct VetScAcc {
  Real ac[2][10] = {}, hp[2][3] = {}, kp[2][VETRAY_NKP] = {}, nl[2] = {};
  Real f2[2][VETRAY_NFS] = {}, fn[2][VETRAY_NFS] = {};
  // ray r of intensity iv in a cell of source s0; g the scalars of its face sums (entry,
  // x2 and x3 interpolation flux, face crossing, exit), each times (1, mu)
  template <class A>
  KOKKOS_INLINE_FUNCTION void Add(const A &ang, const int r, const Real iv, const Real s0,
                                  const Real g[5], const bool hrs, const bool kps,
                                  const bool nls) {
    const Real m1 = ang(r,0), m2 = ang(r,1), m3 = ang(r,2), wr = ang(r,3);
    const bool up = (m1 > 0.0);
    const int u = up ? 0 : 1;
    auto add4 = [&](Real *fx, const Real q) {
      fx[0] += q;
      fx[1] += q*m1;
      fx[2] += q*m2;
      fx[3] += q*m3;
    };
    add4(f2[u] + 4, g[0]);
    add4(fn[u], g[1]);
    add4(fn[u] + 12, g[2]);
    add4(up ? f2[0] + 8 : fn[1] + 8, g[3]);
    add4(f2[u] + 4, g[4]);
    Real *c = ac[u];
    const Real a = wr*iv;
    c[0] += a;
    c[1] += a*m1*m1;
    c[2] += a*m2*m2;
    c[3] += a*m3*m3;
    c[4] += a*m1*m2;
    c[5] += a*m1*m3;
    c[6] += a*m2*m3;
    c[7] += a*m1;
    c[8] += a*m2;
    c[9] += a*m3;
    if (hrs) {
      hp[u][0] += a*fmax(m1, 0.0);
      hp[u][1] += a*fmax(m2, 0.0);
      hp[u][2] += a*fmax(m3, 0.0);
    }
    if (kps) {
      const Real mu[3] = {m1, m2, m3};
      for (int d = 0; d < 3; ++d) {
        if (!(mu[d] > 0.0)) continue;
        for (int b = 0; b < 3; ++b) kp[u][3*d+b] += a*mu[b]*mu[d];
      }
    }
    if (nls) {nl[u] += wr*fabs(iv - s0);}
  }
  template <class V>
  KOKKOS_INLINE_FUNCTION void Write(const V &vc, const int m, const int k, const int j,
                                    const int iu, const int id, const int is,
                                    const int ie, const bool mrlo, const bool mrhi,
                                    const bool hrs, const bool kps, const bool nls,
                                    const bool f2s) const {
    for (int n = 0; n < 10; ++n) {
      vc(m,VETRAY_J+n,k,j,iu) += ac[0][n];
      vc(m,VETRAY_J+n,k,j,id) += ac[1][n];
    }
    if (hrs) {
      for (int n = 0; n < 3; ++n) {
        vc(m,VETRAY_HP1+n,k,j,iu) += hp[0][n];
        vc(m,VETRAY_HP1+n,k,j,id) += hp[1][n];
      }
    }
    if (kps) {
      for (int n = 0; n < VETRAY_NKP; ++n) {
        vc(m,VETRAY_KP+n,k,j,iu) += kp[0][n];
        vc(m,VETRAY_KP+n,k,j,id) += kp[1][n];
      }
    }
    if (nls) {
      vc(m,VETRAY_NL,k,j,iu) += nl[0];
      vc(m,VETRAY_NL,k,j,id) += nl[1];
    }
    if (f2s) {
      const int ium = (iu > is) ? iu - 1 : (mrlo ? is : ie);
      const int idp = (id < ie) ? id + 1 : (mrhi ? ie : is);
      for (int n = 0; n < VETRAY_NFS; ++n) {
        vc(m,VETRAY_F2+n,k,j,iu) += f2[0][n];
        vc(m,VETRAY_F2+n,k,j,id) += f2[1][n];
        vc(m,VETRAY_F2+n,k,j,ium) += fn[0][n];
        vc(m,VETRAY_F2+n,k,j,idp) += fn[1][n];
      }
    }
  }
};

// periodic index in [0, n)
KOKKOS_INLINE_FUNCTION
int VetWrap(const int a, const int n) {
  int r = a % n;
  return (r < 0) ? (r + n) : r;
}

// largest eigenvalue and its unit eigenvector of a symmetric 3x3 matrix (cyclic Jacobi)
KOKKOS_INLINE_FUNCTION
void VetEigMax(Real a[3][3], Real &lmax, Real &e1, Real &e2, Real &e3) {
  Real v[3][3];
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {v[p][q] = (p == q) ? 1.0 : 0.0;}
  }
  for (int sw = 0; sw < 6; ++sw) {
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Real apq = a[p][q];
        if (apq == 0.0) continue;
        Real th = 0.5*(a[q][q] - a[p][p])/apq;
        Real t = ((th >= 0.0) ? 1.0 : -1.0)/(fabs(th) + sqrt(th*th + 1.0));
        Real c = 1.0/sqrt(t*t + 1.0), s = t*c;
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[kk][p], xq = a[kk][q];
          a[kk][p] = c*xp - s*xq;
          a[kk][q] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[p][kk], xq = a[q][kk];
          a[p][kk] = c*xp - s*xq;
          a[q][kk] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = v[kk][p], xq = v[kk][q];
          v[kk][p] = c*xp - s*xq;
          v[kk][q] = s*xp + c*xq;
        }
      }
    }
  }
  int im = 0;
  if (a[1][1] > a[im][im]) im = 1;
  if (a[2][2] > a[im][im]) im = 2;
  lmax = a[im][im];
  e1 = v[0][im];
  e2 = v[1][im];
  e3 = v[2][im];
}

// all eigenvalues lam[p] and unit eigenvectors v[.][p] of a symmetric 3x3 matrix (cyclic
// Jacobi run to round-off); a is destroyed
KOKKOS_INLINE_FUNCTION
void VetEigSym(Real a[3][3], Real lam[3], Real v[3][3]) {
  for (int p = 0; p < 3; ++p) {
    for (int q = 0; q < 3; ++q) {v[p][q] = (p == q) ? 1.0 : 0.0;}
  }
  for (int sw = 0; sw < 12; ++sw) {
    Real off = fabs(a[0][1]) + fabs(a[0][2]) + fabs(a[1][2]);
    if (off <= 1.0e-18*(fabs(a[0][0]) + fabs(a[1][1]) + fabs(a[2][2]))) break;
    for (int p = 0; p < 2; ++p) {
      for (int q = p + 1; q < 3; ++q) {
        Real apq = a[p][q];
        if (apq == 0.0) continue;
        Real th = 0.5*(a[q][q] - a[p][p])/apq;
        Real t = ((th >= 0.0) ? 1.0 : -1.0)/(fabs(th) + sqrt(th*th + 1.0));
        Real c = 1.0/sqrt(t*t + 1.0), s = t*c;
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[kk][p], xq = a[kk][q];
          a[kk][p] = c*xp - s*xq;
          a[kk][q] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = a[p][kk], xq = a[q][kk];
          a[p][kk] = c*xp - s*xq;
          a[q][kk] = s*xp + c*xq;
        }
        for (int kk = 0; kk < 3; ++kk) {
          Real xp = v[kk][p], xq = v[kk][q];
          v[kk][p] = c*xp - s*xq;
          v[kk][q] = s*xp + c*xq;
        }
      }
    }
  }
  for (int p = 0; p < 3; ++p) {lam[p] = a[p][p];}
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetShortChar
//! \brief vet_sc build, with the vet_rotate_avg bookkeeping around VetShortCharSweep

void RadiationVet::VetShortChar() {
  // vet_rotate_avg: the first build of a new block shifts the ring by the old block's
  // Delta (else re-init: sweep the K-1 previous sets); the last build of a block also
  // sweeps the next set at the same state (Delta) and restores vet_cell
  const bool ra = vra_on && vsc_rot;
  const int ilo = pmy_pack->pmesh->mb_indcs.is;
  const int64_t b = VetGdBlock(pmy_pack->pmesh->ncycle);
  if (ra && vra_state == 1 && b != vra_blk) {
    if (vra_pend == 3) {
      VetScRotAvgSetup();
      RotAvgShift(vra_k - 1);
      vra_pend = 0;
      vra_nrot += 1.0;
    } else {
      vra_state = 0;
    }
  }
  if (ra && vra_state == 0) {
    Kokkos::Timer tr;
    vra_head = 0;
    for (int jb = 1; jb < vra_k; ++jb) {
      VetRotSet(VetRotBeta(b - jb));
      vra_raw = true;
      VetShortCharSweep();
      vra_raw = false;
      VetScRotAvgSetup();
      RotAvgSnap(jb - 1, ilo);
      vra_nsweep += 1.0;
    }
    vra_pend = 1;
    vra_state = 1;
    vra_ninit += 1.0;
    vra_time += tr.seconds();
  }
  VetRotUpdate();   // set of this step's rotation block
  vra_last = ra && vra_state == 1 && RotLastStep();
  VetShortCharSweep();
  if (vra_last) {
    Kokkos::Timer tr;
    if (vra_bk_vc.extent(0) != vet_cell.extent(0) ||
        vra_bk_vc.extent(1) != vet_cell.extent(1)) {
      Kokkos::realloc(vra_bk_vc, vet_cell.extent(0), vet_cell.extent(1),
                      vet_cell.extent(2), vet_cell.extent(3), vet_cell.extent(4));
    }
    Kokkos::deep_copy(vra_bk_vc, vet_cell);
    VetRotSet(VetRotBeta(b + 1));
    vra_raw = true;
    VetShortCharSweep();
    vra_raw = false;
    VetScRotAvgSetup();
    RotAvgDelta(vra_k - 1, ilo);
    Kokkos::deep_copy(vet_cell, vra_bk_vc);
    vra_pend = 3;
    vra_nsweep += 1.0;
    vra_time += tr.seconds();
  }
  vra_last = false;
  vra_blk = b;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetShortCharSweep
//! \brief the vet_sc build with the current set

void RadiationVet::VetShortCharSweep() {
  // vet_source_noesrc: source and boundary intensity read the physical E^n = u0, not
  // the solve's old vector (VET_IW_EN)
  const bool noes_ = VsrcNoes();
  auto esv_ = noes_ ? u0 : iw;
  const int esc_ = noes_ ? VET_E : VET_IW_EN;
  // half-range sums sum w I max(mu_a, 0) (VETRAY_HP1..3)
  const bool hrs_ = impl_beam_hr;
  // half-range stress columns K+^d_bd (VETRAY_KP, boundary fluxes)
  const bool kps_ = true;
  // implicit_thin_switch = fs | fskn: non-locality sum w |I - S| (VETRAY_NL, C1)
  const bool nls_ = thsw_c1;
  // clean_consistency: the rays' face fluxes (VETRAY_F2..): x1 flux where the segment
  // crosses the face, the transverse foot interpolation as an x2 (x3) flux, so the cell
  // balance of the face fluxes is the exchange inside the cell
  const bool f2s_ = cln_consist;
  const bool x1per = vet_x1per;
  Kokkos::fence();
  Kokkos::Timer timer;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool thrd = trans_x3;
  const int nray = vet_nray;
  auto iw_ = iw;
  auto vc_ = vet_cell;
  auto ip_ = vet_ipl;
  auto ang_ = vet_ang;
  auto &mbsize = pmy_pack->pmb->mb_size;
  const Real cl = c_light;
  const Real ar = arad;
  const Real efl = e_floor;
  const bool thermal = fl_on && coupling && !opac_zero;
  const bool vsf = vscat;   // vet_scatter
  const Real bblo = VetBathLo();   // vet_bc_inner
  const Real bbhi = (vet_bc_bath && ibc_x1max == VET_IBC_MARSHAK) ? iebath_x1max : 0.0;
  const bool bfc = vet_bc_face;
  const bool x1src = vet_x1per && vet_x1src;
  constexpr bool spl = true;   // (the segment rule: split, emission and jump terms)
  constexpr int spe = 3;
  const bool x2o = vet_x2open, x3o = vet_x3open;   // vacuum x2 (x3) faces
  constexpr int vim = 3;   // (corner at a material edge, log elsewhere)
  const bool plm_ = true;
  // boundary beam rays (ImplBeamInit)
  const int bmn_ = ibm_nr;
  const int bmr0 = ibm_ray[0], bmr1 = ibm_ray[1], bmr2 = ibm_ray[2], bmr3 = ibm_ray[3];
  auto bmi_ = ibm_i;

  // (1) extinction and source, cell by cell
  auto opd_ = opac;
  auto iwd_ = iw;
  auto vcd_ = vet_cell;
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  VetFor("vet_vet_srcc", (nmb1 + 1)*nk*nj*ni, KOKKOS_LAMBDA(const int t) {
    const int i = is + t%ni;
    int q = t/ni;
    const int j = js + q%nj;
    q /= nj;
    const int k = ks + q%nk, m = q/nk;
    Real chx = fmax(opd_(m,VET_OP_T,k,j,i), 1.0e-300);
    Real e = fmax(esv_(m,esc_,k,j,i), efl);
    Real s = e;
    if (thermal) {
      Real tg = iwd_(m,VET_IW_TP,k,j,i);
      Real t2 = tg*tg;
      if (vsf) {   // vet_scatter (VetScatterMix, rad_vet.hpp)
        const Real rkp = fmin(opd_(m,VET_OP_P,k,j,i), chx);
        s = VetScatterMix(rkp, chx, ar*t2*t2, e);
      } else {
        Real eth = fmin(opd_(m,VET_OP_P,k,j,i)/chx, 1.0);
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
    }
    vcd_(m,VETRAY_CHX,k,j,i) = chx;
    vcd_(m,VETRAY_SRC,k,j,i) = s;
    for (int n = VETRAY_J; n < VETRAY_CHI; ++n) {
      vcd_(m,n,k,j,i) = 0.0;
    }
    if (hrs_) {
      for (int n = 0; n < 3; ++n) {vcd_(m,VETRAY_HP1+n,k,j,i) = 0.0;}
    }
    if (kps_) {
      for (int n = 0; n < VETRAY_NKP; ++n) {vcd_(m,VETRAY_KP+n,k,j,i) = 0.0;}
    }
    if (f2s_) {
      for (int n = 0; n < VETRAY_NFS; ++n) {vcd_(m,VETRAY_F2+n,k,j,i) = 0.0;}
    }
    if (nls_) {
      vcd_(m,VETRAY_NL,k,j,i) = 0.0;
      // vacuum: no lagged scattering source, eps = 1
      const Real eps = thermal ? fmin(opd_(m,VET_OP_P,k,j,i)/chx, 1.0) : 0.0;
      vcd_(m,VETRAY_EPS,k,j,i) = !(opd_(m,VET_OP_T,k,j,i) > 0.0) ? 1.0 : eps;
    }
  });

  // (2) the sweep: launch l does layer is+l for the up rays and ie-l for the down rays,
  // reading launch l-1's plane; several MeshBlocks: VetSweepMB instead
  if (vet_mbs != nullptr) {
    VetSweepMB();
  }
  const int nlaunch = (vet_mbs != nullptr) ? 0 : nx1;   // banded sweep ran instead
  // vet_bc_reflect = mirror: at a reflecting x1 end the entering rays carry the leaving
  // intensity of their x1 mirror ray (r +- nray/2), read in a further pass from the
  // boundary layer's plane.  One reflecting end: exact in 2 passes; both (a closed
  // cavity): vet_x1_npass passes, each adding one round trip.
  const bool mlo = vet_bc_mirror && (ibc_x1min == VET_IBC_REFLECT) && !vet_x1per;
  const bool mhi = vet_bc_mirror_hi && (ibc_x1max == VET_IBC_REFLECT) && !vet_x1per;
  if ((mlo || mhi) && vet_mbs != nullptr) {
    VetFatal("<rad_vet>/vet_bc_reflect = mirror (the default with a reflecting x1 end) "
             "needs a single MeshBlock on vet_sc (the banded sweep has no mirror); name "
             "vet_bc_reflect = open");
  }
  const int nhray = nray/2;
  const int npass = vet_x1per ? vet_x1npass : ((mlo && mhi) ? vet_mirror_npass
                                               : ((mlo || mhi) ? 2 : 1));
  // vet_x1_ptol > 0: adaptive periodic passes; each pass accumulates into zeroed moments
  // and its exit plane is compared with the previous one (one reduction per pass)
  const bool adapt = vet_x1per && (vet_x1ptol > 0.0) && (nlaunch > 0);
  Real dprev = -1.0;
  int npdone = 0;
  bool conv = false;
  // an adaptive sweep after the first build starts from the last build's exit plane (its
  // periodic image then), not from the boundary inflow
  const bool warm = adapt && vet_x1_warm;
  if (warm) {
    auto ip_ = vet_ipl;
    auto ent_ = vet_x1ent;
    par_for("vet_vet_x1warm", DevExeSpace(), 0, nmb1, 0, nray - 1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int r, const int k, const int j) {
      ip_(m,1,r,k,j) = ent_(m,r,k,j);
    });
  }
  for (int p = 0; p < npass; ++p) {
  const bool wrap = vet_x1per && (p > 0 || warm);   // layer is (ie) reads the image
  const bool mrlo = mlo && (p > 0), mrhi = mhi && (p > 0);   // or the mirror image
  // only the last pass adds to the moments (adaptive: each pass, on zeroed moments)
  const bool accum = adapt || (p == npass - 1);
  if (adapt && p > 0) {
    par_for("vet_vet_x1zero", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int n = VETRAY_J; n < VETRAY_CHI; ++n) {vc_(m,n,k,j,i) = 0.0;}
      if (hrs_) {for (int n = 0; n < 3; ++n) {vc_(m,VETRAY_HP1+n,k,j,i) = 0.0;}}
      if (kps_) {for (int n = 0; n < VETRAY_NKP; ++n) {vc_(m,VETRAY_KP+n,k,j,i) = 0.0;}}
      if (f2s_) {for (int n = 0; n < VETRAY_NFS; ++n) {vc_(m,VETRAY_F2+n,k,j,i) = 0.0;}}
      if (nls_) {vc_(m,VETRAY_NL,k,j,i) = 0.0;}
    });
  }
  for (int l = 0; l < nlaunch; ++l) {
    const int gl = p*nlaunch + l;
    const int pw = gl & 1, pr = pw ^ 1;
    const int iu = is + l, id = ie - l;
    // ray r of launch l in cell (k, j): its intensity (into the plane buffer) and the
    // scalars g of its face sums, shared by the per-cell and the ray-parallel kernel
    auto rayf = KOKKOS_LAMBDA(const int m, const int k, const int j, const int r,
                              Real g[5]) {
      const Real dx1 = mbsize.d_view(m).dx1;
      const Real dx2 = mbsize.d_view(m).dx2;
      const Real dx3 = mbsize.d_view(m).dx3;
      const Real m1 = ang_(r,0), m2 = ang_(r,1), m3 = ang_(r,2), wr = ang_(r,3);
      const bool up = (m1 > 0.0);
      const int i = up ? iu : id;
      const Real am1 = fabs(m1);
      const Real c0 = vc_(m,VETRAY_CHX,k,j,i);
      const Real s0 = vc_(m,VETRAY_SRC,k,j,i);
      Real iv;
      const bool mimg = (l == 0) && (up ? mrlo : mrhi);   // vet_bc_reflect = mirror
      for (int q = 0; q < 5; ++q) {g[q] = 0.0;}
      Real ifc = 0.0;   // clean_consistency: intensity at the face crossing
      if (l == 0 && !wrap && !mimg) {
        // half a cell from the boundary face; face S linearly extrapolated from the two
        // boundary cells (exact for a deep linear S), linear in tau along the segment
        const int iin = up ? (i + 1) : (i - 1);
        const Real s1 = vc_(m,VETRAY_SRC,k,j,iin);
        const Real su = fmax(1.5*s0 - 0.5*s1, 0.0);
        Real ib = up ? 0.0 : bbhi;
        if (x1src) {
          ib = su;   // vet_x1_periodic_inflow = source: symmetric entry
        } else if (up) {
          if (bblo >= 0.0) {
            ib = bblo;
          } else {
            // diffusion: eps = E + 3 (F . Omega)/c of the bottom cell; vet_bc_diff = face
            // moves E to the face, E_f = E + (S_f - S_0)
            Real e = esv_(m,esc_,k,j,i);
            if (bfc) e += 0.5*(s0 - s1);
            e = fmax(e, efl);
            Real fo = iw_(m,VET_IW_F1,k,j,i)*m1 + iw_(m,VET_IW_F2,k,j,i)*m2
                      + iw_(m,VET_IW_F3,k,j,i)*m3;
            ib = fmax(e + 3.0*fo/cl, 0.0);
          }
          if (bmn_ > 0) {   // boundary beams
            const int q = (r == bmr0) ? 0 : ((r == bmr1) ? 1 : ((r == bmr2) ? 2
                          : ((r == bmr3) ? 3 : -1)));
            if (q >= 0 && q < bmn_) {ib += bmi_(m,q,k-ks,j-js);}
          }
        }
        if (f2s_) g[0] = wr*m1*ib;
        const Real dtau = 0.5*c0*dx1/am1;
        Real ex, w0, wu;
        VetSegW(dtau, ex, w0, wu);
        iv = ib*ex + wu*su + w0*s0;
      } else {
        const int iup = mimg ? i : ((l == 0) ? (up ? ie : is) : (up ? (i - 1) : (i + 1)));
        const int rb = mimg ? ((r < nhray) ? (r + nhray) : (r - nhray)) : r;
        // foot of the characteristic on the upwind plane, in cell units
        const Real sh2 = -dx1*m2/(am1*dx2);
        const Real fl2 = floor(sh2);
        const Real a2 = sh2 - fl2;
        const int o2 = static_cast<int>(fl2);
        int ja = js + VetWrap(j - js + o2, nx2);
        int jb = js + VetWrap(j - js + o2 + 1, nx2);
        // vacuum x2 (x3) faces: a foot outside the box carries I = 0 and the edge cell's
        // chi and S
        Real ga = 1.0, gb = 1.0;
        if (x2o) {
          const int qa = j - js + o2, qb = qa + 1;
          ga = (qa >= 0 && qa < nx2) ? 1.0 : 0.0;
          gb = (qb >= 0 && qb < nx2) ? 1.0 : 0.0;
          ja = js + ((qa < 0) ? 0 : ((qa >= nx2) ? (nx2 - 1) : qa));
          jb = js + ((qb < 0) ? 0 : ((qb >= nx2) ? (nx2 - 1) : qb));
        }
        int ka = k, kb = k, o3 = 0;
        Real a3 = 0.0, ka_g = 1.0, kb_g = 1.0;
        const Real sh3 = thrd ? -dx1*m3/(am1*dx3) : 0.0;
        if (thrd) {
          const Real fl3 = floor(sh3);
          a3 = sh3 - fl3;
          o3 = static_cast<int>(fl3);
          ka = ks + VetWrap(k - ks + o3, nx3);
          kb = ks + VetWrap(k - ks + o3 + 1, nx3);
          if (x3o) {
            const int qa = k - ks + o3, qb = qa + 1;
            ka_g = (qa >= 0 && qa < nx3) ? 1.0 : 0.0;
            kb_g = (qb >= 0 && qb < nx3) ? 1.0 : 0.0;
            ka = ks + ((qa < 0) ? 0 : ((qa >= nx3) ? (nx3 - 1) : qa));
            kb = ks + ((qb < 0) ? 0 : ((qb >= nx3) ? (nx3 - 1) : qb));
          }
        }
        // upwind plane at offsets (dk, dj) from the cell, 0 beyond a vacuum face
        auto ipq = [&](const int dk, const int dj) {
          const int qj = j - js + dj, qk = k - ks + dk;
          if ((x2o && (qj < 0 || qj >= nx2)) || (x3o && (qk < 0 || qk >= nx3))) {
            return 0.0;
          }
          return ip_(m,pr,rb,thrd ? ks + VetWrap(qk, nx3) : k,js + VetWrap(qj, nx2));
        };
        Real dv = 0.0, fl[2] = {0.0, 0.0};
        if (plm_) {
          dv = VetPlmFoot([&](const int dk, const int dj) {
            return ipq(o3 + dk, o2 + dj);
          }, a2, a3, thrd, fl);
        }
        if (f2s_) {
          // Phi: the upwind plane's intensity carried across the low x2 (x3) face by the
          // foot shift.  The x2 and x3 parts S2, S3 commute: the x2 face carries
          // (I + S3 I)/2, the x3 face (I + S2 I)/2, their balance is S3 S2 I - I.
          auto g2 = [&](const int d) {
            return thrd ? 0.5*(ipq(0,d) + (1.0 - a3)*ipq(o3,d) + a3*ipq(o3+1,d))
                        : ipq(0,d);
          };
          auto g3 = [&](const int d) {
            return 0.5*(ipq(d,0) + (1.0 - a2)*ipq(d,o2) + a2*ipq(d,o2+1));
          };
          auto phi = [&](const Real sh, const auto &g) {
            const Real s = fabs(sh), f = s - floor(s);
            const int n = static_cast<int>(floor(s)), sg = (sh < 0.0) ? 1 : -1;
            Real p = f*g(sg > 0 ? -1 - n : n);
            for (int q = 0; q < n; ++q) p += g(sg > 0 ? -1 - q : q);
            return sg*p;
          };
          // mirror: the upwind plane is the wall column's image, so the shift crosses
          // its own x2 faces and the wall face carries no net flux
          g[1] = wr*am1*dx2/dx1*(phi(sh2, g2) + fl[0]);
          if (thrd) {
            g[2] = wr*am1*dx3/dx1*(phi(sh3, g3) + fl[1]);
          }
        }
        const Real waa = (1.0 - a2)*(1.0 - a3), wab = a2*(1.0 - a3);
        const Real wba = (1.0 - a2)*a3, wbb = a2*a3;
        // vacuum-face gates on the mirror-aware ray rb; exactly 1 off the open faces
        const Real gaa = ga*ka_g, gab = gb*ka_g, gba = ga*kb_g, gbb = gb*kb_g;
        Real iup_v;
        if (x2o || x3o) {
          iup_v = waa*gaa*ip_(m,pr,rb,ka,ja) + wab*gab*ip_(m,pr,rb,ka,jb)
                  + wba*gba*ip_(m,pr,rb,kb,ja) + wbb*gbb*ip_(m,pr,rb,kb,jb);
        } else {
          iup_v = waa*ip_(m,pr,rb,ka,ja) + wab*ip_(m,pr,rb,ka,jb)
                  + wba*ip_(m,pr,rb,kb,ja) + wbb*ip_(m,pr,rb,kb,jb);
        }
        iup_v += dv;
        const Real cu = waa*vc_(m,VETRAY_CHX,ka,ja,iup)
                        + wab*vc_(m,VETRAY_CHX,ka,jb,iup)
                        + wba*vc_(m,VETRAY_CHX,kb,ja,iup)
                        + wbb*vc_(m,VETRAY_CHX,kb,jb,iup);
        const Real su = waa*vc_(m,VETRAY_SRC,ka,ja,iup)
                        + wab*vc_(m,VETRAY_SRC,ka,jb,iup)
                        + wba*vc_(m,VETRAY_SRC,kb,ja,iup)
                        + wbb*vc_(m,VETRAY_SRC,kb,jb,iup);
        // jump: the corners' S one x1 cell further upwind, its periodic or mirror image
        // beyond an x1 end; none (q = 1) beyond an open one
        int iuu = up ? (i - 2) : (i + 2);
        if (x1per) {
          iuu = is + VetWrap(iuu - is, nx1);
        } else if (up ? (iuu < is && mrlo) : (iuu > ie && mrhi)) {
          iuu = up ? (2*is - 1 - iuu) : (2*ie + 1 - iuu);
        }
        const bool juu = spl && spe == 2 && iuu >= is && iuu <= ie;
        auto s2 = [&](const int kq, const int jq) {
          return juu ? vc_(m,VETRAY_SRC,kq,jq,iuu) : 0.0;
        };
        const Real squ[4] = {s2(ka,ja), s2(ka,jb), s2(kb,ja), s2(kb,jb)};
        const Real *suq = juu ? squ : nullptr;
        // and the own cell's S one x1 cell downstream (active cells only)
        int idn = up ? (i + 1) : (i - 1);
        if (x1per) idn = is + VetWrap(idn - is, nx1);
        const Real sdd = (spl && spe == 2 && idn >= is && idn <= ie)
                         ? vc_(m,VETRAY_SRC,k,j,idn) : -1.0;
        if (vim != 0) {   // vet_interp: log | corner
          const Real wq[4] = {waa, wab, wba, wbb};
          // a corner outside a vacuum face carries I = 0 and the edge cell's chi and S
          const Real iq[4] = {gaa*ip_(m,pr,rb,ka,ja) + dv, gab*ip_(m,pr,rb,ka,jb) + dv,
                              gba*ip_(m,pr,rb,kb,ja) + dv, gbb*ip_(m,pr,rb,kb,jb) + dv};
          const Real cq[4] = {vc_(m,VETRAY_CHX,ka,ja,iup), vc_(m,VETRAY_CHX,ka,jb,iup),
                              vc_(m,VETRAY_CHX,kb,ja,iup), vc_(m,VETRAY_CHX,kb,jb,iup)};
          const Real sq[4] = {vc_(m,VETRAY_SRC,ka,ja,iup), vc_(m,VETRAY_SRC,ka,jb,iup),
                              vc_(m,VETRAY_SRC,kb,ja,iup), vc_(m,VETRAY_SRC,kb,jb,iup)};
          iv = VetUpCorner(vim, wq, iq, cq, sq, s0, c0, dx1/am1, spl, spe,
                           f2s_ ? &ifc : nullptr, suq, sdd);
        } else if (spl) {
          // split: the two halves at the x1 face
          iv = VetSegSplit(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1, spe,
                           f2s_ ? &ifc : nullptr,
                           juu ? (waa*squ[0] + wab*squ[1] + wba*squ[2] + wbb*squ[3])
                               : -1.0, sdd);
        } else {
          const Real dtau = 0.5*(cu + c0)*dx1/am1;
          Real ex, w0, wu;
          VetSegW(dtau, ex, w0, wu);
          iv = iup_v*ex + wu*su + w0*s0;
          // face crossing at the face (0.5 chi_u hl from u), not at half the tau
          if (f2s_) ifc = VetSegFaceI(iup_v, su, s0, 0.5*cu*dx1/am1, 0.5*c0*dx1/am1,
                                      0);
        }
        // face crossing of the single segment: a face flux
        if (f2s_ && !mimg) g[3] = wr*m1*fmax(ifc, 0.0);
      }
      iv = fmax(iv, 0.0);
      ip_(m,pw,r,k,j) = iv;
      if (f2s_ && !x1per && i == (up ? ie : is)) {
        // ray leaving through the x1 boundary face: the half cell at chi_0, S_0
        const Real ex = exp(-0.5*c0*dx1/am1);
        g[4] = wr*m1*(iv*ex + (1.0 - ex)*s0);
      }
      return iv;
    };
    if (vsc_rayk) {
      auto rg_ = vet_rayg;
      par_for("vet_vet_sweep_ray", DevExeSpace(), 0, nmb1, 0, nray - 1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int r, const int k, const int j) {
        Real g[5];
        rayf(m, k, j, r, g);
        for (int q = 0; q < 5; ++q) {rg_(m,q*nray + r,k,j) = g[q];}
      });
      par_for("vet_vet_sweep_sum", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        VetScAcc a;
        for (int r = 0; r < nray; ++r) {
          const Real g[5] = {rg_(m,r,k,j), rg_(m,nray + r,k,j), rg_(m,2*nray + r,k,j),
                             rg_(m,3*nray + r,k,j), rg_(m,4*nray + r,k,j)};
          const int i = (ang_(r,0) > 0.0) ? iu : id;
          a.Add(ang_, r, ip_(m,pw,r,k,j), vc_(m,VETRAY_SRC,k,j,i), g, hrs_, kps_, nls_);
        }
        if (accum) {
          a.Write(vc_, m, k, j, iu, id, is, ie, mrlo, mrhi, hrs_, kps_, nls_, f2s_);
        }
      });
    } else {
      par_for("vet_vet_sweep", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        VetScAcc a;
        for (int r = 0; r < nray; ++r) {
          Real g[5];
          const Real iv = rayf(m, k, j, r, g);
          const int i = (ang_(r,0) > 0.0) ? iu : id;
          a.Add(ang_, r, iv, vc_(m,VETRAY_SRC,k,j,i), g, hrs_, kps_, nls_);
        }
        if (accum) {
          a.Write(vc_, m, k, j, iu, id, is, ie, mrlo, mrhi, hrs_, kps_, nls_, f2s_);
        }
      });
    }
  }
  ++npdone;
  if (adapt) {
    // exit plane of this pass (= next pass's entry): up rays at ie, down rays at is;
    // d = max|I_p - I_{p-1}|, a = max I
    const int pl = (p*nlaunch + nlaunch - 1) & 1;
    auto ent_ = vet_x1ent;
    const int nj = je - js + 1, nk = ke - ks + 1;
    const int ntot = (nmb1 + 1)*nray*nk*nj;
    const bool cmp = (p > 0);
    Real dmax = 0.0, amax = 0.0;
    Kokkos::parallel_reduce("vet_vet_x1conv",
      Kokkos::RangePolicy<DevExeSpace>(0, ntot),
      KOKKOS_LAMBDA(const int t, Real &dl, Real &al) {
        int q = t;
        const int j = js + q % nj; q /= nj;
        const int k = ks + q % nk; q /= nk;
        const int r = q % nray;
        const int m = q / nray;
        const Real v = ip_(m,pl,r,k,j);
        if (cmp) {dl = fmax(dl, fabs(v - ent_(m,r,k,j)));}
        al = fmax(al, v);
        ent_(m,r,k,j) = v;
      }, Kokkos::Max<Real>(dmax), Kokkos::Max<Real>(amax));
    if (p > 0) {
      const Real d = (amax > 0.0) ? (dmax/amax) : 0.0;
      // error of the entry this pass used: d/(1 - q)
      const Real q = (dprev > 0.0) ? fmin(d/dprev, 0.999) : 0.999;
      conv = (d == 0.0) || (d/(1.0 - q) < vet_x1ptol);
      dprev = d;
      if (conv) {break;}
    }
  }
  }
  if (adapt) {
    vet_x1_nbuild += 1.0;
    vet_x1_warm = true;
    vet_x1_npsum += static_cast<Real>(npdone);
    vet_x1_npmax = std::max(vet_x1_npmax, npdone);
    if (!conv) {
      if (vet_x1_ncap == 0.0 && global_variable::my_rank == 0) {
        std::cout << "<rad_vet> WARNING vet_sc periodic x1 sweep: vet_x1_npass="
                  << vet_x1npass << " reached before vet_x1_ptol=" << vet_x1ptol
                  << " (first time; counted in the final report)" << std::endl;
      }
      vet_x1_ncap += 1.0;
    }
  }

  // vet_rotate_avg: pending offsets from this sweep, then the mean over the last K sets
  // (not on a snapshot)
  if (vra_on && vsc_rot && !vra_raw && vra_state == 1) {VetScRotAvgSetup();}
  // a history loaded from the restart file must fit the table
  if (vra_on && vsc_rot && !vra_raw && vra_state == 1 && RotAvgRstCheck()) {
    Kokkos::Timer tr;
    if (vra_last) {RotAvgSnap(vra_k - 1, is);}   // current set's raw moments (Delta)
    RotAvgUpdate(is);
    RotAvgApply(is);
    Kokkos::fence();
    vra_time += tr.seconds();
  }

  // (3) the uniaxial projection (chi, n) of D = K/J
  const bool axflux = vet_axis_flux;
  par_for("vet_vet_proj", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real jj = vc_(m,VETRAY_J,k,j,i);
    Real chi = 1.0/3.0, n1 = 1.0, n2 = 0.0, n3 = 0.0;
    if (jj > 0.0) {
      Real ij = 1.0/jj;
      Real a[3][3];
      a[0][0] = vc_(m,VETRAY_K11,k,j,i)*ij;
      a[1][1] = vc_(m,VETRAY_K11+1,k,j,i)*ij;
      a[2][2] = vc_(m,VETRAY_K11+2,k,j,i)*ij;
      a[0][1] = a[1][0] = vc_(m,VETRAY_K11+3,k,j,i)*ij;
      a[0][2] = a[2][0] = vc_(m,VETRAY_K11+4,k,j,i)*ij;
      a[1][2] = a[2][1] = vc_(m,VETRAY_K11+5,k,j,i)*ij;
      // axis: the formal solution's own flux direction (vet_axis = flux) or the
      // principal axis of D (eigen)
      Real h1 = vc_(m,VETRAY_H1,k,j,i), h2 = vc_(m,VETRAY_H1+1,k,j,i);
      Real h3 = vc_(m,VETRAY_H1+2,k,j,i);
      Real hm = sqrt(h1*h1 + h2*h2 + h3*h3);
      if (axflux && hm > 1.0e-8*jj) {
        n1 = h1/hm;
        n2 = h2/hm;
        n3 = h3/hm;
        chi = n1*(a[0][0]*n1 + a[0][1]*n2 + a[0][2]*n3)
              + n2*(a[1][0]*n1 + a[1][1]*n2 + a[1][2]*n3)
              + n3*(a[2][0]*n1 + a[2][1]*n2 + a[2][2]*n3);
      } else {
        VetEigMax(a, chi, n1, n2, n3);
      }
      Real nn = sqrt(n1*n1 + n2*n2 + n3*n3);
      if (nn > 0.0) {
        n1 /= nn;
        n2 /= nn;
        n3 /= nn;
      } else {
        n1 = 1.0;
        n2 = 0.0;
        n3 = 0.0;
      }
      chi = fmin(fmax(chi, 1.0/3.0), 1.0);
    }
    vc_(m,VETRAY_CHI,k,j,i) = chi;
    vc_(m,VETRAY_N1,k,j,i) = n1;
    vc_(m,VETRAY_N1+1,k,j,i) = n2;
    vc_(m,VETRAY_N1+2,k,j,i) = n3;
  });

  // (4) vet_tensor = full: the guarded D = K/J the solve reads instead of (chi, n)
  if (vet_full) {VetFullTensor();}
  Kokkos::fence();
  vet_time += timer.seconds();
  vet_ncall += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetFullTensor
//! \brief vet_tensor = full: D_ab = K_ab/J of the last formal solution, realizability
//! guarded, into vet_cell(VETRAY_D11 .. +5), active cells and ghosts.
//! With w_r > 0 and I_r >= 0 the raw D is symmetric, trace 1 to round-off, eigenvalues
//! in [0, 1].  The guard makes that exact:
//!   (1) J <= 0 or non-positive trace: D = delta/3 (counted);
//!   (2) D /= tr D;
//!   (3) eigenvalues below emin = vet_eig_min: lam -> emin + (lam' - emin)
//!       (1 - 3 emin)/(sum lam' - 3 emin), lam' = max(lam, emin) (counted).
//! So D_dd >= emin >= 0: the line rows stay M-matrices.  |D_guarded - D_raw| goes to
//! VETRAY_GD for the run statistics.

void RadiationVet::VetFullTensor() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n1 = indcs.nx1 + 2*(indcs.ng);
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*(indcs.ng)) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*(indcs.ng)) : 1;
  const bool thrd = trans_x3;
  const Real emin = vet_eig_min;
  const bool x2o = vet_x2open, x3o = vet_x3open;
  auto vc_ = vet_cell;

  par_for("vet_vet_full", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jj = vc_(m,VETRAY_J,k,j,i);
    Real d[6];
    Real raw[6];
    for (int n = 0; n < 6; ++n) {raw[n] = 0.0;}
    Real tr = 0.0;
    if (jj > 0.0) {
      const Real ij = 1.0/jj;
      for (int n = 0; n < 6; ++n) {raw[n] = vc_(m,VETRAY_K11+n,k,j,i)*ij;}
      tr = raw[0] + raw[1] + raw[2];
    }
    bool hit = false;
    if (!(jj > 0.0) || !(tr > 0.0)) {
      d[0] = d[1] = d[2] = 1.0/3.0;
      d[3] = d[4] = d[5] = 0.0;
      hit = true;
    } else {
      const Real itr = 1.0/tr;
      for (int n = 0; n < 6; ++n) {d[n] = raw[n]*itr;}
      Real a[3][3], v[3][3], lam[3];
      a[0][0] = d[0];
      a[1][1] = d[1];
      a[2][2] = d[2];
      a[0][1] = a[1][0] = d[3];
      a[0][2] = a[2][0] = d[4];
      a[1][2] = a[2][1] = d[5];
      VetEigSym(a, lam, v);
      if (lam[0] < emin || lam[1] < emin || lam[2] < emin) {
        Real sl = 0.0;
        for (int p = 0; p < 3; ++p) {
          lam[p] = fmax(lam[p], emin);
          sl += lam[p];
        }
        const Real sc = (1.0 - 3.0*emin)/fmax(sl - 3.0*emin, 1.0e-300);
        for (int p = 0; p < 3; ++p) {lam[p] = emin + (lam[p] - emin)*sc;}
        Real b[3][3];
        for (int p = 0; p < 3; ++p) {
          for (int q = 0; q < 3; ++q) {
            b[p][q] = lam[0]*v[p][0]*v[q][0] + lam[1]*v[p][1]*v[q][1]
                      + lam[2]*v[p][2]*v[q][2];
          }
        }
        d[0] = b[0][0];
        d[1] = b[1][1];
        d[2] = b[2][2];
        d[3] = 0.5*(b[0][1] + b[1][0]);
        d[4] = 0.5*(b[0][2] + b[2][0]);
        d[5] = 0.5*(b[1][2] + b[2][1]);
        hit = true;
      }
    }
    Real dev = 0.0;
    for (int n = 0; n < 6; ++n) {
      vc_(m,VETRAY_D11+n,k,j,i) = d[n];
      dev = fmax(dev, fabs(d[n] - raw[n]));
    }
    // <= -1 marks a cell the guard left alone (trace normalisation only): -1 - dev
    vc_(m,VETRAY_GD,k,j,i) = hit ? dev : (-1.0 - dev);
  });

  // ghosts: periodic in x2 (x3), edge copy in x1; several MeshBlocks: the ordinary
  // cell-centred exchange, edge copy only at a physical x1 face
  if (vet_mbs != nullptr) {
    VetMBState &st = *vet_mbs;
    auto ds_ = st.dsc;
    auto lx_ = st.lx1;
    const int nbx1 = st.nbx1;
    par_for("vet_vet_full_dsc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      for (int n = 0; n < 6; ++n) {
        ds_(m,n,k,j,i) = vc_(m,VETRAY_D11+n,k,j,i);
      }
    });
    MeshBoundaryValuesCC *pb = st.pbd;
    while (pb->InitRecv(6) == TaskStatus::incomplete) {}
    while (pb->PackAndSendCC(st.dsc, st.dsc_c) == TaskStatus::incomplete) {}
    while (pb->RecvAndUnpackCC(st.dsc, st.dsc_c) == TaskStatus::incomplete) {}
    while (pb->ClearSend() == TaskStatus::incomplete) {}
    while (pb->ClearRecv() == TaskStatus::incomplete) {}
    par_for("vet_vet_full_ghmb", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      if (act) return;
      const int lb = lx_.d_view(m);
      int ii = i;
      if (i < is && lb == 0) ii = is;
      if (i > ie && lb == nbx1 - 1) ii = ie;
      for (int n = 0; n < 6; ++n) {
        vc_(m,VETRAY_D11+n,k,j,i) = ds_(m,n,k,j,ii);
      }
    });
  } else {
    par_for("vet_vet_full_gh", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      if (act) return;
      const int ii = (i < is) ? is : ((i > ie) ? ie : i);
      // edge copy at a vacuum (non-periodic) x2 (x3) face
      const int jj = x2o ? ((j < js) ? js : ((j > je) ? je : j))
                         : (js + VetWrap(j - js, nx2));
      const int kk = thrd ? (x3o ? ((k < ks) ? ks : ((k > ke) ? ke : k))
                                 : (ks + VetWrap(k - ks, nx3))) : k;
      for (int n = 0; n < 6; ++n) {
        vc_(m,VETRAY_D11+n,k,j,i) = vc_(m,VETRAY_D11+n,kk,jj,ii);
      }
    });
  }

  // run statistics of the guard (host reductions over the active cells)
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  const int nc = (nmb1 + 1)*nk*nj*ni;
  Real nhit = 0.0, gmax = 0.0, dmin = 1.0;
  Kokkos::parallel_reduce("vet_vet_full_st", Kokkos::RangePolicy<>(DevExeSpace(), 0, nc),
  KOKKOS_LAMBDA(const int idx, Real &sh, Real &sg, Real &sd) {
    const int m = idx/(nk*nj*ni);
    int r = idx - m*(nk*nj*ni);
    const int k = ks + r/(nj*ni);
    r -= (k - ks)*(nj*ni);
    const int j = js + r/ni;
    const int i = is + r - (j - js)*ni;
    const Real g = vc_(m,VETRAY_GD,k,j,i);
    if (g >= 0.0) {sh += 1.0;}
    sg = fmax(sg, (g >= 0.0) ? g : (-1.0 - g));
    sd = fmin(sd, vc_(m,VETRAY_D11,k,j,i));
  }, Kokkos::Sum<Real>(nhit), Kokkos::Max<Real>(gmax), Kokkos::Min<Real>(dmin));
  vet_nguard += nhit;
  vet_ncell += static_cast<Real>(nc);
  vet_guard_max = std::max(vet_guard_max, gmax);
  vet_d11_min = std::min(vet_d11_min, dmin);
}

}  // namespace radvet
