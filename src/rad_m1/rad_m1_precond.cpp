//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_precond.cpp
//! \brief implicit M1: preconditioner study (tests_m1/runs_5m_precond).
//!
//! <rad_m1>/implicit_dump_op = N (debug, read only when named): at the first implicit
//! solve of cycle N the stored stencil, b and x0 of the frozen system are written, per
//! rank, to m1op.c<N>.r<rank>.bin for the offline study (runs_5m_precond/study.py).
//!
//! <rad_m1>/implicit_precond = mg: the rbgs_fwd preconditioner (x1 line solves, one
//! forward red-black transverse sweep, block-local) followed by a coarse correction for
//! the smooth transverse error it leaves (runs_5m_precond/README.md):
//!   z0 = S r;  r1 = R (r - A7 z0);  z1 = V(r1);  z = z0 + P z1,
//! with A7 the 7-point part of the row (TA,TB,TC,CJM..CKP), P piecewise constant over
//! 2 x 2 (x2,x3) aggregates of the SAME x1 index inside one MeshBlock (x1 is never
//! coarsened: the line solves are exact there), R = P^T, and the Galerkin coarse rows
//! R A7 P, again 7-point, keeping only the couplings inside the MeshBlock.  V on a
//! coarse level is the same thing recursively (one forward red-black line sweep, then
//! the next level), down to implicit_mg_levels levels (the fine one included); the
//! coarsest level gets the sweep alone.  The only communication is the halo of z0 for
//! the fine residual (implicit_mg_halo = true, default; false = block-local residual).
//! M is a fixed linear map (a valid right preconditioner); it changes the iterates, not
//! the converged answer.
//!
//! <rad_m1>/implicit_precond = mg_gc (tests_m1/runs_5p_coarse2): a GLOBAL coarse space
//! in front of mg.  P_g is piecewise constant over gb2 x gb3 bands of the whole (x2,x3)
//! mesh (implicit_gc_bands2 / 3, default 1 x 1 = the sideways mean of each x1 layer),
//! one unknown per band and GLOBAL x1 index.  Multiplicative, coarse first:
//!   x_g = (P_g^T A P_g)^{-1} P_g^T r;  r' = r - A P_g x_g;  z = P_g x_g + M_mg r',
//! M_mg = mg with implicit_mg_levels (1 = rbgs_fwd alone).  I - A M is then
//! (I - A M_mg)(I - A P_g A_g^{-1} P_g^T), the same spectrum as the coarse-last order
//! (smoother, coarse on its residual), but P_g^T r needs no operator product and
//! A P_g x_g needs no halo (x_g is global).  A is the Krylov operator itself (the
//! stencil of ImplicitStencilOp, incl. edges and folded vimp): a coarse operator from
//! the 7-point preconditioner rows alone misses the vimp terms and made the iteration
//! WORSE with velocities (runs_5p_coarse2/README.md).  Per application: one fine kernel
//! (the p / s update + the partial band sums), one small sum kernel, a host round trip
//! with ONE MPI_Allreduce of N1 x nbands numbers and the banded LU solve on every rank,
//! one fine kernel for r' (5 row sums per cell for one band, the stencil otherwise),
//! and z += P_g x_g in the prolongation kernel.  The coarse matrix (banded, half
//! bandwidth 3 nbands - 1) is assembled with one MPI_Allreduce and factorised once per
//! Picard pass.  Couplings across a non-periodic x2/x3 mesh face and across the x1 mesh
//! faces are left out.  Needs a single-level mesh and implicit_op_stencil.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "bvals/bvals.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitDumpOp
//! \brief file layout (all int32 then float64): nmb, nst, n1, n2, n3; per block gid,
//! lx1, lx2, lx3, level; then per block st(0..nst-1), b, x0, TA, TB, TC, CJM, CJP, CKM,
//! CKP over the active cells (k,j,i order, i fastest)

void RadiationM1::ImplicitDumpOp() {
  auto *pm = pmy_pack->pmesh;
  if (impl_dump_cyc < 0 || pm->ncycle != impl_dump_cyc || impl_dump_done) {return;}
  impl_dump_done = true;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1 = indcs.nx1, n2 = indcs.nx2, n3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nst = static_cast<int>(ost.extent(1));
  auto hst = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ost);
  auto hiw = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iw);
  char fn[128];
  std::snprintf(fn, sizeof(fn), "m1op.c%d.r%d.bin", impl_dump_cyc,
                global_variable::my_rank);
  std::FILE *f = std::fopen(fn, "wb");
  if (f == nullptr) {return;}
  int hd[5] = {nmb, nst, n1, n2, n3};
  std::fwrite(hd, sizeof(int), 5, f);
  pmy_pack->pmb->mb_gid.sync_host();
  for (int m = 0; m < nmb; ++m) {
    int g = pmy_pack->pmb->mb_gid.h_view(m);
    auto &ll = pm->lloc_eachmb[g];
    int b[5] = {g, ll.lx1, ll.lx2, ll.lx3, ll.level};
    std::fwrite(b, sizeof(int), 5, f);
  }
  std::vector<double> buf(static_cast<size_t>(n1)*n2*n3);
  for (int m = 0; m < nmb; ++m) {
    const int cmp[9] = {M1_IW_KB, M1_IW_EP, M1_IW_TA, M1_IW_TB, M1_IW_TC, M1_IW_CJM,
                        M1_IW_CJP, M1_IW_CKM, M1_IW_CKP};
    for (int q = 0; q < nst + 9; ++q) {
      size_t p = 0;
      for (int k = ks; k < ks + n3; ++k) {
        for (int j = js; j < js + n2; ++j) {
          for (int i = is; i < is + n1; ++i) {
            if (q < nst) {
              buf[p++] = hst(m,q,k,j,i);
            } else {
              buf[p++] = hiw(m, cmp[q - nst], k, j, i);
            }
          }
        }
      }
      std::fwrite(buf.data(), sizeof(double), buf.size(), f);
    }
  }
  std::fclose(f);
}

namespace {
// component layout of a coarse level
constexpr int MG_TA = 0, MG_TB = 1, MG_TC = 2, MG_CJM = 3, MG_CJP = 4, MG_CKM = 5,
              MG_CKP = 6, MG_R = 7, MG_Z = 8, MG_NC = 9;

//! one colour of the forward red-black line sweep on a coarse level: columns (k,j) with
//! parity (k+j) = col solve A_line z = r (- the 5-point coupling to z of the other
//! colour when sub), by parallel cyclic reduction in team scratch (as M1PCRX)
template <typename T>
void M1PCRMG(const DvceArray5D<Real> &a_, const int nmb, const int nx, const int nj,
             const int nk, const int col, const bool sub, const bool thrd, int ts) {
  const int njl = (nj + 1)/2;
  const int nkj = nk*njl;
  size_t scr_size = ScrArray1D<T>::shmem_size(8*nx);
  int nround = 0;
  while ((1 << nround) < nx) ++nround;
  Kokkos::TeamPolicy<DevExeSpace> policy;
  if (!std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) {
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, ts);
  } else {
    policy = Kokkos::TeamPolicy<DevExeSpace>(DevExeSpace(), nmb*nkj, Kokkos::AUTO);
  }
  Kokkos::parallel_for("m1_mg_pcr", policy.set_scratch_size(0, Kokkos::PerTeam(scr_size)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    const int m = tm.league_rank()/nkj;
    const int k = (tm.league_rank() - m*nkj)/njl;
    const int jj = (tm.league_rank() - m*nkj)%njl;
    const int j = 2*jj + ((col + k) & 1);
    if (j >= nj) return;   // team-uniform
    ScrArray1D<T> sw(tm.team_scratch(0), 8*nx);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
      sw(i) = static_cast<T>((i == 0) ? 0.0 : a_(m,MG_TA,k,j,i));
      sw(nx + i) = static_cast<T>(a_(m,MG_TB,k,j,i));
      sw(2*nx + i) = static_cast<T>((i == nx-1) ? 0.0 : a_(m,MG_TC,k,j,i));
      Real rr = a_(m,MG_R,k,j,i);
      if (sub) {
        if (j > 0) rr -= a_(m,MG_CJM,k,j,i)*a_(m,MG_Z,k,j-1,i);
        if (j < nj-1) rr -= a_(m,MG_CJP,k,j,i)*a_(m,MG_Z,k,j+1,i);
        if (thrd) {
          if (k > 0) rr -= a_(m,MG_CKM,k,j,i)*a_(m,MG_Z,k-1,j,i);
          if (k < nk-1) rr -= a_(m,MG_CKP,k,j,i)*a_(m,MG_Z,k+1,j,i);
        }
      }
      sw(3*nx + i) = static_cast<T>(rr);
    });
    tm.team_barrier();
    int src = 0;
    for (int rd=0, s=1; rd<nround; ++rd, s*=2) {
      const int o = src*4*nx, d = (1-src)*4*nx;
      Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
        const int im = i - s, ip = i + s;
        T ai = sw(o + i), bi = sw(o + nx + i), ci = sw(o + 2*nx + i);
        T ri = sw(o + 3*nx + i);
        T an = 0.0, cn = 0.0;
        if (im >= 0) {
          T f = -ai/sw(o + nx + im);
          an = f*sw(o + im);
          bi += f*sw(o + 2*nx + im);
          ri += f*sw(o + 3*nx + im);
        }
        if (ip < nx) {
          T g = -ci/sw(o + nx + ip);
          cn = g*sw(o + 2*nx + ip);
          bi += g*sw(o + ip);
          ri += g*sw(o + 3*nx + ip);
        }
        sw(d + i) = an;
        sw(d + nx + i) = bi;
        sw(d + 2*nx + i) = cn;
        sw(d + 3*nx + i) = ri;
      });
      tm.team_barrier();
      src = 1 - src;
    }
    const int o = src*4*nx;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(tm, nx), [&](const int i) {
      a_(m,MG_Z,k,j,i) = static_cast<Real>(sw(o + 3*nx + i)/sw(o + nx + i));
    });
  });
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitMGBuild
//! \brief implicit_precond = mg: the level sizes (once) and the Galerkin coarse rows of
//! this pass, level by level from the fine rows TA..CKP of the work array

void RadiationM1::ImplicitMGBuild() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  if (mgc.empty() || (mgc.size() > 1 && mgc[1].extent_int(0) != nmb)) {
    mgc.clear(); mg_nj.clear(); mg_nk.clear();
    mg_nj.push_back(indcs.nx2);
    mg_nk.push_back(thrd ? indcs.nx3 : 1);
    mgc.push_back(DvceArray5D<Real>());   // level 0 is the work array itself
    for (int l = 1; l < mg_nlev; ++l) {
      const int pj = mg_nj[l-1], pk = mg_nk[l-1];
      if (pj == 1 && pk == 1) break;
      mg_nj.push_back((pj + 1)/2);
      mg_nk.push_back((pk + 1)/2);
      mgc.push_back(DvceArray5D<Real>("m1_mgc", nmb, MG_NC, mg_nk[l], mg_nj[l], nx));
    }
    if (gc_on) {ImplicitGCInit();}
    if (gf_on) {ImplicitGFInit();}
  }
  const int nl = static_cast<int>(mgc.size());
  auto iw_ = iw;
  for (int l = 1; l < nl; ++l) {
    const bool fine = (l == 1);
    auto f_ = mgc[l-1];
    auto c_ = mgc[l];
    const int fj = mg_nj[l-1], fk = mg_nk[l-1];
    const int cj = mg_nj[l], ck = mg_nk[l];
    par_for("m1_mg_build", DevExeSpace(), 0, nmb-1, 0, ck-1, 0, cj-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int i) {
      Real ta = 0.0, tb = 0.0, tc = 0.0, jm = 0.0, jp = 0.0, km = 0.0, kp = 0.0;
      for (int a = 0; a < 2; ++a) {
        const int k = 2*kc + a;
        if (k >= fk) continue;
        for (int b = 0; b < 2; ++b) {
          const int j = 2*jc + b;
          if (j >= fj) continue;
          Real v[7];
          if (fine) {
            v[0] = iw_(m,M1_IW_TA,k+ks,j+js,i+is);
            v[1] = iw_(m,M1_IW_TB,k+ks,j+js,i+is);
            v[2] = iw_(m,M1_IW_TC,k+ks,j+js,i+is);
            v[3] = iw_(m,M1_IW_CJM,k+ks,j+js,i+is);
            v[4] = iw_(m,M1_IW_CJP,k+ks,j+js,i+is);
            v[5] = thrd ? iw_(m,M1_IW_CKM,k+ks,j+js,i+is) : 0.0;
            v[6] = thrd ? iw_(m,M1_IW_CKP,k+ks,j+js,i+is) : 0.0;
          } else {
            for (int q = 0; q < 7; ++q) {v[q] = f_(m,q,k,j,i);}
          }
          ta += v[0]; tb += v[1]; tc += v[2];
          // a coupling inside the aggregate goes onto the diagonal, one to another
          // aggregate of the block stays, one across the block face is dropped
          if (j - 1 >= 2*jc) {tb += v[3];} else if (j - 1 >= 0) {jm += v[3];}
          if (j + 1 <= 2*jc + 1 && j + 1 < fj) {
            tb += v[4];
          } else if (j + 1 < fj) {
            jp += v[4];
          }
          if (thrd) {
            if (k - 1 >= 2*kc) {tb += v[5];} else if (k - 1 >= 0) {km += v[5];}
            if (k + 1 <= 2*kc + 1 && k + 1 < fk) {
              tb += v[6];
            } else if (k + 1 < fk) {
              kp += v[6];
            }
          }
        }
      }
      c_(m,MG_TA,kc,jc,i) = ta; c_(m,MG_TB,kc,jc,i) = tb; c_(m,MG_TC,kc,jc,i) = tc;
      c_(m,MG_CJM,kc,jc,i) = jm; c_(m,MG_CJP,kc,jc,i) = jp;
      c_(m,MG_CKM,kc,jc,i) = km; c_(m,MG_CKP,kc,jc,i) = kp;
    });
  }
  if (gc_on) {ImplicitGCBuild();}
  if (gf_on) {ImplicitGFBuild();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitMGApply
//! \brief implicit_precond = mg: z = M^{-1} r (see the file header); rc / upd / c1 / c2
//! as ImplicitPCRSolveX

void RadiationM1::ImplicitMGApply(int rc, int zc, int upd, Real c1, Real c2) {
  if (gc_on) {   // mg_gc: the global coarse solve first, then mg on r' = r - A P_g x_g
    ImplicitGCPre(rc, upd, c1, c2);
    rc = M1_IW_S1;
    upd = 0;
  }
  const bool gfa = gf_on && gf_ok;
  if (gfa) {     // mg_gf: the Fourier coarse solve first, then mg on r' = r - Pi r
    ImplicitGFPre(rc, upd, c1, c2);
    rc = M1_IW_S1;
    upd = 0;
  }
  // the fine smoother: rbgs_fwd, exactly the implicit_precond = rbgs_fwd map
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 0, -1);
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 1, zc);
  if (mgc.size() < 2) {
    if (gc_on) {ImplicitGCAdd(zc);}
    if (gfa) {ImplicitGFAdd(zc);}
    return;
  }
  const int nl = static_cast<int>(mgc.size());
  const int rr = (upd == 1) ? M1_IW_KP : ((upd == 2) ? M1_IW_KS : rc);
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const bool hal = mg_halo;
  if (hal) {ImplicitKrylovHalo(zc);}
  auto iw_ = iw;
  int ts = impl_pcr_team;
  if (ts == 0) {
    ts = 1;
    while (ts < nx && ts < 256) ts *= 2;
  }
  // the fine residual r - A7 z0, restricted
  {
    auto c_ = mgc[1];
    const int fj = mg_nj[0], fk = mg_nk[0];
    const int cj = mg_nj[1], ck = mg_nk[1];
    par_for("m1_mg_res0", DevExeSpace(), 0, nmb-1, 0, ck-1, 0, cj-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int i) {
      Real sum = 0.0;
      const int ii = i + is;
      for (int a = 0; a < 2; ++a) {
        const int k = 2*kc + a;
        if (k >= fk) continue;
        const int kk = k + ks;
        for (int b = 0; b < 2; ++b) {
          const int j = 2*jc + b;
          if (j >= fj) continue;
          const int jj = j + js;
          Real y = iw_(m,M1_IW_TB,kk,jj,ii)*iw_(m,zc,kk,jj,ii);
          if (hal || i > 0) y += iw_(m,M1_IW_TA,kk,jj,ii)*iw_(m,zc,kk,jj,ii-1);
          if (hal || i < nx-1) y += iw_(m,M1_IW_TC,kk,jj,ii)*iw_(m,zc,kk,jj,ii+1);
          if (hal || j > 0) y += iw_(m,M1_IW_CJM,kk,jj,ii)*iw_(m,zc,kk,jj-1,ii);
          if (hal || j < fj-1) y += iw_(m,M1_IW_CJP,kk,jj,ii)*iw_(m,zc,kk,jj+1,ii);
          if (thrd) {
            if (hal || k > 0) y += iw_(m,M1_IW_CKM,kk,jj,ii)*iw_(m,zc,kk-1,jj,ii);
            if (hal || k < fk-1) y += iw_(m,M1_IW_CKP,kk,jj,ii)*iw_(m,zc,kk+1,jj,ii);
          }
          sum += iw_(m,rr,kk,jj,ii) - y;
        }
      }
      c_(m,MG_R,kc,jc,i) = sum;
    });
  }
  // down the coarse levels: one forward red-black sweep each, then the restricted
  // residual of the next level (block-local)
  for (int l = 1; l < nl; ++l) {
    auto c_ = mgc[l];
    const int cj = mg_nj[l], ck = mg_nk[l];
    const bool th = thrd && (ck > 1);
    if (impl_prec_float) {
      M1PCRMG<float>(c_, nmb, nx, cj, ck, 0, false, th, ts);
      M1PCRMG<float>(c_, nmb, nx, cj, ck, 1, true, th, ts);
    } else {
      M1PCRMG<Real>(c_, nmb, nx, cj, ck, 0, false, th, ts);
      M1PCRMG<Real>(c_, nmb, nx, cj, ck, 1, true, th, ts);
    }
    if (l + 1 < nl) {
      auto n_ = mgc[l+1];
      const int nj2 = mg_nj[l+1], nk2 = mg_nk[l+1];
      par_for("m1_mg_res", DevExeSpace(), 0, nmb-1, 0, nk2-1, 0, nj2-1, 0, nx-1,
      KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int i) {
        Real sum = 0.0;
        for (int a = 0; a < 2; ++a) {
          const int k = 2*kc + a;
          if (k >= ck) continue;
          for (int b = 0; b < 2; ++b) {
            const int j = 2*jc + b;
            if (j >= cj) continue;
            Real y = c_(m,MG_TB,k,j,i)*c_(m,MG_Z,k,j,i);
            if (i > 0) y += c_(m,MG_TA,k,j,i)*c_(m,MG_Z,k,j,i-1);
            if (i < nx-1) y += c_(m,MG_TC,k,j,i)*c_(m,MG_Z,k,j,i+1);
            if (j > 0) y += c_(m,MG_CJM,k,j,i)*c_(m,MG_Z,k,j-1,i);
            if (j < cj-1) y += c_(m,MG_CJP,k,j,i)*c_(m,MG_Z,k,j+1,i);
            if (th) {
              if (k > 0) y += c_(m,MG_CKM,k,j,i)*c_(m,MG_Z,k-1,j,i);
              if (k < ck-1) y += c_(m,MG_CKP,k,j,i)*c_(m,MG_Z,k+1,j,i);
            }
            sum += c_(m,MG_R,k,j,i) - y;
          }
        }
        n_(m,MG_R,kc,jc,i) = sum;
      });
    }
  }
  // back up, in ONE kernel: the aggregates of level l are 2^l x 2^l fine columns, so
  // z += P_1 (z_1 + P_2 (z_2 + ...)) is z(k,j,i) += sum_l z_l(k >> l, j >> l, i)
  auto c1_ = mgc[1];
  auto c2_ = (nl > 2) ? mgc[2] : mgc[1];
  auto c3_ = (nl > 3) ? mgc[3] : mgc[1];
  auto c4_ = (nl > 4) ? mgc[4] : mgc[1];
  if (nl > 5) {
    for (int l = nl - 2; l >= 4; --l) {   // deeper than 4 coarse levels: fold them first
      auto c_ = mgc[l];
      auto n_ = mgc[l+1];
      par_for("m1_mg_pro", DevExeSpace(), 0, nmb-1, 0, mg_nk[l]-1, 0, mg_nj[l]-1, 0, nx-1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        c_(m,MG_Z,k,j,i) += n_(m,MG_Z,k/2,j/2,i);
      });
    }
  }
  const int nlc = std::min(nl - 1, 4);   // coarse levels read by the fine kernel
  if (gfa) {     // the same, plus P x_g of the Fourier coarse space
    auto hx_ = gf_hx;
    auto go_ = gc_off;
    auto v_ = gf_v;
    const int n1 = gc_n1, n2 = gc_n2, n3 = gc_n3, nb3 = gf_nb3;
    par_for("m1_gf_pro0", DevExeSpace(), 0, nmb-1, 0, mg_nk[0]-1, 0, mg_nj[0]-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real s = c1_(m,MG_Z,k >> 1,j >> 1,i);
      if (nlc > 1) s += c2_(m,MG_Z,k >> 2,j >> 2,i);
      if (nlc > 2) s += c3_(m,MG_Z,k >> 3,j >> 3,i);
      if (nlc > 3) s += c4_(m,MG_Z,k >> 4,j >> 4,i);
      const int gi = go_(3*m) + i, gj = go_(3*m+1) + j, gk = go_(3*m+2) + k;
      Real g = 0.0;
      for (int c = 0; c < nb3; ++c) {g += v_(c*n3 + gk)*hx_((c*n2 + gj)*n1 + gi);}
      iw_(m,zc,k+ks,j+js,i+is) += s + g;
    });
    return;
  }
  if (gc_on) {   // the same, plus P_g x_g
    auto gx_ = gc_x;
    auto go_ = gc_off;
    const int n2 = gc_n2, n3 = gc_n3, gb2 = gc_b2, gb3 = gc_b3, nb = gc_nb;
    const bool gok = gc_ok;
    par_for("m1_gc_pro0", DevExeSpace(), 0, nmb-1, 0, mg_nk[0]-1, 0, mg_nj[0]-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real s = 0.0;
      if (nlc > 0) s += c1_(m,MG_Z,k >> 1,j >> 1,i);
      if (nlc > 1) s += c2_(m,MG_Z,k >> 2,j >> 2,i);
      if (nlc > 2) s += c3_(m,MG_Z,k >> 3,j >> 3,i);
      if (nlc > 3) s += c4_(m,MG_Z,k >> 4,j >> 4,i);
      if (gok) {
        const int b2 = ((go_(3*m+1) + j)*gb2)/n2;
        const int b3 = ((go_(3*m+2) + k)*gb3)/n3;
        s += gx_((go_(3*m) + i)*nb + b3*gb2 + b2);
      }
      iw_(m,zc,k+ks,j+js,i+is) += s;
    });
    return;
  }
  par_for("m1_mg_pro0", DevExeSpace(), 0, nmb-1, 0, mg_nk[0]-1, 0, mg_nj[0]-1, 0, nx-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real s = c1_(m,MG_Z,k >> 1,j >> 1,i);
    if (nlc > 1) s += c2_(m,MG_Z,k >> 2,j >> 2,i);
    if (nlc > 2) s += c3_(m,MG_Z,k >> 3,j >> 3,i);
    if (nlc > 3) s += c4_(m,MG_Z,k >> 4,j >> 4,i);
    iw_(m,zc,k+ks,j+js,i+is) += s;
  });
}

namespace {
[[noreturn]] void GCFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

//! the neighbour offset (di, dj, dk) of stencil slot s of ImplicitStencilOp, and whether
//! the slot is part of the operator (edges, 3-D, the folded vimp i/j/k +- 2 slots)
KOKKOS_INLINE_FUNCTION
bool GCSlot(const int s, const bool thrd, const bool edg, const bool vfold, int &di,
            int &dj, int &dk) {
  const int tdi[25] = {0, -1, 1, 0, 0, 0, 0, -1, 1, -1, 1, -1, 1, -1, 1,
                       0, 0, 0, 0, -2, 2, 0, 0, 0, 0};
  const int tdj[25] = {0, 0, 0, -1, 1, 0, 0, -1, -1, 1, 1, 0, 0, 0, 0,
                       -1, 1, -1, 1, 0, 0, -2, 2, 0, 0};
  const int tdk[25] = {0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, -1, -1, 1, 1,
                       -1, -1, 1, 1, 0, 0, 0, 0, -2, 2};
  di = tdi[s]; dj = tdj[s]; dk = tdk[s];
  if (s < 5) return true;
  if (s < 7) return thrd;
  if (s < 11) return edg;
  if (s < 19) return thrd && edg;
  if (s < 23) return vfold;
  return vfold && thrd;
}

//! the band class of a neighbour at global (gj + dj, gk + dk) seen from band (b2, b3):
//! -1 = outside a non-periodic mesh face (the coupling is dropped), else
//! (e3 + 1)*3 + (e2 + 1) with e = -1, 0, +1 the band step (bands are >= 2 cells wide)
KOKKOS_INLINE_FUNCTION
int GCClass(int gj, int gk, const int dj, const int dk, const int n2, const int n3,
            const int gb2, const int gb3, const bool per2, const bool per3) {
  const int b2 = (gj*gb2)/n2, b3 = (gk*gb3)/n3;
  gj += dj; gk += dk;
  if (gj < 0 || gj >= n2) {
    if (!per2) return -1;
    gj = (gj + n2)%n2;
  }
  if (gk < 0 || gk >= n3) {
    if (!per3) return -1;
    gk = (gk + n3)%n3;
  }
  const int e2 = ((gj*gb2)/n2 == b2) ? 0 : ((dj < 0) ? -1 : 1);
  const int e3 = ((gk*gb3)/n3 == b3) ? 0 : ((dk < 0) ? -1 : 1);
  return (e3 + 1)*3 + (e2 + 1);
}

//! gc_part((((m*nlb + lb)*nq + q)*nrc + c)*nx + i) = the sum over the rows (k,j) of
//! band slot lb of block m in row chunk c (rows c, c + nrc, ...) of f(m, q, k, j, i)
//! (active indices from 0); the host adds the nrc chunks.  One team per (m, lb, q, c,
//! 16 x1 cells): 16 row groups x 16 x1 cells, each virtual thread sums its rows, then
//! the 16 row groups are added in a fixed order (deterministic, coalesced in x1; any
//! team size, as the 256 virtual threads are strided over the real ones)
template <class F>
void GCBandReduce(const int nmb, const int nlb, const int nq, const int nrc, const int nx,
                  const DvceArray1D<int> &bi_,
                  const Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace> &p_, F f) {
  constexpr int IB = 16, RG = 16, NV = IB*RG;
  const int nib = (nx + IB - 1)/IB;
  const int nlg = nmb*nlb*nq*nrc*nib;
  int ts = NV;
  if (std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value) ts = 1;
  size_t scr = ScrArray1D<Real>::shmem_size(NV);
  Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), nlg, ts);
  Kokkos::parallel_for("m1_gc_band", pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
  KOKKOS_LAMBDA(TeamMember_t tm) {
    int l = tm.league_rank();
    const int ib = l%nib; l /= nib;
    const int c = l%nrc; l /= nrc;
    const int q = l%nq; l /= nq;
    const int lb = l%nlb;
    const int m = l/nlb;
    const int e = 5*(m*nlb + lb);
    ScrArray1D<Real> sw(tm.team_scratch(0), NV);
    const bool in = (bi_(e) >= 0);
    const int j0 = bi_(e+1), nj = bi_(e+2) - j0, k0 = bi_(e+3);
    const int nrow = in ? nj*(bi_(e+4) - k0) : 0;
    for (int t = tm.team_rank(); t < NV; t += tm.team_size()) {
      const int i = ib*IB + t%IB;
      Real a = 0.0;
      if (i < nx) {
        for (int rr = c + nrc*(t/IB); rr < nrow; rr += nrc*RG) {
          a += f(m, q, k0 + rr/nj, j0 + rr%nj, i);
        }
      }
      sw(t) = a;
    }
    tm.team_barrier();
    for (int t = tm.team_rank(); t < IB; t += tm.team_size()) {
      const int i = ib*IB + t;
      if (i < nx) {
        Real a = 0.0;
        for (int g = 0; g < RG; ++g) {a += sw(g*IB + t);}
        p_((((m*nlb + lb)*nq + q)*nrc + c)*nx + i) = a;
      }
    }
  });
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGCInit
//! \brief implicit_precond = mg_gc: checks, the blocks' global offsets and band slots,
//! the storage (once)

void RadiationM1::ImplicitGCInit() {
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const bool thrd = trans_x3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = indcs.nx1, nx2 = indcs.nx2, nx3 = thrd ? indcs.nx3 : 1;
  if (pm->multilevel) {
    GCFatal("<rad_m1>/implicit_precond = mg_gc needs a single-level mesh");
  }
  if (!impl_stencil) {
    GCFatal("<rad_m1>/implicit_precond = mg_gc needs implicit_op_stencil");
  }
  gc_n1 = pm->mesh_indcs.nx1;
  gc_n2 = pm->mesh_indcs.nx2;
  gc_n3 = thrd ? pm->mesh_indcs.nx3 : 1;
  if (!thrd) {gc_b3 = 1;}
  if (gc_b2 < 1 || gc_b3 < 1 || gc_n2 % gc_b2 != 0 || (gc_b2 > 1 && gc_n2/gc_b2 < 2) ||
      gc_n3 % gc_b3 != 0 || (gc_b3 > 1 && gc_n3/gc_b3 < 2)) {
    GCFatal("<rad_m1>/implicit_gc_bands2/3 must divide the mesh nx2/nx3 into bands of "
            ">= 2 cells");
  }
  gc_nb = gc_b2*gc_b3;
  gc_kl = 3*gc_nb - 1;
  gc_per2 = (pm->mesh_bcs[static_cast<int>(BoundaryFace::inner_x2)]
             == BoundaryFlag::periodic);
  gc_per3 = (pm->mesh_bcs[static_cast<int>(BoundaryFace::inner_x3)]
             == BoundaryFlag::periodic);
  const int w2 = gc_n2/gc_b2, w3 = gc_n3/gc_b3;
  pmy_pack->pmb->mb_gid.sync_host();
  gc_off_h.assign(3*nmb, 0);
  int nlb3 = 1;
  gc_nlb2 = 1;
  for (int m = 0; m < nmb; ++m) {
    auto &ll = pm->lloc_eachmb[pmy_pack->pmb->mb_gid.h_view(m)];
    gc_off_h[3*m] = static_cast<int>(ll.lx1)*nx;
    gc_off_h[3*m+1] = static_cast<int>(ll.lx2)*nx2;
    gc_off_h[3*m+2] = thrd ? static_cast<int>(ll.lx3)*nx3 : 0;
    const int o2 = gc_off_h[3*m+1], o3 = gc_off_h[3*m+2];
    gc_nlb2 = std::max(gc_nlb2, (o2 + nx2 - 1)/w2 - o2/w2 + 1);
    nlb3 = std::max(nlb3, (o3 + nx3 - 1)/w3 - o3/w3 + 1);
  }
  gc_nlb = gc_nlb2*nlb3;
  gc_binfo_h.assign(5*nmb*gc_nlb, 0);
  for (int m = 0; m < nmb; ++m) {
    const int o2 = gc_off_h[3*m+1], o3 = gc_off_h[3*m+2];
    for (int lb = 0; lb < gc_nlb; ++lb) {
      int *e = &gc_binfo_h[5*(m*gc_nlb + lb)];
      const int b2 = o2/w2 + lb%gc_nlb2, b3 = o3/w3 + lb/gc_nlb2;
      if (b2 > (o2 + nx2 - 1)/w2 || b3 > (o3 + nx3 - 1)/w3) {e[0] = -1; continue;}
      e[0] = b3*gc_b2 + b2;
      e[1] = std::max(o2, b2*w2) - o2;
      e[2] = std::min(o2 + nx2, (b2 + 1)*w2) - o2;
      e[3] = std::max(o3, b3*w3) - o3;
      e[4] = std::min(o3 + nx3, (b3 + 1)*w3) - o3;
    }
  }
  gc_off = DvceArray1D<int>("m1_gc_off", 3*nmb);
  gc_binfo = DvceArray1D<int>("m1_gc_binfo", 5*nmb*gc_nlb);
  {
    auto h1 = Kokkos::create_mirror_view(gc_off);
    for (int q = 0; q < 3*nmb; ++q) {h1(q) = gc_off_h[q];}
    Kokkos::deep_copy(gc_off, h1);
    auto h2 = Kokkos::create_mirror_view(gc_binfo);
    for (int q = 0; q < 5*nmb*gc_nlb; ++q) {h2(q) = gc_binfo_h[q];}
    Kokkos::deep_copy(gc_binfo, h2);
  }
  if (gc_nb == 1) {gcw = DvceArray5D<Real>("m1_gcw", nmb, 5, nx3, nx2, nx);}
  if (gc_nb > 1) {
    gc_pd = DvceArray1D<Real>("m1_gc_pd", static_cast<size_t>(nmb)*nx3*gc_nlb2*45*nx);
  }
  // row chunks of the band reductions: ~8 rows per thread of GCBandReduce
  int nrow = 1;
  for (int q = 0; q < nmb*gc_nlb; ++q) {
    const int *e = &gc_binfo_h[5*q];
    if (e[0] >= 0) {nrow = std::max(nrow, (e[2] - e[1])*(e[4] - e[3]));}
  }
  gc_nrc = std::max(1, (nrow + 127)/128);
  gc_part = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_gc_part",
                static_cast<size_t>(nmb)*gc_nlb*45*nx*std::max(gc_nrc, 1));
  const int n = gc_n1*gc_nb;
  gc_x = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_gc_x", n);
  gc_lu.assign(static_cast<size_t>(n)*(3*gc_kl + 1), 0.0);
  gc_piv.assign(n, 0);
  gc_g.assign(n, 0.0);
  gc_ok = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGCSum
//! \brief the x2 partial sums gc_pd(m,k,lb2,q,i) summed over x3 in each band slot of
//! each block, into gc_part(m,lb,q,i) (host-pinned; one thread per output)

void RadiationM1::ImplicitGCSum(int nq) {
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = pmy_pack->pmesh->mb_indcs.nx1;
  const int nk = trans_x3 ? pmy_pack->pmesh->mb_indcs.nx3 : 1;
  const int nlb = gc_nlb, nlb2 = gc_nlb2;
  auto bi_ = gc_binfo;
  auto pd_ = gc_pd;
  auto p_ = gc_part;
  par_for("m1_gc_sum", DevExeSpace(), 0, nmb-1, 0, nlb-1, 0, nq-1, 0, nx-1,
  KOKKOS_LAMBDA(const int m, const int lb, const int q, const int i) {
    const int e = 5*(m*nlb + lb);
    const int l2 = lb%nlb2;
    Real s = 0.0;
    if (bi_(e) >= 0) {
      for (int k = bi_(e+3); k < bi_(e+4); ++k) {
        s += pd_((((m*nk + k)*nlb2 + l2)*nq + q)*nx + i);
      }
    }
    p_(((m*nlb + lb)*nq + q)*nx + i) = s;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGCBuild
//! \brief the coarse operator P_g^T A P_g of this pass from the stencil (the Krylov
//! operator itself), per x1 offset and band step; its banded LU (host, every rank)

void RadiationM1::ImplicitGCBuild() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const int nk = thrd ? indcs.nx3 : 1;
  const int n1 = gc_n1, n2 = gc_n2, n3 = gc_n3, gb2 = gc_b2, gb3 = gc_b3;
  const bool per2 = gc_per2, per3 = gc_per3, one = (gc_nb == 1);
  const bool edg = st_edges, vfold = impl_vfold && vimp_now;
  const int nst = ost.extent_int(1);
  const int nlb = gc_nlb, nlb2 = gc_nlb2;
  auto st_ = ost;
  auto w_ = gcw;
  auto pd_ = gc_pd;
  auto bi_ = gc_binfo;
  auto go_ = gc_off;
  const int nb = gc_nb, n = n1*nb, kl = gc_kl, w = 3*kl + 1;
  if (one) {
    // one band: the row sums of each cell per x1 offset (gcw, read by ImplicitGCPre),
    // then their sums over the (x2,x3) extent of each block
    par_for("m1_gc_w", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, indcs.nx2-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int gi = go_(3*m) + i, gj = go_(3*m+1) + j, gk = go_(3*m+2) + k;
      Real wd[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (int s = 0; s < nst; ++s) {
        int di, dj, dk;
        if (!GCSlot(s, thrd, edg, vfold, di, dj, dk)) continue;
        if (gi + di < 0 || gi + di >= n1) continue;
        if (GCClass(gj, gk, dj, dk, n2, n3, 1, 1, per2, per3) < 0) continue;
        wd[di + 2] += st_(m,s,k+ks,j+js,i+is);
      }
      for (int d = 0; d < 5; ++d) {w_(m,d,k,j,i) = wd[d];}
    });
    GCBandReduce(nmb, nlb, 5, gc_nrc, nx, bi_, gc_part,
    KOKKOS_LAMBDA(const int m, const int q, const int k, const int j, const int i) {
      return w_(m,q,k,j,i);
    });
  } else {
    par_for("m1_gc_rows", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, nlb2-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int k, const int l2, const int i) {
      Real acc[45];
      for (int q = 0; q < 45; ++q) {acc[q] = 0.0;}
      const int e = 5*(m*nlb + l2);
      const int gi = go_(3*m) + i, gk = go_(3*m+2) + k;
      if (bi_(e) >= 0) {
        for (int j = bi_(e+1); j < bi_(e+2); ++j) {
          const int gj = go_(3*m+1) + j;
          for (int s = 0; s < nst; ++s) {
            int di, dj, dk;
            if (!GCSlot(s, thrd, edg, vfold, di, dj, dk)) continue;
            if (gi + di < 0 || gi + di >= n1) continue;
            const int c = GCClass(gj, gk, dj, dk, n2, n3, gb2, gb3, per2, per3);
            if (c < 0) continue;
            acc[(di + 2)*9 + c] += st_(m,s,k+ks,j+js,i+is);
          }
        }
      }
      for (int q = 0; q < 45; ++q) {pd_((((m*nk + k)*nlb2 + l2)*45 + q)*nx + i) = acc[q];}
    });
  ImplicitGCSum(45);
  }
  Kokkos::fence();
  // assemble the banded matrix (fixed order), then over the ranks
  const int nqa = one ? 5 : 45, nrc = gc_nrc;
  std::fill(gc_lu.begin(), gc_lu.end(), 0.0);
  for (int m = 0; m < nmb; ++m) {
    for (int lb = 0; lb < nlb; ++lb) {
      const int b = gc_binfo_h[5*(m*nlb + lb)];
      if (b < 0) continue;
      const int b2 = b%gb2, b3 = b/gb2;
      for (int q = 0; q < nqa; ++q) {
        const int qq = one ? q*9 + 4 : q;   // one band: q = x1 offset q - 2, same band
        const int di = qq/9 - 2, e3 = (qq/3)%3 - 1, e2 = qq%3 - 1;
        const int cb = ((b3 + e3 + gb3)%gb3)*gb2 + (b2 + e2 + gb2)%gb2;
        for (int i = 0; i < nx; ++i) {
          Real a = 0.0;
          if (one) {
            for (int c = 0; c < nrc; ++c) {
              a += gc_part((((m*nlb + lb)*nqa + q)*nrc + c)*nx + i);
            }
          } else {
            a = gc_part(((m*nlb + lb)*nqa + q)*nx + i);
          }
          if (a == 0.0) continue;
          const int gi = gc_off_h[3*m] + i;
          const int r = gi*nb + b, c = (gi + di)*nb + cb;
          gc_lu[static_cast<size_t>(r)*w + (c - r + kl)] += a;
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, gc_lu.data(), static_cast<int>(gc_lu.size()),
                  MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  // banded LU with partial pivoting: row r holds columns r-kl .. r+2kl (the fill);
  // the multipliers stay where they are made (the solve applies the swaps in order)
  auto el = [&](int r, int c) -> Real & {
    return gc_lu[static_cast<size_t>(r)*w + (c - r + kl)];
  };
  gc_ok = true;
  for (int c = 0; c < n && gc_ok; ++c) {
    const int rmax = std::min(n - 1, c + kl), cmax = std::min(n - 1, c + 2*kl);
    int p = c;
    for (int r = c + 1; r <= rmax; ++r) {
      if (std::fabs(el(r, c)) > std::fabs(el(p, c))) p = r;
    }
    gc_piv[c] = p;
    if (!(std::fabs(el(p, c)) > 0.0) || !std::isfinite(el(p, c))) {gc_ok = false; break;}
    if (p != c) {
      for (int q = c; q <= cmax; ++q) {std::swap(el(p, q), el(c, q));}
    }
    const Real ip = 1.0/el(c, c);
    for (int r = c + 1; r <= rmax; ++r) {
      const Real f = el(r, c)*ip;
      el(r, c) = f;
      if (f == 0.0) continue;
      for (int q = c + 1; q <= cmax; ++q) {el(r, q) -= f*el(c, q);}
    }
  }
  if (!gc_ok && global_variable::my_rank == 0) {
    std::cout << "### WARNING rad_m1 mg_gc: singular global coarse operator, the coarse "
              << "correction is skipped this pass" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGCPre
//! \brief mg_gc, before mg: the right-hand side (made as the BiCGStab p / s update when
//! upd > 0, as ImplicitPCRSolveX), its band sums, one MPI_Allreduce, x_g by the banded
//! LU on the host (gc_x, host-pinned), and r' = r - A P_g x_g into M1_IW_S1 (free during
//! the Krylov solve: the Thomas scratch, not used by the pcr path)

void RadiationM1::ImplicitGCPre(int rc, int upd, Real c1, Real c2) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1, nj = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const int nk = thrd ? indcs.nx3 : 1;
  const int nlb = gc_nlb, nrc = gc_nrc;
  const int nb = gc_nb, n1 = gc_n1, n = n1*nb;
  const int up = upd, cr = rc;
  const Real a1 = c1, a2 = c2;
  auto iw_ = iw;
  auto bi_ = gc_binfo;
  GCBandReduce(nmb, nlb, 1, gc_nrc, nx, bi_, gc_part,
  KOKKOS_LAMBDA(const int m, const int, const int k, const int j, const int i) {
    const int kk = k + ks, jj = j + js, ii = i + is;
    Real v;
    if (up == 1) {
      v = iw_(m,M1_IW_KR,kk,jj,ii) + a1*(iw_(m,M1_IW_KP,kk,jj,ii)
                                         - a2*iw_(m,M1_IW_KV,kk,jj,ii));
      iw_(m,M1_IW_KP,kk,jj,ii) = v;
    } else if (up == 2) {
      v = iw_(m,M1_IW_KR,kk,jj,ii) - a1*iw_(m,M1_IW_KV,kk,jj,ii);
      iw_(m,M1_IW_KS,kk,jj,ii) = v;
    } else {
      v = iw_(m,cr,kk,jj,ii);
    }
    return v;
  });
  Kokkos::fence();
  std::fill(gc_g.begin(), gc_g.end(), 0.0);
  for (int m = 0; m < nmb; ++m) {
    for (int lb = 0; lb < nlb; ++lb) {
      const int b = gc_binfo_h[5*(m*nlb + lb)];
      if (b < 0) continue;
      for (int i = 0; i < nx; ++i) {
        Real a = 0.0;
        for (int c = 0; c < nrc; ++c) {a += gc_part(((m*nlb + lb)*nrc + c)*nx + i);}
        gc_g[(gc_off_h[3*m] + i)*nb + b] += a;
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    MPI_Allreduce(MPI_IN_PLACE, gc_g.data(), n, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  if (gc_ok) {
    const int kl = gc_kl, w = 3*kl + 1;
    auto el = [&](int r, int c) {return gc_lu[static_cast<size_t>(r)*w + (c - r + kl)];};
    for (int c = 0; c < n; ++c) {
      if (gc_piv[c] != c) {std::swap(gc_g[c], gc_g[gc_piv[c]]);}
      const int rmax = std::min(n - 1, c + kl);
      for (int r = c + 1; r <= rmax; ++r) {gc_g[r] -= el(r, c)*gc_g[c];}
    }
    for (int r = n - 1; r >= 0; --r) {
      Real s = gc_g[r];
      const int cmax = std::min(n - 1, r + 2*kl);
      for (int c = r + 1; c <= cmax; ++c) {s -= el(r, c)*gc_g[c];}
      gc_g[r] = s/el(r, r);
    }
  } else {
    std::fill(gc_g.begin(), gc_g.end(), 0.0);
  }
  for (int q = 0; q < n; ++q) {gc_x(q) = gc_g[q];}
  // r' = r - A P_g x_g
  const int rr = (upd == 1) ? M1_IW_KP : ((upd == 2) ? M1_IW_KS : rc);
  auto gx_ = gc_x;
  auto go_ = gc_off;
  if (nb == 1) {
    auto w_ = gcw;
    par_for("m1_gc_rp", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, nj-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int gi = go_(3*m) + i;
      Real v = iw_(m,rr,k+ks,j+js,i+is);
      for (int d = -2; d <= 2; ++d) {
        if (gi + d >= 0 && gi + d < n1) v -= w_(m,d+2,k,j,i)*gx_(gi + d);
      }
      iw_(m,M1_IW_S1,k+ks,j+js,i+is) = v;
    });
  } else {
    const int n2 = gc_n2, n3 = gc_n3, gb2 = gc_b2, gb3 = gc_b3;
    const bool per2 = gc_per2, per3 = gc_per3;
    const bool edg = st_edges, vfold = impl_vfold && vimp_now;
    const int nst = ost.extent_int(1);
    auto st_ = ost;
    par_for("m1_gc_rp", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, nj-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const int gi = go_(3*m) + i, gj = go_(3*m+1) + j, gk = go_(3*m+2) + k;
      const int b2 = (gj*gb2)/n2, b3 = (gk*gb3)/n3;
      Real v = iw_(m,rr,k+ks,j+js,i+is);
      for (int s = 0; s < nst; ++s) {
        int di, dj, dk;
        if (!GCSlot(s, thrd, edg, vfold, di, dj, dk)) continue;
        if (gi + di < 0 || gi + di >= n1) continue;
        const int c = GCClass(gj, gk, dj, dk, n2, n3, gb2, gb3, per2, per3);
        if (c < 0) continue;
        const int cb = ((b3 + c/3 - 1 + gb3)%gb3)*gb2 + (b2 + c%3 - 1 + gb2)%gb2;
        v -= st_(m,s,k+ks,j+js,i+is)*gx_((gi + di)*nb + cb);
      }
      iw_(m,M1_IW_S1,k+ks,j+js,i+is) = v;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGCAdd
//! \brief mg_gc without local levels: z += P_g x_g

void RadiationM1::ImplicitGCAdd(int zc) {
  if (!gc_ok) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1, nj = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const int nk = trans_x3 ? indcs.nx3 : 1;
  const int n2 = gc_n2, n3 = gc_n3, gb2 = gc_b2, gb3 = gc_b3, nb = gc_nb;
  auto iw_ = iw;
  auto gx_ = gc_x;
  auto go_ = gc_off;
  par_for("m1_gc_add", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, nj-1, 0, nx-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int b2 = ((go_(3*m+1) + j)*gb2)/n2;
    const int b3 = ((go_(3*m+2) + k)*gb3)/n3;
    iw_(m,zc,k+ks,j+js,i+is) += gx_((go_(3*m) + i)*nb + b3*gb2 + b2);
  });
}


//----------------------------------------------------------------------------------------
// implicit_precond = mg_gf (tests_m1/runs_5t_fast5box): the global Fourier coarse space

namespace {
constexpr int GF_IB = 16, GF_RG = 16, GF_NV = GF_IB*GF_RG;

//! p_(((m*nkc + kc)*NQ + q)*nx + i) = the sum over the k of chunk kc (width kcw) of block
//! m of what f(m, k, g, i, acc) adds to acc[0..NQ-1] for the rows j = g, g + GF_RG, ...
//! (active indices from 0).  One team per (m, kc, 16 x1 cells): 16 row groups g x 16 x1
//! cells, the row groups added in a fixed order (deterministic; the host path, team
//! size 1, sums in the same order)
template <int NQ, class F>
void GFLayerReduce(const char *name, const int nmb, const int nkc, const int kcw,
                   const int nk, const int nx, const DvceArray1D<Real> &p_, F f) {
  constexpr int IB = GF_IB, RG = GF_RG, NV = GF_NV, QB = 8;
  const int nib = (nx + IB - 1)/IB;
  const int nlg = nmb*nkc*nib;
  const bool host = std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value;
  const int ts = host ? 1 : NV;
  const int nsc = host ? NV*NQ : NV*QB;
  size_t scr = ScrArray1D<Real>::shmem_size(nsc);
  using Pol = Kokkos::TeamPolicy<DevExeSpace, Kokkos::LaunchBounds<NV,1>>;
  using TM = typename Pol::member_type;
  Pol pol(DevExeSpace(), nlg, ts);
  Kokkos::parallel_for(name, pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
  KOKKOS_LAMBDA(const TM &tm) {
    int l = tm.league_rank();
    const int ib = l%nib; l /= nib;
    const int kc = l%nkc;
    const int m = l/nkc;
    const int k0 = kc*kcw;
    const int k1 = (k0 + kcw < nk) ? (k0 + kcw) : nk;
    const size_t ob = (static_cast<size_t>(m)*nkc + kc)*NQ;
    ScrArray1D<Real> sw(tm.team_scratch(0), nsc);
    if (tm.team_size() == NV) {
      const int t = tm.team_rank();
      const int il = t%IB, g = t/IB, i = ib*IB + il;
      Real acc[NQ];
      for (int q = 0; q < NQ; ++q) {acc[q] = 0.0;}
      if (i < nx) {
        for (int k = k0; k < k1; ++k) {f(m, k, g, i, acc);}
      }
      for (int q0 = 0; q0 < NQ; q0 += QB) {
        for (int qq = 0; qq < QB; ++qq) {
          if (q0 + qq < NQ) {sw(qq*NV + t) = acc[q0 + qq];}
        }
        tm.team_barrier();
        if (t < QB*IB) {
          const int qq = t/IB, i2 = ib*IB + t%IB;
          if (q0 + qq < NQ && i2 < nx) {
            Real a = 0.0;
            for (int g2 = 0; g2 < RG; ++g2) {a += sw(qq*NV + g2*IB + t%IB);}
            p_((ob + q0 + qq)*nx + i2) = a;
          }
        }
        tm.team_barrier();
      }
    } else {
      for (int t = 0; t < NV; ++t) {
        const int il = t%IB, g = t/IB, i = ib*IB + il;
        Real acc[NQ];
        for (int q = 0; q < NQ; ++q) {acc[q] = 0.0;}
        if (i < nx) {
          for (int k = k0; k < k1; ++k) {f(m, k, g, i, acc);}
        }
        for (int q = 0; q < NQ; ++q) {sw(q*NV + t) = acc[q];}
      }
      for (int q = 0; q < NQ; ++q) {
        for (int il = 0; il < IB; ++il) {
          const int i2 = ib*IB + il;
          if (i2 >= nx) continue;
          Real a = 0.0;
          for (int g2 = 0; g2 < RG; ++g2) {a += sw(q*NV + g2*IB + il);}
          p_((ob + q)*nx + i2) = a;
        }
      }
    }
  });
}

//! out_(gi*nq + q) = the sum over the blocks m holding layer gi (in order) and the k
//! chunks kc (in order) of p_(((m*nkc + kc)*nq + q)*nx + gi - x1 offset of m)
void GFLayerSum(const int nmb, const int nkc, const int nq, const int nx, const int n1,
                const DvceArray1D<int> &go_, const DvceArray1D<Real> &p_,
                const DvceArray1D<Real> &out_) {
  Kokkos::parallel_for("m1_gf_sum", Kokkos::RangePolicy<DevExeSpace>(0, n1*nq),
  KOKKOS_LAMBDA(const int idx) {
    const int q = idx/n1, gi = idx%n1;
    Real s = 0.0;
    for (int m = 0; m < nmb; ++m) {
      const int li = gi - go_(3*m);
      if (li < 0 || li >= nx) continue;
      const size_t o = (static_cast<size_t>(m)*nkc*nq + q)*nx + li;
      const size_t st = static_cast<size_t>(nq)*nx;
      Real a[4] = {0.0, 0.0, 0.0, 0.0};   // independent partial sums, fixed order
      int kc = 0;
      for (; kc + 3 < nkc; kc += 4) {
        a[0] += p_(o + kc*st);
        a[1] += p_(o + (kc + 1)*st);
        a[2] += p_(o + (kc + 2)*st);
        a[3] += p_(o + (kc + 3)*st);
      }
      for (; kc < nkc; ++kc) {a[0] += p_(o + kc*st);}
      s += (a[0] + a[1]) + (a[2] + a[3]);
    }
    out_(gi*nq + q) = s;
  });
}

//! the layer projections g(gi, c*NB2 + a) = sum_k v_c(gk) sum_j u_a(gj) r(gi,j,k), with
//! the BiCGStab p / s update when up > 0 (as ImplicitPCRSolveX)
template <int NB2, int NB3>
void GFProject(const int nmb, const int nkc, const int kcw, const int nk, const int nj,
               const int nx, const int is, const int js, const int ks, const int n2,
               const int n3, const int up, const int cr, const Real a1, const Real a2,
               const DvceArray5D<Real> &iw_, const DvceArray1D<int> &go_,
               const DvceArray1D<Real> &u_, const DvceArray1D<Real> &v_,
               const DvceArray1D<Real> &p_) {
  GFLayerReduce<NB2*NB3>("m1_gf_red", nmb, nkc, kcw, nk, nx, p_,
  KOKKOS_LAMBDA(const int m, const int k, const int g, const int i, Real *acc) {
    const int kk = k + ks, ii = i + is;
    const int o2 = go_(3*m+1), gk = go_(3*m+2) + k;
    Real b[NB2];
    for (int a = 0; a < NB2; ++a) {b[a] = 0.0;}
    for (int j = g; j < nj; j += GF_RG) {
      const int jj = j + js;
      Real v;
      if (up == 1) {
        v = iw_(m,M1_IW_KR,kk,jj,ii) + a1*(iw_(m,M1_IW_KP,kk,jj,ii)
                                           - a2*iw_(m,M1_IW_KV,kk,jj,ii));
        iw_(m,M1_IW_KP,kk,jj,ii) = v;
      } else if (up == 2) {
        v = iw_(m,M1_IW_KR,kk,jj,ii) - a1*iw_(m,M1_IW_KV,kk,jj,ii);
        iw_(m,M1_IW_KS,kk,jj,ii) = v;
      } else {
        v = iw_(m,cr,kk,jj,ii);
      }
      for (int a = 0; a < NB2; ++a) {b[a] += u_(a*n2 + o2 + j)*v;}
    }
    for (int c = 0; c < NB3; ++c) {
      const Real vc = v_(c*n3 + gk);
      for (int a = 0; a < NB2; ++a) {acc[c*NB2 + a] += vc*b[a];}
    }
  });
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGFInit
//! \brief implicit_precond = mg_gf: checks, the blocks' global offsets, the 1-D Fourier
//! bases and their correlation sums, the storage (once)

void RadiationM1::ImplicitGFInit() {
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const bool thrd = trans_x3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx = indcs.nx1, nx2 = indcs.nx2, nx3 = thrd ? indcs.nx3 : 1;
  if (pm->multilevel) {
    GCFatal("<rad_m1>/implicit_precond = mg_gf needs a single-level mesh");
  }
  if (!impl_stencil) {
    GCFatal("<rad_m1>/implicit_precond = mg_gf needs implicit_op_stencil");
  }
  gc_n1 = pm->mesh_indcs.nx1;
  gc_n2 = pm->mesh_indcs.nx2;
  gc_n3 = thrd ? pm->mesh_indcs.nx3 : 1;
  gc_per2 = (pm->mesh_bcs[static_cast<int>(BoundaryFace::inner_x2)]
             == BoundaryFlag::periodic);
  gc_per3 = (pm->mesh_bcs[static_cast<int>(BoundaryFace::inner_x3)]
             == BoundaryFlag::periodic);
  if (!gc_per2 || (thrd && !gc_per3)) {
    GCFatal("<rad_m1>/implicit_precond = mg_gf needs periodic x2 (and x3) boundaries");
  }
  const int K = gf_k;
  if (2*K >= gc_n2 || (thrd && 2*K >= gc_n3)) {
    GCFatal("<rad_m1>/implicit_gf_modes: 2 K must be below the mesh nx2 (and nx3)");
  }
  if (gc_n1 > 1024) {
    GCFatal("<rad_m1>/implicit_precond = mg_gf: at most 1024 x1 cells");
  }
  gf_nb2 = 2*K + 1;
  gf_nb3 = thrd ? 2*K + 1 : 1;
  gf_nm = gf_nb2*gf_nb3;
  pmy_pack->pmb->mb_gid.sync_host();
  gc_off_h.assign(3*nmb, 0);
  for (int m = 0; m < nmb; ++m) {
    auto &ll = pm->lloc_eachmb[pmy_pack->pmb->mb_gid.h_view(m)];
    gc_off_h[3*m] = static_cast<int>(ll.lx1)*nx;
    gc_off_h[3*m+1] = static_cast<int>(ll.lx2)*nx2;
    gc_off_h[3*m+2] = thrd ? static_cast<int>(ll.lx3)*nx3 : 0;
  }
  gc_off = DvceArray1D<int>("m1_gf_off", 3*nmb);
  {
    auto h1 = Kokkos::create_mirror_view(gc_off);
    for (int q = 0; q < 3*nmb; ++q) {h1(q) = gc_off_h[q];}
    Kokkos::deep_copy(gc_off, h1);
  }
  // 1-D bases at the cell centres: 1, cos(2 pi k x), sin(2 pi k x), k = 1..K, and
  // their correlation sums C(a, d) = sum_j u_a(j) u_a(j + d) (periodic), d = -2..2
  std::vector<Real> cu, cv;
  auto basis = [&](const int n, const int nb, DvceArray1D<Real> &d, DvceArray1D<Real> &dc,
                   std::vector<Real> &cc) {
    std::vector<Real> b(static_cast<size_t>(nb)*n, 1.0);
    const Real tp = 8.0*std::atan(1.0);
    for (int a = 1; a < nb; ++a) {
      const int kk = (a + 1)/2;
      for (int j = 0; j < n; ++j) {
        const Real ph = tp*kk*(j + 0.5)/n;
        b[a*n + j] = (a % 2 == 1) ? std::cos(ph) : std::sin(ph);
      }
    }
    cc.assign(5*nb, 0.0);
    for (int a = 0; a < nb; ++a) {
      for (int dd = -2; dd <= 2; ++dd) {
        Real s = 0.0;
        for (int j = 0; j < n; ++j) {s += b[a*n + j]*b[a*n + ((j + dd)%n + n)%n];}
        cc[5*a + dd + 2] = s;
      }
    }
    d = DvceArray1D<Real>("m1_gf_basis", nb*n);
    auto h = Kokkos::create_mirror_view(d);
    for (int q = 0; q < nb*n; ++q) {h(q) = b[q];}
    Kokkos::deep_copy(d, h);
    dc = DvceArray1D<Real>("m1_gf_corr", 5*nb);
    auto hc = Kokkos::create_mirror_view(dc);
    for (int q = 0; q < 5*nb; ++q) {hc(q) = cc[q];}
    Kokkos::deep_copy(dc, hc);
  };
  basis(gc_n2, gf_nb2, gf_u, gf_cud, cu);
  basis(gc_n3, gf_nb3, gf_v, gf_cvd, cv);
  gf_nrm = DvceArray1D<Real>("m1_gf_nrm", gf_nm);
  {
    auto h = Kokkos::create_mirror_view(gf_nrm);
    for (int c = 0; c < gf_nb3; ++c) {
      for (int a = 0; a < gf_nb2; ++a) {h(c*gf_nb2 + a) = 1.0/(cu[5*a + 2]*cv[5*c + 2]);}
    }
    Kokkos::deep_copy(gf_nrm, h);
  }
  // k chunks of the layer reductions: ~192 rows per team
  gf_kcw = std::max(1, std::min(nx3, 192/std::max(nx2, 1)));
  gf_nkc = (nx3 + gf_kcw - 1)/gf_kcw;
  const int nqm = std::max(gf_nm, 25), n1 = gc_n1;
  const size_t npart = static_cast<size_t>(nmb)*std::max(gf_nkc*gf_nm, nx3*25)*nx;
  gf_part = DvceArray1D<Real>("m1_gf_part", npart);
  gf_gh = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("m1_gf_gh", nqm*gc_n1);
  gf_g = DvceArray1D<Real>("m1_gf_g", nqm*n1);
  gf_gn = DvceArray1D<Real>("m1_gf_gn", gf_nm*n1);
  gf_x = DvceArray1D<Real>("m1_gf_x", gf_nm*n1);
  gf_hg = DvceArray1D<Real>("m1_gf_hg", static_cast<size_t>(gf_nb3)*gc_n2*n1);
  gf_hx = DvceArray1D<Real>("m1_gf_hx", static_cast<size_t>(gf_nb3)*gc_n2*n1);
  gf_ainv = DvceArray1D<Real>("m1_gf_ainv", static_cast<size_t>(gf_nm)*n1*n1);
  gf_okd = DvceArray1D<int>("m1_gf_okd", gf_nm);
  gf_ok = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGFBuild
//! \brief the mode-diagonal coarse matrices of this pass, A_q(gi, gi+di) =
//! sum_s Sbar_s(gi)/(n2 n3) C2(a, dj_s) C3(c, dk_s), from the layer sums Sbar of the
//! stencil slots (the coefficients of the Krylov operator: the Galerkin operator of the
//! layer-mean stencil without the cos-sin couplings), and their banded LU with partial
//! pivoting (device, one team per mode; a singular mode is left out, gf_okd = 0)

void RadiationM1::ImplicitGFBuild() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1, nj = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  const int nk = thrd ? indcs.nx3 : 1;
  const int n1 = gc_n1;
  const bool edg = st_edges, vfold = impl_vfold && vimp_now;
  const int nst = ost.extent_int(1);
  auto st_ = ost;
  auto go_ = gc_off;
  GFLayerReduce<25>("m1_gf_sbar", nmb, nk, 1, nk, nx, gf_part,
  KOKKOS_LAMBDA(const int m, const int k, const int g, const int i, Real *acc) {
    const int gi = go_(3*m) + i;
    for (int s = 0; s < 25; ++s) {
      int di, dj, dk;
      if (s >= nst || !GCSlot(s, thrd, edg, vfold, di, dj, dk)) continue;
      if (gi + di < 0 || gi + di >= n1) continue;
      Real a = 0.0;
      for (int j = g; j < nj; j += GF_RG) {a += st_(m,s,k+ks,j+js,i+is);}
      acc[s] += a;
    }
  });
  GFLayerSum(nmb, nk, 25, nx, n1, go_, gf_part, gf_g);
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    auto sb = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), gf_g);
    MPI_Allreduce(MPI_IN_PLACE, sb.data(), 25*n1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    Kokkos::deep_copy(gf_g, sb);
  }
#endif
  const int nb2 = gf_nb2;
  const Real inv = 1.0/(static_cast<Real>(gc_n2)*gc_n3);
  auto sb_ = gf_g;
  auto cu_ = gf_cud;
  auto cv_ = gf_cvd;
  auto ok_ = gf_okd;
  const int ts = std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value
                 ? 1 : 128;
  size_t scr = ScrArray1D<Real>::shmem_size(7*n1) + ScrArray1D<int>::shmem_size(n1);
  auto ai_ = gf_ainv;
  Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), gf_nm, ts);
  Kokkos::parallel_for("m1_gf_fac", pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
  KOKKOS_LAMBDA(const TeamMember_t &tm) {
    const int q = tm.league_rank();
    const int a = q%nb2, c = q/nb2;
    ScrArray1D<Real> el(tm.team_scratch(0), 7*n1);   // (r, cc) at r*7 + cc - r + 2
    ScrArray1D<int> pvs(tm.team_scratch(0), n1);
    for (int r = tm.team_rank(); r < n1; r += tm.team_size()) {
      Real w[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (int s = 0; s < 25; ++s) {
        int di, dj, dk;
        if (s >= nst || !GCSlot(s, thrd, edg, vfold, di, dj, dk)) continue;
        if (r + di < 0 || r + di >= n1) continue;
        w[di + 2] += sb_(r*25 + s)*inv*cu_(5*a + dj + 2)*cv_(5*c + dk + 2);
      }
      for (int e = 0; e < 7; ++e) {el(r*7 + e) = (e < 5) ? w[e] : 0.0;}
    }
    tm.team_barrier();
    // banded LU with partial pivoting (the multipliers stay where they are made, as
    // the host LU of mg_gc), by one thread with the three active rows in registers:
    // w[t][e] = element (cc + t, cc + e), e = 0..4
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      auto band = [&](const int r, const int cc) -> Real {   // original row r, column cc
        const int o = cc - r + 2;
        return (r < n1 && cc < n1 && o >= 0 && o < 5) ? el(r*7 + o) : 0.0;
      };
      Real w0[5], w1[5], w2[5];
      for (int e = 0; e < 5; ++e) {
        w0[e] = band(0, e);
        w1[e] = band(1, e);
        w2[e] = band(2, e);
      }
      int good = 1;
      for (int cc = 0; cc < n1; ++cc) {
        const int nr = (n1 - 1 - cc < 2) ? n1 - 1 - cc : 2;   // rows below the pivot
        int p = 0;
        if (nr >= 1 && fabs(w1[0]) > fabs(w0[0])) p = 1;
        if (nr >= 2 && fabs(w2[0]) > fabs((p == 1) ? w1[0] : w0[0])) p = 2;
        pvs(cc) = cc + p;
        for (int e = 0; e < 5; ++e) {
          const Real t = w0[e];
          if (p == 1) {w0[e] = w1[e]; w1[e] = t;}
          if (p == 2) {w0[e] = w2[e]; w2[e] = t;}
        }
        if (!(fabs(w0[0]) > 0.0) || !Kokkos::isfinite(w0[0])) {good = 0; break;}
        const Real ip = 1.0/w0[0];
        const Real f1 = (nr >= 1) ? w1[0]*ip : 0.0;
        const Real f2 = (nr >= 2) ? w2[0]*ip : 0.0;
        for (int e = 1; e < 5; ++e) {
          if (f1 != 0.0) w1[e] -= f1*w0[e];
          if (f2 != 0.0) w2[e] -= f2*w0[e];
        }
        for (int e = 0; e < 5; ++e) {el(cc*7 + 2 + e) = w0[e];}   // U row cc
        if (nr >= 1) el((cc + 1)*7 + 1) = f1;                     // multipliers
        if (nr >= 2) el((cc + 2)*7 + 0) = f2;
        for (int e = 0; e < 4; ++e) {w0[e] = w1[e + 1]; w1[e] = w2[e + 1];}
        w0[4] = 0.0; w1[4] = 0.0;
        for (int e = 0; e < 5; ++e) {w2[e] = band(cc + 3, cc + 1 + e);}
      }
      ok_(q) = good;
    });
    tm.team_barrier();
    // the inverse, one column per thread: ai_((q*n1 + col)*n1 + r) = (A_q^{-1})(r, col),
    // the forward / back substitution of ImplicitGFPre's former per-application solve
    const bool okq = (ok_(q) != 0);
    for (int col = tm.team_rank(); col < n1; col += tm.team_size()) {
      const size_t oa = (static_cast<size_t>(q)*n1 + col)*n1;
      if (!okq) {
        for (int r = 0; r < n1; ++r) {ai_(oa + r) = 0.0;}
        continue;
      }
      Real w0 = (col == 0) ? 1.0 : 0.0;
      Real w1 = (col == 1) ? 1.0 : 0.0;
      Real w2 = (col == 2) ? 1.0 : 0.0;
      for (int cc = 0; cc < n1; ++cc) {
        const int p = pvs(cc) - cc;
        if (p == 1) {const Real t = w0; w0 = w1; w1 = t;}
        if (p == 2) {const Real t = w0; w0 = w2; w2 = t;}
        if (cc + 1 < n1) w1 -= el((cc + 1)*7 + 1)*w0;
        if (cc + 2 < n1) w2 -= el((cc + 2)*7 + 0)*w0;
        ai_(oa + cc) = w0;
        w0 = w1;
        w1 = w2;
        w2 = (cc + 3 == col) ? 1.0 : 0.0;
      }
      Real x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0;
      for (int r = n1 - 1; r >= 0; --r) {
        const int o = r*7 + 2;
        const Real xr = (ai_(oa + r) - el(o + 1)*x1 - el(o + 2)*x2 - el(o + 3)*x3
                         - el(o + 4)*x4)/el(o);
        ai_(oa + r) = xr;
        x4 = x3; x3 = x2; x2 = x1; x1 = xr;
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGFPre
//! \brief mg_gf, before mg: the right-hand side (made as the BiCGStab p / s update when
//! upd > 0), its layer projections g = P^T r (one MPI_Allreduce on several ranks),
//! x_g = A_q^{-1} g per mode (one team per mode, the banded solve in scratch), the
//! tables sum_a u_a (g/|phi|^2, x_g), and r' = r - P diag(1/|phi|^2) g into M1_IW_S1

void RadiationM1::ImplicitGFPre(int rc, int upd, Real c1, Real c2) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1, nj = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const int nk = trans_x3 ? indcs.nx3 : 1;
  const int n1 = gc_n1, n2 = gc_n2, n3 = gc_n3, nm = gf_nm, nb2 = gf_nb2, nb3 = gf_nb3;
  const int nkc = gf_nkc, kcw = gf_kcw;
  auto iw_ = iw;
  auto go_ = gc_off;
  auto u_ = gf_u;
  auto v_ = gf_v;
  auto p_ = gf_part;
  const int key = nb2*10 + nb3;
#define M1_GF_PROJ(A, B) GFProject<A, B>(nmb, nkc, kcw, nk, nj, nx, is, js, ks, n2, n3, \
                                         upd, rc, c1, c2, iw_, go_, u_, v_, p_)
  switch (key) {
    case 11: M1_GF_PROJ(1, 1); break;
    case 31: M1_GF_PROJ(3, 1); break;
    case 33: M1_GF_PROJ(3, 3); break;
    case 51: M1_GF_PROJ(5, 1); break;
    case 55: M1_GF_PROJ(5, 5); break;
    case 71: M1_GF_PROJ(7, 1); break;
    default: M1_GF_PROJ(7, 7); break;
  }
#undef M1_GF_PROJ
  bool mpi = false;
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    mpi = true;
    auto gh_ = gf_gh;
    Kokkos::parallel_for("m1_gf_sum", Kokkos::RangePolicy<DevExeSpace>(0, n1*nm),
    KOKKOS_LAMBDA(const int idx) {
      const int q = idx/n1, gi = idx%n1;
      Real s = 0.0;
      for (int m = 0; m < nmb; ++m) {
        const int li = gi - go_(3*m);
        if (li < 0 || li >= nx) continue;
        for (int kc = 0; kc < nkc; ++kc) {
          s += p_(((static_cast<size_t>(m)*nkc + kc)*nm + q)*nx + li);
        }
      }
      gh_(gi*nm + q) = s;
    });
    Kokkos::fence();
    MPI_Allreduce(MPI_IN_PLACE, gf_gh.data(), nm*n1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
  }
#endif
  // x_g per mode: one team per mode; the layer sums (in the order of GFLayerSum), then
  // the banded forward / back substitution by one thread in scratch
  {
    auto g_ = gf_gh;
    auto gn_ = gf_gn;
    auto x_ = gf_x;
    auto nr_ = gf_nrm;
    auto ok_ = gf_okd;
    const int ts = std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value
                   ? 1 : 128;
    auto ai_ = gf_ainv;
    size_t scr = ScrArray1D<Real>::shmem_size(n1);
    Kokkos::TeamPolicy<DevExeSpace> pol(DevExeSpace(), nm, ts);
    Kokkos::parallel_for("m1_gf_solve", pol.set_scratch_size(0, Kokkos::PerTeam(scr)),
    KOKKOS_LAMBDA(const TeamMember_t &tm) {
      const int q = tm.league_rank();
      ScrArray1D<Real> y(tm.team_scratch(0), n1);
      const bool ok = (ok_(q) != 0);
      const Real nq = nr_(q);
      for (int gi = tm.team_rank(); gi < n1; gi += tm.team_size()) {
        Real s = 0.0;
        if (mpi) {
          s = g_(gi*nm + q);
        } else {
          // independent partial sums (latency: the partials are L2 reads), fixed order
          for (int m = 0; m < nmb; ++m) {
            const int li = gi - go_(3*m);
            if (li < 0 || li >= nx) continue;
            const size_t o = (static_cast<size_t>(m)*nkc*nm + q)*nx + li;
            const size_t st = static_cast<size_t>(nm)*nx;
            Real a[4] = {0.0, 0.0, 0.0, 0.0};
            int kc = 0;
            for (; kc + 3 < nkc; kc += 4) {
              a[0] += p_(o + kc*st);
              a[1] += p_(o + (kc + 1)*st);
              a[2] += p_(o + (kc + 2)*st);
              a[3] += p_(o + (kc + 3)*st);
            }
            for (; kc < nkc; ++kc) {a[0] += p_(o + kc*st);}
            s += (a[0] + a[1]) + (a[2] + a[3]);
          }
        }
        y(gi) = s;
        gn_(gi*nm + q) = ok ? s*nq : 0.0;
      }
      tm.team_barrier();
      for (int gi = tm.team_rank(); gi < n1; gi += tm.team_size()) {
        Real s = 0.0;
        const size_t oa = static_cast<size_t>(q)*n1*n1 + gi;
        Real a[4] = {0.0, 0.0, 0.0, 0.0};
        int c = 0;
        for (; c + 3 < n1; c += 4) {
          a[0] += ai_(oa + static_cast<size_t>(c)*n1)*y(c);
          a[1] += ai_(oa + static_cast<size_t>(c + 1)*n1)*y(c + 1);
          a[2] += ai_(oa + static_cast<size_t>(c + 2)*n1)*y(c + 2);
          a[3] += ai_(oa + static_cast<size_t>(c + 3)*n1)*y(c + 3);
        }
        for (; c < n1; ++c) {a[0] += ai_(oa + static_cast<size_t>(c)*n1)*y(c);}
        s = (a[0] + a[1]) + (a[2] + a[3]);
        x_(gi*nm + q) = s;
      }
    });
  }
  // the tables hg / hx ((c*n2 + gj)*n1 + gi) = sum_a u_a(gj) (gn, x)(gi, c*nb2 + a)
  {
    auto gn_ = gf_gn;
    auto x_ = gf_x;
    auto hg_ = gf_hg;
    auto hx_ = gf_hx;
    Kokkos::parallel_for("m1_gf_h", Kokkos::RangePolicy<DevExeSpace>(0, nb3*n2*n1),
    KOKKOS_LAMBDA(const int idx) {
      const int gi = idx%n1, cj = idx/n1;
      const int gj = cj%n2, c = cj/n2;
      Real sg = 0.0, sx = 0.0;
      for (int a = 0; a < nb2; ++a) {
        const Real ua = u_(a*n2 + gj);
        sg += ua*gn_(gi*nm + c*nb2 + a);
        sx += ua*x_(gi*nm + c*nb2 + a);
      }
      hg_(idx) = sg;
      hx_(idx) = sx;
    });
  }
  // r' = r - P diag(1/|phi_q|^2) P^T r
  const int rr = (upd == 1) ? M1_IW_KP : ((upd == 2) ? M1_IW_KS : rc);
  auto hg_ = gf_hg;
  par_for("m1_gf_rp", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, nj-1, 0, nx-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int gi = go_(3*m) + i, gj = go_(3*m+1) + j, gk = go_(3*m+2) + k;
    Real g = 0.0;
    for (int c = 0; c < nb3; ++c) {g += v_(c*n3 + gk)*hg_((c*n2 + gj)*n1 + gi);}
    iw_(m,M1_IW_S1,k+ks,j+js,i+is) = iw_(m,rr,k+ks,j+js,i+is) - g;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGFAdd
//! \brief mg_gf without local levels: z += P x_g

void RadiationM1::ImplicitGFAdd(int zc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1, nj = indcs.nx2;
  const int nmb = pmy_pack->nmb_thispack;
  const int nk = trans_x3 ? indcs.nx3 : 1;
  const int n1 = gc_n1, n2 = gc_n2, n3 = gc_n3, nb3 = gf_nb3;
  auto iw_ = iw;
  auto hx_ = gf_hx;
  auto go_ = gc_off;
  auto v_ = gf_v;
  par_for("m1_gf_add", DevExeSpace(), 0, nmb-1, 0, nk-1, 0, nj-1, 0, nx-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const int gi = go_(3*m) + i, gj = go_(3*m+1) + j, gk = go_(3*m+2) + k;
    Real g = 0.0;
    for (int c = 0; c < nb3; ++c) {g += v_(c*n3 + gk)*hx_((c*n2 + gj)*n1 + gi);}
    iw_(m,zc,k+ks,j+js,i+is) += g;
  });
}

} // namespace radm1
