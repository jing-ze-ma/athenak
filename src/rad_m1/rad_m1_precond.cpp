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
//! implicit_mg_fuse = true (default; tests_m1/runs_5n_mgfuse): the restricted residual
//! of level l is formed in the load phase of level l's own line sweeps instead of in a
//! separate kernel after level l-1 -- the same numbers, 1 launch fewer per level.


#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

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

void MGFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

//! the restricted fine residual R (r - A7 z0) of level-1 cell (kc,jc,i): the 2 x 2 fine
//! columns of the aggregate, with (hal) or without the halo of z0
KOKKOS_INLINE_FUNCTION
Real MGRes0(const DvceArray5D<Real> &iw_, const int m, const int kc, const int jc,
            const int i, const int is, const int js, const int ks, const int nx,
            const int fj, const int fk, const bool thrd, const bool hal, const int rr,
            const int zc) {
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
  return sum;
}

//! the restricted residual of level l-1 (array c_, size cj x ck, x3 couplings th) at
//! level-l cell (kc,jc,i)
KOKKOS_INLINE_FUNCTION
Real MGResL(const DvceArray5D<Real> &c_, const int m, const int kc, const int jc,
            const int i, const int nx, const int cj, const int ck, const bool th) {
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
  return sum;
}

//! where the right-hand side of a coarse line sweep comes from (implicit_mg_fuse):
//! mode 0 = MG_R as stored by a separate restriction kernel; mode 1 = restricted from the
//! fine rows (MGRes0) and mode 2 = from level l-1 (MGResL), both in the load phase of the
//! sweep and stored to MG_R for the next level's restriction.  The value is the very
//! number the separate kernel stores, so the fused and unfused maps agree.
struct MGRes {
  int mode = 0;
  DvceArray5D<Real> iw;
  int rr = 0, zc = 0, is = 0, js = 0, ks = 0;
  bool hal = false, fthrd = false;
  DvceArray5D<Real> f;
  int fj = 1, fk = 1;
  bool fth = false;
};

//! one colour of the forward red-black line sweep on a coarse level: columns (k,j) with
//! parity (k+j) = col solve A_line z = r (- the 5-point coupling to z of the other
//! colour when sub), by parallel cyclic reduction in team scratch (as M1PCRX)
template <typename T>
void M1PCRMG(const DvceArray5D<Real> &a_, const int nmb, const int nx, const int nj,
             const int nk, const int col, const bool sub, const bool thrd, int ts,
             const MGRes &rs) {
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
  const int mode = rs.mode;
  auto riw = rs.iw;
  auto rf = rs.f;
  const int rrc = rs.rr, rzc = rs.zc, ris = rs.is, rjs = rs.js, rks = rs.ks;
  const bool rhal = rs.hal, rfthrd = rs.fthrd, rfth = rs.fth;
  const int rfj = rs.fj, rfk = rs.fk;
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
      Real rr;
      if (mode == 1) {
        rr = MGRes0(riw, m, k, j, i, ris, rjs, rks, nx, rfj, rfk, rfthrd, rhal, rrc,
                    rzc);
        a_(m,MG_R,k,j,i) = rr;
      } else if (mode == 2) {
        rr = MGResL(rf, m, k, j, i, nx, rfj, rfk, rfth);
        a_(m,MG_R,k,j,i) = rr;
      } else {
        rr = a_(m,MG_R,k,j,i);
      }
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
//! \fn void RadiationM1::ImplicitMGInit
//! \brief the implicit_precond = mg keys, read only under mg (every other configuration
//! untouched)

void RadiationM1::ImplicitMGInit(ParameterInput *pin) {
  mg_nlev = 0;
  mg_halo = true;
  mg_fuse = true;
  if (impl_prec != 3) {return;}
  mg_nlev = pin->GetOrAddInteger("rad_m1","implicit_mg_levels",2);
  mg_halo = pin->GetOrAddBoolean("rad_m1","implicit_mg_halo",true);
  mg_fuse = pin->GetOrAddBoolean("rad_m1","implicit_mg_fuse",true);
  if (mg_nlev < 2) {MGFatal("<rad_m1>/implicit_mg_levels must be >= 2");}
}

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
          if (j + 1 <= 2*jc + 1 && j + 1 < fj) {tb += v[4];}
          else if (j + 1 < fj) {jp += v[4];}
          if (thrd) {
            if (k - 1 >= 2*kc) {tb += v[5];} else if (k - 1 >= 0) {km += v[5];}
            if (k + 1 <= 2*kc + 1 && k + 1 < fk) {tb += v[6];}
            else if (k + 1 < fk) {kp += v[6];}
          }
        }
      }
      c_(m,MG_TA,kc,jc,i) = ta; c_(m,MG_TB,kc,jc,i) = tb; c_(m,MG_TC,kc,jc,i) = tc;
      c_(m,MG_CJM,kc,jc,i) = jm; c_(m,MG_CJP,kc,jc,i) = jp;
      c_(m,MG_CKM,kc,jc,i) = km; c_(m,MG_CKP,kc,jc,i) = kp;
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitMGApply
//! \brief implicit_precond = mg: z = M^{-1} r (see the file header); rc / upd / c1 / c2
//! as ImplicitPCRSolveX

void RadiationM1::ImplicitMGApply(int rc, int zc, int upd, Real c1, Real c2) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int nx = indcs.nx1;
  const int nmb = pmy_pack->nmb_thispack;
  const bool thrd = trans_x3;
  auto iw_ = iw;
  // the fine smoother: rbgs_fwd, exactly the implicit_precond = rbgs_fwd map
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 0, -1);
  ImplicitPCRSolveX(rc, zc, upd, c1, c2, 1, zc);
  const int nl = static_cast<int>(mgc.size());
  if (nl < 2) {return;}
  const int rr = (upd == 1) ? M1_IW_KP : ((upd == 2) ? M1_IW_KS : rc);
  const bool hal = mg_halo;
  if (hal) {ImplicitKrylovHalo(zc);}
  int ts = impl_pcr_team;
  if (ts == 0) {
    ts = 1;
    while (ts < nx && ts < 256) ts *= 2;
  }
  const bool fuse = mg_fuse;
  // the fine residual r - A7 z0, restricted (unfused: its own kernel)
  if (!fuse) {
    auto c_ = mgc[1];
    const int fj = mg_nj[0], fk = mg_nk[0];
    const int cj = mg_nj[1], ck = mg_nk[1];
    par_for("m1_mg_res0", DevExeSpace(), 0, nmb-1, 0, ck-1, 0, cj-1, 0, nx-1,
    KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int i) {
      c_(m,MG_R,kc,jc,i) = MGRes0(iw_, m, kc, jc, i, is, js, ks, nx, fj, fk, thrd, hal,
                                  rr, zc);
    });
  }
  // down the coarse levels: one forward red-black sweep each, the restricted residual
  // of the next level either in the next level's sweep (fused) or its own kernel
  for (int l = 1; l < nl; ++l) {
    auto c_ = mgc[l];
    const int cj = mg_nj[l], ck = mg_nk[l];
    const bool th = thrd && (ck > 1);
    MGRes rs;
    if (fuse && l == 1) {
      rs.mode = 1; rs.iw = iw_; rs.rr = rr; rs.zc = zc; rs.is = is; rs.js = js;
      rs.ks = ks; rs.hal = hal; rs.fthrd = thrd; rs.fj = mg_nj[0]; rs.fk = mg_nk[0];
    } else if (fuse) {
      rs.mode = 2; rs.f = mgc[l-1]; rs.fj = mg_nj[l-1]; rs.fk = mg_nk[l-1];
      rs.fth = thrd && (mg_nk[l-1] > 1);
    }
    if (impl_prec_float) {
      M1PCRMG<float>(c_, nmb, nx, cj, ck, 0, false, th, ts, rs);
      M1PCRMG<float>(c_, nmb, nx, cj, ck, 1, true, th, ts, rs);
    } else {
      M1PCRMG<Real>(c_, nmb, nx, cj, ck, 0, false, th, ts, rs);
      M1PCRMG<Real>(c_, nmb, nx, cj, ck, 1, true, th, ts, rs);
    }
    if (!fuse && l + 1 < nl) {
      auto n_ = mgc[l+1];
      const int nj2 = mg_nj[l+1], nk2 = mg_nk[l+1];
      par_for("m1_mg_res", DevExeSpace(), 0, nmb-1, 0, nk2-1, 0, nj2-1, 0, nx-1,
      KOKKOS_LAMBDA(const int m, const int kc, const int jc, const int i) {
        n_(m,MG_R,kc,jc,i) = MGResL(c_, m, kc, jc, i, nx, cj, ck, th);
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
  par_for("m1_mg_pro0", DevExeSpace(), 0, nmb-1, 0, mg_nk[0]-1, 0, mg_nj[0]-1, 0, nx-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real s = c1_(m,MG_Z,k >> 1,j >> 1,i);
    if (nlc > 1) s += c2_(m,MG_Z,k >> 2,j >> 2,i);
    if (nlc > 2) s += c3_(m,MG_Z,k >> 3,j >> 3,i);
    if (nlc > 3) s += c4_(m,MG_Z,k >> 4,j >> 4,i);
    iw_(m,zc,k+ks,j+js,i+is) += s;
  });
}

} // namespace radm1
