//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_rt_driver.cpp
//! \brief the two-stream RT drivers: picket_fence_two_stream_RT (+ _body), the ck
//! implicit column step (CkImplStep, CkWarmSeed), the flux conservation /
//! non-convergence diagnostics and the ck_impl_every cadence (CkCadStore, CkCadLin,
//! CkCadStep).
//! tusplit-1102: moved verbatim out of the two-stream headers (which keep the
//! declarations and the doc comments) so that they are compiled once, in their
//! own translation unit.

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>

#include "utils/two_stream_rt.hpp"
#include "utils/ck_2nd.hpp"     // ck2nd-1102 (keys read only when named)
#include "pgen/pgen.hpp"

namespace two_stream_rt {

void CkWarmSeed(Mesh *pm, DvceArray5D<Real> u0) {
  if (ck_dep_ptr == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  auto dep_ = *ck_dep_ptr;
  auto sd_ = *ck_seed_ptr;
  auto ei_ = *ck_ei_ptr;
  // ck2nd-1102 problem/ck_t2_seed: the stage predictor instead (full calls only)
  const bool ext = ck_seed_ext_on && !ck_cad_partial && (ck_seed_ext_ptr != nullptr);
  const bool warm = ck_impl_warm || ext;
  auto sx_ = ext ? *ck_seed_ext_ptr : dep_;
  const Real cap = ck_impl_dtmax;
  par_for("ck_warm_seed", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real s = 0.0;
    if (warm) {
      const Real e = ei_(m,k,j,i);          // the previous call's last iterate
      s = ext ? sx_(m,k,j,i) : dep_(m,k,j,i);
      if (e > 0.0) {
        const Real lim = cap*e;
        if (s > lim) s = lim;
        if (s < -lim) s = -lim;
      } else {
        s = 0.0;
      }
      if (!(s == s)) s = 0.0;
      u0(m,IEN,k,j,i) += s;
    }
    sd_(m,k,j,i) = s;
    dep_(m,k,j,i) = s;     // the accumulator starts from what has already been applied
  });
}

int CkImplStep(Mesh *pm, DvceArray5D<Real> u0, DvceArray3D<int> icut_,
               DvceArray4D<Real> T_, DvceArray4D<Real> dx1_, const Real bdt) {
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  auto jac_ = *ck_jac_ptr;
  auto src_ = *ck_src_ptr;
  auto ei_ = *ck_ei_ptr;
  auto est_ = *ck_estar_ptr;
  auto cnv_ = *ck_conv_ptr;
  auto cp_ = *ck_cp_ptr;
  auto dp_ = *ck_dp_ptr;
  auto em_ = *ck_em_ptr;
  auto thk_ = *ck_thu_ptr;
  auto done_ = *ck_done_ptr;
  // problem/ck_impl_warm: the running total this call has applied, carried to the next
  auto dep_ = *ck_dep_ptr;
  const bool skipcol = ck_impl_colskip;
  if (ck_impl_nosync) {
    Kokkos::deep_copy(DevExeSpace(), cnv_, 0.0);   // stream-ordered, no fence
  } else {
    Kokkos::deep_copy(cnv_, 0.0);
  }
  const Real tol = ck_impl_tol;
  const Real dtol = ck_impl_dtol;
  const Real eps = ck_impl_norm_eps;
  const Real dcap = ck_impl_dtmax;
  const Real detot = ck_impl_demax;
  // ---- problem/ck_impl_seed / ck_impl_reuse_jac ----------------------------------
  // seedp_: this pass takes the per-cell initial guess instead of a Newton step.
  // t0rec_: this pass is the one that BUILT the reused Jacobian, so the temperature it
  // was built at is latched here for the rescaling of mode 2.
  const int seedm_ = ck_impl_seed;
  const bool seedp_ = (seedm_ > 0) && (ck_impl_pass == 0);
  const int jscl_ = (ck_impl_reuse_jac == 2) ? 1 : 0;
  const bool t0rec_ = (jscl_ > 0) && (ck_impl_pass == CkImplJacPass());
  auto t0_ = *ck_t0_ptr;

  // ---- problem/ck_impl_fuse: the residual test and the step in ONE team kernel -----
  // One team per column.  Everything per cell -- the residual, the row, the cap, the
  // fallback -- runs in parallel over the cells; the Thomas recurrences run on one lane
  // out of team scratch.  The arithmetic of every per-cell quantity and of the two
  // recurrences is the unfused kernels' in the same order, so the state is bitwise
  // theirs.  The two diagnostic sums (slots 4, 5) are team sums and differ at round-off.
  if (ck_impl_fuse && ck_impl_debug <= 0) {
    const int nk = ke - ks + 1;
    const int nj = je - js + 1;
    const int nkj = nk*nj;
    const int nlg = (nmb1 + 1)*nkj;
    const int ns = ie + 1;                   // scratch rows are indexed by q = i - ic
    // problem/ck_impl_aa: the Gram matrix, right-hand side and coefficients of the
    // least squares go to team scratch too (one element when off)
    const int naa_ = ck_impl_aa;
    const int nsg = (naa_ > 0) ? (CK_AA_MAX*CK_AA_MAX + 2*CK_AA_MAX) : 1;
    const size_t scr = 6*ScrArray1D<Real>::shmem_size(ns)
                       + ScrArray1D<Real>::shmem_size(nsg);
    // problem/ck_impl_cvsec: the secant heat capacity (1-element dummies when off)
    const bool cvs_ = ck_impl_cvsec;
    const bool cv0_ = (ck_impl_pass <= 0);
    const bool cvk_ = cv0_ && ck_impl_cvkeep;
    auto ep_ = cvs_ ? *ck_ep_ptr : CkDum<DvceArray4D<Real>>("ck_ep_d");
    // problem/ck_impl_rsec (1-element dummies when off)
    const Real rsec_ = ck_impl_rsec;
    const bool rsc_ = (rsec_ > 0.0);
    const bool rsu_ = rsc_ && (ck_impl_pass > 0);
    auto rsr_ = rsc_ ? *ck_rsr_ptr : CkDum<DvceArray4D<Real>>("ck_rsr_d");
    auto rse_ = rsc_ ? *ck_rse_ptr : CkDum<DvceArray4D<Real>>("ck_rse_d");
    auto tp_ = cvs_ ? *ck_tp_ptr : CkDum<DvceArray4D<Real>>("ck_tp_d");
    auto cv_ = cvs_ ? *ck_cv_ptr : CkDum<DvceArray4D<Real>>("ck_cv_d");
    // problem/ck_impl_glob: line search and sub-steps (1-element dummies when off)
    const bool glb_ = (ck_impl_glob > 0);
    const bool sub_ = (ck_impl_glob == 2);
    auto lsc_ = glb_ ? *ck_lsc_ptr : CkDum<DvceArray4D<Real>>("ck_lsc_d");
    auto lsd_ = glb_ ? *ck_lsd_ptr : CkDum<DvceArray4D<Real>>("ck_lsd_d");
    auto sacc_ = sub_ ? *ck_sacc_ptr : CkDum<DvceArray4D<Real>>("ck_sacc_d");
    const int ntry_ = ck_impl_ls_ntry;
    const Real lsca_ = ck_impl_ls_c;
    const int maxit_ = ck_impl_maxit;
    const int submax_ = ck_impl_sub_max;
    // problem/ck_impl_esc and ck_impl_aa (1-element dummies when off)
    const bool esc_ = (ck_impl_esc > 0) && sub_;
    const Real escr_ = ck_impl_esc_rho;
    const int escx_ = ck_impl_esc_extra;
    const bool pz_ = (ck_impl_pass <= 0);
    // ck_impl_pred (1-element dummy when off)
    const bool pred_ = ck_impl_pred;
    const Real pfac_ = ck_impl_pred_fac;
    const bool evo_ = ck_impl_evalonly;
    // problem/ck_impl_warm_step: no column is accepted on the seeded state itself
    const bool wfs_ = ck_impl_warm && ck_impl_warm_step && pz_;
    auto rprev_ = pred_ ? *ck_rprev_ptr : CkDum<DvceArray3D<Real>>("ck_rprev_d");
    // ck_impl_floorbound / ck_impl_kkt_demax (see their notes)
    const bool fb_ = ck_impl_floorbound;
    const bool kd_ = ck_impl_kkt_demax && (detot > 0.0);
    EOS_Data eosfb_ = (pm->pmb_pack->pmhd != nullptr)
                      ? pm->pmb_pack->pmhd->peos->eos_data
                      : pm->pmb_pack->phydro->peos->eos_data;
    const bool aarst_ = ck_impl_aa_rst;
    // ck_impl_stalldbg / ck_impl_kkt_row (see their notes)
    const bool sdbg_ = (ck_impl_stalldbg > 0) && (ck_impl_pass >= ck_impl_stalldbg)
                       && !glb_;
    const bool krow_ = ck_impl_kkt_row && (fb_ || kd_) && !glb_;
    // ck_impl_osc (1-element dummy when off)
    const bool osc_ = (ck_impl_osc > 0) && !glb_ && (naa_ == 0);
    const bool oscp_ = osc_ && (ck_impl_pass >= ck_impl_osc);
    auto osa_ = osc_ ? *ck_osc_ptr : CkDum<DvceArray4D<Real>>("ck_osc_d");
    const int pass_ = ck_impl_pass;
    const bool csph_ = pm->use_cubed_sphere;
    // ck_impl_negdbg (1-element counter when off)
    const int ngd_ = ck_impl_negdbg;
    const bool rfb_ = ck_impl_rowfb;
    const bool npv_ = ck_impl_negpiv;
    const bool fstp_ = ck_force_step;       // ck2nd-1102 ck_impl_xpass
    const Real tinc_ = ck_tol_inc;          // ck2nd-1102 (0 = off)
    const Real tflr_ = ck_tol_floor;
    const Real txfl_ = ck_tol_xfloor;
    const Real ptp_ = (ck_pbar_ptr != nullptr) ? ck_impl_tol_ptop : 0.0;
    const Real ept_ = ck_impl_norm_eps_top;
    auto pbr_ = (ptp_ > 0.0) ? *ck_pbar_ptr : CkDum<DvceArray4D<Real>>("ck_pbar_d");
    DvceArray1D<int> ngc_("ck_negdbg_cnt", 1);
    auto &mbpan_ = pm->pmb_pack->pmb->mb_panel;
    auto &x2v_ = pm->pmb_pack->pcoord->x2v;
    auto &x3v_ = pm->pmb_pack->pcoord->x3v;
    // ck_impl_every_thr: a column masked out of this call (ck_done = 2) is not touched
    const bool msk_ = ck_cad_partial;
    auto aah_ = (naa_ > 0) ? *ck_aah_ptr : CkDum<DvceArray5D<Real>>("ck_aah_d");
    auto aac_ = (naa_ > 0) ? *ck_aac_ptr : CkDum<DvceArray4D<Real>>("ck_aac_d");
    // one wavefront per column on a device; the host backends take their own size
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
    Kokkos::TeamPolicy<> tpol(DevExeSpace(), nlg, 64);
#else
    Kokkos::TeamPolicy<> tpol(DevExeSpace(), nlg, Kokkos::AUTO);
#endif
    // level-0 team scratch unless the device cannot fit it (6 ns reals grow with nx1)
    const int cklv_ = TeamScratchLevel(scr, 0, tpol.team_size());
    tpol.set_scratch_size(cklv_, Kokkos::PerTeam(scr));
    auto body = KOKKOS_LAMBDA(TeamMember_t tm) {
      const int lr = tm.league_rank();
      const int m = lr/nkj;
      const int k = (lr - m*nkj)/nj + ks;
      const int j = (lr - m*nkj) % nj + js;
      if (msk_ && done_(m,k,j) > 1.5) return;
      const int ic = icut_(m,k,j);
      // ---- the residual of the column (ck_impl_res) ----
      Real emax = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
      [&](const int i, Real &mx) {
        const Real e = ei_(m,k,j,i);
        if (e > mx) mx = e;
      }, Kokkos::Max<Real>(emax));
      if (!(emax > 0.0)) {
        Kokkos::single(Kokkos::PerTeam(tm), [&]() { done_(m,k,j) = 1.0; });
        return;
      }
      Real rn = 0.0, sg = 0.0, ss = 0.0;
      // this column's implicit interval: bdt, or bdt/2^L on sub-step level L (ls_sub)
      Real hh = bdt;
      // ck_impl_glob: the merit of this iterate and the sub-step bookkeeping to store
      Real phi = 0.0;
      int sst = 0, cnt = 0, lvl = 0;
      // ck_impl_aa: the largest e* of the column (the merit weights)
      Real esmax = 0.0;
      // ck_impl_pred: this column's step is predicted to converge it
      bool pwill = false;
      // ck_impl_pred: a column marked converged by the prediction on an earlier pass
      // of this call (its arrays hold the state BEFORE that last step, so its residual
      // here is stale and must not steer anything)
      const bool pdone = pred_ && (done_(m,k,j) > 0.0);
      int iwd = -1;                     // ck_impl_stalldbg: the worst cell
      if (!glb_) {
        // KKT (ck_impl_floorbound / ck_impl_kkt_demax): a cell on a bound whose
        // residual points through it is at its constrained solution (see the notes)
        auto kktcell = [&](const int i, const Real r) -> bool {
          const Real e = ei_(m,k,j,i);
          bool kkt = false;
          if (fb_ && r > 0.0) {
            const Real rho = u0(m,IDN,k,j,i);
            Real efl = eosfb_.EnergyFromTemperature(rho, eosfb_.tfloor);
            if (e < eosfb_.EnergyFloorBound(rho)) {
              const Real ep = eosfb_.EnergyFromPressure(rho, eosfb_.pfloor);
              if (ep > efl) efl = ep;
            }
            kkt = (e <= efl*(1.0 + 1.0e-10));
          }
          const Real es = est_(m,k,j,i);
          if (kd_ && !kkt && es > 0.0) {
            kkt = (r > 0.0 && e <= es*(1.0 - detot)*(1.0 + 1.0e-10))
                  || (r < 0.0 && e >= es*(1.0 + detot)*(1.0 - 1.0e-10));
          }
          return kkt;
        };
        const bool kk_ = fb_ || kd_;
        // ck_impl_tol_ptop: the norm's share eps of a top cell
        auto epsc = [&](const int i) -> Real {
          return (ptp_ > 0.0 && pbr_(m,k,j,i) < ptp_) ? ept_ : eps;
        };
        // ck2nd-1102 problem/ck_impl_tol_inc: the norm relative to the step (s <= tol
        // <=> |R| <= tau |e - e*| + floor |e| + xfloor bdt (|S| + 2 E): the cell's OWN
        // e and its gross exchange A + E, the round-off scale of S = A - E; no eps e_max
        // share, which in a thin top cell is 1e5-1e9 x its e); off: today's norm
        auto snrm = [&](const int i, const Real r) -> Real {
          // (ck_impl_tol_ptop > 0: the top cells p < ptop keep today's test)
          if (tinc_ > 0.0 && !(ptp_ > 0.0 && pbr_(m,k,j,i) < ptp_)) {
            const Real ec = ei_(m,k,j,i);
            const Real den = tinc_*fabs(ec - est_(m,k,j,i)) + tflr_*fabs(ec)
                             + txfl_*bdt*(fabs(src_(m,k,j,i)) + 2.0*fabs(em_(m,k,j,i)));
            return (den > 0.0) ? tol*fabs(r)/den : ((r == 0.0) ? 0.0 : 1.0e300);
          }
          return fabs(r)/(ei_(m,k,j,i) + epsc(i)*emax);
        };
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &mx) {
          const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
          Real s = snrm(i, r);
          if (kk_ && s > tol && kktcell(i, r)) s = 0.0;
          if (s > mx) mx = s;
          if (t0rec_) t0_(m,k,j,i) = T_(m,k,j,i);
        }, Kokkos::Max<Real>(rn));
        if (!(rn > 0.0)) rn = 0.0;             // the unfused max starts from 0
        if (sdbg_) {
          using MLc = Kokkos::MaxLoc<Real, int>;
          typename MLc::value_type vw;
          Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
          [&](const int i, typename MLc::value_type &u) {
            const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
            Real s = snrm(i, r);
            if (kk_ && s > tol && kktcell(i, r)) s = 0.0;
            if (s > u.val) {
              u.val = s;
              u.loc = i;
            }
          }, MLc(vw));
          iwd = vw.loc;
        }
        if (kk_) {
          Real nkk = 0.0;
          Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
          [&](const int i, Real &sm) {
            const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
            const Real s = snrm(i, r);
            if (s > tol && kktcell(i, r)) sm += 1.0;
          }, nkk);
          if (nkk > 0.0) {
            Kokkos::single(Kokkos::PerTeam(tm), [&]() {
              Kokkos::atomic_add(&cnv_(17), nkk);
            });
          }
        }
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &sm) {
          sm += (ei_(m,k,j,i) - est_(m,k,j,i))*dx1_(m,k,j,i);
        }, sg);
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &sm) {
          sm += bdt*src_(m,k,j,i)*dx1_(m,k,j,i);
        }, ss);
        // ck_impl_warm_step: pass 0 of a warm-started call always takes a step
        const bool cvd = (rn <= tol) && !wfs_ && !fstp_;
        if (!pred_) {
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            Kokkos::atomic_max(&cnv_(0), rn);
            Kokkos::atomic_add(&cnv_(4), sg);
            Kokkos::atomic_add(&cnv_(5), ss);
            if (cvd) {
              done_(m,k,j) = 1.0;
            } else {
              Kokkos::atomic_add(&cnv_(6), 1.0);
            }
          });
          if (cvd) return;
        } else {
          // ck_impl_pred: the same bookkeeping, minus the stale residual of a column
          // already predicted converged; the contraction estimate for the rest
          const Real rp = pz_ ? 0.0 : rprev_(m,k,j);
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            if (!pdone) Kokkos::atomic_max(&cnv_(0), rn);
            Kokkos::atomic_add(&cnv_(4), sg);
            Kokkos::atomic_add(&cnv_(5), ss);
            if (cvd) {
              done_(m,k,j) = 1.0;
            } else if (!pdone) {
              Kokkos::atomic_add(&cnv_(6), 1.0);
              rprev_(m,k,j) = rn;
            }
          });
          if (cvd || pdone || evo_) return;
          pwill = (rp > 0.0) && (rn*rn <= pfac_*tol*rp);
        }
      } else {
        // ---- problem/ck_impl_glob: evaluate the pending trial, then decide ----
        const Real phia = lsc_(m,0,k,j);
        const Real alp = lsc_(m,1,k,j);
        const int ntr = static_cast<int>(lsc_(m,2,k,j));
        lvl = static_cast<int>(lsc_(m,3,k,j));
        sst = static_cast<int>(lsc_(m,4,k,j));
        cnt = static_cast<int>(lsc_(m,5,k,j)) + 1;
        hh = bdt/static_cast<Real>(1 << lvl);
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &mx) {
          const Real e = est_(m,k,j,i);
          if (e > mx) mx = e;
        }, Kokkos::Max<Real>(esmax));
        if (!(esmax > 0.0)) esmax = emax;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &mx) {
          const Real sac = sub_ ? sacc_(m,k,j,i) : 0.0;
          const Real r = ((ei_(m,k,j,i) - est_(m,k,j,i)) - sac) - hh*src_(m,k,j,i);
          const Real s = fabs(r)/(ei_(m,k,j,i) + eps*emax);
          if (s > mx) mx = s;
          if (t0rec_) t0_(m,k,j,i) = T_(m,k,j,i);
        }, Kokkos::Max<Real>(rn));
        if (!(rn > 0.0)) rn = 0.0;
        Real ph2 = 0.0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &sm) {
          const Real sac = sub_ ? sacc_(m,k,j,i) : 0.0;
          const Real r = ((ei_(m,k,j,i) - est_(m,k,j,i)) - sac) - hh*src_(m,k,j,i);
          const Real es = est_(m,k,j,i);
          const Real s = r/(((es > 0.0) ? es : 0.0) + eps*esmax);
          sm += s*s;
        }, ph2);
        phi = sqrt(ph2);
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &sm) {
          sm += (ei_(m,k,j,i) - est_(m,k,j,i))*dx1_(m,k,j,i);
        }, sg);
        // the source the gas has been handed: the finished sub-steps plus this one's
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &sm) {
          const Real sac = sub_ ? sacc_(m,k,j,i) : 0.0;
          sm += (sac + hh*src_(m,k,j,i))*dx1_(m,k,j,i);
        }, ss);
        const bool isdone = (done_(m,k,j) > 0.0);
        tm.team_barrier();
        Kokkos::single(Kokkos::PerTeam(tm), [&]() {
          Kokkos::atomic_max(&cnv_(0), rn);
          Kokkos::atomic_add(&cnv_(4), sg);
          Kokkos::atomic_add(&cnv_(5), ss);
        });
        if (isdone) return;               // converged, or failed at the finest level
        const bool cvd = (rn <= tol);
        // (a) the pending trial: not enough decrease -> halve it, no step this pass
        if (alp > 0.0 && !cvd && ntr < ntry_ && !(phi <= (1.0 - lsca_*alp)*phia)) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
            const Real h = 0.5*lsd_(m,k,j,i);
            u0(m,IEN,k,j,i) -= h;
            dep_(m,k,j,i) -= h;
            lsd_(m,k,j,i) = h;
          });
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            lsc_(m,1,k,j) = 0.5*alp;
            lsc_(m,2,k,j) = static_cast<Real>(ntr + 1);
            lsc_(m,5,k,j) = static_cast<Real>(cnt);
            Kokkos::atomic_add(&cnv_(6), 1.0);
            Kokkos::atomic_add(&cnv_(8), 1.0);
          });
          return;
        }
        bool adv = false;
        if (cvd) {
          if (sub_ && sst < (1 << lvl) - 1) {
            adv = true;                   // this sub-step is done; start the next one
          } else {
            Kokkos::single(Kokkos::PerTeam(tm), [&]() { done_(m,k,j) = 1.0; });
            return;
          }
        }
        // ck_impl_esc: a (sub-)step past its maxit evaluations that is still contracting
        // fast enough to reach tol within ck_impl_esc_extra more passes keeps going
        bool keep = false;
        if (esc_ && !cvd && cnt > maxit_ && cnt <= maxit_ + escx_ && phia > 0.0
            && rn > tol && tol > 0.0) {
          const Real rt = phi/phia;
          if (rt > 0.0 && rt <= escr_) {
            const Real np = log(tol/rn)/log(rt);
            if (np <= static_cast<Real>(maxit_ + escx_ + 1 - cnt)) keep = true;
          }
        }
        if (keep) {
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            Kokkos::atomic_add(&cnv_(12), 1.0);
          });
        } else if (!cvd && sub_ && cnt > maxit_) {
          // (b) the sub-step did not converge in maxit passes: redo the column from e^n
          // with twice the sub-steps, or give up at the finest level
          if ((2 << lvl) <= submax_ && esc_) {
            // ck_impl_esc: redo ONLY this sub-step, from its own start e^n + sacc, as two
            // sub-steps of half its length
            Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
              u0(m,IEN,k,j,i) -= dep_(m,k,j,i) - sacc_(m,k,j,i);
              dep_(m,k,j,i) = sacc_(m,k,j,i);
            });
            Kokkos::single(Kokkos::PerTeam(tm), [&]() {
              lsc_(m,0,k,j) = 0.0;
              lsc_(m,1,k,j) = 0.0;
              lsc_(m,2,k,j) = 0.0;
              lsc_(m,3,k,j) = static_cast<Real>(lvl + 1);
              lsc_(m,4,k,j) = static_cast<Real>(2*sst);
              lsc_(m,5,k,j) = 0.0;
              Kokkos::atomic_add(&cnv_(6), 1.0);
              Kokkos::atomic_add(&cnv_(9), 1.0);
            });
          } else if ((2 << lvl) <= submax_) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
              u0(m,IEN,k,j,i) -= dep_(m,k,j,i);
              dep_(m,k,j,i) = 0.0;
              sacc_(m,k,j,i) = 0.0;
            });
            Kokkos::single(Kokkos::PerTeam(tm), [&]() {
              lsc_(m,0,k,j) = 0.0;
              lsc_(m,1,k,j) = 0.0;
              lsc_(m,2,k,j) = 0.0;
              lsc_(m,3,k,j) = static_cast<Real>(lvl + 1);
              lsc_(m,4,k,j) = 0.0;
              lsc_(m,5,k,j) = 0.0;
              Kokkos::atomic_add(&cnv_(6), 1.0);
              Kokkos::atomic_add(&cnv_(9), 1.0);
            });
          } else {
            Kokkos::single(Kokkos::PerTeam(tm), [&]() {
              done_(m,k,j) = 1.0;
              Kokkos::atomic_add(&cnv_(10), 1.0);
            });
          }
          return;
        }
        if (adv) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
            sacc_(m,k,j,i) += hh*src_(m,k,j,i);
          });
          tm.team_barrier();
          sst += 1;
          // ck_impl_esc: a sub-step that converged quickly hands the column back to the
          // coarser level wherever the finished sub-steps end on its grid
          if (esc_ && lvl > 0 && (sst % 2) == 0 && 2*(cnt - 1) <= maxit_) {
            lvl -= 1;
            sst /= 2;
            hh = bdt/static_cast<Real>(1 << lvl);
            Kokkos::single(Kokkos::PerTeam(tm), [&]() {
              Kokkos::atomic_add(&cnv_(13), 1.0);
            });
          }
          cnt = 1;
          ph2 = 0.0;
          Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
          [&](const int i, Real &sm) {
            const Real r = ((ei_(m,k,j,i) - est_(m,k,j,i)) - sacc_(m,k,j,i))
                           - hh*src_(m,k,j,i);
            const Real es = est_(m,k,j,i);
            const Real s = r/(((es > 0.0) ? es : 0.0) + eps*esmax);
            sm += s*s;
          }, ph2);
          phi = sqrt(ph2);
        }
        Kokkos::single(Kokkos::PerTeam(tm), [&]() { Kokkos::atomic_add(&cnv_(6), 1.0); });
      }
      // ls_sub: the first pass of every sub-step (after a restart from e^n, or after
      // the previous sub-step converged) takes the ck_impl_seed = 2 guess too
      const bool sdc = seedp_ || (sub_ && seedm_ == 2 && cnt == 1);
      // ---- the step (ck_impl_tri) ----
      if (ic > ie) return;
      if (done_(m,k,j) > 0.0) return;
      const int n = ie - ic + 1;
      ScrArray1D<Real> sa(tm.team_scratch(cklv_), ns);
      ScrArray1D<Real> sb(tm.team_scratch(cklv_), ns);
      ScrArray1D<Real> sc(tm.team_scratch(cklv_), ns);
      ScrArray1D<Real> sd(tm.team_scratch(cklv_), ns);
      ScrArray1D<Real> sx(tm.team_scratch(cklv_), ns);
      ScrArray1D<Real> sv(tm.team_scratch(cklv_), ns);
      ScrArray1D<Real> sgm(tm.team_scratch(cklv_), nsg);
      constexpr int g0 = CK_AA_MAX*CK_AA_MAX + CK_AA_MAX;   // gamma's slots in sg
      // ck_impl_cvsec: this pass's cv of every cell, before any row reads a neighbour's
      if (cvs_) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, 0, n), [&](const int q) {
          const int i = ic + q;
          const Real ei = ei_(m,k,j,i);
          const Real Ti = T_(m,k,j,i);
          Real cv = (ei > 0.0 && Ti > 0.0) ? ei/Ti : 1.0;
          // ck_impl_cvkeep: pass 0 starts from the secant cv the previous call ended on
          if (cvk_) {
            const Real cvo = cv_(m,k,j,i);
            if (cvo > 0.02*cv && cvo < 50.0*cv) cv = cvo;
          }
          if (!cv0_) {
            const Real cvo = cv_(m,k,j,i);
            const Real dT = Ti - tp_(m,k,j,i);
            const Real de = ei - ep_(m,k,j,i);
            Real cs = (cvo > 0.0) ? cvo : cv;
            if (fabs(dT) > 1.0e-7*Ti) {
              const Real r = de/dT;
              if (r > 0.02*cv && r < 50.0*cv) cs = r;
            }
            cv = cs;
          }
          sv(q) = cv;
        });
        tm.team_barrier();
      }
      // the rows, in parallel; a row that the unfused sweep would call bad sets the flag.
      // One int carries both counts: bad rows + 65536 x thin rows (n < 65536).
      int nbt = 0;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, 0, n),
      [&](const int q, int &nb) {
        const int i = ic + q;
        const Real ei = ei_(m,k,j,i);
        const Real Ti = T_(m,k,j,i);
        if (!(ei > 0.0) || !(Ti > 0.0)) {
          nb += 1;
          return;
        }
        // a DEAD cell (see CkDeadE): identity row, no step
        if (!(ei > CkDeadE)) {
          sa(q) = 0.0;
          sb(q) = 1.0;
          sc(q) = 0.0;
          sd(q) = 0.0;
          return;
        }
        Real cvi = ei/Ti;
        const Real eim = (q > 0) ? ei_(m,k,j,i-1) : 0.0;
        const Real Tim = (q > 0) ? T_(m,k,j,i-1) : 1.0;
        const Real eip = (q < n-1) ? ei_(m,k,j,i+1) : 0.0;
        const Real Tip = (q < n-1) ? T_(m,k,j,i+1) : 1.0;
        Real cvm = (q > 0 && Tim > 0.0 && eim > 0.0) ? (eim/Tim) : 1.0;
        Real cvp = (q < n-1 && Tip > 0.0 && eip > 0.0) ? (eip/Tip) : 1.0;
        if (cvs_) {
          cvi = sv(q);
          if (q > 0) cvm = sv(q-1);
          if (q < n-1) cvp = sv(q+1);
        }
        const bool thin = !(thk_(m,k,j,i) > 0.0);
        Real a, b, c, d;
        if (sdc) {
          const Real emc = em_(m,k,j,i);
          const Real sc0 = src_(m,k,j,i);
          a = 0.0;
          b = 1.0;
          c = 0.0;
          if (sub_ && seedm_ == 2) {
            d = CkThinSolve(ei, est_(m,k,j,i) + sacc_(m,k,j,i), sc0, emc, hh) - ei;
          } else if (seedm_ == 2) {
            d = CkThinSolve(ei, est_(m,k,j,i), sc0, emc, bdt) - ei;
          } else {
            const Real de0 = est_(m,k,j,i) - ei;
            if (emc > 0.0) {
              const Real lam = 4.0*emc/ei;
              const Real x = lam*bdt;
              d = de0 + ((x > 1.0e-4) ? (sc0/lam)*(-expm1(-x)) : sc0*bdt);
            } else {
              d = de0 + sc0*bdt;
            }
          }
        } else if (thin) {
          a = 0.0;
          b = 1.0;
          c = 0.0;
          if (sub_) {
            d = CkThinSolve(ei, est_(m,k,j,i) + sacc_(m,k,j,i), src_(m,k,j,i),
                            em_(m,k,j,i), hh) - ei;
          } else {
            d = CkThinSolve(ei, est_(m,k,j,i), src_(m,k,j,i), em_(m,k,j,i), hh) - ei;
          }
          nb += 65536;
        } else {
          Real s0 = 1.0, s1 = 1.0, s2 = 1.0;
          if (jscl_ > 0) {
            const Real tb = t0_(m,k,j,i);
            if (Ti > 0.0 && tb > 0.0) {
              const Real r = Ti/tb;
              s1 = r*r*r;
            }
            if (q > 0) {
              const Real tbm = t0_(m,k,j,i-1);
              if (Tim > 0.0 && tbm > 0.0) {
                const Real r = Tim/tbm;
                s0 = r*r*r;
              }
            }
            if (q < n-1) {
              const Real tbp = t0_(m,k,j,i+1);
              if (Tip > 0.0 && tbp > 0.0) {
                const Real r = Tip/tbp;
                s2 = r*r*r;
              }
            }
          }
          a = (q > 0) ? (-hh*s0*jac_(m,0,k,j,i)/cvm) : 0.0;
          b = 1.0 - hh*s1*jac_(m,1,k,j,i)/cvi;
          c = (q < n-1) ? (-hh*s2*jac_(m,2,k,j,i)/cvp) : 0.0;
          if (!(eim > CkDeadE)) a = 0.0;          // no coupling to a dead neighbour
          if (!(eip > CkDeadE)) c = 0.0;
          if (sub_) {
            d = -(((ei - est_(m,k,j,i)) - sacc_(m,k,j,i)) - hh*src_(m,k,j,i));
          } else {
            d = -(ei - est_(m,k,j,i) - hh*src_(m,k,j,i));
          }
          // problem/ck_impl_rsec: the secant bound on the diagonal (see its note)
          if (rsc_ && !sub_) {
            if (rsu_ && b > 0.0) {
              const Real dei = ei - rse_(m,k,j,i);
              if (fabs(dei) > 1.0e-9*ei) {
                const Real js = (-d - rsr_(m,k,j,i))/dei;
                if (js > b) b = fmin(js, rsec_*b);
              }
            }
            rsr_(m,k,j,i) = -d;
            rse_(m,k,j,i) = ei;
          }
        }
        // problem/ck_impl_kkt_row: a thick KKT cell holds (identity row, no step); after
        // rsec, so the secant history keeps the cell's true residual
        if (krow_ && !sdc && !thin && CkKktCell(eosfb_, fb_, kd_, detot,
                                                u0(m,IDN,k,j,i), ei, est_(m,k,j,i), -d)) {
          a = 0.0;
          b = 1.0;
          c = 0.0;
          d = 0.0;
        }
        if (ngd_ > 0 && !(b > 0.0)) {
          if (Kokkos::atomic_fetch_add(&ngc_(0), 1) < ngd_) {
            Kokkos::printf("### ck_negdbg pass=%d m=%d k=%d j=%d i=%d ic=%d thin=%d "
                           "T=%.4e e=%.4e cv=%.4e h=%.4e J=%.4e %.4e %.4e a=%.4e b=%.4e "
                           "c=%.4e S/e=%.4e\n", pass_, m, k, j, i, ic, thin ? 1 : 0, Ti,
                           ei, cvi, hh, jac_(m,0,k,j,i), jac_(m,1,k,j,i), jac_(m,2,k,j,i),
                           a, b, c, hh*src_(m,k,j,i)/ei);
          }
        }
        // ck_impl_rowfb: a row with a non-positive diagonal takes the bounded per-cell
        // step of a thin cell (identity row) instead of sending its column to the
        // fallback
        if (rfb_ && !npv_ && !sdc && !thin && !(b > 0.0)) {
          a = 0.0;
          b = 1.0;
          c = 0.0;
          if (sub_) {
            d = CkThinSolve(ei, est_(m,k,j,i) + sacc_(m,k,j,i), src_(m,k,j,i),
                            em_(m,k,j,i), hh) - ei;
          } else {
            d = CkThinSolve(ei, est_(m,k,j,i), src_(m,k,j,i), em_(m,k,j,i), hh) - ei;
          }
          nb += 65536;
        }
        if (npv_ && !thin && !sdc) {
          // ck_impl_negpiv: the exact Newton row stands; only a vanishing pivot is bad
          if (!(fabs(b) > 1.0e-12)) nb += 1;
        } else if (!(b > 0.0)) {
          nb += 1;
        }
        sa(q) = a;
        sb(q) = b;
        sc(q) = c;
        sd(q) = d;
      }, nbt);
      const int nbad = nbt % 65536;
      const int nthn = nbt/65536;
      tm.team_barrier();
      // ---- problem/ck_impl_aa: the Anderson history and the least squares ----
      bool acc = false;
      int nh = 0;
      if (naa_ > 0 && nbad == 0) {
        const int na = naa_;
        if (!glb_) {
          Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
          [&](const int i, Real &mx) {
            const Real e = est_(m,k,j,i);
            if (e > mx) mx = e;
          }, Kokkos::Max<Real>(esmax));
          if (!(esmax > 0.0)) esmax = emax;
        }
        // a new backward-Euler residual (pass 0, a new sub-step) restarts the history
        const bool rs = sdc || (glb_ ? (cnt == 1) : pz_);
        // the weighted residual norm^2 of this iterate (thick rows)
        Real rr = 0.0;
        Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
        [&](const int i, Real &sm) {
          if (!(thk_(m,k,j,i) > 0.0)) return;
          const Real sac = sub_ ? sacc_(m,k,j,i) : 0.0;
          const Real r = ((ei_(m,k,j,i) - est_(m,k,j,i)) - sac) - hh*src_(m,k,j,i);
          const Real es = est_(m,k,j,i);
          const Real w = 1.0/(((es > 0.0) ? es : 0.0) + eps*esmax);
          sm += w*w*r*r;
        }, rr);
        // restart: a new residual, or the norm grew since the last iterate
        const bool grw = aarst_ && !rs && (aac_(m,2,k,j) > 0.0)
                         && !(rr <= aac_(m,3,k,j));
        nh = (rs || grw) ? 0 : static_cast<int>(aac_(m,0,k,j));
        int hd = (rs || grw) ? 0 : static_cast<int>(aac_(m,1,k,j));
        const bool pv = !rs && !grw && (aac_(m,2,k,j) > 0.0);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
          const Real e = ei_(m,k,j,i);
          const Real sac = sub_ ? sacc_(m,k,j,i) : 0.0;
          const Real r = ((e - est_(m,k,j,i)) - sac) - hh*src_(m,k,j,i);
          if (pv) {
            aah_(m,hd,k,j,i) = e - aah_(m,2*na,k,j,i);
            aah_(m,na+hd,k,j,i) = r - aah_(m,2*na+1,k,j,i);
          }
          aah_(m,2*na,k,j,i) = e;
          aah_(m,2*na+1,k,j,i) = r;
        });
        if (pv) {
          hd = (hd + 1) % na;
          if (nh < na) nh += 1;
        }
        tm.team_barrier();
        Kokkos::single(Kokkos::PerTeam(tm), [&]() {
          aac_(m,0,k,j) = static_cast<Real>(nh);
          aac_(m,1,k,j) = static_cast<Real>(hd);
          aac_(m,2,k,j) = sdc ? 0.0 : 1.0;      // a seed step's pair is never stored
          aac_(m,3,k,j) = rr;
        });
        if (!sdc && nh > 0) {
          // the weighted Gram matrix of dR over the thick rows, and dR^T W^2 r
          for (int a=0; a<nh; ++a) {
            for (int b=0; b<=a+1; ++b) {
              Real g = 0.0;
              Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, ic, ie+1),
              [&](const int i, Real &sm) {
                if (!(thk_(m,k,j,i) > 0.0)) return;
                const Real es = est_(m,k,j,i);
                const Real w = 1.0/(((es > 0.0) ? es : 0.0) + eps*esmax);
                const Real v = (b <= a) ? aah_(m,na+b,k,j,i) : aah_(m,2*na+1,k,j,i);
                sm += w*w*aah_(m,na+a,k,j,i)*v;
              }, g);
              Kokkos::single(Kokkos::PerTeam(tm), [&]() {
                if (b <= a) {
                  sgm(a*CK_AA_MAX + b) = g;
                } else {
                  sgm(CK_AA_MAX*CK_AA_MAX + a) = g;
                }
              });
            }
          }
          tm.team_barrier();
          // regularised Cholesky of the Gram matrix, in place, and gamma
          int ok = 0;
          Kokkos::single(Kokkos::PerTeam(tm), [&](int &okl) {
            Real dmx = 0.0;
            for (int a=0; a<nh; ++a) {
              if (sgm(a*CK_AA_MAX + a) > dmx) dmx = sgm(a*CK_AA_MAX + a);
            }
            okl = (dmx > 0.0) ? 1 : 0;
            const Real lam = 1.0e-10*dmx;
            for (int a=0; a<nh && okl == 1; ++a) {
              for (int b=0; b<=a; ++b) {
                Real v = sgm(a*CK_AA_MAX + b) + ((a == b) ? lam : 0.0);
                for (int c=0; c<b; ++c) v -= sgm(a*CK_AA_MAX + c)*sgm(b*CK_AA_MAX + c);
                if (a == b) {
                  if (!(v > 1.0e-14*dmx)) {
                    okl = 0;
                    break;
                  }
                  sgm(a*CK_AA_MAX + a) = sqrt(v);
                } else {
                  sgm(a*CK_AA_MAX + b) = v/sgm(b*CK_AA_MAX + b);
                }
              }
            }
            if (okl == 1) {
              for (int a=0; a<nh; ++a) {
                Real y = sgm(CK_AA_MAX*CK_AA_MAX + a);
                for (int c=0; c<a; ++c) y -= sgm(a*CK_AA_MAX + c)*sgm(g0 + c);
                sgm(g0 + a) = y/sgm(a*CK_AA_MAX + a);
              }
              for (int a=nh-1; a>=0; --a) {
                Real x = sgm(g0 + a);
                for (int c=a+1; c<nh; ++c) x -= sgm(c*CK_AA_MAX + a)*sgm(g0 + c);
                sgm(g0 + a) = x/sgm(a*CK_AA_MAX + a);
                if (!(x == x)) okl = 0;
              }
            }
          }, ok);
          acc = (ok == 1);
          tm.team_barrier();
        }
        if (acc) {
          // the thick rows solve for the minimised residual r - dR gamma
          Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, 0, n), [&](const int q) {
            const int i = ic + q;
            if (!(thk_(m,k,j,i) > 0.0)) return;
            Real sm = 0.0;
            for (int a=0; a<nh; ++a) sm += sgm(g0 + a)*aah_(m,na+a,k,j,i);
            sd(q) += sm;
          });
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            Kokkos::atomic_add(&cnv_(14), 1.0);
          });
          tm.team_barrier();
        }
      }
      // the Thomas recurrences on one lane: cp overwrites c, dp overwrites d, and the
      // uncapped solution (NaN-guarded, as the unfused back substitution carries it)
      // goes to x
      int bad = (nbad > 0) ? 1 : 0;
      Kokkos::single(Kokkos::PerTeam(tm), [&](int &bd) {
        if (bd == 0) {
          for (int q=0; q<n; ++q) {
            const Real a = sa(q);
            const Real den = sb(q) - a*((q > 0) ? sc(q-1) : 0.0);
            if (!(fabs(den) > 0.0)) {
              bd = 1;
              break;
            }
            const Real cq = sc(q);
            sc(q) = cq/den;
            sd(q) = (sd(q) - a*((q > 0) ? sd(q-1) : 0.0))/den;
          }
        }
        if (bd == 0) {
          Real prev = 0.0;
          for (int q=n-1; q>=0; --q) {
            Real de = sd(q) - sc(q)*prev;
            if (!(de == de)) de = 0.0;
            prev = de;
            sx(q) = de;
          }
        }
        if (nthn > 0) Kokkos::atomic_add(&cnv_(7), static_cast<Real>(nthn));
      }, bad);
      tm.team_barrier();
      if (bad) {
        Kokkos::single(Kokkos::PerTeam(tm), [&]() { Kokkos::atomic_add(&cnv_(3), 1.0); });
        Kokkos::parallel_for(Kokkos::TeamThreadRange(tm, ic, ie+1), [&](const int i) {
          const Real ei = ei_(m,k,j,i);
          if (!(ei > 0.0)) return;
          const Real esb = sub_ ? (est_(m,k,j,i) + sacc_(m,k,j,i)) : est_(m,k,j,i);
          Real de = CkThinSolve(ei, esb, src_(m,k,j,i), em_(m,k,j,i), hh) - ei;
          const Real lim = dcap*ei;
          if (de > lim) de = lim;
          if (de < -lim) de = -lim;
          u0(m,IEN,k,j,i) += de;
          dep_(m,k,j,i) += de;
          if (glb_) lsd_(m,k,j,i) = 0.0;
        });
        if (naa_ > 0) {
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            aac_(m,0,k,j) = 0.0;          // a fallback step restarts the history
            aac_(m,1,k,j) = 0.0;
            aac_(m,2,k,j) = 0.0;
          });
        }
        if (glb_) {
          Kokkos::single(Kokkos::PerTeam(tm), [&]() {
            lsc_(m,0,k,j) = phi;
            lsc_(m,1,k,j) = 0.0;          // a fallback step is not line-searched
            lsc_(m,2,k,j) = 0.0;
            lsc_(m,4,k,j) = static_cast<Real>(sst);
            if (esc_) lsc_(m,3,k,j) = static_cast<Real>(lvl);
            lsc_(m,5,k,j) = static_cast<Real>(cnt);
          });
        }
        return;
      }
      // the caps and the apply, in parallel; the cap flag of row q goes to sa(q)
      Real dn = 0.0;
      int ncp = 0;
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, 0, n),
      [&](const int q, Real &mx) {
        const int i = ic + q;
        Real de = sx(q);
        // ck_impl_aa: from the combination of the past iterates
        if (acc && thk_(m,k,j,i) > 0.0) {
          Real sm = 0.0;
          for (int a=0; a<nh; ++a) sm += sgm(g0 + a)*aah_(m,a,k,j,i);
          de -= sm;
        }
        // ck_impl_osc: Aitken damping of a cell whose step reversed sign (see the note)
        if (oscp_ && thk_(m,k,j,i) > 0.0) {
          const Real sp = osa_(m,k,j,i);
          if (sp*de < 0.0) {
            de = de*(sp/(sp - de));
            Kokkos::atomic_add(&cnv_(18), 1.0);
          }
        }
        bool cap = false;
        const Real ei = ei_(m,k,j,i);
        const Real lim = dcap*ei;
        if (de > lim) {
          de = lim;
          cap = true;
        }
        if (de < -lim) {
          de = -lim;
          cap = true;
        }
        const Real es = est_(m,k,j,i);
        if (detot > 0.0 && es > 0.0) {
          Real en = ei + de;
          const Real ehi = es*(1.0 + detot);
          const Real elo = es*(1.0 - detot);
          if (en > ehi) {
            en = ehi;
            cap = true;
          }
          if (en < elo) {
            en = elo;
            cap = true;
          }
          de = en - ei;
        }
        // ck_impl_floorbound: a cooling step stops at e_floor (rt_floor_consistent rule)
        if (fb_ && de < 0.0) {
          const Real rho = u0(m,IDN,k,j,i);
          Real efl = eosfb_.EnergyFromTemperature(rho, eosfb_.tfloor);
          if (ei + de < eosfb_.EnergyFloorBound(rho)) {
            const Real ep = eosfb_.EnergyFromPressure(rho, eosfb_.pfloor);
            if (ep > efl) efl = ep;
          }
          if (ei + de < efl) {
            de = (ei > efl) ? (efl - ei) : 0.0;
            Kokkos::atomic_add(&cnv_(16), 1.0);
          }
        }
        if (sdbg_) {
          const bool clp = (de != sx(q));
          if (i == iwd || i == iwd - 1 || i == iwd + 1 || clp) {
            const Real rho = u0(m,IDN,k,j,i);
            const Real r = ei - es - bdt*src_(m,k,j,i);
            const Real Ti = T_(m,k,j,i);
            Real efl = eosfb_.EnergyFromTemperature(rho, eosfb_.tfloor);
            const Real ep = eosfb_.EnergyFromPressure(rho, eosfb_.pfloor);
            if (ep > efl) efl = ep;
            const bool kk = CkKktCell(eosfb_, fb_, kd_, detot, rho, ei, es, r);
            Real pth = 0.0, plat = 0.0, plon = 0.0;
            if (csph_) {
              atm_column::CSCellAngles(mbpan_.d_view(m), x2v_(m,j), x3v_(m,k), pth, plat,
                                       plon);
            }
            Kokkos::printf("### ckstall pass=%d m=%d k=%d j=%d i=%d w=%d lat=%.2f "
                           "lon=%.2f "
                           "T=%.5e rho=%.4e p=%.4e r/e=%.4e e=%.10e est=%.6e hS/e=%.4e "
                           "hEm/e=%.4e Etot/e=%.4e efl/e=%.6e cvT=%.4e cvEOS=%.4e "
                           "thick=%d kkt=%d b=%.4e dx/e=%.4e de/e=%.4e cap=%d\n",
                           pass_, m, k, j, i, (i == iwd) ? 1 : 0,
                           plat*57.29577951308232, plon*57.29577951308232, Ti, rho,
                           eosfb_.Pressure(rho, ei, Ti), r/ei, ei, es,
                           bdt*src_(m,k,j,i)/ei, bdt*em_(m,k,j,i)/ei,
                           u0(m,IEN,k,j,i)/ei, efl/ei, (Ti > 0.0) ? ei/Ti : 0.0,
                           eosfb_.SpecificHeatCv(rho, ei, Ti)*rho,
                           (thk_(m,k,j,i) > 0.0) ? 1 : 0, kk ? 1 : 0, sb(q),
                           sx(q)/ei, de/ei, cap ? 1 : 0);
          }
        }
        u0(m,IEN,k,j,i) += de;
        dep_(m,k,j,i) += de;
        if (osc_) osa_(m,k,j,i) = de;
        if (glb_) lsd_(m,k,j,i) = de;       // the trial the next pass evaluates
        const Real s = (ei > 0.0) ? fabs(de)/ei : 0.0;
        if (s > mx) mx = s;
        sa(q) = cap ? 1.0 : 0.0;
        if (cvs_) {
          ep_(m,k,j,i) = ei;
          tp_(m,k,j,i) = T_(m,k,j,i);
          cv_(m,k,j,i) = sv(q);
        }
      }, Kokkos::Max<Real>(dn));
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, 0, n),
      [&](const int q, int &nc) {
        if (sa(q) > 0.0) nc += 1;
      }, ncp);
      if (!(dn > 0.0)) dn = 0.0;
      Kokkos::single(Kokkos::PerTeam(tm), [&]() {
        Kokkos::atomic_max(&cnv_(1), dn);
        if (ncp > 0) Kokkos::atomic_add(&cnv_(2), 1.0);
        // ck_impl_pred: an uncapped Newton step predicted to land below tol ends the
        // column here, without a confirmation sweep
        if (pwill && ncp == 0) {
          done_(m,k,j) = 1.0;
          Kokkos::atomic_add(&cnv_(15), 1.0);
        }
        if (glb_) {
          lsc_(m,0,k,j) = phi;
          lsc_(m,1,k,j) = sdc ? 0.0 : 1.0;      // a seed step is not line-searched
          // a sub-step seed: the next pass rebuilds the chord Jacobian at its result
          if (sub_ && sdc && !seedp_) Kokkos::atomic_add(&cnv_(11), 1.0);
          lsc_(m,2,k,j) = 0.0;
          lsc_(m,4,k,j) = static_cast<Real>(sst);
          if (esc_) lsc_(m,3,k,j) = static_cast<Real>(lvl);
          lsc_(m,5,k,j) = static_cast<Real>(cnt);
        }
      });
    };
    if (ck_impl_lw) {
      Kokkos::parallel_for("ck_impl_fused", Kokkos::Experimental::require(tpol,
                           Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                           body);
    } else {
      Kokkos::parallel_for("ck_impl_fused", tpol, body);
    }
    // ck_impl_nosync: one cached host mirror instead of an allocation per pass
    static DvceArray1D<Real>::HostMirror *hcf_c = nullptr;
    if (ck_impl_nosync && hcf_c == nullptr) {
      hcf_c = new DvceArray1D<Real>::HostMirror(Kokkos::create_mirror_view(cnv_));
    }
    auto hcf = ck_impl_nosync ? *hcf_c : Kokkos::create_mirror_view(cnv_);
    Kokkos::deep_copy(hcf, cnv_);
    ck_impl_last_res = hcf(0);
    ck_impl_last_gap = (hcf(5) != 0.0) ? (hcf(4)/hcf(5) - 1.0) : 0.0;
    ck_impl_nactive = static_cast<int>(hcf(6));
    if (ck_impl_floorbound || ck_impl_kkt_demax) {
      ck_impl_nfloor = static_cast<int>(hcf(16));
      ck_impl_nkkt = static_cast<int>(hcf(17));
    }
    if (ck_impl_osc > 0) {
      ck_impl_nosc = static_cast<int>(hcf(18));
    }
    if (ck_impl_glob > 0) {
      ck_impl_nrej += static_cast<int>(hcf(8));
      ck_impl_nsub += static_cast<int>(hcf(9));
      ck_impl_nsubfail += static_cast<int>(hcf(10));
      ck_impl_jac_again = (ck_impl_glob == 2) && (hcf(11) > 0.0);
    }
    if (ck_impl_esc > 0) {
      ck_impl_nesc += static_cast<int>(hcf(12));
      ck_impl_ncoarse += static_cast<int>(hcf(13));
    }
    if (ck_impl_aa > 0) ck_impl_naa += static_cast<int>(hcf(14));
    if (hcf(6) == 0.0) return 0;
    ck_impl_last_dst = hcf(1);
    ck_impl_ncap = static_cast<int>(hcf(2));
    ck_impl_nfall = static_cast<int>(hcf(3));
    ck_impl_nthin = static_cast<int>(hcf(7));
    if (hcf(1) <= dtol && hcf(0) <= tol) return 0;
    // ck_impl_pred: every column still active has been predicted converged by its step
    ck_impl_npred = static_cast<int>(hcf(15));
    if (ck_impl_pred && hcf(15) >= hcf(6)) return 0;
    return 1;
  }

  // ---- pass 1: the residual norm of the CURRENT iterate --------------------------
  par_for("ck_impl_res", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    Real emax = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real e = ei_(m,k,j,i);
      if (e > emax) emax = e;
    }
    if (!(emax > 0.0)) {
      done_(m,k,j) = 1.0;      // nothing usable here; it is not part of the system
      return;
    }
    Real rn = 0.0;
    Real sg = 0.0, ss = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
      const Real s = fabs(r)/(ei_(m,k,j,i) + eps*emax);
      if (s > rn) rn = s;
      const Real dx = dx1_(m,k,j,i);
      sg += (ei_(m,k,j,i) - est_(m,k,j,i))*dx;
      ss += bdt*src_(m,k,j,i)*dx;
      // ck_impl_reuse_jac = 2: latch the T the Jacobian of this pass was built at
      if (t0rec_) t0_(m,k,j,i) = T_(m,k,j,i);
    }
    Kokkos::atomic_max(&cnv_(0), rn);
    Kokkos::atomic_add(&cnv_(4), sg);
    Kokkos::atomic_add(&cnv_(5), ss);
    // ---- ck_impl_colskip: this COLUMN's own verdict.  Once set it is never cleared:
    // the sweep no longer touches the column, so neither its source nor its energy can
    // change and the residual it was retired on is the one it keeps.
    if (rn <= tol) {
      done_(m,k,j) = 1.0;
    } else {
      Kokkos::atomic_add(&cnv_(6), 1.0);
    }
  });
  auto hc = Kokkos::create_mirror_view(cnv_);
  Kokkos::deep_copy(hc, cnv_);
  ck_impl_last_res = hc(0);
  ck_impl_last_gap = (hc(5) != 0.0) ? (hc(4)/hc(5) - 1.0) : 0.0;
  ck_impl_nactive = static_cast<int>(hc(6));
  // ---- problem/ck_impl_debug: print the WORST cell of the pass, with everything the
  // row is built from.  Diagnosis only; nothing reads it.
  if (ck_impl_debug > 0) {
    const bool dbg2_ = (ck_impl_debug >= 2);
    const Real rmax = hc(0);
    const Real emref = ck_impl_norm_eps;
    par_for("ck_impl_dbg", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      const int ic = icut_(m,k,j);
      Real emax = 0.0;
      for (int i=ic; i<ie+1; ++i) {
        if (ei_(m,k,j,i) > emax) emax = ei_(m,k,j,i);
      }
      if (!(emax > 0.0)) return;
      for (int i=ic; i<ie+1; ++i) {
        const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
        const Real s = fabs(r)/(ei_(m,k,j,i) + emref*emax);
        if (s >= 0.999*rmax) {
          const Real emc = em_(m,k,j,i);
          Kokkos::printf("### ckdbg m=%d k=%d j=%d i=%d r=%.4e e=%.5e est=%.5e "
                         "src=%.5e em=%.5e AoverE=%.4e thick=%d T=%.4e bdt=%.3e\n",
                         m, k, j, i, s, ei_(m,k,j,i), est_(m,k,j,i), src_(m,k,j,i),
                         emc, (emc > 0.0) ? (src_(m,k,j,i) + emc)/emc : -1.0,
                         (thk_(m,k,j,i) > 0.0) ? 1 : 0, T_(m,k,j,i), bdt);
          // ck_impl_debug >= 2: every cell of this column the call has moved by more
          // than 20 % so far, or whose Jacobian diagonal is not negative, with its row
          if (dbg2_) {
            for (int i2=ic; i2<ie+1; ++i2) {
              const Real e2 = ei_(m,k,j,i2), d2 = dep_(m,k,j,i2);
              if ((e2 > 0.0 && fabs(d2) > 0.2*e2) || !(jac_(m,1,k,j,i2) < 0.0)) {
                Kokkos::printf("### ckcap i=%d e=%.5e dep=%.5e est=%.5e T=%.4e "
                               "src=%.5e em=%.5e thick=%d jac=%.4e %.4e %.4e\n",
                               i2, e2, d2, est_(m,k,j,i2), T_(m,k,j,i2),
                               src_(m,k,j,i2), em_(m,k,j,i2),
                               (thk_(m,k,j,i2) > 0.0) ? 1 : 0, jac_(m,0,k,j,i2),
                               jac_(m,1,k,j,i2), jac_(m,2,k,j,i2));
              }
            }
          }
          return;
        }
      }
    });
  }
  if (hc(6) == 0.0) return 0;

  // ---- pass 2: assemble, solve, apply --------------------------------------------
  par_for("ck_impl_tri", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    if (done_(m,k,j) > 0.0) return;      // converged; see ck_impl_colskip
    // the Thomas sweep.  One thread owns the column, exactly as the sweeps do; the two
    // work rows live in device arrays rather than on the stack so the kernel carries no
    // per-thread footprint that would keep it off an accelerator.
    const int n = ie - ic + 1;
    if (n < 1) return;
    bool bad = false;
    for (int q=0; q<n; ++q) {
      const int i = ic + q;
      const Real ei = ei_(m,k,j,i);
      const Real Ti = T_(m,k,j,i);
      if (!(ei > 0.0) || !(Ti > 0.0)) {
        bad = true;
        break;
      }
      const Real cvi = ei/Ti;
      const Real eim = (q > 0) ? ei_(m,k,j,i-1) : 0.0;
      const Real Tim = (q > 0) ? T_(m,k,j,i-1) : 1.0;
      const Real eip = (q < n-1) ? ei_(m,k,j,i+1) : 0.0;
      const Real Tip = (q < n-1) ? T_(m,k,j,i+1) : 1.0;
      const Real cvm = (q > 0 && Tim > 0.0 && eim > 0.0) ? (eim/Tim) : 1.0;
      const Real cvp = (q < n-1 && Tip > 0.0 && eip > 0.0) ? (eip/Tip) : 1.0;
      // ---- THE TWO-LEVEL SPLIT, see ck_impl_tau_min ------------------------------
      // A THIN row is an identity row whose right-hand side is the cell's own
      // bracketed implicit step.  It is decoupled from the tridiagonal in BOTH
      // directions: from the thin cell outwards because its row is the identity, and
      // from a thick neighbour inwards because the sweep zeroed its dB/dT, which is
      // the factor every off-diagonal entry pointing at it carries (see the jck block
      // in rt_chain_ck).  So the matrix stays the same M-matrix, one size smaller.
      const bool thin = !(thk_(m,k,j,i) > 0.0);
      Real a, b, c, d;
      if (seedp_) {
        // ---- problem/ck_impl_seed: THE INITIAL GUESS, as an identity row ----------
        // a = c = 0 and b = 1 make the Thomas sweep hand `d` back unchanged (cp = 0,
        // dp = d), so the guess reaches the gas through the same caps and the same
        // accounting as a Newton step, and the column stays decoupled.
        const Real emc = em_(m,k,j,i);
        const Real sc = src_(m,k,j,i);
        a = 0.0;
        b = 1.0;
        c = 0.0;
        if (seedm_ == 2) {
          d = CkThinSolve(ei, est_(m,k,j,i), sc, emc, bdt) - ei;
        } else {
          // the rt_semi_lin step of the default scheme: relax the cell's own emission
          // about the current state at lambda = 4 E/e, absorbed field lagged
          const Real de0 = est_(m,k,j,i) - ei;   // the step already applied this stage
          if (emc > 0.0) {
            const Real lam = 4.0*emc/ei;
            const Real x = lam*bdt;
            d = de0 + ((x > 1.0e-4) ? (sc/lam)*(-expm1(-x)) : sc*bdt);
          } else {
            d = de0 + sc*bdt;
          }
        }
      } else if (thin) {
        a = 0.0;
        b = 1.0;
        c = 0.0;
        d = CkThinSolve(ei, est_(m,k,j,i), src_(m,k,j,i), em_(m,k,j,i), bdt) - ei;
        Kokkos::atomic_add(&cnv_(7), 1.0);
      } else {
        // ---- problem/ck_impl_reuse_jac = 2: the SCALED CHORD.  Each column of the row
        // is rescaled by (T/T_build)^3 of the cell that column points at, which is the
        // dominant T dependence of the dB_b/dT the frozen entry carries.  Positive, so
        // the M-matrix property survives.  Mode 1 leaves the entries alone (the
        // rescaling is skipped by jscl_ = 0, not by a second branch, so the two modes
        // share one code path).
        Real s0 = 1.0, s1 = 1.0, s2 = 1.0;
        if (jscl_ > 0) {
          const Real tb = t0_(m,k,j,i);
          if (Ti > 0.0 && tb > 0.0) {
            const Real r = Ti/tb;
            s1 = r*r*r;
          }
          if (q > 0) {
            const Real tbm = t0_(m,k,j,i-1);
            if (Tim > 0.0 && tbm > 0.0) {
              const Real r = Tim/tbm;
              s0 = r*r*r;
            }
          }
          if (q < n-1) {
            const Real tbp = t0_(m,k,j,i+1);
            if (Tip > 0.0 && tbp > 0.0) {
              const Real r = Tip/tbp;
              s2 = r*r*r;
            }
          }
        }
        a = (q > 0) ? (-bdt*s0*jac_(m,0,k,j,i)/cvm) : 0.0;
        b = 1.0 - bdt*s1*jac_(m,1,k,j,i)/cvi;
        c = (q < n-1) ? (-bdt*s2*jac_(m,2,k,j,i)/cvp) : 0.0;
        d = -(ei - est_(m,k,j,i) - bdt*src_(m,k,j,i));
        if (!(eim > CkDeadE)) a = 0.0;            // no coupling to a dead neighbour
        if (!(eip > CkDeadE)) c = 0.0;
      }
      if (!(ei > CkDeadE)) {                      // a DEAD cell: identity row, no step
        a = 0.0;
        b = 1.0;
        c = 0.0;
        d = 0.0;
      }
      if (!(b > 0.0)) {
        bad = true;
        break;
      }
      const Real den = b - a*((q > 0) ? cp_(m,k,j,i-1) : 0.0);
      if (!(fabs(den) > 0.0)) {
        bad = true;
        break;
      }
      cp_(m,k,j,i) = c/den;
      dp_(m,k,j,i) = (d - a*((q > 0) ? dp_(m,k,j,i-1) : 0.0))/den;
    }
    if (bad) {
      Kokkos::atomic_add(&cnv_(3), 1.0);
      // the bounded fallback for the whole column: every cell takes the per-cell
      // bracketed implicit step of the thin branch, which is positivity-preserving and
      // has the same fixed point as the tridiagonal, only with the neighbours lagged.
      for (int i=ic; i<ie+1; ++i) {
        const Real ei = ei_(m,k,j,i);
        if (!(ei > 0.0)) continue;
        Real de = CkThinSolve(ei, est_(m,k,j,i), src_(m,k,j,i), em_(m,k,j,i), bdt) - ei;
        const Real lim = dcap*ei;
        if (de > lim) de = lim;
        if (de < -lim) de = -lim;
        u0(m,IEN,k,j,i) += de;
        dep_(m,k,j,i) += de;
      }
      return;
    }
    Real dn = 0.0;
    Real prev = 0.0;
    bool cap = false;
    for (int q=n-1; q>=0; --q) {
      const int i = ic + q;
      Real de = dp_(m,k,j,i) - cp_(m,k,j,i)*prev;
      // NaN guard; <cmath> is not available on the device side
      if (!(de == de)) de = 0.0;
      prev = de;
      const Real ei = ei_(m,k,j,i);
      const Real lim = dcap*ei;
      if (de > lim) {
        de = lim;
        cap = true;
      }
      if (de < -lim) {
        de = -lim;
        cap = true;
      }
      // the total-excursion bound, see ck_impl_demax
      const Real es = est_(m,k,j,i);
      if (detot > 0.0 && es > 0.0) {
        Real en = ei + de;
        const Real ehi = es*(1.0 + detot);
        const Real elo = es*(1.0 - detot);
        if (en > ehi) {
          en = ehi;
          cap = true;
        }
        if (en < elo) {
          en = elo;
          cap = true;
        }
        de = en - ei;
      }
      u0(m,IEN,k,j,i) += de;
      dep_(m,k,j,i) += de;
      const Real s = (ei > 0.0) ? fabs(de)/ei : 0.0;
      if (s > dn) dn = s;
    }
    Kokkos::atomic_max(&cnv_(1), dn);
    if (cap) Kokkos::atomic_add(&cnv_(2), 1.0);
  });
  Kokkos::deep_copy(hc, cnv_);
  ck_impl_last_dst = hc(1);
  ck_impl_ncap = static_cast<int>(hc(2));
  ck_impl_nfall = static_cast<int>(hc(3));
  ck_impl_nthin = static_cast<int>(hc(7));
  // the step test.  The residual of the NEW iterate is measured at the top of the next
  // pass, which is the only honest place to measure it.
  if (hc(1) <= dtol && hc(0) <= tol) return 0;
  return 1;
}

void CkConserveFlux(Mesh *pm, const Real bdt) {
  if (ck_impl_conserve <= 0 || rt_Fb_ptr == nullptr || rt_icut_ptr == nullptr ||
      ck_dep_ptr == nullptr || ck_src_ptr == nullptr || rt_Em_ptr == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  MeshBlockPack *pp = pm->pmb_pack;
  const int nmb1 = pp->nmb_thispack - 1;
  auto fb_ = *rt_Fb_ptr;
  auto em_ = *rt_Em_ptr;
  auto dep_ = *ck_dep_ptr;
  auto src_ = *ck_src_ptr;
  auto est_ = *ck_estar_ptr;
  auto icut_ = *rt_icut_ptr;
  auto done_ = (ck_done_ptr != nullptr) ? *ck_done_ptr : CkDum<DvceArray3D<Real>>("d");
  const bool hasdone = (ck_done_ptr != nullptr);
  const bool msk = ck_cad_partial && hasdone;
  DvceArray5D<Real> u0 = (pp->pmhd != nullptr) ? pp->pmhd->u0 : pp->phydro->u0;
  EOS_Data eos = (pp->pmhd != nullptr) ? pp->pmhd->peos->eos_data
                                       : pp->phydro->peos->eos_data;
  const bool pp_ = !ck_spherical;          // the ck's own divergence (see rt_apply)
  auto area1 = pp->pcoord->area.x1f;
  auto vol = pp->pcoord->volume;
  auto dx1 = pp->pcoord->dx1;
  const int nblk = fb_.extent_int(1);
  const bool all = (ck_impl_conserve >= 2);
  const bool fbnd = ck_impl_floorbound;
  const bool kd = ck_impl_kkt_demax && (ck_impl_demax > 0.0);
  const Real detot = ck_impl_demax;
  const Real h = bdt;
  if (ck_cons_stat_ptr == nullptr) {
    ck_cons_stat_ptr = new DvceArray1D<Real>("ck_cons_stat", 2);
  }
  auto st_ = *ck_cons_stat_ptr;
  Kokkos::deep_copy(DevExeSpace(), st_, 0.0);
  // pass 1, per cell: the retained power r V/h of a corrected cell (0 elsewhere) and
  // the cell's total emission (the block weights), into ck_cons_buf
  const int n1 = fb_.extent_int(2);
  if (ck_cons_buf_ptr != nullptr &&
      (ck_cons_buf_ptr->extent_int(0) != nmb1 + 1 ||
       ck_cons_buf_ptr->extent_int(4) != n1)) {
    delete ck_cons_buf_ptr;                    // the pack changed size
    ck_cons_buf_ptr = nullptr;
  }
  if (ck_cons_buf_ptr == nullptr) {
    ck_cons_buf_ptr = new DvceArray5D<Real>("ck_cons_buf", nmb1 + 1, 2,
                                            fb_.extent_int(3), fb_.extent_int(4), n1);
  }
  auto cb_ = *ck_cons_buf_ptr;
  par_for("ck_conserve", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (msk && done_(m,k,j) > 1.5) return;       // masked out of a partial call
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    Real csum = 0.0, nc = 0.0;
    for (int i=ic; i<=ie; ++i) {
      const Real v = pp_ ? dx1(m,k,j,i) : vol(m,k,j,i);
      const Real r = dep_(m,k,j,i) - h*src_(m,k,j,i);
      bool use = all;
      if (!use && r != 0.0) {
        const Real es = est_(m,k,j,i);
        const Real e = es + dep_(m,k,j,i);
        use = CkKktCell(eos, fbnd, kd, detot, u0(m,IDN,k,j,i), e, es, r);
      }
      Real rv = 0.0;
      if (use && r == r) {
        rv = r*v/h;
        csum += r*v;
        nc += 1.0;
        Real et = 0.0;
        for (int b=0; b<nblk; ++b) et += em_(m,b,i,k,j);
        cb_(m,1,k,j,i) = et;
      }
      cb_(m,0,k,j,i) = rv;
    }
    if (nc > 0.0) {
      Kokkos::atomic_add(&st_(0), nc);
      Kokkos::atomic_add(&st_(1), csum);
    }
    // whether pass 2 has anything to do in this column: slot i = 0, a ghost cell
    // (ic >= is > 0), which pass 1 never writes
    cb_(m,0,k,j,0) = (nc > 0.0) ? 1.0 : 0.0;
  });
  // pass 2, per (block, column): carry the block's share up the column
  par_for("ck_conserve_fb", DevExeSpace(), 0, nmb1, 0, nblk-1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int b, const int k, const int j) {
    if (msk && done_(m,k,j) > 1.5) return;
    const int ic = icut_(m,k,j);
    if (ic > ie || !(cb_(m,0,k,j,0) > 0.0)) return;
    Real c = 0.0;
    for (int i=ic; i<=ie; ++i) {
      const Real rv = cb_(m,0,k,j,i);
      if (rv != 0.0) {
        const Real et = cb_(m,1,k,j,i);
        const Real w = (et > 0.0) ? em_(m,b,i,k,j)/et : ((b == 0) ? 1.0 : 0.0);
        c += w*rv;
      }
      if (c != 0.0) fb_(m,b,i+1,k,j) -= c/(pp_ ? 1.0 : area1(m,k,j,i+1));
    }
  });
  if (ck_impl_verbose) {
    auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), st_);
    ck_cons_ncell = static_cast<int64_t>(hs(0));
    ck_cons_de = hs(1);
  }
}

void CkNonconvLoc(Mesh *pm, const Real bdt) {
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  if (ck_done_ptr == nullptr || rt_pb_ptr == nullptr || rt_icut_ptr == nullptr) return;
  auto done_ = *ck_done_ptr;
  auto ei_ = *ck_ei_ptr;
  auto est_ = *ck_estar_ptr;
  auto src_ = *ck_src_ptr;
  auto pb_ = *rt_pb_ptr;
  auto T_ = *rt_T_ptr;
  auto icut_ = *rt_icut_ptr;
  auto gid_ = pm->pmb_pack->pmb->mb_gid.d_view;
  const Real eps = ck_impl_norm_eps;
  const int nmax = ck_impl_ncloc;
  const int rank = global_variable::my_rank;
  const int ncyc = static_cast<int>(pm->ncycle);
  const bool csph_ = pm->use_cubed_sphere;
  auto &mbpan_ = pm->pmb_pack->pmb->mb_panel;
  auto &x2v_ = pm->pmb_pack->pcoord->x2v;
  auto &x3v_ = pm->pmb_pack->pcoord->x3v;
  DvceArray1D<int> cnt("ck_ncloc_cnt", 1);
  // the cells the convergence test counts as converged (KKT: on the floor or the demax
  // bound with the residual pointing through it) are skipped, as the test skips them
  const Real detot = ck_impl_demax;
  const bool fbk = ck_impl_floorbound;
  const bool kdk = ck_impl_kkt_demax && (detot > 0.0);
  EOS_Data eosk = (pm->pmb_pack->pmhd != nullptr)
                  ? pm->pmb_pack->pmhd->peos->eos_data
                  : pm->pmb_pack->phydro->peos->eos_data;
  auto u0k = (pm->pmb_pack->pmhd != nullptr) ? pm->pmb_pack->pmhd->u0
                                             : pm->pmb_pack->phydro->u0;
  par_for("ck_ncloc", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    if (done_(m,k,j) != 0.0) return;
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    Real emax = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      if (ei_(m,k,j,i) > emax) emax = ei_(m,k,j,i);
    }
    Real rn = -1.0;
    int iw = ic;
    for (int i=ic; i<ie+1; ++i) {
      const Real r = ei_(m,k,j,i) - est_(m,k,j,i) - bdt*src_(m,k,j,i);
      const Real s = fabs(r)/(ei_(m,k,j,i) + eps*emax);
      if ((fbk || kdk) && CkKktCell(eosk, fbk, kdk, detot, u0k(m,IDN,k,j,i),
                                    ei_(m,k,j,i), est_(m,k,j,i), r)) continue;
      if (s > rn) {
        rn = s;
        iw = i;
      }
    }
    const int n = Kokkos::atomic_fetch_add(&cnt(0), 1);
    if (n < nmax) {
      Real pth = 0.0, plat = 0.0, plon = 0.0;
      if (csph_) {
        atm_column::CSCellAngles(mbpan_.d_view(m), x2v_(m,j), x3v_(m,k), pth, plat, plon);
      }
      Kokkos::printf("### ck_ncloc ncycle=%d rank=%d gid=%d m=%d k=%d j=%d i=%d ie=%d "
                     "icut=%d res=%.3e p_bar=%.4e T=%.1f ptop_bar=%.4e lat=%.2f "
                     "lon=%.2f\n", ncyc, rank,
                     gid_(m), m, k, j, iw, ie, ic, rn, pb_(m,k,j,iw), T_(m,k,j,iw),
                     pb_(m,k,j,ie), plat*57.29577951308232, plon*57.29577951308232);
    }
  });
  Kokkos::fence();
}

void CkNcRefSnap(Mesh *pm) {
  if (ck_dep_ptr == nullptr) return;
  if (ck_ncref_ptr == nullptr) {
    auto d = *ck_dep_ptr;
    ck_ncref_ptr = new DvceArray4D<Real>("ck_ncref", d.extent(0), d.extent(1),
                                         d.extent(2), d.extent(3));
  }
  Kokkos::deep_copy(*ck_ncref_ptr, *ck_dep_ptr);
}

void CkNcRefCmp(Mesh *pm, const int nsnap, const int npass, const bool conv) {
  if (ck_dep_ptr == nullptr || ck_ncref_ptr == nullptr || rt_icut_ptr == nullptr) return;
  auto &indcs = pm->mb_indcs;
  const int ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  auto dep_ = *ck_dep_ptr;
  auto snp_ = *ck_ncref_ptr;
  auto ei_ = *ck_ei_ptr;
  auto pb_ = *rt_pb_ptr;
  auto T_ = *rt_T_ptr;
  auto icut_ = *rt_icut_ptr;
  auto dx1_ = pm->pmb_pack->pcoord->dx1;
  DvceArray1D<Real> acc("ck_ncref_acc", 8);
  par_for("ck_ncref", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    Real mx = 0.0, sd = 0.0, sr = 0.0, mxd = 0.0, mxd4 = 0.0;
    for (int i=ic; i<ie+1; ++i) {
      const Real e = ei_(m,k,j,i);
      const Real d = fabs(dep_(m,k,j,i) - snp_(m,k,j,i));
      if (e > 0.0 && d/e > mx) mx = d/e;
      // the region that counts (p > 1e-6 bar), and p > 1e-4 bar
      if (e > 0.0 && pb_(m,k,j,i) > 1.0e-6 && d/e > mxd) mxd = d/e;
      if (e > 0.0 && pb_(m,k,j,i) > 1.0e-4 && d/e > mxd4) mxd4 = d/e;
      sd += d*dx1_(m,k,j,i);
      sr += fabs(dep_(m,k,j,i))*dx1_(m,k,j,i);
    }
    Kokkos::atomic_max(&acc(0), mx);
    if (sr > 0.0) Kokkos::atomic_max(&acc(1), sd/sr);
    Kokkos::atomic_add(&acc(2), sd);
    Kokkos::atomic_add(&acc(3), sr);
    if (mx > 0.0) Kokkos::atomic_add(&acc(4), 1.0);
    if (mx > 1.0e-6) Kokkos::atomic_add(&acc(5), 1.0);
    Kokkos::atomic_max(&acc(6), mxd);
    Kokkos::atomic_max(&acc(7), mxd4);
  });
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), acc);
  const Real thr = 0.5*h(0);
  const int rank = global_variable::my_rank;
  const int ncyc = static_cast<int>(pm->ncycle);
  std::cout << "### ck_ncref ncycle=" << ncyc << " rank=" << rank << " snap=" << nsnap
            << " passes=" << npass << " conv=" << (conv ? 1 : 0) << " maxde_e=" << h(0)
            << " colheat_err_max=" << h(1) << " glob_heat_err="
            << ((h(3) > 0.0) ? h(2)/h(3) : 0.0) << " ncol_diff=" << h(4)
            << " ncol_gt1e-6=" << h(5) << " maxde_e_p1e-6=" << h(6)
            << " maxde_e_p1e-4=" << h(7) << std::endl;
  if (!(h(0) > 1.0e-6)) return;
  DvceArray1D<int> cnt("ck_ncref_cnt", 1);
  const bool csph_ = pm->use_cubed_sphere;
  auto &mbpan_ = pm->pmb_pack->pmb->mb_panel;
  auto &x2v_ = pm->pmb_pack->pcoord->x2v;
  auto &x3v_ = pm->pmb_pack->pcoord->x3v;
  par_for("ck_ncref_loc", DevExeSpace(), 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const int ic = icut_(m,k,j);
    if (ic > ie) return;
    for (int i=ic; i<ie+1; ++i) {
      const Real e = ei_(m,k,j,i);
      const Real d = dep_(m,k,j,i) - snp_(m,k,j,i);
      if (e > 0.0 && fabs(d)/e >= thr) {
        const int n = Kokkos::atomic_fetch_add(&cnt(0), 1);
        if (n < 6) {
          Real pth = 0.0, plat = 0.0, plon = 0.0;
          if (csph_) {
            atm_column::CSCellAngles(mbpan_.d_view(m), x2v_(m,j), x3v_(m,k), pth, plat,
                                     plon);
          }
          Kokkos::printf("### ck_ncref_cell ncycle=%d m=%d k=%d j=%d i=%d ic=%d "
                         "lat=%.2f lon=%.2f p_bar=%.4e T=%.1f de/e=%.3e dep/e=%.3e "
                         "snap/e=%.3e\n", ncyc, m, k, j, i, ic,
                         plat*57.29577951308232, plon*57.29577951308232, pb_(m,k,j,i),
                         T_(m,k,j,i), d/e, dep_(m,k,j,i)/e, snp_(m,k,j,i)/e);
        }
      }
    }
  });
  Kokkos::fence();
}

void picket_fence_two_stream_RT_body(Mesh *pm, Real bdt) {
  // ---- problem/rt_implicit_column: the merged tridiagonal Newton --------------------
  // Requirements, checked once.  The column solve linearises the GREY SPLIT sweep, reads
  // the running state out of u0 (rt_use_cons) between passes, and replaces -- not
  // supplements -- the semi-implicit per-cell relaxation, which therefore has to be the
  // step it is replacing.
  Conduction *pc = (pm->pmb_pack->pmhd != nullptr) ? pm->pmb_pack->pmhd->pcond
                                                   : pm->pmb_pack->phydro->pcond;
  if (rt_implicit_column == 3) {
    // ---- mode 3, the exact block-tridiagonal column solve ---------------------------
    // It needs the grey split sweep (it IS that sweep, differentiated), the conserved
    // state (it applies de to u0 and must read the post-RK state), the direct per-cell
    // source (its rows ARE the per-half-layer absorbed-minus-emitted balance), and the
    // centre-to-centre layers.  It does NOT need rad_implicit_x1: the radial conduction
    // stays a separate implicit tridiagonal, applied by its own task.
    static bool c3checked = false;
    if (!c3checked) {
      c3checked = true;
      if (!(rt_grey && rt_split) || rt_ck || !rt_use_cons ||
          !rt_semi_implicit || !rt_src_direct || rt_layer_legacy || rt_top_re) {
        std::cout << "### FATAL ERROR in two_stream_rt: problem/rt_implicit_column = 3 "
                  << "needs the GREY SPLIT sweep (rt_grey + rt_split, and NOT "
                  << "correlated-k or picket fence), problem/rt_use_cons, "
                  << "rt_src_direct, the centre-to-centre layers (rt_layer_legacy = "
                  << "false), rt_top_re = false and the semi-implicit "
                  << "apply it replaces (rt_semi_implicit).  Got rt_grey="
                  << rt_grey << " rt_split=" << rt_split << " rt_ck=" << rt_ck
                  << " rt_use_cons=" << rt_use_cons
                  << " rt_semi_implicit=" << rt_semi_implicit
                  << " rt_src_direct=" << rt_src_direct
                  << " rt_layer_legacy=" << rt_layer_legacy
                  << " rt_top_re=" << rt_top_re << "." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  // ---- problem/ck_implicit: THE BACKWARD-EULER CORRELATED-K COLUMN ------------------
  // See utils/two_stream_column_ck.hpp for the balance, the Jacobian and what is left in
  // the residual.  Each pass is ONE ORDINARY ck sweep, which re-evaluates the residual at
  // the running state (hence rt_use_cons); the tridiagonal accumulated inside that same
  // sweep then takes the Newton step.  The apply block deposits NOTHING under this switch
  // -- it only records S_i and e_i -- so the gas receives exactly the converged flux
  // divergence, with no relaxation factor and no rt_de_max clip.
  if (ck_implicit) {
    if (!rt_ck || !rt_use_cons || rt_layer_legacy || rt_implicit_column != 0) {
      std::cout << "### FATAL ERROR in two_stream_rt: problem/ck_implicit needs the "
                << "correlated-k sweep on the centre-to-centre layers (rt_ck, "
                << "!rt_layer_legacy), problem/rt_use_cons (the passes read the running "
                << "state out of u0), and the older column iteration off "
                << "(rt_implicit_column = 0).  Got rt_ck="
                << rt_ck << " rt_use_cons=" << rt_use_cons
                << " rt_layer_legacy=" << rt_layer_legacy
                << " rt_implicit_column=" << rt_implicit_column << "." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // the first call of a restarted run: the cross-call state the file carried, and the
    // ck_impl_xstep operator rebuilt when this call re-applies it (two_stream_ck_rst.hpp)
    CkRstBeginCall(pm);
    // problem/ck_impl_glob = ls_sub: the passes a column may spend on its sub-steps
    // (maxit per sub-step at every level 1, 2, 4, .. sub_max, plus the restarts)
    int nitmax = ck_impl_maxit;
    if (ck_impl_glob == 2) {
      int nsl = 0;
      for (int l=1; l<=ck_impl_sub_max; l*=2) nsl += l;
      nitmax = ck_impl_maxit*nsl + 2*nsl + 2;
      // ck_impl_esc: the extra passes a contracting (sub-)step may take
      if (ck_impl_esc > 0) nitmax += ck_impl_esc_extra*nsl;
    }
    int npass = nitmax;
    bool conv = false;
    // ck_impl_colskip: every column is live again at the start of a call
    // (ck_impl_every_thr: except the columns this partial call masks out, ck_done = 2)
    if (ck_done_ptr != nullptr) {
      if (ck_cad_partial) {
        Kokkos::deep_copy(DevExeSpace(), *ck_done_ptr, *ck_cad_mask_ptr);
      } else if (ck_impl_nosync) {
        Kokkos::deep_copy(DevExeSpace(), *ck_done_ptr, 0.0);   // stream-ordered, no fence
      } else {
        Kokkos::deep_copy(*ck_done_ptr, 0.0);
      }
    }
    // ck_impl_jreuse: the first pass of a call always builds the Jacobian
    ck_impl_prev_res = -1.0;
    // ck_impl_glob: every column starts the call on level 0 with no pending trial
    ck_impl_nrej = 0;
    ck_impl_nsub = 0;
    ck_impl_nsubfail = 0;
    ck_impl_jac_again = false;
    ck_impl_nesc = 0;
    ck_impl_ncoarse = 0;
    ck_impl_naa = 0;
    if (ck_impl_glob > 0 && ck_lsc_ptr != nullptr) {
      Kokkos::deep_copy(*ck_lsc_ptr, 0.0);
      if (ck_sacc_ptr != nullptr) Kokkos::deep_copy(*ck_sacc_ptr, 0.0);
    }
    // problem/ck_impl_warm: reset the per-call increment accumulator and, when the
    // switch is on, seed the gas with the previous call's converged increment
    {
      MeshBlockPack *pp0 = pm->pmb_pack;
      DvceArray5D<Real> u0w = (pp0->pmhd != nullptr) ? pp0->pmhd->u0 : pp0->phydro->u0;
      CkWarmSeed(pm, u0w);
    }
    // ck_impl_debug <= -2: the residual and the step of every pass, on the report line
    std::string hist;
    int ncsnap = -1;   // ck_impl_ncref: active columns at the latch, -1 = no latch
    int xleft = -1;    // ck2nd-1102 ck_impl_xpass: forced passes left (-1 = not started)
    for (int it=0; it<nitmax; ++it) {
      ck_impl_pass = it;
      // ck2nd-1102 problem/ck_top_live: the ghosts (the top slot's column above) re-formed
      // from the CURRENT iterate before every pass, not left at the end of the last stage
      if (ck_top_live && pm->pgen->user_bcs) (pm->pgen->user_bcs_func)(pm);
      picket_fence_two_stream_RT_pass(pm, bdt);
      MeshBlockPack *pp = pm->pmb_pack;
      DvceArray5D<Real> u0c = (pp->pmhd != nullptr) ? pp->pmhd->u0 : pp->phydro->u0;
      ++ck_impl_nsweep;
      const int stp = CkImplStep(pm, u0c, *rt_icut_ptr, *rt_T_ptr, pp->pcoord->dx1, bdt);
      // problem/ck_impl_jreuse: rebuild the chord Jacobian on the next pass only when the
      // max residual contracted by less than ck_impl_jreuse over this pass
      if (ck_impl_jreuse > 0.0 && ck_impl_glob == 0) {
        const Real rr = ck_impl_last_res;
        ck_impl_jac_again = (ck_impl_prev_res > 0.0)
                            && !(rr <= ck_impl_jreuse*ck_impl_prev_res);
        ck_impl_prev_res = rr;
        // ck_impl_jreuse_act: once few columns are left the rebuild is cheap (the
        // Jacobian kernels skip converged columns), and the few slow columns are the
        // ones a stale chord holds back, so rebuild whatever the max residual did
        const int ncol_ = pp->nmb_thispack*pm->mb_indcs.nx2*pm->mb_indcs.nx3;
        if (ck_impl_jreuse_act > 0.0 &&
            static_cast<Real>(ck_impl_nactive) < ck_impl_jreuse_act*ncol_) {
          ck_impl_jac_again = true;
        }
      }
      if (ck_impl_debug <= -2) {
        char hb[64];
        std::snprintf(hb, sizeof(hb), "%s%.2e/%.2e/%d", (it > 0) ? "," : "",
                      ck_impl_last_res, ck_impl_last_dst, ck_impl_nactive);
        hist += hb;
      }
      // ck2nd-1102 problem/ck_impl_xpass: n forced passes after the first convergence
      if (stp == 0 && ck_impl_xpass > 0 && !ck_force_step && xleft < 0) {
        xleft = ck_impl_xpass;
      }
      ck_force_step = false;
      if (stp == 0 && xleft > 0 && it + 1 < nitmax) {
        --xleft;
        ++ck_impl_nxpass;
        if (ck_cad_partial) {
          Kokkos::deep_copy(DevExeSpace(), *ck_done_ptr, *ck_cad_mask_ptr);
        } else {
          Kokkos::deep_copy(DevExeSpace(), *ck_done_ptr, 0.0);
        }
        ck_force_step = true;
        continue;
      }
      if (stp == 0) {
        npass = it + 1;
        conv = true;
        break;
      }
      // problem/ck_impl_ncref: latch the state a maxit = N call would have accepted
      if (ck_impl_ncref > 0 && it + 1 == ck_impl_ncref) {
        CkNcRefSnap(pm);
        ncsnap = ck_impl_nactive;
      }
    }
    ck_impl_last_it = npass;
    if (ncsnap >= 0) CkNcRefCmp(pm, ncsnap, npass, conv);
    // problem/ck_impl_pred_chk: the final state's own residual and gap (one more sweep,
    // no step; diagnostic only -- the gas is not touched)
    int pchk_act = -1;
    if (ck_impl_pred && ck_impl_pred_chk && conv) {
      if (ck_cad_partial) {
        Kokkos::deep_copy(*ck_done_ptr, *ck_cad_mask_ptr);
      } else {
        Kokkos::deep_copy(*ck_done_ptr, 0.0);
      }
      ck_impl_pass = npass;
      picket_fence_two_stream_RT_pass(pm, bdt);
      MeshBlockPack *pp = pm->pmb_pack;
      DvceArray5D<Real> u0c = (pp->pmhd != nullptr) ? pp->pmhd->u0 : pp->phydro->u0;
      ck_impl_evalonly = true;
      CkImplStep(pm, u0c, *rt_icut_ptr, *rt_T_ptr, pp->pcoord->dx1, bdt);
      ck_impl_evalonly = false;
      pchk_act = ck_impl_nactive;
    }
    // problem/ck_impl_conserve: the reported fluxes of the cells the call left on a bound
    CkConserveFlux(pm, bdt);
    if (ck_impl_nsubfail > 0) conv = false;   // ls_sub: a column failed at the finest L
    ck_impl_nonconv = conv ? 0 : 1;
    if (!conv && ck_impl_ncloc > 0) CkNonconvLoc(pm, bdt);
    ck_impl_pass = -1;
    // ck_impl_debug = -1: the report line from EVERY rank, tagged; the Newton loop is
    // rank-local, so rank 0's verdict is not the call's
    if (ck_impl_verbose && (global_variable::my_rank == 0 || ck_impl_debug < 0)) {
      std::cout << "### ck_implicit ncycle=" << pm->ncycle << " passes=" << npass
                << " rank=" << global_variable::my_rank
                << " res=" << ck_impl_last_res << " dstep=" << ck_impl_last_dst
                << " ckdesum=" << ck_impl_last_gap
                << " capped=" << ck_impl_ncap << " fallback=" << ck_impl_nfall
                << " thin=" << ck_impl_nthin << " active=" << ck_impl_nactive
                << " nsweep=" << ck_impl_nsweep
                << " rj=" << ck_impl_reuse_jac << " sd=" << ck_impl_seed
                << ((ck_impl_glob > 0) ? (" lsrej=" + std::to_string(ck_impl_nrej)
                    + " subrst=" + std::to_string(ck_impl_nsub)
                    + " subfail=" + std::to_string(ck_impl_nsubfail)) : std::string(""))
                << ((ck_impl_esc > 0) ? (" escx=" + std::to_string(ck_impl_nesc)
                    + " coarse=" + std::to_string(ck_impl_ncoarse)) : std::string(""))
                << ((ck_impl_aa > 0) ? (" aa=" + std::to_string(ck_impl_naa))
                    : std::string(""))
                << ((ck_impl_xstep > 0) ? (std::string(" xs=")
                    + (ck_impl_reuse_op ? "R" : "S") + " xsd="
                    + std::to_string(ck_impl_xs_dmax)) : std::string(""))
                << (ck_impl_pred ? (" pred=" + std::to_string(ck_impl_npred))
                    : std::string(""))
                << ((ck_dif_dtau > 0.0) ? (" difcol=" + std::to_string(ck_dif_ncol)
                    + " difdep=" + std::to_string((ck_dif_ncol > 0)
                        ? ck_dif_fsum/static_cast<double>(ck_dif_ncol) : 0.0))
                    : std::string(""))
                << ((pchk_act >= 0) ? (" pchk_active=" + std::to_string(pchk_act))
                    : std::string(""))
                << ((ck_impl_floorbound || ck_impl_kkt_demax)
                    ? (" floor=" + std::to_string(ck_impl_nfloor) + " kkt="
                       + std::to_string(ck_impl_nkkt)) : std::string(""))
                << ((ck_impl_conserve > 0) ? (" cons=" + std::to_string(ck_cons_ncell)
                    + " consde=" + std::to_string(ck_cons_de)) : std::string(""))
                << ((ck_impl_osc > 0) ? (" osc=" + std::to_string(ck_impl_nosc))
                    : std::string(""))
                << (conv ? "" : " NOT-CONVERGED")
                << ((ck_impl_debug <= -2) ? (" hist=" + hist) : std::string(""))
                << std::endl;
    }
    return;
  }
  // mode 3 solves and applies the column inside the pass; the radial conduction stays
  // a separate operator, applied by the ImplicitConduction task as usual.
  picket_fence_two_stream_RT_pass(pm, bdt);
}

void picket_fence_two_stream_RT(Mesh *pm, Real bdt) {
  picket_fence_two_stream_RT_body(pm, bdt);
  RTSourceLimiterReduce();
}

void CkCadStore(Mesh *pm, const Real dt, const bool partial, const bool diag) {
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  const CkCadState st = CkCadMakeState(pm);
  auto cad_ = *ck_cad_ptr;
  auto msk_ = *ck_cad_mask_ptr;
  auto dep_ = *ck_dep_ptr;
  auto jac_ = *ck_jac_ptr;
  auto thu_ = *ck_thu_ptr;
  auto em_ = *ck_em_ptr;
  auto trt_ = *rt_T_ptr;
  auto icut_ = *rt_icut_ptr;
  auto dx1_ = pm->pmb_pack->pcoord->dx1;
  auto sta_ = *ck_cad_stat_ptr;
  const Real idt = 1.0/dt;
  // ck2nd-1102 problem/ck_cad_drift: G and t0 per cell (a 1-element dummy when off)
  const bool drf = ck_cad_drift;
  if (drf && ck_cad2_ptr == nullptr) {
    // ck2ndt-1102 ck_cad_drift2: slots 6 (the last slope), 7 (its midpoint time)
    const int ns2 = (ck_cad_drift2 || ck_cad_cont) ? 8 : 6;
    ck_cad2_ptr = new DvceArray5D<Real>("ck_cad2", cad_.extent(0), ns2, cad_.extent(2),
                                        cad_.extent(3), cad_.extent(4));
    // t0 = "none yet" (a restart's Q0 has no measured drift)
    Kokkos::deep_copy(Kokkos::subview(*ck_cad2_ptr, Kokkos::ALL, 1, Kokkos::ALL,
                                      Kokkos::ALL, Kokkos::ALL), -1.0e300);
    if (ck_cad_drift2) {
      Kokkos::deep_copy(Kokkos::subview(*ck_cad2_ptr, Kokkos::ALL, 7, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL), -1.0e300);
    }
  }
  auto cd2_ = drf ? *ck_cad2_ptr : cad_;    // (not read when off)
  const Real tev = ck_cad_teval;
  if (drf) ++ck_cad_ndrift;
  // ck2nd-1102 problem/ck_cad_tri: the drift of the TRIDIAGONAL model (the neighbours'
  // own Jacobian terms J0 (i-1), J2 (i+1) taken out too); slot 5 = current T, slot 4 = G
  const bool tri = drf && ck_cad_tri;
  const bool dr2 = drf && ck_cad_drift2 && !tri;   // ck2ndt-1102 (ck_2nd.hpp (3c))
  // ck2ndt-1102 problem/ck_cad_cont (ck_2nd.hpp (8b)): the continuous probe store
  if (drf && ck_cad_cont && ck_cad_un.size() > 0 && !partial) {
    CkCadState sn = st;
    sn.u0 = ck_cad_un;
    const Real pest = ((ck_cad_dt > 0.0) ? ck_cad_dt
                       : static_cast<Real>(ck_impl_every)*pm->dt)*ck_cad_cont_len;
    ck_cad_cont_l = pest;
    par_for("ck_cad_cont", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = st.u0(m,IDN,k,j,i);
      Real tk = 0.0, cv = 0.0;
      const bool ok = st.TCv(m, k, j, i, d, st.Eint(m, k, j, i), tk, cv);
      Real tn = 0.0, cvn = 0.0;
      const bool okn = sn.TCv(m, k, j, i, sn.u0(m,IDN,k,j,i), sn.Eint(m, k, j, i), tn, cvn);
      Real q = 0.0, dd = 0.0;
      if (ok && okn && i >= icut_(m,k,j)) {
        q = dep_(m,k,j,i)*idt;
        if (thu_(m,k,j,i) > 0.0) {
          dd = jac_(m,1,k,j,i);
        } else {
          const Real tr = trt_(m,k,j,i);
          dd = (tr > 0.0) ? -4.0*em_(m,k,j,i)/tr : 0.0;
        }
        if (!(dd <= 0.0)) dd = 0.0;
        if (!(q == q)) q = 0.0;
      }
      const Real sn_ = (ok && okn) ? q + dd*(tn - tk) : 0.0;   // the rate at U^n
      const Real q0 = cad_(m,0,k,j,i), d0 = cad_(m,1,k,j,i), t00 = cad_(m,2,k,j,i);
      const Real g0 = cd2_(m,0,k,j,i), tp = cd2_(m,1,k,j,i), so = cd2_(m,6,k,j,i);
      const bool hav = okn && (t00 > 0.0) && (q0 != 0.0 || d0 != 0.0);
      const bool tv = (tp > -1.0e299) && (tev > tp);
      Real aold = sn_, gm = 0.0;
      if (hav) {
        aold = q0 + d0*(tn - t00) + (tv ? g0*(tev - tp) : 0.0);
        if (tv) gm = (sn_ - so - d0*(tn - t00))/(tev - tp);
      }
      if (!(aold == aold)) aold = sn_;
      Real g = gm;
      if (!(g == g)) g = 0.0;
      // the model rate S_n + G tau + C max(0, 1 - tau/L), C = M_old - S_n (slot 7)
      cd2_(m,7,k,j,i) = (ok && okn) ? aold - sn_ : 0.0;
      cad_(m,0,k,j,i) = (ok && okn) ? sn_ : 0.0;
      cad_(m,1,k,j,i) = dd;
      cad_(m,2,k,j,i) = (ok && okn) ? tn : 0.0;
      cad_(m,3,k,j,i) = d;
      cd2_(m,0,k,j,i) = (ok && okn) ? g : 0.0;
      cd2_(m,1,k,j,i) = tev;
      cd2_(m,6,k,j,i) = sn_;
    });
    return;
  }
  if (tri) {
    par_for("ck_cad_tri_t", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = st.u0(m,IDN,k,j,i);
      const Real e = st.Eint(m, k, j, i);
      Real tk = 0.0, cv = 0.0;
      const bool ok = st.TCv(m, k, j, i, d, e, tk, cv);
      cd2_(m,5,k,j,i) = ok ? tk : 0.0;
    });
    par_for("ck_cad_tri_g", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (partial && msk_(m,k,j) > 0.5) return;
      Real g = 0.0;
      const Real tk = cd2_(m,5,k,j,i);
      const Real q0 = cad_(m,0,k,j,i), d0 = cad_(m,1,k,j,i), t00 = cad_(m,2,k,j,i);
      const Real tp = cd2_(m,1,k,j,i);
      const Real tl = tev - tp;
      const Real q = (tk > 0.0 && i >= icut_(m,k,j)) ? dep_(m,k,j,i)*idt : 0.0;
      if (tk > 0.0 && t00 > 0.0 && tp > -1.0e299 && tl > 0.0 && (q0 != 0.0 || d0 != 0.0)
          && q == q && q != 0.0) {
        Real r = q - q0 - d0*(tk - t00);
        if (i > is && cad_(m,2,k,j,i-1) > 0.0 && cd2_(m,5,k,j,i-1) > 0.0) {
          r -= cd2_(m,2,k,j,i)*(cd2_(m,5,k,j,i-1) - cad_(m,2,k,j,i-1));
        }
        if (i < ie && cad_(m,2,k,j,i+1) > 0.0 && cd2_(m,5,k,j,i+1) > 0.0) {
          r -= cd2_(m,3,k,j,i)*(cd2_(m,5,k,j,i+1) - cad_(m,2,k,j,i+1));
        }
        g = r/tl;
        if (!(g == g)) g = 0.0;
      }
      cd2_(m,4,k,j,i) = g;
    });
  }
  par_for("ck_cad_store", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (partial && msk_(m,k,j) > 0.5) return;
    const Real d = st.u0(m,IDN,k,j,i);
    const Real e = st.Eint(m, k, j, i);
    Real tk = 0.0, cv = 0.0;
    const bool ok = st.TCv(m, k, j, i, d, e, tk, cv);
    Real q = 0.0, dd = 0.0;
    if (ok && i >= icut_(m,k,j)) {
      q = dep_(m,k,j,i)*idt;
      if (thu_(m,k,j,i) > 0.0) {
        dd = jac_(m,1,k,j,i);
      } else {
        const Real tr = trt_(m,k,j,i);
        dd = (tr > 0.0) ? -4.0*em_(m,k,j,i)/tr : 0.0;
      }
      if (!(dd <= 0.0)) dd = 0.0;          // positive or NaN: no implicit term
      if (!(q == q)) q = 0.0;
      if (diag && cad_(m,2,k,j,i) > 0.0) {
        const Real qp = cad_(m,0,k,j,i) + cad_(m,1,k,j,i)*(tk - cad_(m,2,k,j,i));
        const Real w = dx1_(m,k,j,i);
        Kokkos::atomic_add(&sta_(6), w*fabs(q - qp));
        Kokkos::atomic_add(&sta_(7), w*fabs(q));
      }
    }
    if (drf) {
      // the drift of the part of the rate the diagonal model missed since the last store
      Real g = 0.0;
      const Real q0 = cad_(m,0,k,j,i), d0 = cad_(m,1,k,j,i), t00 = cad_(m,2,k,j,i);
      const Real tp = cd2_(m,1,k,j,i);
      const Real tl = tev - tp;
      if (tri) {
        g = cd2_(m,4,k,j,i);
      } else if (ok && t00 > 0.0 && tp > -1.0e299 && tl > 0.0 &&
                 (q0 != 0.0 || d0 != 0.0) && (q != 0.0 || dd != 0.0)) {
        g = (q - q0 - d0*(tk - t00))/tl;
        if (!(g == g)) g = 0.0;
      }
      if (dr2) {
        // g is the mean slope over [tp, tev], i.e. the slope at its midpoint; with the
        // previous mean slope (slot 6, midpoint slot 7) the curvature H, the slope at tev
        // G = g + H tl/2, and the model rate Q0 + G tau + H tau^2/2
        Real hh = 0.0;
        const Real mid = tev - 0.5*tl;
        const Real sp = cd2_(m,6,k,j,i), mp = cd2_(m,7,k,j,i);
        const bool gv = (g != 0.0) && tl > 0.0 && tp > -1.0e299;
        if (gv && mp > -1.0e299 && mid > mp) {
          hh = (g - sp)/(mid - mp);
          if (!(hh == hh)) hh = 0.0;
        }
        cd2_(m,6,k,j,i) = g;
        cd2_(m,7,k,j,i) = gv ? mid : -1.0e300;
        cd2_(m,4,k,j,i) = hh;
        g += 0.5*hh*tl;
      }
      cd2_(m,0,k,j,i) = g;
      cd2_(m,1,k,j,i) = tev;
      if (tri) {
        // the neighbour couplings of the Newton tridiagonal (thick cells; >= 0)
        Real j0 = 0.0, j2 = 0.0;
        if (ok && i >= icut_(m,k,j) && thu_(m,k,j,i) > 0.0) {
          if (i > icut_(m,k,j)) j0 = jac_(m,0,k,j,i);
          if (i < ie) j2 = jac_(m,2,k,j,i);
        }
        cd2_(m,2,k,j,i) = (j0 > 0.0) ? j0 : 0.0;
        cd2_(m,3,k,j,i) = (j2 > 0.0) ? j2 : 0.0;
      }
    }
    cad_(m,0,k,j,i) = q;
    cad_(m,1,k,j,i) = dd;
    cad_(m,2,k,j,i) = ok ? tk : 0.0;
    cad_(m,3,k,j,i) = d;
  });
}

int CkCadLin(Mesh *pm, const Real dt, const bool diag) {
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pm->pmb_pack->nmb_thispack - 1;
  const CkCadState st = CkCadMakeState(pm);
  auto cad_ = *ck_cad_ptr;
  auto msk_ = *ck_cad_mask_ptr;
  auto dx1_ = pm->pmb_pack->pcoord->dx1;
  auto sta_ = *ck_cad_stat_ptr;
  const Real thr = ck_impl_every_thr;
  const bool thr_on = (thr > 0.0);
  // ck2nd-1102 problem/ck_cad_drift (G = 0 until a store has measured one)
  const bool drf = ck_cad_drift && (ck_cad2_ptr != nullptr);
  auto cd2_ = drf ? *ck_cad2_ptr : cad_;    // (not read when off)
  const Real tev = ck_cad_teval;
  const bool tri = drf && ck_cad_tri;       // ck2nd-1102 ck_cad_tri
  const bool dr2 = drf && ck_cad_drift2 && !tri;   // ck2ndt-1102 ck_cad_drift2
  const int nlit = ck_cad_lin_iter;          // ck2ndt-1102 ck_cad_lin_iter
  const bool cnt = drf && ck_cad_cont && ck_cad_probe && !tri;   // ck2ndt-1102
  const Real cnl = ck_cad_cont_l;
  Real capf = ck_impl_dtmax;
  if (ck_impl_demax > 0.0 && ck_impl_demax < capf) capf = ck_impl_demax;
  const int nk = ke - ks + 1;
  const int nj = je - js + 1;
  const int nkj = nk*nj;
  const int nlg = (nmb1 + 1)*nkj;
#if defined(KOKKOS_ENABLE_HIP) || defined(KOKKOS_ENABLE_CUDA)
  Kokkos::TeamPolicy<> tpol(DevExeSpace(), nlg, 64);
#else
  Kokkos::TeamPolicy<> tpol(DevExeSpace(), nlg, Kokkos::AUTO);
#endif
  Kokkos::parallel_for("ck_cad_lin", tpol, KOKKOS_LAMBDA(TeamMember_t tm) {
    const int lr = tm.league_rank();
    const int m = lr/nkj;
    const int k = (lr - m*nkj)/nj + ks;
    const int j = (lr - m*nkj) % nj + js;
    // the current T of every cell (slot 4) and, with the guard, the column's deviation
    Real dv = 0.0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, is, ie+1),
    [&](const int i, Real &mx) {
      const Real d = st.u0(m,IDN,k,j,i);
      const Real e = st.Eint(m, k, j, i);
      Real tk = 0.0, cv = 0.0;
      const bool ok = st.TCv(m, k, j, i, d, e, tk, cv);
      cad_(m,4,k,j,i) = ok ? tk : 0.0;
      cad_(m,5,k,j,i) = ok ? cv : 0.0;
      const Real t0 = cad_(m,2,k,j,i);
      if (!(cad_(m,0,k,j,i) != 0.0 || cad_(m,1,k,j,i) != 0.0)) return;
      Real r = 0.0;
      if (!ok || !(t0 > 0.0)) {
        r = 1.0e30;                          // a cell went bad: re-solve the column
      } else {
        r = fabs(tk/t0 - 1.0);
        const Real d0 = cad_(m,3,k,j,i);
        if (d0 > 0.0) {
          const Real rd = fabs(d/d0 - 1.0);
          if (rd > r) r = rd;
        }
      }
      if (r > mx) mx = r;
    }, Kokkos::Max<Real>(dv));
    if (!(dv > 0.0)) dv = 0.0;
    if (thr_on && dv > thr) {
      Kokkos::single(Kokkos::PerTeam(tm), [&]() {
        msk_(m,k,j) = 0.0;
        Kokkos::atomic_add(&sta_(0), 1.0);
        Kokkos::atomic_max(&sta_(4), dv);
      });
      return;
    }
    tm.team_barrier();
    Real sde = 0.0, sq = 0.0, saq = 0.0;
    int ncp = 0;
    if (tri) {
      // ck2nd-1102 ck_cad_tri: the linearised BE step with the Newton tridiagonal
      //   (cv - dt D) dT_i - dt J0 dT_{i-1} - dt J2 dT_{i+1}
      //     = dt (Q0 + G (t - t0) + D (T_i - T0_i) + J0 (T_{i-1} - T0_{i-1})
      //           + J2 (T_{i+1} - T0_{i+1})),
      // Thomas on one lane (slots 4/5 of ck_cad2 = c', r'); caps as the diagonal step
      Kokkos::single(Kokkos::PerTeam(tm), [&]() {
        auto act = [&](const int i) -> bool {
          return (cad_(m,0,k,j,i) != 0.0 || cad_(m,1,k,j,i) != 0.0) &&
                 cad_(m,4,k,j,i) > 0.0 && cad_(m,5,k,j,i) > 0.0;
        };
        Real cp = 0.0, rp = 0.0;
        for (int i=is; i<=ie; ++i) {
          Real a = 0.0, b = 1.0, c = 0.0, r = 0.0;
          if (act(i)) {
            const Real dd = cad_(m,1,k,j,i);
            const Real tk = cad_(m,4,k,j,i);
            const Real cv = cad_(m,5,k,j,i);
            Real q0 = cad_(m,0,k,j,i) + cd2_(m,0,k,j,i)*(tev - cd2_(m,1,k,j,i));
            r = q0 + dd*(tk - cad_(m,2,k,j,i));
            b = cv - dd*dt;
            if (i > is && act(i-1)) {
              const Real j0 = cd2_(m,2,k,j,i);
              r += j0*(cad_(m,4,k,j,i-1) - cad_(m,2,k,j,i-1));
              a = -dt*j0;
            }
            if (i < ie && act(i+1)) {
              const Real j2 = cd2_(m,3,k,j,i);
              r += j2*(cad_(m,4,k,j,i+1) - cad_(m,2,k,j,i+1));
              c = -dt*j2;
            }
            r *= dt;
          }
          const Real den = b - a*cp;
          cp = c/den;
          rp = (r - a*rp)/den;
          cd2_(m,4,k,j,i) = cp;
          cd2_(m,5,k,j,i) = rp;
        }
        Real dTn = 0.0;
        for (int i=ie; i>=is; --i) {
          const Real dT = cd2_(m,5,k,j,i) - cd2_(m,4,k,j,i)*dTn;
          dTn = dT;
          Real de = 0.0;
          if (act(i)) {
            de = cad_(m,5,k,j,i)*dT;
            const Real e = st.Eint(m, k, j, i);
            const Real lim = capf*e;
            if (de > lim) {
              de = lim;
              ncp += 1;
            } else if (de < -lim) {
              de = -lim;
              ncp += 1;
            }
            if (!(de == de)) de = 0.0;
          }
          cd2_(m,5,k,j,i) = de;    // applied below, after every row has read its T
        }
        for (int i=is; i<=ie; ++i) {
          const Real de = cd2_(m,5,k,j,i);
          if (de != 0.0) st.u0(m,IEN,k,j,i) += de;
          cad_(m,5,k,j,i) = de;
        }
      });
    } else {
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, is, ie+1),
    [&](const int i, int &nc) {
      Real q0 = cad_(m,0,k,j,i);
      const Real dd = cad_(m,1,k,j,i);
      if (drf) {
        const Real tau = tev - cd2_(m,1,k,j,i);
        q0 += cd2_(m,0,k,j,i)*tau;
        if (cnt && tau < cnl) q0 += cd2_(m,7,k,j,i)*(1.0 - tau/cnl);   // ck_cad_cont
        if (dr2) q0 += 0.5*cd2_(m,4,k,j,i)*tau*tau;   // ck2ndt-1102 ck_cad_drift2
      }
      const Real tk = cad_(m,4,k,j,i);
      const Real cv = cad_(m,5,k,j,i);
      Real de = 0.0;
      if ((q0 != 0.0 || dd != 0.0) && tk > 0.0 && cv > 0.0) {
        const Real dT = (q0 + dd*(tk - cad_(m,2,k,j,i)))*dt/(cv - dd*dt);
        de = cv*dT;
        const Real e = st.Eint(m, k, j, i);
        // ck2ndt-1102 problem/ck_cad_lin_iter = n: the BE step of the model solved by n
        // more Newton iterations in e (T(e) of the EOS) instead of one linearisation in T
        // (the unconverged stage leaves an O(de^2) error that a stiff cell does not damp)
        for (int it=0; it<nlit; ++it) {
          const Real e1 = e + de;
          Real t1 = 0.0, c1 = 0.0;
          if (!st.TCv(m, k, j, i, st.u0(m,IDN,k,j,i), e1, t1, c1)) break;
          const Real r = de - dt*(q0 + dd*(t1 - cad_(m,2,k,j,i)));
          const Real ddn = r/(1.0 - dt*dd/c1);
          if (!(ddn == ddn)) break;
          de -= ddn;
        }
        const Real lim = capf*e;
        if (de > lim) {
          de = lim;
          nc += 1;
        } else if (de < -lim) {
          de = -lim;
          nc += 1;
        }
        if (!(de == de)) de = 0.0;
        st.u0(m,IEN,k,j,i) += de;
      }
      cad_(m,5,k,j,i) = de;
    }, ncp);
    }
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      msk_(m,k,j) = 2.0;
      Kokkos::atomic_max(&sta_(4), dv);
    });
    if (!diag) return;
    tm.team_barrier();
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, is, ie+1),
    [&](const int i, Real &sm) {
      sm += dx1_(m,k,j,i)*cad_(m,5,k,j,i);
    }, sde);
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, is, ie+1),
    [&](const int i, Real &sm) {
      sm += dx1_(m,k,j,i)*cad_(m,0,k,j,i)*dt;
    }, sq);
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(tm, is, ie+1),
    [&](const int i, Real &sm) {
      sm += dx1_(m,k,j,i)*fabs(cad_(m,0,k,j,i))*dt;
    }, saq);
    Kokkos::single(Kokkos::PerTeam(tm), [&]() {
      Kokkos::atomic_add(&sta_(1), sde);
      Kokkos::atomic_add(&sta_(2), sq);
      Kokkos::atomic_add(&sta_(3), static_cast<Real>(ncp));
      Kokkos::atomic_add(&sta_(5), saq);
    });
  });
  if (!thr_on && !diag) return 0;
  auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), sta_);
  return static_cast<int>(hs(0));
}

// ck2ndt-1102 problem/ck_cad_probe (ck_2nd.hpp (8)): the probe call of pdt on U^n has
// just run; store its rate as the cadence model at time tev (allocating on first use)
void CkCadProbeStore(Mesh *pm, const Real pdt, const Real tev) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  if (ck_cad_stat_ptr == nullptr) {
    ck_cad_stat_ptr = new DvceArray1D<Real>("ck_cad_stat", 8);
  }
  if (ck_cad_ptr == nullptr) {
    const int n1 = indcs.nx1 + 2*indcs.ng;
    const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
    const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
    ck_cad_ptr = new DvceArray5D<Real>("ck_cad", pmbp->nmb_thispack, 6, n3, n2, n1);
    ck_cad_mask_ptr = new DvceArray3D<Real>("ck_cad_mask", pmbp->nmb_thispack, n3, n2);
    Kokkos::deep_copy(*ck_cad_mask_ptr, 2.0);
  }
  ck_cad_teval = tev;
  CkCadStore(pm, pdt, false, false);
  ++ck_cad_nfull;
  ++ck_cad_nprobe;
}

void CkCadStep(Mesh *pm, const Real dt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ncol = pmbp->nmb_thispack*indcs.nx2*indcs.nx3;
  const bool diag = ck_impl_verbose;
  CkRstCadRestore(pm);     // restart: Q0, D, T0, rho0 (two_stream_ck_rst.hpp)
  ck_cad_teval = pm->time + pm->dt;   // ck2nd-1102 ck_cad_drift (both T2 stages: c = 1)
  // rkgfuse-1102 ck_strang: the second half of a cycle is never a scheduled full call
  // ck2ndt-1102 problem/ck_cad_dt: the schedule in simulated time (ck_2nd.hpp (7))
  const bool due = (ck_cad_dt > 0.0)
                   ? CkSchedCross(pm->time, pm->dt, ck_cad_dt)
                   : (pm->ncycle % static_cast<int64_t>(ck_impl_every) == 0);
  const bool full = (ck_cad_ptr == nullptr) || (due && ck_cad_half != 2);
  if (ck_cad_stat_ptr == nullptr) {
    ck_cad_stat_ptr = new DvceArray1D<Real>("ck_cad_stat", 8);
  }
  Kokkos::deep_copy(DevExeSpace(), *ck_cad_stat_ptr, 0.0);
  const bool had = (ck_cad_ptr != nullptr);
  int nref = 0;
  if (full) {
    ck_cad_partial = false;
    picket_fence_two_stream_RT(pm, dt);
    if (ck_cad_ptr == nullptr) {
      const int n1 = indcs.nx1 + 2*indcs.ng;
      const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
      const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
      ck_cad_ptr = new DvceArray5D<Real>("ck_cad", pmbp->nmb_thispack, 6, n3, n2, n1);
      ck_cad_mask_ptr = new DvceArray3D<Real>("ck_cad_mask", pmbp->nmb_thispack, n3, n2);
      Kokkos::deep_copy(*ck_cad_mask_ptr, 2.0);
    }
    CkCadStore(pm, dt, false, diag && had);
    ++ck_cad_nfull;
  } else {
    nref = CkCadLin(pm, dt, diag);
    ++ck_cad_nlin;
    if (nref > 0) {
      ck_cad_partial = true;
      picket_fence_two_stream_RT(pm, dt);
      CkCadStore(pm, dt, true, diag);
      ck_cad_partial = false;
      ck_impl_xs_cyc = -1;          // the operator this call stored is partial: drop it
      ++ck_cad_nguard;
    }
    ck_cad_fref += static_cast<double>(nref)/static_cast<double>(ncol);
  }
  if (diag && (global_variable::my_rank == 0 || ck_impl_debug < 0)) {
    auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), *ck_cad_stat_ptr);
    std::printf("### ck_cadence ncycle=%d rank=%d mode=%s nref=%d ncol=%d fref=%.4e "
                "dvmax=%.3e edef=%.4e sde=%.6e sq=%.6e ncap=%d linerr=%.4e "
                "nfull=%lld nlin=%lld nguard=%lld\n",
                pm->ncycle, global_variable::my_rank, full ? "F" : (nref > 0 ? "G" : "L"),
                nref, ncol, static_cast<double>(nref)/ncol, hs(4),
                (hs(5) > 0.0) ? (hs(1) - hs(2))/hs(5) : 0.0, hs(1), hs(2),
                static_cast<int>(hs(3)), (hs(7) > 0.0) ? hs(6)/hs(7) : -1.0,
                static_cast<long long>(ck_cad_nfull), static_cast<long long>(ck_cad_nlin),
                static_cast<long long>(ck_cad_nguard));
  }
}

//----------------------------------------------------------------------------------------
// ck2nd-1102 problem/ck_vef_tint (utils/ck_2nd.hpp (4)): the VEF factors of a storing
// pass extrapolated linearly in time from the last two refreshes
namespace {
DvceArray5D<Real> *vt_sk = nullptr, *vt_gk = nullptr, *vt_sp = nullptr, *vt_gp = nullptr;
Real vt_tk = 0.0, vt_tp = 0.0;
int vt_n = 0;
}  // namespace

void CkVefTint(Mesh *pm, const bool refreshed, const bool storing) {
  if (ck_vef_s_ptr == nullptr || ck_vef_g_ptr == nullptr) return;
  auto vs_ = *ck_vef_s_ptr;
  auto vg_ = *ck_vef_g_ptr;
  const Real t = pm->time + pm->dt;
  if (vt_sk == nullptr || vt_sk->extent(0) != vs_.extent(0)) {
    auto mk = [&](const char *lab) {
      return new DvceArray5D<Real>(std::string(lab), vs_.extent(0), vs_.extent(1), vs_.extent(2),
                                   vs_.extent(3), vs_.extent(4));
    };
    vt_sk = mk("ck_vt_sk"); vt_gk = mk("ck_vt_gk");
    vt_sp = mk("ck_vt_sp"); vt_gp = mk("ck_vt_gp");
    vt_n = 0;
  }
  if (refreshed) {
    if (vt_n > 0) {
      Kokkos::deep_copy(*vt_sp, *vt_sk);
      Kokkos::deep_copy(*vt_gp, *vt_gk);
      vt_tp = vt_tk;
    }
    Kokkos::deep_copy(*vt_sk, vs_);
    Kokkos::deep_copy(*vt_gk, vg_);
    vt_tk = t;
    ++vt_n;
    return;
  }
  if (!storing || vt_n < 2 || !(vt_tk > vt_tp)) return;
  const Real w = (t - vt_tk)/(vt_tk - vt_tp);
  auto sk = *vt_sk, gk = *vt_gk, sp = *vt_sp, gp = *vt_gp;
  const int n0 = static_cast<int>(vs_.extent(0)), n1 = static_cast<int>(vs_.extent(1));
  const int n2 = static_cast<int>(vs_.extent(2)), n3 = static_cast<int>(vs_.extent(3));
  const int n4 = static_cast<int>(vs_.extent(4));
  par_for("ck_vef_tint", DevExeSpace(), 0, n0 - 1, 0, n1 - 1, 0, n2 - 1, 0, n3 - 1,
  KOKKOS_LAMBDA(const int m, const int c, const int q, const int k) {
    for (int j=0; j<n4; ++j) {
      const Real a = sk(m,c,q,k,j), b = sp(m,c,q,k,j);
      Real x = a;
      if (a*b > 0.0) {
        x = a + (a - b)*w;
        if (!(x*a > 0.0)) x = a;
      }
      vs_(m,c,q,k,j) = x;
      const Real ga = gk(m,c,q,k,j), gb = gp(m,c,q,k,j);
      Real y = ga;
      if (ga*gb > 0.0) {
        y = ga + (ga - gb)*w;
        if (!(y*ga > 0.0)) y = ga;
      }
      vg_(m,c,q,k,j) = y;
    }
  });
  ++ck_vef_ntint;
}

}  // namespace two_stream_rt
