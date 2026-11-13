//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_fluxes.cpp
//! \brief Calculate fluxes of the conserved variables, and area-averaged electric fields
//! E = - (v X B) on cell faces for mhd.  Fluxes are stored in face-centered vector
//! 'uflx', while electric fields are stored in individual arrays: e2x1,e3x1 on x1-faces;
//! e1x2,e3x2 on x2-faces; e1x3,e2x3 on x3-faces.

#include <iostream>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd.hpp"
#include "hydro/hydro_wball.hpp"
#include "reconstruct/wb_wall.hpp"
#include "eos/eos.hpp"
#include "reconstruct/dc.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/plm_c1.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
#include "reconstruct/mignone_curvilinear.hpp"
#include "mhd/rsolvers/advect_mhd.hpp"
#include "mhd/rsolvers/cs_lowbeta_fallback.hpp"
#include "mhd/rsolvers/boris_hlle_mhd.hpp"
#include "mhd/rsolvers/llf_mhd.hpp"
#include "mhd/rsolvers/hlle_mhd.hpp"
#include "mhd/rsolvers/hlld_mhd.hpp"
#include "mhd/rsolvers/lhlld_mhd.hpp"
#include "mhd/rsolvers/llf_srmhd.hpp"
#include "mhd/rsolvers/hlle_srmhd.hpp"
#include "mhd/rsolvers/llf_grmhd.hpp"
#include "mhd/rsolvers/hlle_grmhd.hpp"
#include "bvals/cs_seam_states.hpp"
// #include "mhd/rsolvers/roe_mhd.hpp"

namespace mhd {

//----------------------------------------------------------------------------------------
//! \fn void UctFaceStore
//! \brief <mhd>/emf_scheme = uct | uct_hlld (lhlld2-1102, ucthlld-1102): store, per face,
//! what the UCT corner EMF in MHD::CornerE needs (Mignone & Del Zanna 2021, JCP 424,
//! 109748, "MD21").  The Riemann solver's induction flux is written as (MD21 eq. 29)
//!   F^[B_t] = a^L F^L + a^R F^R - (d^R B_t^R - d^L B_t^L),   F = v_n B_t - v_t B_n,
//! one (a, d) pair for every face-parallel component B_t.  ONE store for every solver:
//!  hll = true (emf_scheme = uct, rsolver = hlle, a spherical-polar polar-row face where
//!   polar_hlle_rows replaced HLLD/LHLLD by HLLE (spdef-1103), or a cubed-sphere face where
//!   CSLowBetaFallback replaced the flux by HLLE: same test (p_L + p_R) < lb_beta <B^2>):
//!   MD21 eq. 30, a^L = aR/(aR+aL), d^L = d^R = aR aL/(aR+aL), aR = max(0,S_R),
//!   aL = max(0,-S_L), S_L/S_R the Davis fast speeds of every solver here;
//!  hll = false (uct_hlld, rsolver hlld/lhlld): MD21 eqs. 40, 44-46,
//!   a^L = (1+nu*)/2, d^s = (nu^s - nu*) (lam*s - lam^s) chi^s/2 + (|lam*s| - nu* lam*s)/2,
//!   nu^s = (lam*s + lam^s)/(|lam*s| + |lam^s|),  nu* = (lam*R + lam*L)/(|lam*R| + |lam*L|)
//!   (nu* = 0 when |lam*R - lam*L| <= 1e-9 |S_R - S_L|, MD21 eq. 46), every speed and
//!   chi^s = B_t*s/B_t^s - 1 recomputed with the formulas and degenerate branches of
//!   hlld_mhd.hpp / lhlld_mhd.hpp (chi = 0 where B_t* = B_t; lam*L = lam*R = S_M, nu* = 0,
//!   the HLLC limit, where the solver drops the ** states).  LHLLD's low-Mach fix enters
//!   only p_T* in the momentum/energy flux, so its a/d are the HLLD ones; its Step-8
//!   resistivity (lhlld_emf_eta_pos, eta_pos here) adds ceta/2 to d^L and d^R.
//! Slots: 0 a^L, 1 d^L, 2 d^R, 3/4/5 a^L v^L + a^R v^R of the face-frame components
//! ivx, ivx+1, ivx+2 (cyclic; on the cubed sphere the gnomonic face frame the solver
//! used, see gnomonic_kernels.hpp).  Ideal or general EOS (isothermal refused).

KOKKOS_INLINE_FUNCTION
void UctFaceStore(const EOS_Data &eos, const ScrArray2D<Real> &wl,
                  const ScrArray2D<Real> &wr, const ScrArray2D<Real> &bl,
                  const ScrArray2D<Real> &br, const ScrArray2D<Real> &dl,
                  const ScrArray2D<Real> &dr, const int ivx, const int i,
                  const Real bxi, const bool hll_in, const bool eta_pos,
                  const Real lb_beta, const DvceArray5D<Real> &u,
                  const int m, const int k, const int j,
                  const DvceArray5D<Real> &w0) {
  const int ivy = IVX + ((ivx-IVX)+1)%3;
  const int ivz = IVX + ((ivx-IVX)+2)%3;
  const int iby = ((ivx-IVX) + 1)%3;
  const int ibz = ((ivx-IVX) + 2)%3;
  Real pl, pr, g1l, g1r;
  if (eos.IsGeneral()) {
    pl = dl(IDPR,i); pr = dr(IDPR,i); g1l = dl(IDG1,i); g1r = dr(IDG1,i);
  } else {
    pl = eos.IdealGasPressure(wl(IEN,i)); pr = eos.IdealGasPressure(wr(IEN,i));
    g1l = eos.gamma; g1r = eos.gamma;
  }
  const Real rl = wl(IDN,i), rr = wr(IDN,i);
  const Real vl = wl(ivx,i), vr = wr(ivx,i);
  const Real bxsq = bxi*bxi;
  const Real btl = SQR(bl(iby,i)) + SQR(bl(ibz,i));
  const Real btr = SQR(br(iby,i)) + SQR(br(ibz,i));
  const Real ptl = pl + 0.5*(bxsq + btl);
  const Real ptr = pr + 0.5*(bxsq + btr);
  Real cfl = eos.FastSpeedFromP(rl, pl, g1l, bxi, bl(iby,i), bl(ibz,i));
  Real cfr = eos.FastSpeedFromP(rr, pr, g1r, bxi, br(iby,i), br(ibz,i));
  // boris-1102 <mhd>/boris_c: the faces BorisHLLEFace took over carry an HLLE flux with
  // the limited (B -> sqrt(f) B) fast speeds: same test, same speeds, HLL weights
  // (uctaudit-1103: the one test BorisFaceTest, same |B|^2 / |v|^2 expressions)
  bool boris_face = false;
  if (eos.boris_c2inv > 0.0) {
    Real ffl, ffr;
    boris_face = BorisFaceTest(eos,
        rl, bxi*bxi + bl(iby,i)*bl(iby,i) + bl(ibz,i)*bl(ibz,i), pl,
        SQR(wl(IVX,i)) + SQR(wl(IVY,i)) + SQR(wl(IVZ,i)),
        rr, bxi*bxi + br(iby,i)*br(iby,i) + br(ibz,i)*br(ibz,i), pr,
        SQR(wr(IVX,i)) + SQR(wr(IVY,i)) + SQR(wr(IVZ,i)), ffl, ffr);
    if (boris_face) {
      const Real sfl = sqrt(ffl), sfr = sqrt(ffr);
      cfl = eos.FastSpeedFromP(rl, pl, g1l, sfl*bxi, sfl*bl(iby,i), sfl*bl(ibz,i));
      cfr = eos.FastSpeedFromP(rr, pr, g1r, sfr*bxi, sfr*br(iby,i), sfr*br(ibz,i));
    }
  }
  const Real sl = fmin(vl - cfl, vr - cfr);
  const Real sr = fmax(vl + cfl, vr + cfr);
  // CSLowBetaFallback's test, the shared CSLowBetaFires (lb_beta = 0 off the cubed sphere)
  const bool hll = hll_in || boris_face ||
                   ((lb_beta > 0.0) && CSLowBetaFires(pl, pr, bxi, bl(iby,i), bl(ibz,i),
                                                      br(iby,i), br(ibz,i), lb_beta));
  Real al, dlc, drc;
  if (sl >= 0.0) {          // F = F^L
    al = 1.0; dlc = 0.0; drc = 0.0;
  } else if (sr <= 0.0) {   // F = F^R
    al = 0.0; dlc = 0.0; drc = 0.0;
  } else if (hll) {         // MD21 eq. 30
    al  = sr/(sr - sl);
    dlc = -sr*sl/(sr - sl);
    drc = dlc;
  } else {
    const Real sdl = sl - vl, sdr = sr - vr;
    const Real sm = (sdr*(rr*vr) - sdl*(rl*vl) + (ptl - ptr))/(sdr*rr - sdl*rl);
    const Real sdml = sl - sm, sdmr = sr - sm;
    const Real ptst = 0.5*((ptr + rr*sdr*(sm - vr)) + (ptl + rl*sdl*(sm - vl)));
    const Real denl = rl*sdl*sdml - bxsq;
    const Real denr = rr*sdr*sdmr - bxsq;
    const Real chil = (fabs(denl) < (HLLD_SMALL_NUMBER)*ptst) ? 0.0 :
                      rl*sdl*(sm - vl)/denl;
    const Real chir = (fabs(denr) < (HLLD_SMALL_NUMBER)*ptst) ? 0.0 :
                      rr*sdr*(sm - vr)/denr;
    Real lsl, lsr, nus;
    if (0.5*bxsq < eos.hlld_bx_zero_tol*ptst) {   // solver: ** = *, HLLC limit
      lsl = sm; lsr = sm; nus = 0.0;
    } else {
      lsl = sm - fabs(bxi)/sqrt(rl*sdl/sdml);
      lsr = sm + fabs(bxi)/sqrt(rr*sdr/sdmr);
      const Real den = fabs(lsr) + fabs(lsl);
      nus = (fabs(lsr - lsl) > 1.0e-9*fabs(sr - sl) && den > 0.0) ?
            (lsr + lsl)/den : 0.0;
    }
    const Real dnl = fabs(lsl) + fabs(sl);
    const Real dnr = fabs(lsr) + fabs(sr);
    const Real nul = (dnl > 0.0) ? (lsl + sl)/dnl : 0.0;
    const Real nur = (dnr > 0.0) ? (lsr + sr)/dnr : 0.0;
    al  = 0.5*(1.0 + nus);
    dlc = 0.5*(nul - nus)*(lsl - sl)*chil + 0.5*(fabs(lsl) - nus*lsl);
    drc = 0.5*(nur - nus)*(lsr - sr)*chir + 0.5*(fabs(lsr) - nus*lsr);
    // alfmach-1102 <mhd>/lhlld_alfmach_rot (LHLLD only, refused for other solvers): the
    // solver's B_t flux is F_ref + psi_A (F_HLLD - F_ref), F_ref = (F*L + F*R)/2 =
    // (F^L + F^R)/2 + (S_L chi^L B^L + S_R chi^R B^R)/2, i.e. (a, d^L, d^R)_ref =
    // (1/2, S_L chi^L/2, -S_R chi^R/2); blend the coefficients with the same psi_A
    if (eos.alfmach_rot && !(0.5*bxsq < eos.hlld_bx_zero_tol*ptst) && lsl < 0.0 &&
        lsr > 0.0) {
      Real usql = SQR(vl) + (SQR(wl(ivy,i)) + SQR(wl(ivz,i)));
      Real usqr = SQR(vr) + (SQR(wr(ivy,i)) + SQR(wr(ivz,i)));
      if (eos.chi_normal) { usql = SQR(vl); usqr = SQR(vr); }
      const Real caxs = fabs(bxi)*fmax(1.0/sqrt(rl*sdl/sdml), 1.0/sqrt(rr*sdr/sdmr));
      const Real psia = LhlldAlfRotPsi(w0, m, k, j, i, ivx, ivy, ivz, fmax(usql, usqr),
                                       caxs, eos.cb_sensor, eos.cb_cont,
                                       eos.alfmach_rot_min);
      al  = 0.5 + psia*(al - 0.5);
      dlc = 0.5*sl*chil + psia*(dlc - 0.5*sl*chil);
      drc = -0.5*sr*chir + psia*(drc + 0.5*sr*chir);
    }
    if (eta_pos) {   // LHLLD Step 8 (lhlld_emf_eta_pos)
      const Real ceta = fmax(0.0, -(fabs(sm) + 0.5*(vr - vl)));
      dlc += 0.5*ceta;
      drc += 0.5*ceta;
    }
  }
  u(m,0,k,j,i) = al;
  u(m,1,k,j,i) = dlc;
  u(m,2,k,j,i) = drc;
  u(m,3,k,j,i) = al*vl + (1.0 - al)*vr;
  u(m,4,k,j,i) = al*wl(ivy,i) + (1.0 - al)*wr(ivy,i);
  u(m,5,k,j,i) = al*wl(ivz,i) + (1.0 - al)*wr(ivz,i);
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::CalculateFlux
//! \brief Calculate fluxes of conserved variables, and face-centered area-averaged EMFs
//! for evolution of magnetic field
//! Note this function is templated over RS for better performance on GPUs.

template <MHD_RSolver rsolver_method_>
void MHD::CalculateFluxes(Driver *pdriver, int stage) {
  RegionIndcs &indcs_ = pmy_pack->pmesh->mb_indcs;
  int is = indcs_.is, ie = indcs_.ie;
  int js = indcs_.js, je = indcs_.je;
  int ks = indcs_.ks, ke = indcs_.ke;
  int ncells1 = indcs_.nx1 + 2*(indcs_.ng);

  int &nmhd_ = nmhd;
  int nvars = nmhd + nscalars;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const auto recon_method_ = recon_method;
  bool extrema = false;
  if (recon_method == ReconstructionMethod::ppmx) {
    extrema = true;
  }

  auto &eos_ = peos->eos_data;
  auto &size_ = pmy_pack->pmb->mb_size;
  auto &coord_ = pmy_pack->pcoord->coord_data;
  auto &w0_ = w0;
  auto &b0_ = bcc0;

  // For a general EOS, pressure and Gamma_1 were evaluated once per cell in ConsToPrim
  // and are reconstructed to interfaces alongside the primitives, so the Riemann solvers
  // consume pressure directly instead of calling the (expensive) EOS. nder is 0 for an
  // ideal gas, leaving the original code path untouched.
  auto &wder_ = wder;
  const int nder = eos_.IsGeneral() ? NDERIVED : 0;
    
  auto &pwb_ = pwb;
  auto pfacewb_x1f = pfacewb.x1f;
  auto w0facewb_x1f = w0facewb.x1f;
  auto pfacewb_x2f = pfacewb.x2f;
  auto w0facewb_x2f = w0facewb.x2f;
  auto pfacewb_x3f = pfacewb.x3f;
  auto w0facewb_x3f = w0facewb.x3f;
  auto &phicc0_ = phicc0;
  auto &wbq0_ = wbq0;     // the per-cell well-balanced background, built just below
  // the x1 well-balanced walk uses the EFFECTIVE potential (MHD::phicc_wb); it IS
  // phicc0 / phi0.x1f unless a problem generator enabled its own (bitwise default)
  auto &phiwbcc_ = phicc_wb;
  auto phiwb_x1f_ = phi_wb_x1f;
  auto phi0_x2f = phi0.x2f;
  auto phi0_x3f = phi0.x3f;
    
  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  const Real cs_lb_beta = cs_lowbeta_fallback;
  // boris-1102 <mhd>/boris_c: limited-speed HLLE on the faces the limiter acts on
  const bool boris_on_ = (peos->eos_data.boris_c2inv > 0.0);
  const bool polar_hlle_x1_ = (polar_hlle_rows & 1) != 0;
  const bool polar_hlle_x2_ = (polar_hlle_rows & 2) != 0;
  const bool polar_hlle_x3_ = (polar_hlle_rows & 4) != 0;
  const auto wb_option_ = wb_option;
  const bool use_wb_rho_ = use_wb_rho;
  const bool use_wb_x1_ = use_wb_x1;
  const bool use_wb_x2_ = use_wb_x2;
  const bool wb_all_ = wb_all;   // wball-1104
  auto wbq2_ = wbq2;
  auto wbq3_ = wbq3;
  auto wbw_ = wbw;
  const bool use_wellbalance_dynamic_ = use_wellbalance_dynamic;
  const Real wb_rmax_ = wb_rmax;
  const Real wb_phimax_ = wb_phimax;   // dhjwball-1106
  const Real wb_rmin_ = wb_rmin;
  const bool use_wellbalance_static_reconst_perturb_ =
      use_wellbalance_static_reconst_perturb;
  // The gnomonic rotations need one (sin, cos) pair each, so hand the kernels that pair
  // rather than the Coordinates object: reaching them through `pmy_pack->pcoord->` made
  // the lambda capture `this` and dereference a host pointer on the device. ONE pair per
  // sweep, never all six -- see the note in gnomonic_kernels.hpp.
  auto gtrig_cell = pmy_pack->pcoord->GnomonicTrigCell();
  auto gtrig_xi   = pmy_pack->pcoord->GnomonicTrigFaceXi();
  auto gtrig_eta  = pmy_pack->pcoord->GnomonicTrigFaceEta();
  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  // A STRETCHED radial grid on the cubed sphere takes the same position-aware x1
  // reconstruction as spherical polar; on a uniform radial grid nothing changes.
  // (x1v is the volume centroid on both spherical grids, so the plain uniform stencil
  // is off-centre even without a stretch; `reconstruct` governs the ANGULAR sweeps)
  const bool str_r1_ = pmy_pack->pmesh->use_cubed_sphere;
  const bool wall_ix1_ = wall_closed_ix1;
  const bool diode_ox1_ = diode_ox1;   // topbc-1101 (mhd.hpp)
  const bool wall_sf_ = wall_closed_ix1 && wall_stress_free;
  // the dhj conducting wall is a closed wall on every geometry too (legacy mirror)
  const bool wall_cd_ = wall_closed_ix1 && wall_conducting_ix1;
  // wallwb-1105 <mhd>/wall_wb_onesided (mhd.hpp, reconstruct/wb_wall.hpp, shared with
  // Hydro::Fluxes): every geometry, on the x1 dynamic WB paths that read wbq0 (PLM)
  const bool curv1_ = pmy_pack->pmesh->use_spherical_polar || str_r1_;
  const bool wwb_ = wall_wb_onesided && use_wellbalance_dynamic && use_wb_x1 &&
                    (wb_all || recon_method_x1 == ReconstructionMethod::plm || !curv1_);
  // <mhd>/reconstruct_x1: separate knob for the radial sweep on cs/sp (see mhd.hpp).
  const auto recon_method_x1_ = recon_method_x1;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool pquad_ = pmy_pack->pmesh->use_polar_quadratic_recon;
  const int px3_ = pmy_pack->pmesh->polar_x3_shift;

  //--------------------------------------------------------------------------------------
  // i-direction

  size_t scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
                     ScrArray2D<Real>::shmem_size(3, ncells1) +
                     ScrArray2D<Real>::shmem_size(nder, ncells1)) * 2;
  int scr_level = 0;
  auto &flx1_ = uflx.x1f;
  auto &e31_ = e3x1;
  auto &e21_ = e2x1;
  auto &bx_ = b0.x1f;
  // <mhd>/emf_scheme = uct | uct_hlld: UCT corner-EMF face store (UctFaceStore)
  const bool uct_on_ = emf_uct || emf_uct_hlld;
  const bool uct_hll_ = emf_uct || (rsolver_method_ == MHD_RSolver::hlle);
  const bool uct_eta_ = (rsolver_method_ == MHD_RSolver::lhlld) &&
                        peos->eos_data.lhlld_emf_eta_pos;
  const Real uct_lb_ = cs_lowbeta_fallback;   // share-1104 G8: every grid (0 off cs by default)
  auto &uct1_ = uct1;
  auto &uct2_ = uct2;
  auto &uct3_ = uct3;
  // ucthshare-1102 <mhd>/cs_seam_mhd = faces | shared: the neighbour's own seam face
  // states (MHD::CsSeamFaceStates, mhd_cs_seam.cpp), substituted on the ghost side of
  // every panel-seam face (rows is-1..ie+1) after the gnomonic face rotation
  const bool cs_shr_ = (cs_seam_mhd > 0);
  if (wb_all) BuildWBCacheAll();   // wball-1104 (before the seam face states)
  if (wb_all && wb_poison) {   // dhjwball-1106 TEST ONLY
    hydro::WbPoison(pmy_pack, wbq2, phicc_wb, wb_phimax);
    if (wbq3.extent_int(2) > 1) hydro::WbPoison(pmy_pack, wbq3, phicc_wb, wb_phimax);
  }
  if (cs_shr_) {CsSeamFaceStates();}
  DvceArray5D<Real> sfs2_, sfs3_;
  if (cs_shr_) {sfs2_ = pseamfs->x2f; sfs3_ = pseamfs->x3f;}
  const int sfob_ = nvars + NDERIVED;   // record slots of (B.r, B.e_xi, B.e_eta)
    
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;

  // set the loop limits for 1D/2D/3D problems
  int jl,ju,kl,ku;
  if (pmy_pack->pmesh->one_d) {
    jl = js, ju = je, kl = ks, ku = ke;
  } else if (pmy_pack->pmesh->two_d) {
    jl = js-1, ju = je+1, kl = ks, ku = ke;
  } else {
    jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
  }
  int il = is, iu = ie+1;
  if (use_fofc) { il = is-1, iu = ie+2; }
  // TRANSVERSE WIDENING.  FOFC replaces the fluxes on all six faces of every flagged
  // cell, and cells are flagged over [js-1,je+1] x [ks-1,ke+1] (MHD::FOFC), so the x1
  // fluxes are needed one row beyond the active zone in x2 and x3 as well.  Unlike
  // hydro, MHD already computes this sweep over exactly those transverse ranges
  // UNCONDITIONALLY -- the corner-EMF stencil needs them -- so no fofc widening is
  // required here (nor in the x2/x3 sweeps, which run over is-1,ie+1 in x1 and
  // js-1,je+1 in x2 for the same reason).  There is no MHD counterpart of the
  // spherical-polar face-average (sp_favg) widening, so there is nothing to union with.

  // the well-balanced background of every cell in this sweep, walked once

  if (use_wellbalance_dynamic && use_wb_x1) {
    const int ncyc = pmy_pack->pmesh->ncycle;
    if (wb_cache_every <= 0 || !wb_cache_built ||
        (stage == 1 && (ncyc % wb_cache_every) == 0)) {
      BuildWBCache(jl, ju, kl, ku);
      if (wb_poison) hydro::WbPoison(pmy_pack, wbq0, phicc_wb, wb_phimax);   // dhjwball-1106 TEST
      wb_cache_built = true;
    }
  }


  scr_level = TeamScratchLevel(scr_size, 0);
  par_for_outer("mhd_flux1",DevExeSpace(), scr_size, scr_level, 0, nmb1, kl, ku, jl, ju,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray2D<Real> wl(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> wr(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> bl(member.team_scratch(scr_level), 3, ncells1);
    ScrArray2D<Real> br(member.team_scratch(scr_level), 3, ncells1);
    ScrArray2D<Real> dl(member.team_scratch(scr_level), nder, ncells1);
    ScrArray2D<Real> dr(member.team_scratch(scr_level), nder, ncells1);

    if (use_spherical_polar || str_r1_) {
      if (recon_method_x1_ == ReconstructionMethod::plm) {
        GridPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
                              use_wellbalance_dynamic_,
                              use_wb_x1_, wb_rmax_, wb_rmin_,
                              m, k, j, il-1, iu, w0_, x1v_, x1f_, phiwbcc_,
                              phiwb_x1f_, wbq0_, true, wl, wr, wb_all_, wb_phimax_);
        GridPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
                              use_wellbalance_dynamic_,
                              use_wb_x1_, wb_rmax_, wb_rmin_,
                              m, k, j, il-1, iu, b0_, x1v_, x1f_, phiwbcc_,
                              phiwb_x1f_, wbq0_, false, bl, br);
      } else {
        // ppm4/wenoz on the RADIAL sweep: the published (Mignone 2014) curvilinear
        // volume-coordinate reconstruction -- see reconstruct/mignone_curvilinear.hpp
        // for the full derivation and equation references. ppmx is refused at setup
        // (no curvilinear Colella-Sekora limiter in that paper). Setup also refused
        // wellbalance_dynamic+wb_x1 with this choice, so there is no background to
        // subtract here.
        switch (recon_method_x1_) {
          case ReconstructionMethod::ppm4:
            MignonePPM4X1(member,eos_,true, m,k,j,il-1,iu,w0_,x1v_,x1f_,wl,wr);
            MignonePPM4X1(member,eos_,false,m,k,j,il-1,iu,b0_,x1v_,x1f_,bl,br);
            break;
          case ReconstructionMethod::wenoz:
            MignoneWENO3X1(member,eos_,true, m,k,j,il-1,iu,w0_,x1v_,x1f_,wl,wr);
            MignoneWENO3X1(member,eos_,false,m,k,j,il-1,iu,b0_,x1v_,x1f_,bl,br);
            break;
          default:
            break;
        }
      }
    } else {
        
    if (use_wellbalance_dynamic_ && use_wb_x1_)
    {
      if (wb_all_) {   // wball-1104: x1 through the shared routine (w = 1)
        hydro::WbAllPLM(member, use_wb_rho_, 1, m, k, j, il-1, iu, w0_, wbq0_, wbw_, wl, wr);
      } else {
      WbLocalPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
          m, k, j, il-1, iu, w0_, phiwbcc_, phiwb_x1f_, wbq0_, wl, wr);
      }
      if (eos_.plm_c1) {
        PiecewiseLinearC1X1(member, m, k, j, il-1, iu, b0_, bl, br);
      } else {
        PiecewiseLinearX1(member, m, k, j, il-1, iu, b0_, bl, br);
      }
    } else {

    // Reconstruct qR[i] and qL[i+1], for both W and Bcc
    switch (recon_method_) {
      case ReconstructionMethod::dc:
        DonorCellX1(member, m, k, j, il-1, iu, w0_, wl, wr);
        DonorCellX1(member, m, k, j, il-1, iu, b0_, bl, br);
        break;
      case ReconstructionMethod::plm:
        if (eos_.plm_c1) {
          PiecewiseLinearC1X1(member, m, k, j, il-1, iu, w0_, wl, wr);
        } else {
          PiecewiseLinearX1(member, m, k, j, il-1, iu, w0_, wl, wr);
        }
        if (eos_.plm_c1) {
          PiecewiseLinearC1X1(member, m, k, j, il-1, iu, b0_, bl, br);
        } else {
          PiecewiseLinearX1(member, m, k, j, il-1, iu, b0_, bl, br);
        }
        break;
      case ReconstructionMethod::ppm4:
      case ReconstructionMethod::ppmx:
        PiecewiseParabolicX1(member,eos_,extrema,true,  m, k, j, il-1, iu, w0_, wl, wr);
        PiecewiseParabolicX1(member,eos_,extrema,false, m, k, j, il-1, iu, b0_, bl, br);
        break;
      case ReconstructionMethod::wenoz:
        WENOZX1(member, eos_, true,  m, k, j, il-1, iu, w0_, wl, wr);
        WENOZX1(member, eos_, false, m, k, j, il-1, iu, b0_, bl, br);
        break;
      default:
        break;
    }
        
    }
    }
      
    // Reconstruct the derived thermodynamic variables (general EOS only). apply_floors is
    // false: that logic keys off n==IDN/n==IEN, whose values collide with DerivedIndex.
    if (nder > 0) {
      // the well-balanced scheme must treat pressure in deviation form too: a general EOS
      // hands the Riemann solver this reconstructed pressure rather than recomputing it
      // from the reconstructed (d,e), so leaving it to plain PLM would put the entire
      // hydrostatic gradient back into the solver's pressure and unbalance the scheme.
      if ((use_spherical_polar || str_r1_) &&
          recon_method_x1_ == ReconstructionMethod::plm) {
        GridPiecewiseLinearDerX1(member, eos_, wb_option_, use_wellbalance_dynamic_,
                                 use_wb_x1_, wb_rmax_, wb_rmin_,
                                 use_wellbalance_static_reconst_perturb_,
                                 pwb_, pfacewb_x1f,
                                 m, k, j, il-1, iu, w0_, wder_,
                                 x1v_, x1f_, phiwbcc_, phiwb_x1f_, wbq0_, dl, dr, wb_phimax_);
      } else if (use_spherical_polar || str_r1_) {
        // reconstruct_x1 != plm: same Mignone curvilinear reconstruction as w0_/b0_
        // above, applied to the derived thermodynamic variables (no WB background,
        // refused at setup for this combination).
        switch (recon_method_x1_) {
          case ReconstructionMethod::ppm4:
            MignonePPM4X1(member,eos_,false,m,k,j,il-1,iu,wder_,x1v_,x1f_,dl,dr);
            break;
          case ReconstructionMethod::wenoz:
            MignoneWENO3X1(member,eos_,false,m,k,j,il-1,iu,wder_,x1v_,x1f_,dl,dr);
            break;
          default:
            break;
        }
      } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX1(member, m, k, j, il-1, iu,
                                     pwb_, pfacewb_x1f,
                                    wder_, dl, dr);
      } else if (use_wellbalance_dynamic_ && use_wb_x1_) {
        WbPiecewiseLinearDerX1(member, eos_, wb_option_,
            m, k, j, il-1, iu, w0_, wder_,
                               phiwbcc_, phiwb_x1f_, wbq0_, dl, dr);
      } else {
      switch (recon_method_) {
        case ReconstructionMethod::dc:
          DonorCellX1(member, m, k, j, il-1, iu, wder_, dl, dr);
          break;
        case ReconstructionMethod::plm:
          if (eos_.plm_c1) {
            PiecewiseLinearC1X1(member, m, k, j, il-1, iu, wder_, dl, dr);
          } else {
            PiecewiseLinearX1(member, m, k, j, il-1, iu, wder_, dl, dr);
          }
          break;
        case ReconstructionMethod::ppm4:
        case ReconstructionMethod::ppmx:
          PiecewiseParabolicX1(member,eos_,extrema,false,m,k,j,il-1,iu,wder_,dl,dr);
          break;
        case ReconstructionMethod::wenoz:
          WENOZX1(member, eos_, false, m, k, j, il-1, iu, wder_, dl, dr);
          break;
        default:
          break;
      }
      }
      // reconstruction can undershoot; keep the interface pressure positive
      // The x1 reconstruction writes dl(n,i+1) from the thread that owns i, so the
      // floor below reads a slot ANOTHER thread wrote.  Successive par_for_inner
      // loops are not implicitly synchronised, so without this barrier whether the
      // floor sees the reconstructed value or a stale one is scheduling-dependent,
      // and the run is not reproducible.  x2/x3 need no barrier: those
      // reconstructions write index i from thread i.
      member.team_barrier();
      par_for_inner(member, il, iu, [&](const int i) {
        dl(IDPR,i) = fmax(dl(IDPR,i), eos_.pfloor);
        dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
      });
    }

      if (use_wellbalance_static_reconst_perturb_) {
        AddWbPrimFaceX1(member,m,k,j,il-1,iu,w0facewb_x1f,wl,wr);
      }

      // CUBED SPHERE. The gnomonic tangent basis is not orthogonal, so the primitive
      // velocity (contravariant) has to be put into a locally ORTHONORMAL frame before
      // the Riemann solver sees it, and the returned momentum flux has to be rotated
      // back and its index lowered. hydro_fluxes.cpp has done this since the grid was
      // added; the MHD path never did, which left the fluxes in a different frame from
      // the geometric source terms that SrcTermsGnomonicEquiangle adds to them.
      if (use_cubed_sphere) {
        // the field needs no rotation in this sweep: bcc is already stored in this
        // frame. See the note at GnomonicEquiangleFaceBX2.
        GnomonicEquianglePrimFaceX1(gtrig_cell,member,m,k,j,il-1,iu,wl,wr);
      }

    // A REFLECTING radial wall passes exactly zero mass flux only if the two states the
    // Riemann solver sees at the wall face are mirror images.  The position-aware x1
    // reconstruction (spherical polar and, since the radial unification, the cubed
    // sphere) is not mirror-symmetric about the wall -- the ghost centroids do not mirror
    // the interior ones -- so the ghost-side state differed from the interior one at
    // O(dx^2) and a strong blast leaked 2-3e-4 of the mass through a "closed" wall.  Make
    // the ghost-side state the exact mirror of the interior-side one: any solver then
    // returns a zero mass flux (S_M = 0 by symmetry) to round-off.
    // A Cartesian mesh takes the user-wall mirror only under the stress-free or the
    // conducting wall (the legacy closed wall was never applied there; reflect is
    // untouched on Cartesian).
    if (use_spherical_polar || str_r1_ || wall_sf_ || wall_cd_ || wwb_) {
      // wallwb-1105: under <mhd>/wall_wb_onesided a Cartesian reflect wall takes the
      // mirror as the curvilinear one (key off: unchanged)
      const bool curv_ = use_spherical_polar || str_r1_ || wwb_;
      member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        // (user inner boundary under wall_closed_ix1 too -- see hydro_fluxes.cpp)
        const BoundaryFlag bix1 = mb_bcs.d_view(m,BoundaryFace::inner_x1);
        const BoundaryFlag box1 = mb_bcs.d_view(m,BoundaryFace::outer_x1);
        // wallwb-1105: the wall cells' one-sided WB deviation (reconstruct/wb_wall.hpp)
        if (wwb_) {
          auto inr = [&](const int i) {
            return !curv1_ || (!(wb_rmax_ > 0.0 && x1v_(m,i) > wb_rmax_) &&
                               !(wb_rmin_ > 0.0 && x1v_(m,i) < wb_rmin_) &&
                               !(phiwbcc_(m,k,j,i) > wb_phimax_));
          };
          if ((bix1 == BoundaryFlag::reflect || (wall_ix1_ && bix1 == BoundaryFlag::user))
              && inr(is)) {
            wbwall::WallCellX1(true, m, k, j, is, w0_, wder_, wbq0_, x1v_, x1f_,
                               use_wb_rho_, nder, wl, wr, dl, dr);
          }
          if (box1 == BoundaryFlag::reflect && inr(ie)) {
            wbwall::WallCellX1(false, m, k, j, ie, w0_, wder_, wbq0_, x1v_, x1f_,
                               use_wb_rho_, nder, wl, wr, dl, dr);
          }
        }
        if ((curv_ && bix1 == BoundaryFlag::reflect) ||
            (wall_ix1_ && bix1 == BoundaryFlag::user)) {
          for (int n=0; n<nvars; ++n) { wl(n,is) = wr(n,is); }
          wl(IVX,is) = -wr(IVX,is);
          bl(1,is) = br(1,is);  bl(2,is) = br(2,is);
          if (bix1 == BoundaryFlag::user || wwb_) {  // reflect stays bitwise (key off)
            for (int n=0; n<nder; ++n) {
              dl(n,is) = dr(n,is);
            }
          }
          if (bix1 == BoundaryFlag::user && wall_sf_) {  // stress-free: B_t odd
            bl(1,is) = -br(1,is);  bl(2,is) = -br(2,is);
          }
        }
        if (curv_ && box1 == BoundaryFlag::reflect) {
          for (int n=0; n<nvars; ++n) { wr(n,ie+1) = wl(n,ie+1); }
          wr(IVX,ie+1) = -wl(IVX,ie+1);
          br(1,ie+1) = bl(1,ie+1);  br(2,ie+1) = bl(2,ie+1);
          if (wwb_) {
            for (int n=0; n<nder; ++n) { dr(n,ie+1) = dl(n,ie+1); }
          }
        }
      });
    }
    // topbc-1101 diode_ox1: a one-way outer radial face (see mhd.hpp)
    if (diode_ox1_) {
      member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user &&
            wl(IVX,ie+1) < 0.0) {
          for (int n=0; n<nvars; ++n) { wr(n,ie+1) = wl(n,ie+1); }
          wr(IVX,ie+1) = -wl(IVX,ie+1);
          br(1,ie+1) = bl(1,ie+1);  br(2,ie+1) = bl(2,ie+1);
          for (int n=0; n<nder; ++n) { dr(n,ie+1) = dl(n,ie+1); }
        }
      });
    }
    // Sync all threads in the team so that scratch memory is consistent
    member.team_barrier();

    // compute fluxes over [is,ie+1].  MHD RS also computes electric fields, where
    // (IBY) component of flx = E_{z} = -(v x B)_{z} = -(v1*b2 - v2*b1)
    // (IBZ) component of flx = E_{y} = -(v x B)_{y} =  (v1*b3 - v3*b1)
    // NOTE(@pdmullen): Capture variables prior to if constexpr.  Required for cuda 11.6+.
    // dhjcuda2-1113: as const references to the captured copies, not new copies: on
    // CUDA the body sits in device memory, and a reference reads each field where it
    // is used, while a copy holds every View pointer and stride (and EOS_Data) in
    // registers through the whole solver (A100 ptxas: 255 registers, spills).  Same
    // values: bitwise.
    const auto &eos = eos_;
    const auto &indcs = indcs_;
    const auto &size = size_;
    const auto &coord = coord_;
    const auto &bx = bx_;
    const auto &flx1 = flx1_;
    const auto &e31 = e31_;
    const auto &e21 = e21_;
    const bool do_pole = polar_hlle_x1_ &&
      ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) ||
       (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je));
    if constexpr (rsolver_method_ == MHD_RSolver::advect) {
      Advect(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::llf) {
      LLF(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::hlle) {
      HLLE(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::hlld) {
      if (do_pole) {
        HLLE(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,
             bx,flx1,e31,e21);
      } else {
      HLLD(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,
           bx,flx1,e31,e21);
      }
    } else if constexpr (rsolver_method_ == MHD_RSolver::lhlld) {
      if (do_pole) {
        HLLE(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,
             bx,flx1,e31,e21);
      } else {
        LHLLD(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,
              bx,flx1,e31,e21,w0_);
      }
    } else if constexpr (rsolver_method_ == MHD_RSolver::llf_sr) {
      LLF_SR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_sr) {
      HLLE_SR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::llf_gr) {
      LLF_GR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_gr) {
      HLLE_GR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    }

    // share-1104 G8: the per-face low-beta HLLD -> HLLE fallback on every grid (default
    // 0 = off except on the cubed sphere, so bitwise unless named); same order as on cs
    if (!use_cubed_sphere && cs_lb_beta > 0.0) {
      member.team_barrier();
      CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,il,iu,IVX,cs_lb_beta,
                        wl,wr,bl,br,dl,dr,bx,flx1,e31,e21);
      member.team_barrier();
    }
    if (boris_on_ && !use_cubed_sphere) {
      member.team_barrier();
      BorisHLLEFace(member,eos,eos.IsGeneral(),m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,
                    bx,flx1,e31,e21);
    }
    if (use_cubed_sphere) {
      // low-beta dissipation fallback, BEFORE the rotation into covariant slots
      if (cs_lb_beta > 0.0) {
        CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,il,iu,IVX,cs_lb_beta,
                     wl,wr,bl,br,dl,dr,bx,flx1,e31,e21);
        member.team_barrier();
      }
      if (boris_on_) {
        member.team_barrier();
        BorisHLLEFace(member,eos,eos.IsGeneral(),m,k,j,il,iu,IVX,wl,wr,bl,br,dl,dr,
                      bx,flx1,e31,e21);
        member.team_barrier();
      }
      GnomonicEquiangleFluxX1(gtrig_cell,member,m,k,j,il,iu,flx1);
      GnomonicEquiangleEmfX1(gtrig_cell,member,m,k,j,il,iu,e31,e21);
    }
    member.team_barrier();
    if (uct_on_) {   // UCT corner EMF inputs (x1 faces)
      par_for_inner(member, il, iu, [&](const int i) {
        UctFaceStore(eos, wl, wr, bl, br, dl, dr, IVX, i, bx(m,k,j,i),
                     uct_hll_ || do_pole, uct_eta_,
                     uct_lb_, uct1_, m, k, j, w0_);
      });
      member.team_barrier();
    }

    // calculate fluxes of scalars (if any)
    if (nvars > nmhd_) {
      for (int n=nmhd_; n<nvars; ++n) {
        par_for_inner(member, is, ie+1, [&](const int i) {
          if (flx1_(m,IDN,k,j,i) >= 0.0) {
            flx1_(m,n,k,j,i) = flx1_(m,IDN,k,j,i)*wl(n,i);
          } else {
            flx1_(m,n,k,j,i) = flx1_(m,IDN,k,j,i)*wr(n,i);
          }
        });
      }
    }
  });

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmy_pack->pmesh->multi_d) {
    scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
                ScrArray2D<Real>::shmem_size(3, ncells1) +
                ScrArray2D<Real>::shmem_size(nder, ncells1)) * 3;
    auto &flx2_ = uflx.x2f;
    auto &by_ = b0.x2f;
    auto &e12_ = e1x2;
    auto &e32_ = e3x2;
      
    auto &x2v_ = pmy_pack->pcoord->x2v_rec;  // == x2v unless sp_x2_periodic_image
    auto &x2f_ = pmy_pack->pcoord->xx2f;

    // set the loop limits for 2D/3D problems
    if (pmy_pack->pmesh->two_d) {
      kl = ks, ku = ke;
    } else { // 3D
      kl = ks-1, ku = ke+1;
    }
    jl = js-1, ju = je+1;
    if (use_fofc) { jl = js-2, ju = je+2; }

    scr_level = TeamScratchLevel(scr_size, 0);
    par_for_outer("mhd_flux2",DevExeSpace(),scr_size,scr_level,0,nmb1, kl, ku,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k) {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr4(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr5(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr6(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> dscr1(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr2(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr3(member.team_scratch(scr_level), nder, ncells1);

      for (int j=jl; j<=ju; ++j) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_jp1 = scr2;
        auto wr     = scr3;
        auto bl     = scr4;
        auto bl_jp1 = scr5;
        auto br     = scr6;
        auto dl     = dscr1;
        auto dl_jp1 = dscr2;
        auto dr     = dscr3;
        if ((j%2) == 0) {
          wl     = scr2;
          wl_jp1 = scr1;
          bl     = scr5;
          bl_jp1 = scr4;
          dl     = dscr2;
          dl_jp1 = dscr1;
        }
          
        const bool pq = pquad_ &&
          ((j == js && mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar) ||
           (j == je && mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar));
        if (use_spherical_polar) {
          GridPiecewiseLinearX2(member, m, k, j, is-1, ie+1, w0_, x2v_, x2f_, wl_jp1, wr,
                                pq);
          GridPiecewiseLinearX2(member, m, k, j, is-1, ie+1, b0_, x2v_, x2f_, bl_jp1, br,
                                pq);
          if (wb_all_) {
            hydro::WbAllBlendGrid(member, use_wb_rho_, 2, false, m, k, j, is-1, ie+1, w0_,
                                  x2v_, x2f_, wbq2_, wbw_, wl_jp1, wr);
          }
        } else {
            
        if (use_wellbalance_dynamic_ && use_wb_x2_)
        {
          // is-1,ie+1, not il,iu -- the b0_ line below, the spherical-polar branch
          // above, the whole default switch, and the Riemann solver all use
          // is-1,ie+1.  With il,iu (the x1 sweep's limits, never reset for x2/x3)
          // this leaves wl_jp1/wr unwritten at i = is-1, which the solver reads.
          WbLocalPiecewiseLinearX2(member, eos_, wb_option_, use_wb_rho_,
              m, k, j, is-1, ie+1, w0_,
                                   phicc0_, phi0_x2f, wl_jp1, wr);
          if (eos_.plm_c1) {
            PiecewiseLinearC1X2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
          } else {
            PiecewiseLinearX2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
          }
        } else {

        // Reconstruct qR[j] and qL[j+1], for both W and Bcc
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX2(member, m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            DonorCellX2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
            break;
          case ReconstructionMethod::plm:
            if (eos_.plm_c1) {
              PiecewiseLinearC1X2(member, m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            } else {
              PiecewiseLinearX2(member, m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            }
            if (eos_.plm_c1) {
              PiecewiseLinearC1X2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
            } else {
              PiecewiseLinearX2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
            }
            break;
          case ReconstructionMethod::ppm4:
          case ReconstructionMethod::ppmx:
            PiecewiseParabolicX2(member,eos_,extrema,true, m,k,j,is-1,ie+1,w0_,wl_jp1,wr);
            PiecewiseParabolicX2(member,eos_,extrema,false,m,k,j,is-1,ie+1,b0_,bl_jp1,br);
            break;
          case ReconstructionMethod::wenoz:
            WENOZX2(member, eos_, true,  m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            WENOZX2(member, eos_, false, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
            break;
          default:
            break;
        }
            
        }
        if (wb_all_) {   // wball-1104: the WB deviation blended into the plain faces
          hydro::WbAllBlendUni(member, use_wb_rho_, 2, false, m, k, j, is-1, ie+1, w0_,
                               wbq2_, wbw_, wl_jp1, wr);
        }
        }
          
        // Reconstruct the derived thermodynamic variables (general EOS only). NOTE the
        // range is is-1,ie+1 to match what the Riemann solver below actually reads: il,iu
        // are only ever set for the i-direction sweep and are not reset here.
        if (nder > 0) {
          if (use_spherical_polar) {
            GridPiecewiseLinearX2(member, m, k, j, is-1, ie+1, wder_, x2v_, x2f_,
                                  dl_jp1, dr,
                                  pq);
          } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX2(member, m, k, j, is-1, ie+1,
                                     pwb_, pfacewb_x2f,
                                    wder_, dl_jp1, dr);
      } else if (use_wellbalance_dynamic_ && use_wb_x2_) {
            WbPiecewiseLinearDerX2(member, eos_, wb_option_,
                m, k, j, is-1, ie+1, w0_, wder_,
                                   phicc0_, phi0_x2f, dl_jp1, dr);
          } else {
          switch (recon_method_) {
            case ReconstructionMethod::dc:
              DonorCellX2(member, m, k, j, is-1, ie+1, wder_, dl_jp1, dr);
              break;
            case ReconstructionMethod::plm:
              if (eos_.plm_c1) {
                PiecewiseLinearC1X2(member, m, k, j, is-1, ie+1, wder_, dl_jp1, dr);
              } else {
                PiecewiseLinearX2(member, m, k, j, is-1, ie+1, wder_, dl_jp1, dr);
              }
              break;
            case ReconstructionMethod::ppm4:
            case ReconstructionMethod::ppmx:
              PiecewiseParabolicX2(member,eos_,extrema,false,m,k,j,is-1,ie+1,
                                   wder_,dl_jp1,dr);
              break;
            case ReconstructionMethod::wenoz:
              WENOZX2(member, eos_, false, m, k, j, is-1, ie+1, wder_, dl_jp1, dr);
              break;
            default:
              break;
          }
          }
          if (wb_all_) {
            if (use_spherical_polar) {
              hydro::WbAllBlendGrid(member, false, 2, true, m, k, j, is-1, ie+1, wder_,
                                    x2v_, x2f_, wbq2_, wbw_, dl_jp1, dr);
            } else {
              hydro::WbAllBlendUni(member, false, 2, true, m, k, j, is-1, ie+1, wder_,
                                   wbq2_, wbw_, dl_jp1, dr);
            }
          }
          par_for_inner(member, is-1, ie+1, [&](const int i) {
            dl_jp1(IDPR,i) = fmax(dl_jp1(IDPR,i), eos_.pfloor);
            dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
          });
        }

          if (use_wellbalance_static_reconst_perturb_) {
            // is-1,ie+1 -- NOT il,iu.  il/iu are the x1 sweep's limits (set once at
            // the top of this function) and are never reset for x2/x3, where every
            // other limit here -- the reconstruction, the pressure floor, the
            // Riemann solver -- is is-1,ie+1.  Measured to change no answer today:
            // the only index it adds is the ghost i = is-1, and w0facewb is zero
            // there because problem generators fill it over active cells only.  It
            // is fixed because the mismatch is a trap, not because it is a live bug:
            // a shifted range also remaps par_for_inner's lane->index assignment
            // relative to the reconstruction above, which would make this
            // read-modify-write touch a scratch slot another lane wrote with no
            // team_barrier in between.  Neither symptom reproduced (5 replicates x
            // 600 cycles, bitwise identical with and without this change).
            AddWbPrimFaceX2(member,m,k,j,is-1,ie+1,w0facewb_x2f,wl_jp1,wr);
          }

          // CUBED SPHERE -- see the note in the x1 sweep. is-1,ie+1 (not il,iu) is the
          // range this sweep actually reconstructs and solves over.
          if (use_cubed_sphere) {
            GnomonicEquianglePrimFaceX2(gtrig_xi,member,m,k,j,is-1,ie+1,
                                                          wl_jp1,wr);
            GnomonicEquiangleFaceBX2(gtrig_xi,member,m,k,j,is-1,ie+1,
                                                       bl_jp1,br);
          }

        member.team_barrier();

        // condwall-1028: a REFLECTING theta wall on spherical polar.  The position-aware
        // x2 reconstruction is not mirror-symmetric about the wall (as the x1 one about a
        // radial wall, see the x1 sweep), so the ghost-side state is not the mirror of
        // the interior one and the wall passed mass (6e-6 of M in t = 1, cw_test sp).
        // Make it the exact mirror (v_theta flipped, B_t as the reflect ghosts: even).
        if (use_spherical_polar &&
            ((j == js && mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::reflect)
             || (j == je+1 &&
                 mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::reflect))) {
          const bool in2 = (j == js);
          par_for_inner(member, is-1, ie+1, [&](const int i) {
            for (int n=0; n<nvars; ++n) {
              if (in2) {wl(n,i) = wr(n,i);} else {wr(n,i) = wl(n,i);}
            }
            if (in2) {wl(IVY,i) = -wr(IVY,i);} else {wr(IVY,i) = -wl(IVY,i);}
            for (int n=0; n<3; ++n) {   // bcc in x1,x2,x3 order (slot 1 unused: by)
              if (in2) {bl(n,i) = br(n,i);} else {br(n,i) = bl(n,i);}
            }
            for (int n=0; n<nder; ++n) {
              if (in2) {dl(n,i) = dr(n,i);} else {dr(n,i) = dl(n,i);}
            }
          });
          member.team_barrier();
        }

        // compute fluxes over [js,je+1].  MHD RS also computes electric fields, where
        // (IBY) component of flx = E_{x} = -(v x B)_{x} = -(v2*b3 - v3*b2)
        // (IBZ) component of flx = E_{z} = -(v x B)_{z} =  (v2*b1 - v1*b2)
        if (j>jl) {
          // NOTE(@pdmullen): Capture variables prior to if constexpr.
          // dhjcuda2-1113: const references, see the x1 sweep.
          const auto &eos = eos_;
          const auto &indcs = indcs_;
          const auto &size = size_;
          const auto &coord = coord_;
          const auto &by = by_;
          const auto &flx2 = flx2_;
          const auto &e12 = e12_;
          const auto &e32 = e32_;
          if (cs_shr_ && k >= ks && k <= ke &&
              ((j == js && size.d_view(m).x2min <= -1.0 + 1.0e-12) ||
               (j == je+1 && size.d_view(m).x2max >= 1.0 - 1.0e-12))) {
            member.team_barrier();
            const bool lo = (j == js);
            const Real sf = gtrig_xi.sn(m,k,j), cf = gtrig_xi.cs(m,k,j);
            const auto &sfs = sfs2_;
            cs_seam::CsSeamGhostPrim(member, sfs, m, k, j, is-1, ie+1, lo, sf, cf, nvars,
                                     nder, IVY, IVZ, wl, wr, dl, dr);
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              // x2 frame {n_xi, e_eta, rhat}: slot 0 = B.r, slot 2 = B.e_eta
              if (lo) {
                bl(0,i) = sfs(m,sfob_,k,j,i); bl(2,i) = sfs(m,sfob_+2,k,j,i);
              } else {
                br(0,i) = sfs(m,sfob_,k,j,i); br(2,i) = sfs(m,sfob_+2,k,j,i);
              }
            });
            member.team_barrier();
          }
          const bool do_pole = polar_hlle_x2_ &&
            ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar &&
              (j == js || j == js+1)) ||
             (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar &&
              (j == je+1 || j == je)));
          if constexpr (rsolver_method_ == MHD_RSolver::advect) {
            Advect(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,by,flx2,e12,e32);
          } else if constexpr (rsolver_method_ == MHD_RSolver::llf) {
            LLF(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlle) {
            HLLE(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlld) {
            if (do_pole) {
              HLLE(member,eos,indcs,size,coord,
                          m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
            } else {
            HLLD(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
            }
          } else if constexpr (rsolver_method_ == MHD_RSolver::lhlld) {
            if (do_pole) {
              HLLE(member,eos,indcs,size,coord,
                          m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
            } else {
              LHLLD(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,by,flx2,e12,e32,
                    w0_);
            }
          } else if constexpr (rsolver_method_ == MHD_RSolver::llf_sr) {
            LLF_SR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,by,flx2,e12,e32);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_sr) {
            HLLE_SR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,by,flx2,e12,e32);
          } else if constexpr (rsolver_method_ == MHD_RSolver::llf_gr) {
            LLF_GR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,by,flx2,e12,e32);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_gr) {
            HLLE_GR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,by,flx2,e12,e32);
          }

          if (!use_cubed_sphere && cs_lb_beta > 0.0) {   // share-1104 G8 (see x1)
            member.team_barrier();
            CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVY,cs_lb_beta,
                              wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
            member.team_barrier();
          }
          if (boris_on_ && !use_cubed_sphere) {
            member.team_barrier();
            BorisHLLEFace(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,
                          by,flx2,e12,e32);
          }
          if (use_cubed_sphere) {
            if (cs_lb_beta > 0.0) {
              CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVY,cs_lb_beta,
                           wl,wr,bl,br,dl,dr,by,flx2,e12,e32);
              member.team_barrier();
            }
            if (boris_on_) {
              member.team_barrier();
              BorisHLLEFace(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVY,wl,wr,bl,br,dl,dr,
                            by,flx2,e12,e32);
              member.team_barrier();
            }
            GnomonicEquiangleFluxX2(gtrig_xi,member,m,k,j,is-1,ie+1,flx2);
          }
          member.team_barrier();
          if (uct_on_) {   // UCT corner EMF inputs (x2 faces)
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              UctFaceStore(eos, wl, wr, bl, br, dl, dr, IVY, i, by(m,k,j,i),
                           uct_hll_ || do_pole,
                           uct_eta_, uct_lb_, uct2_, m, k, j, w0_);
            });
            member.team_barrier();
          }
        }

        // calculate fluxes of scalars (if any)
        if (nvars > nmhd_) {
          for (int n=nmhd_; n<nvars; ++n) {
            par_for_inner(member, is, ie, [&](const int i) {
              if (flx2_(m,IDN,k,j,i) >= 0.0) {
                flx2_(m,n,k,j,i) = flx2_(m,IDN,k,j,i)*wl(n,i);
              } else {
                flx2_(m,n,k,j,i) = flx2_(m,IDN,k,j,i)*wr(n,i);
              }
            });
          }
        }
      } // end of loop over j
    });
  }

  //--------------------------------------------------------------------------------------
  // k-direction. Note order of k,j loops switched

  if (pmy_pack->pmesh->three_d) {
    scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
                ScrArray2D<Real>::shmem_size(3, ncells1) +
                ScrArray2D<Real>::shmem_size(nder, ncells1)) * 3;
    auto &flx3_ = uflx.x3f;
    auto &bz_ = b0.x3f;
    auto &e23_ = e2x3;
    auto &e13_ = e1x3;
      
    auto &x3v_ = pmy_pack->pcoord->x3v;
    auto &x3f_ = pmy_pack->pcoord->xx3f;
    auto &x2v3_ = pmy_pack->pcoord->x2v;
    auto &x2f3_ = pmy_pack->pcoord->xx2f;

    // set the loop limits
    kl = ks-1, ku = ke+1;
    if (use_fofc) { kl = ks-2, ku = ke+2; }

    scr_level = TeamScratchLevel(scr_size, 0);
    par_for_outer("mhd_flux3",DevExeSpace(), scr_size, scr_level, 0, nmb1, js-1, je+1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j) {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr4(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr5(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> scr6(member.team_scratch(scr_level), 3, ncells1);
      ScrArray2D<Real> dscr1(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr2(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr3(member.team_scratch(scr_level), nder, ncells1);

      for (int k=kl; k<=ku; ++k) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_kp1 = scr2;
        auto wr     = scr3;
        auto bl     = scr4;
        auto bl_kp1 = scr5;
        auto br     = scr6;
        auto dl     = dscr1;
        auto dl_kp1 = dscr2;
        auto dr     = dscr3;
        if ((k%2) == 0) {
          wl     = scr2;
          wl_kp1 = scr1;
          bl     = scr5;
          bl_kp1 = scr4;
          dl     = dscr2;
          dl_kp1 = dscr1;
        }
          
        if (use_spherical_polar) {
          GridPiecewiseLinearX3(member, m, k, j, is-1, ie+1, w0_, x3v_, x3f_, wl_kp1, wr);
          GridPiecewiseLinearX3(member, m, k, j, is-1, ie+1, b0_, x3v_, x3f_, bl_kp1, br);
          if (wb_all_) {
            hydro::WbAllBlendGrid(member, use_wb_rho_, 3, false, m, k, j, is-1, ie+1, w0_,
                                  x3v_, x3f_, wbq3_, wbw_, wl_kp1, wr);
          }
        } else {
        // Reconstruct qR[k] and qL[k+1], for both W and Bcc
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX3(member, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            DonorCellX3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
            break;
          case ReconstructionMethod::plm:
            if (eos_.plm_c1) {
              PiecewiseLinearC1X3(member, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            } else {
              PiecewiseLinearX3(member, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            }
            if (eos_.plm_c1) {
              PiecewiseLinearC1X3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
            } else {
              PiecewiseLinearX3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
            }
            break;
          case ReconstructionMethod::ppm4:
          case ReconstructionMethod::ppmx:
            PiecewiseParabolicX3(member,eos_,extrema,true, m,k,j,is-1,ie+1,w0_,wl_kp1,wr);
            PiecewiseParabolicX3(member,eos_,extrema,false,m,k,j,is-1,ie+1,b0_,bl_kp1,br);
            break;
          case ReconstructionMethod::wenoz:
            WENOZX3(member, eos_, true,  m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            WENOZX3(member, eos_, false, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
            break;
          default:
            break;
        }
        if (wb_all_) {
          hydro::WbAllBlendUni(member, use_wb_rho_, 3, false, m, k, j, is-1, ie+1, w0_,
                               wbq3_, wbw_, wl_kp1, wr);
        }
        }

        // Reconstruct the derived thermodynamic variables (general EOS only); range
        // matches the Riemann solver below, as in the j-direction sweep above.
        if (nder > 0) {
          if (use_spherical_polar) {
            GridPiecewiseLinearX3(member, m, k, j, is-1, ie+1, wder_, x3v_, x3f_,
                                  dl_kp1, dr);
          } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX3(member, m, k, j, is-1, ie+1,
                                     pwb_, pfacewb_x3f,
                                    wder_, dl_kp1, dr);
      } else {
          switch (recon_method_) {
            case ReconstructionMethod::dc:
              DonorCellX3(member, m, k, j, is-1, ie+1, wder_, dl_kp1, dr);
              break;
            case ReconstructionMethod::plm:
              if (eos_.plm_c1) {
                PiecewiseLinearC1X3(member, m, k, j, is-1, ie+1, wder_, dl_kp1, dr);
              } else {
                PiecewiseLinearX3(member, m, k, j, is-1, ie+1, wder_, dl_kp1, dr);
              }
              break;
            case ReconstructionMethod::ppm4:
            case ReconstructionMethod::ppmx:
              PiecewiseParabolicX3(member,eos_,extrema,false,m,k,j,is-1,ie+1,
                                   wder_,dl_kp1,dr);
              break;
            case ReconstructionMethod::wenoz:
              WENOZX3(member, eos_, false, m, k, j, is-1, ie+1, wder_, dl_kp1, dr);
              break;
            default:
              break;
          }
          }
          if (wb_all_) {
            if (use_spherical_polar) {
              hydro::WbAllBlendGrid(member, false, 3, true, m, k, j, is-1, ie+1, wder_,
                                    x3v_, x3f_, wbq3_, wbw_, dl_kp1, dr);
            } else {
              hydro::WbAllBlendUni(member, false, 3, true, m, k, j, is-1, ie+1, wder_,
                                   wbq3_, wbw_, dl_kp1, dr);
            }
          }
          par_for_inner(member, is-1, ie+1, [&](const int i) {
            dl_kp1(IDPR,i) = fmax(dl_kp1(IDPR,i), eos_.pfloor);
            dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
          });
        }

          if (use_wellbalance_static_reconst_perturb_) {
            // is-1,ie+1, not il,iu -- see the note on AddWbPrimFaceX2 above.
            AddWbPrimFaceX3(member,m,k,j,is-1,ie+1,w0facewb_x3f,wl_kp1,wr);
          }

          // CUBED SPHERE -- see the note in the x1 sweep.
          if (use_cubed_sphere) {
            // the field needs no rotation in this sweep either.
            GnomonicEquianglePrimFaceX3(gtrig_eta,member,m,k,j,is-1,ie+1,
                                                          wl_kp1,wr);
          }

        // SPHERICAL POLAR: the states were reconstructed along phi at the cell's volume
        // centroid theta = x2v, but an x3 face (area element r dr dtheta) is centred at
        // the theta MIDPOINT of the cell.  The two differ by O(dtheta^2 cot theta),
        // negligible in the interior, but in the polar row the midpoint is dtheta/2
        // against a centroid of 2 dtheta/3 and the phi-face fluxes enter the balance
        // with weight ~2/(r dtheta dphi): the r-components of v and B (~ sin theta) were
        // 4/3 too large -- a resolution-independent spurious radial force of ~0.5
        // (B^2/2)/r (sp_test iprob=11) -- and any scalar with a gradient across the pole
        // (p = p0 + g x, say) was read dtheta/6 off, an O(1) error in the polar cell's
        // phi-force.  Shift EVERY reconstructed variable to the midpoint with the
        // centred theta-derivative of the cell values; the polar ghosts hold the analytic
        // continuation through the axis (the exchange flips the tangential components),
        // so the centred difference is valid in the polar row too.  For a vector this
        // shift includes the basis rotation rhat' = rhat cos d + thhat sin d.
        // <mesh>/polar_x3_shift selects the form: 1 = "rotate", the exact basis rotation
        // of the (r,theta) pairs of v and B alone (a8cfb83e); 2 = "shift", the unlimited
        // centred-difference shift of every variable above (863e8337) -- BISECTED to the
        // polar-row radial-field blow-up of the sp MHD production run (2026-09-10);
        // 3 = "scalar", rotate the vectors and shift only the remaining variables;
        // 0 = "off".
        if (use_spherical_polar) {
          if (px3_ == 1 || px3_ == 3) {
            const Real dang = 0.5*(x2f3_(m,j) + x2f3_(m,j+1)) - x2v3_(m,j);
            const Real cd = cos(dang), sd0 = sin(dang);
            // A REFLECTING radial wall: the ghost row holds the MIRROR image (v_r, B_r
            // flipped) of the interior row, and its phi-face EMFs enter the wall-edge EMF
            // (CornerE), which vanishes only if they are minus the interior ones.
            // Rotating the mirrored (r,theta) pair by +dang is not the mirror of the
            // rotated interior pair; rotating it by -dang is.  With +dang the wall-edge
            // E_theta was off by O(sin(dang) v_theta B_phi): B.rhat grew at the closed
            // wall wherever the flow is not axisymmetric, and the Poynting flux through
            // it drained energy (equator MHD blast, -5e-8 at nx2 = 32; now round-off).
            const bool gin_ = mb_bcs.d_view(m,BoundaryFace::inner_x1) ==
                              BoundaryFlag::reflect;
            const bool gout_ = mb_bcs.d_view(m,BoundaryFace::outer_x1) ==
                               BoundaryFlag::reflect;
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              const Real sd = ((gin_ && i < is) || (gout_ && i > ie)) ? -sd0 : sd0;
              Real a = wl_kp1(IVX,i), b = wl_kp1(IVY,i);
              wl_kp1(IVX,i) = a*cd + b*sd;  wl_kp1(IVY,i) = -a*sd + b*cd;
              a = wr(IVX,i); b = wr(IVY,i);
              wr(IVX,i) = a*cd + b*sd;  wr(IVY,i) = -a*sd + b*cd;
              a = bl_kp1(0,i); b = bl_kp1(1,i);
              bl_kp1(0,i) = a*cd + b*sd;  bl_kp1(1,i) = -a*sd + b*cd;
              a = br(0,i); b = br(1,i);
              br(0,i) = a*cd + b*sd;  br(1,i) = -a*sd + b*cd;
            });
          }
          if (px3_ == 2 || px3_ == 3) {
            const Real fac = (0.5*(x2f3_(m,j) + x2f3_(m,j+1)) - x2v3_(m,j))
                             /(x2v3_(m,j+1) - x2v3_(m,j-1));
            for (int n=0; n<nvars; ++n) {
              if (px3_ == 3 && (n == IVX || n == IVY)) continue;
              par_for_inner(member, is-1, ie+1, [&](const int i) {
                const Real sh = fac*(w0_(m,n,k,j+1,i) - w0_(m,n,k,j-1,i));
                wl_kp1(n,i) += sh;
                wr(n,i) += sh;
              });
            }
            for (int n=0; n<3; ++n) {
              if (px3_ == 3 && n != 2) continue;
              par_for_inner(member, is-1, ie+1, [&](const int i) {
                const Real sh = fac*(b0_(m,n,k,j+1,i) - b0_(m,n,k,j-1,i));
                bl_kp1(n,i) += sh;
                br(n,i) += sh;
              });
            }
          }
        }

        member.team_barrier();

        // compute fluxes over [ks,ke+1].  MHD RS also computes electric fields, where
        // (IBY) component of flx = E_{y} = -(v x B)_{y} = -(v3*b1 - v1*b3)
        // (IBZ) component of flx = E_{x} = -(v x B)_{x} =  (v3*b2 - v2*b3)
        if (k>kl) {
          // NOTE(@pdmullen): Capture variables prior to if constexpr.
          // dhjcuda2-1113: const references, see the x1 sweep.
          const auto &eos = eos_;
          const auto &indcs = indcs_;
          const auto &size = size_;
          const auto &coord = coord_;
          const auto &bz = bz_;
          const auto &flx3 = flx3_;
          const auto &e23 = e23_;
          const auto &e13 = e13_;
          if (cs_shr_ && j >= js && j <= je &&
              ((k == ks && size.d_view(m).x3min <= -1.0 + 1.0e-12) ||
               (k == ke+1 && size.d_view(m).x3max >= 1.0 - 1.0e-12))) {
            member.team_barrier();
            const bool lo = (k == ks);
            const Real sf = gtrig_eta.sn(m,k,j), cf = gtrig_eta.cs(m,k,j);
            const auto &sfs = sfs3_;
            cs_seam::CsSeamGhostPrim(member, sfs, m, k, j, is-1, ie+1, lo, sf, cf, nvars,
                                     nder, IVZ, IVY, wl, wr, dl, dr);
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              // x3 frame {n_eta, rhat, e_xi}: slot 0 = B.r, slot 1 = B.e_xi
              if (lo) {
                bl(0,i) = sfs(m,sfob_,k,j,i); bl(1,i) = sfs(m,sfob_+1,k,j,i);
              } else {
                br(0,i) = sfs(m,sfob_,k,j,i); br(1,i) = sfs(m,sfob_+1,k,j,i);
              }
            });
            member.team_barrier();
          }
          const bool do_pole = polar_hlle_x3_ &&
            ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar &&
              j == js) ||
             (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je));
          if constexpr (rsolver_method_ == MHD_RSolver::advect) {
            Advect(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,bz,flx3,e23,e13);
          } else if constexpr (rsolver_method_ == MHD_RSolver::llf) {
            LLF(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlle) {
            HLLE(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlld) {
            if (do_pole) {
              HLLE(member,eos,indcs,size,coord,
                        m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
            } else {
            HLLD(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
            }
          } else if constexpr (rsolver_method_ == MHD_RSolver::lhlld) {
            if (do_pole) {
              HLLE(member,eos,indcs,size,coord,
                        m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
            } else {
              LHLLD(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,bz,flx3,e23,e13,
                    w0_);
            }
          } else if constexpr (rsolver_method_ == MHD_RSolver::llf_sr) {
            LLF_SR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,bz,flx3,e23,e13);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_sr) {
            HLLE_SR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,bz,flx3,e23,e13);
          } else if constexpr (rsolver_method_ == MHD_RSolver::llf_gr) {
            LLF_GR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,bz,flx3,e23,e13);
          } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_gr) {
            HLLE_GR(member,eos,indcs,size,coord,
                    m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,bz,flx3,e23,e13);
          }

          if (!use_cubed_sphere && cs_lb_beta > 0.0) {   // share-1104 G8 (see x1)
            member.team_barrier();
            CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVZ,cs_lb_beta,
                              wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
            member.team_barrier();
          }
          if (boris_on_ && !use_cubed_sphere) {
            member.team_barrier();
            BorisHLLEFace(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,
                          bz,flx3,e23,e13);
          }
          if (use_cubed_sphere) {
            if (cs_lb_beta > 0.0) {
              CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVZ,cs_lb_beta,
                           wl,wr,bl,br,dl,dr,bz,flx3,e23,e13);
              member.team_barrier();
            }
            if (boris_on_) {
              member.team_barrier();
              BorisHLLEFace(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVZ,wl,wr,bl,br,dl,dr,
                            bz,flx3,e23,e13);
              member.team_barrier();
            }
            GnomonicEquiangleFluxX3(gtrig_eta,member,m,k,j,is-1,ie+1,flx3);
          }
          member.team_barrier();
          if (uct_on_) {   // UCT corner EMF inputs (x3 faces)
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              UctFaceStore(eos, wl, wr, bl, br, dl, dr, IVZ, i, bz(m,k,j,i),
                           uct_hll_ || do_pole,
                           uct_eta_, uct_lb_, uct3_, m, k, j, w0_);
            });
            member.team_barrier();
          }
        }

        // calculate fluxes of scalars (if any)
        if (nvars > nmhd_) {
          for (int n=nmhd_; n<nvars; ++n) {
            par_for_inner(member, is, ie, [&](const int i) {
              if (flx3_(m,IDN,k,j,i) >= 0.0) {
                flx3_(m,n,k,j,i) = flx3_(m,IDN,k,j,i)*wl(n,i);
              } else {
                flx3_(m,n,k,j,i) = flx3_(m,IDN,k,j,i)*wr(n,i);
              }
            });
          }
        }
      } // end loop over k
    });
  }

  // whole-mesh Cartesian momentum on spherical polar: fourth-order face averages of the
  // hydrodynamic fluxes (see Coordinates::SphericalPolarFaceAverageFluxes)
  if (pmy_pack->pcoord->sp_cart_all_momentum) {
    pmy_pack->pcoord->SphericalPolarFaceAverageFluxes(uflx, IEN+1);
  }

  return;
}

// function definitions for each template parameter
template void MHD::CalculateFluxes<MHD_RSolver::advect>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::llf>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlle>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlld>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::lhlld>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::llf_sr>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlle_sr>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::llf_gr>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlle_gr>(Driver *pdriver, int stage);

} // namespace mhd
