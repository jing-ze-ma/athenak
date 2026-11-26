//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_closure.cpp
//! \brief Closure of the clean implicit radiation solve: the formal solution on the
//! predicted end-of-step state gives each cell's Eddington tensor D, half-range factors
//! h+- and non-locality C1.
//!
//!  called from: rad_vet_tasks.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_closure_faces.cpp, rad_vet_closure_rays.cpp,
//!    rad_vet_closure_ring.cpp, rad_vet_init.cpp, rays/sc_build.cpp,
//!    rays/sc_cart_sweep.cpp, rays/sc_directions.cpp, rays/sc_source.cpp,
//!    rays/sc_sph_walls.cpp

#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "rad_vet/rad_vet_momentum.hpp"
#include "rad_vet/rad_vet_moments.hpp"
#include "utils/imex_sdirk2.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::BuildClosure
//! \brief builds the closure for t^{n+1}, once per step and dt (all stage solves end
//! there); closure_form = eddington: isotropic

void RadiationVet::BuildClosure() {
  if (cln_sphere) {
    if (!cln_havec) ClosureSphere();
    cln_havec = true;
    return;
  }
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  if (cln_havec && time == cln_tc && dt == cln_dtb) return;

  if (cln_edd) {
    ClosureEddington();
  } else {
    PredictState(true);
    TraceRays();
    if (cln_opsaved) Kokkos::deep_copy(opac, cln_opb);
    ClosureFill();
  }
  GhostFill(cln_s, 0, cln_s.extent_int(1));

  cln_tc = time;
  cln_dtb = dt;
  cln_havec = true;
  cln_cip = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureEddington
//! \brief closure_form = eddington: D = I/3, h+- = +-1/4, k+ = I/6, C1 = 0

void RadiationVet::ClosureEddington() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto cs = cln_s;
  SolveFor("cln_edd", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    for (int a = 0; a < 3; ++a) {
      cs(m,CLN_D11+a,k,j,i) = 1.0/3.0;
      cs(m,CLN_D11+3+a,k,j,i) = 0.0;
      cs(m,CLN_HP+a,k,j,i) = 0.25;
      cs(m,CLN_HM+a,k,j,i) = -0.25;
    }
    for (int a = 0; a < 9 && CLN_KP + a < cs.extent_int(1); ++a) {
      cs(m,CLN_KP+a,k,j,i) = (a%4 == 0) ? 1.0/6.0 : 0.0;
    }
    cs(m,CLN_C1,k,j,i) = 0.0;
    cs(m,CLN_J,k,j,i) = 0.0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureIterate
//! \brief closure_form = iterate: closure from the formal solution of the current
//! iterate's E, T and F

void RadiationVet::ClosureIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto u0_ = u0;
  auto w = cln_w;
  par_for("cln_clositer", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,VET_IW_EN,k,j,i) = fmax(w(m,CLW_EK,k,j,i), 0.0);
    iw_(m,VET_IW_TP,k,j,i) = w(m,CLW_T,k,j,i);
    iw_(m,VET_IW_F1,k,j,i) = u0_(m,VET_F1,k,j,i);
    iw_(m,VET_IW_F2,k,j,i) = u0_(m,VET_F1+1,k,j,i);
    iw_(m,VET_IW_F3,k,j,i) = u0_(m,VET_F1+2,k,j,i);
  });
  TraceRays();
  ClosureFill();
  GhostFill(cln_s, 0, cln_s.extent_int(1));
  cln_havec = false;
  cln_cip = true;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RayMoments
//! \brief ray moments J, K, H+-, <|I - S|> and face fluxes (CLN_*, not normalised) into
//! dst: of a new formal solution of the current state, or (solve = false) of the closure
//! in place

void RadiationVet::RayMoments(DvceArray5D<Real> dst, const bool solve) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  if (solve) {
    PredictState(false);
    TraceRays();
    if (cln_opsaved) Kokkos::deep_copy(opac, cln_opb);
  }

  auto cs = cln_s;
  SolveFor("cln_raymom", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jr = cs(m,CLN_J,k,j,i);
    for (int n = 0; n < cs.extent_int(1); ++n) {
      dst(m,n,k-ks,j-js,i-is) = ((CellMoment(n) && n != CLN_J) ? jr : 1.0)*
                                cs(m,n,k,j,i);
    }
  });
  cln_havec = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureAfterSolve
//! \brief closure_form = extrap(2): store the ray moments of a converged solution: the
//! step's last solve at t^{n+1}, extrap2 also stage 1 at its time.  extrap1: the stage-1
//! closure's moments at t^n + (a21 + g) dt; only backward Euler solves here

void RadiationVet::ClosureAfterSolve() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  if (cln_extrap == 3 && t2_solve == VET_T2S_STAGE2) return;
  if (cln_extrap == 3 && t2_solve == VET_T2S_STAGE1) {
    RayMoments(cln_mc, false);
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    cln_tm = time + ((t2_tab == 0) ? dt : (tb.a21 + tb.g)*dt);
    cln_have_c = true;
    cln_cand_be = false;
    RingJump();
    return;
  }
  if (t2_solve == VET_T2S_STAGE1) {
    if (cln_extrap != 2) return;
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    RayMoments(cln_m1);
    cln_t1 = time + (tb.a21 + tb.g)*dt;
    cln_have_1 = true;
    return;
  }
  RayMoments(cln_mc);
  cln_tm = time + dt;
  cln_have_c = true;
  cln_cand_be = (t2_solve != VET_T2S_STAGE2);
  RingJump();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingSolveSet
//! \brief ring = jump: the face half ranges use the swept set of the step's
//! rotation block

void RadiationVet::RingSolveSet() {
  if (cln_cssc) return;   // (cubed-sphere sc rays: RingUnion's two sets)
  const int64_t bs = RotBlock(pmy_pack->pmesh->ncycle);
  const int ns = static_cast<int>(cln_fvdm.size());
  for (int q = 0; q < ns && ns > 1; ++q) {
    if (cln_fvblk[q] == bs) Kokkos::deep_copy(cln_fvd, cln_fvdm[q]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RingJump
//! \brief ring = jump: at a rotation-block end the moments jump to the next
//! block's; the next step does not extrapolate through it

void RadiationVet::RingJump() {
  if (!cln_rjump) return;
  Kokkos::deep_copy(cln_mc, cln_rnew);
  cln_cand_be = true;
  cln_rjump = false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureExtrap
//! \brief closure_form = extrap(2): the closure at the stage time from ray moments of
//! converged solutions, extrapolated linearly in time (exact for S and inflow linear in
//! t).  The line runs through t^{n-1} and t^n, or (stage 2, extrap2) t^n and stage 1.
//! Cells whose extrapolate is not realisable keep the latest moments (counted); first
//! step, restart and backward Euler use t^n alone (or build it).  extrap1 with stage 1
//! before t^{n+1} (sdirk2, trbdf2): stage 1 uses the last step's moments, stage 2 the
//! line through the last and this step's stage-1 sets, weighted 1/(1 + z), z = c dt/dx
//! (no over-relaxation where the field is quasi-static)

void RadiationVet::ClosureExtrap() {
  const Real time = pmy_pack->pmesh->time;
  const Real dt = pmy_pack->pmesh->dt;
  cln_nguard = 0;
  const bool sx = ExtrapSx();   // (stage-1 moments before t^n+1)
  // (stage 2 of sx: the line through two stage-1 sets, the last not from backward
  // Euler, a restart or a ring jump)
  const bool usec = sx && t2_solve == VET_T2S_STAGE2 && cln_have_c && !cln_cand_be &&
                    cln_have_n && cln_tn < time && cln_tm > cln_tn;
  if (cln_extrap == 3 && t2_solve == VET_T2S_STAGE2 && !usec) return;

  // a new step: the end of the last step (sx: its stage-1 time) becomes t^n
  if (cln_have_c && (cln_tm == time || (sx && cln_tm < time))) {
    if (!sx) std::swap(cln_mp, cln_mn);   // (sx: two sets, the line never reads t^{n-1})
    std::swap(cln_mn, cln_mc);
    cln_tp = cln_tn;
    cln_tn = cln_tm;
    cln_have_p = cln_have_n && !cln_cand_be;
    cln_have_n = true;
    cln_have_c = false;
    cln_have_1 = false;
  }
  if (!cln_have_n || !(cln_tn == time || (sx && cln_tn < time))) {
    BuildClosure();
    return;
  }

  // the line M(t) = Mb + (t - tb)/(tb - ta) (Mb - Ma) at the stage time ts
  const bool stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  Real ts = time + dt;
  if (t2_solve == VET_T2S_STAGE1) {
    const imex_sdirk2::Tableau tb = imex_sdirk2::Make(t2_tab);
    ts = time + (tb.a21 + tb.g)*dt;
  }
  const bool use1 = (cln_extrap == 2) && (t2_solve == VET_T2S_STAGE2) && cln_have_1 &&
                    (cln_t1 > cln_tn);
  const bool line = stage && (use1 || usec || (cln_have_p && !sx));
  auto ma = (use1 || usec || sx) ? cln_mn : cln_mp;
  auto mb = use1 ? cln_m1 : (usec ? cln_mc : cln_mn);
  const Real tb = use1 ? cln_t1 : (usec ? cln_tm : cln_tn);
  const Real ta = (use1 || usec) ? cln_tn : cln_tp;
  const Real wx = line ? (ts - tb)/(tb - ta) : 0.0;
  auto mn = usec ? cln_mc : cln_mn;   // (a guarded cell: the latest moments)

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ni = ie - is + 1, nj = je - js + 1, nk = ke - ks + 1;
  auto cs = cln_s;
  const int ns = cs.extent_int(1);
  auto valid = cln_v;
  SolveArgs x = MakeArgs();
  int nguard = 0;
  Kokkos::parallel_reduce("cln_extrap",
                          Kokkos::RangePolicy<DevExeSpace>(0, (nmb1 + 1)*nk*nj*ni),
  KOKKOS_LAMBDA(const int t, int &ng) {
    const int i = is + t%ni;
    const int j = js + (t/ni)%nj;
    const int k = ks + (t/(ni*nj))%nk;
    const int m = t/(ni*nj*nk);
    // (stage 2 of sx: the line times 1/(1 + z), z = c dt/dx of the thinnest direction)
    Real z = 0.0;
    for (int d = 0; d < x.nf && usec; ++d) {
      z = fmax(z, x.c*dt*DivF(x, d, m, k, j, i, 0.0, 1.0));
    }
    const Real wc = wx/(1.0 + z);
    Real q[CLN_N];
    const int kr = k - ks, jr = j - js, ir = i - is;
    for (int n = 0; n < CLN_N; ++n) {
      q[n] = (n >= ns) ? 0.0 :
             mb(m,n,kr,jr,ir) + wc*(mb(m,n,kr,jr,ir) - ma(m,n,kr,jr,ir));
    }
    // the extrapolated moments must be those of a non-negative intensity
    const bool ok = line && ma(m,CLN_J,kr,jr,ir) > 0.0 && q[CLN_J] > 0.0 &&
                    Realizable(x, q, m, k, j, i);
    if (!ok) {
      if (line && mn(m,CLN_J,kr,jr,ir) > 0.0) ++ng;
      // (cubed sphere: face sums, held on both sides, stay extrapolated alike)
      for (int n = 0; n < ns; ++n) {
        if (!x.cub || CellMoment(n)) q[n] = mn(m,n,kr,jr,ir);
      }
    }
    StoreMoments(cs, valid, q, m, k, j, i);
  }, Kokkos::Sum<int>(nguard));
  cln_nguard = static_cast<int>(GlobalReduce(nguard, false));
  cln_nguard_sum += cln_nguard;
  if (line) cln_next_sum += 1.0;

  ClosureFill();
  GhostFill(cln_s, 0, cln_s.extent_int(1));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::PredictState
//! \brief predict the t^{n+1} state for the formal solution (iw: E, T, F; opac if
//! T-dependent), O(dt^2), from the hydro's explicit stage and K1 = f(U^n) (FSAL):
//!  1. momentum: hydro + dt K1_M (cancels the explicit gravity kick where supported);
//!  2. transport: dt (K1_E + K1_EN), losses beyond the content at constant relative rate;
//!  3. exchange: implicit local relaxation of (eint, E) over dt.
//! No valid K1 (first step, backward Euler): only 3.  No feedback: gas fixed.
//! use_k1 = false: the current state itself

void RadiationVet::PredictState(const bool use_k1) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int n1 = indcs.nx1 + 2*indcs.ng;
  int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  Real dt = pmy_pack->pmesh->dt;

  FluidRef fl = FluidRef::Get(pmy_pack);
  const bool fluid = fl.on;
  const bool exch = use_k1;
  const bool rates = exch && t2_ok && (t2k1.extent_int(0) > nmb1);

  // opacities at the predicted state; the solve itself uses those of the stage
  cln_opsaved = fluid && !opac_zero && (opacity_type != VET_OPAC_CONST);
  if (cln_opsaved) {
    if (cln_opb.extent_int(1) != opac.extent_int(1)) {
      Kokkos::realloc(cln_opb, opac.extent_int(0), opac.extent_int(1),
                      opac.extent_int(2), opac.extent_int(3), opac.extent_int(4));
    }
    Kokkos::deep_copy(cln_opb, opac);
  }
  const bool new_opac = cln_opsaved;

  auto iw_ = iw;
  auto u0_ = u0;
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto k1 = t2k1;
  auto opac_ = opac;
  int otype = opacity_type;
  Real kp = kappa_p, ke_ = kappa_e, kf = kappa_f, ksc = kappa_s;
  Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  VetOpacTab otab_ = otab;
  bool scat = vscat;
  Real arad_ = arad, cl = c_light;
  const bool feedback = cln_feedback;
  SolveArgs x = MakeArgs();

  SolveFor("cln_predict", nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // K1 is stored on active cells only; ghost cells take the explicit state
    const bool active = (k >= ks && k <= ke && j >= js && j <= je && i >= is && i <= ie);
    const bool use_rates = rates && active;
    const bool relax = exch && feedback && active;
    const Real e0 = u0_(m,VET_E,k,j,i);
    Real e = e0, d = 1.0, t = 0.0;
    if (fluid) {
      d = uh(m,IDN,k,j,i);
      Real mom[3] = {uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i)};
      const Real ekin = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
      const Real eint_h = fmax(GasEint(x, uh, ekin, m, k, j, i), 0.0);
      t = eos.Temperature(d, eint_h);
      if (use_rates && !feedback) {
        // gas held fixed: E advances at its own rate, which includes the exchange
        const Real de = dt*k1(m,VET_T2_E,k,j,i);
        e = PredictE(e0, de);
      } else if (relax) {
        // transport over dt (valid K1) less the radiation force's work, then the
        // implicit local exchange of (eint, E)
        if (use_rates) {
          for (int n = 0; n < 3; ++n) mom[n] += dt*k1(m,VET_T2_M1+n,k,j,i);
          const Real ekin_new = 0.5*(SQR(mom[0]) + SQR(mom[1]) + SQR(mom[2]))/d;
          const Real de = dt*(k1(m,VET_T2_E,k,j,i) + k1(m,VET_T2_EN,k,j,i))
                          - (ekin_new - ekin);
          e = PredictE(e0, de);
        }
        // E (1 + h c kE) = e + h c kP a T^4,  eint(T) + E = eint_h + e
        const Real hce = dt*cl*opac_(m,VET_OP_E,k,j,i);
        const Real hcp = dt*cl*opac_(m,VET_OP_P,k,j,i)*arad_;
        t = ExchangeTemperature(eos, d, hcp/(1.0 + hce), eint_h + e*hce/(1.0 + hce),
                                t);
        e = (e + hcp*SQR(SQR(t)))/(1.0 + hce);
      }
    } else if (use_rates) {
      const Real de = dt*k1(m,VET_T2_E,k,j,i);
      e = PredictE(e0, de);
    }

    if (active) {
      iw_(m,VET_IW_EN,k,j,i) = e;
      iw_(m,VET_IW_TP,k,j,i) = t;
      iw_(m,VET_IW_F1,k,j,i) = u0_(m,VET_F1,k,j,i);
      iw_(m,VET_IW_F2,k,j,i) = u0_(m,VET_F1+1,k,j,i);
      iw_(m,VET_IW_F3,k,j,i) = u0_(m,VET_F1+2,k,j,i);
    }
    if (new_opac) {
      Real op, oe, of, os;
      if (otype == VET_OPAC_TABLE) {
        VetTableOpacities(otab_, d, t, op, oe, of, os);
        VetScatterEos(otab_, eos, d, t, of, os);
      } else {
        VetOpacities(otype, d, t, kp, ke_, kf, ksc, rref, tref, aa, bb, op, oe, of, os);
      }
      opac_(m,VET_OP_P,k,j,i) = d*op;
      opac_(m,VET_OP_E,k,j,i) = d*oe;
      opac_(m,VET_OP_T,k,j,i) = d*(of + os);
      if (scat) opac_(m,VET_OP_S,k,j,i) = d*os;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureSphere
//! \brief exact closure of a uniformly bright sphere of radius R = closure_sphere_radius
//! in vacuum.  Outside, the cone mu > mu_c = sqrt(1 - R^2/r^2): D_rr = (1 + mu_c +
//! mu_c^2)/3, D_tt = D_pp = (1 - D_rr)/2, h+_r = (1 + mu_c)/2, h-_r = 0, h+-_t = h+-_p =
//! +-S/(pi (1 - mu_c)), S = int_{mu_c}^1 sqrt(1 - mu^2) dmu, C1 = 1; k+^r = (D_rr, 0, 0),
//! k+^t = (k, D_tt/2, 0), k+^p = (k, 0, D_pp/2), k = (R/r)(1 + mu_c)/(3 pi).  Inside:
//! isotropic.  Cubed sphere: Cartesian about r_hat and the lateral family normals n_d

void RadiationVet::ClosureSphere() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const Real rs = cln_sph_r;
  const bool cub = pmy_pack->pmesh->use_cubed_sphere;
  auto x1v = pmy_pack->pcoord->x1v;
  auto cs = cln_s;
  SolveArgs x = MakeArgs();

  SolveFor("cln_sphere", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    Real drr = 1.0/3.0, hpr = 0.25, hmr = -0.25, ht = 0.25, c1 = 0.0, kr = 0.0,
         krr = 1.0/6.0;
    if (r > rs) {
      const Real q = rs/r;
      const Real mc = std::sqrt(fmax(1.0 - q*q, 0.0));
      drr = (1.0 + mc + mc*mc)/3.0;
      hpr = 0.5*(1.0 + mc);
      hmr = 0.0;
      // S/(pi (1 - mu_c)) with 1 - mu_c = q^2/(1 + mu_c); S ~ 2 q^3/3 for small q
      const Real s = 0.5*(std::asin(q) - mc*q);
      ht = (q > 1.0e-4) ? s*(1.0 + mc)/(M_PI*q*q) : 2.0*q/(3.0*M_PI);
      c1 = 1.0;
      kr = q*(1.0 + mc)/(3.0*M_PI);
      krr = drr;
    }
    cs(m,CLN_D11,k,j,i) = drr;
    cs(m,CLN_D11+1,k,j,i) = 0.5*(1.0 - drr);
    cs(m,CLN_D11+2,k,j,i) = 0.5*(1.0 - drr);
    for (int a = 3; a < 6; ++a) cs(m,CLN_D11+a,k,j,i) = 0.0;
    cs(m,CLN_HP,k,j,i) = hpr;
    cs(m,CLN_HM,k,j,i) = hmr;
    for (int a = 1; a < 3; ++a) {
      cs(m,CLN_HP+a,k,j,i) = ht;
      cs(m,CLN_HM+a,k,j,i) = -ht;
    }
    cs(m,CLN_C1,k,j,i) = c1;
    cs(m,CLN_J,k,j,i) = 1.0;
    const bool kp = cs.extent_int(1) > CLN_KP + 8;   // (k+, not CLN_HO)
    if (kp) {
      for (int a = 0; a < 9; ++a) cs(m,CLN_KP+a,k,j,i) = 0.0;
      cs(m,CLN_KP,k,j,i) = krr;
      cs(m,CLN_KP+3,k,j,i) = kr;
      cs(m,CLN_KP+4,k,j,i) = 0.25*(1.0 - drr);
      cs(m,CLN_KP+6,k,j,i) = kr;
      cs(m,CLN_KP+8,k,j,i) = 0.25*(1.0 - drr);
    }
    if (!cub) return;
    const Real dt = 0.5*(1.0 - drr);
    Real nd[3][3];
    CsFamily(x, m, k, j, i, nd);
    const int ab[6][2] = {{0, 0}, {1, 1}, {2, 2}, {0, 1}, {0, 2}, {1, 2}};
    for (int a = 0; a < 6; ++a) {
      cs(m,CLN_D11+a,k,j,i) = (drr - dt)*nd[0][ab[a][0]]*nd[0][ab[a][1]] +
                              ((a < 3) ? dt : 0.0);
    }
    for (int b = 0; b < 3 && kp; ++b) {
      cs(m,CLN_KP+b,k,j,i) = krr*nd[0][b];
      for (int d = 1; d < 3; ++d) {
        cs(m,CLN_KP+3*d+b,k,j,i) = kr*nd[0][b] + 0.5*dt*nd[d][b];
      }
    }
  });
  GhostFill(cln_s, 0, cln_s.extent_int(1));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::TraceRays
//! \brief formal solution of the state in iw (chi from opac) and its closure: vet_sc,
//! or vet_gd on curved meshes; the state's own source (no dt-dependent relaxation)

void RadiationVet::TraceRays() {
  // (the closure's time; a step of at least one rotation block takes its start's sets)
  const Real dt = pmy_pack->pmesh->dt, tc = pmy_pack->pmesh->time +
                  ((vgd_rot_dt > 0.0 && dt >= vgd_rot_dt) ? 0.0 : dt);
  if (!vgd_on && vgd_slab < 0) ScSlabSelect();
  if (!vgd_on && vgd_slab == 0 && cln_ring) {   // (the rotating vet_sc set: RingMoments)
    RingSets(tc, pmy_pack->pmesh->ncycle);
    RingMoments();
    return;
  }
  if (!vgd_on && vgd_slab == 0) {
    ScCartBuild();
    ClosureFromRays();
    return;
  }
  if (!vray_ready) VetRayInit();
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto iw_ = iw;
  auto opac_ = opac;
  par_for("cln_gd_chi", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    iw_(m,VET_IW_KT,k,j,i) = opac_(m,VET_OP_T,k,j,i);
  });
  if (cln_ring) {   // (the sc rays: each set's sweep in RingMoments)
    ScSource();
    RingSets(tc, pmy_pack->pmesh->ncycle);
  } else {
    ScBuild();
    if (cln_scmem) ScTopCopy();
  }
  if (cln_ring) {
    RingMoments();
  } else {
    ClosureFromSc();
    // (cubed-sphere sc rays: also the face half ranges from cln_fvi, FaceHalf)
    // (Cartesian slab: the face sums are made in the sweep)
    if ((cln_consist && vgd_slab != 1) || cln_cssc) RayFaceSums();
    if (cln_cssc && !cln_rfh) FaceHalf(vgd_dir, 0, true, 0);
    if (cln_cssc) cln_fvi = DvceArray5D<Real>();   // (build scratch)
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ClosureFill
//! \brief cells without ray radiation take the mean closure (D, h+-, k+) of their valid
//! neighbours, one layer per pass (clean_fill_passes passes).  face_thin: the
//! first layer's, next to the rays, weighted by their E, so that a neighbour the rays
//! hardly light does not count as much as a beam; the deeper layers hold no radiation

void RadiationVet::ClosureFill() {
  cln_sout = true;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const int dj = pmy_pack->pmesh->multi_d ? 1 : 0;
  const int dk = pmy_pack->pmesh->three_d ? 1 : 0;
  auto cs = cln_s;
  auto valid = cln_v;
  // (D, h+- and, where held, k+: the slots t < CLN_C1, then CLN_KP + t - CLN_C1)
  const int nfl = CLN_C1 + std::max(std::min(cln_s.extent_int(1) - CLN_KP, 9), 0);
  auto iw_ = iw;
  const bool thf = cln_thinf;
  const Real efl = fmax(e_floor, 1.0e-300);
  // (no cell without rays: the passes would change nothing)
  const int ni = ie - is + 1, nji = (je - js + 1)*ni, nkji = (ke - ks + 1)*nji;
  Real nbad = 0.0;
  Kokkos::parallel_reduce("cln_fill_n", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          (nmb1 + 1)*nkji),
  KOKKOS_LAMBDA(const int idx, Real &a) {
    const int m = idx/nkji, k = (idx - m*nkji)/nji, j = (idx - m*nkji - k*nji)/ni;
    if (!(valid(m,0,ks+k,js+j,is+idx-m*nkji-k*nji-j*ni) > 0.0)) a += 1.0;
  }, Kokkos::Sum<Real>(nbad));
  if (GlobalReduce(nbad, false) == 0.0) return;
  if (thf) {   // the weight, E
    par_for("cln_fill_w", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,2,k,j,i) = fmax(iw_(m,VET_IW_EN,k,j,i), efl);
    });
    GhostFill(cln_v, 2, 1);
  }

  constexpr int npass = 8;   // fill layers
  for (int pass = 0; pass < npass; ++pass) {
    // (only what cln_fill reads: D, h+- (and the held k+ or HO), valid(0), at depth 1;
    // callers exchange all of cln_s afterwards; the last pass also valid(1))
    GhostDefer(cln_s, 0, CLN_C1, 1);
    if (nfl > CLN_C1) GhostDefer(cln_s, CLN_KP, nfl - CLN_C1, 1);
    GhostFill(cln_v, 0, (pass + 1 < npass) ? 1 : 2, true, 0, 1);
    // valid(1) marks the cells filled in this pass; only cells valid before it are read
    par_for("cln_fill", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,1,k,j,i) = valid(m,0,k,j,i);
      if (valid(m,0,k,j,i) > 0.0) return;
      Real sum[CLN_C1 + 9] = {0.0}, nv = 0.0;
      for (int kk = k-dk; kk <= k+dk; ++kk) {
        for (int jj = j-dj; jj <= j+dj; ++jj) {
          for (int ii = i-1; ii <= i+1; ++ii) {
            if (valid(m,0,kk,jj,ii) <= 0.0) continue;
            const Real wt = (thf && pass == 0) ? valid(m,2,kk,jj,ii) : 1.0;
            nv += wt;
            for (int t = 0; t < nfl; ++t) {
              sum[t] += wt*cs(m,(t < CLN_C1) ? t : CLN_KP + t - CLN_C1,kk,jj,ii);
            }
          }
        }
      }
      if (nv == 0.0) return;
      valid(m,1,k,j,i) = 1.0;
      for (int t = 0; t < nfl; ++t) {
        cs(m,(t < CLN_C1) ? t : CLN_KP + t - CLN_C1,k,j,i) = sum[t]/nv;
      }
    });
    par_for("cln_fill_mark", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      valid(m,0,k,j,i) = valid(m,1,k,j,i);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::StageClosure
//! \brief the closure a stage solve starts with

void RadiationVet::StageClosure() {
  // from the gas before the stage increments; closure_form = iterate with no closure
  // in place takes the predicted one
  if (cln_extrap > 0) {
    ClosureExtrap();
  } else if (!cln_closit || !cln_cip) {
    BuildClosure();
  }
  if (cln_ring == 2) RingSolveSet();
}

}  // namespace radvet
