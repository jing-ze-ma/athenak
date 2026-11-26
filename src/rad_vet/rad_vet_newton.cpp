//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_newton.cpp
//! \brief the Newton / Picard passes on the gas coupling: first iterate, gas temperature
//! and velocity updates, opacity of the iterate, face bracket, counts of the solution.
//!
//!  called from: rad_vet_tasks.cpp
//!  calls into: rad_vet_bvals.cpp, rad_vet_flux.cpp, rad_vet_flux_seam.cpp,
//!    rad_vet_init.cpp, rad_vet_solver.cpp

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "eos/eos.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "rad_vet/rad_vet_momentum.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void ReduceMax
//! \brief the maxima of a, b, c over all ranks, one reduction (GlobalReduce)

static void ReduceMax(Real &a, Real &b, Real &c) {
  Real v[3] = {a, b, c};
  MaxOverRanks(v, 3);
  a = v[0];
  b = v[1];
  c = v[2];
}

//----------------------------------------------------------------------------------------
//! \fn void VelocityCell
//! \brief VelocityUpdate's cell (k,j,i): the momentum change (CLW_DM), the advective
//! coefficients of the new velocity if fa (CLW_A1), and the maxima of |dv|, |v| and the
//! force scale; kFs: F known from the F solve (xh.fs, xh.fa)

template <bool kFs>
KOKKOS_INLINE_FUNCTION
static void VelocityCell(const SolveArgs &x, const SolveArgs &xh,
                         const SolveArgs &xa,
                         const DvceArray5D<Real> &uh,
                         const int fref, const DvceArray5D<Real> &fr,
                         const DvceArray4D<Real> &aref, const bool cap,
                         const DvceArray5D<Real> &op, const Real arad_,
                         const bool fa,
                         const int m, const int k, const int j, const int i,
                         Real &dv_max, Real &v_max, Real &fv_max) {
  const Real d = uh(m,IDN,k,j,i);
  Real v2 = 0.0, dmv[3], mg[3], pg[3], rf[3];
  CellMomentumChange<kFs>(xh, CLW_EK, m, k, j, i, dmv, mg);
  FaceFrameGasMomentum(x, uh, m, k, j, i, pg);
  RefForce(x, fref, fr, aref, d, m, k, j, i, rf);
  for (int dir = 0; dir < 3; ++dir) {
    Real dm = 0.0;
    if (dir == 0 || (dir == 1 && x.multi_d) || (dir == 2 && x.three_d)) {
      const Real mag = mg[dir];
      dm = dmv[dir];
      v_max = fmax(v_max, mag/d);   // |v| floor: size of the stage's force terms
      // force_reference = clean: the first coupled solve's force is the reference,
      // withheld from the gas (as wb_arad's)
      if (cap) fr(m,dir,k,j,i) = dm/x.h;
      if (fref != VET_FREF_NONE) dm -= cap ? dm : x.h*rf[dir];
      if (x.vzero) dm = 0.0;
      // force scale relative to E in the stress
      const Real e = fmax(fabs(x.w(m,CLW_EK,k,j,i)), 1.0e-300), hc = x.h*x.c;
      fv_max = fmax(fv_max, mag/d*(2.0 + (e*hc*op(m,VET_OP_E,k,j,i) + fabs(x.w(m,
               CLW_EO,k,j,i)) + hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(x.w(m,CLW_T,k,j,
               i))))/e));
    }
    dv_max = fmax(dv_max, fabs(dm - x.w(m,CLW_DM+dir,k,j,i))/d);
    x.w(m,CLW_DM+dir,k,j,i) = dm;
    v2 += SQR((pg[dir] + dm)/d);
  }
  v_max = fmax(v_max, std::sqrt(v2));
  if (!fa) return;
  Real v[3], a[3];
  for (int n = 0; n < 3; ++n) v[n] = (pg[n] + x.w(m,CLW_DM+n,k,j,i))/d;
  AdvectCoeff(xa, m, k, j, i, v, a);
  for (int n = 0; n < 3; ++n) x.w(m,CLW_A1+n,k,j,i) = a[n];
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::WorstCell
//! \brief CLEAN WORST diagnostic of an unconverged backward-Euler solve: its kind from
//! the last two changes (STALL: fell less than 2x) and the cell of the largest |dE|/S

void RadiationVet::WorstCell(const Real dnl, const Real dnl0) {
  const bool ok = std::isfinite(dnl) && dnl0 >= 0.0 && cln_bon;
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int ni = id.nx1, nji = id.nx2*ni, nkji = id.nx3*nji, is = id.is, js = id.js;
  const int ks = id.ks;
  auto w = cln_w;
  auto ci = cln_wci;
  auto op = opac;
  Kokkos::MaxLoc<Real,int>::value_type ml;
  Kokkos::parallel_reduce("cln_worst", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          ok ? pmy_pack->nmb_thispack*nkji : 0),
  KOKKOS_LAMBDA(const int q, Kokkos::MaxLoc<Real,int>::value_type &v) {
    const int m = q/nkji, k = ks + (q%nkji)/nji, j = js + (q%nji)/ni, i = is + q%ni;
    const Real d = fabs(w(m,CLW_EK,k,j,i) - ci(m,0,k,j,i))/w(m,CLW_SE,k,j,i);
    if (d > v.val) {v.val = d; v.loc = q;}
  }, Kokkos::MaxLoc<Real,int>(ml));
  const Real gmx = GlobalReduce(ml.val, true);
  if (ok ? !(ml.loc >= 0 && ml.val >= gmx) : global_variable::my_rank != 0) return;
  const int q = std::max(ml.loc, 0), m = q/nkji, k = ks + (q%nkji)/nji;
  const int j = js + (q%nji)/ni, i = is + q%ni;
  DvceArray1D<Real> cv("cln_worst", 5);
  if (ok) {
    par_for("cln_worst_cell", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      cv(0) = w(m,CLW_EK,k,j,i);
      cv(1) = w(m,CLW_EO,k,j,i);
      cv(2) = w(m,CLW_T,k,j,i);
      cv(3) = w(m,CLW_SE,k,j,i);
      cv(4) = op(m,VET_OP_T,k,j,i);
    });
  }
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cv);
  std::printf("CLEAN WORST cycle %d %s last change %.3e (before %.3e) gid %d k %d j %d "
              "i %d E %.6e E_old %.6e T %.6e S %.3e chi_T %.3e\n",
              pmy_pack->pmesh->ncycle,
              !std::isfinite(dnl) ? "BREAKDOWN" : (dnl0 < 0.0 ? "?" : (dnl > 2.0*dnl0 ?
              "DIVERGE" : (dnl > 0.5*dnl0 ? "STALL" : "SLOW"))), dnl, dnl0,
              ok ? pmy_pack->pmb->mb_gid.h_view(m) : -1, k, j, i, h(0), h(1), h(2), h(3),
              h(4));
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::AdvectiveCoefficients
//! \brief a_d = v_d + D_db v_b, the advective flux per unit E, v the gas velocity of the
//! old vector plus the iterate's momentum change

void RadiationVet::AdvectiveCoefficients() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto w = cln_w;
  SolveArgs x = MakeArgs(), xa = x;
  // staggered on the cubed sphere: Cartesian (FaceFrameGasMomentum)
  const bool sc = x.cub;
  xa.cub = false;

  SolveFor("cln_advcoef", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = uh(m,IDN,k,j,i);
    Real v[3], a[3];
    if (sc) {
      FaceFrameGasMomentum(x, uh, m, k, j, i, v);
    } else {
      GasMomentum(x, uh, m, k, j, i, v);
    }
    for (int n = 0; n < 3; ++n) v[n] = (v[n] + w(m,CLW_DM+n,k,j,i))/d;
    AdvectCoeff(xa, m, k, j, i, v, a);
    for (int n = 0; n < 3; ++n) w(m,CLW_A1+n,k,j,i) = a[n];
  });
  // their ghosts are first read after the next exchange (the RHS's)
  GhostDefer(cln_w, CLW_A1, 3);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::StageOpacity
//! \brief clean_stag_opac_update: the active cells' opacities at the Newton iterate's gas
//! temperature (stage density), then the face coefficients and bracket that use them

void RadiationVet::StageOpacity() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;
  const int otype = opacity_type;
  const VetOpacTab ot = otab;
  const Real kp = kappa_p, kev = kappa_e, kf = kappa_f, ks_ = kappa_s;
  const Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  const bool osc = vscat;
  par_for("cln_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = uh(m,IDN,k,j,i);
    Real o[4];
    VetCellOpacities(otype, ot, eos, d, w(m,CLW_T,k,j,i), kp, kev, kf, ks_, rref, tref,
                     aa,
                    bb, o[0], o[1], o[2], o[3]);
    op(m,VET_OP_P,k,j,i) = d*o[0];
    op(m,VET_OP_E,k,j,i) = d*o[1];
    op(m,VET_OP_T,k,j,i) = d*(o[2] + o[3]);
    if (osc) op(m,VET_OP_S,k,j,i) = d*o[3];
  });
  // (cubed sphere: FaceCoefficients' own exchange)
  if (!pmy_pack->pmesh->use_cubed_sphere) {
    GhostFill(opac, 0, opac.extent_int(1), false);
  }
  FaceCoefficients();
  if (cln_consist) FaceBracket();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FirstIterate
//! \brief first Newton iterate: E of the old vector, no momentum change, gas T in local
//! exchange equilibrium with that E.  implicit_predictor: E + h K1_E instead (K1 the
//! slope of U^n), where that changes the radiation temperature by less than
//! implicit_predictor_tcap.  Only the starting point: the solution does not move.

void RadiationVet::FirstIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const bool feedback = cln_feedback;
  const Real hc = dt_sub*c_light, arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;
  auto k1 = t2k1;
  // sdirk2's stage 2 from (K1 + z K2^{n-1})/(1 + z), z = h (ch kE +
  // 4 c kP a T^4/e) the exchange stiffness: stiff cells follow the last step's stage-2
  // slope (the slow manifold's drift), slow ones the latest slope
  const bool pxc = t2_tab == 2 && t2_solve == VET_T2S_STAGE2;
  auto ks1 = t2k2, ks2 = t2k1;   // sdirk2 stage 2: K1 of this step, K2 of the last
  const Real ch = chat;
  const bool pred = (t2_ok || cln_k1ok) &&
                    t2_solve != VET_T2S_NONE && !cln_vzero && t2k1.extent_int(0) > nmb1;
  const Real h = dt_sub, ecap = SQR(SQR(kPredTcap)), rtol = cln_trtol;
  const Real cl = c_light;

  SolveFor("cln_iter0", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real eo = fmax(w(m,CLW_EO,k,j,i), 0.0);
    Real ep = pred ? eo + h*k1(m,VET_T2_E,k,j,i) : eo;
    if (pred && pxc) {
      const Real ei = feedback ? w(m,CLW_EI,k,j,i) : 0.0;
      const Real tk = (ei > 0.0) ? eos.Temperature(uh(m,IDN,k,j,i), ei) : 0.0;
      const Real z = h*ch*op(m,VET_OP_E,k,j,i) +
                     ((ei > 0.0) ? 4.0*hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(tk))/ei
                                 : 0.0);
      ep = eo + h*(ks1(m,VET_T2_E,k,j,i) + z*ks2(m,VET_T2_E,k,j,i))/(1.0 + z);
    }
    const Real e = (ep > eo/ecap && ep < eo*ecap) ? ep : eo;
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_RR,k,j,i) = 0.0;
    for (int n = 0; n < 3; ++n) w(m,CLW_DM+n,k,j,i) = 0.0;
    if (feedback) {
      w(m,CLW_EP,k,j,i) = 0.0;   // beta terms of this T (RightHandSide)
      const Real gm = 1.0 - 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                                 SQR(uh(m,IM3,k,j,i)))/SQR(uh(m,IDN,k,j,i)*cl);
      w(m,CLW_T,k,j,i) = ExchangeTemperature(eos, uh(m,IDN,k,j,i),
                             hc*op(m,VET_OP_P,k,j,i)*arad_*gm,
                             w(m,CLW_EI,k,j,i) + hc*op(m,VET_OP_E,k,j,i)*e,
                             w(m,CLW_T,k,j,i), rtol);
    }
  });
  // ghosts of E_k for the gas terms, with the next exchange
  if (feedback) GhostDefer(cln_w, CLW_EK, 1);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::NewtonUpdate
//! \brief new iterate E_k = X and the gas temperature in exchange equilibrium with it
//! (absorption of max(X, 0); also the iterate without clean_stag_gflux).  Returns
//! max |dE|/E, |dT|/T.

Real RadiationVet::NewtonUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  const Real hc = dt_sub*c_light, arad_ = arad;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  auto op = opac;
  auto w = cln_w;
  const bool poslim = cln_gflux;
  const bool gfd = cln_gflux && cln_feedback;
  // clean_conv = budget: changes over the rows' gross terms, and the time error
  // h/2 |f(Y) - f(U^n)| of the iterate (Heun against Euler, K1 = f(U^n))
  const bool bud = cln_bon;
  const bool kok = bud && (t2_ok || cln_k1ok) && t2_solve != VET_T2S_NONE &&
                   !cln_vzero && t2k1.extent_int(0) >= nmb;
  const Real h = dt_sub, rtol = cln_trtol;
  auto k1 = t2k1;

  Real change = 0.0, femax = 0.0, ermax = 0.0;
  // CUDA: 4 blocks of 256 per SM; the reduction keeps its blocks and grid (bitwise)
  using NewtonRp =
      Kokkos::RangePolicy<DevExeSpace, Kokkos::LaunchBounds<256, SolveLbMb(4)>>;
  AthReduce("cln_newton", NewtonRp(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &dmax, Real &fe_max, Real &er_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    k += ks;
    j += js;
    i += is;
    const Real e_old = w(m,CLW_EK,k,j,i);
    const Real t_old = w(m,CLW_T,k,j,i);
    Real x = w(m,CLW_X,k,j,i);
    const Real xr = poslim ? x : fmax(x, 0.0);   // the change: undamped
    const Real ex = fmax(x, 0.0);
    const Real e = poslim ? x : ex, g = w(m,CLW_GM,k,j,i);
    // gflux: the gas takes the Doppler part hc kE dop E of sb at X, as the E row in r
    // (lagged at E_k, T and E drifted apart pass by pass in moving opaque cells)
    if (gfd) w(m,CLW_EP,k,j,i) -= w(m,CLW_JV+1,k,j,i)*(x - e_old);
    const Real hcke = hc*op(m,VET_OP_E,k,j,i);
    const Real t = ExchangeTemperature(eos, uh(m,IDN,k,j,i),
        hc*op(m,VET_OP_P,k,j,i)*g*arad_,
        w(m,CLW_EI,k,j,i) + hcke*ex + w(m,CLW_EP,k,j,i), t_old, rtol);
    w(m,CLW_EK,k,j,i) = e;
    w(m,CLW_T,k,j,i) = t;
    const Real de = fabs(xr - e_old)/fmax(fmax(fabs(xr), fabs(e_old)), 1.0e-300);
    const Real dt = fabs(t - t_old)/fmax(fmax(t, t_old), 1.0e-300);
    if (bud) {
      const Real d = uh(m,IDN,k,j,i), se = w(m,CLW_SE,k,j,i);
      Real e1, e0, p, chi_rho, chi_t, cv;
      eos.ThermoAt(d, t, e1, p, chi_rho, chi_t, cv);
      eos.ThermoAt(d, t_old, e0, p, chi_rho, chi_t, cv);
      dmax = fmax(dmax, fmax(fabs(xr - e_old)/se, fabs(e1 - e0)/w(m,CLW_SG,k,j,i)));
      if (kok) er_max = fmax(er_max, 0.5*fabs(e - w(m,CLW_EO,k,j,i) -
                                              h*k1(m,VET_T2_E,k,j,i))/se);
    } else {
      dmax = fmax(dmax, fmax(de, dt));
    }
    // round-off level of E (its terms over E)
    fe_max = fmax(fe_max, (fabs(e)*(1.0 + hc*op(m,VET_OP_E,k,j,i)) +
        fabs(w(m,CLW_EO,k,j,i)) + hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(t)))/
        fmax(fabs(e), 1.0e-300));
  }, Kokkos::Max<Real>(change), Kokkos::Max<Real>(femax), Kokkos::Max<Real>(ermax));
  // one reduction over the ranks
  Real mx[3] = {change, femax, ermax};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, mx, 3, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  change = mx[0]; femax = mx[1]; ermax = mx[2];
  cln_fe = std::numeric_limits<Real>::epsilon()*femax;
  if (bud) {   // (the tolerance of the next test and solve; for settled, the level of a
    // change the linear solve fixes at the floor, its 0.1 tfl in the budget's norm)
    // (a step dt_retry redoes, whose state just failed, is solved strictly too)
    // (clean_clos1_budget, clean_bud_maxnorm: the time error down to
    // clean_tol_nl, as each relaxes the passes to it; at 10 clean_tol a step whose time
    // error is below, a 1e-6 wave's 1e-11 at a change of 1e-9, lost the wave's phase)
    const Real tfl = fmax(cln_tol_nl, 10.0*cln_tol);
    const Real tlo = (cln_c1b || cln_bmn) ? cln_tol_nl : tfl;
    cln_tolb = (kok && !dtr_redo) ?
               fmin(fmax(cln_bkap*ermax, tlo), cln_bmax) : tfl;
    cln_fe = 0.1*tfl;
  }
  // fwork: emission hc kP a T^4 of the new T for the face beta terms of the velocity
  // update's F solve
  if (cln_fwork && cln_feedback) {
    par_for("cln_sbre", DevExeSpace(), 0, nmb-1, ks, ks+nx3-1, js, js+nx2-1, is,
            is+nx1-1, KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      w(m,CLW_RE,k,j,i) = hc*op(m,VET_OP_P,k,j,i)*arad_*SQR(SQR(w(m,CLW_T,k,j,i)));
    });
    GhostDefer(cln_w, CLW_RE, 1);
  }
  // staggered: the ghosts of the new E_k, for the velocity update and the next RHS
  if (cln_feedback) GhostFill(cln_w, CLW_EK, 1);
  return change;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VelocityUpdate
//! \brief momentum change of the new iterate less force_reference's reference force, and
//! the new velocity's advective coefficients.  Returns max |dv| over max |v|.

Real RadiationVet::VelocityUpdate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  CsLateralStress(CLW_EK);   // E_k's ghosts: NewtonUpdate
  SolveArgs x = MakeArgs(), xh = x;
  xh.hi = cln_split;
  xh.mcl = true;
  if (cln_mclc.extent_int(0) != 1) Kokkos::realloc(cln_mclc, 1);
  Kokkos::deep_copy(cln_mclc, 0.0);
  xh.mclc = cln_mclc;
  // the F solve's face tie after the drag (same face T: CLF_PT)
  const bool td = cln_fwork;
  // clean_stag_fwork: the F solve's momentum, which the write-back gives the gas (work W)
  if (cln_fwork) {
    FluxSolve(CLW_EK, !td, false);
    xh.fs = true;
    cln_fsk = true;
  }
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto op = opac;
  const Real arad_ = arad;
  const int fref = force_ref;
  const bool cap = (fref == VET_FREF_CLEAN) && !cln_fref_set && !cln_vzero;
  auto fr = cln_fref;
  auto aref = arad_ref;

  // with the F solve's F: the faces' terms and drag first, once per face, then the cells'
  // momentum and advective coefficients
  const bool fa = xh.fs;
  if (fa) FaceDrag(false, true, false);   // lateral sums of E_k: above
  if (td) FluxTie(CLW_EK, (t2_solve == VET_T2S_NONE) ? 0 : ((t2_solve ==
                  VET_T2S_STAGE1) ? CLF_KD + 1 : CLF_KD), true);
  xh.fa = fa;
  auto xa = x;
  xa.cub = false;
  Real dvmax = 0.0, vmax = 0.0, fvmax = 0.0;
  CLEAN_HIPX(const SolveArgs *px = DevArgs(x, 0), *pxh = DevArgs(xh, 1);
             const SolveArgs *pxa = DevArgs(xa, 4);)
  // fa: CellMomentumChange<true> reads the faces' terms, in a kernel of its own (HIP: 256
  // VGPRs, 2 waves per SIMD)
#if defined(KOKKOS_ENABLE_HIP)
  constexpr int vmb = 2;
#else
  constexpr int vmb = SolveLbMb(2);
#endif
  auto cell = KOKKOS_LAMBDA(const int idx, int &m, int &k, int &j, int &i) {
    m = idx/nkji;
    k = (idx - m*nkji)/nji + ks;
    j = (idx - m*nkji - (k - ks)*nji)/nx1 + js;
    i = idx - m*nkji - (k - ks)*nji - (j - js)*nx1 + is;
  };
  if (fa) {
    AthReduce("cln_velocity", Kokkos::RangePolicy<DevExeSpace,
              Kokkos::LaunchBounds<256, vmb>>(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
      CLEAN_HIPX(const SolveArgs &x = *px, &xh = *pxh, &xa = *pxa;)   // device memory
      int m, k, j, i;
      cell(idx, m, k, j, i);
      VelocityCell<true>(x, xh, xa, uh, fref, fr, aref, cap, op, arad_, true, m,
                         k, j, i, dv_max, v_max, fv_max);
    }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  } else {
    AthReduce("cln_velocity", SolveRpLb2(DevExeSpace(), 0, nmb*nkji),
    KOKKOS_LAMBDA(const int idx, Real &dv_max, Real &v_max, Real &fv_max) {
      CLEAN_HIPX(const SolveArgs &x = *px, &xh = *pxh, &xa = *pxa;)   // device memory
      int m, k, j, i;
      cell(idx, m, k, j, i);
      VelocityCell<false>(x, xh, xa, uh, fref, fr, aref, cap, op, arad_, false, m,
                          k, j, i, dv_max, v_max, fv_max);
    }, Kokkos::Max<Real>(dvmax), Kokkos::Max<Real>(vmax), Kokkos::Max<Real>(fvmax));
  }
  ReduceMax(dvmax, vmax, fvmax);
  if (fa) {
    GhostDefer(cln_w, CLW_A1, 3);   // as AdvectiveCoefficients
  } else {
    FaceDrag(false, false, false);   // this dm's coefficients
    AdvectiveCoefficients();
  }
  vmax = fmax(vmax, 1.0e-300);
  cln_fv = std::numeric_limits<Real>::epsilon()*fvmax/vmax;
  return dvmax/vmax;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CopyIterate
//! \brief the solution is the last Newton iterate

void RadiationVet::CopyIterate() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto w = cln_w;
  SolveFor("cln_copy", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    w(m,CLW_X,k,j,i) = w(m,CLW_EK,k,j,i);
  });
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::CountNegative
//! \brief the number of cells with E < 0 beyond the solver tolerance (10 tol max|B|);
//! WriteBack raises those within it to zero, counted

int RadiationVet::CountNegative() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto w = cln_w;

  Real bmax = 0.0;
  Kokkos::parallel_reduce("cln_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(w(m,CLW_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));

  const Real zero = 10.0*cln_tol*GlobalReduce(bmax, true);
  Real nneg = 0.0;
  Kokkos::parallel_reduce("cln_neg", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &n) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    if (!(w(m,CLW_X,k+ks,j+js,i+is) >= -zero)) n += 1.0;
  }, Kokkos::Sum<Real>(nneg));
  return static_cast<int>(GlobalReduce(nneg, false));
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::CountUnrealizable
//! \brief staggered: realizability of the solution X beyond the solver tolerance (10
//! clean_tol max|B|).  R1: cells with E < 0 (cln_nr1); R2: interior faces whose flux
//! exceeds c E of the donor (cln_nr2).  Returns cln_nr1 + cln_nr2.

int RadiationVet::CountUnrealizable() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, js = indcs.js, ks = indcs.ks;
  int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  int nmb = pmy_pack->nmb_thispack;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  GhostFill(cln_w, CLW_X, 1);
  CsLateralStress(CLW_X);
  SolveArgs x = MakeArgs();
  Real bmax = 0.0;
  Kokkos::parallel_reduce("cln_bmax", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
  KOKKOS_LAMBDA(const int idx, Real &b_max) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = idx - m*nkji - k*nji - j*nx1;
    b_max = fmax(b_max, fabs(x.w(m,CLW_B,k+ks,j+js,i+is)));
  }, Kokkos::Max<Real>(bmax));
  const Real zero = 10.0*cln_tol*GlobalReduce(bmax, true);
  Real n1 = 0.0, n2 = 0.0;
  CLEAN_HIPX(const SolveArgs *px = DevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    const int f1 = nx1 + di, f2 = nx2 + dj, f3 = nx3 + dk, fji = f2*f1, fkji = f3*fji;
    Real r1 = 0.0, r2 = 0.0;
    Kokkos::parallel_reduce("cln_real", SolveRpLb2(DevExeSpace(), 0, nmb*fkji),
    KOKKOS_LAMBDA(const int idx, Real &q1, Real &q2) {
      CLEAN_HIPX(const SolveArgs &x = *px;)   // in device memory: no per-thread copy
      int m = idx/fkji;
      int k = (idx - m*fkji)/fji;
      int j = (idx - m*fkji - k*fji)/f1;
      int i = idx - m*fkji - k*fji - j*f1;
      const bool cell = (i < nx1 && j < nx2 && k < nx3);   // its low face
      k += ks;
      j += js;
      i += is;
      const Real el = x.w(m,CLW_X,k-dk,j-dj,i-di), er = x.w(m,CLW_X,k,j,i);
      if (dir == 0 && cell && !(er >= -zero)) q1 += 1.0;
      if (IsBoundaryFace(x, dir, m, k, j, i)) return;
      const Real f = FaceFlux(x, CLW_X, dir, m, k, j, i, true);
      x.fc[dir](m,CLF_PT,k,j,i) = f;   // for WriteBack, next (QL is free)
      const Real ed = (f > 0.0) ? el : er;
      const bool r2 = !(fabs(f) <= (1.0 + 1.0e-12)*x.c*ed + x.c*zero);
      if (r2) q2 += 1.0;
    }, Kokkos::Sum<Real>(r1), Kokkos::Sum<Real>(r2));
    n1 += r1;
    n2 += r2;
  }
  Real nr[2] = {n1, n2};
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, nr, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  cln_nr1 = static_cast<int>(nr[0]);
  cln_nr2 = static_cast<int>(nr[1]);
  return cln_nr1 + cln_nr2;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FaceBracket
//! \brief staggered clean_consistency: the bracket C_f = w (F_ray,f - Phi_f(J)) of every
//! interior face, F_ray the rays' energy flux (CLN_FX, CLN_F2, CLN_F3 of the upper cell),
//! Phi_f(J) the flux of the rays' own J, so Phi(J) + C = (1 - w) Phi(J) + w F_ray.
//! Lagged (fixed in the solve), energy flux only.  It enters the streaming share as flux
//! factors, Phi_w = c (g+ E_L + g- E_R) with
//!   c (g+ J_L + g- J_R) w = c (h+ J_L + h- J_R) w + C_f  (CLF_HP, CLF_HM):
//! the half ranges H+-_f = h+-_f J tilted by dd = C_f/(w c) in proportion to their sizes,
//!   H+-' = H+-_f (1 +- dd/s),  s = H+_f - H-_f,
//! the shrinking one at most to 0.  Linear in E with upwind signs, so the M-matrix holds
//! for any E.  g+- are not capped at 1 (a beam edge's face flux can exceed c J_donor).
//! Cubed sphere: a panel seam's high x2, x3 faces take the rays' flux from their last cell
//! (CLN_FH, RayFace).  clean_stag_thinf takes the flux of the cells' half ranges there
//! (the face sum with fixed sets): the two MeshBlocks of a face mix or extrapolate its
//! sum with their own cells' weights, and the write-back's face tie would move the
//! difference into one side's E.

void RadiationVet::FaceBracket() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  SolveArgs x = MakeArgs();
  x.feedback = false;   // the gas terms read the last iterate's work vectors
  const bool thf = cln_thinf;
  const bool lat = thf && pmy_pack->pmesh->use_cubed_sphere;
  const Real rmn = cln_raymin;
  // clean_stag_thinf: R, what the rays hold of the cell: its J and its faces' fluxes
  SolveFor("cln_rayj", nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x.w(m,CLW_Z,k,j,i) = x.cs(m,CLN_J,k,j,i);
    if (!thf) return;
    Real r = fmax(x.cs(m,CLN_JR,k,j,i), 0.0);
    for (int n = 0; n < x.nf; ++n) {
      for (int sd = 0; sd < 2; ++sd) {
        const int kf = k + sd*(n == 2), jf = j + sd*(n == 1), iff = i + sd*(n == 0);
        if (sd > 0 && IsBoundaryFace(x, n, m, kf, jf, iff)) continue;
        r = fmax(r, fabs(RayFace(x, n, m, kf, jf, iff)));
      }
    }
    x.w(m,CLW_Y,k,j,i) = r;
  });
  GhostFill(cln_w, thf ? CLW_Y : CLW_Z, thf ? 2 : 1);
  CsLateralStress(CLW_Z);
  CLEAN_HIPX(const SolveArgs *px = DevArgs(x, 0);)
  for (int dir = 0; dir < x.nf; ++dir) {
    const int di = (dir == 0), dj = (dir == 1), dk = (dir == 2);
    SolveFor("cln_sbracket", nmb1, ks, ke+dk, js, je+dj, is, ie+di,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      CLEAN_HIPX(const SolveArgs &x = *px;)   // in device memory: no per-thread copy
      // Phi(J) at its steady value, -(c/chi_f) div(D J), in the share 1 - w, so the
      // bracket (and the scheme's steady state) does not depend on h or G.  Where the
      // diffusion flux exceeds free streaming, lambda = |div(D J)|/(chi_f J_f) > 1, it
      // takes c J_f along -div(D J) (Larsen's max limiter)
      const Real w = x.fc[dir](m,CLF_W,k,j,i);
      const Real chi = FaceChi(x.op, x.x1f, x.mbc, dir, m, k, j, i, is, ie);
      const Real st = FaceStress(x, CLW_Z, dir, m, k, j, i, true);
      const Real jf = 0.5*(x.w(m,CLW_Z,k,j,i) + x.w(m,CLW_Z,k-dk,j-dj,i-di));
      const Real cd = (jf > 0.0 && fabs(st) > chi*jf) ? x.c*jf/fabs(st) :
                      ((chi > 0.0) ? x.c/chi : x.fc[dir](m,CLF_TH,k,j,i)*x.h*x.c*x.c);
      // cells' half ranges realisable, h+ in [0, 1], h- in [-1, 0]; the rest is DG below
      x.fc[dir](m,CLF_HP,k,j,i) = fmin(fmax(x.fc[dir](m,CLF_HP,k,j,i), 0.0), 1.0);
      x.fc[dir](m,CLF_HM,k,j,i) = fmax(fmin(x.fc[dir](m,CLF_HM,k,j,i), 0.0), -1.0);
      Real gl = 0.0, gr = 0.0;
      AreaStream(x, dir, m, k, j, i, w, gl, gr);
      const Real pj = -(1.0 - w)*cd*st +
          w*x.c*(x.fc[dir](m,CLF_HP,k,j,i)*x.w(m,CLW_Z,k-dk,j-dj,i-di) +
                 x.fc[dir](m,CLF_HM,k,j,i)*x.w(m,CLW_Z,k,j,i)) +
          gl*x.w(m,CLW_Z,k-dk,j-dj,i-di) + gr*x.w(m,CLW_Z,k,j,i);
      const Real fr0 = RayFace(x, dir, m, k, j, i);
      const bool bf = IsBoundaryFace(x, dir, m, k, j, i);
      if (bf) return;
      // the rays' J also below clean_ray_min E (a diagonal beam's side cells);
      const Real jl = x.cs(m,CLN_JR,k-dk,j-dj,i-di), jr = x.cs(m,CLN_JR,k,j,i);
      const Real jd = (fr0 > 0.0) ? x.cs(m,CLN_JR,k-dk,j-dj,i-di) : x.cs(m,CLN_JR,k,j,i);
      // the rays' face flux interpolated linearly in optical depth,
      // H_f = (tau_R H_L + tau_L H_R)/(tau_L + tau_R), tau the half-cell depths (H r^2 on
      // sp / cs radial faces; x1 faces there, every face on Cartesian).  F_ray (x
      // (r_d/r_f)^2 on sp / cs) keeps the share e = (1 - min/max of the cells' r^2 H)
      // e^-tau_d, so a jump in H (beam edge, shadow) keeps the transport solution and
      // smooth or opaque regions the interpolation; a convex blend, single-valued per
      // face.
      Real fy = fr0;
      if (lat && dir > 0) {
        fy = x.fc[dir](m,CLF_HP,k,j,i)*jl + x.fc[dir](m,CLF_HM,k,j,i)*jr;
      }
      if (dir == 0 || !(x.sph || x.cub)) {
        const int il = i - di;
        auto hc = [&](const int kk, const int jj, const int ii) {
          return (x.cs(m,CLN_HP+dir,kk,jj,ii) + x.cs(m,CLN_HM+dir,kk,jj,ii))*
                 x.cs(m,CLN_JR,kk,jj,ii);
        };
        Real wl, wr, tg[4];
        HcTauW(x, dir, m, k, j, i, wl, wr, tg);
        const Real hl = hc(k-dk, j-dj, il), hr = hc(k, j, i);
        const Real al = fabs(tg[2]*hl), ar = fabs(tg[3]*hr);
        const Real e = ((hl*hr > 0.0) ? 1.0 - fmin(al, ar)/fmax(al, ar) : 1.0)*
                       exp(-tg[(fr0 > 0.0) ? 0 : 1]);
        const Real ga = (tg[(fr0 > 0.0) ? 2 : 3] - 1.0)*fy;
        fy = (1.0 - e)*(wl*hl + wr*hr) + e*(fy + ga);
      }
      // clean_stag_thinf: the rest of the rays' face flux, dd = F_y - h+ J_L - h- J_R,
      // tilted into the face half ranges in proportion to their sizes (the shrinking one
      // at most to 0), per J of the donor (at least min(clean_ray_min E_old, R)), g+ in
      // [0, 1], g- in [-1, 0]: the thin share w c (g+ E_L + g- E_R) is within [-c E_R,
      // c E_L], realisable by construction, no DG.  A flux the cap cuts is carried by a
      // higher E (a diagonal beam's side cells: J ~ 0, R ~ E).  In a ray hole (R <<
      // clean_ray_min E_old) E streams along the rays' own direction: per clean_ray_min
      // E_old, g -> 0 as J -> 0, so E would stay, its flux independent of it.
      if (thf) {
        const Real hl = x.fc[dir](m,CLF_HP,k,j,i)*jl, hr = x.fc[dir](m,CLF_HM,k,j,i)*jr;
        const Real dd = fy - (hl + hr), fr = hl + hr + dd, s = hl - hr;
        const Real tl = hl + dd*hl/s, tr = hr - dd*hr/s;
        const bool cut = (tl < 0.0) || (tr > 0.0);
        const Real hp = !(s > 0.0) ? hl : (cut ? ((tl < 0.0) ? 0.0 : fmax(fr, 0.0)) : tl);
        const Real hm = !(s > 0.0) ? hr : (cut ? ((tl < 0.0) ? fmin(fr, 0.0) : 0.0) : tr);
        const Real el = fmax(x.w(m,CLW_EO,k-dk,j-dj,i-di), 0.0);
        const Real er = fmax(x.w(m,CLW_EO,k,j,i), 0.0);
        const Real dl = fmax(jl, fmin(rmn*el, x.w(m,CLW_Y,k-dk,j-dj,i-di)));
        const Real dr = fmax(jr, fmin(rmn*er, x.w(m,CLW_Y,k,j,i)));
        const Real gp = fmin(hp/dl, 1.0), gm = fmax(hm/dr, -1.0);
        if (jl > 0.0) x.fc[dir](m,CLF_HP,k,j,i) = gp;
        if (jr > 0.0) x.fc[dir](m,CLF_HM,k,j,i) = gm;
        x.fc[dir](m,x.fdg,k,j,i) = 0.0;
        return;
      }
      // default form: half ranges h+- realisable, the rest the lagged flux DG/w = q c
      // (F_y - h+ J_L - h- J_R) (ThinFaceDG), F_y clamped to [-J_R, J_L], q = r^2
      // the quasi-steadiness gate, r the worse of the two cells' min(E_old/J, J/E_old)
      auto rq = [&](const Real a, const Real b) {
        return (fmax(a, b) > 0.0) ? fmin(a, b)/fmax(a, b) : 1.0;
      };
      const Real r = fmin(rq(fmax(x.w(m,CLW_EO,k-dk,j-dj,i-di), 0.0), fmax(jl, 0.0)),
                          rq(fmax(x.w(m,CLW_EO,k,j,i), 0.0), fmax(jr, 0.0)));
      const Real hup = x.fc[dir](m,CLF_HP,k,j,i)*jl + x.fc[dir](m,CLF_HM,k,j,i)*jr;
      x.fc[dir](m,x.fdg,k,j,i) = r*r*x.c*(fmin(fmax(fy, -jr), jl) - hup);
    });
  }
}

//----------------------------------------------------------------------------------------
//! \fn SolvePass RadiationVet::PassSetup
//! \brief what the passes of this stage solve iterate on, their cap and tolerances

SolvePass RadiationVet::PassSetup() {
  SolvePass p;
  // Newton iteration on the gas coupling, Picard passes on the lagged slopes and (with
  // clean_closure = iterate, at most clean_maxpic) the closure; one linear solve
  // without these.  clean_closure = extrap(2): the stage-1 solve takes the closure of
  // its own first iterate once (pass 1), else the step is first order in time.
  p.clos1 = (cln_extrap > 0) && (t2_solve == VET_T2S_STAGE1);
  p.iterate = cln_feedback || cln_closit || p.clos1;
  p.t2stage = (t2_solve == VET_T2S_STAGE1 || t2_solve == VET_T2S_STAGE2);
  p.converged = !p.iterate;
  p.nmax = p.iterate ? cln_maxnl : 1;
  // lagged slopes settle only to the linear tolerance.  clean_conv = budget: the pass's
  // tolerance (NewtonUpdate), the first solve's from the floor (restarts bitwise)
  p.tol_nl = fmax(cln_tol_nl, 10.0*cln_tol);
  cln_tolb = p.tol_nl;
  // budget norm and time error only in gas-coupled Newton passes with K1 = f(U^n) known;
  // otherwise the strict test
  cln_bon = cln_budget && p.iterate && cln_feedback && (t2_ok || cln_k1ok) &&
            t2_solve != VET_T2S_NONE && !cln_vzero &&
            t2k1.extent_int(0) >= pmy_pack->nmb_thispack;
  cln_ewr = -1.0;
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::FirstPass
//! \brief the first iterate: E, T, opacity, beta, drag, donor fluxes; resets of the
//! solve's passes

void RadiationVet::FirstPass() {
  const bool opu = cln_opt && cln_feedback;
  FirstIterate();
  cln_fsk = false;
  cln_csyok = false;   // lateral sums' coefficients once per closure
  if (opu) StageOpacity();
  if (cln_feedback) AdvectiveCoefficients();
  // the drag's face flux from the first iterate's beta and emission (a restart, failed
  // stage or discarded attempt does not carry the last solve's)
  if (cln_feedback) RightHandSide(true);
  if (cln_feedback) FaceDrag(true);
  // clean_stag_xdonor: q of the first iterate with its drag share b (CLW_L holds E_k's
  // lateral sums after the first drag)
  FaceDonor(true, !cln_feedback);
  cln_gpic = false;   // clean_gas_newton_switch: the solve's history
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::PassClosure
//! \brief pass p.n: closure, opacities and donor fluxes of the iterate, in the passes
//! that take them

void RadiationVet::PassClosure(SolvePass &p) {
  const bool opu = cln_opt && cln_feedback;
  p.clos = (cln_closit && p.n < cln_maxpic) || (p.clos1 && p.n == 1);
  // opacities of the iterate's T, in passes < clean_stag_opac_passes
  p.opr = opu && p.n > 0 && p.n < cln_optn;
  if (p.clos) {
    cln_csyok = false;
    ClosureIterate();
    if (!p.opr) FaceCoefficients();   // (else in StageOpacity)
    if (cln_consist && !p.opr) FaceBracket();
    if (cln_feedback) AdvectiveCoefficients();
  }
  if (p.opr) {
    StageOpacity();
    FaceDrag(false);
    AdvectiveCoefficients();
  }
  // clean_stag_xdonor: q of the iterate; CLW_L holds E_k's lateral sums unless the
  // closure, bracket or opacity refresh rewrote them (and on the cubed sphere)
  if (p.n > 0) FaceDonor(false, p.clos || p.opr || !cln_feedback,
                         cln_xdfz > 0 && p.n >= cln_xdfz);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::PassRows
//! \brief a pass's linear system: right-hand side, operator rows, preconditioner

void RadiationVet::PassRows(const SolvePass &p) {
  // after the xdonor: the lateral sums of E_k are its
  RightHandSide(false, !(p.n > 0 && cln_xd_on));
  // the first iterate's T from the pass's gas row, then its RHS
  if (p.n == 0 && cln_feedback) {
    NewtonUpdate();
    RightHandSide();
  }
  PreconditionerSetup();
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::PassSolve
//! \brief a pass's linear solve for E; false: the passes end here (no iteration, or a
//! failed linear solve)

bool RadiationVet::PassSolve(SolvePass &p) {
  // budget: a pass whose iterate gives the next closure is solved in the strict norm
  const bool cnext = (cln_closit && p.n + 1 < cln_maxpic) || (p.clos1 && p.n == 0);
  if (cnext) cln_ewr = (cln_c1b && !dtr_redo) ? -3.0 : -2.0;
  p.nit = SolveBiCGStab();
  if (cnext) cln_ewr = -1.0;
  p.nit_sum += std::abs(p.nit);
  ++p.nsolve;
  cln_nsolve += 1.0;
  cln_nit += std::abs(p.nit);
  if (p.nit < 0) p.dnl = std::numeric_limits<Real>::infinity();
  // implicit_lin_maxit: a linear solve that ends at its cap leaves its iterate to the
  // pass, which is then not the last one
  p.lp = p.nit < 0 && p.n + 1 < p.nmax;
  return p.iterate && (p.nit >= 0 || p.lp);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::PassUpdate
//! \brief a pass's updates: E and gas temperature (change p.de), then face fluxes, gas
//! momentum and drag (change p.dv)

void RadiationVet::PassUpdate(SolvePass &p) {
  if (!p.t2stage && p.n + 1 == p.nmax) {   // for WorstCell
    auto sv = Kokkos::subview(cln_w.s, Kokkos::ALL, std::make_pair(CLW_EK,
                              CLW_EK + 1), Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    Kokkos::realloc(cln_wci, sv.extent(0), 1, sv.extent(2), sv.extent(3),
                    sv.extent(4));
    Kokkos::deep_copy(cln_wci, sv);
  }
  p.de = NewtonUpdate();
  if (cln_bon) p.tol_nl = cln_tolb;
  p.dv = 0.0;
  if (cln_feedback) p.dv = VelocityUpdate();
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationVet::PassSettled
//! \brief the stop test of the passes; also the stall switch clean_gas_newton_switch

bool RadiationVet::PassSettled(SolvePass &p) {
  p.dnl0 = (p.n + 1 == p.nmax) ? p.dnl : -1.0;
  p.dnl = fmax(p.de, p.dv);
  // clean_gas_newton_switch = W: a stalled solve (pass + 1 >= _min and > W, the
  // minimum over the last W passes >= clean_stall_fac x the one before, up at least
  // twice in the window) drops the gas Newton slope for its remaining passes: the gas
  // T takes the bracketed root of its row, a contraction (Newton cycles where c_v
  // jumps, e.g. H ionisation).  Deterministic in the passes: restarts bitwise.
  if (!cln_gpic) {
    p.rhist.push_back(p.dnl);
    constexpr int nw = 4;
    if (p.n + 1 >= 6) {
      Real mwin = 1.0e300, mpre = 1.0e300;
      int nup = 0;
      for (int q = 0; q <= p.n; ++q) {
        if (q > p.n - nw) {
          mwin = fmin(mwin, p.rhist[q]);
          if (p.rhist[q] > p.rhist[q-1]) ++nup;
        } else {
          mpre = fmin(mpre, p.rhist[q]);
        }
      }
      if (mwin >= 0.5*mpre && nup >= 2) {
        cln_gpic = true;
        cln_ngnsw += 1.0;
        if (global_variable::my_rank == 0) {
          std::printf("CLEAN GNSW cycle %d t2s %d pass %d: stalled, the gas Newton "
                      "slope off (%.0f so far)\n", pmy_pack->pmesh->ncycle, t2_solve,
                      p.n, cln_ngnsw);
        }
      }
    }
  }
  // strict test: stop at 0.1 x the first pass's change, at most 1e-4
  if (!cln_bon && p.n == 0) p.trel = fmin(0.1*p.dnl, 1.0e-4);
  p.tol_nl = fmax(p.tol_nl, p.trel);
  // clean_verbose >= 2: each pass's changes against the tolerance it is tested with
  if (cln_verbose > 1 && global_variable::my_rank == 0) {
    std::printf("CLEAN PASS cycle %d t2s %d pass %d it %d de %.3e dv %.3e tol %.3e "
                "fe %.3e fv %.3e%s%s\n", pmy_pack->pmesh->ncycle, t2_solve, p.n,
                p.nit, p.de, p.dv, p.tol_nl, cln_fe, cln_fv, p.clos ? " clos" : "",
                p.opr ? " opac" : "");
  }
  if (p.lp) {
    p.de0 = p.de;
    p.dv0 = p.dv;
    return false;
  }
  if (p.Settled(p.de, p.de0, cln_fe) && p.Settled(p.dv, p.dv0, cln_fv) &&
      !(p.clos1 && p.n == 0)) {
    p.converged = true;
    return true;
  }
  p.de0 = p.de;
  p.dv0 = p.dv;
  return false;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SolveEnd
//! \brief the solution of the passes: its iterate, realizability and positivity counts

void RadiationVet::SolveEnd(SolvePass &p) {
  if (p.iterate) CopyIterate();
  // R1, R2 of the solution, counted
  CountUnrealizable();
  p.nfirst = {cln_nr1, cln_nr2};
  if (!p.converged) {
    p.nit = -std::abs(p.nit) - 1;
  }

  // positivity of the solution
  p.nneg = CountNegative();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VerdictSolve
//! \brief what an unconverged solve means (a failed stage, dt_retry, FATAL or a warning),
//! and the CLEAN line of the solve

void RadiationVet::VerdictSolve(const SolvePass &p) {
  if (p.nit < 0) cln_nfail += 1.0;
  if (p.nneg > 0) cln_nneg += 1.0;
  // an unconverged backward-Euler solve fails the cycle: with <time>/dt_retry it is
  // retaken at a smaller dt (a failed stage is redone by backward Euler via t2_fail);
  // without, the run stops if the last change exceeds implicit_resid_fatal
  if (p.nit < 0 && !p.t2stage) {
    const bool fatal = !dtr_on && impl_res_fatal > 0.0 && !(p.dnl <= impl_res_fatal);
    if (dtr_on) {
      if (!dtr_fail) {
        dtr_why = 1;
        dtr_resid = -1.0;
      }
      dtr_fail = true;
    }
    if (global_variable::my_rank == 0) {
      std::printf("CLEAN WARNING cycle %d: backward-Euler solve not converged (%s, %d "
                  "so far; last change %.3e)\n", pmy_pack->pmesh->ncycle,
                  dtr_on ? "dt_retry" : (fatal ? "no dt_retry: > implicit_resid_fatal" :
                  "no dt_retry: accepted"), static_cast<int>(cln_nfail), p.dnl);
    }
    WorstCell(p.dnl, p.dnl0);
    if (fatal) RadVetFatal("<rad_vet> clean: a backward-Euler solve did not converge and "
                           "its last change exceeds <rad_vet>/implicit_resid_fatal");
  }
  if (global_variable::my_rank == 0) {
    std::printf("CLEAN cycle %d t2s %d it %d (sum %d) pos 0 mark 0 all 0 neg %d "
                "r1 %d r2 %d first %d %d guard %d solves %d\n",
                pmy_pack->pmesh->ncycle, t2_solve, p.nit, p.nit_sum, p.nneg, cln_nr1,
                cln_nr2, p.nfirst[0], p.nfirst[1], cln_nguard, p.nsolve);
  }
  if (cln_mclc.extent_int(0) == 1) {   // last pass's face momentum terms at the cap
    auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cln_mclc);
    const Real n = GlobalReduce(h(0), false);
    if (global_variable::my_rank == 0 && cln_verbose > 0 && n > 0.0) {
      std::printf(
          "CLEAN MCAP cycle %d t2s %d: %.0f face momentum terms at the cap h chi_f "
          "E in the last pass\n", pmy_pack->pmesh->ncycle, t2_solve, n);
    }
  }
}

}  // namespace radvet
