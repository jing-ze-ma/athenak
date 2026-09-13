//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction_transverse.cpp
//! \brief the IMPLICIT (unconditionally stable) TRANSVERSE radiative diffusion operator,
//! <hydro>/ or <mhd>/rad_implicit_ang.  Operator-split: runs after the RK update and the
//! source terms, and before the ghost exchange, so what it writes is what is
//! communicated.  Only u0(IEN) is touched.
//!
//! WHY.  The transverse (x2/x3) radiative operator is limited by
//! dt ~ dx^2 rho c_v/kappa_rad, and in the layers a radiative box actually cares about
//! that limit is orders of magnitude below the hydrodynamic step.  rad_cap_ang keeps the
//! run alive by scaling every transverse face down to what an explicit step can carry,
//! which is stable and exactly conservative but physically wrong: it throttles the
//! horizontal radiative exchange itself.  This operator removes the constraint instead
//! of the physics.
//!
//! WHAT IS SOLVED.  With e the internal energy density and T linearised about the frozen
//! pre-update state w0,  T_i = T*_i + alpha_i (e_i - e*_i),  alpha_i = 1/(rho_i c_v,i),
//! and the frozen face coefficient C_f = A_f K_f/dl_f of BuildAngularCoeffs (the SAME
//! RadFaceKappa, flux limiter and tau-blend weight the explicit x2/x3 face flux uses),
//! the increment y_i = e_i - e*_i obeys the linear ODE
//!     dy_i/dt = M(y)_i = (1/V_i) sum_f s_f C_f (T_j - T_i),   y_i(0) = 0,
//! integrated over the stage, t = 0 -> beta_dt.  M is affine, and the matrix
//!     dM_i/dy_j = C_f alpha_j/V_i  (off-diagonal),  -alpha_i sum_f C_f/V_i (diagonal)
//! is an M-matrix, so the exact solution of the ODE cannot overshoot.
//!
//! HOW IT IS SOLVED (step 1).  By an RKL1 super-time-stepping loop (Meyer, Balsara &
//! Aslam 2012) rather than a linear solve: s substages of the same 5-point stencil, with
//! the shifted-Legendre coefficients
//!     w1 = 2/(s^2+s),  mu_j = (2j-1)/j,  nu_j = (1-j)/j,  mu~_j = mu_j w1,
//!     Y_0 = 0,  Y_1 = Y_0 + mu~_1 tau M(Y_0),
//!     Y_j = mu_j Y_{j-1} + nu_j Y_{j-2} + mu~_j tau M(Y_{j-1}),   y = Y_s,
//! whose amplification polynomial is L_s(1 + w1 z) and is therefore stable for
//! |lambda| tau <= 2/w1 = s^2+s, i.e. for a super-step of (s^2+s)/2 explicit steps.  s is
//! chosen per call from a GLOBAL reduction of the Gershgorin row radius
//!     z_i = tau (alpha_i sum_f C_f + sum_f C_f alpha_j)/V_i  >=  lambda_max tau,
//! as s = ceil((sqrt(1 + 8R) - 1)/2) with R = 0.55 max_i z_i (the 1.1 is round-off
//! margin), so the cost is sqrt of the stiffness ratio, not the ratio.  A direct
//! (PCG) solve of the backward-Euler system is the intended step 2; the machinery here
//! -- the frozen coefficients, the dedicated one-variable halo exchange, the flux-form
//! stencil and the conservation check -- is the same either way.
//!
//! CONSERVATION.  The stencil is written in FLUX form: the two cells sharing a face
//! evaluate the same product C_f (Th_j - Th_i) from the same operands, so they are
//! bitwise equal and opposite and sum_i V_i y_i = 0 to round-off over the whole mesh
//! (the report below measures exactly that).  A face on a physical (non-periodic) x2/x3
//! boundary carries no flux, since the increment is not defined outside the mesh.

#include <float.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "utils/eint_from_cons.hpp"
#include "conduction.hpp"

//----------------------------------------------------------------------------------------
//! \fn void Conduction::ImplicitTransverseUpdate
//! \brief advance u0(IEN) by the transverse radiative diffusion over one stage; see the
//! file comment for the system and the scheme.  A no-op unless rad_implicit_ang is set.

void Conduction::ImplicitTransverseUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                                          const Real beta_dt) {
  if (!rad_implicit_ang) return;
  if (rad_sts_all) return;    // the unified operator owns the transverse faces instead
  RklConductionUpdate(u0, eos, beta_dt, false);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::StsConductionUpdate
//! \brief advance u0(IEN) by the UNIFIED radiative conduction operator -- x1, x2 and x3
//! in ONE RKL1 loop -- over one stage (<hydro>/ or <mhd>/rad_sts_all).
//!
//! WHY ONE LOOP.  The split of an implicit tridiagonal x1 solve (ImplicitRadialUpdate)
//! and the RKL1 transverse operator applies two diffusion operators to the same energy
//! one after the other, each linearised about a state the other has already moved.  That
//! is a first-order-in-dt error, and it is not small where both directions carry flux:
//! it is what grew max|de| by ~1.5x per cycle in the He-star FeCZ box, at the first
//! active row above the bottom wall where the radial solve does its largest work.  With
//! every face of the 7-point stencil in one loop, linearised about one state, there is no
//! splitting error between the directions at all -- and because the x1 faces of an
//! isotropic box carry the same conductance as the x2/x3 ones, they raise the Gershgorin
//! row radius by only ~3/2 and the substage count by only ~20 %.
//!
//! WHAT IT IS NOT FOR: a radially stiff grid.  See the note on rad_sts_all in
//! conduction.hpp -- one global substage count is chosen from the stiffest row of the
//! whole mesh, so a stretched spherical grid would pay the radial stiffness ratio's
//! square root on every transverse face, where the column-local tridiagonal solve pays
//! nothing.  The constructor refuses curvilinear meshes outright.
//!
//! WHAT THE x1 FACES CARRY.  The interior x1 faces carry the frozen face conductance
//! cap_c1 = A_f K_f/dl_f that BuildAngularCoeffs forms with the SAME face_kcode (and so
//! the same RadFaceKappa, flux limiter, tau-blend weight and density gate) as the x2/x3
//! faces and as the explicit x1 face flux.  The two PHYSICAL x1 faces are NOT in the
//! stencil: they keep the explicit treatment they have under rad_implicit_x1 -- the
//! imposed internal flux rad_flux_inner through the bottom wall and the ghost-based
//! gradient at the top, both added to flx1(IEN) in AddIsotropicHeatFluxRadiative before
//! the RK update -- so the stencil's interior fluxes telescope and sum_i V_i y_i is the
//! imposed wall flux alone, to round-off.  The two-stream's explicit handover deposit is
//! untouched.

void Conduction::StsConductionUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                                     const Real beta_dt) {
  if (!rad_sts_all) return;
  RklConductionUpdate(u0, eos, beta_dt, true);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::RklConductionUpdate
//! \brief the RKL1 loop itself, over the transverse faces alone (`with_x1` false) or
//! over all three directions (`with_x1` true).  Every added term is guarded by that flag
//! and appended after the existing ones, so the transverse-only path is bitwise what it
//! was.

void Conduction::RklConductionUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                                     const Real beta_dt, const bool with_x1) {
  if (!(beta_dt > 0.0)) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool three_d = pmy_pack->pmesh->three_d;
  auto c1 = cap_c1;
  auto c2 = cap_c2;
  auto c3 = cap_c3;
  const bool sts1 = with_x1;
  auto st = tr_st;
  const int it_ = TRST, ia_ = TRSA;
  const Real tau = beta_dt;

  // ---- THE LINEARISATION POINT IS REFRESHED FROM THE CURRENT CONSERVED ENERGY.
  // BuildAngularCoeffs forms T*_i and alpha_i = 1/(rho_i c_v,i) from w0 in Hydro::Fluxes,
  // i.e. at the START of the stage -- before the RK update, the source terms and the
  // implicit RADIAL solve.  The increment y this operator computes is then added to the
  // energy those three have already moved.  So the operator relaxes the STALE horizontal
  // temperature differences and writes the result on top of a state in which the radial
  // solve has already changed them: the two split operators are one stage out of phase,
  // and a lagged diffusion operator is ANTI-diffusive with respect to the state it acts
  // on.  The error re-enters through K(T) on the next cycle, which is why the He-star
  // FeCZ box grew max|de| by ~1.5x per cycle from cycle 1, at the first active row above
  // the bottom wall where ImplicitRadialUpdate does its largest work.  The Gaussian test
  // cannot see it: there the state is uniform in x1 and the radial operator moves nothing.
  //
  // Re-evaluating T* and alpha here -- from u0, through the same EintFromCons and the
  // same EOS the radial solve uses -- puts the linearisation point back on the state the
  // increment is applied to.  The face conductances C_f stay frozen at w0: K is a far
  // weaker function of the state than the temperature difference it multiplies, and
  // keeping them fixed is what keeps the system linear and the stencil exactly
  // conservative.
  //
  // The refreshed values are needed on the x2/x3 ghost ring as well -- this operator runs
  // BEFORE the stage's ghost exchange, so u0 there is stale -- and the two cells sharing
  // a block face must use IDENTICAL numbers or the flux form stops cancelling.  They are
  // therefore built on the active cells and exchanged through the module's own one-
  // variable boundary object, exactly as the substages exchange Y.
  Real tshift = 0.0;    // max |T_new - T*|/T*, the size of the lag; reported below
  {
    auto eos_ = eos;
    const bool gen = eos.IsGeneral();
    const bool ismhd = (my_block.compare("mhd") == 0);
    auto &wtemp_ = ismhd ? pmy_pack->pmhd->wtemp : pmy_pack->phydro->wtemp;
    auto &phicc_ = ismhd ? pmy_pack->pmhd->phicc0 : pmy_pack->phydro->phicc0;
    const bool etg = ismhd ? pmy_pack->pmhd->use_etotgrav
                           : pmy_pack->phydro->use_etotgrav;
    const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
    auto &cosc_ = pmy_pack->pcoord->cos_cell;
    DvceArray5D<Real> bcc_("trs_bcc_dummy", 1, 1, 1, 1, 1);
    if (ismhd) bcc_ = pmy_pack->pmhd->bcc0;
    auto tnew = tr_ya;      // scratch: the registers are zeroed after this block
    auto anew = tr_yb;
    // start from the w0 values everywhere, so a ghost the exchange does not reach (a
    // PHYSICAL x2/x3 boundary, whose faces are closed anyway) still holds a sane number
    par_for("radtrseed", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      tnew(m,0,k,j,i) = st(m,it_,k,j,i);
      anew(m,0,k,j,i) = st(m,ia_,k,j,i);
    });
    par_for("radtrrefr", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = u0(m,IDN,k,j,i);
      const Real ei = EintFromCons(u0, m, k, j, i, cs_ ? cosc_(m,k,j) : 0.0, cs_, etg,
                                   etg ? phicc_(m,k,j,i) : 0.0,
                                   ismhd ? MagEnergyCC(bcc_,m,k,j,i) : 0.0);
      // a cell the EOS cannot invert is DROPPED from the system (alpha = 0 closes every
      // one of its faces, and the neighbour closes the same face by the same test), which
      // is what the stencil already does with a non-positive alpha
      if (!(ei > 0.0) || !(d > 0.0) || !isfinite(ei)) { anew(m,0,k,j,i) = 0.0; return; }
      const Real tt = eos_.Temperature(d, ei, gen ? wtemp_(m,k,j,i) : -1.0);
      if (!(tt > 0.0) || !isfinite(tt)) { anew(m,0,k,j,i) = 0.0; return; }
      const Real cv = eos_.SpecificHeatCv(d, ei, tt);
      if (!(cv > 0.0) || !isfinite(cv)) { anew(m,0,k,j,i) = 0.0; return; }
      tnew(m,0,k,j,i) = tt;
      anew(m,0,k,j,i) = 1.0/(d*cv);
    });
    // how far the state moved between the two operators, on the cells that carry flux
    {
      const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
      const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;
      Kokkos::parallel_reduce("radtrlag",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
      KOKKOS_LAMBDA(const int &idx, Real &sres) {
        const int m = idx/nkji_;
        const int k = (idx - m*nkji_)/nji_ + ks;
        const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
        const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
        const Real t0 = st(m,it_,k,j,i);
        if (!(t0 > 0.0) || !(anew(m,0,k,j,i) > 0.0)) return;
        const Real r = fabs(tnew(m,0,k,j,i) - t0)/t0;
        if (isfinite(r) && r > sres) sres = r;
      }, Kokkos::Max<Real>(tshift));
    }
    // the ghost ring, through the module's own exchange (one variable at a time)
    pbval_tr->InitRecv(1);
    pbval_tr->PackAndSendCC(tnew, tr_ycoar);
    while (pbval_tr->RecvAndUnpackCC(tnew, tr_ycoar) != TaskStatus::complete) {}
    while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
    while (pbval_tr->ClearSend() != TaskStatus::complete) {}
    pbval_tr->InitRecv(1);
    pbval_tr->PackAndSendCC(anew, tr_ycoar);
    while (pbval_tr->RecvAndUnpackCC(anew, tr_ycoar) != TaskStatus::complete) {}
    while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
    while (pbval_tr->ClearSend() != TaskStatus::complete) {}
    par_for("radtrcopy", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      st(m,it_,k,j,i) = tnew(m,0,k,j,i);
      st(m,ia_,k,j,i) = anew(m,0,k,j,i);
    });
  }

  // A face on a physical x2/x3 boundary is CLOSED: the cell on its far side is a ghost
  // whose increment nothing computes.  block and periodic neighbours are real and are
  // filled by the exchange below.  Both cells of a face test the same face index, so
  // they cannot disagree.
  auto open2 = [=] (const int m, const int jf) {
    if (jf == js) {
      const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::inner_x2);
      return (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
    }
    if (jf == je+1) {
      const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::outer_x2);
      return (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
    }
    return true;
  };
  // An x1 face is in the stencil only if it is INTERIOR.  The two physical x1 faces are
  // left to the explicit path (see the file comment), and the whole x1 extent is in one
  // MeshBlock by construction, so there is no x1 block face to exchange: the refreshed
  // T*/alpha of the x1 ghosts are never read, and their registers stay at the zero
  // Kokkos allocated them with.
  auto open1 = [=] (const int i) {
    return (i > is && i < ie+1);
  };
  auto open3 = [=] (const int m, const int kf) {
    if (kf == ks) {
      const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::inner_x3);
      return (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
    }
    if (kf == ke+1) {
      const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::outer_x3);
      return (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
    }
    return true;
  };

  // ---- the stiffness of the step: the largest Gershgorin row radius of tau dM/dy.  A
  // global maximum, because every rank must take the same number of substages.
  Real zmax = 0.0;
  {
    const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
    const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;
    Kokkos::parallel_reduce("radtrz",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
    KOKKOS_LAMBDA(const int &idx, Real &zres) {
      const int m = idx/nkji_;
      const int k = (idx - m*nkji_)/nji_ + ks;
      const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
      const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
      const Real ai = st(m,ia_,k,j,i);
      const Real cl2 = open2(m,j) ? c2(m,k,j,i) : 0.0;
      const Real cr2 = open2(m,j+1) ? c2(m,k,j+1,i) : 0.0;
      Real sumc = cl2 + cr2;
      Real sumca = cl2*st(m,ia_,k,j-1,i) + cr2*st(m,ia_,k,j+1,i);
      if (three_d) {
        const Real cl3 = open3(m,k) ? c3(m,k,j,i) : 0.0;
        const Real cr3 = open3(m,k+1) ? c3(m,k+1,j,i) : 0.0;
        sumc += cl3 + cr3;
        sumca += cl3*st(m,ia_,k-1,j,i) + cr3*st(m,ia_,k+1,j,i);
      }
      // rad_sts_all: the x1 faces are part of the same row, so they must be part of the
      // same bound -- the substage count has to cover the stiffest direction
      if (sts1) {
        const Real cl1 = open1(i) ? c1(m,k,j,i) : 0.0;
        const Real cr1 = open1(i+1) ? c1(m,k,j,i+1) : 0.0;
        sumc += cl1 + cr1;
        sumca += cl1*st(m,ia_,k,j,i-1) + cr1*st(m,ia_,k,j,i+1);
      }
      // V_i = 1 on a Cartesian mesh: the face coefficients already carry the 1/dx that
      // the flux-divergence form of the RK update applies (see BuildAngularCoeffs)
      const Real zi = tau*(ai*sumc + sumca);
      if (isfinite(zi) && zi > zres) zres = zi;
    }, Kokkos::Max<Real>(zmax));
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &zmax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  if (!(zmax > 0.0)) return;    // no face carries any flux this stage: nothing to do

  // ---- the substage count.  R is the super-step in units of the explicit limit
  // 2/lambda_max, with 10 % of round-off margin; RKL1 with s stages covers (s^2+s)/2.
  // rad_sts_margin: the round-off margin on the RKL1 stability limit, 0.10 by default
  // (0.5 x 1.1 = 0.55, the value this loop has always used).  z_i is a Gershgorin radius
  // and therefore already an over-estimate of |lambda|, so 0.02 is safe and buys a few
  // per cent of the substage count.
  const Real rstiff = 0.5*(1.0 + rad_sts_margin)*zmax;
  int nsub = static_cast<int>(std::ceil(0.5*(std::sqrt(1.0 + 8.0*rstiff) - 1.0)));
  if (nsub < 1) nsub = 1;
  bool clamped = false;
  if (nsub > rad_ang_maxit) {
    nsub = rad_ang_maxit;
    clamped = true;
  }
  const Real w1 = 2.0/(static_cast<Real>(nsub)*static_cast<Real>(nsub) + nsub);

  // ---- the RKL1 loop on the increment.  Y_0 = 0 everywhere INCLUDING the ghosts, so
  // the first substage needs no exchange; every later one exchanges Y_{j-1} through the
  // module's own one-variable boundary object before the stencil reads it.
  auto ycur = tr_ya;   // Y_{j-1}
  auto yold = tr_yb;   // Y_{j-2}
  auto ynew = tr_yc;   // Y_j
  Kokkos::deep_copy(ycur, 0.0);
  Kokkos::deep_copy(yold, 0.0);
  for (int js_ = 1; js_ <= nsub; ++js_) {
    if (js_ > 1) {
      // post the receives, send, and spin on the unpack: RecvAndUnpackCC is the only
      // one of these that can legitimately come back incomplete (the MPI traffic of
      // this substage is all there is to overlap it with)
      pbval_tr->InitRecv(1);
      pbval_tr->PackAndSendCC(ycur, tr_ycoar);
      while (pbval_tr->RecvAndUnpackCC(ycur, tr_ycoar) != TaskStatus::complete) {}
      while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
      while (pbval_tr->ClearSend() != TaskStatus::complete) {}
    }
    const Real muj = (js_ == 1) ? 1.0 : (2.0*js_ - 1.0)/js_;
    const Real nuj = (js_ == 1) ? 0.0 : (1.0 - js_)/js_;
    const Real mut = muj*w1*tau;
    auto yc_ = ycur;
    auto yo_ = yold;
    auto yn_ = ynew;
    par_for("radtrsub", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real ai = st(m,ia_,k,j,i);
      const Real thc = st(m,it_,k,j,i) + ai*yc_(m,0,k,j,i);
      // the FLUX through each face, written so that the neighbour forms the identical
      // expression from the identical operands
      const Real alm2 = st(m,ia_,k,j-1,i), alp2 = st(m,ia_,k,j+1,i);
      Real fl = (open2(m,j) && ai > 0.0 && alm2 > 0.0)
                ? c2(m,k,j,i)*(thc - (st(m,it_,k,j-1,i) + alm2*yc_(m,0,k,j-1,i))) : 0.0;
      Real fr = (open2(m,j+1) && ai > 0.0 && alp2 > 0.0)
                ? c2(m,k,j+1,i)*((st(m,it_,k,j+1,i) + alp2*yc_(m,0,k,j+1,i)) - thc) : 0.0;
      Real mi = fr - fl;
      if (three_d) {
        const Real alm3 = st(m,ia_,k-1,j,i), alp3 = st(m,ia_,k+1,j,i);
        const Real gl = (open3(m,k) && ai > 0.0 && alm3 > 0.0)
                ? c3(m,k,j,i)*(thc - (st(m,it_,k-1,j,i) + alm3*yc_(m,0,k-1,j,i))) : 0.0;
        const Real gr = (open3(m,k+1) && ai > 0.0 && alp3 > 0.0)
                ? c3(m,k+1,j,i)*((st(m,it_,k+1,j,i) + alp3*yc_(m,0,k+1,j,i)) - thc) : 0.0;
        mi += gr - gl;
      }
      if (sts1) {
        const Real alm1 = st(m,ia_,k,j,i-1), alp1 = st(m,ia_,k,j,i+1);
        const Real hl = (open1(i) && ai > 0.0 && alm1 > 0.0)
                ? c1(m,k,j,i)*(thc - (st(m,it_,k,j,i-1) + alm1*yc_(m,0,k,j,i-1))) : 0.0;
        const Real hr = (open1(i+1) && ai > 0.0 && alp1 > 0.0)
                ? c1(m,k,j,i+1)*((st(m,it_,k,j,i+1) + alp1*yc_(m,0,k,j,i+1)) - thc) : 0.0;
        mi += hr - hl;
      }
      if (!isfinite(mi)) mi = 0.0;
      yn_(m,0,k,j,i) = muj*yc_(m,0,k,j,i) + nuj*yo_(m,0,k,j,i) + mut*mi;
    });
    auto tmp = yold;
    yold = ycur;
    ycur = ynew;
    ynew = tmp;
  }

  // ---- write the increment into the energy, and measure what it did to the total
  Real esum = 0.0, eabs = 0.0, emax = 0.0;
  {
    const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
    const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;
    auto yfin = ycur;
    Kokkos::parallel_reduce("radtrend",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
    KOKKOS_LAMBDA(const int &idx, Real &ssum, Real &sabs, Real &smax) {
      const int m = idx/nkji_;
      const int k = (idx - m*nkji_)/nji_ + ks;
      const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
      const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
      Real y = yfin(m,0,k,j,i);
      if (!isfinite(y)) y = 0.0;
      u0(m,IEN,k,j,i) += y;
      const Real dv = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
      ssum += dv*y;
      sabs += dv*fabs(y);
      smax = fmax(smax, fabs(y));
    }, Kokkos::Sum<Real>(esum), Kokkos::Sum<Real>(eabs), Kokkos::Max<Real>(emax));
  }
#if MPI_PARALLEL_ENABLED
  {
    Real buf[2] = {esum, eabs};
    MPI_Allreduce(MPI_IN_PLACE, buf, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    esum = buf[0];
    eabs = buf[1];
    MPI_Allreduce(MPI_IN_PLACE, &emax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
  }
#endif
  const Real viol = (eabs > 0.0) ? fabs(esum)/eabs : 0.0;
  sts_nsub_tot += nsub;
  ++sts_ncall;

  // the report: the substage count is the cost of the operator and the residual is the
  // one thing that can silently go wrong, so both are printed.  A clamped substage count
  // means the step was NOT covered and is always worth a line.
  if (ang_lines < 400 && global_variable::my_rank == 0) {
    const bool tell = rad_ang_verbose && (ang_lines < 20 || viol > 1.0e-10);
    if (clamped || tell) {
      ++ang_lines;
      std::cout << (sts1 ? "### rad_sts_all cycle " : "### rad_implicit_ang cycle ")
                << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time
                << ": max z_i = " << zmax << ", substages = " << nsub
                << ", substages_total = " << sts_nsub_tot
                << " in " << sts_ncall << " calls"
                << (clamped ? " (CLAMPED at rad_ang_maxit -- the step is NOT covered)"
                            : "")
                << ", max |de| = " << emax
                << ", max |dT*|/T* = " << tshift
                << ", |sum V de|/sum V|de| = " << viol << std::endl;
    }
  }
  return;
}
