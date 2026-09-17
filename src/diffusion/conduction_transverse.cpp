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
#include <cstddef>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/nghbr_index.hpp"
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
  // ---- CUBED-SPHERE GEOMETRY.  Everything here is inert off a cubed sphere: curv gates
  // the CELL VOLUME (the Cartesian face coefficients already carry the 1/dx that the
  // flux-divergence form applies, so V_i = 1 there and the division is by an exact 1.0),
  // and csx gates the METRIC CROSS TERM.  See docs/dev/cs_implicit_transverse.md.
  // Spherical polar is refused by the constructor, so use_cubed_sphere IS "curvilinear"
  // here.
  const bool curv = pmy_pack->pmesh->use_cubed_sphere;
  const bool csx = curv && rad_cs_exact && three_d;
  auto vol_ = pmy_pack->pcoord->volume;
  auto dx2c_ = pmy_pack->pcoord->dx2;
  auto dx3c_ = pmy_pack->pcoord->dx3;
  auto g2 = csx ? cap_g2 : DvceArray4D<Real>("radtrg2dummy", 1, 1, 1, 1);
  auto g3 = csx ? cap_g3 : DvceArray4D<Real>("radtrg3dummy", 1, 1, 1, 1);
  // ---- THE GHOST-SKIN BOOKKEEPING (rad_tr_halo_every).  N = 1 is the old loop: the
  // increment is exchanged before every substage and computed on the ACTIVE cells only.
  //
  // With N > 1 the exchange happens only every N substages and the intermediate
  // substages compute on a GHOST SKIN as well, so that the ghost values they need are
  // there without a message.  What the skin costs is arithmetic on at most ng-1 extra
  // ghost layers; what it saves is (N-1)/N of the halo traffic.
  //
  // THE VALID REGION.  Call d(y) the number of x2/x3 ghost layers on which y is correct
  // (a FULL frame -- the x2x3 edge/corner ghosts included, which is why this switch
  // needs the diagonal buffers and rad_tr_halo_faces_only is refused with it).  An
  // exchange sets d = ng.  The recurrence Y_j = mu Y_{j-1} + nu Y_{j-2} + mu~ tau
  // M(Y_{j-1}) reads Y_{j-1} through the 5-POINT CROSS (one layer of reach) and Y_{j-2}
  // POINTWISE, so
  //     d(Y_j) = min( d(Y_{j-1}) - 1, d(Y_{j-2}) ),
  // which is tracked exactly below in dcur/dold and must never go negative on the active
  // cells (d >= 0).  Feeding only Y_{j-1} to a refill gives the steady state
  // d = (ng-1, ng-2, ...) and survives N = 2 but NOT N >= 3 -- Y_{j-2} then runs out
  // first -- so a refill at N >= 3 exchanges the second register too, and the traffic is
  // 2 arrays per N substages instead of N.  N = 2 is therefore the cheapest cadence
  // (1 array per 2 substages) and N = 3 the deepest one nghost = 3 allows.
  const int ng_ = indcs.ng;
  int hevery = rad_tr_halo_every;
  if (hevery < 1) hevery = 1;
  if (hevery > ng_) hevery = ng_;
  const bool skin_ = (hevery > 1);
  const bool xch2_ = (hevery >= 3);   // a refill also exchanges Y_{j-2}
  auto c1 = cap_c1;
  auto c2 = cap_c2;
  auto c3 = cap_c3;
  const bool sts1 = with_x1;
  // rad_sts_split: cap_c1/2/3 hold C_sts, the part of each face conductance the explicit
  // face fluxes of the stage could NOT carry, and a MeshBlock on which nothing is left
  // is skipped by the stencil below.  Its increment is identically zero, so the skip
  // changes no result -- it only stops the kernel from reading the block's 7-point
  // neighbourhood.  Nothing here is aware of the split otherwise: the loop integrates
  // whatever conductances it is handed.
  const bool blkon = rad_sts_split && sts_blk_used;
  auto blk = blkon ? sts_blk : DvceArray1D<int>("rklblkdummy", 1);
  auto st = tr_st;
  const int it_ = TRST, ia_ = TRSA;
  const Real tau = beta_dt;
  // PER-PLANE SUBSTAGE COUNTS (rad_sts_perplane).  Only for the 5-point transverse
  // stencil: with_x1 puts x1 faces in the row, which couple the planes, and the substage
  // count must then be global.  See conduction.hpp.
  // rad_ang_solver = adi replaces the RKL1 loop below with one alternating-direction
  // implicit step; it is transverse-only, and it needs the per-plane stiffness radii (for
  // the plane skip), so it turns the per-plane path on for itself.
  const bool adi = rad_ang_adi && !sts1;
  const bool perpl = (rad_sts_perplane || adi) && !sts1;
  const int nplane = ie - is + 1;
  const int isv = is;
  auto spl = perpl ? tr_spl : DvceArray1D<int>("rklspldummy", 1);
  auto w1pl = perpl ? tr_w1pl : DvceArray1D<Real>("rklw1dummy", 1);

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
  // cannot see it: there the state is uniform in x1 and the radial operator moves
  // nothing.
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
    // the SKIN needs T* and alpha over the whole ng-deep frame, not one layer: the
    // exchange below already fills tnew/anew that deep (a same-level face buffer is ng
    // cells deep), so this only has to copy them out that far.  At hevery = 1 the extra
    // layers are never read and sdg is 1, i.e. exactly the old ranges.
    const int sdg = skin_ ? ng_ : 1;
    par_for("radtrseed", DevExeSpace(), 0, nmb1, ks-sdg, ke+sdg, js-sdg, je+sdg, is, ie,
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
    // ONLY the ACTIVE x1 range is exchanged.  Every read of the T*/alpha ghost ring is
    // st(...,k+-1,j+-1,i) with i an ACTIVE plane (the stencil, the Gershgorin reduction
    // and the copy-back all run over is..ie), so the x1-ghost columns of the transverse
    // halo are dead weight -- radtrseed/radtrcopy do not even write them.
    // rad_tr_window = false exchanges the FULL x1 range instead of the window
    const int sl_ = rad_tr_window ? is : -1, su_ = rad_tr_window ? ie : -1;
    pbval_tr->InitRecv(1);
    pbval_tr->PackAndSendCC(tnew, tr_ycoar, sl_, su_);
    while (pbval_tr->RecvAndUnpackCC(tnew, tr_ycoar, sl_, su_) != TaskStatus::complete) {}
    // RecvAndUnpackCC launches the unpack kernel ASYNCHRONOUSLY and returns.  The next
    // InitRecv re-posts MPI_Irecv into the SAME device recv buffer, so with GPU-aware
    // MPI a neighbour that is ahead can have the NIC overwrite that buffer while the
    // unpack kernel is still reading it: fence before ClearRecv/ClearSend and InitRecv.
    Kokkos::fence();
    while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
    while (pbval_tr->ClearSend() != TaskStatus::complete) {}
    pbval_tr->InitRecv(1);
    pbval_tr->PackAndSendCC(anew, tr_ycoar, sl_, su_);
    while (pbval_tr->RecvAndUnpackCC(anew, tr_ycoar, sl_, su_) != TaskStatus::complete) {}
    Kokkos::fence();   // see the comment above: the recv buffer is re-posted below
    while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
    while (pbval_tr->ClearSend() != TaskStatus::complete) {}
    par_for("radtrcopy", DevExeSpace(), 0, nmb1, ks-sdg, ke+sdg, js-sdg, je+sdg, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      st(m,it_,k,j,i) = tnew(m,0,k,j,i);
      st(m,ia_,k,j,i) = anew(m,0,k,j,i);
    });
  }

  // A face on a physical x2/x3 boundary is CLOSED: the cell on its far side is a ghost
  // whose increment nothing computes.  block and periodic neighbours are real and are
  // filled by the exchange below.  Both cells of a face test the same face index, so
  // they cannot disagree.
  // A PANEL SEAM is a real, flux-carrying face: the cell on its far side is a live cell
  // of another panel and the module's own exchange fills it (MeshBoundaryValuesCC does
  // the index permutation and the along-seam resample for a cell-centred scalar).  It is
  // OPEN but not LINKED: no tridiagonal line can cross it, because the neighbouring
  // panel's x2 may be this panel's x3 and the two charts' cells do not coincide along
  // the seam.  The ADI therefore leaves seam faces out of its sweeps and applies them in
  // its own pair-implicit sub-step; RKL1 has no line structure and treats them as
  // ordinary open faces.  On a Cartesian mesh `panel` never occurs, so open == linked
  // and every truth table below is what it was.
  auto bflag2 = [=] (const int m, const int jf) {
    if (jf == js) return mb_bcs.d_view(m,BoundaryFace::inner_x2);
    if (jf == je+1) return mb_bcs.d_view(m,BoundaryFace::outer_x2);
    return BoundaryFlag::block;
  };
  auto open2 = [=] (const int m, const int jf) {
    const BoundaryFlag f = bflag2(m,jf);
    return (f == BoundaryFlag::block || f == BoundaryFlag::periodic ||
            f == BoundaryFlag::panel);
  };
  auto lnk2 = [=] (const int m, const int jf) {
    const BoundaryFlag f = bflag2(m,jf);
    return (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
  };
  auto seam2 = [=] (const int m, const int jf) {
    return (bflag2(m,jf) == BoundaryFlag::panel);
  };
  // An x1 face is in the stencil only if it is INTERIOR.  The two physical x1 faces are
  // left to the explicit path (see the file comment), and the whole x1 extent is in one
  // MeshBlock by construction, so there is no x1 block face to exchange: the refreshed
  // T*/alpha of the x1 ghosts are never read, and their registers stay at the zero
  // Kokkos allocated them with.
  auto open1 = [=] (const int i) {
    return (i > is && i < ie+1);
  };
  auto bflag3 = [=] (const int m, const int kf) {
    if (kf == ks) return mb_bcs.d_view(m,BoundaryFace::inner_x3);
    if (kf == ke+1) return mb_bcs.d_view(m,BoundaryFace::outer_x3);
    return BoundaryFlag::block;
  };
  auto open3 = [=] (const int m, const int kf) {
    const BoundaryFlag f = bflag3(m,kf);
    return (f == BoundaryFlag::block || f == BoundaryFlag::periodic ||
            f == BoundaryFlag::panel);
  };
  auto lnk3 = [=] (const int m, const int kf) {
    const BoundaryFlag f = bflag3(m,kf);
    return (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
  };
  auto seam3 = [=] (const int m, const int kf) {
    return (bflag3(m,kf) == BoundaryFlag::panel);
  };

  // ---- THE METRIC CROSS TERM.  On the gnomonic panel the xi and eta lines meet at an
  // angle alpha, so the total flux through an x2 face is
  //     Phi_f = -C_f (T_j - T_i) + G_f ge_f,   G_f = C_f dl_f cos(alpha),
  // with ge_f the face-tangential (eta) derivative -- the SAME two-cell average the
  // explicit cubed-sphere face flux of AddIsotropicHeatFluxRadiative forms, so the two
  // operators cannot disagree about the geometry.  The second term couples j to k and is
  // therefore not part of the 5-point/tridiagonal structure: it is carried explicitly.
  // Th = T* + alpha y when a register is supplied (usey), T* alone otherwise.
  auto thof = [=] (const DvceArray5D<Real> &yv, const bool usey,
                   const int m, const int k, const int j, const int i) {
    return usey ? (st(m,it_,k,j,i) + st(m,ia_,k,j,i)*yv(m,0,k,j,i))
                : st(m,it_,k,j,i);
  };
  auto cross2 = [=] (const DvceArray5D<Real> &yv, const bool usey,
                     const int m, const int k, const int j, const int i) {
    const Real gg = g2(m,k,j,i);
    const Real dzm = 0.5*dx3c_(m,k-1,j-1,i) + dx3c_(m,k,j-1,i) + 0.5*dx3c_(m,k+1,j-1,i);
    const Real dzp = 0.5*dx3c_(m,k-1,j,i) + dx3c_(m,k,j,i) + 0.5*dx3c_(m,k+1,j,i);
    const Real ge = 0.5*((thof(yv,usey,m,k+1,j-1,i) - thof(yv,usey,m,k-1,j-1,i))/dzm
                       + (thof(yv,usey,m,k+1,j,i) - thof(yv,usey,m,k-1,j,i))/dzp);
    const Real q = gg*ge;
    return isfinite(q) ? q : static_cast<Real>(0.0);
  };
  auto cross3 = [=] (const DvceArray5D<Real> &yv, const bool usey,
                     const int m, const int k, const int j, const int i) {
    const Real gg = g3(m,k,j,i);
    const Real dym = 0.5*dx2c_(m,k-1,j-1,i) + dx2c_(m,k-1,j,i) + 0.5*dx2c_(m,k-1,j+1,i);
    const Real dyp = 0.5*dx2c_(m,k,j-1,i) + dx2c_(m,k,j,i) + 0.5*dx2c_(m,k,j+1,i);
    const Real gx = 0.5*((thof(yv,usey,m,k-1,j+1,i) - thof(yv,usey,m,k-1,j-1,i))/dym
                       + (thof(yv,usey,m,k,j+1,i) - thof(yv,usey,m,k,j-1,i))/dyp);
    const Real q = gg*gx;
    return isfinite(q) ? q : static_cast<Real>(0.0);
  };
  // the cross term's contribution to the Gershgorin ROW RADIUS of the face: the sum of
  // |entry| over the four cells its stencil reads, each weighted by that cell's alpha
  auto crad2 = [=] (const int m, const int k, const int j, const int i) {
    const Real gg = fabs(g2(m,k,j,i));
    if (!(gg > 0.0)) return static_cast<Real>(0.0);
    const Real dzm = 0.5*dx3c_(m,k-1,j-1,i) + dx3c_(m,k,j-1,i) + 0.5*dx3c_(m,k+1,j-1,i);
    const Real dzp = 0.5*dx3c_(m,k-1,j,i) + dx3c_(m,k,j,i) + 0.5*dx3c_(m,k+1,j,i);
    const Real r = 0.5*gg*((st(m,ia_,k+1,j-1,i) + st(m,ia_,k-1,j-1,i))/dzm
                         + (st(m,ia_,k+1,j,i) + st(m,ia_,k-1,j,i))/dzp);
    return isfinite(r) ? r : static_cast<Real>(0.0);
  };
  auto crad3 = [=] (const int m, const int k, const int j, const int i) {
    const Real gg = fabs(g3(m,k,j,i));
    if (!(gg > 0.0)) return static_cast<Real>(0.0);
    const Real dym = 0.5*dx2c_(m,k-1,j-1,i) + dx2c_(m,k-1,j,i) + 0.5*dx2c_(m,k-1,j+1,i);
    const Real dyp = 0.5*dx2c_(m,k,j-1,i) + dx2c_(m,k,j,i) + 0.5*dx2c_(m,k,j+1,i);
    const Real r = 0.5*gg*((st(m,ia_,k-1,j+1,i) + st(m,ia_,k-1,j-1,i))/dym
                         + (st(m,ia_,k,j+1,i) + st(m,ia_,k,j-1,i))/dyp);
    return isfinite(r) ? r : static_cast<Real>(0.0);
  };
  // V_i.  Exactly 1.0 off a cubed sphere, so the divisions below are bitwise no-ops
  // there (see the note on the Cartesian face coefficients above).
  auto vcell = [=] (const int m, const int k, const int j, const int i) {
    return curv ? vol_(m,k,j,i) : static_cast<Real>(1.0);
  };

  // ---- the stiffness of the step: the largest Gershgorin row radius of tau dM/dy.  A
  // global maximum, because every rank must take the same number of substages -- or, with
  // rad_sts_perplane, one global maximum PER x1 PLANE, because the planes do not talk to
  // each other and each may take its own.
  Real zmax = 0.0;
  if (perpl) {
    auto zpl = tr_zpl;
    Kokkos::deep_copy(zpl, 0.0);
    par_for("radtrzpl", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
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
      // the cubed-sphere cross term adds four entries per face to the row; V_i = 1 off
      // a cubed sphere, so the Cartesian radius is bitwise what it was
      // the CARTESIAN expression is kept VERBATIM in its own branch: adding a
      // zero cross radius and dividing by an exact 1.0 are algebraic no-ops, but they
      // change how the compiler contracts the product into an FMA, and this radius
      // decides the substage count and the active-plane bracket
      Real zi;
      if (curv) {
        Real xrad = 0.0;
        if (csx) {
          xrad = crad2(m,k,j,i) + crad2(m,k,j+1,i)
               + crad3(m,k,j,i) + crad3(m,k+1,j,i);
        }
        zi = tau*(ai*sumc + sumca + xrad)/vol_(m,k,j,i);
      } else {
        zi = tau*(ai*sumc + sumca);
      }
      if (isfinite(zi) && zi > 0.0) Kokkos::atomic_max(&zpl(i - isv), zi);
    });
    Kokkos::deep_copy(tr_zpl_h, zpl);
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, tr_zpl_h.data(), nplane, MPI_ATHENA_REAL, MPI_MAX,
                  MPI_COMM_WORLD);
#endif
    for (int p=0; p<nplane; ++p) {
      if (tr_zpl_h(p) > zmax) zmax = tr_zpl_h(p);
    }
  } else {
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
      // the flux-divergence form of the RK update applies (see BuildAngularCoeffs).  On
      // a cubed sphere they carry the face AREA instead and the volume divides here, and
      // the metric cross term adds four entries per face to the row.
      Real zi;
      if (curv) {
        Real xrad = 0.0;
        if (csx) {
          xrad = crad2(m,k,j,i) + crad2(m,k,j+1,i)
               + crad3(m,k,j,i) + crad3(m,k+1,j,i);
        }
        zi = tau*(ai*sumc + sumca + xrad)/vol_(m,k,j,i);
      } else {
        zi = tau*(ai*sumc + sumca);
      }
      if (isfinite(zi) && zi > zres) zres = zi;
    }, Kokkos::Max<Real>(zmax));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &zmax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  }
  if (!(zmax > 0.0)) return;    // no face carries any flux this stage: nothing to do

  // ====================================================================================
  // rad_ang_solver = adi: ONE ALTERNATING-DIRECTION IMPLICIT STEP over the whole stage,
  // in place of the RKL1 loop below.
  //
  // THE SYSTEM is the one the file comment defines and the RKL1 substage discretises.
  // M(y) is affine, M(y) = b + A y with
  //     b_i     = (1/V_i) sum_f s_f C_f (T*_j - T*_i)                       = M(0)_i,
  //     (A y)_i = (1/V_i) sum_f s_f C_f (alpha_j y_j - alpha_i y_i),
  // over the SAME open faces, with the SAME alpha > 0 gates and the SAME flux-form
  // expression, so ADI and RKL1 are two time discretisations of ONE operator.  Split
  // A = A2 + A3 by direction and integrate y(0) = 0 to tau = beta_dt by the DOUGLAS
  // scheme, theta-weighted:
  //     (I - theta tau A2) Y1 = tau b,      (I - theta tau A3) y = Y1,
  // whose product form is (I - theta tau A2)(I - theta tau A3) y = tau b.  theta = 1 is
  // backward Euler plus the O(tau^2 A2 A3) splitting term -- unconditionally stable and
  // damping in the stiff limit, which is what a stiff transverse operator needs -- and
  // theta = 0.5 is the Peaceman-Rachford variant, second order in tau.  One step is taken
  // however stiff the row is: there is no substage count.
  //
  // CONSERVATION.  sum_i V_i (A v)_i = 0 identically for any v (the flux form telescopes,
  // and a face on a physical boundary is closed), and sum_i V_i b_i = 0 for the same
  // reason.  Summing the first sweep gives sum V Y1 = tau sum V b + theta tau sum V A2 Y1
  // = 0, and the second gives sum V y = sum V Y1 = 0.  What is left is the round-off of
  // the line solves, which the report measures exactly as before.
  //
  // THE LINE SOLVES CROSS MeshBlocks AND RANKS, and are solved by PARTITION.  See the
  // rad_ang_solver note in conduction.hpp for the pattern; here is the algebra.  Order
  // the cells of one line inside one block s = 1..n.  The local rows are tridiagonal
  // except that row 1 couples to the left neighbour's row n through a_1 and row n to the
  // right neighbour's row 1 through c_n, so with x_L, x_R those two outside unknowns,
  //     T x = d - a_1 x_L e_1 - c_n x_R e_n,   x = u - x_L v' - x_R w',
  //     u = T^-1 d,   v' = a_1 T^-1 e_1,   w' = c_n T^-1 e_n,
  // i.e. ONE Thomas factorisation with THREE right-hand sides.  Its first and last rows,
  //     p_b + A_b q_{b-1} + B_b p_{b+1} = U_b,   q_b + C_b q_{b-1} + D_b p_{b+1} = Q_b,
  //     (A,B,U,C,D,Q)_b = (v'_1, w'_1, u_1, v'_n, w'_n, u_n),
  // are the REDUCED SYSTEM: 2 unknowns per block, 2*nblk in all, cyclic because the line
  // is periodic.  The six numbers are gathered around the ring of blocks, the reduced
  // system is assembled in ABSOLUTE block order (so every block of the ring does bitwise
  // identical arithmetic and they cannot disagree about an interface value) and solved
  // redundantly, and each block back-substitutes x = u - x_L v' - x_R w' locally.  A
  // closed face gives a_1 = 0 or c_n = 0 and the chain simply decouples there, so a
  // single-block periodic direction (nblk = 1, the ring is the block itself) and a
  // physical boundary are the same code.
  //
  // PLANE SKIP, as in the RKL1 path: a plane whose Gershgorin radius z_i <= 1 is not
  // stiff, RKL1 gives it s_i = 1 and hence y = tau M(0) = tau b, and that is what it is
  // given here; the sweeps run over the contiguous bracket of ACTIVE planes only.
  if (adi) {
    auto yy = tr_ya;      // the working register: right-hand side in, answer out
    auto uu = tr_yb;      // T^-1 d
    auto vv = tr_yc;      // the left spike, scaled by a_1
    auto ww = tr_aw;      // the right spike, scaled by c_n
    auto cpv = tr_acp;    // the Thomas upper-diagonal scratch
    auto yf = tr_yf;      // lod2 only: the full-step LOD answer, kept for Richardson
    auto rd_ = tr_ared;
    auto act = tr_aact;
    auto ab0 = tr_ab0;
    const bool adilod = rad_adi_lod;
    const int scm = rad_adi_scm;
    // a cubed sphere ALWAYS has panel seams (every panel edge is one), and a Cartesian
    // mesh never has any.  A collective flag, because the seam sub-step exchanges.
    const bool anyseam = curv;
    // A cubed-sphere line is an OPEN CHAIN (it ends at the two panel seams); a Cartesian
    // multi-block line is a CLOSED RING (it is required to be periodic).  Only the
    // GATHER differs -- see (2) below; the reduced system itself is the same cyclic
    // assembly, which degenerates to the chain because a seam face gives a_1 = c_n = 0.
    const bool openchain = curv;

    // ---- the active plane bracket, from the per-plane Gershgorin radii computed above
    int alo = nplane, ahi = -1;
    for (int p=0; p<nplane; ++p) {
      const int a = (tr_zpl_h(p) > 1.0) ? 1 : 0;
      tr_aact_h(p) = a;
      if (a) {
        if (p < alo) alo = p;
        ahi = p;
      }
    }
    Kokkos::deep_copy(act, tr_aact_h);
    const bool anyact = (ahi >= alo);
    const int ilo_ = anyact ? (is + alo) : is;
    const int ihi_ = anyact ? (is + ahi) : is;
    const int nplact = ihi_ - ilo_ + 1;
    const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
    const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;

    // ---- the CHAIN of MeshBlocks each line crosses.  One entry per (direction, local
    // block): this block's index along the direction WITHIN ITS OWN PANEL, and the rank
    // and local id of its two face neighbours, from the neighbour table (no SMR, so one
    // neighbour per face and the same-level slot is the only one filled).
    //
    // A neighbour counts here only if the face is LINKED, i.e. `block` or `periodic` --
    // the same truth table lnk2/lnk3 use, so the chain of the gather and the chain of
    // the tridiagonal system cannot disagree.  A PANEL SEAM is therefore an END of the
    // chain: the line stops there and the seam face is applied by the pair-implicit
    // sub-step below.  Each panel tree is its own root grid, so ll.lx2/ll.lx3 already
    // count blocks WITHIN the panel and adi_nb2/adi_nb3 are blocks PER PANEL.  On a
    // Cartesian mesh nothing here changes: every multi-block direction is periodic, so
    // every face of it is linked and the chain is the old closed ring.
    const int nmbl = nmb1 + 1;
    std::vector<int> ab0h(2*nmbl, 0), lrk(2*nmbl, -1), rrk(2*nmbl, -1);
    std::vector<int> llid(2*nmbl, -1), rlid(2*nmbl, -1);
    {
      auto &nghbr = pmy_pack->pmb->nghbr;
      auto &gidh = pmy_pack->pmb->mb_gid;
      const int nl2 = NeighborIndex(0,-1,0,0,0), nr2 = NeighborIndex(0,1,0,0,0);
      const int nl3 = NeighborIndex(0,0,-1,0,0), nr3 = NeighborIndex(0,0,1,0,0);
      const int bface[4] = {BoundaryFace::inner_x2, BoundaryFace::outer_x2,
                            BoundaryFace::inner_x3, BoundaryFace::outer_x3};
      for (int m=0; m<nmbl; ++m) {
        const LogicalLocation &ll = pmy_pack->pmesh->lloc_eachmb[gidh.h_view(m)];
        ab0h[m] = static_cast<int>(ll.lx2);
        ab0h[nmbl + m] = static_cast<int>(ll.lx3);
        const int slot[4] = {nl2, nr2, nl3, nr3};
        for (int q=0; q<4; ++q) {
          if (q >= 2 && !three_d) continue;
          const BoundaryFlag bf = mb_bcs.h_view(m,bface[q]);
          if (!(bf == BoundaryFlag::block || bf == BoundaryFlag::periodic)) continue;
          const int gg = nghbr.h_view(m,slot[q]).gid;
          const int rr = nghbr.h_view(m,slot[q]).rank;
          if (gg < 0) continue;
          const int lid = gg - pmy_pack->pmesh->gids_eachrank[rr];
          const int o = (q/2)*nmbl + m;
          if ((q%2) == 0) {
            lrk[o] = rr;
            llid[o] = lid;
          } else {
            rrk[o] = rr;
            rlid[o] = lid;
          }
        }
      }
    }

    // ---- THE SCHEMES.  A SUB-STEP is a pair of directional sweeps over an interval
    // dts, starting from a given y_in.  Written out, with b = b2 + b3 the direction-split
    // source and A = A2 + A3 the direction-split matrix:
    //
    // DOUGLAS (rad_adi_scheme = douglas) puts the WHOLE source in the first sweep,
    //     (I - th dts A2) Y1 = y_in + dts b,    (I - th dts A3) y = Y1,
    // product form (I - th dts A2)(I - th dts A3) y = y_in + dts b.  It is the textbook
    // scheme and second order at th = 0.5, but it is NOT STIFFLY ACCURATE: for one
    // Fourier mode with dts lambda2 = -a, dts lambda3 = -b it returns
    //     y = -(dT/alpha) (a + b)/((1+a)(1+b))
    // against the exact -(dT/alpha)(1 - e^-(a+b)).  A mode stiff in BOTH directions --
    // the horizontal checkerboard, which is what this operator exists to damp -- has
    // a = b >> 1 and gets ~2/a of the relaxation it needs: at the B-star box's z ~ 1500
    // that is 0.3 %, i.e. the operator does essentially NOTHING to the checkerboard.
    // MEASURED: the B-star 1-rank gate arm dies at cycle 3411 with the transverse kinetic
    // energy 2e4 times the RKL1 arm's.  Kept only as a switch, for the record.
    //
    // LOD (the default), the sequential (Lie / locally-one-dimensional) splitting: one
    // BACKWARD-EULER step of each sub-problem in turn over the whole interval, each
    // carrying ITS OWN part of the source,
    //     (I - dts A2) Y1 = y_in + dts b2,     (I - dts A3) y = Y1 + dts b3.
    // Same two line solves, same cost, same exact conservation (each sweep's flux form
    // telescopes on its own), first order in dts -- but the stiff limit is
    //     y = -(dT/alpha) [a/(1+a) + b]/(1+b)  ->  -(dT/alpha),
    // i.e. a mode stiff in either or both directions is relaxed essentially completely,
    // which is what backward Euler on the unsplit operator would do.  L-stable, and that
    // is why it, not Douglas, is what a stiff transverse operator can use.
    //
    // WHAT LOD STILL COSTS: it is only FIRST order, and its stiff damping factor is the
    // backward-Euler 1/(1+z), which for a moderately stiff mode is far above the exact
    // e^-z -- the mode is UNDER-DAMPED.  MEASURED on the Gaussian test: L1 5.78e-7 (32^2)
    // and 1.61e-7 (64^2) against RKL1's 1.04e-7 and 2.79e-8, and on the B-star gate box
    // dT_tau10 and dT_tau1 come out 7-8x the RKL1 arm's.  The two cures below both keep
    // the sweeps, the plane skip and the partitioned line solve exactly as they are and
    // only change WHICH sub-steps are taken:
    //
    // LODN (rad_adi_scheme = lodn, rad_adi_nsub = N): N sequential LOD sub-steps of
    // dts = tau/N, with the SWEEP ORDER ALTERNATED between sub-steps (x2-x3, then x3-x2,
    // ...), which cancels the leading Lie splitting error over each pair (Strang-like)
    // without a half-step.  Still first order in the backward-Euler sense, amplification
    // 1/(1+z/N)^N, so it converges to the exact e^-z only as N grows: the honest
    // brute-force fallback, and the comparison the extrapolation has to beat.
    //
    // LOD2 (rad_adi_scheme = lod2): RICHARDSON EXTRAPOLATION of LOD in the step size,
    //     y_full = LOD(tau) from 0,
    //     y_half = LOD(tau/2) from 0, then LOD(tau/2) again from y_half,
    //     y      = 2 y_half - y_full.
    // ALL THREE SUB-STEPS SWEEP IN THE SAME ORDER (x2 then x3), and that is not a detail:
    // Richardson needs y_full and y_half to be the SAME method at two step sizes, so that
    // their leading errors differ only by the factor 2.  Alternating the order between
    // the two halves (rad_adi_scheme = lod2a, kept as a switch) symmetrises the half
    // sequence and cancels its Lie SPLITTING error on its own -- which sounds better but
    // breaks the extrapolation: the full step's splitting error then has nothing to
    // cancel against and survives with the WRONG SIGN, so 2 y_half - y_full carries as
    // much splitting error as plain lod does.  MEASURED on the B-star gate box: lod2a
    // leaves dT_tau10 = 2.6e-6 and dT_tau1 = 6.0e-5, i.e. the lod arm's 2.3e-6/5.6e-5,
    // while RKL1 has 2.8e-7/7.5e-6.  The B-star residual is SPLITTING error, not
    // backward-Euler damping error.
    // LOD's error is O(dts) with a fixed leading coefficient, so the combination cancels
    // it and the scheme is SECOND order.  It is also still L-stable: one Fourier mode
    // with tau lambda = -z sees
    //     R(z) = 2/(1 + z/2)^2 - 1/(1 + z)  ->  0   as z -> infinity,
    // and |R| <= 1 on z >= 0 (R(z) - 1 = -z^2(3 + z + z^2/4)/((1+z/2)^2 (1+z)) <= 0 and
    // R(z) >= -1/8 at its minimum), so the stiff checkerboard is still damped, not
    // amplified -- unlike Douglas, whose R does not go to zero.  THE COST is three sweep
    // pairs per stage instead of one, with TWO factorisation sets (tau and tau/2): the
    // tau/2 tridiagonal differs from the tau one in every entry, so nothing is shared.
    //
    // CONSERVATION survives all of them: every sweep is a flux-form solve that conserves
    // sum_i V_i y_i exactly (see the note above), each sub-step therefore preserves it,
    // and 2 y_half - y_full is a linear combination with weights summing to 1 of two
    // vectors that each have sum V y = 0.  The report below measures exactly that.
    //
    // THE SCHEDULE.  One pass loop over (sub-step, direction), so pass p belongs to
    // sub-step p/nsweep and runs direction p%nsweep of it.
    const int nbm_ = tr_ared.extent_int(1);
    const int nlmax_ = tr_ared.extent_int(2);
    const int nsweep = three_d ? 2 : 1;
    int nsb = 1;
    if (scm == ADISCM_LOD2 || scm == ADISCM_LOD2A) {
      nsb = 3;
    } else if (scm == ADISCM_LODN) {
      nsb = (rad_adi_nsub > 1) ? rad_adi_nsub : 1;
    }
    std::vector<Real> sdt(nsb, tau);
    std::vector<int> smode(nsb, 1), sfrst(nsb, 1);
    // smode: 0 = this sub-step starts from y = 0, 1 = from what is in yy.
    // sfrst: 1 = x2 sweeps first in this sub-step, 0 = x3 first.
    if (scm == ADISCM_LOD2 || scm == ADISCM_LOD2A) {
      // the full step
      sdt[0] = tau;
      smode[0] = 0;
      sfrst[0] = 1;
      // the first half, from zero again
      sdt[1] = 0.5*tau;
      smode[1] = 0;
      sfrst[1] = 1;
      // the second half, continuing.  lod2 keeps the SAME sweep order as the full step,
      // which is what makes the extrapolation valid (see the note above); lod2a swaps it,
      // symmetrising the half sequence at the price of the cancellation.
      sdt[2] = 0.5*tau;
      smode[2] = 1;
      sfrst[2] = (scm == ADISCM_LOD2A) ? 0 : 1;
    } else {
      const Real dsb = tau/static_cast<Real>(nsb);
      for (int s=0; s<nsb; ++s) {
        sdt[s] = dsb;
        smode[s] = (s == 0) ? 0 : 1;
        sfrst[s] = ((s % 2) == 0) ? 1 : 0;
      }
    }
    // ---- THE OUTER CROSS-TERM ITERATION (cubed sphere only).  The metric cross term
    // is not in the tridiagonal, so one ADI step carries it frozen at T*.  That is
    // unconditionally stable (the cross/diagonal symbol ratio is <= 1/2, and -> 0 at the
    // Nyquist checkerboard), but a mid-band mode is left with amplification ~ -1/4
    // instead of ~0.  rad_adi_cross_iter > 1 re-solves with the cross term evaluated at
    // the previous answer, a fixed-point iteration with that same contraction factor.
    // ncrit = 1 off the cubed sphere, and then `lag` is never read.
    const int ncrit = (csx && rad_adi_cross_iter > 1) ? rad_adi_cross_iter : 1;
    auto lag = (ncrit > 1) ? tr_ylg : tr_ya;
    bool uselg = false;
    // the x1 window of every exchange this section makes, as in the refresh above
    const int sla_ = rad_tr_window ? is : -1, sua_ = rad_tr_window ? ie : -1;
    for (int cit=0; cit<ncrit; ++cit) {
    if (cit > 0) {
      Kokkos::deep_copy(lag, yy);
      pbval_tr->InitRecv(1);
      pbval_tr->PackAndSendCC(lag, tr_ycoar, sla_, sua_);
      while (pbval_tr->RecvAndUnpackCC(lag, tr_ycoar, sla_, sua_)
             != TaskStatus::complete) {}
      Kokkos::fence();
      while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
      while (pbval_tr->ClearSend() != TaskStatus::complete) {}
      uselg = true;
    }
    for (int pass=0; pass<nsb*nsweep; ++pass) {
      const int sb = pass/nsweep;
      const int dpass = pass - sb*nsweep;
      const Real dts = sdt[sb];
      const bool dfrst = three_d ? (sfrst[sb] != 0) : true;
      // LOD is a sequence of BACKWARD-EULER sub-steps: the implicit weight is 1 by
      // construction (rad_adi_theta is refused with it, see conduction.cpp)
      const Real thtau = adilod ? dts : (rad_adi_theta*dts);
      const bool d2 = ((dpass == 0) == dfrst);
      // LOD2: the full-step answer is stashed before the half-step sequence restarts
      if (dpass == 0 && sb == 1 && (scm == ADISCM_LOD2 ||
                                    scm == ADISCM_LOD2A)) {
        Kokkos::deep_copy(yf, yy);
      }
      // ---- THE RIGHT-HAND SIDE of this sweep.  The first sweep of a sub-step loads the
      // sub-step's starting value plus the source of the direction it integrates
      // (Douglas: the whole source); the second loads its own direction's source on top
      // of the first sub-problem's answer.  Over the FULL plane range, so an INACTIVE
      // plane -- which no sweep touches -- still accumulates the full tau b, the answer
      // the RKL1 path's single substage gives it, in every scheme (and the Richardson
      // combination of 2 x (tau/2) b and tau b is again tau b).
      if (dpass == 0 || adilod) {
        const bool ub2 = (dpass == 0) ? (adilod ? d2 : true) : d2;
        const bool ub3 = three_d && ((dpass == 0) ? (adilod ? !d2 : true) : !d2);
        const int lmode = (dpass == 0) ? smode[sb] : 2;    // 2 = add in place
        auto o_ = yy;
        auto lg_ = lag;
        const bool ulg_ = uselg;
        par_for("radtradirhs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
          const Real ai = st(m,ia_,k,j,i);
          const Real phc = st(m,it_,k,j,i);
          Real mi = 0.0;
          if (ub2) {
            const Real alm2 = st(m,ia_,k,j-1,i), alp2 = st(m,ia_,k,j+1,i);
            Real fl = (lnk2(m,j) && ai > 0.0 && alm2 > 0.0)
                      ? c2(m,k,j,i)*(phc - st(m,it_,k,j-1,i)) : 0.0;
            Real fr = (lnk2(m,j+1) && ai > 0.0 && alp2 > 0.0)
                      ? c2(m,k,j+1,i)*(st(m,it_,k,j+1,i) - phc) : 0.0;
            // the metric cross term of the SAME faces, explicit: it is the part of
            // M(0) the tridiagonal cannot carry (rad_adi_cross_iter re-evaluates it at
            // the current answer instead of at T*)
            if (csx) {
              if (lnk2(m,j) && ai > 0.0 && alm2 > 0.0) {
                fl -= cross2(lg_, ulg_, m, k, j, i);
              }
              if (lnk2(m,j+1) && ai > 0.0 && alp2 > 0.0) {
                fr -= cross2(lg_, ulg_, m, k, j+1, i);
              }
            }
            mi = fr - fl;
          }
          if (ub3) {
            const Real alm3 = st(m,ia_,k-1,j,i), alp3 = st(m,ia_,k+1,j,i);
            Real gl = (lnk3(m,k) && ai > 0.0 && alm3 > 0.0)
                    ? c3(m,k,j,i)*(phc - st(m,it_,k-1,j,i)) : 0.0;
            Real gr = (lnk3(m,k+1) && ai > 0.0 && alp3 > 0.0)
                    ? c3(m,k+1,j,i)*(st(m,it_,k+1,j,i) - phc) : 0.0;
            if (csx) {
              if (lnk3(m,k) && ai > 0.0 && alm3 > 0.0) {
                gl -= cross3(lg_, ulg_, m, k, j, i);
              }
              if (lnk3(m,k+1) && ai > 0.0 && alp3 > 0.0) {
                gr -= cross3(lg_, ulg_, m, k+1, j, i);
              }
            }
            const Real m3 = gr - gl;
            mi = ub2 ? (mi + m3) : m3;
          }
          if (curv) mi /= vcell(m,k,j,i);
          if (!isfinite(mi)) mi = 0.0;
          if (lmode == 0) {
            o_(m,0,k,j,i) = dts*mi;
          } else {
            o_(m,0,k,j,i) += dts*mi;   // lmode 1 (from y_in) and 2 (add) coincide
          }
        });
      }
      const int nb = d2 ? adi_nb2 : adi_nb3;
      const int ns = d2 ? js : ks;
      const int ne = d2 ? je : ke;
      const int t1s = d2 ? ks : js;
      const int t1e = d2 ? ke : je;
      const int nline = (t1e - t1s + 1)*nplact;
      const int doff = d2 ? 0 : nmbl;
      for (int m=0; m<nmbl; ++m) tr_ab0_h(m) = ab0h[doff + m];
      Kokkos::deep_copy(ab0, tr_ab0_h);
      auto y_ = yy;
      auto u_ = uu;
      auto v_ = vv;
      auto w_ = ww;
      auto cp_ = cpv;
      const int nplact_ = nplact, ilo2_ = ilo_, t1sv_ = t1s;

      // (1) the LOCAL Thomas factorisation, three right-hand sides at once, and the six
      // interface coefficients of this block's piece of every line
      par_for("radtradiln", DevExeSpace(), 0, nmb1, t1s, t1e, ilo_, ihi_,
      KOKKOS_LAMBDA(const int m, const int t, const int i) {
        const int l = (t - t1sv_)*nplact_ + (i - ilo2_);
        const int b0 = ab0(m);
        if (act(i - isv) == 0) {
          for (int c=0; c<6; ++c) rd_(m,b0,l,c) = 0.0;
          return;
        }
        Real a1 = 0.0, cn = 0.0;
        for (int s=ns; s<=ne; ++s) {
          const int k = d2 ? t : s;
          const int j = d2 ? s : t;
          const Real ai = st(m,ia_,k,j,i);
          Real cl, cr, alm, alp;
          if (d2) {
            alm = st(m,ia_,k,j-1,i);
            alp = st(m,ia_,k,j+1,i);
            cl = (lnk2(m,j) && ai > 0.0 && alm > 0.0) ? c2(m,k,j,i) : 0.0;
            cr = (lnk2(m,j+1) && ai > 0.0 && alp > 0.0) ? c2(m,k,j+1,i) : 0.0;
          } else {
            alm = st(m,ia_,k-1,j,i);
            alp = st(m,ia_,k+1,j,i);
            cl = (lnk3(m,k) && ai > 0.0 && alm > 0.0) ? c3(m,k,j,i) : 0.0;
            cr = (lnk3(m,k+1) && ai > 0.0 && alp > 0.0) ? c3(m,k+1,j,i) : 0.0;
          }
          // The row of (I - thtau A_d).  The three Cartesian statements are kept
          // EXACTLY as they were, in exactly this order, and the curvilinear form is
          // appended as an overwrite: A_d carries the 1/V_i of the flux divergence and
          // V_i is exactly 1.0 off a cubed sphere, but neither dividing by that 1.0 nor
          // hoisting the same expressions into an if/else is a compile no-op -- both
          // change how the product contracts into an FMA, and MEASURED on the He-star
          // FeCZ box (the production Cartesian configuration) that moves the last bit of
          // the tridiagonal row and, through it, the solution.  This is the one hunk of
          // this branch that broke the Cartesian bitwise regression; see tests_adi.
          Real dj = 1.0 + thtau*ai*(cl + cr);
          if (!(dj > 0.0) || !isfinite(dj)) dj = 1.0;
          Real aj = -thtau*cl*alm;
          Real cj = -thtau*cr*alp;
          if (curv) {
            const Real thv = thtau/vol_(m,k,j,i);
            dj = 1.0 + thv*ai*(cl + cr);
            if (!(dj > 0.0) || !isfinite(dj)) dj = 1.0;
            aj = -thv*cl*alm;
            cj = -thv*cr*alp;
          }
          if (!isfinite(aj)) aj = 0.0;
          if (!isfinite(cj)) cj = 0.0;
          Real e1 = 0.0, en = 0.0;
          // the two rows that couple OUT of this block: the coupling is remembered and
          // the row is truncated, which is what turns it into a spike right-hand side
          if (s == ns) {
            a1 = aj;
            aj = 0.0;
            e1 = 1.0;
          }
          if (s == ne) {
            cn = cj;
            cj = 0.0;
            en = 1.0;
          }
          const int kp = d2 ? k : (k - 1);
          const int jp = d2 ? (j - 1) : j;
          const Real cpm = (s == ns) ? 0.0 : cp_(m,0,kp,jp,i);
          Real bet = dj - aj*cpm;
          if (!(fabs(bet) > 0.0) || !isfinite(bet)) bet = 1.0;
          const Real pu = (s == ns) ? 0.0 : u_(m,0,kp,jp,i);
          const Real pv = (s == ns) ? 0.0 : v_(m,0,kp,jp,i);
          const Real pw = (s == ns) ? 0.0 : w_(m,0,kp,jp,i);
          u_(m,0,k,j,i) = (y_(m,0,k,j,i) - aj*pu)/bet;
          v_(m,0,k,j,i) = (e1 - aj*pv)/bet;
          w_(m,0,k,j,i) = (en - aj*pw)/bet;
          cp_(m,0,k,j,i) = cj/bet;
        }
        // the back substitution, with the two spikes scaled by their couplings as they
        // are written (v'_s = a_1 v_s - cp_s v'_{s+1} is the same recurrence)
        {
          const int k = d2 ? t : ne;
          const int j = d2 ? ne : t;
          v_(m,0,k,j,i) *= a1;
          w_(m,0,k,j,i) *= cn;
        }
        for (int s=ne-1; s>=ns; --s) {
          const int k = d2 ? t : s;
          const int j = d2 ? s : t;
          const int kn = d2 ? k : (k + 1);
          const int jn = d2 ? (j + 1) : j;
          const Real cc = cp_(m,0,k,j,i);
          u_(m,0,k,j,i) -= cc*u_(m,0,kn,jn,i);
          v_(m,0,k,j,i) = a1*v_(m,0,k,j,i) - cc*v_(m,0,kn,jn,i);
          w_(m,0,k,j,i) = cn*w_(m,0,k,j,i) - cc*w_(m,0,kn,jn,i);
        }
        const int k1 = d2 ? t : ns, j1 = d2 ? ns : t;
        const int k2 = d2 ? t : ne, j2 = d2 ? ne : t;
        rd_(m,b0,l,0) = v_(m,0,k1,j1,i);
        rd_(m,b0,l,1) = w_(m,0,k1,j1,i);
        rd_(m,b0,l,2) = u_(m,0,k1,j1,i);
        rd_(m,b0,l,3) = v_(m,0,k2,j2,i);
        rd_(m,b0,l,4) = w_(m,0,k2,j2,i);
        rd_(m,b0,l,5) = u_(m,0,k2,j2,i);
      });

      // (2) the RING GATHER.  Round r shifts one slab one block to the right, so after
      // nb-1 rounds every block of the ring holds all nb slabs, indexed by ABSOLUTE block
      // index.  A neighbour on this rank is a device copy; one on another rank is a
      // message on the module's own communicator, tagged with the RECEIVER's local id (a
      // rank can receive one slab per local block per round, so that is unique).
      //
      // ON THE CUBED SPHERE the chain is OPEN and a cyclic shift has nothing to walk, so
      // the same shift is run TWICE, once in each direction: pass 0 carries slabs
      // rightward (block b receives slab b-r from its left neighbour, which holds it
      // after round r-1), pass 1 carries them leftward (slab b+r from the right
      // neighbour).  A round is skipped where the slab index leaves [0,nb) or the face is
      // not linked, which is what makes the two ends of the chain ends.  After nb-1
      // rounds each way every block of the chain again holds all nb slabs, indexed by the
      // block's index WITHIN THE PANEL.  The extra tag bit (p << 4) is free: r < NADIB=8
      // uses bits 1..3 only.
      Real *rbase = rd_.data();
      const int chunk = nlmax_*6;
      for (int p=0; openchain && p<2; ++p) {
        const bool rt = (p == 0);
        for (int r=1; r<nb; ++r) {
          Kokkos::fence();
#if MPI_PARALLEL_ENABLED
          std::vector<MPI_Request> reqs;
#endif
          // what I receive, and from which side
          for (int m=0; m<nmbl; ++m) {
            const int b0 = ab0h[doff + m];
            const int br = rt ? (b0 - r) : (b0 + r);
            if (br < 0 || br >= nb) continue;
            const int srk = rt ? lrk[doff + m] : rrk[doff + m];
            const int slid = rt ? llid[doff + m] : rlid[doff + m];
            if (srk < 0) continue;
            if (srk == global_variable::my_rank) {
              Kokkos::deep_copy(
                Kokkos::subview(rd_, m, br, Kokkos::make_pair(0, nline), Kokkos::ALL),
                Kokkos::subview(rd_, slid, br, Kokkos::make_pair(0, nline), Kokkos::ALL));
            } else {
#if MPI_PARALLEL_ENABLED
              MPI_Request rq;
              const int tg = (m << 5) | (p << 4) | (r << 1) | (d2 ? 0 : 1);
              MPI_Irecv(rbase + (static_cast<std::size_t>(m)*nbm_ + br)*chunk,
                        6*nline, MPI_ATHENA_REAL, srk, tg, adi_comm, &rq);
              reqs.push_back(rq);
#endif
            }
          }
          // what I send, and to which side: the slab my downstream neighbour needs
          for (int m=0; m<nmbl; ++m) {
            const int b0 = ab0h[doff + m];
            const int bs = rt ? (b0 - r + 1) : (b0 + r - 1);
            if (bs < 0 || bs >= nb) continue;
            const int drk = rt ? rrk[doff + m] : lrk[doff + m];
            const int dlid = rt ? rlid[doff + m] : llid[doff + m];
            if (drk < 0 || drk == global_variable::my_rank) continue;
#if MPI_PARALLEL_ENABLED
            MPI_Request rq;
            const int tg = (dlid << 5) | (p << 4) | (r << 1) | (d2 ? 0 : 1);
            MPI_Isend(rbase + (static_cast<std::size_t>(m)*nbm_ + bs)*chunk,
                      6*nline, MPI_ATHENA_REAL, drk, tg, adi_comm, &rq);
            reqs.push_back(rq);
#endif
          }
#if MPI_PARALLEL_ENABLED
          if (!reqs.empty()) {
            MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
          }
#endif
        }
      }
      for (int r=1; !openchain && r<nb; ++r) {
        Kokkos::fence();
#if MPI_PARALLEL_ENABLED
        std::vector<MPI_Request> reqs;
#endif
        for (int m=0; m<nmbl; ++m) {
          const int b0 = ab0h[doff + m];
          const int br = ((b0 - r) % nb + nb) % nb;          // the slab I receive
          if (lrk[doff + m] == global_variable::my_rank) {
            const int ml = llid[doff + m];
            Kokkos::deep_copy(
              Kokkos::subview(rd_, m, br, Kokkos::make_pair(0, nline), Kokkos::ALL),
              Kokkos::subview(rd_, ml, br, Kokkos::make_pair(0, nline), Kokkos::ALL));
          } else {
#if MPI_PARALLEL_ENABLED
            MPI_Request rq;
            const int tg = (m << 5) | (r << 1) | (d2 ? 0 : 1);
            MPI_Irecv(rbase + (static_cast<std::size_t>(m)*nbm_ + br)*chunk,
                      6*nline, MPI_ATHENA_REAL, lrk[doff + m], tg, adi_comm, &rq);
            reqs.push_back(rq);
#endif
          }
        }
        for (int m=0; m<nmbl; ++m) {
          if (rrk[doff + m] == global_variable::my_rank) continue;
#if MPI_PARALLEL_ENABLED
          const int b0 = ab0h[doff + m];
          const int bs = ((b0 - r + 1) % nb + nb) % nb;
          MPI_Request rq;
          const int tg = (rlid[doff + m] << 5) | (r << 1) | (d2 ? 0 : 1);
          MPI_Isend(rbase + (static_cast<std::size_t>(m)*nbm_ + bs)*chunk,
                    6*nline, MPI_ATHENA_REAL, rrk[doff + m], tg, adi_comm, &rq);
          reqs.push_back(rq);
#endif
        }
#if MPI_PARALLEL_ENABLED
        if (!reqs.empty()) {
          MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
        }
#endif
      }
      Kokkos::fence();

      // (3) the REDUCED SYSTEM, assembled in absolute block order and solved redundantly
      // on every block of the ring, then the local back substitution
      par_for("radtradird", DevExeSpace(), 0, nmb1, t1s, t1e, ilo_, ihi_,
      KOKKOS_LAMBDA(const int m, const int t, const int i) {
        if (act(i - isv) == 0) return;
        const int l = (t - t1sv_)*nplact_ + (i - ilo2_);
        const int nu = 2*nb;
        Real mm[2*NADIB*(2*NADIB + 1)];
        for (int r=0; r<nu; ++r) {
          for (int c=0; c<=nu; ++c) mm[r*(nu+1) + c] = 0.0;
        }
        // accumulate, so that nb = 1 (the block is its own two neighbours) and nb = 2
        // (one neighbour on both sides) fall out of the same cyclic assembly
        for (int b=0; b<nb; ++b) {
          const int bm = (b + nb - 1)%nb, bp = (b + 1)%nb;
          const int r0 = 2*b, r1 = 2*b + 1;
          mm[r0*(nu+1) + 2*b]      += 1.0;
          mm[r0*(nu+1) + 2*bm + 1] += rd_(m,b,l,0);
          mm[r0*(nu+1) + 2*bp]     += rd_(m,b,l,1);
          mm[r0*(nu+1) + nu]        = rd_(m,b,l,2);
          mm[r1*(nu+1) + 2*b + 1]  += 1.0;
          mm[r1*(nu+1) + 2*bm + 1] += rd_(m,b,l,3);
          mm[r1*(nu+1) + 2*bp]     += rd_(m,b,l,4);
          mm[r1*(nu+1) + nu]        = rd_(m,b,l,5);
        }
        for (int c=0; c<nu; ++c) {
          int piv = c;
          Real best = fabs(mm[c*(nu+1) + c]);
          for (int r=c+1; r<nu; ++r) {
            const Real vr = fabs(mm[r*(nu+1) + c]);
            if (vr > best) {
              best = vr;
              piv = r;
            }
          }
          if (piv != c) {
            for (int q=c; q<=nu; ++q) {
              const Real tmp = mm[c*(nu+1) + q];
              mm[c*(nu+1) + q] = mm[piv*(nu+1) + q];
              mm[piv*(nu+1) + q] = tmp;
            }
          }
          Real dg = mm[c*(nu+1) + c];
          if (!(fabs(dg) > 0.0) || !isfinite(dg)) dg = 1.0;
          for (int r=c+1; r<nu; ++r) {
            const Real f = mm[r*(nu+1) + c]/dg;
            if (f == 0.0) continue;
            for (int q=c; q<=nu; ++q) mm[r*(nu+1) + q] -= f*mm[c*(nu+1) + q];
          }
        }
        Real zz[2*NADIB];
        for (int r=nu-1; r>=0; --r) {
          Real acc = mm[r*(nu+1) + nu];
          for (int q=r+1; q<nu; ++q) acc -= mm[r*(nu+1) + q]*zz[q];
          Real dg = mm[r*(nu+1) + r];
          if (!(fabs(dg) > 0.0) || !isfinite(dg)) dg = 1.0;
          zz[r] = acc/dg;
        }
        const int b0 = ab0(m);
        const int bm = (b0 + nb - 1)%nb, bp = (b0 + 1)%nb;
        Real xl = zz[2*bm + 1], xr = zz[2*bp];
        if (!isfinite(xl)) xl = 0.0;
        if (!isfinite(xr)) xr = 0.0;
        for (int s=ns; s<=ne; ++s) {
          const int k = d2 ? t : s;
          const int j = d2 ? s : t;
          Real yv = u_(m,0,k,j,i) - xl*v_(m,0,k,j,i) - xr*w_(m,0,k,j,i);
          if (!isfinite(yv)) yv = 0.0;
          y_(m,0,k,j,i) = yv;
        }
      });
    }
    // ---- LOD2: the RICHARDSON COMBINATION, y = 2 y_half - y_full.  A linear combination
    // with weights 2 - 1 = 1, so sum V y = 0 survives it exactly (both operands have it).
    if (scm == ADISCM_LOD2 || scm == ADISCM_LOD2A) {
      auto y_ = yy;
      auto f_ = yf;
      par_for("radtradirx", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real yv = 2.0*y_(m,0,k,j,i) - f_(m,0,k,j,i);
        if (!isfinite(yv)) yv = 0.0;
        y_(m,0,k,j,i) = yv;
      });
    }

    // ---- THE PANEL-SEAM SUB-STEP.  A tridiagonal line cannot cross a panel seam, so
    // the sweeps above left every panel face out of the system.  Adding it as a plain
    // explicit source would put an EXPLICIT face inside an implicit operator: for the
    // mode that is uniform on each side and jumps at the seam, every other face carries
    // nothing, the amplification is 1 - 2 z, and at the production z ~ 1e3 that is a
    // blow-up.  Each seam face is therefore solved as an ISOLATED TWO-CELL backward
    // Euler problem,
    //     delta = tau [ C_f (Th_B - Th_A) + s_f Q_f ]
    //             / ( 1 + tau C_f (alpha_A/V_A + alpha_B/V_B) ),    y_A += delta/V_A,
    // which is L-stable -- as tau -> infinity delta saturates at the pair's equilibrium
    // exchange and never overshoots -- and exactly conservative: the cell on the far
    // side forms the same expression with the roles swapped and gets -delta, up to the
    // along-seam RESAMPLE of the halo, which is the same error the explicit cubed-sphere
    // operator already makes at a seam.  (rad_ang_solver = sts needs none of this: RKL1
    // has no line structure and a panel face is an ordinary open face to it.)
    if (anyseam) {
      pbval_tr->InitRecv(1);
      pbval_tr->PackAndSendCC(yy, tr_ycoar, sla_, sua_);
      while (pbval_tr->RecvAndUnpackCC(yy, tr_ycoar, sla_, sua_)
             != TaskStatus::complete) {}
      Kokkos::fence();
      while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
      while (pbval_tr->ClearSend() != TaskStatus::complete) {}
      // a SNAPSHOT, so that a cell with two seam faces (a cube-vertex block) reads the
      // pre-step value on both of them and the exchange with each neighbour stays
      // antisymmetric
      Kokkos::deep_copy(uu, yy);
      auto y_ = yy;
      auto sn_ = uu;
      // WHERE THE SEAM FLUX IS EVALUATED.  The sub-step runs AFTER the sweeps, so
      // reading the post-sweep Th (w = 1) makes the seam flux one sub-step out of phase
      // with the interior -- a Lie splitting error of relative size ~ z, which shows up
      // as a 3x local residual at the CUBE VERTICES, where a cell has two seam faces.
      // w = 1/2 is the trapezoidal evaluation between the state the sweeps started from
      // (y = 0) and the one they ended at, i.e. second order in the sub-step, and it
      // leaves the ISOLATED pair (the sweeps move nothing) exactly backward Euler, so
      // the L-stability argument above is untouched.
      const Real wsm = rad_adi_seam_w;
      par_for("radtradiseam", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real ai = st(m,ia_,k,j,i);
        const Real va = vcell(m,k,j,i);
        if (!(ai > 0.0) || !(va > 0.0)) return;
        const Real tha = st(m,it_,k,j,i) + ai*wsm*sn_(m,0,k,j,i);
        Real dy = 0.0;
        for (int f=0; f<4; ++f) {
          if (f >= 2 && !three_d) break;
          int kn = k, jn = j;
          Real cf = 0.0, qf = 0.0, sgn = 1.0;
          if (f == 0) {
            if (!seam2(m,j)) continue;
            jn = j - 1;
            cf = c2(m,k,j,i);
            if (csx) qf = cross2(sn_, true, m, k, j, i);
          } else if (f == 1) {
            if (!seam2(m,j+1)) continue;
            jn = j + 1;
            cf = c2(m,k,j+1,i);
            sgn = -1.0;
            if (csx) qf = cross2(sn_, true, m, k, j+1, i);
          } else if (f == 2) {
            if (!seam3(m,k)) continue;
            kn = k - 1;
            cf = c3(m,k,j,i);
            if (csx) qf = cross3(sn_, true, m, k, j, i);
          } else {
            if (!seam3(m,k+1)) continue;
            kn = k + 1;
            cf = c3(m,k+1,j,i);
            sgn = -1.0;
            if (csx) qf = cross3(sn_, true, m, k+1, j, i);
          }
          const Real ab = st(m,ia_,kn,jn,i);
          const Real vb = vcell(m,kn,jn,i);
          if (!(cf > 0.0) || !(ab > 0.0) || !(vb > 0.0)) continue;
          const Real thb = st(m,it_,kn,jn,i) + ab*wsm*sn_(m,0,kn,jn,i);
          const Real den = 1.0 + tau*cf*(ai/va + ab/vb);
          const Real dlt = tau*(cf*(thb - tha) + sgn*qf)/den;
          if (isfinite(dlt)) dy += dlt/va;
        }
        if (dy != 0.0) y_(m,0,k,j,i) += dy;
      });
    }
    }
    ++sts_ncall;

    // ---- write the increment into the energy and measure what it did to the total
    Real esum = 0.0, eabs = 0.0, emax = 0.0;
    {
      auto y_ = yy;
      Kokkos::parallel_reduce("radtradiend",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
      KOKKOS_LAMBDA(const int &idx, Real &ssum, Real &sabs, Real &smax) {
        const int m = idx/nkji_;
        const int k = (idx - m*nkji_)/nji_ + ks;
        const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
        const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
        Real y = y_(m,0,k,j,i);
        if (!isfinite(y)) y = 0.0;
        u0(m,IEN,k,j,i) += y;
        const Real dv = curv ? vol_(m,k,j,i)
                           : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
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
    if (ang_lines < 400 && global_variable::my_rank == 0) {
      const bool tell = rad_ang_verbose && (ang_lines < 20 || viol > 1.0e-10);
      if (tell) {
        ++ang_lines;
        std::cout << "### rad_implicit_ang cycle " << pmy_pack->pmesh->ncycle
                  << " t = " << pmy_pack->pmesh->time
                  << ": adi " << ((scm == ADISCM_DOUGLAS) ? "douglas" :
                                  (scm == ADISCM_LOD2) ? "lod2" :
                                  (scm == ADISCM_LOD2A) ? "lod2a" :
                                  (scm == ADISCM_LODN) ? "lodn" : "lod")
                  << ", sweep pairs = " << nsb
                  << ", theta = " << (adilod ? 1.0 : rad_adi_theta)
                  << ", max z_i = " << zmax
                  << ", planes active/total = " << (anyact ? (ahi - alo + 1) : 0)
                  << "/" << nplane
                  << ", blocks/line = " << adi_nb2 << " x " << adi_nb3
                  << ", max |de| = " << emax
                  << ", max |dT*|/T* = " << tshift
                  << ", |sum V de|/sum V|de| = " << viol << std::endl;
      }
    }
    return;
  }

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
  // ---- the PER-PLANE substage counts.  ceil((sqrt(1+8R)-1)/2) is monotone in R, so
  // max_i s_i is exactly the nsub the global radius above gives, i.e. the loop below
  // still runs nsub substages -- but only the planes that need them do any work.
  int splmin = nsub, splsum = nsub*nplane;
  if (perpl) {
    splmin = nsub;
    splsum = 0;
    for (int p=0; p<nplane; ++p) {
      const Real rp = 0.5*(1.0 + rad_sts_margin)*tr_zpl_h(p);
      int sp = static_cast<int>(std::ceil(0.5*(std::sqrt(1.0 + 8.0*rp) - 1.0)));
      if (sp < 1) sp = 1;
      if (sp > rad_ang_maxit) sp = rad_ang_maxit;
      tr_spl_h(p) = sp;
      tr_w1pl_h(p) = 2.0/(static_cast<Real>(sp)*static_cast<Real>(sp) + sp);
      if (sp < splmin) splmin = sp;
      splsum += sp;
    }
    Kokkos::deep_copy(tr_spl, tr_spl_h);
    Kokkos::deep_copy(tr_w1pl, tr_w1pl_h);
  }

  // ---- THE ACTIVE PLANE RANGE OF EACH SUBSTAGE.  The loop below still runs s_max
  // substages, but at substage j only the planes with s_i >= j do any arithmetic.  Under
  // the transverse stencil the planes never talk to each other, so those planes can be
  // dropped from the KERNEL as well as from the arithmetic: the x1 index is the outermost
  // of the par_for range, and the active planes of a stiffness profile that varies
  // smoothly with depth are contiguous, so one min/max bracket per substage is enough.
  // (A bracket, not a compact list: a plane inside the bracket with s_i < j still enters
  // the kernel, but returns immediately.)
  //
  // What the finished planes used to pay is the register copy y_j = y_{j-1}, needed only
  // because the three RKL1 registers rotate.  It is removed by NOT writing them at all
  // and remembering where each plane's final value was left: the rotation is a fixed
  // 3-cycle, so the register written at substage j is always index 2 - (j-1) % 3 of
  // (tr_yc, tr_yb, tr_ya) -- see the write-out below.  No value changes, so this is
  // bitwise what the copies gave.
  std::vector<int> plo(nsub + 1, 0), phi(nsub + 1, nplane - 1);
  int pswept = nsub*nplane;
  if (perpl) {
    pswept = 0;
    for (int j=1; j<=nsub; ++j) {
      int lo = nplane, hi = -1;
      for (int p=0; p<nplane; ++p) {
        if (tr_spl_h(p) >= j) {
          if (p < lo) lo = p;
          hi = p;
        }
      }
      if (hi < lo) { lo = 0; hi = 0; }     // cannot happen: max_i s_i = nsub
      plo[j] = lo;
      phi[j] = hi;
      pswept += hi - lo + 1;
    }
  }

  // ---- the RKL1 loop on the increment.  Y_0 = 0 everywhere INCLUDING the ghosts, so
  // the first substage needs no exchange; every later one exchanges Y_{j-1} through the
  // module's own one-variable boundary object before the stencil reads it.
  auto ycur = tr_ya;   // Y_{j-1}
  auto yold = tr_yb;   // Y_{j-2}
  auto ynew = tr_yc;   // Y_j
  Kokkos::deep_copy(ycur, 0.0);
  Kokkos::deep_copy(yold, 0.0);
  // Y_0 = Y_{-1} = 0 over the WHOLE array, ghosts included: both start fully valid.
  int dcur = ng_, dold = ng_;
  for (int js_ = 1; js_ <= nsub; ++js_) {
    // the refill cadence.  hevery = 1 reproduces "exchange before every substage but the
    // first" exactly; hevery = N refills at js_ = 1 + N, 1 + 2N, ...
    const bool refill = (js_ > 1) && (((js_ - 1) % hevery) == 0);
    if (refill) {
      // post the receives, send, and spin on the unpack: RecvAndUnpackCC is the only
      // one of these that can legitimately come back incomplete (the MPI traffic of
      // this substage is all there is to overlap it with)
      // ... and only over THIS SUBSTAGE'S ACTIVE PLANE BRACKET.  The stencil below
      // reads the transverse ghosts of Y_{j-1} at i in [ilo_,ihi_] and nowhere else:
      // the x1 term (rad_sts_all only, which turns rad_sts_perplane off) reads i+-1 on
      // the block's OWN row, never a transverse ghost, and the write-out reads active
      // cells only.  The brackets NEST as j grows (the set {i : s_i >= j} shrinks), so a
      // plane dropped here is never read again from this register either.
      // rad_tr_window = false: (-1,-1) = the full x1 range, no window
      const int iwl_ = rad_tr_window ? (is + plo[js_]) : -1;
      const int iwu_ = rad_tr_window ? (is + phi[js_]) : -1;
      pbval_tr->InitRecv(1);
      pbval_tr->PackAndSendCC(ycur, tr_ycoar, iwl_, iwu_);
      while (pbval_tr->RecvAndUnpackCC(ycur, tr_ycoar, iwl_, iwu_)
             != TaskStatus::complete) {}
      // RecvAndUnpackCC launches the unpack kernel ASYNCHRONOUSLY and returns.  The
      // next InitRecv re-posts MPI_Irecv into the SAME device recv buffer, so with
      // GPU-aware MPI a neighbour that is ahead can have the NIC overwrite it while
      // the unpack kernel still reads it: fence before ClearRecv/ClearSend/InitRecv.
      Kokkos::fence();
      while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
      while (pbval_tr->ClearSend() != TaskStatus::complete) {}
      dcur = ng_;
      // hevery >= 3: Y_{j-2} as well, or the pointwise term of the recurrence runs the
      // skin out before the next refill (see the bookkeeping note above)
      if (xch2_) {
        pbval_tr->InitRecv(1);
        pbval_tr->PackAndSendCC(yold, tr_ycoar, iwl_, iwu_);
        while (pbval_tr->RecvAndUnpackCC(yold, tr_ycoar, iwl_, iwu_)
               != TaskStatus::complete) {}
        Kokkos::fence();
        while (pbval_tr->ClearRecv() != TaskStatus::complete) {}
        while (pbval_tr->ClearSend() != TaskStatus::complete) {}
        dold = ng_;
      }
    }
    // the depth this substage can legitimately write (0 = the active cells alone, which
    // is what hevery = 1 always gives)
    int dnew = 0;
    if (skin_) {
      const int dlim = (dcur - 1 < dold) ? (dcur - 1) : dold;
      dnew = (dlim < ng_ - 1) ? dlim : (ng_ - 1);
      if (dnew < 0) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "rad_tr_halo_every: the ghost skin ran out"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    const Real muj = (js_ == 1) ? 1.0 : (2.0*js_ - 1.0)/js_;
    const Real nuj = (js_ == 1) ? 0.0 : (1.0 - js_)/js_;
    const Real mut = muj*w1*tau;
    const int jsub = js_;
    auto yc_ = ycur;
    auto yo_ = yold;
    auto yn_ = ynew;
    const int ilo_ = is + plo[js_], ihi_ = is + phi[js_];
    // the ghost SKIN: a full frame of depth dnew in x2/x3 (0 = the active cells only).
    // Its own stencil inputs are valid because the frame is FULL -- a cell at
    // (j = je+d, k = ke+d) reads (je+d+1, ke+d) and (je+d, ke+d+1), which lie in the
    // depth-(d+1) frame the previous substage wrote.  The x1 window is untouched.
    const int jsk = js - dnew, jek = je + dnew;
    const int ksk = three_d ? (ks - dnew) : ks, kek = three_d ? (ke + dnew) : ke;
    par_for("radtrsub", DevExeSpace(), 0, nmb1, ksk, kek, jsk, jek, ilo_, ihi_,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      // the stiffness split left this block nothing to do: its increment stays zero,
      // but the register still has to be written, because the three registers rotate
      if (blkon && blk(m) == 0) { yn_(m,0,k,j,i) = 0.0; return; }
      // rad_sts_perplane: this plane's own s_i-stage scheme is already finished, so its
      // final increment is already sitting in the register substage s_i wrote it to, and
      // it is left there -- nothing is copied and nothing is written (the write-out below
      // picks the register up).  mu_j and nu_j do not depend on s; only mu~_j = mu_j w1
      // does, through w1 = 2/(s^2+s).
      Real mutp = mut;
      if (perpl) {
        if (jsub > spl(i - isv)) return;
        mutp = muj*w1pl(i - isv)*tau;
      }
      const Real ai = st(m,ia_,k,j,i);
      const Real thc = st(m,it_,k,j,i) + ai*yc_(m,0,k,j,i);
      // the FLUX through each face, written so that the neighbour forms the identical
      // expression from the identical operands
      const Real alm2 = st(m,ia_,k,j-1,i), alp2 = st(m,ia_,k,j+1,i);
      Real fl = (open2(m,j) && ai > 0.0 && alm2 > 0.0)
                ? c2(m,k,j,i)*(thc - (st(m,it_,k,j-1,i) + alm2*yc_(m,0,k,j-1,i))) : 0.0;
      Real fr = (open2(m,j+1) && ai > 0.0 && alp2 > 0.0)
                ? c2(m,k,j+1,i)*((st(m,it_,k,j+1,i) + alp2*yc_(m,0,k,j+1,i)) - thc) : 0.0;
      // the CUBED-SPHERE cross flux of the same two faces, at the SAME Th the diagonal
      // part uses -- so RKL1 integrates the exact cubed-sphere operator, cross term and
      // all, and the Gershgorin radius above covers the extra entries
      if (csx) {
        if (open2(m,j) && ai > 0.0 && alm2 > 0.0) fl -= cross2(yc_, true, m, k, j, i);
        if (open2(m,j+1) && ai > 0.0 && alp2 > 0.0) {
          fr -= cross2(yc_, true, m, k, j+1, i);
        }
      }
      Real mi = fr - fl;
      if (three_d) {
        const Real alm3 = st(m,ia_,k-1,j,i), alp3 = st(m,ia_,k+1,j,i);
        Real gl = (open3(m,k) && ai > 0.0 && alm3 > 0.0)
                ? c3(m,k,j,i)*(thc - (st(m,it_,k-1,j,i) + alm3*yc_(m,0,k-1,j,i))) : 0.0;
        Real gr = (open3(m,k+1) && ai > 0.0 && alp3 > 0.0)
                ? c3(m,k+1,j,i)*((st(m,it_,k+1,j,i) + alp3*yc_(m,0,k+1,j,i)) - thc) : 0.0;
        if (csx) {
          if (open3(m,k) && ai > 0.0 && alm3 > 0.0) gl -= cross3(yc_, true, m, k, j, i);
          if (open3(m,k+1) && ai > 0.0 && alp3 > 0.0) {
            gr -= cross3(yc_, true, m, k+1, j, i);
          }
        }
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
      // the face sum is a sum of TOTAL face fluxes; the cell volume divides it here,
      // outside the sum, so the flux form still cancels between the two cells of a face
      if (curv) mi /= vcell(m,k,j,i);
      if (!isfinite(mi)) mi = 0.0;
      yn_(m,0,k,j,i) = muj*yc_(m,0,k,j,i) + nuj*yo_(m,0,k,j,i) + mutp*mi;
    });
    auto tmp = yold;
    yold = ycur;
    ycur = ynew;
    ynew = tmp;
    dold = dcur;
    dcur = dnew;
  }

  // ---- write the increment into the energy, and measure what it did to the total
  Real esum = 0.0, eabs = 0.0, emax = 0.0;
  {
    const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
    const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;
    auto yfin = ycur;
    // rad_sts_perplane: each plane's final increment is in the register its OWN last
    // substage wrote.  The rotation (old <- cur, cur <- new, new <- old) has period 3
    // starting from new = tr_yc, so substage j wrote register 2 - (j-1) % 3 of
    // (tr_ya, tr_yb, tr_yc).  For s_i = nsub this reduces to ycur, which is what the
    // global path uses.
    auto yf0_ = tr_ya;
    auto yf1_ = tr_yb;
    auto yf2_ = tr_yc;
    Kokkos::parallel_reduce("radtrend",
    Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
    KOKKOS_LAMBDA(const int &idx, Real &ssum, Real &sabs, Real &smax) {
      const int m = idx/nkji_;
      const int k = (idx - m*nkji_)/nji_ + ks;
      const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
      const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
      Real y;
      if (perpl) {
        const int b = 2 - (spl(i - isv) - 1)%3;
        y = (b == 2) ? yf2_(m,0,k,j,i) : ((b == 1) ? yf1_(m,0,k,j,i)
                                                   : yf0_(m,0,k,j,i));
      } else {
        y = yfin(m,0,k,j,i);
      }
      if (!isfinite(y)) y = 0.0;
      u0(m,IEN,k,j,i) += y;
      const Real dv = curv ? vol_(m,k,j,i)
                           : size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
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
      // rad_sts_perplane: nsub above is then s_max, and these are the rest of the
      // distribution -- the whole point of the switch is that they sit well below it
      std::string plstr;
      if (perpl) {
        plstr = ", per-plane s min/mean = " + std::to_string(splmin) + "/"
                + std::to_string(static_cast<double>(splsum)
                                 /static_cast<double>(nplane))
                + ", planes swept/needed/full = " + std::to_string(pswept) + "/"
                + std::to_string(splsum) + "/" + std::to_string(nsub*nplane);
      }
      std::cout << (sts1 ? "### rad_sts_all cycle " : "### rad_implicit_ang cycle ")
                << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time
                << ": max z_i = " << zmax << ", substages = " << nsub
                << ", substages_total = " << sts_nsub_tot
                << " in " << sts_ncall << " calls"
                << (clamped ? " (CLAMPED at rad_ang_maxit -- the step is NOT covered)"
                            : "")
                << plstr
                << ", max |de| = " << emax
                << ", max |dT*|/T* = " << tshift
                << ", |sum V de|/sum V|de| = " << viol << std::endl;
    }
  }
  // the per-plane substage PROFILE, once: it is the map of where the operator is stiff,
  // and it is what says whether a taper is buying anything.  One block of nx1 integers.
  if (perpl && rad_ang_verbose && sts_ncall == 1 && global_variable::my_rank == 0) {
    std::cout << "### rad_sts_perplane s_i (i = is .. ie):";
    for (int p=0; p<nplane; ++p) std::cout << " " << tr_spl_h(p);
    std::cout << std::endl;
  }
  return;
}
