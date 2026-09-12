//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_fofc.cpp
//! \brief Implements functions for first-order flux correction (FOFC) algorithm.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/gnomonic_kernels.hpp"
#include "coordinates/gnomonic_raisevel_mhd.hpp"
#include "eos/eos.hpp"
#include "mhd/rsolvers/llf_mhd_singlestate.hpp"
#include "mhd/rsolvers/hlle_mhd_singlestate.hpp"
#include "mhd.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn void MHD::FOFC
//! \brief Implements first-order flux-correction (FOFC) algorithm for MHD.  First an
//! estimate of the updated conserved variables is made. This estimate is then used to
//! flag any cell where floors will be required during the conversion to primitives. Then
//! the fluxes on the faces of flagged cells are replaced with first-order LLF fluxes.
//! Often this is enough to prevent floors from being needed.  The FOFC infrastructure is
//! also exploited for BH excision.  If a cell is about the horizon, FOFC is automatically
//! triggered (without estimating updated conserved variables).

void MHD::FOFC(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
  int js = indcs.js, je = indcs.je, nx2 = indcs.nx2;
  int ks = indcs.ks, ke = indcs.ke, nx3 = indcs.nx3;

  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  int nmb = pmy_pack->nmb_thispack;
  auto flx1 = uflx.x1f;
  auto flx2 = uflx.x2f;
  auto flx3 = uflx.x3f;
  auto &size = pmy_pack->pmb->mb_size;

  auto &bcc0_ = bcc0;
  auto &e3x1_ = e3x1;
  auto &e2x1_ = e2x1;
  auto &e1x2_ = e1x2;
  auto &e3x2_ = e3x2;
  auto &e2x3_ = e2x3;
  auto &e1x3_ = e1x3;

  // CURVILINEAR GRIDS.  Both spherical grids store the flux divergence in conservative
  // form with FACE AREAS and CELL VOLUMES from Coordinates -- see MHD::RKUpdate, which
  // the trial update below mirrors exactly.  Dividing by the Cartesian mb_size.dx1/2/3
  // instead is silently wrong on spherical polar (the error is a factor ~r).
  //
  // GEOMETRIC SOURCE TERMS ARE DELIBERATELY LEFT OUT, exactly as RKUpdate leaves them
  // out: on these grids they are added afterwards by Coordinates::CoordSrcTerms, from
  // MHDSrcTerms.  The published FOFC semantics are "trial update from the Riemann
  // fluxes alone", and the trial state is used only to DECIDE whether a cell needs
  // floors, never stored.  Including the geometric source would make the test state
  // differ from the one RKUpdate produces at that point in the stage and, near the
  // polar axis where the cot(theta) term is largest, would flag cells on the strength
  // of a term the flux replacement cannot influence.
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  const bool sp_fofc_ = pmy_pack->pmesh->use_spherical_polar;
  const bool curv_fofc_ = sp_fofc_ || cs_;
  // The gnomonic rotations the single-state solves below need.  Three (sin,cos) pairs:
  // the x1 sweep uses the CELL-centred pair, the x2 sweep the xi-FACE pair (staggered in
  // j), the x3 sweep the eta-FACE pair (staggered in k) -- the same staggering the
  // high-order sweeps in mhd_fluxes.cpp use.  See the note in gnomonic_kernels.hpp about
  // the constant-memory budget if this kernel ever grows.
  const GnomonicTrig gtc_  = pmy_pack->pcoord->GnomonicTrigCell();
  const GnomonicTrig gtxi_ = pmy_pack->pcoord->GnomonicTrigFaceXi();
  const GnomonicTrig gtet_ = pmy_pack->pcoord->GnomonicTrigFaceEta();
  auto &vol_ = pmy_pack->pcoord->volume;
  auto &ar1_ = pmy_pack->pcoord->area.x1f;
  auto &ar2_ = pmy_pack->pcoord->area.x2f;
  auto &ar3_ = pmy_pack->pcoord->area.x3f;
  auto &spx1v_ = pmy_pack->pcoord->x1v;
  auto &spx2v_ = pmy_pack->pcoord->x2v;
  auto &spx3v_ = pmy_pack->pcoord->x3v;
  auto &spx1f_ = pmy_pack->pcoord->xx1f;
  auto &spx2f_ = pmy_pack->pcoord->xx2f;
  auto &spx3f_ = pmy_pack->pcoord->xx3f;

  if (use_fofc) {
    Real &gam0 = pdriver->gam0[stage-1];
    Real &gam1 = pdriver->gam1[stage-1];
    Real beta_dt = (pdriver->beta[stage-1])*(pmy_pack->pmesh->dt);

    int &nmhd_ = nmhd;
    auto &u0_ = u0;
    auto &u1_ = u1;
    auto &utest_ = utest;
    auto &bcctest_ = bcctest;
    auto &b1_ = b1;
    const bool cs_fofc_ = cs_;
    auto &dx2_fofc_ = pmy_pack->pcoord->dx2;
    auto &dx3_fofc_ = pmy_pack->pcoord->dx3;
    auto &ccs_fofc_ = pmy_pack->pcoord->cos_cell;
    auto &scs_fofc_ = pmy_pack->pcoord->sin_cell;
    auto &x1v_fofc_ = pmy_pack->pcoord->x1v;
    auto &x1f_fofc_ = pmy_pack->pcoord->xx1f;

    // Index bounds
    int il = is-1, iu = ie+1, jl = js, ju = je, kl = ks, ku = ke;
    if (multi_d) { jl = js-1, ju = je+1; }
    if (three_d) { kl = ks-1, ku = ke+1; }

    // Estimate updated conserved variables and cell-centered fields
    par_for("FOFC-newu", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dtodx1 = beta_dt/size.d_view(m).dx1;
      Real dtodx2 = beta_dt/size.d_view(m).dx2;
      Real dtodx3 = beta_dt/size.d_view(m).dx3;

      // Estimate conserved variables
      if (curv_fofc_) {
        // area/volume divergence: the same form as MHD::RKUpdate.  At a POLAR axis the
        // x2 face area is identically zero, so the pole face contributes nothing here,
        // exactly as in the real update.
        const Real dtodv = beta_dt/vol_(m,k,j,i);
        for (int n=0; n<nmhd_; ++n) {
          Real divf = dtodv*(flx1(m,n,k,j,i+1)*ar1_(m,k,j,i+1) -
                             flx1(m,n,k,j,i  )*ar1_(m,k,j,i  ));
          if (multi_d) {
            divf += dtodv*(flx2(m,n,k,j+1,i)*ar2_(m,k,j+1,i) -
                           flx2(m,n,k,j  ,i)*ar2_(m,k,j  ,i));
          }
          if (three_d) {
            divf += dtodv*(flx3(m,n,k+1,j,i)*ar3_(m,k+1,j,i) -
                           flx3(m,n,k  ,j,i)*ar3_(m,k  ,j,i));
          }
          utest_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - divf;
        }
      } else {
      for (int n=0; n<nmhd_; ++n) {
        Real divf = dtodx1*(flx1(m,n,k,j,i+1) - flx1(m,n,k,j,i));
        if (multi_d) {
          divf += dtodx2*(flx2(m,n,k,j+1,i) - flx2(m,n,k,j,i));
        }
        if (three_d) {
          divf += dtodx3*(flx3(m,n,k+1,j,i) - flx3(m,n,k,j,i));
        }
        utest_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - divf;
      }
      }

      // Estimate updated cell-centered fields
      Real b1old = 0.5*(b1_.x1f(m,k,j,i) + b1_.x1f(m,k,j,i+1));
      Real b2old = 0.5*(b1_.x2f(m,k,j,i) + b1_.x2f(m,k,j+1,i));
      Real b3old = 0.5*(b1_.x3f(m,k,j,i) + b1_.x3f(m,k+1,j,i));
      if (sp_fofc_) {
        // SPHERICAL POLAR: the radial, polar and azimuthal cell centres are not the
        // midpoints of their faces on a stretched grid, and the real ConsToPrim
        // (ideal_mhd.cpp / general_mhd.cpp) reconstructs bcc with the linear weights
        // below.  The plain 0.5 averages disagree with it, and the difference is a
        // magnetic energy the floor test would then see but the real pass would not.
        Real lw = (spx1f_(m,i+1)-spx1v_(m,i))/(spx1f_(m,i+1)-spx1f_(m,i));
        Real rw = (spx1v_(m,i)-spx1f_(m,i))/(spx1f_(m,i+1)-spx1f_(m,i));
        b1old = lw*b1_.x1f(m,k,j,i) + rw*b1_.x1f(m,k,j,i+1);
        lw = (spx2f_(m,j+1)-spx2v_(m,j))/(spx2f_(m,j+1)-spx2f_(m,j));
        rw = (spx2v_(m,j)-spx2f_(m,j))/(spx2f_(m,j+1)-spx2f_(m,j));
        b2old = lw*b1_.x2f(m,k,j,i) + rw*b1_.x2f(m,k,j+1,i);
        lw = (spx3f_(m,k+1)-spx3v_(m,k))/(spx3f_(m,k+1)-spx3f_(m,k));
        rw = (spx3v_(m,k)-spx3f_(m,k))/(spx3f_(m,k+1)-spx3f_(m,k));
        b3old = lw*b1_.x3f(m,k,j,i) + rw*b1_.x3f(m,k+1,j,i);
      }
      if (cs_fofc_) {
        // CUBED SPHERE: bcc0 is stored in the orthonormal frame with the radial
        // component at the centroid (GnomonicEquiangleRaiseVelMHD); the plain averages
        // are components on the non-orthogonal tangent pair.  Mixing the two in bcctest
        // mis-states the test state's magnetic energy by the cross term.
        b1old = CellCenteredRadialFld(b1_.x1f(m,k,j,i), b1_.x1f(m,k,j,i+1),
                                      x1f_fofc_(m,i), x1f_fofc_(m,i+1), x1v_fofc_(m,i));
        const Real c = ccs_fofc_(m,k,j), sn = scs_fofc_(m,k,j);
        b2old = (b2old + c*b3old)/sn;
      }

      bcctest_(m,IBX,k,j,i) = gam0*bcc0_(m,IBX,k,j,i) + gam1*b1old;
      bcctest_(m,IBY,k,j,i) = gam0*bcc0_(m,IBY,k,j,i) + gam1*b2old;
      bcctest_(m,IBZ,k,j,i) = gam0*bcc0_(m,IBZ,k,j,i) + gam1*b3old;

      if (sp_fofc_) {
        // SPHERICAL POLAR.  dB/dt = -curl E with the orthogonal scale factors h1 = 1,
        // h2 = r, h3 = r sin(theta), i.e. the same operator the Cartesian lines below
        // apply, generalised.  This is the CELL-CENTRED counterpart of the face-centred
        // CT update in mhd_ct.cpp (which carries the same weights as area/dxedge); the
        // Cartesian form used here before was low by a factor ~r in every transverse
        // term and by the sin(theta) weighting in the x2 term.
        //   dB1/dt = -(1/(r sin))d_th(sin E3) + (1/(r sin))d_ph E2
        //   dB2/dt = -(1/(r sin))d_ph E1      + (1/r)d_r(r E3)
        //   dB3/dt = -(1/r)d_r(r E2)          + (1/r)d_th E1
        const Real rc = spx1v_(m,i);
        const Real snc = sin(spx2v_(m,j));
        const Real rfl = spx1f_(m,i), rfr = spx1f_(m,i+1);
        const Real dtordr = beta_dt/(rc*(rfr - rfl));
        bcctest_(m,IBY,k,j,i) += dtordr*(rfr*e3x1_(m,k,j,i+1) - rfl*e3x1_(m,k,j,i));
        bcctest_(m,IBZ,k,j,i) -= dtordr*(rfr*e2x1_(m,k,j,i+1) - rfl*e2x1_(m,k,j,i));
        if (multi_d) {
          const Real dth = spx2f_(m,j+1) - spx2f_(m,j);
          const Real snl = sin(spx2f_(m,j)), snr = sin(spx2f_(m,j+1));
          bcctest_(m,IBX,k,j,i) -= (beta_dt/(rc*snc*dth))*
              (snr*e3x2_(m,k,j+1,i) - snl*e3x2_(m,k,j,i));
          bcctest_(m,IBZ,k,j,i) += (beta_dt/(rc*dth))*
              (e1x2_(m,k,j+1,i) - e1x2_(m,k,j,i));
        }
        if (three_d) {
          const Real dtordp = beta_dt/(rc*snc*(spx3f_(m,k+1) - spx3f_(m,k)));
          bcctest_(m,IBX,k,j,i) += dtordp*(e2x3_(m,k+1,j,i) - e2x3_(m,k,j,i));
          bcctest_(m,IBY,k,j,i) -= dtordp*(e1x3_(m,k+1,j,i) - e1x3_(m,k,j,i));
        }
      } else if (cs_fofc_) {
        // CUBED SPHERE.  There is no cell-centred curl on this grid that MHD::CT can be
        // read off: CT updates the FACE fields from the EDGE (corner-averaged) EMFs and
        // GnomonicEquiangleRaiseVelMHD then rebuilds bcc from those faces, and the
        // corner EMFs do not exist yet when FOFC runs.  This is the cheapest estimate
        // that is CONSISTENT with the real update, and it is used only to decide whether
        // a cell needs floors -- never stored:
        //   * the two angular basis vectors both have length r, so the RADIAL terms take
        //     exactly the spherical-polar form (1/r) d(r E)/dr, which is exact here;
        //   * the two ANGULAR terms use the PHYSICAL edge lengths dx2/dx3 = r*dtheta
        //     that Coordinates stores (the Cartesian branch above divides by the
        //     coordinate widths mb_size.dx2/dx3, which on this grid are the raw
        //     panel-angle deltas and are wrong by a factor ~r).  What is NEGLECTED is
        //     the variation of each edge length along the transverse angle -- the
        //     analogue of the sin(theta) weighting in the spherical-polar branch -- and
        //     the O(cos_cell) cross term of the non-orthogonal tangent pair.  Both are
        //     corrections to a term that is itself one stage's change in B.
        const Real rc = x1v_fofc_(m,i);
        const Real rfl = x1f_fofc_(m,i), rfr = x1f_fofc_(m,i+1);
        const Real dtordr = beta_dt/(rc*(rfr - rfl));
        bcctest_(m,IBY,k,j,i) += dtordr*(rfr*e3x1_(m,k,j,i+1) - rfl*e3x1_(m,k,j,i));
        bcctest_(m,IBZ,k,j,i) -= dtordr*(rfr*e2x1_(m,k,j,i+1) - rfl*e2x1_(m,k,j,i));
        if (multi_d) {
          const Real dtod2 = beta_dt/dx2_fofc_(m,k,j,i);
          bcctest_(m,IBX,k,j,i) -= dtod2*(e3x2_(m,k,j+1,i) - e3x2_(m,k,j,i));
          bcctest_(m,IBZ,k,j,i) += dtod2*(e1x2_(m,k,j+1,i) - e1x2_(m,k,j,i));
        }
        if (three_d) {
          const Real dtod3 = beta_dt/dx3_fofc_(m,k,j,i);
          bcctest_(m,IBX,k,j,i) += dtod3*(e2x3_(m,k+1,j,i) - e2x3_(m,k,j,i));
          bcctest_(m,IBY,k,j,i) -= dtod3*(e1x3_(m,k+1,j,i) - e1x3_(m,k,j,i));
        }
      } else {
      bcctest_(m,IBY,k,j,i) += dtodx1*(e3x1_(m,k,j,i+1) - e3x1_(m,k,j,i));
      bcctest_(m,IBZ,k,j,i) -= dtodx1*(e2x1_(m,k,j,i+1) - e2x1_(m,k,j,i));
      if (multi_d) {
        bcctest_(m,IBX,k,j,i) -= dtodx2*(e3x2_(m,k,j+1,i) - e3x2_(m,k,j,i));
        bcctest_(m,IBZ,k,j,i) += dtodx2*(e1x2_(m,k,j+1,i) - e1x2_(m,k,j,i));
      }
      if (three_d) {
        bcctest_(m,IBX,k,j,i) += dtodx3*(e2x3_(m,k+1,j,i) - e2x3_(m,k,j,i));
        bcctest_(m,IBY,k,j,i) -= dtodx3*(e1x3_(m,k+1,j,i) - e1x3_(m,k,j,i));
      }
      }
    });

    // Test whether conversion to primitives requires floors
    // Note b0 and w0 passed to function, but not used/changed.
    peos->ConsToPrim(utest_, b0, w0, bcctest_, true, il, iu, jl, ju, kl, ku);

    // CUBED SPHERE: the DEFERRED half of the floors, exactly as in Hydro::FOFC.
    // ConsToPrim above ran under EOS_Data::defer_cons_floors, i.e. it tested the
    // pressure/temperature/energy floors and the velocity ceiling against an ORTHONORMAL
    // kinetic energy, which is NOT the kinetic energy on this grid (and against a
    // magnetic energy that is a sum of squares on a non-orthogonal triple).  The real
    // decision is taken afterwards by Coordinates::GnomonicEquiangleRaiseVelMHD, and the
    // floor-TEST pass is followed by no such call, so a cell that needs one of those
    // floors was invisible to FOFC.  Run exactly that arithmetic on the trial state here
    // -- the same free function the real pass uses,
    // coordinates/gnomonic_raisevel_mhd.hpp -- and flag whatever it would floor.
    // Nothing is written back: the trial state is scratch.  bcctest_ already holds the
    // ORTHONORMAL cell-centred field (built above the way RaiseVelMHD builds it), so the
    // magnetic energy it forms is the right one.  Cells the DENSITY floor already
    // flagged are skipped.
    if (cs_ && peos->eos_data.defer_cons_floors &&
        !pmy_pack->pcoord->is_special_relativistic &&
        !pmy_pack->pcoord->is_general_relativistic) {
      auto eos_ = peos->eos_data;
      const bool gen_ = eos_.IsGeneral();
      const Real vceil_ = eos_.vceil;
      auto fofcc_ = fofc;
      auto cosc_ = pmy_pack->pcoord->cos_cell;
      auto wtemp_ = wtemp;
      const int ni   = (iu - il + 1);
      const int nji  = (ju - jl + 1)*ni;
      const int nkji = (ku - kl + 1)*nji;
      int nflag_ = 0;
      Kokkos::parallel_reduce("FOFC-csfloor",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
      KOKKOS_LAMBDA(const int &idx, int &sumf) {
        const int m = idx/nkji;
        const int k = (idx - m*nkji)/nji + kl;
        const int j = (idx - m*nkji - (k - kl)*nji)/ni + jl;
        const int i = (idx - m*nkji - (k - kl)*nji - (j - jl)*ni) + il;
        if (fofcc_(m,k,j,i)) { return; }
        const Real d = utest_(m,IDN,k,j,i);
        if (!(d > 0.0)) { return; }
        Real m1 = utest_(m,IM1,k,j,i);
        Real m2 = utest_(m,IM2,k,j,i);
        Real m3 = utest_(m,IM3,k,j,i);
        Real etot = utest_(m,IEN,k,j,i);
        Real v1, v2, v3, eint, pnew, g1new, temp;
        bool ceil_used, floored;
        GnomonicRaiseVelMHDFloors(cosc_(m,k,j), eos_, gen_, vceil_, true,
                                  gen_ ? wtemp_(m,k,j,i) : 0.0, d,
                                  bcctest_(m,IBX,k,j,i), bcctest_(m,IBY,k,j,i),
                                  bcctest_(m,IBZ,k,j,i),
                                  m1, m2, m3, etot, v1, v2, v3, eint, pnew, g1new, temp,
                                  ceil_used, floored);
        if (ceil_used || floored) {
          fofcc_(m,k,j,i) = true;
          sumf++;
        }
      }, Kokkos::Sum<int>(nflag_));
      pmy_pack->pmesh->ecounter.nfofc += nflag_;
    }
    // Accumulate the PER-CELL FOFC flag count for the `mhd_fofc` output variable.
    // Placed after every writer of the flag array and before the flux kernel that reads
    // and clears it, so it sees exactly the set of cells this stage flagged -- the same
    // set the scalar `fofc` column of the event-log output counts.  It reads the flags
    // and writes only its own array: where and how the correction is applied is
    // untouched.  Reset to zero when the output variable is loaded (derived_variables).
    {
      auto fofcf_ = fofc;
      auto fcnt_ = fofc_cnt;
      par_for("FOFC-count", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        if (fofcf_(m,k,j,i)) { fcnt_(m,k,j,i) += 1.0; }
      });
    }
  }

  auto &coord = pmy_pack->pcoord->coord_data;
  bool &is_sr = pmy_pack->pcoord->is_special_relativistic;
  bool &is_gr = pmy_pack->pcoord->is_general_relativistic;
  auto &eos = peos->eos_data;
  // derived thermodynamic variables (general EOS only; empty if ideal)
  auto &wder_ = wder;
  auto &use_fofc_ = use_fofc;
  auto fofc_ = fofc;
  auto &use_excise_ = pmy_pack->pcoord->coord_data.bh_excise;
  auto &excision_flux_ = pmy_pack->pcoord->excision_flux;
  auto &w0_ = w0;
  // During MHD::Fluxes the well-balanced PERTURBATION reconstruction has already
  // subtracted the static background from w0 (RemoveWbVar in MHD::Fluxes), so w0 is NOT
  // a full state in this window.  The first-order fallback is a physical Riemann solve
  // and needs the full state, so the background is added back here, for the few flagged
  // cells only.  With the perturbation scheme off this is the identity on w0, so every
  // non-well-balanced run is bit-for-bit unchanged.  w0wb carries the nmhd hydrodynamic
  // channels only -- the magnetic field in bcc0 is never perturbed -- so the accessor
  // falls through for the passive scalars exactly as the hydro one does.
  auto &w0wb_ = w0wb;
  const bool wbpert_ = use_wellbalance_static_reconst_perturb;
  // <mhd>/fofc_rsolver: see the note in mhd.hpp.
  const bool fofc_hlle_ = fofc_hlle;
  const int nmhd_f = nmhd;
  const int nvar_f = nmhd + nscalars;
  auto &b0_ = b0;

  // Index bounds
  int il = is-1, iu = ie+1, jl = js, ju = je, kl = ks, ku = ke;
  if (multi_d) { jl = js-1, ju = je+1; }
  if (three_d) { kl = ks-1, ku = ke+1; }

  // Replace fluxes with first-order LLF fluxes at i,j,k faces for any cell where FOFC
  // and/or excision is used (if GR+excising)
  par_for("FOFC-flx", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // full-state primitive accessor; see the note on the WB perturbation above
    auto wfull = [&](const int n, const int kk, const int jj, const int ii) -> Real {
      return (wbpert_ && n < nmhd_f) ? (w0_(m,n,kk,jj,ii) + w0wb_(m,n,kk,jj,ii))
                                     : w0_(m,n,kk,jj,ii);
    };

    // CUBED SPHERE.  The two panel-tangential basis vectors are not orthogonal, so a
    // Riemann solve must be handed the velocity AND the field in the sweep's locally
    // ORTHONORMAL frame, and the momentum flux it returns must be rotated back and its
    // index LOWERED -- exactly what GnomonicEquianglePrimFaceX* / GnomonicEquiangle-
    // FaceBX2 / GnomonicEquiangleFluxX* / GnomonicEquiangleEmfX1 do for the high-order
    // sweeps (coordinates/gnomonic_kernels.hpp).  Without this the first-order fallback
    // writes momentum fluxes and face EMFs in the wrong basis, which is an inconsistency
    // and not merely a loss of accuracy.  The radial component is orthogonal to both
    // angles and is never touched.  Off the cubed sphere every one of these is the
    // identity, so all other grids are bit-for-bit unchanged.
    //
    // bcc0 is stored in the orthonormal frame {rhat, e_xi, (e_eta - c e_xi)/s} (see
    // Coordinates::GnomonicEquiangleRaiseVelMHD), and each sweep's frame is {nhat, then
    // its two face-parallel axes} with the NORMAL component taken from b0.x*f: so x1 and
    // x3 need no field rotation at all, and only the x2 sweep's "by" slot does --
    // B.e_eta = c*(B.e_xi) + s*(B.f2).
    //
    // The face index `f` selects the trig: the x1 sweep uses the CELL-centred pair (it
    // depends on (k,j) only, so both x1 faces of a cell share it), the x2 sweep the
    // xi-FACE pair at the face's own j, the x3 sweep the eta-FACE pair at the face's own
    // k -- the same staggering the high-order sweeps use.
    auto ldx1 = [&](const int ii, MHDPrim1D &q) {
      q.d = wfull(IDN,k,j,ii);
      Real q1 = wfull(IVX,k,j,ii), q2 = wfull(IVY,k,j,ii), q3 = wfull(IVZ,k,j,ii);
      if (cs_) {
        const Real cc = gtc_.cs(m,k,j), ss = gtc_.sn(m,k,j);
        q2 += cc*q3;  q3 *= ss;
      }
      q.vx = q1;  q.vy = q2;  q.vz = q3;
      if (eos.is_ideal) { q.e = wfull(IEN,k,j,ii); }
      q.by = bcc0_(m,IBY,k,j,ii);
      q.bz = bcc0_(m,IBZ,k,j,ii);
    };
    // called AFTER both the momentum flux and the two face EMFs have been written
    auto rotflx1 = [&](const int f) {
      if (cs_) {
        const Real cc = gtc_.cs(m,k,j), ss = gtc_.sn(m,k,j);
        const Real fb = flx1(m,IM3,k,j,f)/ss;
        const Real fa = flx1(m,IM2,k,j,f) - fb*cc;
        flx1(m,IM2,k,j,f) = fa + fb*cc;
        flx1(m,IM3,k,j,f) = fb + fa*cc;
        // CT integrates E around a face, so e3 must be E on the ETA EDGE, whose length
        // dxedge.x3e it is multiplied by.  The x1 frame's third axis is (e_eta - c
        // e_xi)/s, not that edge, so E.e_eta = c*(E.e_xi) + s*(E.f2) = c*e2 + s*e3.
        e3x1_(m,k,j,f) = cc*e2x1_(m,k,j,f) + ss*e3x1_(m,k,j,f);
      }
    };
    // x2/x3 load the state ALREADY PERMUTED so the face normal comes first, exactly as
    // the Cartesian code did inline.
    auto ldx2 = [&](const int jj, const int f, MHDPrim1D &q) {
      q.d = wfull(IDN,k,jj,i);
      Real q1 = wfull(IVX,k,jj,i), q2 = wfull(IVY,k,jj,i), q3 = wfull(IVZ,k,jj,i);
      Real b1 = bcc0_(m,IBX,k,jj,i);
      Real b2 = bcc0_(m,IBY,k,jj,i), b3 = bcc0_(m,IBZ,k,jj,i);
      if (cs_) {
        const Real cc = gtxi_.cs(m,k,f), ss = gtxi_.sn(m,k,f);
        q3 += cc*q2;  q2 *= ss;
        b3 = cc*b2 + ss*b3;
      }
      q.vx = q2;  q.vy = q3;  q.vz = q1;
      if (eos.is_ideal) { q.e = wfull(IEN,k,jj,i); }
      q.by = b3;  q.bz = b1;
    };
    auto rotflx2 = [&](const int f) {
      if (cs_) {
        const Real cc = gtxi_.cs(m,k,f), ss = gtxi_.sn(m,k,f);
        const Real fb = flx2(m,IM3,k,f,i) - (cc/ss)*flx2(m,IM2,k,f,i);
        const Real fa = flx2(m,IM2,k,f,i)/ss;
        flx2(m,IM2,k,f,i) = fa + fb*cc;
        flx2(m,IM3,k,f,i) = fb + fa*cc;
      }
    };
    auto ldx3 = [&](const int kk, const int f, MHDPrim1D &q) {
      q.d = wfull(IDN,kk,j,i);
      Real q1 = wfull(IVX,kk,j,i), q2 = wfull(IVY,kk,j,i), q3 = wfull(IVZ,kk,j,i);
      if (cs_) {
        const Real cc = gtet_.cs(m,f,j), ss = gtet_.sn(m,f,j);
        q2 += cc*q3;  q3 *= ss;
      }
      q.vx = q3;  q.vy = q1;  q.vz = q2;
      if (eos.is_ideal) { q.e = wfull(IEN,kk,j,i); }
      q.by = bcc0_(m,IBX,kk,j,i);
      q.bz = bcc0_(m,IBY,kk,j,i);
    };
    auto rotflx3 = [&](const int f) {
      if (cs_) {
        const Real cc = gtet_.cs(m,f,j), ss = gtet_.sn(m,f,j);
        const Real fa = flx3(m,IM2,f,j,i) - (cc/ss)*flx3(m,IM3,f,j,i);
        const Real fb = flx3(m,IM3,f,j,i)/ss;
        flx3(m,IM2,f,j,i) = fa + fb*cc;
        flx3(m,IM3,f,j,i) = fb + fa*cc;
      }
    };

    // Check for FOFC flag
    bool fofc_flag = false;
    if (use_fofc_) { fofc_flag = fofc_(m,k,j,i); }

    // Check for GR + excision
    bool fofc_excision = false;
    if (is_gr) {
      if (use_excise_) { fofc_excision = excision_flux_(m,k,j,i); }
    }

    // Apply FOFC
    if (fofc_flag || fofc_excision) {
      // load W_{i-1} state
      MHDPrim1D wim1;
      ldx1(i-1, wim1);

      // load W_{i} state
      MHDPrim1D wi;
      ldx1(i, wi);

      // compute new 1st-order LLF flux at i-face
      {
        Real bxi = b0_.x1f(m,k,j,i);
        MHDCons1D flux;
        if (is_gr) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = LeftEdgeX(i-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          SingleStateLLF_GRMHD(wim1, wi, bxi, x1v, x2v, x3v, IVX, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRMHD(wim1, wi, bxi, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            if (fofc_hlle_) {
              SingleStateHLLE_GenMHD(wim1, wi, bxi,
                wder_(m,IDPR,k,j,i-1), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k,j,i-1), wder_(m,IDG1,k,j,i), flux);
            } else {
              SingleStateLLF_GenMHD(wim1, wi, bxi,
                wder_(m,IDPR,k,j,i-1), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k,j,i-1), wder_(m,IDG1,k,j,i), flux);
            }
          } else {
            if (fofc_hlle_) {
              SingleStateHLLE_MHD(wim1, wi, bxi, eos, flux);
            } else {
              SingleStateLLF_MHD(wim1, wi, bxi, eos, flux);
            }
          }
        }

        // store 1st-order fluxes.
        flx1(m,IDN,k,j,i) = flux.d;
        flx1(m,IM1,k,j,i) = flux.mx;
        flx1(m,IM2,k,j,i) = flux.my;
        flx1(m,IM3,k,j,i) = flux.mz;
        if (eos.is_ideal) {flx1(m,IEN,k,j,i) = flux.e;}
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in mhd_fluxes.cpp used the flux just replaced
        for (int n=nmhd_f; n<nvar_f; ++n) {
          flx1(m,n,k,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i-1)
                              : flux.d*wfull(n,k,j,i);
        }
        e3x1_(m,k,j,i) = flux.by;
        e2x1_(m,k,j,i) = flux.bz;
        rotflx1(i);
      }

      if (multi_d) {
        // load W_{j-1} state, permutting components of vectors
        MHDPrim1D wjm1;
        ldx2(j-1, j, wjm1);

        // load W_{j} state, permutting components of vectors
        MHDPrim1D wj;
        ldx2(j, j, wj);

        // compute new first-order flux at j-face
        Real bxi = b0_.x2f(m,k,j,i);
        MHDCons1D flux;
        if (is_gr) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = LeftEdgeX(j-js, nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          SingleStateLLF_GRMHD(wjm1, wj, bxi, x1v, x2v, x3v, IVY, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRMHD(wjm1, wj, bxi, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            if (fofc_hlle_) {
              SingleStateHLLE_GenMHD(wjm1, wj, bxi,
                wder_(m,IDPR,k,j-1,i), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k,j-1,i), wder_(m,IDG1,k,j,i), flux);
            } else {
              SingleStateLLF_GenMHD(wjm1, wj, bxi,
                wder_(m,IDPR,k,j-1,i), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k,j-1,i), wder_(m,IDG1,k,j,i), flux);
            }
          } else {
            if (fofc_hlle_) {
              SingleStateHLLE_MHD(wjm1, wj, bxi, eos, flux);
            } else {
              SingleStateLLF_MHD(wjm1, wj, bxi, eos, flux);
            }
          }
        }

        // store 1st-order fluxes, permutting indices.
        flx2(m,IDN,k,j,i) = flux.d;
        flx2(m,IM2,k,j,i) = flux.mx;
        flx2(m,IM3,k,j,i) = flux.my;
        flx2(m,IM1,k,j,i) = flux.mz;
        if (eos.is_ideal) {flx2(m,IEN,k,j,i) = flux.e;}
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in mhd_fluxes.cpp used the flux just replaced
        for (int n=nmhd_f; n<nvar_f; ++n) {
          flx2(m,n,k,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j-1,i)
                              : flux.d*wfull(n,k,j,i);
        }
        e1x2_(m,k,j,i) = flux.by;
        e3x2_(m,k,j,i) = flux.bz;
        rotflx2(j);
      }

      if (three_d) {
        // load W_{k-1} state, permutting components of vectors
        MHDPrim1D wkm1;
        ldx3(k-1, k, wkm1);

        // load W_{k} state, permutting components of vectors
        MHDPrim1D wk;
        ldx3(k, k, wk);

        // compute new first-order flux at k-face
        Real bxi = b0_.x3f(m,k,j,i);
        MHDCons1D flux;
        if (is_gr) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          SingleStateLLF_GRMHD(wkm1, wk, bxi, x1v, x2v, x3v, IVZ, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRMHD(wkm1, wk, bxi, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            if (fofc_hlle_) {
              SingleStateHLLE_GenMHD(wkm1, wk, bxi,
                wder_(m,IDPR,k-1,j,i), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k-1,j,i), wder_(m,IDG1,k,j,i), flux);
            } else {
              SingleStateLLF_GenMHD(wkm1, wk, bxi,
                wder_(m,IDPR,k-1,j,i), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k-1,j,i), wder_(m,IDG1,k,j,i), flux);
            }
          } else {
            if (fofc_hlle_) {
              SingleStateHLLE_MHD(wkm1, wk, bxi, eos, flux);
            } else {
              SingleStateLLF_MHD(wkm1, wk, bxi, eos, flux);
            }
          }
        }

        // store 1st-order fluxes, permutting indices.
        flx3(m,IDN,k,j,i) = flux.d;
        flx3(m,IM3,k,j,i) = flux.mx;
        flx3(m,IM1,k,j,i) = flux.my;
        flx3(m,IM2,k,j,i) = flux.mz;
        if (eos.is_ideal) {flx3(m,IEN,k,j,i) = flux.e;}
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in mhd_fluxes.cpp used the flux just replaced
        for (int n=nmhd_f; n<nvar_f; ++n) {
          flx3(m,n,k,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k-1,j,i)
                              : flux.d*wfull(n,k,j,i);
        }
        e2x3_(m,k,j,i) = flux.by;
        e1x3_(m,k,j,i) = flux.bz;
        rotflx3(k);
      }
    }
  });

  // Replace fluxes with first-order LLF fluxes at i+1,j+1,k+1 faces for any cell where
  // FOFC and/or excision is used (if GR+excising)
  par_for("FOFC-flx", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // full-state primitive accessor; see the note on the WB perturbation above
    auto wfull = [&](const int n, const int kk, const int jj, const int ii) -> Real {
      return (wbpert_ && n < nmhd_f) ? (w0_(m,n,kk,jj,ii) + w0wb_(m,n,kk,jj,ii))
                                     : w0_(m,n,kk,jj,ii);
    };

    // CUBED SPHERE.  The two panel-tangential basis vectors are not orthogonal, so a
    // Riemann solve must be handed the velocity AND the field in the sweep's locally
    // ORTHONORMAL frame, and the momentum flux it returns must be rotated back and its
    // index LOWERED -- exactly what GnomonicEquianglePrimFaceX* / GnomonicEquiangle-
    // FaceBX2 / GnomonicEquiangleFluxX* / GnomonicEquiangleEmfX1 do for the high-order
    // sweeps (coordinates/gnomonic_kernels.hpp).  Without this the first-order fallback
    // writes momentum fluxes and face EMFs in the wrong basis, which is an inconsistency
    // and not merely a loss of accuracy.  The radial component is orthogonal to both
    // angles and is never touched.  Off the cubed sphere every one of these is the
    // identity, so all other grids are bit-for-bit unchanged.
    //
    // bcc0 is stored in the orthonormal frame {rhat, e_xi, (e_eta - c e_xi)/s} (see
    // Coordinates::GnomonicEquiangleRaiseVelMHD), and each sweep's frame is {nhat, then
    // its two face-parallel axes} with the NORMAL component taken from b0.x*f: so x1 and
    // x3 need no field rotation at all, and only the x2 sweep's "by" slot does --
    // B.e_eta = c*(B.e_xi) + s*(B.f2).
    //
    // The face index `f` selects the trig: the x1 sweep uses the CELL-centred pair (it
    // depends on (k,j) only, so both x1 faces of a cell share it), the x2 sweep the
    // xi-FACE pair at the face's own j, the x3 sweep the eta-FACE pair at the face's own
    // k -- the same staggering the high-order sweeps use.
    auto ldx1 = [&](const int ii, MHDPrim1D &q) {
      q.d = wfull(IDN,k,j,ii);
      Real q1 = wfull(IVX,k,j,ii), q2 = wfull(IVY,k,j,ii), q3 = wfull(IVZ,k,j,ii);
      if (cs_) {
        const Real cc = gtc_.cs(m,k,j), ss = gtc_.sn(m,k,j);
        q2 += cc*q3;  q3 *= ss;
      }
      q.vx = q1;  q.vy = q2;  q.vz = q3;
      if (eos.is_ideal) { q.e = wfull(IEN,k,j,ii); }
      q.by = bcc0_(m,IBY,k,j,ii);
      q.bz = bcc0_(m,IBZ,k,j,ii);
    };
    // called AFTER both the momentum flux and the two face EMFs have been written
    auto rotflx1 = [&](const int f) {
      if (cs_) {
        const Real cc = gtc_.cs(m,k,j), ss = gtc_.sn(m,k,j);
        const Real fb = flx1(m,IM3,k,j,f)/ss;
        const Real fa = flx1(m,IM2,k,j,f) - fb*cc;
        flx1(m,IM2,k,j,f) = fa + fb*cc;
        flx1(m,IM3,k,j,f) = fb + fa*cc;
        // CT integrates E around a face, so e3 must be E on the ETA EDGE, whose length
        // dxedge.x3e it is multiplied by.  The x1 frame's third axis is (e_eta - c
        // e_xi)/s, not that edge, so E.e_eta = c*(E.e_xi) + s*(E.f2) = c*e2 + s*e3.
        e3x1_(m,k,j,f) = cc*e2x1_(m,k,j,f) + ss*e3x1_(m,k,j,f);
      }
    };
    // x2/x3 load the state ALREADY PERMUTED so the face normal comes first, exactly as
    // the Cartesian code did inline.
    auto ldx2 = [&](const int jj, const int f, MHDPrim1D &q) {
      q.d = wfull(IDN,k,jj,i);
      Real q1 = wfull(IVX,k,jj,i), q2 = wfull(IVY,k,jj,i), q3 = wfull(IVZ,k,jj,i);
      Real b1 = bcc0_(m,IBX,k,jj,i);
      Real b2 = bcc0_(m,IBY,k,jj,i), b3 = bcc0_(m,IBZ,k,jj,i);
      if (cs_) {
        const Real cc = gtxi_.cs(m,k,f), ss = gtxi_.sn(m,k,f);
        q3 += cc*q2;  q2 *= ss;
        b3 = cc*b2 + ss*b3;
      }
      q.vx = q2;  q.vy = q3;  q.vz = q1;
      if (eos.is_ideal) { q.e = wfull(IEN,k,jj,i); }
      q.by = b3;  q.bz = b1;
    };
    auto rotflx2 = [&](const int f) {
      if (cs_) {
        const Real cc = gtxi_.cs(m,k,f), ss = gtxi_.sn(m,k,f);
        const Real fb = flx2(m,IM3,k,f,i) - (cc/ss)*flx2(m,IM2,k,f,i);
        const Real fa = flx2(m,IM2,k,f,i)/ss;
        flx2(m,IM2,k,f,i) = fa + fb*cc;
        flx2(m,IM3,k,f,i) = fb + fa*cc;
      }
    };
    auto ldx3 = [&](const int kk, const int f, MHDPrim1D &q) {
      q.d = wfull(IDN,kk,j,i);
      Real q1 = wfull(IVX,kk,j,i), q2 = wfull(IVY,kk,j,i), q3 = wfull(IVZ,kk,j,i);
      if (cs_) {
        const Real cc = gtet_.cs(m,f,j), ss = gtet_.sn(m,f,j);
        q2 += cc*q3;  q3 *= ss;
      }
      q.vx = q3;  q.vy = q1;  q.vz = q2;
      if (eos.is_ideal) { q.e = wfull(IEN,kk,j,i); }
      q.by = bcc0_(m,IBX,kk,j,i);
      q.bz = bcc0_(m,IBY,kk,j,i);
    };
    auto rotflx3 = [&](const int f) {
      if (cs_) {
        const Real cc = gtet_.cs(m,f,j), ss = gtet_.sn(m,f,j);
        const Real fa = flx3(m,IM2,f,j,i) - (cc/ss)*flx3(m,IM3,f,j,i);
        const Real fb = flx3(m,IM3,f,j,i)/ss;
        flx3(m,IM2,f,j,i) = fa + fb*cc;
        flx3(m,IM3,f,j,i) = fb + fa*cc;
      }
    };

    // Check for FOFC flag
    bool fofc_flag = false;
    if (use_fofc_) { fofc_flag = fofc_(m,k,j,i); }

    // Check for GR + excision
    bool fofc_excision = false;
    if (is_gr) {
      if (use_excise_) { fofc_excision = excision_flux_(m,k,j,i); }
    }

    // Apply FOFC
    if (fofc_flag || fofc_excision) {
      // load W_{i} state
      MHDPrim1D wi;
      ldx1(i, wi);

      // load W_{i+1} state
      MHDPrim1D wip1;
      ldx1(i+1, wip1);

      // compute new 1st-order LLF flux at (i+1)-face
      {
        Real bxi = b0_.x1f(m,k,j,i+1);
        MHDCons1D flux;
        if (is_gr) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          SingleStateLLF_GRMHD(wi, wip1, bxi, x1v, x2v, x3v, IVX, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRMHD(wi, wip1, bxi, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            if (fofc_hlle_) {
              SingleStateHLLE_GenMHD(wi, wip1, bxi,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k,j,i+1),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k,j,i+1), flux);
            } else {
              SingleStateLLF_GenMHD(wi, wip1, bxi,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k,j,i+1),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k,j,i+1), flux);
            }
          } else {
            if (fofc_hlle_) {
              SingleStateHLLE_MHD(wi, wip1, bxi, eos, flux);
            } else {
              SingleStateLLF_MHD(wi, wip1, bxi, eos, flux);
            }
          }
        }

        // store 1st-order fluxes.
        flx1(m,IDN,k,j,i+1) = flux.d;
        flx1(m,IM1,k,j,i+1) = flux.mx;
        flx1(m,IM2,k,j,i+1) = flux.my;
        flx1(m,IM3,k,j,i+1) = flux.mz;
        if (eos.is_ideal) {flx1(m,IEN,k,j,i+1) = flux.e;}
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in mhd_fluxes.cpp used the flux just replaced
        for (int n=nmhd_f; n<nvar_f; ++n) {
          flx1(m,n,k,j,i+1) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i)
                              : flux.d*wfull(n,k,j,i+1);
        }
        e3x1_(m,k,j,i+1) = flux.by;
        e2x1_(m,k,j,i+1) = flux.bz;
        rotflx1(i+1);
      }

      if (multi_d) {
        // load W_{j} state, permutting components of vectors
        MHDPrim1D wj;
        ldx2(j, j+1, wj);

        // load W_{j+1} state, permutting components of vectors
        MHDPrim1D wjp1;
        ldx2(j+1, j+1, wjp1);

        // compute new first-order flux at (j+1)-face
        Real bxi = b0_.x2f(m,k,j+1,i);
        MHDCons1D flux;
        if (is_gr) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          SingleStateLLF_GRMHD(wj, wjp1, bxi, x1v, x2v, x3v, IVY, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRMHD(wj, wjp1, bxi, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            if (fofc_hlle_) {
              SingleStateHLLE_GenMHD(wj, wjp1, bxi,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k,j+1,i),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k,j+1,i), flux);
            } else {
              SingleStateLLF_GenMHD(wj, wjp1, bxi,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k,j+1,i),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k,j+1,i), flux);
            }
          } else {
            if (fofc_hlle_) {
              SingleStateHLLE_MHD(wj, wjp1, bxi, eos, flux);
            } else {
              SingleStateLLF_MHD(wj, wjp1, bxi, eos, flux);
            }
          }
        }

        // store 1st-order fluxes, permutting indices.
        flx2(m,IDN,k,j+1,i) = flux.d;
        flx2(m,IM2,k,j+1,i) = flux.mx;
        flx2(m,IM3,k,j+1,i) = flux.my;
        flx2(m,IM1,k,j+1,i) = flux.mz;
        if (eos.is_ideal) {flx2(m,IEN,k,j+1,i) = flux.e;}
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in mhd_fluxes.cpp used the flux just replaced
        for (int n=nmhd_f; n<nvar_f; ++n) {
          flx2(m,n,k,j+1,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i)
                              : flux.d*wfull(n,k,j+1,i);
        }
        e1x2_(m,k,j+1,i) = flux.by;
        e3x2_(m,k,j+1,i) = flux.bz;
        rotflx2(j+1);
      }

      if (three_d) {
        // load W_{k} state, permutting components of vectors
        MHDPrim1D wk;
        ldx3(k, k+1, wk);

        // load W_{k+1} state, permutting components of vectors
        MHDPrim1D wkp1;
        ldx3(k+1, k+1, wkp1);

        // compute new first-order flux at (k+1)-face
        Real bxi = b0_.x3f(m,k+1,j,i);
        MHDCons1D flux;
        if (is_gr) {
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
          SingleStateLLF_GRMHD(wk, wkp1, bxi, x1v, x2v, x3v, IVZ, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRMHD(wk, wkp1, bxi, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            if (fofc_hlle_) {
              SingleStateHLLE_GenMHD(wk, wkp1, bxi,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k+1,j,i),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k+1,j,i), flux);
            } else {
              SingleStateLLF_GenMHD(wk, wkp1, bxi,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k+1,j,i),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k+1,j,i), flux);
            }
          } else {
            if (fofc_hlle_) {
              SingleStateHLLE_MHD(wk, wkp1, bxi, eos, flux);
            } else {
              SingleStateLLF_MHD(wk, wkp1, bxi, eos, flux);
            }
          }
        }

        // store 1st-order fluxes, permutting indices.
        flx3(m,IDN,k+1,j,i) = flux.d;
        flx3(m,IM3,k+1,j,i) = flux.mx;
        flx3(m,IM1,k+1,j,i) = flux.my;
        flx3(m,IM2,k+1,j,i) = flux.mz;
        if (eos.is_ideal) {flx3(m,IEN,k+1,j,i) = flux.e;}
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in mhd_fluxes.cpp used the flux just replaced
        for (int n=nmhd_f; n<nvar_f; ++n) {
          flx3(m,n,k+1,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i)
                              : flux.d*wfull(n,k+1,j,i);
        }
        e2x3_(m,k+1,j,i) = flux.by;
        e1x3_(m,k+1,j,i) = flux.bz;
        rotflx3(k+1);
      }
    }
  });

  // reset FOFC flag (do not reset excision flag)
  if (use_fofc_) {
    Kokkos::deep_copy(fofc, false);
  }

  return;
}

} // namespace mhd
