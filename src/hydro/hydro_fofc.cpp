//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_fofc.cpp
//! \brief Implements functions for first-order flux correction (FOFC) algorithm.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/gnomonic_kernels.hpp"
#include "coordinates/gnomonic_raisevel.hpp"
#include "eos/eos.hpp"
#include "hydro/rsolvers/llf_hyd_singlestate.hpp"
#include "hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "hydro/fofc_etotgrav.hpp"
#include "bvals/bvals.hpp"

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void Hydro::FOFC
//! \brief Implements first-order flux-correction (FOFC) algorithm for Hydro.  First an
//! estimate of the updated conserved variables is made. This estimate is then used to
//! flag any cell where floors will be required during the conversion to primitives. Then
//! the fluxes on the faces of flagged cells are replaced with first-order LLF fluxes.
//! Often this is enough to prevent floors from being needed. The FOFC infrastructure is
//! also exploited for BH excision. If a cell is about the horizon, FOFC is automatically
//! triggered (without estimating updated conserved variables).

void Hydro::FOFC(Driver *pdriver, int stage) {
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

  // CURVILINEAR GRIDS.  Both spherical grids store the flux divergence in conservative
  // form with FACE AREAS and CELL VOLUMES from Coordinates -- see Hydro::RKUpdate, which
  // this mirrors exactly.  Dividing by the Cartesian mb_size.dx1/2/3 instead is silently
  // wrong on spherical polar and was a startup fatal on the cubed sphere.
  //
  // GEOMETRIC SOURCE TERMS ARE DELIBERATELY LEFT OUT, exactly as RKUpdate leaves them
  // out: on these grids they are added afterwards by Coordinates::CoordSrcTerms, from
  // HydroSrcTerms.  The published FOFC semantics are "trial update from the Riemann
  // fluxes alone", and the trial state is used only to DECIDE whether a cell needs
  // floors, never stored.  Including the geometric source would make the test state
  // differ from the one RKUpdate produces at that point in the stage and, near the
  // polar axis where the cot(theta) term is largest, would flag cells on the strength
  // of a term the flux replacement cannot influence.
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  const bool curv_ = cs_ || pmy_pack->pmesh->use_spherical_polar;
  auto &vol_  = pmy_pack->pcoord->volume;
  auto &ar1_  = pmy_pack->pcoord->area.x1f;
  auto &ar2_  = pmy_pack->pcoord->area.x2f;
  auto &ar3_  = pmy_pack->pcoord->area.x3f;
  // The gnomonic rotations the single-state solves below need.  Three (sin,cos) pairs,
  // i.e. six Views: on a GPU everything a kernel captures travels in the launch's
  // constant-memory buffer, so if this kernel ever grows, split the three sweeps into
  // three launches rather than adding a seventh -- see the note in gnomonic_kernels.hpp.
  const GnomonicTrig gtc_  = pmy_pack->pcoord->GnomonicTrigCell();
  const GnomonicTrig gtxi_ = pmy_pack->pcoord->GnomonicTrigFaceXi();
  const GnomonicTrig gtet_ = pmy_pack->pcoord->GnomonicTrigFaceEta();

  if (use_fofc) {
    Real &gam0 = pdriver->gam0[stage-1];
    Real &gam1 = pdriver->gam1[stage-1];
    Real beta_dt = (pdriver->beta[stage-1])*(pmy_pack->pmesh->dt);

    int &nhyd_ = nhydro;
    // under <hydro>/dual_energy the trial state includes the scalars
    const int ntv_ = nhydro + (dual.on ? nscalars : 0);
    auto &u0_ = u0;
    auto &u1_ = u1;
    auto &utest_ = utest;

    // Index bounds
    int il = is-1, iu = ie+1, jl = js, ju = je, kl = ks, ku = ke;
    if (multi_d) { jl = js-1, ju = je+1; }
    if (three_d) { kl = ks-1, ku = ke+1; }

    // <mesh>/cs_vertex_3cell: at a CUBE VERTEX the ng x ng corner ghost
    // block has no physical counterpart.  The trial state of a seam GHOST cell within ng
    // of a vertex corner along its row took the along-seam flux divergence from faces
    // whose reconstruction stencils reach into that block; those terms are left out of
    // the trial (the test then sees the radial and seam-normal fluxes only), and the
    // corner-block cells themselves are never flagged.  Flags of the active cells and of
    // every other ghost are unchanged.
    const bool v3_ = cs_ && pmy_pack->pmesh->cs_vertex_3cell && three_d;
    auto &vnb_ = pmy_pack->pmb->nghbr;
    auto &vpn_ = pmy_pack->pmb->mb_panel;
    const int ng_ = indcs.ng;

    // Estimate updated conserved variables and cell-centered fields
    par_for("FOFC-newu", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (curv_) {
        bool drop2 = false, drop3 = false;   // (v3_ only)
        if (v3_) {
          CsVertexFofcDrop(vnb_.d_view, vpn_.d_view, m, k, j, js, je, ks, ke, ng_,
                           drop2, drop3);
        }
        // area/volume divergence: the same form as Hydro::RKUpdate.  At a POLAR axis
        // the x2 face area is identically zero, so the pole face contributes nothing
        // here, exactly as in the real update.
        const Real dtodv = beta_dt/vol_(m,k,j,i);
        for (int n=0; n<ntv_; ++n) {
          Real divf = dtodv*(flx1(m,n,k,j,i+1)*ar1_(m,k,j,i+1) -
                             flx1(m,n,k,j,i  )*ar1_(m,k,j,i  ));
          if (multi_d && !drop2) {
            divf += dtodv*(flx2(m,n,k,j+1,i)*ar2_(m,k,j+1,i) -
                           flx2(m,n,k,j  ,i)*ar2_(m,k,j  ,i));
          }
          if (three_d && !drop3) {
            divf += dtodv*(flx3(m,n,k+1,j,i)*ar3_(m,k+1,j,i) -
                           flx3(m,n,k  ,j,i)*ar3_(m,k  ,j,i));
          }
          utest_(m,n,k,j,i) = gam0*u0_(m,n,k,j,i) + gam1*u1_(m,n,k,j,i) - divf;
        }
      } else {
      Real dtodx1 = beta_dt/size.d_view(m).dx1;
      Real dtodx2 = beta_dt/size.d_view(m).dx2;
      Real dtodx3 = beta_dt/size.d_view(m).dx3;

      // Estimate conserved variables
      for (int n=0; n<ntv_; ++n) {
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
    });

    // <hydro>/etotgrav: the trial energy still carries rho*Phi and lacks the
    // gravitational energy flux; make it e + KE of the trial state before the floor test
    // (and the cubed-sphere GnomonicRaiseVelFloors test below).  See fofc_etotgrav.hpp.
    if (use_etotgrav) {
      FofcTrialRemoveGrav(pmy_pack, utest_, uflx, phicc0, phi0, beta_dt,
                          il, iu, jl, ju, kl, ku);
    }

    // Test whether conversion to primitives requires floors
    // Note b0 and w0 passed to function, but not used/changed.
    peos->ConsToPrim(utest_, w0, true, il, iu, jl, ju, kl, ku);

    // CUBED SPHERE: the DEFERRED half of the floors.  ConsToPrim above ran under
    // EOS_Data::defer_cons_floors, i.e. it tested the pressure/temperature/energy floors
    // and the velocity ceiling against an ORTHONORMAL kinetic energy, which is NOT the
    // kinetic energy on this grid; the real decision is taken afterwards by
    // Coordinates::GnomonicEquiangleRaiseVel, and the floor-TEST pass is followed by no
    // such call.  A cell that needs one of those deferred floors was therefore invisible
    // to FOFC.  Run exactly that arithmetic on the trial state here -- the same free
    // function the real pass uses, coordinates/gnomonic_raisevel.hpp, which raises the
    // covariant momentum with the metric BEFORE testing anything -- and flag whatever it
    // would floor.  Nothing is written back: the trial state is scratch.
    // Cells the DENSITY floor already flagged are skipped; they are flagged either way,
    // and their trial density is not the density the raise would be handed.
    // (isothermal: no energy slot and no energy floor to defer)
    if (cs_ && peos->eos_data.defer_cons_floors && peos->eos_data.is_ideal &&
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
        Real v1, v2, v3, eint, pnew, g1new, temp, de, vde;
        bool ceil_used, floored, tset_used;
        // the floor-TEST pass: `vde` is discarded, since this state is a trial one and
        // nothing is written back (see EventCounters::vceil_de)
        GnomonicRaiseVelFloors(cosc_(m,k,j), eos_, gen_, false, vceil_, true,
                               1.0, gen_ ? wtemp_(m,k,j,i) : 0.0, d,
                               m1, m2, m3, etot, v1, v2, v3, eint, pnew, g1new, temp,
                               ceil_used, floored, tset_used, de, vde);
        if (ceil_used || floored) {
          fofcc_(m,k,j,i) = true;
          sumf++;
        }
      }, Kokkos::Sum<int>(nflag_));
      pmy_pack->pmesh->ecounter.nfofc += nflag_;
    }
    // <hydro>/fofc_retest (and dual_energy): RE-TEST every flagged trial cell
    // with the metric-correct internal energy (0.5 m_i v^i on the cubed sphere; the
    // ConsToPrim test above used the orthonormal m^2/2rho there) and, under dual_energy,
    // the dual selection of DualEnergy::Post; a cell no floor and no ceiling would touch
    // is unflagged.  Density-floor and non-finite trial cells stay flagged.
    if (dual.retest || dual.on) {
      auto eos_ = peos->eos_data;
      auto fofcr_ = fofc;
      auto wtempr_ = wtemp;
      auto pstr_ = dual.pst;
      const bool sel_ = dual.on && dual.valid;
      const Real eta_ = dual.eta;
      const int iden_ = dual.iden, ivol_ = dual.ivol;
      const Real vceil_ = eos_.vceil;
      auto coscr_ = cs_ ? pmy_pack->pcoord->cos_cell : DvceArray3D<Real>("dr_c",1,1,1);
      const int ni   = (iu - il + 1);
      const int nji  = (ju - jl + 1)*ni;
      const int nkji = (ku - kl + 1)*nji;
      Real nret_ = 0.0, nunf_ = 0.0;
      Kokkos::parallel_reduce("FOFC-retest",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
      KOKKOS_LAMBDA(const int &idx, Real &sr, Real &su) {
        const int m = idx/nkji;
        const int k = (idx - m*nkji)/nji + kl;
        const int j = (idx - m*nkji - (k - kl)*nji)/ni + jl;
        const int i = (idx - m*nkji - (k - kl)*nji - (j - jl)*ni) + il;
        if (!fofcr_(m,k,j,i)) return;
        sr += 1.0;
        const Real d = utest_(m,IDN,k,j,i);
        const Real m1 = utest_(m,IM1,k,j,i), m2 = utest_(m,IM2,k,j,i);
        const Real m3 = utest_(m,IM3,k,j,i), et = utest_(m,IEN,k,j,i);
        if (!Kokkos::isfinite(d) || !Kokkos::isfinite(m1) || !Kokkos::isfinite(m2) ||
            !Kokkos::isfinite(m3) || !Kokkos::isfinite(et)) return;
        if (!(d >= eos_.dfloor) || !(d > 0.0)) return;
        Real ekin;
        if (cs_) {
          const Real c = coscr_(m,k,j), det = 1.0 - c*c;
          ekin = 0.5*(m1*m1/d + m2*(m2 - c*m3)/(d*det) + m3*(m3 - c*m2)/(d*det));
        } else {
          ekin = 0.5*(m1*m1 + m2*m2 + m3*m3)/d;
        }
        if (vceil_ > 0.0 && 2.0*ekin/d > vceil_*vceil_) return;
        const Real econs = et - ekin;
        Real esel = econs;
        if (sel_ && (!(econs > 0.0) || econs < eta_*et)) {
          const Real se = utest_(m,iden_,k,j,i), sv = utest_(m,ivol_,k,j,i);
          if (se > 0.0 && sv > 0.0) {
            const Real ed = se*exp(pstr_(m,k,j,i)*(sv - 1.0)/se);
            if (ed > 0.0 && Kokkos::isfinite(ed)) esel = ed;
          }
        }
        if (!(esel > 0.0)) return;
        Real t, p, g1;
        bool tcl = false;
        eos_.TemperaturePressureGamma1(d, esel, wtempr_(m,k,j,i), t, p, g1, tcl);
        if (!(p >= eos_.pfloor) || !(t >= eos_.tfloor)) return;
        if (eos_.efloor_as_tfloor && tcl) return;
        Real etry = esel;
        if (eos_.ApplyEntropyFloor(d, 1.0/d, etry)) return;
        fofcr_(m,k,j,i) = false;
        su += 1.0;
      }, Kokkos::Sum<Real>(nret_), Kokkos::Sum<Real>(nunf_));
      pmy_pack->pmesh->ecounter.nfofc -= static_cast<int>(nunf_);
      dual.hc[fluid::DC_RETEST] += nret_;
      dual.hc[fluid::DC_UNFLAG] += nunf_;
    }
    // the cube-vertex corner-block cells are never flagged (v3_ only)
    if (v3_) {
      auto fofcv_ = fofc;
      par_for("FOFC-vtxclr", DevExeSpace(), 0, nmb-1, 0, 3, 0, ng_-1, 0, ng_-1, il, iu,
      KOKKOS_LAMBDA(const int m, const int c, const int gk, const int gj, const int i) {
        const int sj = c & 1, sk = (c >> 1) & 1;
        if (!IsCubeVertexCorner(vnb_.d_view, vpn_.d_view, m, 40 + 2*(sj + 2*sk))) return;
        const int j = sj ? (je + 1 + gj) : (js - 1 - gj);
        const int k = sk ? (ke + 1 + gk) : (ks - 1 - gk);
        if (j < jl || j > ju || k < kl || k > ku) return;
        fofcv_(m,k,j,i) = false;
      });
    }
    // Accumulate the PER-CELL FOFC flag count for the `hydro_fofc` output variable.
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
    // <hydro>/fofc_report (diagnostic, read only when named): reads the flags only
    fofc_rep.Count(pmy_pack, fofc, w0, peos->eos_data.dfloor);
  }

  auto &coord = pmy_pack->pcoord->coord_data;
  bool &is_sr = pmy_pack->pcoord->is_special_relativistic;
  bool &is_gr = pmy_pack->pcoord->is_general_relativistic;
  auto &eos = peos->eos_data;
  // derived thermodynamic variables (general EOS only; empty if ideal)
  auto &wder_ = wder;
  auto &use_fofc_ = use_fofc;
  auto &fofc_ = fofc;
  auto &use_excise = pmy_pack->pcoord->coord_data.bh_excise;
  auto &excision_flux_ = pmy_pack->pcoord->excision_flux;
  auto &w0_ = w0;
  // During Hydro::Fluxes the well-balanced PERTURBATION reconstruction has already
  // subtracted the static background from w0 (RemoveWbVar in Hydro::Fluxes), so w0 is
  // NOT a full state in this window.  The first-order fallback is a physical Riemann
  // solve and needs the full state, so the background is added back here, for the few
  // flagged cells only.  With the perturbation scheme off this is the identity on w0,
  // so every non-well-balanced run is bit-for-bit unchanged.
  auto &w0wb_ = w0wb;
  const bool wbpert_ = use_wellbalance_static_reconst_perturb;
  const int nhyd_f = nhydro;
  const int nvar_f = nhydro + nscalars;

  // Index bounds
  int il = is-1, iu = ie+1, jl = js, ju = je, kl = ks, ku = ke;
  if (multi_d) { jl = js-1, ju = je+1; }
  if (three_d) { kl = ks-1, ku = ke+1; }

  // Now replace fluxes with first-order LLF fluxes for any cell where floors needed (if
  // using FOFC) and/or for any cell about the excision (if GR+excising)
  par_for("FOFC-flx", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // full-state primitive accessor; see the note on the WB perturbation above
    auto wfull = [&](const int n, const int kk, const int jj, const int ii) -> Real {
      return (wbpert_ && n < nhyd_f) ? (w0_(m,n,kk,jj,ii) + w0wb_(m,n,kk,jj,ii))
                                     : w0_(m,n,kk,jj,ii);
    };

    // CUBED SPHERE.  The two panel-tangential basis vectors are not orthogonal, so a
    // Riemann solve must be handed the velocity in the sweep's locally ORTHONORMAL
    // frame, and the momentum flux it returns must be rotated back and its index
    // LOWERED -- exactly what GnomonicEquianglePrimFaceX*/GnomonicEquiangleFluxX* do
    // for the high-order sweeps (coordinates/gnomonic_kernels.hpp).  Without this the
    // first-order fallback writes momentum fluxes in the wrong basis, which is an
    // inconsistency and not merely a loss of accuracy.  The radial component is
    // orthogonal to both angles and is never touched.  Off the cubed sphere every one
    // of these is the identity, so all other grids are bit-for-bit unchanged.
    // The face index `f` selects the trig: the x1 sweep uses the CELL-centred pair, the
    // x2 sweep the xi-FACE pair at the face's own j, the x3 sweep the eta-FACE pair at
    // the face's own k -- the same staggering the high-order sweeps use.
    auto ldx1 = [&](const int ii, HydPrim1D &q) {
      q.d = wfull(IDN,k,j,ii);
      Real q1 = wfull(IVX,k,j,ii), q2 = wfull(IVY,k,j,ii), q3 = wfull(IVZ,k,j,ii);
      if (cs_) {
        const Real cc = gtc_.cs(m,k,j), ss = gtc_.sn(m,k,j);
        q2 += cc*q3;  q3 *= ss;
      }
      q.vx = q1;  q.vy = q2;  q.vz = q3;
      if (eos.is_ideal) { q.e = wfull(IEN,k,j,ii); }
    };
    auto rotflx1 = [&](const int f) {
      if (cs_) {
        const Real cc = gtc_.cs(m,k,j), ss = gtc_.sn(m,k,j);
        const Real fb = flx1(m,IM3,k,j,f)/ss;
        const Real fa = flx1(m,IM2,k,j,f) - fb*cc;
        flx1(m,IM2,k,j,f) = fa + fb*cc;
        flx1(m,IM3,k,j,f) = fb + fa*cc;
      }
    };
    // x2/x3 load the state ALREADY PERMUTED so the face normal comes first, exactly as
    // the Cartesian code below did inline.
    auto ldx2 = [&](const int jj, const int f, HydPrim1D &q) {
      q.d = wfull(IDN,k,jj,i);
      Real q1 = wfull(IVX,k,jj,i), q2 = wfull(IVY,k,jj,i), q3 = wfull(IVZ,k,jj,i);
      if (cs_) {
        const Real cc = gtxi_.cs(m,k,f), ss = gtxi_.sn(m,k,f);
        q3 += cc*q2;  q2 *= ss;
      }
      q.vx = q2;  q.vy = q3;  q.vz = q1;
      if (eos.is_ideal) { q.e = wfull(IEN,k,jj,i); }
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
    auto ldx3 = [&](const int kk, const int f, HydPrim1D &q) {
      q.d = wfull(IDN,kk,j,i);
      Real q1 = wfull(IVX,kk,j,i), q2 = wfull(IVY,kk,j,i), q3 = wfull(IVZ,kk,j,i);
      if (cs_) {
        const Real cc = gtet_.cs(m,f,j), ss = gtet_.sn(m,f,j);
        q2 += cc*q3;  q3 *= ss;
      }
      q.vx = q3;  q.vy = q1;  q.vz = q2;
      if (eos.is_ideal) { q.e = wfull(IEN,kk,j,i); }
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

    // ONE WRITER PER FACE.  A face between two flagged cells is replaced by both; the two
    // first-order fluxes are the same algebra on the same states but are compiled as
    // two inlined copies, which the GPU compiler may contract (FMA) differently, so they
    // can differ in the last bit -- and which store lands last is up to the scheduler:
    // GPU runs were not bitwise reproducible.  So a cell
    // replaces its UPPER face only when the cell above is not itself replacing it (it
    // is outside this kernel's range or not flagged).  The cell above writes the face
    // as its LOWER face, which is also the store that lands last in the serial loop
    // order: a host (Serial) build is bitwise unchanged.  The flags are therefore no
    // longer cleared here, but by a separate kernel below.
    auto nbr = [&](const int kk, const int jj, const int ii) -> bool {
      if (ii > iu || jj > ju || kk > ku) {return false;}
      bool f = false;
      if (use_fofc_) {f = fofc_(m,kk,jj,ii);}
      if (is_gr && use_excise) {f = f || excision_flux_(m,kk,jj,ii);}
      return f;
    };

    // Check for FOFC flag
    bool fofc_flag = false;
    if (use_fofc_) { fofc_flag = fofc_(m,k,j,i); }

    // Check for GR + excision
    bool fofc_excision = false;
    if (is_gr) {
      if (use_excise) { fofc_excision = excision_flux_(m,k,j,i); }
    }

    // Apply FOFC
    if (fofc_flag || fofc_excision) {
      // replace x1-flux at i
      // load left state
      HydPrim1D wim1;
      ldx1(i-1, wim1);

      // load right state
      HydPrim1D wi;
      ldx1(i, wi);

      // compute new 1st-order LLF flux
      HydCons1D flux;
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
        SingleStateLLF_GRHyd(wim1, wi, x1v, x2v, x3v, IVX, coord, eos, flux);
      } else if (is_sr) {
        SingleStateLLF_SRHyd(wim1, wi, eos, flux);
      } else {
        if (eos.IsGeneral()) {
          SingleStateLLF_GenHyd(wim1, wi,
              wder_(m,IDPR,k,j,i-1), wder_(m,IDPR,k,j,i),
              wder_(m,IDG1,k,j,i-1), wder_(m,IDG1,k,j,i), flux);
        } else {
          SingleStateLLF_Hyd(wim1, wi, eos, flux);
        }
      }

      // store 1st-order fluxes
      flx1(m,IDN,k,j,i) = flux.d;
      flx1(m,IM1,k,j,i) = flux.mx;
      flx1(m,IM2,k,j,i) = flux.my;
      flx1(m,IM3,k,j,i) = flux.mz;
      if (eos.is_ideal) {flx1(m,IEN,k,j,i) = flux.e;}
      rotflx1(i);
      // passive scalars ride the NEW mass flux: the higher-order scalar
      // flux formed in hydro_fluxes.cpp used the flux just replaced
      for (int n=nhyd_f; n<nvar_f; ++n) {
        flx1(m,n,k,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i-1)
                              : flux.d*wfull(n,k,j,i);
      }

      // replace x1-flux at i+1 (unless the cell at i+1 replaces it)
      if (!nbr(k, j, i+1)) {
      // load right state (left state just wi from above)
      HydPrim1D wip1;
      ldx1(i+1, wip1);

      // compute new 1st-order LLF flux
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
        SingleStateLLF_GRHyd(wi, wip1, x1v, x2v, x3v, IVX, coord, eos, flux);
      } else if (is_sr) {
        SingleStateLLF_SRHyd(wi, wip1, eos, flux);
      } else {
        if (eos.IsGeneral()) {
          SingleStateLLF_GenHyd(wi, wip1,
              wder_(m,IDPR,k,j,i), wder_(m,IDPR,k,j,i+1),
              wder_(m,IDG1,k,j,i), wder_(m,IDG1,k,j,i+1), flux);
        } else {
          SingleStateLLF_Hyd(wi, wip1, eos, flux);
        }
      }

      // store 1st-order fluxes
      flx1(m,IDN,k,j,i+1) = flux.d;
      flx1(m,IM1,k,j,i+1) = flux.mx;
      flx1(m,IM2,k,j,i+1) = flux.my;
      flx1(m,IM3,k,j,i+1) = flux.mz;
      if (eos.is_ideal) {flx1(m,IEN,k,j,i+1) = flux.e;}
      rotflx1(i+1);
      // passive scalars ride the NEW mass flux: the higher-order scalar
      // flux formed in hydro_fluxes.cpp used the flux just replaced
      for (int n=nhyd_f; n<nvar_f; ++n) {
        flx1(m,n,k,j,i+1) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i)
                              : flux.d*wfull(n,k,j,i+1);
      }
      }

      if (multi_d) {
        // replace x2-flux at j
        // load left state, permutting components of vectors
        HydPrim1D wjm1;
        ldx2(j-1, j, wjm1);

        // load right state, permutting components of vectors
        HydPrim1D wj;
        ldx2(j, j, wj);

        // compute new first-order flux
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
          SingleStateLLF_GRHyd(wjm1, wj, x1v, x2v, x3v, IVY, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRHyd(wjm1, wj, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            SingleStateLLF_GenHyd(wjm1, wj,
                wder_(m,IDPR,k,j-1,i), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k,j-1,i), wder_(m,IDG1,k,j,i), flux);
          } else {
            SingleStateLLF_Hyd(wjm1, wj, eos, flux);
          }
        }

        // store 1st-order fluxes, permutting indices
        flx2(m,IDN,k,j,i) = flux.d;
        flx2(m,IM2,k,j,i) = flux.mx;
        flx2(m,IM3,k,j,i) = flux.my;
        flx2(m,IM1,k,j,i) = flux.mz;
        if (eos.is_ideal) {flx2(m,IEN,k,j,i) = flux.e;}
        rotflx2(j);
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in hydro_fluxes.cpp used the flux just replaced
        for (int n=nhyd_f; n<nvar_f; ++n) {
          flx2(m,n,k,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j-1,i)
                                : flux.d*wfull(n,k,j,i);
        }

        // replace x2-flux at j+1 (unless the cell at j+1 replaces it)
        if (!nbr(k, j+1, i)) {
        // load left state, permutting components of vectors (just wj from above)
        // load right state, permutting components of vectors
        HydPrim1D wjp1;
        ldx2(j+1, j+1, wjp1);
        // the LEFT state at face j+1 is cell j rotated with the face-(j+1) trig, not
        // the wj built above for face j
        HydPrim1D wjl;
        ldx2(j, j+1, wjl);

        // compute new first-order flux
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
          SingleStateLLF_GRHyd(wjl, wjp1, x1v, x2v, x3v, IVY, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRHyd(wjl, wjp1, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            SingleStateLLF_GenHyd(wjl, wjp1,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k,j+1,i),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k,j+1,i), flux);
          } else {
            SingleStateLLF_Hyd(wjl, wjp1, eos, flux);
          }
        }

        // store 1st-order fluxes, permutting indices
        flx2(m,IDN,k,j+1,i) = flux.d;
        flx2(m,IM2,k,j+1,i) = flux.mx;
        flx2(m,IM3,k,j+1,i) = flux.my;
        flx2(m,IM1,k,j+1,i) = flux.mz;
        if (eos.is_ideal) {flx2(m,IEN,k,j+1,i) = flux.e;}
        rotflx2(j+1);
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in hydro_fluxes.cpp used the flux just replaced
        for (int n=nhyd_f; n<nvar_f; ++n) {
          flx2(m,n,k,j+1,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i)
                                : flux.d*wfull(n,k,j+1,i);
        }
        }
      }

      if (three_d) {
        // replace x3-flux at k
        // load left state, permutting components of vectors
        HydPrim1D wkm1;
        ldx3(k-1, k, wkm1);

        // load right state, permutting components of vectors
        HydPrim1D wk;
        ldx3(k, k, wk);

        // compute new first-order flux
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
          SingleStateLLF_GRHyd(wkm1, wk, x1v, x2v, x3v, IVZ, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRHyd(wkm1, wk, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            SingleStateLLF_GenHyd(wkm1, wk,
                wder_(m,IDPR,k-1,j,i), wder_(m,IDPR,k,j,i),
                wder_(m,IDG1,k-1,j,i), wder_(m,IDG1,k,j,i), flux);
          } else {
            SingleStateLLF_Hyd(wkm1, wk, eos, flux);
          }
        }

        // store 1st-order fluxes, permutting indices
        flx3(m,IDN,k,j,i) = flux.d;
        flx3(m,IM3,k,j,i) = flux.mx;
        flx3(m,IM1,k,j,i) = flux.my;
        flx3(m,IM2,k,j,i) = flux.mz;
        if (eos.is_ideal) {flx3(m,IEN,k,j,i) = flux.e;}
        rotflx3(k);
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in hydro_fluxes.cpp used the flux just replaced
        for (int n=nhyd_f; n<nvar_f; ++n) {
          flx3(m,n,k,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k-1,j,i)
                                : flux.d*wfull(n,k,j,i);
        }

        // replace x3-flux at k+1 (unless the cell at k+1 replaces it)
        if (!nbr(k+1, j, i)) {
        // load left state, permutting components of vectors (just wk from above)
        // load right state, permutting components of vectors
        HydPrim1D wkp1;
        ldx3(k+1, k+1, wkp1);
        // the LEFT state at face k+1 is cell k rotated with the face-(k+1) trig
        HydPrim1D wkl;
        ldx3(k, k+1, wkl);

        // compute new first-order flux
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
          SingleStateLLF_GRHyd(wkl, wkp1, x1v, x2v, x3v, IVZ, coord, eos, flux);
        } else if (is_sr) {
          SingleStateLLF_SRHyd(wkl, wkp1, eos, flux);
        } else {
          if (eos.IsGeneral()) {
            SingleStateLLF_GenHyd(wkl, wkp1,
                wder_(m,IDPR,k,j,i), wder_(m,IDPR,k+1,j,i),
                wder_(m,IDG1,k,j,i), wder_(m,IDG1,k+1,j,i), flux);
          } else {
            SingleStateLLF_Hyd(wkl, wkp1, eos, flux);
          }
        }

        // store 1st-order fluxes, permutting indices
        flx3(m,IDN,k+1,j,i) = flux.d;
        flx3(m,IM3,k+1,j,i) = flux.mx;
        flx3(m,IM1,k+1,j,i) = flux.my;
        flx3(m,IM2,k+1,j,i) = flux.mz;
        if (eos.is_ideal) {flx3(m,IEN,k+1,j,i) = flux.e;}
        rotflx3(k+1);
        // passive scalars ride the NEW mass flux: the higher-order scalar
        // flux formed in hydro_fluxes.cpp used the flux just replaced
        for (int n=nhyd_f; n<nvar_f; ++n) {
          flx3(m,n,k+1,j,i) = (flux.d >= 0.0) ? flux.d*wfull(n,k,j,i)
                                : flux.d*wfull(n,k+1,j,i);
        }
        }
      }
    }
  });
  // reset the FOFC flags (not the excision flags) once every face is written
  if (use_fofc_) {
    par_for("FOFC-reset", DevExeSpace(), 0, nmb-1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fofc_(m,k,j,i) = false;
    });
  }

  return;
}

namespace {
//! the outflow limiter's face correction: the mass flux fd at a face whose donor has the
//! factor t < 1 is scaled by t, and the momentum (mr: the donor's m/rho, covariant on the
//! cubed sphere), energy (h: (e + p + m.v/2)/rho, + the face potential under etotgrav)
//! and scalars (sr) that mass carries leave with it; the pressure part of the flux is
//! left alone (well-balanced)
KOKKOS_INLINE_FUNCTION
void DflFace(const DvceArray5D<Real> &f, const int m, const int k, const int j,
             const int i, const Real t, const Real m1, const Real m2, const Real m3,
             const Real h, const DvceArray5D<Real> &w, const int nvar, const int kd,
             const int jd, const int id) {
  const Real df = (1.0 - t)*f(m,IDN,k,j,i);
  f(m,IDN,k,j,i) -= df;
  f(m,IM1,k,j,i) -= df*m1;
  f(m,IM2,k,j,i) -= df*m2;
  f(m,IM3,k,j,i) -= df*m3;
  f(m,IEN,k,j,i) -= df*h;
  for (int n = IEN + 1; n < nvar; ++n) f(m,n,k,j,i) -= df*w(m,n,kd,jd,id);
}
//! the donor's specific values: 0 th, 1-3 m/rho, 4 (e + p + m.v/2)/rho
KOKKOS_INLINE_FUNCTION
Real DflSpec(const int q, const DvceArray5D<Real> &u, const DvceArray5D<Real> &w,
             const DvceArray5D<Real> &wd, const DvceArray5D<Real> &th, const bool gen,
             const Real gm1, const int m, const int k, const int j, const int i) {
  if (q == 0) return th(m,0,k,j,i);
  const Real d = u(m,IDN,k,j,i);
  if (q < 4) return u(m,q,k,j,i)/d;
  const Real e = w(m,IEN,k,j,i), p = gen ? wd(m,IDPR,k,j,i) : gm1*e;
  return (e + p + 0.5*(u(m,IM1,k,j,i)*w(m,IVX,k,j,i) + u(m,IM2,k,j,i)*w(m,IVY,k,j,i) +
                       u(m,IM3,k,j,i)*w(m,IVZ,k,j,i)))/d;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void Hydro::DfloorFluxLimit
//! \brief <hydro>/dfloor_flux_limit: the stage's outflows of every cell are scaled so
//! that the RK update cannot take rho below dfloor (inflows ignored, so the bound holds
//! for any inflow; convex RK combinations keep it).  Per cell,
//! th = min(1, (gam0 rho0 + gam1 rho1 - dfloor)/(beta dt sum_out A F/V)); the factors
//! go to the ghosts by one exchange, and every face takes its donor cell's factor and
//! specific values, so the two blocks of a block face hold the same flux and mass,
//! momentum and energy stay exact (no mass is created).  On the cubed sphere the donor
//! across a panel seam (resampled ghosts) comes from the seam record pdfl_fs (seam_raw,
//! the m/rho pair re-expressed on the receiver's basis).  Called by RKUpdate, on the
//! final fluxes (the seam fluxes already averaged).

void Hydro::DfloorFluxLimit(Driver *pdrive, int stage) {
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, ng = indcs.ng;
  const int nmb = pmy_pack->nmb_thispack;
  const bool md = pm->multi_d, td = pm->three_d, cs = pm->use_cubed_sphere;
  const int n1 = indcs.nx1 + 2*ng, n2 = md ? indcs.nx2 + 2*ng : 1;
  const int n3 = td ? indcs.nx3 + 2*ng : 1;
  if (dfl_th.extent_int(0) != nmb) {
    if (pm->multilevel || nhydro != 5) {
      std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << "<hydro>/"
                << "dfloor_flux_limit: uniform mesh, adiabatic hydro only" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Kokkos::realloc(dfl_th, nmb, 1, n3, n2, n1);
    Kokkos::realloc(dfl_thc, nmb, 1, 1, 1, 1);
    pbval_dfl->InitializeBuffers(1);
    if (cs) {   // the panel-seam faces (bits: x2 lo / hi, x3 lo / hi)
      pbval_dfs->InitializeBuffers(5);
      pdfl_fs = new DvceFaceFld5D<Real>("dfl_fs", nmb, 5, n3, n2, n1);
      Kokkos::realloc(pdfl_fs->x1f, nmb, 5, 1, 1, 1);
      Kokkos::realloc(dfl_sm, nmb);
      auto hm = Kokkos::create_mirror_view(dfl_sm);
      const int nb_[4] = {8, 12, 24, 28};
      for (int m = 0; m < nmb; ++m) {
        hm(m) = 0;
        for (int q = 0; q < 4; ++q) {
          const auto nb = pmy_pack->pmb->nghbr.h_view(m, nb_[q]);
          if (nb.gid >= 0 && nb.panel != pmy_pack->pmb->mb_panel.h_view(m)) {
            hm(m) |= (1 << q);
          }
        }
      }
      Kokkos::deep_copy(dfl_sm, hm);
    }
  }
  const Real g0 = pdrive->gam0[stage-1], g1 = pdrive->gam1[stage-1];
  const Real bdt = pdrive->beta[stage-1]*pm->dt;
  const Real dfl = peos->eos_data.dfloor;
  const bool curv = cs || pm->use_spherical_polar;
  auto &vol = pmy_pack->pcoord->volume;
  auto &a1 = pmy_pack->pcoord->area.x1f;
  auto &a2 = pmy_pack->pcoord->area.x2f;
  auto &a3 = pmy_pack->pcoord->area.x3f;
  auto &size = pmy_pack->pmb->mb_size;
  auto f1 = uflx.x1f, f2 = uflx.x2f, f3 = uflx.x3f;
  auto u0_ = u0, u1_ = u1, th = dfl_th;
  Kokkos::deep_copy(th, 1.0);
  Real nlim = 0.0;
  Kokkos::parallel_reduce("dfl_th", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(),
                          {0, ks, js, is}, {nmb, ke + 1, je + 1, ie + 1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &sn) {
    const auto &sz = size.d_view(m);
    const Real w1 = curv ? 1.0 : 1.0/sz.dx1, w2 = curv ? 1.0 : 1.0/sz.dx2;
    const Real w3 = curv ? 1.0 : 1.0/sz.dx3;
    Real q = (fmax(f1(m,IDN,k,j,i+1), 0.0)*(curv ? a1(m,k,j,i+1) : w1) +
              fmax(-f1(m,IDN,k,j,i), 0.0)*(curv ? a1(m,k,j,i) : w1));
    if (md) q += (fmax(f2(m,IDN,k,j+1,i), 0.0)*(curv ? a2(m,k,j+1,i) : w2) +
                  fmax(-f2(m,IDN,k,j,i), 0.0)*(curv ? a2(m,k,j,i) : w2));
    if (td) q += (fmax(f3(m,IDN,k+1,j,i), 0.0)*(curv ? a3(m,k+1,j,i) : w3) +
                  fmax(-f3(m,IDN,k,j,i), 0.0)*(curv ? a3(m,k,j,i) : w3));
    q *= bdt/(curv ? vol(m,k,j,i) : 1.0);
    const Real room = g0*u0_(m,IDN,k,j,i) + g1*u1_(m,IDN,k,j,i) - dfl;
    if (q > room) {
      th(m,0,k,j,i) = (room > 0.0) ? room/q : 0.0;
      sn += 1.0;
    }
  }, Kokkos::Sum<Real>(nlim));
  dfl_n += nlim;
  while (pbval_dfl->InitRecv(1) == TaskStatus::incomplete) {}
  while (pbval_dfl->PackAndSendCC(dfl_th, dfl_thc) == TaskStatus::incomplete) {}
  while (pbval_dfl->RecvAndUnpackCC(dfl_th, dfl_thc) == TaskStatus::incomplete) {}
  while (pbval_dfl->ClearSend() == TaskStatus::incomplete) {}
  while (pbval_dfl->ClearRecv() == TaskStatus::incomplete) {}
  auto w0_ = w0, wd = wder;
  const int nvar = nhydro + nscalars;
  const bool gen = peos->eos_data.IsGeneral(), eg = use_etotgrav;
  const Real gm1 = peos->eos_data.gamma - 1.0;
  auto p1 = phi0.x1f, p2 = phi0.x2f, p3 = phi0.x3f;
  DvceArray5D<Real> r2 = f2, r3 = f3;   // (the seam record's faces on the cubed sphere)
  DvceArray1D<int> sm;
  if (cs) {
    r2 = pdfl_fs->x2f;
    r3 = pdfl_fs->x3f;
    sm = dfl_sm;
    par_for("dfl_rec2", DevExeSpace(), 0, nmb - 1, ks, ke, 0, 1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int sd, const int i) {
      const int jf = sd ? je + 1 : js, jc = sd ? je : js;
      for (int q = 0; q < 5; ++q) {
        r2(m,q,k,jf,i) = DflSpec(q, u0_, w0_, wd, th, gen, gm1, m, k, jc, i);
      }
    });
    par_for("dfl_rec3", DevExeSpace(), 0, nmb - 1, 0, 1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int sd, const int j, const int i) {
      const int kf = sd ? ke + 1 : ks, kc = sd ? ke : ks;
      for (int q = 0; q < 5; ++q) {
        r3(m,q,kf,j,i) = DflSpec(q, u0_, w0_, wd, th, gen, gm1, m, kc, j, i);
      }
    });
    pbval_dfs->InitFluxSeamRecv(5);
    pbval_dfs->PackAndSendFluxSeamCC(*pdfl_fs);
    while (pbval_dfs->RecvAndUnpackFluxSeamCC(*pdfl_fs) != TaskStatus::complete) {}
    pbval_dfs->ClearFluxSend();
    pbval_dfs->ClearFluxRecv();
  }
  // dr = 1, 2, 3: the face's direction; the donor is (kd, jd, id); a seam face's outer
  // donor (the ghost side) comes from the record
  auto face = [=] KOKKOS_FUNCTION (const DvceArray5D<Real> &f, const DvceArray5D<Real> &r,
                                   const int dr, const int m, const int k, const int j,
                                   const int i, const Real phif) {
    const bool up = (f(m,IDN,k,j,i) > 0.0);   // the donor is the lower cell
    const int kd = (dr == 3 && up) ? k - 1 : k, jd = (dr == 2 && up) ? j - 1 : j;
    const int id = (dr == 1 && up) ? i - 1 : i;
    const int bit = (dr == 2) ? ((j == js) ? 1 : ((j == je + 1) ? 2 : 0))
                              : ((dr == 3) ? ((k == ks) ? 4 : ((k == ke + 1) ? 8 : 0)) : 0);
    const bool srec = cs && bit != 0 && (sm(m) & bit) &&
                      (((bit & 5) && up) || ((bit & 10) && !up));
    Real v[5];
    for (int q = 0; q < 5; ++q) {
      v[q] = srec ? r(m,q,k,j,i) : DflSpec(q, u0_, w0_, wd, th, gen, gm1, m, kd, jd, id);
    }
    if (!(v[0] < 1.0)) return;
    DflFace(f, m, k, j, i, v[0], v[1], v[2], v[3], v[4] + phif, w0_, nvar, kd, jd, id);
  };
  par_for("dfl_f1", DevExeSpace(), 0, nmb - 1, ks, ke, js, je, is, ie + 1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    face(f1, f1, 1, m, k, j, i, eg ? p1(m,k,j,i) : 0.0);
  });
  if (md) {
    par_for("dfl_f2", DevExeSpace(), 0, nmb - 1, ks, ke, js, je + 1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      face(f2, r2, 2, m, k, j, i, eg ? p2(m,k,j,i) : 0.0);
    });
  }
  if (td) {
    par_for("dfl_f3", DevExeSpace(), 0, nmb - 1, ks, ke + 1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      face(f3, r3, 3, m, k, j, i, eg ? p3(m,k,j,i) : 0.0);
    });
  }
}

} // namespace hydro
