//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_fluxes.cpp
//! \brief Calculate 3D fluxes for hydro

#include <iostream>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro.hpp"
#include "eos/eos.hpp"
#include "reconstruct/dc.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
#include "hydro/rsolvers/advect_hyd.hpp"
#include "hydro/rsolvers/llf_hyd.hpp"
#include "hydro/rsolvers/hlle_hyd.hpp"
#include "hydro/rsolvers/hllc_hyd.hpp"
#include "hydro/rsolvers/lhllc_hyd.hpp"
#include "hydro/rsolvers/hllclm_hyd.hpp"
#include "hydro/rsolvers/ausmpup_hyd.hpp"
#include "hydro/rsolvers/roe_hyd.hpp"
#include "hydro/rsolvers/llf_srhyd.hpp"
#include "hydro/rsolvers/hlle_srhyd.hpp"
#include "hydro/rsolvers/hllc_srhyd.hpp"
#include "hydro/rsolvers/llf_grhyd.hpp"
#include "hydro/rsolvers/hlle_grhyd.hpp"

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void Hydro::CalculateFluxes
//! \brief Calls reconstruction and Riemann solver functions to compute hydro fluxes
//! Note this function is templated over RS for better performance on GPUs.

template <Hydro_RSolver rsolver_method_>
void Hydro::CalculateFluxes(Driver *pdriver, int stage) {
  RegionIndcs &indcs_ = pmy_pack->pmesh->mb_indcs;
  int is = indcs_.is, ie = indcs_.ie;
  int js = indcs_.js, je = indcs_.je;
  int ks = indcs_.ks, ke = indcs_.ke;
  int ncells1 = indcs_.nx1 + 2*(indcs_.ng);

  int &nhyd_  = nhydro;
  int nvars = nhydro + nscalars;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const auto recon_method_ = recon_method;
  bool extrema = false;
  if (recon_method == ReconstructionMethod::ppmx) {
    extrema = true;
  }

  // --- <problem>/nan_report: the HLLC in-kernel capture (see hllc_hyd.hpp).  Zeroed once
  // per Fluxes call; the first face with a non-finite energy flux wins the record.  When
  // the switch is off these are 1-element dummies and the solver never touches them.
  const bool nanrep_r_ = nan_report;
  if (nanrep_r_ && wbrec_nanrep_cnt == nullptr) {
    wbrec_nanrep_cnt = new DvceArray1D<int>("wbrec_nanrep_cnt", 1);
    wbrec_nanrep_rec = new DvceArray1D<Real>("wbrec_nanrep_rec", 48);
  }
  auto wbcnt_ = nanrep_r_ ? *wbrec_nanrep_cnt : DvceArray1D<int>("d", 1);
  auto wbrec_ = nanrep_r_ ? *wbrec_nanrep_rec : DvceArray1D<Real>("d", 1);
  if (nanrep_r_) {
    Kokkos::deep_copy(wbcnt_, 0);
    Kokkos::deep_copy(wbrec_, 0.0);
  }
  const Real nanrep_stage_ = static_cast<Real>(stage);
  auto wtemp_nr_ = wtemp;
  if (nanrep_r_ && rsolv_nanrep_cnt == nullptr) {
    rsolv_nanrep_cnt = new DvceArray1D<int>("rsolv_nanrep_cnt", 1);
    rsolv_nanrep_rec = new DvceArray1D<Real>("rsolv_nanrep_rec", 40);
  }
  auto rscnt_ = nanrep_r_ ? *rsolv_nanrep_cnt : DvceArray1D<int>("d", 1);
  auto rsrec_ = nanrep_r_ ? *rsolv_nanrep_rec : DvceArray1D<Real>("d", 1);
  if (nanrep_r_) {
    Kokkos::deep_copy(rscnt_, 0);
    Kokkos::deep_copy(rsrec_, 0.0);
  }

  auto &eos_ = peos->eos_data;
  auto &size_ = pmy_pack->pmb->mb_size;
  auto &coord_ = pmy_pack->pcoord->coord_data;
  auto &w0_ = w0;

  // For a general EOS, pressure and Gamma_1 were evaluated once per cell in ConsToPrim
  // and are reconstructed to interfaces alongside the primitives, so that the Riemann
  // solvers consume pressure directly instead of calling the (expensive) EOS. For an
  // ideal gas nder is 0: no extra scratch is allocated and no extra reconstruction is
  // performed, leaving the original code path untouched.
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
  auto phi0_x1f = phi0.x1f;
  auto phi0_x2f = phi0.x2f;
  auto phi0_x3f = phi0.x3f;
    
  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  const auto wb_option_ = wb_option;
  const bool use_wb_rho_ = use_wb_rho;
  const bool use_wb_x1_ = use_wb_x1;
  const bool use_wb_x2_ = use_wb_x2;
  const bool use_wb_x3_ = use_wb_x3;
  const bool use_wellbalance_dynamic_ = use_wellbalance_dynamic;
  const Real wb_rmax_ = wb_rmax;
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
  auto &mb_bcs_pq = pmy_pack->pmb->mb_bcs;
  const bool pquad_ = pmy_pack->pmesh->use_polar_quadratic_recon;
  const int px3_ = pmy_pack->pmesh->polar_x3_shift;

  //--------------------------------------------------------------------------------------
  // i-direction

  size_t scr_size = (ScrArray2D<Real>::shmem_size(nvars, ncells1) +
                     ScrArray2D<Real>::shmem_size(nder, ncells1)) * 2;
  int scr_level = 0;
  auto &flx1_ = uflx.x1f;
    
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;

  const bool sp_favg_ = pmy_pack->pcoord->sp_cart_all_momentum &&
                        pmy_pack->pcoord->sp_face_avg;
  // set the loop limits for 1D/2D/3D problems
  int il = is, iu = ie+1, jl = js, ju = je, kl = ks, ku = ke;
  if (use_fofc) {
    il = is-1, iu = ie+2;
    if (pmy_pack->pmesh->one_d) {
      // 1D: j,k have no ghost zones to widen into (the arrays are one cell wide)
      jl = js, ju = je, kl = ks, ku = ke;
    } else if (pmy_pack->pmesh->two_d) {
      jl = js-1, ju = je+1, kl = ks, ku = ke;
    } else {
      jl = js-1, ju = je+1, kl = ks-1, ku = ke+1;
    }
  }

  // the well-balanced background of every cell in this sweep, walked once

  bool wbreb_h = false;
  if (use_wellbalance_dynamic && use_wb_x1) {
    const int ncyc = pmy_pack->pmesh->ncycle;
    if (wb_cache_every <= 0 || !wb_cache_built ||
        (stage == 1 && (ncyc % wb_cache_every) == 0)) {
      BuildWBCache(jl, ju, kl, ku);
      wb_cache_built = true;
      wbreb_h = true;
    }
  }


  const Real wbreb_ = wbreb_h ? 1.0 : 0.0;

  par_for_outer("hflux_x1",DevExeSpace(), scr_size, scr_level, 0, nmb1, kl, ku, jl, ju,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    auto wbcnt = wbcnt_;
    auto wbrec = wbrec_;
    const Real nanrep_stage = nanrep_stage_;
    const Real wbrebuilt = wbreb_;
    auto wbq0 = wbq0_;
    auto w0 = w0_;
    auto wder = wder_;
    auto wtemp = wtemp_nr_;
    auto phicc = phicc0_;
    auto phif = phi0_x1f;
    const bool nanrep_r = nanrep_r_;
    ScrArray2D<Real> wl(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> wr(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> dl(member.team_scratch(scr_level), nder, ncells1);
    ScrArray2D<Real> dr(member.team_scratch(scr_level), nder, ncells1);

    if (use_spherical_polar || str_r1_)
      {
        GridPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
                              use_wellbalance_dynamic_,
                              use_wb_x1_, wb_rmax_,
                              m, k, j, il-1, iu, w0_, x1v_, x1f_, phicc0_,
                              phi0_x1f, wbq0_, wl, wr);
      } else {

    if (use_wellbalance_dynamic_ && use_wb_x1_)
    {
      WbLocalPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
          m, k, j, il-1, iu, w0_, phicc0_, phi0_x1f, wbq0_, wl, wr);
    } else {
          
    // Reconstruct qR[i] and qL[i+1]
    switch (recon_method_) {
      case ReconstructionMethod::dc:
        DonorCellX1(member, m, k, j, il-1, iu, w0_, wl, wr);
        break;
      case ReconstructionMethod::plm:
        PiecewiseLinearX1(member, m, k, j, il-1, iu, w0_, wl, wr);
        break;
      case ReconstructionMethod::ppm4:
      case ReconstructionMethod::ppmx:
        PiecewiseParabolicX1(member,eos_,extrema,true, m, k, j, il-1, iu, w0_, wl, wr);
        break;
      case ReconstructionMethod::wenoz:
        WENOZX1(member, eos_, true, m, k, j, il-1, iu, w0_, wl, wr);
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
      if (use_spherical_polar || str_r1_) {
        GridPiecewiseLinearDerX1(member, eos_, wb_option_, use_wellbalance_dynamic_,
                                 use_wb_x1_, wb_rmax_,
                                 m, k, j, il-1, iu, w0_, wder_,
                                 x1v_, x1f_, phicc0_, phi0_x1f, wbq0_, dl, dr);
      } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX1(member, m, k, j, il-1, iu,
                                     pwb_, pfacewb_x1f,
                                    wder_, dl, dr);
      } else if (use_wellbalance_dynamic_ && use_wb_x1_) {
        WbPiecewiseLinearDerX1(member, eos_, wb_option_,
            m, k, j, il-1, iu, w0_, wder_,
                               phicc0_, phi0_x1f, wbq0_, dl, dr);
      } else {
      switch (recon_method_) {
        case ReconstructionMethod::dc:
          DonorCellX1(member, m, k, j, il-1, iu, wder_, dl, dr);
          break;
        case ReconstructionMethod::plm:
          PiecewiseLinearX1(member, m, k, j, il-1, iu, wder_, dl, dr);
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
      // --- <problem>/nan_report: look at the reconstruction BEFORE the fmax mask below.
      // fmax(NaN, pfloor) returns pfloor, so a NaN interface pressure is invisible after
      // it; the interface ENERGY has no such mask and reaches the solver as NaN.
      if (nanrep_r) {
        par_for_inner(member, il, iu, [&](const int i) {
          if (isfinite(dl(IDPR,i)) && isfinite(dr(IDPR,i)) &&
              isfinite(wl(IEN,i)) && isfinite(wr(IEN,i))) return;
          if (Kokkos::atomic_fetch_add(&wbcnt(0), 1) != 0) return;
          Real bd[5], be[5], bp[5];
          WBReadCache(wbq0, WBVar::wb_dens, m, k, j, i,
                      bd[0], bd[1], bd[2], bd[3], bd[4]);
          WBReadCache(wbq0, WBVar::wb_eint, m, k, j, i,
                      be[0], be[1], be[2], be[3], be[4]);
          WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i,
                      bp[0], bp[1], bp[2], bp[3], bp[4]);
          wbrec(0) = static_cast<Real>(m);
          wbrec(1) = static_cast<Real>(k);
          wbrec(2) = static_cast<Real>(j);
          wbrec(3) = static_cast<Real>(i);
          wbrec(4) = nanrep_stage;
          wbrec(5) = wbrebuilt;
          for (int q = 0; q < 3; ++q) {
            wbrec(6+q)  = w0(m,IDN,k,j,i-1+q);
            wbrec(9+q)  = w0(m,IEN,k,j,i-1+q);
            wbrec(12+q) = wtemp(m,k,j,i-1+q);
            wbrec(15+q) = wder(m,IDPR,k,j,i-1+q);
            wbrec(18+q) = phicc(m,k,j,i-1+q);
          }
          for (int q = 0; q < 5; ++q) {
            wbrec(21+q) = bd[q];
            wbrec(26+q) = be[q];
            wbrec(31+q) = bp[q];
          }
          wbrec(36) = dl(IDPR,i);
          wbrec(37) = dr(IDPR,i);
          wbrec(38) = wl(IEN,i);
          wbrec(39) = wr(IEN,i);
          wbrec(40) = dl(IDG1,i);
          wbrec(41) = dr(IDG1,i);
          wbrec(42) = wl(IDN,i);
          wbrec(43) = wr(IDN,i);
          wbrec(44) = phif(m,k,j,i);
        });
        member.team_barrier();
      }
      par_for_inner(member, il, iu, [&](const int i) {
        // fmax(NaN, pfloor) returns pfloor, so the floor below does not just floor an
        // undershoot -- it MASKS a non-finite reconstructed pressure and leaves the
        // interface energy wl/wr(IEN), which has no floor of its own, NaN; the HLLC
        // general-EOS branch then poisons flx(IEN) through 0*NaN while every momentum
        // flux stays finite (the signature of the wb_recon_x1 death).  Repair the pair
        // consistently instead: floored pressure AND the internal energy that belongs
        // to it, pfloor/(Gamma_1 - 1) with this interface's reconstructed Gamma_1
        // (gamma-law fallback if that is unusable too).  The nan_report block above has
        // already counted and recorded the event, so nothing is hidden any more.
        if (!isfinite(dl(IDPR,i)) || !isfinite(wl(IEN,i))) {
          Real g1 = (isfinite(dl(IDG1,i)) && dl(IDG1,i) > 1.0) ? dl(IDG1,i) : eos_.gamma;
          dl(IDPR,i) = eos_.pfloor;
          wl(IEN,i) = eos_.pfloor/(g1 - 1.0);
        }
        if (!isfinite(dr(IDPR,i)) || !isfinite(wr(IEN,i))) {
          Real g1 = (isfinite(dr(IDG1,i)) && dr(IDG1,i) > 1.0) ? dr(IDG1,i) : eos_.gamma;
          dr(IDPR,i) = eos_.pfloor;
          wr(IEN,i) = eos_.pfloor/(g1 - 1.0);
        }
        dl(IDPR,i) = fmax(dl(IDPR,i), eos_.pfloor);
        dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
      });
    }

      if (use_wellbalance_static_reconst_perturb_) {
        AddWbPrimFaceX1(member,m,k,j,il-1,iu,w0facewb_x1f,wl,wr);
      }

      if (use_cubed_sphere) {
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
    if (use_spherical_polar || str_r1_) {
      member.team_barrier();
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        if (mb_bcs_pq.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect) {
          for (int n=0; n<nvars; ++n) { wl(n,is) = wr(n,is); }
          wl(IVX,is) = -wr(IVX,is);
        }
        if (mb_bcs_pq.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::reflect) {
          for (int n=0; n<nvars; ++n) { wr(n,ie+1) = wl(n,ie+1); }
          wr(IVX,ie+1) = -wl(IVX,ie+1);
        }
      });
    }
    // Sync all threads in the team so that scratch memory is consistent
    member.team_barrier();

    // compute fluxes over [is,ie+1]
    // NOTE(@pdmullen): Capture variables prior to if constexpr.  Required for cuda 11.6+.
    auto eos = eos_;
    auto rscnt = rscnt_;
    auto rsrec = rsrec_;
    auto indcs = indcs_;
    auto size = size_;
    auto coord = coord_;
    auto flx1 = flx1_;
    if constexpr (rsolver_method_ == Hydro_RSolver::advect) {
      Advect(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::llf) {
      LLF(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, dl, dr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle) {
      HLLE(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, dl, dr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::hllc) {
      HLLC(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,dl,dr,flx1,
           nanrep_r,rscnt,rsrec);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::lhllc) {
      LHLLC(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,dl,dr,flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::hllclm) {
      HLLCLM(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,dl,dr,flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::ausmpup) {
      AUSMPUP(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,dl,dr,flx1);
      // hybrid: HLLC on the physical x1 wall faces (zero wall mass flux by construction)
      if (eos.ausm_wall_hllc) {
        const BoundaryFlag bi = mb_bcs_pq.d_view(m,BoundaryFace::inner_x1);
        const BoundaryFlag bo = mb_bcs_pq.d_view(m,BoundaryFace::outer_x1);
        if (bi == BoundaryFlag::reflect || bi == BoundaryFlag::user) {
          member.team_barrier();
          HLLC(member,eos,indcs,size,coord,m,k,j,is,is,IVX,wl,wr,dl,dr,flx1,
               nanrep_r,rscnt,rsrec);
        }
        if (bo == BoundaryFlag::reflect || bo == BoundaryFlag::user) {
          member.team_barrier();
          HLLC(member,eos,indcs,size,coord,m,k,j,ie+1,ie+1,IVX,wl,wr,dl,dr,flx1,
               nanrep_r,rscnt,rsrec);
        }
      }
    } else if constexpr (rsolver_method_ == Hydro_RSolver::roe) {
      Roe(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::llf_sr) {
      LLF_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle_sr) {
      HLLE_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::hllc_sr) {
      HLLC_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::llf_gr) {
      LLF_GR(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle_gr) {
      HLLE_GR(member, eos, indcs, size, coord, m, k, j, il, iu, IVX, wl, wr, flx1);
    }
      
      if (use_cubed_sphere) {
          GnomonicEquiangleFluxX1(gtrig_cell,member,m,k,j,il,iu,flx1);
      }
    member.team_barrier();

    // calculate fluxes of scalars (if any)
    if (nvars > nhyd_) {
      for (int n=nhyd_; n<nvars; ++n) {
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
                ScrArray2D<Real>::shmem_size(nder, ncells1)) * 3;
    auto &flx2_ = uflx.x2f;
      
    auto &x2v_ = pmy_pack->pcoord->x2v;
    auto &x2f_ = pmy_pack->pcoord->xx2f;

    // set the loop limits for 1D/2D/3D problems
    il = is, iu = ie, jl = js-1, ju = je+1, kl = ks, ku = ke;
    if (use_fofc) {
      jl = js-2, ju = je+2;
      if (pmy_pack->pmesh->two_d) {
        il = is-1, iu = ie+1, kl = ks, ku = ke;
      } else {
        il = is-1, iu = ie+1, kl = ks-1, ku = ke+1;
      }
    } else if (sp_favg_) {
      // the face-average stencils read the x2 fluxes at i +- 1 and k +- 1
      il = is-1, iu = ie+1;
      if (!pmy_pack->pmesh->two_d) { kl = ks-1, ku = ke+1; }
    }

    par_for_outer("hflux_x2",DevExeSpace(), scr_size, scr_level, 0, nmb1, kl, ku,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k) {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> dscr1(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr2(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr3(member.team_scratch(scr_level), nder, ncells1);

      for (int j=jl; j<=ju; ++j) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_jp1 = scr2;
        auto wr     = scr3;
        auto dl     = dscr1;
        auto dl_jp1 = dscr2;
        auto dr     = dscr3;
        if ((j%2) == 0) {
          wl     = scr2;
          wl_jp1 = scr1;
          dl     = dscr2;
          dl_jp1 = dscr1;
        }

        const bool pq = pquad_ &&
          ((j == js &&
            mb_bcs_pq.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar) ||
           (j == je &&
            mb_bcs_pq.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar));
        if (use_spherical_polar) {
            GridPiecewiseLinearX2(member, m, k, j, is-1, ie+1, w0_, x2v_, x2f_,
                                  wl_jp1, wr, pq);
        } else {
            
        if (use_wellbalance_dynamic_ && use_wb_x2_)
        {
          WbLocalPiecewiseLinearX2(member, eos_, wb_option_, use_wb_rho_,
              m, k, j, il, iu, w0_, phicc0_, phi0_x2f, wl_jp1, wr);
        } else {

        // Reconstruct qR[j] and qL[j+1]
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX2(member, m, k, j, il, iu, w0_, wl_jp1, wr);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX2(member, m, k, j, il, iu, w0_, wl_jp1, wr);
            break;
          case ReconstructionMethod::ppm4:
          case ReconstructionMethod::ppmx:
            PiecewiseParabolicX2(member,eos_,extrema,true,m,k,j,il,iu, w0_, wl_jp1, wr);
            break;
          case ReconstructionMethod::wenoz:
            WENOZX2(member, eos_, true, m, k, j, il, iu, w0_, wl_jp1, wr);
            break;
          default:
            break;
        }
            
        }
        }
          
        // Reconstruct the derived thermodynamic variables (general EOS only)
        if (nder > 0) {
          if (use_spherical_polar) {
            GridPiecewiseLinearX2(member, m, k, j, il, iu, wder_, x2v_, x2f_, dl_jp1, dr,
                                  pq);
          } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX2(member, m, k, j, il, iu,
                                     pwb_, pfacewb_x2f,
                                    wder_, dl_jp1, dr);
      } else if (use_wellbalance_dynamic_ && use_wb_x2_) {
            WbPiecewiseLinearDerX2(member, eos_, wb_option_,
                m, k, j, il, iu, w0_, wder_,
                                   phicc0_, phi0_x2f, dl_jp1, dr);
          } else {
          switch (recon_method_) {
            case ReconstructionMethod::dc:
              DonorCellX2(member, m, k, j, il, iu, wder_, dl_jp1, dr);
              break;
            case ReconstructionMethod::plm:
              PiecewiseLinearX2(member, m, k, j, il, iu, wder_, dl_jp1, dr);
              break;
            case ReconstructionMethod::ppm4:
            case ReconstructionMethod::ppmx:
              PiecewiseParabolicX2(member,eos_,extrema,false,m,k,j,il,iu,wder_,dl_jp1,dr);
              break;
            case ReconstructionMethod::wenoz:
              WENOZX2(member, eos_, false, m, k, j, il, iu, wder_, dl_jp1, dr);
              break;
            default:
              break;
          }
          }
          par_for_inner(member, il, iu, [&](const int i) {
            dl_jp1(IDPR,i) = fmax(dl_jp1(IDPR,i), eos_.pfloor);
            dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
          });
        }

          if (use_wellbalance_static_reconst_perturb_) {
            AddWbPrimFaceX2(member,m,k,j,il,iu,w0facewb_x2f,wl_jp1,wr);
          }

          if (use_cubed_sphere) {
              GnomonicEquianglePrimFaceX2(gtrig_xi,member,m,k,j,il,iu,wl_jp1,wr);
          }
          
        member.team_barrier();

        // compute fluxes over [js,je+1].  RS returns flux in input wr array
        if (j>jl) {
          // NOTE(@pdmullen): Capture variables prior to if constexpr.
          auto eos = eos_;
          const bool nanrep_r = nanrep_r_;
          auto rscnt = rscnt_;
          auto rsrec = rsrec_;
          auto indcs = indcs_;
          auto size = size_;
          auto coord = coord_;
          auto flx2 = flx2_;
          if constexpr (rsolver_method_ == Hydro_RSolver::advect) {
            Advect(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::llf) {
            LLF(member,eos,indcs,size,coord,m,k,j,il,iu,IVY,wl,wr,dl,dr,flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle) {
            HLLE(member,eos,indcs,size,coord,m,k,j,il,iu,IVY,wl,wr,dl,dr,flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hllc) {
            HLLC(member,eos,indcs,size,coord,m,k,j,il,iu,IVY,wl,wr,dl,dr,flx2,
                 nanrep_r,rscnt,rsrec);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::lhllc) {
            LHLLC(member,eos,indcs,size,coord,m,k,j,il,iu,IVY,wl,wr,dl,dr,flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hllclm) {
            HLLCLM(member,eos,indcs,size,coord,m,k,j,il,iu,IVY,wl,wr,dl,dr,flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::ausmpup) {
            AUSMPUP(member,eos,indcs,size,coord,m,k,j,il,iu,IVY,wl,wr,dl,dr,flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::roe) {
            Roe(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::llf_sr) {
            LLF_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle_sr) {
            HLLE_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hllc_sr) {
            HLLC_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::llf_gr) {
            LLF_GR(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle_gr) {
            HLLE_GR(member, eos, indcs, size, coord, m, k, j, il, iu, IVY, wl, wr, flx2);
          }
            if (use_cubed_sphere) {
                GnomonicEquiangleFluxX2(gtrig_xi,member,m,k,j,il,iu,flx2);
            }
          member.team_barrier();
        }

        // calculate fluxes of scalars (if any)
        if (nvars > nhyd_) {
          for (int n=nhyd_; n<nvars; ++n) {
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
                ScrArray2D<Real>::shmem_size(nder, ncells1)) * 3;
    auto &flx3_ = uflx.x3f;
      
    auto &x3v_ = pmy_pack->pcoord->x3v;
    auto &x3f_ = pmy_pack->pcoord->xx3f;
    auto &x2v3_ = pmy_pack->pcoord->x2v;
    auto &x2f3_ = pmy_pack->pcoord->xx2f;

    // set the loop limits
    il = is, iu = ie, jl = js, ju = je, kl = ks-1, ku = ke+1;
    if (use_fofc) { il = is-1, iu = ie+1, jl = js-1, ju = je+1, kl = ks-2, ku = ke+2; }
    // the face-average stencils read the x3 fluxes at i +- 1 and j +- 1 (the polar
    // ghost row included)
    else if (sp_favg_) { il = is-1, iu = ie+1, jl = js-1, ju = je+1; }

    par_for_outer("hflux_x3",DevExeSpace(), scr_size, scr_level, 0, nmb1, jl, ju,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j) {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> dscr1(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr2(member.team_scratch(scr_level), nder, ncells1);
      ScrArray2D<Real> dscr3(member.team_scratch(scr_level), nder, ncells1);

      for (int k=kl; k<=ku; ++k) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_kp1 = scr2;
        auto wr     = scr3;
        auto dl     = dscr1;
        auto dl_kp1 = dscr2;
        auto dr     = dscr3;
        if ((k%2) == 0) {
          wl     = scr2;
          wl_kp1 = scr1;
          dl     = dscr2;
          dl_kp1 = dscr1;
        }

        if (use_spherical_polar) {
          GridPiecewiseLinearX3(member, m, k, j, is-1, ie+1, w0_, x3v_, x3f_, wl_kp1, wr);
        } else {
          
        if (use_wellbalance_dynamic_ && use_wb_x3_)
        {
          WbLocalPiecewiseLinearX3(member, eos_, wb_option_, use_wb_rho_,
              m, k, j, il, iu, w0_, phicc0_, phi0_x3f, wl_kp1, wr);
        } else {

        // Reconstruct qR[k] and qL[k+1]
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX3(member, m, k, j, il, iu, w0_, wl_kp1, wr);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX3(member, m, k, j, il, iu, w0_, wl_kp1, wr);
            break;
          case ReconstructionMethod::ppm4:
          case ReconstructionMethod::ppmx:
            PiecewiseParabolicX3(member,eos_,extrema,true,m,k,j,il,iu, w0_, wl_kp1, wr);
            break;
          case ReconstructionMethod::wenoz:
            WENOZX3(member, eos_, true, m, k, j, il, iu, w0_, wl_kp1, wr);
            break;
          default:
            break;
        }
              
        }
        }
          
        // Reconstruct the derived thermodynamic variables (general EOS only)
        if (nder > 0) {
          if (use_spherical_polar) {
            GridPiecewiseLinearX3(member, m, k, j, il, iu, wder_, x3v_, x3f_, dl_kp1, dr);
          } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX3(member, m, k, j, il, iu,
                                     pwb_, pfacewb_x3f,
                                    wder_, dl_kp1, dr);
      } else if (use_wellbalance_dynamic_ && use_wb_x3_) {
            WbPiecewiseLinearDerX3(member, eos_, wb_option_,
                m, k, j, il, iu, w0_, wder_,
                                   phicc0_, phi0_x3f, dl_kp1, dr);
          } else {
          switch (recon_method_) {
            case ReconstructionMethod::dc:
              DonorCellX3(member, m, k, j, il, iu, wder_, dl_kp1, dr);
              break;
            case ReconstructionMethod::plm:
              PiecewiseLinearX3(member, m, k, j, il, iu, wder_, dl_kp1, dr);
              break;
            case ReconstructionMethod::ppm4:
            case ReconstructionMethod::ppmx:
              PiecewiseParabolicX3(member,eos_,extrema,false,m,k,j,il,iu,wder_,dl_kp1,dr);
              break;
            case ReconstructionMethod::wenoz:
              WENOZX3(member, eos_, false, m, k, j, il, iu, wder_, dl_kp1, dr);
              break;
            default:
              break;
          }
          }
          par_for_inner(member, il, iu, [&](const int i) {
            dl_kp1(IDPR,i) = fmax(dl_kp1(IDPR,i), eos_.pfloor);
            dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
          });
        }

          if (use_wellbalance_static_reconst_perturb_) {
            AddWbPrimFaceX3(member,m,k,j,il,iu,w0facewb_x3f,wl_kp1,wr);
          }

          if (use_cubed_sphere) {
              GnomonicEquianglePrimFaceX3(gtrig_eta,member,m,k,j,il,iu,wl_kp1,wr);
          }
          
        // SPHERICAL POLAR: correct the reconstructed x3-face states from the cell's
        // volume centroid theta to the x3 face's midpoint theta (see the long note in
        // mhd_fluxes.cpp).  <mesh>/polar_x3_shift selects the form: 1 = "rotate", the
        // exact basis rotation of the (r,theta) vector pair; 2 = "shift", the unlimited
        // centred-difference theta shift of every variable, which also carries the
        // scalar gradient but was BISECTED to the polar-row radial-field blow-up of the
        // sp MHD production run (2026-09-10); 3 = "scalar", rotate the vectors and shift
        // only the remaining variables; 0 = "off".
        if (use_spherical_polar) {
          if (px3_ == 1 || px3_ == 3) {
            const Real dang = 0.5*(x2f3_(m,j) + x2f3_(m,j+1)) - x2v3_(m,j);
            const Real cd = cos(dang), sd = sin(dang);
            par_for_inner(member, il, iu, [&](const int i) {
              Real a = wl_kp1(IVX,i), b = wl_kp1(IVY,i);
              wl_kp1(IVX,i) = a*cd + b*sd;  wl_kp1(IVY,i) = -a*sd + b*cd;
              a = wr(IVX,i); b = wr(IVY,i);
              wr(IVX,i) = a*cd + b*sd;  wr(IVY,i) = -a*sd + b*cd;
            });
          }
          if (px3_ == 2 || px3_ == 3) {
            const Real fac = (0.5*(x2f3_(m,j) + x2f3_(m,j+1)) - x2v3_(m,j))
                             /(x2v3_(m,j+1) - x2v3_(m,j-1));
            for (int n=0; n<nvars; ++n) {
              if (px3_ == 3 && (n == IVX || n == IVY)) continue;
              par_for_inner(member, il, iu, [&](const int i) {
                const Real sh = fac*(w0_(m,n,k,j+1,i) - w0_(m,n,k,j-1,i));
                wl_kp1(n,i) += sh;
                wr(n,i) += sh;
              });
            }
          }
        }

        member.team_barrier();

        // compute fluxes over [ks,ke+1].  RS returns flux in input wr array
        if (k>kl) {
          // NOTE(@pdmullen): Capture variables prior to if constexpr.
          auto eos = eos_;
          const bool nanrep_r = nanrep_r_;
          auto rscnt = rscnt_;
          auto rsrec = rsrec_;
          auto indcs = indcs_;
          auto size = size_;
          auto coord = coord_;
          auto flx3 = flx3_;
          if constexpr (rsolver_method_ == Hydro_RSolver::advect) {
            Advect(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::llf) {
            LLF(member,eos,indcs,size,coord,m,k,j,il,iu,IVZ,wl,wr,dl,dr,flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle) {
            HLLE(member,eos,indcs,size,coord,m,k,j,il,iu,IVZ,wl,wr,dl,dr,flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hllc) {
            HLLC(member,eos,indcs,size,coord,m,k,j,il,iu,IVZ,wl,wr,dl,dr,flx3,
                 nanrep_r,rscnt,rsrec);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::lhllc) {
            LHLLC(member,eos,indcs,size,coord,m,k,j,il,iu,IVZ,wl,wr,dl,dr,flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hllclm) {
            HLLCLM(member,eos,indcs,size,coord,m,k,j,il,iu,IVZ,wl,wr,dl,dr,flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::ausmpup) {
            AUSMPUP(member,eos,indcs,size,coord,m,k,j,il,iu,IVZ,wl,wr,dl,dr,flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::roe) {
            Roe(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::llf_sr) {
            LLF_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle_sr) {
            HLLE_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hllc_sr) {
            HLLC_SR(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::llf_gr) {
            LLF_GR(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          } else if constexpr (rsolver_method_ == Hydro_RSolver::hlle_gr) {
            HLLE_GR(member, eos, indcs, size, coord, m, k, j, il, iu, IVZ, wl, wr, flx3);
          }
            if (use_cubed_sphere) {
                GnomonicEquiangleFluxX3(gtrig_eta,member,m,k,j,il,iu,flx3);
            }
          member.team_barrier();
        }

        // calculate fluxes of scalars (if any)
        if (nvars > nhyd_) {
          for (int n=nhyd_; n<nvars; ++n) {
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
  // fluxes (see Coordinates::SphericalPolarFaceAverageFluxes)
  if (pmy_pack->pcoord->sp_cart_all_momentum) {
    pmy_pack->pcoord->SphericalPolarFaceAverageFluxes(uflx, nhydro);
  }

  return;
}

// function definitions for each template parameter
template void Hydro::CalculateFluxes<Hydro_RSolver::advect>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::llf>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::hlle>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::hllc>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::lhllc>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::hllclm>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::ausmpup>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::roe>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::llf_sr>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::hlle_sr>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::hllc_sr>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::llf_gr>(Driver *pdriver, int stage);
template void Hydro::CalculateFluxes<Hydro_RSolver::hlle_gr>(Driver *pdriver, int stage);

} // namespace hydro
