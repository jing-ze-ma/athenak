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
#include "eos/eos.hpp"
#include "reconstruct/dc.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
#include "mhd/rsolvers/advect_mhd.hpp"
#include "mhd/rsolvers/cs_lowbeta_fallback.hpp"
#include "mhd/rsolvers/llf_mhd.hpp"
#include "mhd/rsolvers/hlle_mhd.hpp"
#include "mhd/rsolvers/hlld_mhd.hpp"
#include "mhd/rsolvers/llf_srmhd.hpp"
#include "mhd/rsolvers/hlle_srmhd.hpp"
#include "mhd/rsolvers/llf_grmhd.hpp"
#include "mhd/rsolvers/hlle_grmhd.hpp"
// #include "mhd/rsolvers/roe_mhd.hpp"

namespace mhd {

void ReportLowBetaDiag(const DvceArray2D<Real> &lbd, const int nbin, const Real thresh);

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
  auto phi0_x1f = phi0.x1f;
  auto phi0_x2f = phi0.x2f;
  auto phi0_x3f = phi0.x3f;
    
  auto &use_cubed_sphere = pmy_pack->pmesh->use_cubed_sphere;
  const Real cs_lb_beta = cs_lowbeta_fallback;
  // MEASUREMENT ONLY: count the fallback's hits on ONE stage-1 sweep every N cycles, so
  // the printed profile is a snapshot of a single sweep rather than a running sum whose
  // denominator depends on how long the run has been going.
  const int lbd_on = (cs_lowbeta_diag > 0 && stage == 1 &&
                      (pmy_pack->pmesh->ncycle % cs_lowbeta_diag) == 0) ? 1 : 0;
  const int lbd_nbin = indcs_.nx1;
  auto lbd_ = lb_diag;
  if (lbd_on) {
    const Real big = std::numeric_limits<Real>::max();
    par_for("lbdiag_zero", DevExeSpace(), 0, lbd_nbin-1, KOKKOS_LAMBDA(const int b) {
      lbd_(b,0) = 0.0; lbd_(b,1) = 0.0; lbd_(b,2) = big;
    });
  }
  const auto wb_option_ = wb_option;
  const bool use_wb_rho_ = use_wb_rho;
  const bool use_wb_x1_ = use_wb_x1;
  const bool use_wb_x2_ = use_wb_x2;
  const bool use_wb_x3_ = use_wb_x3;
  const bool use_wellbalance_dynamic_ = use_wellbalance_dynamic;
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
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool pquad_ = pmy_pack->pmesh->use_polar_quadratic_recon;

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

  par_for_outer("mhd_flux1",DevExeSpace(), scr_size, scr_level, 0, nmb1, kl, ku, jl, ju,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray2D<Real> wl(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> wr(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> bl(member.team_scratch(scr_level), 3, ncells1);
    ScrArray2D<Real> br(member.team_scratch(scr_level), 3, ncells1);
    ScrArray2D<Real> dl(member.team_scratch(scr_level), nder, ncells1);
    ScrArray2D<Real> dr(member.team_scratch(scr_level), nder, ncells1);

    if (use_spherical_polar || str_r1_) {
      GridPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
                            use_wellbalance_dynamic_,
                            use_wb_x1_, m, k, j, il-1, iu, w0_, x1v_, x1f_, phicc0_,
                            phi0_x1f, true, wl, wr);
      GridPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
                            use_wellbalance_dynamic_,
                            use_wb_x1_, m, k, j, il-1, iu, b0_, x1v_, x1f_, phicc0_,
                            phi0_x1f, false, bl, br);
    } else {
        
    if (use_wellbalance_dynamic_ && use_wb_x1_)
    {
      WbLocalPiecewiseLinearX1(member, eos_, wb_option_, use_wb_rho_,
          m, k, j, il-1, iu, w0_, phicc0_, phi0_x1f, wl, wr);
      PiecewiseLinearX1(member, m, k, j, il-1, iu, b0_, bl, br);
    } else {

    // Reconstruct qR[i] and qL[i+1], for both W and Bcc
    switch (recon_method_) {
      case ReconstructionMethod::dc:
        DonorCellX1(member, m, k, j, il-1, iu, w0_, wl, wr);
        DonorCellX1(member, m, k, j, il-1, iu, b0_, bl, br);
        break;
      case ReconstructionMethod::plm:
        PiecewiseLinearX1(member, m, k, j, il-1, iu, w0_, wl, wr);
        PiecewiseLinearX1(member, m, k, j, il-1, iu, b0_, bl, br);
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
      if (use_spherical_polar || str_r1_) {
        GridPiecewiseLinearDerX1(member, eos_, wb_option_, use_wellbalance_dynamic_,
                                 use_wb_x1_, m, k, j, il-1, iu, w0_, wder_,
                                 x1v_, x1f_, phicc0_, phi0_x1f, dl, dr);
      } else if (use_wellbalance_static_reconst_perturb_) {
        WbStaticPiecewiseLinearDerX1(member, m, k, j, il-1, iu,
                                     pwb_, pfacewb_x1f,
                                    wder_, dl, dr);
      } else if (use_wellbalance_dynamic_ && use_wb_x1_) {
        WbPiecewiseLinearDerX1(member, eos_, wb_option_,
            m, k, j, il-1, iu, w0_, wder_,
                               phicc0_, phi0_x1f, dl, dr);
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

    // Sync all threads in the team so that scratch memory is consistent
    member.team_barrier();

    // compute fluxes over [is,ie+1].  MHD RS also computes electric fields, where
    // (IBY) component of flx = E_{z} = -(v x B)_{z} = -(v1*b2 - v2*b1)
    // (IBZ) component of flx = E_{y} = -(v x B)_{y} =  (v1*b3 - v3*b1)
    // NOTE(@pdmullen): Capture variables prior to if constexpr.  Required for cuda 11.6+.
    auto eos = eos_;
    auto indcs = indcs_;
    auto size = size_;
    auto coord = coord_;
    auto bx = bx_;
    auto flx1 = flx1_;
    auto e31 = e31_;
    auto e21 = e21_;
    const bool do_pole = (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je);
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
    } else if constexpr (rsolver_method_ == MHD_RSolver::llf_sr) {
      LLF_SR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_sr) {
      HLLE_SR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::llf_gr) {
      LLF_GR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    } else if constexpr (rsolver_method_ == MHD_RSolver::hlle_gr) {
      HLLE_GR(member,eos,indcs,size,coord,m,k,j,il,iu,IVX,wl,wr,bl,br,bx,flx1,e31,e21);
    }

    if (use_cubed_sphere) {
      // low-beta dissipation fallback, BEFORE the rotation into covariant slots
      if (cs_lb_beta > 0.0) {
        CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,il,iu,IVX,cs_lb_beta,
                     wl,wr,bl,br,dl,dr,bx,flx1,e31,e21,lbd_on,is,lbd_nbin,lbd_);
        member.team_barrier();
      }
      GnomonicEquiangleFluxX1(gtrig_cell,member,m,k,j,il,iu,flx1);
      GnomonicEquiangleEmfX1(gtrig_cell,member,m,k,j,il,iu,e31,e21);
    }
    member.team_barrier();

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
      
    auto &x2v_ = pmy_pack->pcoord->x2v;
    auto &x2f_ = pmy_pack->pcoord->xx2f;

    // set the loop limits for 2D/3D problems
    if (pmy_pack->pmesh->two_d) {
      kl = ks, ku = ke;
    } else { // 3D
      kl = ks-1, ku = ke+1;
    }
    jl = js-1, ju = je+1;
    if (use_fofc) { jl = js-2, ju = je+2; }

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
          PiecewiseLinearX2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
        } else {

        // Reconstruct qR[j] and qL[j+1], for both W and Bcc
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX2(member, m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            DonorCellX2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX2(member, m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            PiecewiseLinearX2(member, m, k, j, is-1, ie+1, b0_, bl_jp1, br);
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
              PiecewiseLinearX2(member, m, k, j, is-1, ie+1, wder_, dl_jp1, dr);
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

        // compute fluxes over [js,je+1].  MHD RS also computes electric fields, where
        // (IBY) component of flx = E_{x} = -(v x B)_{x} = -(v2*b3 - v3*b2)
        // (IBZ) component of flx = E_{z} = -(v x B)_{z} =  (v2*b1 - v1*b2)
        if (j>jl) {
          // NOTE(@pdmullen): Capture variables prior to if constexpr.
          auto eos = eos_;
          auto indcs = indcs_;
          auto size = size_;
          auto coord = coord_;
          auto by = by_;
          auto flx2 = flx2_;
          auto e12 = e12_;
          auto e32 = e32_;
          const bool do_pole = (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && (j == js || j == js+1)) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && (j == je+1 || j == je));
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

          if (use_cubed_sphere) {
            if (cs_lb_beta > 0.0) {
              CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVY,cs_lb_beta,
                           wl,wr,bl,br,dl,dr,by,flx2,e12,e32,lbd_on,is,lbd_nbin,lbd_);
              member.team_barrier();
            }
            GnomonicEquiangleFluxX2(gtrig_xi,member,m,k,j,is-1,ie+1,flx2);
          }
          member.team_barrier();
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
        } else {
            
        if (use_wellbalance_dynamic_ && use_wb_x3_)
        {
          // is-1,ie+1, not il,iu -- the b0_ line below, the spherical-polar branch
          // above, the whole default switch, and the Riemann solver all use
          // is-1,ie+1.  With il,iu (the x1 sweep's limits, never reset for x2/x3)
          // this leaves wl_kp1/wr unwritten at i = is-1, which the solver reads.
          WbLocalPiecewiseLinearX3(member, eos_, wb_option_, use_wb_rho_,
              m, k, j, is-1, ie+1, w0_,
                                   phicc0_, phi0_x3f, wl_kp1, wr);
          PiecewiseLinearX3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
        } else {

        // Reconstruct qR[k] and qL[k+1], for both W and Bcc
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX3(member, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            DonorCellX3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX3(member, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            PiecewiseLinearX3(member, m, k, j, is-1, ie+1, b0_, bl_kp1, br);
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
      } else if (use_wellbalance_dynamic_ && use_wb_x3_) {
            WbPiecewiseLinearDerX3(member, eos_, wb_option_,
                m, k, j, is-1, ie+1, w0_, wder_,
                                   phicc0_, phi0_x3f, dl_kp1, dr);
          } else {
          switch (recon_method_) {
            case ReconstructionMethod::dc:
              DonorCellX3(member, m, k, j, is-1, ie+1, wder_, dl_kp1, dr);
              break;
            case ReconstructionMethod::plm:
              PiecewiseLinearX3(member, m, k, j, is-1, ie+1, wder_, dl_kp1, dr);
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
        if (use_spherical_polar) {
          const Real fac = (0.5*(x2f3_(m,j) + x2f3_(m,j+1)) - x2v3_(m,j))
                           /(x2v3_(m,j+1) - x2v3_(m,j-1));
          for (int n=0; n<nvars; ++n) {
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              const Real sh = fac*(w0_(m,n,k,j+1,i) - w0_(m,n,k,j-1,i));
              wl_kp1(n,i) += sh;
              wr(n,i) += sh;
            });
          }
          for (int n=0; n<3; ++n) {
            par_for_inner(member, is-1, ie+1, [&](const int i) {
              const Real sh = fac*(b0_(m,n,k,j+1,i) - b0_(m,n,k,j-1,i));
              bl_kp1(n,i) += sh;
              br(n,i) += sh;
            });
          }
        }

        member.team_barrier();

        // compute fluxes over [ks,ke+1].  MHD RS also computes electric fields, where
        // (IBY) component of flx = E_{y} = -(v x B)_{y} = -(v3*b1 - v1*b3)
        // (IBZ) component of flx = E_{x} = -(v x B)_{x} =  (v3*b2 - v2*b3)
        if (k>kl) {
          // NOTE(@pdmullen): Capture variables prior to if constexpr.
          auto eos = eos_;
          auto indcs = indcs_;
          auto size = size_;
          auto coord = coord_;
          auto bz = bz_;
          auto flx3 = flx3_;
          auto e23 = e23_;
          auto e13 = e13_;
          const bool do_pole = (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je);
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

          if (use_cubed_sphere) {
            if (cs_lb_beta > 0.0) {
              CSLowBetaFallback(member,eos,eos.IsGeneral(),m,k,j,is-1,ie+1,IVZ,cs_lb_beta,
                           wl,wr,bl,br,dl,dr,bz,flx3,e23,e13,lbd_on,is,lbd_nbin,lbd_);
              member.team_barrier();
            }
            GnomonicEquiangleFluxX3(gtrig_eta,member,m,k,j,is-1,ie+1,flx3);
          }
          member.team_barrier();
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

  if (lbd_on) { ReportLowBetaDiag(lbd_, lbd_nbin, cs_lb_beta); }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ReportLowBetaDiag
//! \brief print the low-beta fallback's hit rate as a profile in radius.
//!
//! MEASUREMENT ONLY; see mhd.hpp and rsolvers/cs_lowbeta_fallback.hpp.  The counters are
//! summed over every MeshBlock and all three flux directions on this sweep, so "faces" is
//! the number of faces the switch examined, not the number of cells.  Radial bins are
//! grouped for printing: nx1 = 128 rows would bury the answer.

void ReportLowBetaDiag(const DvceArray2D<Real> &lbd, const int nbin, const Real thresh) {
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), lbd);
  std::vector<double> nface(nbin), nsw(nbin), bmin(nbin);
  for (int b=0; b<nbin; ++b) {
    nface[b] = h(b,0); nsw[b] = h(b,1); bmin[b] = h(b,2);
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, nface.data(), nbin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, nsw.data(),   nbin, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, bmin.data(),  nbin, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif
  if (global_variable::my_rank != 0) return;

  const int ngroup = (nbin < 16) ? nbin : 16;
  const Real big = 0.5*static_cast<Real>(std::numeric_limits<Real>::max());
  double tface = 0.0, tsw = 0.0;
  std::cout << "## cs_lowbeta_diag (threshold beta < " << thresh
            << "): faces examined / switched to HLLE, by radial index" << std::endl;
  std::cout << "##   i-range        faces      switched   frac     min beta" << std::endl;
  for (int g=0; g<ngroup; ++g) {
    const int lo = (g*nbin)/ngroup, hi = ((g+1)*nbin)/ngroup;
    double f = 0.0, w = 0.0, bm = std::numeric_limits<double>::max();
    for (int b=lo; b<hi; ++b) {
      f += nface[b]; w += nsw[b];
      if (bmin[b] < bm) { bm = bmin[b]; }
    }
    tface += f; tsw += w;
    std::cout << "##  " << std::setw(4) << lo << "-" << std::setw(4) << hi
              << std::setw(14) << static_cast<std::int64_t>(f)
              << std::setw(14) << static_cast<std::int64_t>(w)
              << std::setw(9) << std::fixed << std::setprecision(4)
              << ((f > 0.0) ? w/f : 0.0);
    // a bin with no field anywhere never recorded a beta at all
    if (bm < big) {
      std::cout << "  " << std::scientific << std::setprecision(3) << bm << std::endl;
    } else {
      std::cout << "        --" << std::endl;
    }
  }
  std::cout << "##  total" << std::setw(18) << static_cast<std::int64_t>(tface)
            << std::setw(14) << static_cast<std::int64_t>(tsw)
            << std::setw(9) << std::fixed << std::setprecision(4)
            << ((tface > 0.0) ? tsw/tface : 0.0) << std::endl;
  std::cout.unsetf(std::ios_base::floatfield);
}

// function definitions for each template parameter
template void MHD::CalculateFluxes<MHD_RSolver::advect>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::llf>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlle>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlld>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::llf_sr>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlle_sr>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::llf_gr>(Driver *pdriver, int stage);
template void MHD::CalculateFluxes<MHD_RSolver::hlle_gr>(Driver *pdriver, int stage);

} // namespace mhd
