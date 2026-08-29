//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_seam_cc.cpp
//! \brief makes the flux of cell-centered variables through a cubed-sphere PANEL SEAM
//! single-valued, so that the seam is exactly conservative.
//!
//! WHY THIS IS NEEDED. A finite-volume update telescopes -- what leaves one cell enters
//! its neighbour, to round-off -- only because both cells multiply the SAME stored flux
//! by the SAME stored area. Inside a panel that is true by construction, and it is true
//! between MeshBlocks of one panel as well, because their ghost zones are COPIES and so
//! reproduce the neighbour's Riemann solve bit for bit. At a panel seam it is false: the
//! ghosts there are INTERPOLATED (see the along-seam resample in bvals_cc.cpp), so the
//! two panels reconstruct slightly different states either side of the one physical face
//! and their fluxes disagree at the order of the halo's interpolation error. The
//! difference is mass, momentum and energy created out of nothing, and it accumulates
//! secularly. Measured before this exchange existed, it was the dominant term in the
//! global mass and energy drift of the rigid-rotation test, and it was worst at the cube
//! vertices. Note this is exactly the role `SendFlux`/`RecvFlux` already play at a
//! fine/coarse boundary -- but that path runs only when `multilevel`, and a uniform
//! cubed-sphere grid is not multilevel, so nothing reconciled a seam.
//!
//! WHAT IS EXCHANGED. Each panel packs the OUTWARD flux through its copy of the shared
//! face, `sgn*flx`, and each replaces its own with the average of the two outward values,
//!     F_out <- 0.5*(F_out - F_out_neighbour).
//! Both sides evaluate the same expression, so afterwards what leaves one panel is
//! exactly minus what leaves the other and the sum over the sphere telescopes to the
//! physical boundaries alone. No ordering or single-owner rule is needed.
//!
//! THE MOMENTUM COMPONENTS NEED THE BASIS TRANSFORM, the scalars do not. Mass and total
//! energy flux densities are true scalars, and so is the radial momentum flux, since rhat
//! is common to both charts. The two TANGENTIAL momentum fluxes are components on the
//! local gnomonic basis, and the two charts carry different tangent bases at the same
//! physical point, so they go through `cubed_sphere::TransformMomentum` -- the same
//! transform the halo applies to the momentum itself.
//!
//! NO ALONG-SEAM RESAMPLE IS NEEDED, unlike the halo. The halo interpolates because a
//! ghost CELL CENTRE of one panel does not sit at a cell centre of the other. The seam
//! FACES do coincide: both charts parameterise the shared cube edge identically, so the
//! faces pair one-to-one and a signed index map is exact. (Verified geometrically by the
//! `CS SEAM FLUX MISMATCH` gate in `pgen/cs_test.cpp`, which pairs faces by physical
//! position and reports a max area mismatch at round-off.)
//!
//! It is safe to reuse the `flux` buffers, `comm_flux` and `flux_req` here: mesh
//! refinement and the cubed sphere are mutually exclusive (a startup FATAL in
//! `mesh.cpp`), so the fine/coarse path in `flux_correct_cc.cpp` can never be live at the
//! same time as this one.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"

namespace {
//! Classify a boundary buffer from its own stored index range, rather than from the
//! buffer NUMBER. The numbering is not what it looks like -- x1 faces are 0-7, x2 faces
//! 8-15, then the x1x2 EDGES are 16-23 and the x3 faces do not start until 24 -- so a
//! plain "n < 24 is a face" test silently drops every x3 seam. The index range says it
//! unambiguously: a face buffer is one cell thick along its own normal and spans the full
//! active range along the other two axes.
//!
//! Returns false unless this is an x2 or x3 FACE. Only those can be a seam: x1 is radial
//! and no seam crosses it, and an edge or corner buffer carries no face flux at all.
//! `x2face` says which of the two it is, and `sgn` is the sign that turns the stored flux
//! (always the flux in the +x2 or +x3 direction) into the OUTWARD flux.
KOKKOS_INLINE_FUNCTION
bool SeamFaceGeom(const MeshBufferIndcs &ix, const int is, const int ie,
                  const int js, const int ks, bool &x2face, Real &sgn) {
  // A face buffer spans the FULL active range in x1. Testing only "thin in j XOR thin in
  // k" is not enough: an x1x2 edge buffer is thin in j and full in k, and so looks
  // exactly like an x2 face. That costs nothing while x1 is a single MeshBlock -- the
  // radial direction is bounded by physical boundaries, so those edge buffers have no
  // neighbour -- but it would hand an edge buffer to the face packer the moment x1 is
  // split, with the wrong extents and the wrong flux array.
  if (ix.bis != is || ix.bie != ie) {return false;}
  const bool jthin = (ix.bjs == ix.bje);
  const bool kthin = (ix.bks == ix.bke);
  if (jthin == kthin) {return false;}     // x1 face, or an x2x3 edge
  x2face = jthin;
  sgn = (x2face ? (ix.bjs == js) : (ix.bks == ks)) ? -1.0 : 1.0;
  return true;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::InitFluxSeamRecv
//! \brief Post non-blocking receives for the seam flux exchange.

TaskStatus MeshBoundaryValuesCC::InitFluxSeamRecv(const int nvars) {
#if MPI_PARALLEL_ENABLED
  int &nmb = pmy_pack->nmb_thispack;
  int &nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const int is_ = pmy_pack->pmesh->mb_indcs.is;
  const int ie_ = pmy_pack->pmesh->mb_indcs.ie;
  const int js_ = pmy_pack->pmesh->mb_indcs.js;
  const int ks_ = pmy_pack->pmesh->mb_indcs.ks;

  bool no_errors=true;
  bool x2f; Real sg;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ((nghbr.h_view(m,n).gid >= 0) &&
          SeamFaceGeom(recvbuf[n].iflux_same[0], is_, ie_, js_, ks_, x2f, sg) &&
          (nghbr.h_view(m,n).panel != mbpanel.h_view(m))) {
        int drank = nghbr.h_view(m,n).rank;
        if (drank != global_variable::my_rank) {
          // tag uses the local ID and buffer index of the *receiving* MeshBlock
          int tag = CreateBvals_MPI_Tag(m, n);
          int data_size = nvars*(recvbuf[n].iflxs_ndat);
          auto recv_ptr = Kokkos::subview(recvbuf[n].flux, m, Kokkos::ALL);
          int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_flux, &(recvbuf[n].flux_req[m]));
          if (ierr != MPI_SUCCESS) {no_errors=false;}
        }
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::PackAndSendFluxSeamCC
//! \brief Pack the OUTWARD flux through each panel-seam face and send it to the panel on
//! the other side.
//!
//! The index map across a seam is the same signed permutation the halo uses -- the
//! along-seam index may reverse, and across a "swap" seam the two tangential axes
//! exchange roles, which is absorbed here by transposing the buffer strides so that the
//! receiver's generic unpack lands the data on its own j,k.

TaskStatus MeshBoundaryValuesCC::PackAndSendFluxSeamCC(DvceFaceFld5D<Real> &flx) {
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  int nvar = flx.x1f.extent_int(1);

  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;

  Kokkos::TeamPolicy<> policy(DevExeSpace(), (nmb*nnghbr*nvar), Kokkos::AUTO);
  Kokkos::parallel_for("SeamFluxSend", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    if (nghbr.d_view(m,n).gid < 0) {return;}
    bool x2face; Real sgn;
    if (!SeamFaceGeom(sbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                      cs_indcs.js, cs_indcs.ks, x2face, sgn)) {
      return;
    }
    const int my_panel = mbpanel.d_view(m);
    const int dst_panel = nghbr.d_view(m,n).panel;
    if (dst_panel == my_panel) {return;}

    const int il = sbuf[n].iflux_same[0].bis;
    const int iu = sbuf[n].iflux_same[0].bie;
    const int jl = sbuf[n].iflux_same[0].bjs;
    const int ju = sbuf[n].iflux_same[0].bje;
    const int kl = sbuf[n].iflux_same[0].bks;
    const int ku = sbuf[n].iflux_same[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nkj = nk*nj;

    const int dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
    const int dn = nghbr.d_view(m,n).dest;

    // the seam index map, exactly as in MeshBoundaryValuesCC::PackAndSendCC
    const PanelBoundaries pb = GetPanelBoundary(my_panel, dst_panel);
    int aj = 1, bj = 0, ak = 1, bk = 0, sj = 1, sk = nj;
    const int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
    const int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
    if (pb.swap_ax == 1) {sj = nk; sk = 1;}
    if (rev_a_preswap) {aj = -1; bj = jl + ju;}
    if (rev_b_preswap) {ak = -1; bk = kl + ku;}

    // Mass and energy fluxes are true scalars, and so is the RADIAL momentum flux since
    // rhat is common to both charts. Only the two tangential momentum fluxes are
    // components on a basis the seam changes.
    const bool xform = (v == IM2) || (v == IM3);

    Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
      int k = idx / nj;
      int j = (idx - k * nj) + jl;
      k += kl;
      const int jj = aj*j + bj;
      const int kk = ak*k + bk;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1), [&](const int i) {
        Real val;
        if (!xform) {
          val = x2face ? flx.x2f(m,v,kk,jj,i) : flx.x3f(m,v,kk,jj,i);
        } else {
          // (xi,eta) of the FACE CENTRE: the seam-normal coordinate is a cell EDGE.
          const int js_ = cs_indcs.js, ks_ = cs_indcs.ks;
          const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
          const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
          const Real xi = 0.25*M_PI*(x2face
              ? LeftEdgeX(jj-js_, cs_indcs.nx2, x2mn, x2mx)
              : CellCenterX(jj-js_, cs_indcs.nx2, x2mn, x2mx));
          const Real eta = 0.25*M_PI*(x2face
              ? CellCenterX(kk-ks_, cs_indcs.nx3, x3mn, x3mx)
              : LeftEdgeX(kk-ks_, cs_indcs.nx3, x3mn, x3mx));
          const Real f2 = x2face ? flx.x2f(m,IM2,kk,jj,i) : flx.x3f(m,IM2,kk,jj,i);
          const Real f3 = x2face ? flx.x2f(m,IM3,kk,jj,i) : flx.x3f(m,IM3,kk,jj,i);
          Real f2o, f3o;
          cubed_sphere::TransformMomentum(my_panel, dst_panel, xi, eta, f2, f3, f2o, f3o);
          val = (v == IM2) ? f2o : f3o;
        }
        const int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*v);
        if (nghbr.d_view(m,n).rank == my_rank) {
          rbuf[dn].flux(dm, index) = sgn*val;
        } else {
          sbuf[n].flux(m, index) = sgn*val;
        }
      });
    });
  });

#if MPI_PARALLEL_ENABLED
  Kokkos::fence();
  bool no_errors=true;
  auto &nghbr_h = pmy_pack->pmb->nghbr;
  auto &mbpanel_h = pmy_pack->pmb->mb_panel;
  const int is_h = pmy_pack->pmesh->mb_indcs.is;
  const int ie_h = pmy_pack->pmesh->mb_indcs.ie;
  const int js_h = pmy_pack->pmesh->mb_indcs.js;
  const int ks_h = pmy_pack->pmesh->mb_indcs.ks;
  bool x2f_h; Real sg_h;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ((nghbr_h.h_view(m,n).gid >= 0) &&
          SeamFaceGeom(sendbuf[n].iflux_same[0], is_h, ie_h, js_h, ks_h, x2f_h, sg_h) &&
          (nghbr_h.h_view(m,n).panel != mbpanel_h.h_view(m))) {
        int dn = nghbr_h.h_view(m,n).dest;
        int drank = nghbr_h.h_view(m,n).rank;
        if (drank != my_rank) {
          int lid = nghbr_h.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);
          int data_size = nvar*(sendbuf[n].iflxs_ndat);
          auto send_ptr = Kokkos::subview(sendbuf[n].flux, m, Kokkos::ALL);
          int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_flux, &(sendbuf[n].flux_req[m]));
          if (ierr != MPI_SUCCESS) {no_errors=false;}
        }
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::RecvAndUnpackFluxSeamCC
//! \brief Replace the flux on each seam face by the average of the two panels' OUTWARD
//! values, so that what leaves one panel is exactly minus what enters the other.
//!
//! Each side computes F_out <- 0.5*(F_out - F_out_neighbour) from its own stored flux and
//! the buffer, so the two results are exact negatives of each other by construction and
//! no side has to own the face.

TaskStatus MeshBoundaryValuesCC::RecvAndUnpackFluxSeamCC(DvceFaceFld5D<Real> &flx) {
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &rbuf = recvbuf;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;

#if MPI_PARALLEL_ENABLED
  bool bflag = false;
  bool no_errors=true;
  bool x2f_h; Real sg_h;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ((nghbr.h_view(m,n).gid >= 0) &&
          SeamFaceGeom(recvbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                       cs_indcs.js, cs_indcs.ks, x2f_h, sg_h) &&
          (nghbr.h_view(m,n).panel != mbpanel.h_view(m)) &&
          (nghbr.h_view(m,n).rank != global_variable::my_rank)) {
        int test;
        int ierr = MPI_Test(&(rbuf[n].flux_req[m]), &test, MPI_STATUS_IGNORE);
        if (ierr != MPI_SUCCESS) {no_errors=false;}
        if (!(static_cast<bool>(test))) {bflag = true;}
      }
    }
  }
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (bflag) {return TaskStatus::incomplete;}
#endif

  int nvar = flx.x1f.extent_int(1);
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (nmb*nnghbr*nvar), Kokkos::AUTO);
  Kokkos::parallel_for("SeamFluxRecv", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    if (nghbr.d_view(m,n).gid < 0) {return;}
    bool x2face; Real sgn;
    if (!SeamFaceGeom(rbuf[n].iflux_same[0], cs_indcs.is, cs_indcs.ie,
                      cs_indcs.js, cs_indcs.ks, x2face, sgn)) {
      return;
    }
    if (nghbr.d_view(m,n).panel == mbpanel.d_view(m)) {return;}

    const int il = rbuf[n].iflux_same[0].bis;
    const int iu = rbuf[n].iflux_same[0].bie;
    const int jl = rbuf[n].iflux_same[0].bjs;
    const int ju = rbuf[n].iflux_same[0].bje;
    const int kl = rbuf[n].iflux_same[0].bks;
    const int ku = rbuf[n].iflux_same[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nkj = nk*nj;

    Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
      int k = idx / nj;
      int j = (idx - k * nj) + jl;
      k += kl;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1), [&](const int i) {
        const int index = i-il + ni*((j-jl) + nj*(k-kl) + nk*nj*v);
        // rbuf holds the neighbour's OUTWARD flux; ours is sgn*flx. Averaging the two
        // outward values and writing back gives 0.5*flx - 0.5*sgn*recv.
        const Real recv = rbuf[n].flux(m, index);
        if (x2face) {
          flx.x2f(m,v,k,j,i) = 0.5*flx.x2f(m,v,k,j,i) - 0.5*sgn*recv;
        } else {
          flx.x3f(m,v,k,j,i) = 0.5*flx.x3f(m,v,k,j,i) - 0.5*sgn*recv;
        }
      });
    });
  });
  return TaskStatus::complete;
}
