//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for face-centered (FC)
//! Mesh variables.
//! Prolongation of FC variables  occurs in ProlongateFC() function called from task list

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <iomanip>    // std::setprecision()

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"
#include "bvals/cs_seam_fc3.hpp"

//----------------------------------------------------------------------------------------
// BValFC constructor:

MeshBoundaryValuesFC::MeshBoundaryValuesFC(MeshBlockPack *pp, ParameterInput *pin) :
  MeshBoundaryValues(pp, pin, false) {
}

//----------------------------------------------------------------------------------------
//! \!fn void MeshBoundaryValuesFC::BuildCubedSphereSeamTables()
//! \brief Precompute everything the cubed-sphere seam packer derives from (xi,eta).
//!
//! x1 is RADIAL on the cubed sphere and no seam crosses it, so every geometric quantity
//! the seam branch of PackAndSendFC evaluates is a function of (kk,jj) only and is
//! recomputed identically for all nx1 cells of the radial loop.  Profiled on a 300-cycle
//! dhj MHD run that made "SendBuff" cost 1432 ms, a quarter of the GPU time, almost all
//! of it trigonometry: three TransformFieldToDstNormals per output cell (the resample
//! needs one per stencil cell), each with four PanelFrame switches, four tan/sqrt pairs,
//! a PanelToCart and a CartToPanel with two atan, plus the resample's own atan/tan and
//! the shear term's two sin.
//!
//! HOISTING IT INTO A NESTED LOOP -- which is what bvals_cc.cpp does, 8.9x there -- does
//! NOT work here: the face-centred pack needs THREE transforms per column rather than
//! one, and turning the flat (k,j,i) loop into TeamThreadRange(k,j) x
//! ThreadVectorRange(i) both serialises the radial direction (the TeamPolicy's vector
//! length is 1) and pushes the per-thread private segment from 908 to 2224 bytes.
//! Measured 1432 -> 8596 ms.
//! So the geometry goes into a TABLE instead and the loop structure is left alone.  The
//! table is indexed by (block, destination-panel slot, staggering kind, kk, jj); the
//! packer reads it with the (kk,jj) of each sample.  Consecutive threads of a team differ
//! only in i, so they all read the SAME table row and it broadcasts out of cache.
//!
//! BITWISE IDENTITY is the point: the fill below evaluates the SAME expressions, in the
//! same order, with the same inputs as the kernel used to inline, and the kernel then
//! calls ApplyFieldXform on the stored SeamFieldXform.  Nothing is refactored into a
//! smaller form (the 2x2 matrix ApplyFieldXform effectively builds would re-associate the
//! sums and change the last bits), so all 15 fields are stored.
//!
//! WHEN IT IS BUILT.  The table depends on the neighbour table, the panel of each block
//! and mb_size -- none of which change after startup on a cubed sphere: ADAPTIVE
//! refinement is refused outright (mesh.cpp), and static refinement is refused wherever a
//! coarse/fine interface would lie on a seam (build_tree.cpp's
//! CheckCubedSphereRefinement), so no regrid can alter it.  It is therefore built lazily
//! on the first PackAndSendFC call, once, and never invalidated.  Cross-level SEAM
//! buffers still exist on the diagonals, so the coarse staggering variants are built
//! whenever the mesh is multilevel.

void MeshBoundaryValuesFC::BuildCubedSphereSeamTables() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  // TWO KIND AXES.  The transform is needed only for the two ANGULAR components, the
  // resample for all three -- b.x1f crosses a seam as a scalar but is still resampled
  // along it -- so the geometry table carries 2 staggerings and the stencil table 3,
  // doubled when the mesh is multilevel and cross-level seam buffers exist.
  const int nkg = ml_ ? 4 : 2;
  const int nks = ml_ ? 6 : 3;
  // The tables are indexed by the ABSOLUTE (kk,jj) of the sample, GHOSTS INCLUDED.  On a
  // multilevel mesh a same-level send buffer is widened by ng in its transverse
  // directions -- that is what feeds the neighbour's prolongation -- so the packer reads
  // source cells outside the active range, and indexing the table from the active corner
  // ran it off the end.  One entry per face index of the block array covers everything.
  const int nk2 = indcs.nx2 + 2*indcs.ng + 1;
  const int nk3 = indcs.nx3 + 2*indcs.ng + 1;
  // fcorder3-1103: <mesh>/cs_seam_fc_order = 3 needs one level, ng >= 2 and >= 4 cells
  // per block side (3 cells / 4 faces next to the seam, 3-4 along it)
  if (pmy_pack->pmesh->cs_seam_fc_order == 3 &&
      (ml_ || indcs.ng < 2 || indcs.nx2 < 4 || indcs.nx3 < 4)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mesh>/cs_seam_fc_order = 3 needs a single-level cubed sphere, "
              << "nghost >= 2 and MeshBlocks of >= 4 cells in x2 and x3" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Which destination panels does each block actually send to?  At most four (the four
  // edge neighbours of a panel; the cube-vertex diagonals are skipped on both sides), so
  // a slot map keeps the table a quarter of the size of one indexed by panel.
  Kokkos::realloc(cs_seam_slot, nmb, 6);
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  int nslot = 1;
  for (int m=0; m<nmb; ++m) {
    for (int p=0; p<6; ++p) {
      cs_seam_slot.h_view(m,p) = -1;
    }
    int ns = 0;
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid < 0) continue;
      if (IsCubeVertexCorner(nghbr.h_view, mbpanel.h_view, m, n)) continue;
      const int dp = nghbr.h_view(m,n).panel;
      if (dp == mbpanel.h_view(m)) continue;
      if (cs_seam_slot.h_view(m,dp) < 0) {
        cs_seam_slot.h_view(m,dp) = ns;
        ++ns;
      }
    }
    if (ns > nslot) nslot = ns;
  }
  cs_seam_nslot = nslot;
  // modify_host()/sync_device(), not the explicit-space templates: those fail the
  // DualView static_assert under hipcc (same trap as cdd7d2a5, and the same fix the
  // radiation tetrads already carry).
  cs_seam_slot.modify_host();
  cs_seam_slot.sync_device();

  Kokkos::realloc(cs_seam_geom, nmb, nslot, nkg, nk3, nk2, 15);
  Kokkos::realloc(cs_seam_stnc, nmb, nks, 2, nk3, nk2, 7);
  // layoutbug-1103: the stencil of the second (extended) pass, <mesh>/cs_seam_blockend
  if (cs_blockend_ext_) {
    Kokkos::realloc(cs_seam_stnc_x, nmb, nks, 2, nk3, nk2, 7);
    Kokkos::realloc(cs_ext_jk, nmb, 4);
    auto hx = Kokkos::create_mirror_view(cs_ext_jk);
    auto &msz = pmy_pack->pmb->mb_size;
    for (int m=0; m<nmb; ++m) {
      const Real w2 = msz.h_view(m).x2max - msz.h_view(m).x2min;
      const Real w3 = msz.h_view(m).x3max - msz.h_view(m).x3min;
      const int nb2 = static_cast<int>(std::floor(2.0/w2 + 0.5));
      const int nb3 = static_cast<int>(std::floor(2.0/w3 + 0.5));
      const int ib2 = static_cast<int>(std::floor((msz.h_view(m).x2min + 1.0)/w2 + 0.5));
      const int ib3 = static_cast<int>(std::floor((msz.h_view(m).x3min + 1.0)/w3 + 0.5));
      hx(m,0) = (ib2 > 0) ? indcs.ng : 0;  hx(m,1) = (ib2 < nb2-1) ? indcs.ng : 0;
      hx(m,2) = (ib3 > 0) ? indcs.ng : 0;  hx(m,3) = (ib3 < nb3-1) ? indcs.ng : 0;
    }
    Kokkos::deep_copy(cs_ext_jk, hx);
    // seamexact-1103: with > 2 blocks per panel side a block end lies off the panel
    // midline, the along-seam resample shifts there and reaches the outermost along-seam
    // ghost, whose partner interpolation needs one face more than this block holds: those
    // interpolants come from the owning same-panel neighbour (CsSeamExtPartner).
    auto &pmi = pmy_pack->pmesh->mesh_indcs;
    cs_ext_q_ = (pmi.nx2 > 2*indcs.nx2) || (pmi.nx3 > 2*indcs.nx3);
    // fcorder3-1103: order 3 keeps every stencil within ng of the sample (srcval3), so
    // it needs no partner interpolants from the neighbour
    if (pmy_pack->pmesh->cs_seam_fc_order == 3) cs_ext_q_ = false;
    if (cs_ext_q_) {
      Kokkos::realloc(cs_ext_q, nmb, 4, indcs.nx3 + 2*indcs.ng, indcs.nx2 + 2*indcs.ng,
                      indcs.nx1 + 2*indcs.ng);
      Kokkos::deep_copy(cs_ext_q, 0.0);
      pbq_ = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
      pbq_->InitializeBuffers(4);
      pbq_->SetVectorPairs(4, {});     // plain scalars; only same-panel ghosts are read
      pbq_->cs_noresample = true;
    }
  }
  const bool bext_ = cs_blockend_ext_;
  auto stnx_ = cs_seam_stnc_x;
  const int ngx_ = indcs.ng;

  auto &mbsize = pmy_pack->pmb->mb_size;
  auto slot_ = cs_seam_slot;
  auto geom_ = cs_seam_geom;
  auto stnc_ = cs_seam_stnc;
  auto &cs_indcs = indcs;

  par_for("cs_seam_tbl", DevExeSpace(), 0, nmb-1, 0, nks-1, 0, nk3-1, 0, nk2-1,
  KOKKOS_LAMBDA(const int m, const int kd, const int kk, const int jj) {
    // kind = (coarse ? 3 : 0) + vv: the staggering of the SOURCE component vv the packer
    // reads, which is what fixes whether xi/eta are a face or a centre.  vv = 0 is a real
    // kind and not a spare: b.x1f crosses a seam as a SCALAR and needs no transform, but
    // it IS resampled along the seam, and it lives at BOTH cell centres.
    const int vv = kd % 3;
    const bool coar = (kd >= 3);
    const int js_ = coar ? cs_indcs.cjs : cs_indcs.js;
    const int ks_ = coar ? cs_indcs.cks : cs_indcs.ks;
    const int je_ = coar ? cs_indcs.cje : cs_indcs.je;
    const int ke_ = coar ? cs_indcs.cke : cs_indcs.ke;
    const int nx2_ = coar ? cs_indcs.cnx2 : cs_indcs.nx2;
    const int nx3_ = coar ? cs_indcs.cnx3 : cs_indcs.nx3;
    // the coarse array is half as wide; its rows sit in the same table, the rest unused
    if (coar && (jj > nx2_ + 2*cs_indcs.ng || kk > nx3_ + 2*cs_indcs.ng)) return;
    const int jq = jj - js_;
    const int kq = kk - ks_;

    const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
    const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
    // exactly the packer's `angles` lambda
    const Real xi  = 0.25*M_PI*((vv == 1) ? LeftEdgeX(jq, nx2_, x2mn, x2mx)
                                          : CellCenterX(jq, nx2_, x2mn, x2mx));
    const Real eta = 0.25*M_PI*((vv == 2) ? LeftEdgeX(kq, nx3_, x3mn, x3mx)
                                          : CellCenterX(kq, nx3_, x3mn, x3mx));

    // the field transform's geometry, one entry per destination panel this block uses
    const int psrc = mbpanel.d_view(m);
    if (vv > 0) {
      const int gkd = (coar ? 2 : 0) + (vv - 1);
      for (int p=0; p<6; ++p) {
        const int sl = slot_.d_view(m,p);
        if (sl < 0) continue;
        cubed_sphere::SeamFieldXform t;
        cubed_sphere::SeamFieldXformAt(psrc, p, xi, eta, t);
        for (int c=0; c<3; ++c) {
          geom_(m,sl,gkd,kk,jj,c)    = t.s1[c];
          geom_(m,sl,gkd,kk,jj,3+c)  = t.s2[c];
          geom_(m,sl,gkd,kk,jj,6+c)  = t.d1[c];
          geom_(m,sl,gkd,kk,jj,9+c)  = t.d2[c];
        }
        geom_(m,sl,gkd,kk,jj,12) = t.sns;
        geom_(m,sl,gkd,kk,jj,13) = t.csd;
        geom_(m,sl,gkd,kk,jj,14) = t.snd;
      }
    }

    // the along-seam resample's stencil and the shear term's sigma, for both seam
    // orientations (cs_seam = 2 stored at sidx 0, cs_seam = 3 at sidx 1).  Neither
    // depends on the destination panel.
    for (int sidx=0; sidx<2; ++sidx) {
      Real ang, nrm, dang;
      int sc, blo, bhi;
      if (sidx == 0) {
        ang = eta; nrm = xi;
        dang = 0.25*M_PI*(x3mx - x3mn)/static_cast<Real>(nx3_);
        sc = kk; blo = ks_; bhi = ke_ + ((vv == 2) ? 1 : 0) - 2;
      } else {
        ang = xi; nrm = eta;
        dang = 0.25*M_PI*(x2mx - x2mn)/static_cast<Real>(nx2_);
        sc = jj; blo = js_; bhi = je_ + ((vv == 1) ? 1 : 0) - 2;
      }
      const Real pos = sc + (atan(tan(ang)*tan(fabs(nrm))) - ang)/dang;
      int bs = static_cast<int>(floor(pos + 0.5)) - 1;
      bs = (bs < blo) ? blo : ((bs > bhi) ? bhi : bs);
      const Real u = pos - static_cast<Real>(bs + 1);
      stnc_(m,kd,sidx,kk,jj,0) = static_cast<Real>(bs);
      stnc_(m,kd,sidx,kk,jj,1) = 0.5*u*(u - 1.0);
      stnc_(m,kd,sidx,kk,jj,2) = 1.0 - u*u;
      stnc_(m,kd,sidx,kk,jj,3) = 0.5*u*(u + 1.0);
      stnc_(m,kd,sidx,kk,jj,4) = pos - static_cast<Real>(bs);
      // SHEAR CORRECTION geometry.  `along_eta` is (cs_seam == 2), i.e. sidx == 0, and
      // its (ang_,nrm_) are exactly the (ang,nrm) the resample already picked, so one
      // pair serves both.  sigma is stored unconditionally; the packer keeps the
      // |sin(2 nrm)| > 1e-8 guard, so a division by zero here is never read back.
      const Real sden = sin(2.0*nrm);
      stnc_(m,kd,sidx,kk,jj,5) = sden;
      stnc_(m,kd,sidx,kk,jj,6) = sin(2.0*ang)/sden;
      if (bext_) {
        // the same stencil clamped to the PANEL (up to ng cells into this block's own
        // along-seam ghosts at an interior block end); bitwise the entry above elsewhere
        int blx = blo, bhx = bhi;
        if (!coar) {
          const Real w_ = (sidx == 0) ? (x3mx - x3mn) : (x2mx - x2mn);
          const Real x0_ = (sidx == 0) ? x3mn : x2mn;
          const int nbp_ = static_cast<int>(floor(2.0/w_ + 0.5));
          const int ib_ = static_cast<int>(floor((x0_ + 1.0)/w_ + 0.5));
          if (ib_ > 0) {blx -= ngx_;}
          if (ib_ < nbp_ - 1) {bhx += ngx_;}
        }
        int bx = static_cast<int>(floor(pos + 0.5)) - 1;
        bx = (bx < blx) ? blx : ((bx > bhx) ? bhx : bx);
        const Real ux = pos - static_cast<Real>(bx + 1);
        stnx_(m,kd,sidx,kk,jj,0) = static_cast<Real>(bx);
        stnx_(m,kd,sidx,kk,jj,1) = 0.5*ux*(ux - 1.0);
        stnx_(m,kd,sidx,kk,jj,2) = 1.0 - ux*ux;
        stnx_(m,kd,sidx,kk,jj,3) = 0.5*ux*(ux + 1.0);
        stnx_(m,kd,sidx,kk,jj,4) = pos - static_cast<Real>(bx);
        stnx_(m,kd,sidx,kk,jj,5) = sden;
        stnx_(m,kd,sidx,kk,jj,6) = sin(2.0*ang)/sden;
      }
    }
  });

  cs_seam_tbl_ready = true;
  return;
}

//----------------------------------------------------------------------------------------
//! \!fn void MeshBoundaryValuesFC::BuildPackSegmentsFC()
//! \brief halofc-1102: resolve the per-(m,n,v) setup of the cubed-sphere face-field pack
//! ONCE.
//!
//! WHERE THE TIME WENT (r1m, 1 H200, developed state, ncu on the old team kernel): 1.36-
//! 1.50 ms per call, 4032 teams of 128 threads, 128 registers and a 2720-byte STACK per
//! thread.  The nested by-reference lambdas of the seam path (seamval -> srcval ->
//! x3f_at_etacentre -> x3f_at_xiface, each capturing the team's locals by reference)
//! kept every setup variable addressable, so the kernel issued 15.1 M local loads and
//! 6.1 M local stores (60 M + 24 M sectors, the stores mostly L1 misses) against 1.9 M
//! global load requests; achieved occupancy 13 % (theoretical 25 %, register bound),
//! SMs active 47 % of the elapsed cycles (a few teams carry ~3600 seam samples each,
//! most of the 4032 are empty).  Nothing of it is bandwidth: 50 MB of DRAM per call.
//!
//! So: (1) the setup -- which depends only on the neighbour table, the block levels and
//! panels and the buffer index ranges, none of which change on a cubed sphere (adaptive
//! refinement is refused, see BuildCubedSphereSeamTables) -- is evaluated here by the
//! SAME code the team kernel ran, into one PackSegFC per active triple; (2) the samples
//! are walked flat, one thread each, in chunks of FC_PACK_CHUNK samples of one triple, so
//! the load balances and no thread repeats any setup; (3) the seam arithmetic is a
//! functor with by-value arguments (FCSeamPack below), every expression unchanged in
//! form and order, so the buffers are bitwise those of the team kernel.

namespace {
constexpr int FC_PACK_CHUNK = 128;
constexpr int kFc3Max = 48;   // fcorder3-1103: max distinct faces of one order-3 sample
}  // namespace

void MeshBoundaryValuesFC::BuildPackSegmentsFC() {
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  const int nmnv = 3*nmb*nnghbr;
  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool use_pole = pmy_pack->pmesh->use_polar_boundary;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto cs_sltb = cs_seam_slot;
  auto sbuf = SendBufDv();

  DvceArray1D<PackSegFC> all("fc_seg_all", nmnv);
  par_for("fc_seg_build", DevExeSpace(), 0, nmnv-1, KOKKOS_LAMBDA(const int lr) {
    const int m = lr/(3*nnghbr);
    const int n = (lr/3)%nnghbr;
    const int v = lr%3;
    PackSegFC s;
    s.m = m; s.n = n; s.v = v; s.vv = v;
    s.il = s.jl = s.kl = 0; s.ni = s.nj = s.nk = 0; s.ndat = 0;
    s.dm = s.dn = 0; s.local = 0; s.mode = -1;
    s.aj = 1; s.bj = 0; s.ak = 1; s.bk = 0; s.sj = 1; s.sk = 0; s.signvar = 1;
    s.cs_xform = 0; s.cs_coar = 0; s.cs_slotidx = 0; s.cs_gkind = 0; s.cs_skind = 0;
    s.cs_seam = 0;
    s.cs_psrc = mbpanel.d_view(m); s.cs_pdst = s.cs_psrc;
    if (nghbr.d_view(m,n).gid >= 0 &&
      !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
      int il, iu, jl, ju, kl, ku, ndat;
      int vsrc = v;
      if (use_cs &&
          nghbr.d_view(m,n).panel != mbpanel.d_view(m) && v > 0) {
        PanelBoundaries pb0;
        pb0 = GetPanelBoundary(mbpanel.d_view(m),
                                                nghbr.d_view(m,n).panel);
        if (pb0.swap_ax == 1) vsrc = (v == 1) ? 2 : 1;
      }
      if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
        il = sbuf[n].icoar[vsrc].bis;
        iu = sbuf[n].icoar[vsrc].bie;
        jl = sbuf[n].icoar[vsrc].bjs;
        ju = sbuf[n].icoar[vsrc].bje;
        kl = sbuf[n].icoar[vsrc].bks;
        ku = sbuf[n].icoar[vsrc].bke;
        ndat = sbuf[n].icoar_ndat;
      } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
        il = sbuf[n].isame[vsrc].bis;
        iu = sbuf[n].isame[vsrc].bie;
        jl = sbuf[n].isame[vsrc].bjs;
        ju = sbuf[n].isame[vsrc].bje;
        kl = sbuf[n].isame[vsrc].bks;
        ku = sbuf[n].isame[vsrc].bke;
        ndat = sbuf[n].isame_ndat;
      } else {
        il = sbuf[n].ifine[vsrc].bis;
        iu = sbuf[n].ifine[vsrc].bie;
        jl = sbuf[n].ifine[vsrc].bjs;
        ju = sbuf[n].ifine[vsrc].bje;
        kl = sbuf[n].ifine[vsrc].bks;
        ku = sbuf[n].ifine[vsrc].bke;
        ndat = sbuf[n].ifine_ndat;
      }
      const int ni = iu - il + 1;
      const int nj = ju - jl + 1;
      const int nk = ku - kl + 1;
      s.il = il; s.jl = jl; s.kl = kl; s.ni = ni; s.nj = nj; s.nk = nk; s.ndat = ndat;
      s.dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
      s.dn = nghbr.d_view(m,n).dest;
      s.local = (nghbr.d_view(m,n).rank == my_rank) ? 1 : 0;

      const bool do_cs = use_cs &&
                         (nghbr.d_view(m,n).panel != mbpanel.d_view(m));
      const bool do_pole = use_pole &&
                           (nghbr.d_view(m,n).polar > 0);
      if (do_cs || do_pole) {
        int aj = 1, bj = 0;
        int ak = 1, bk = 0;
        int signvar = 1;
        int sj = 1, sk = nj;
        int vv = v;
        bool cs_xform = false;
        const bool cs_coar = (nghbr.d_view(m,n).lev < mblev.d_view(m));
        int cs_slotidx = 0, cs_gkind = 0, cs_skind = 0;
        int cs_seam = 0;
        if (do_cs) {
          const auto ngh = nghbr.d_view(m,n);
          const int my_panel = mbpanel.d_view(m);
          PanelBoundaries pb;
          pb = GetPanelBoundary(my_panel, ngh.panel);
          int map_vy = 1;
          int map_vz = 2;
          int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
          int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
          if (pb.swap_ax == 1) {
            map_vy = 2;
            map_vz = 1;
            sj = nk;
            sk = 1;
          }
          if (rev_a_preswap) {
            aj = -1;
            bj = jl + ju;
          }
          if (rev_b_preswap) {
            ak = -1;
            bk = kl + ku;
          }
          if (v == 1) vv = map_vy;
          if (v == 2) vv = map_vz;
          cs_xform = (v == 1) || (v == 2);
          cs_slotidx = cs_sltb.d_view(m, ngh.panel);
          cs_gkind = (cs_coar ? 2 : 0) + ((vv > 0) ? (vv - 1) : 0);
          cs_skind = (cs_coar ? 3 : 0) + vv;
          if (n >= 8 && n < 24) {
            cs_seam = 2;
          } else if (n >= 24 && n < 40) {
            cs_seam = 3;
          }
          if (n >= 40 && n < 56) {
            int sjc, skc;
            if (n < 48) {
              const int cq = (n - 40)/2;
              sjc = cq & 1;  skc = (cq >> 1) & 1;
            } else {
              const int cq = n - 48;
              sjc = (cq >> 1) & 1;  skc = (cq >> 2) & 1;
            }
            const int njid = sjc ? 12 : 8;
            const int nkid = skc ? 28 : 24;
            const bool sj_ = (nghbr.d_view(m,njid).gid >= 0 &&
                              nghbr.d_view(m,njid).panel != my_panel);
            const bool sk_ = (nghbr.d_view(m,nkid).gid >= 0 &&
                              nghbr.d_view(m,nkid).panel != my_panel);
            if (sj_ && !sk_) {
              cs_seam = 2;
            } else if (sk_ && !sj_) {
              cs_seam = 3;
            }
          }
          const int seam_extent = cs_coar
              ? ((cs_seam == 2) ? cs_indcs.cnx3 : cs_indcs.cnx2)
              : ((cs_seam == 2) ? cs_indcs.nx3 : cs_indcs.nx2);
          if (cs_seam != 0 && seam_extent < 3) { cs_seam = 0; }
        } else if (do_pole) {
          aj = -1;
          bj = jl + ju;
          if (v == 1) signvar = -1;
          if (v == 2) signvar = -1;
        }
        s.mode = 2;
        s.aj = aj; s.bj = bj; s.ak = ak; s.bk = bk; s.sj = sj; s.sk = sk;
        s.signvar = signvar; s.vv = vv;
        s.cs_xform = cs_xform ? 1 : 0; s.cs_coar = cs_coar ? 1 : 0;
        s.cs_slotidx = cs_slotidx; s.cs_gkind = cs_gkind; s.cs_skind = cs_skind;
        s.cs_seam = cs_seam;
        if (do_cs) s.cs_pdst = nghbr.d_view(m,n).panel;
      } else {
        s.mode = (nghbr.d_view(m,n).lev >= mblev.d_view(m)) ? 0 : 1;
      }
    }
    all(lr) = s;
  });

  // compact on the host: plain triples first, then the seam/pole ones, each cut into
  // chunks of FC_PACK_CHUNK samples
  auto h_all = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), all);
  std::vector<int> act;
  for (int pass=0; pass<2; ++pass) {
    for (int lr=0; lr<nmnv; ++lr) {
      const int md = h_all(lr).mode;
      if (md < 0) continue;
      if ((pass == 0) != (md != 2)) continue;
      act.push_back(lr);
    }
  }
  const int nseg = static_cast<int>(act.size());
  Kokkos::realloc(fc_seg, (nseg > 0) ? nseg : 1);
  auto h_seg = Kokkos::create_mirror_view(fc_seg);
  std::vector<int> cseg, cbeg;
  int nplain = 0;
  for (int q=0; q<nseg; ++q) {
    h_seg(q) = h_all(act[q]);
    const int nkji = h_seg(q).nk*h_seg(q).nj*h_seg(q).ni;
    for (int c0=0; c0<nkji; c0+=FC_PACK_CHUNK) {
      cseg.push_back(q);
      cbeg.push_back(c0);
      if (h_seg(q).mode != 2) ++nplain;
    }
  }
  fc_nchunk_plain = nplain;
  fc_nchunk_seam = static_cast<int>(cseg.size()) - nplain;
  const int nch = static_cast<int>(cseg.size());
  Kokkos::realloc(fc_chunk_seg, (nch > 0) ? nch : 1);
  Kokkos::realloc(fc_chunk_beg, (nch > 0) ? nch : 1);
  auto h_cs = Kokkos::create_mirror_view(fc_chunk_seg);
  auto h_cb = Kokkos::create_mirror_view(fc_chunk_beg);
  for (int c=0; c<nch; ++c) { h_cs(c) = cseg[c]; h_cb(c) = cbeg[c]; }
  Kokkos::deep_copy(fc_seg, h_seg);
  Kokkos::deep_copy(fc_chunk_seg, h_cs);
  Kokkos::deep_copy(fc_chunk_beg, h_cb);

  // r1mspeed-1102: the two walk sets (all layers / depth 1) with the seam srcval boxes
  {
    const int js_ = cs_indcs.js, je_ = cs_indcs.je, ks_ = cs_indcs.ks, ke_ = cs_indcs.ke;
    // fcorder3-1103: order 3 evaluates every seam sample itself (srcval3, no resample box)
    const bool fc3_ = (pmy_pack->pmesh->cs_seam_fc_order == 3);
    int64_t npv_max = 0;
    for (int d=0; d<2; ++d) {
      std::vector<SubSegFC> sub;
      std::vector<int> cs_p, cb_p, cs_s, cb_s, cs_r, cb_r;
      int64_t pv = 0;
      for (int q=0; q<nseg; ++q) {
        const PackSegFC &g = h_seg(q);
        SubSegFC u;
        u.q = q; u.il = g.il; u.ni = g.ni; u.jl = g.jl; u.nj = g.nj; u.kl = g.kl;
        u.nk = g.nk;
        u.pvoff = -1; u.ptlo = 0; u.pnt = 0; u.pnlo = 0; u.pnn = 0; u.pil = 0; u.pni = 0;
        if (d == 1) {
          // direction of buffer n (nghbr_index.hpp: 0-7 x1 faces, 8-15 x2 faces, 16-23
          // x1x2 edges, 24-31 x3 faces, 32-39 x3x1 edges, 40-47 x2x3 edges, 48-55
          // corners); in each direction the buffer is ghost in, keep the source layer
          // whose READ index (k,j after the seam map) is adjacent to the boundary
          const int n = g.n;
          int o1 = 0, o2 = 0, o3 = 0;
          if (n < 8) {
            o1 = (n < 4) ? -1 : 1;
          } else if (n < 16) {
            o2 = (n < 12) ? -1 : 1;
          } else if (n < 24) {
            o1 = (((n-16) >> 1) & 1) ? 1 : -1;  o2 = (((n-16) >> 2) & 1) ? 1 : -1;
          } else if (n < 32) {
            o3 = (n < 28) ? -1 : 1;
          } else if (n < 40) {
            o1 = (((n-32) >> 1) & 1) ? 1 : -1;  o3 = (((n-32) >> 2) & 1) ? 1 : -1;
          } else if (n < 48) {
            o2 = (((n-40) >> 1) & 1) ? 1 : -1;  o3 = (((n-40) >> 2) & 1) ? 1 : -1;
          } else {
            o1 = ((n-48) & 1) ? 1 : -1;  o2 = (((n-48) >> 1) & 1) ? 1 : -1;
            o3 = (((n-48) >> 2) & 1) ? 1 : -1;
          }
          if (o1 != 0) {u.il = (o1 < 0) ? g.il : g.il + g.ni - 1; u.ni = 1;}
          if (o2 != 0) {
            const int jt = (o2 < 0) ? g.jl : g.jl + g.nj - 1;
            u.jl = (g.aj == 1) ? jt : (g.bj - jt); u.nj = 1;
          }
          if (o3 != 0) {
            const int kt = (o3 < 0) ? g.kl : g.kl + g.nk - 1;
            u.kl = (g.ak == 1) ? kt : (g.bk - kt); u.nk = 1;
          }
        }
        if (g.mode == 2 && g.cs_seam != 0 && !g.cs_coar && !fc3_) {
          // the samples' resample reads srcval at (bs..bs+2) along the seam, bs within
          // [lo, hi-2] (BuildCubedSphereSeamTables), at the sample's own normal index
          u.pil = u.il; u.pni = u.ni;
          // layoutbug-1103: with cs_seam_blockend = extend the box also covers ng cells
          // on each side (the second pass's stencil reaches the along-seam ghosts)
          const int bxg_ = cs_blockend_ext_ ? cs_indcs.ng : 0;
          if (g.cs_seam == 2) {   // along eta (k); normal index jj = aj*j + bj
            u.ptlo = ks_ - bxg_; u.pnt = ke_ + ((g.vv == 2) ? 1 : 0) - ks_ + 1 + 2*bxg_;
            const int a0 = g.aj*u.jl + g.bj, a1 = g.aj*(u.jl + u.nj - 1) + g.bj;
            u.pnlo = (a0 < a1) ? a0 : a1; u.pnn = u.nj;
          } else {                 // along xi (j); normal index kk = ak*k + bk
            u.ptlo = js_ - bxg_; u.pnt = je_ + ((g.vv == 1) ? 1 : 0) - js_ + 1 + 2*bxg_;
            const int a0 = g.ak*u.kl + g.bk, a1 = g.ak*(u.kl + u.nk - 1) + g.bk;
            u.pnlo = (a0 < a1) ? a0 : a1; u.pnn = u.nk;
          }
          u.pvoff = static_cast<int>(pv);
          const int npts = u.pnn*u.pnt*u.pni;
          const int us = static_cast<int>(sub.size());
          for (int c0=0; c0<npts; c0+=FC_PACK_CHUNK) {cs_r.push_back(us); cb_r.push_back(c0);}
          pv += npts;
        } else if (g.mode == 2) {
          // every other seam / pole segment: srcval at each sample itself, indexed by the
          // sample (pnt = 0), so the sample pass never evaluates srcval
          u.pvoff = static_cast<int>(pv); u.pnt = 0;
          const int npts = u.nk*u.nj*u.ni;
          const int us = static_cast<int>(sub.size());
          for (int c0=0; c0<npts; c0+=FC_PACK_CHUNK) {cs_r.push_back(us); cb_r.push_back(c0);}
          pv += npts;
        }
        const int us = static_cast<int>(sub.size());
        const int nkji = u.nk*u.nj*u.ni;
        for (int c0=0; c0<nkji; c0+=FC_PACK_CHUNK) {
          if (g.mode != 2) {cs_p.push_back(us); cb_p.push_back(c0);}
          else {cs_s.push_back(us); cb_s.push_back(c0);}
        }
        sub.push_back(u);
      }
      npv_max = (pv > npv_max) ? pv : npv_max;
      fc_np[d] = static_cast<int>(cs_p.size());
      fc_ns[d] = static_cast<int>(cs_s.size());
      fc_npre[d] = static_cast<int>(cs_r.size());
      // chunk lists: plain, then seam, then precompute
      std::vector<int> ca(cs_p), cb(cb_p);
      ca.insert(ca.end(), cs_s.begin(), cs_s.end()); cb.insert(cb.end(), cb_s.begin(), cb_s.end());
      ca.insert(ca.end(), cs_r.begin(), cs_r.end()); cb.insert(cb.end(), cb_r.begin(), cb_r.end());
      const int nsub = static_cast<int>(sub.size());
      Kokkos::realloc(fc_sub[d], (nsub > 0) ? nsub : 1);
      auto h_sub = Kokkos::create_mirror_view(fc_sub[d]);
      for (int u=0; u<nsub; ++u) {h_sub(u) = sub[u];}
      Kokkos::deep_copy(fc_sub[d], h_sub);
      const int nc = static_cast<int>(ca.size());
      Kokkos::realloc(fc_cs[d], (nc > 0) ? nc : 1);
      Kokkos::realloc(fc_cb[d], (nc > 0) ? nc : 1);
      auto h_a = Kokkos::create_mirror_view(fc_cs[d]);
      auto h_b = Kokkos::create_mirror_view(fc_cb[d]);
      for (int c=0; c<nc; ++c) {h_a(c) = ca[c]; h_b(c) = cb[c];}
      Kokkos::deep_copy(fc_cs[d], h_a);
      Kokkos::deep_copy(fc_cb[d], h_b);
    }
    Kokkos::realloc(fc_preval, (npv_max > 0) ? npv_max : 1);
  }
  fc_seg_ready = true;
}

namespace {

//----------------------------------------------------------------------------------------
//! \struct FCPlainPack
//! \brief halofc-1102: one thread per sample of a plain (same-panel) face-field buffer.

struct FCPlainPack {
  DvceFaceFld4D<Real> b, cb;
  DvceArray1D<MeshBufferDv> sbuf, rbuf;
  DvceArray1D<PackSegFC> seg;
  DvceArray1D<SubSegFC> sub;   // r1mspeed-1102: chunk -> sample box (all / depth 1)
  DvceArray1D<int> cseg, cbeg;
  int c0;   // first chunk of this launch

  FCPlainPack(DvceFaceFld4D<Real> &b_, DvceFaceFld4D<Real> &cb_,
              DvceArray1D<MeshBufferDv> sb_, DvceArray1D<MeshBufferDv> rb_,
              DvceArray1D<PackSegFC> sg_, DvceArray1D<SubSegFC> su_,
              DvceArray1D<int> cs_, DvceArray1D<int> cb2_,
              int c0_) : b(b_), cb(cb_), sbuf(sb_), rbuf(rb_), seg(sg_), sub(su_),
              cseg(cs_), cbeg(cb2_), c0(c0_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int t) const {
    const int c = c0 + t/FC_PACK_CHUNK;
    const SubSegFC u = sub(cseg(c));
    const int q = u.q;
    const int idx = cbeg(c) + t%FC_PACK_CHUNK;
    if (idx >= u.nk*u.nj*u.ni) return;
    const int ni = seg(q).ni, nj = seg(q).nj;
    const int il = seg(q).il, jl = seg(q).jl, kl = seg(q).kl;
    const int m = seg(q).m, v = seg(q).v, ndat = seg(q).ndat;
    const int uji = u.nj*u.ni;
    int k = (idx)/uji;
    int j = (idx - k*uji)/u.ni;
    int i = (idx - k*uji - j*u.ni) + u.il;
    k += u.kl;
    j += u.jl;
    Real val;
    if (seg(q).mode == 0) {
      val = (v == 0) ? b.x1f(m,k,j,i) : ((v == 1) ? b.x2f(m,k,j,i) : b.x3f(m,k,j,i));
    } else {
      val = (v == 0) ? cb.x1f(m,k,j,i) : ((v == 1) ? cb.x2f(m,k,j,i) : cb.x3f(m,k,j,i));
    }
    const int ix = ndat*v + i-il + ni*(j-jl + nj*(k-kl));
    if (seg(q).local) {
      rbuf(seg(q).dn).vars(seg(q).dm, ix) = val;
    } else {
      sbuf(seg(q).n).vars(m, ix) = val;
    }
  }
};

//----------------------------------------------------------------------------------------
//! \struct FCSeamPack
//! \brief halofc-1102: one thread per sample of a cross-panel (or polar) face-field
//! buffer.  The member functions are the old kernel's lambdas (bx1..bx3, lag4, the four
//! partner interpolations, srcval, seamval) with every captured setup variable passed BY
//! VALUE in a SeamArgs; each expression is copied unchanged in form and order.  See the
//! long comments at the old kernel (PackAndSendFC, non-cubed-sphere path) for the
//! physics of every step.

struct SeamArgs {
  int m, v, vv, signvar, cs_xform, cs_coar, cs_slotidx, cs_gkind, cs_skind, cs_seam;
  int js_, ks_, je_, ke_;
  bool wide2, wide3;
  int pvoff, ptlo, pnt, pnlo, pil, pni;   // r1mspeed-1102 (SubSegFC); pvoff < 0: none
};

struct FCSeamPre {};   // r1mspeed-1102: tag of the srcval precompute pass
struct FCSeamQ {};     // seamexact-1103: tag of the partner-interpolant pass
struct FCSeamPre3 {};  // fcorder3-1103: tag of the order-3 precompute pass
struct FCSeamTab {};   // fcorder3-1103: tag of the order-3 coefficient-table build

// fcorder3-1103: the order-3 coefficient tables, one row per (seam segment, sample k,j),
// [0] with the first pass's bounds, [1] with the extended pass's (cs_seam_blockend =
// extend).  Kept here, keyed by the boundary object, so bvals.hpp is untouched.
struct Fc3Tables {
  bool ready[2] = {false, false};
  int nrow = 0;
  DvceArray1D<int> gbase, rq;
  DvceArray1D<int> cnt[2];
  DvceArray2D<int> code[2];
  DvceArray2D<Real> w[2];
};
// never destroyed: its Views must not outlive Kokkos::finalize
std::map<const void*, Fc3Tables> *fc3_tables_ = new std::map<const void*, Fc3Tables>;

struct FCSeamPack {
  DvceFaceFld4D<Real> b, cb;
  DvceArray1D<MeshBufferDv> sbuf, rbuf;
  DvceArray1D<PackSegFC> seg;
  DvceArray1D<SubSegFC> sub;
  DvceArray1D<int> cseg, cbeg;
  DvceArray6D<Real> cs_gtbl, cs_stbl;
  DvceArray1D<Real> pre;   // r1mspeed-1102: srcval scratch (SubSegFC pv* boxes)
  int c0;
  RegionIndcs cs_indcs;
  Real cs_shear_sgn;
  bool cs_nolim_;
  // layoutbug-1103: second pass of <mesh>/cs_seam_blockend = extend -- the source index
  // bounds are widened by (j lo, j hi, k lo, k hi) = ext_jk(m,:) at interior block ends
  bool ext_ = false;
  DvceArray2D<int> ext_jk;
  // seamexact-1103: the extended pass reads the partner interpolants at the along-seam
  // ghost positions of an interior block end from q (CsSeamExtPartner) instead of
  // evaluating them on this block's own faces (one-sided at the outermost ghost)
  bool useq_ = false;
  DvceArray5D<Real> q;
  // fcorder3-1103 <mesh>/cs_seam_fc_order = 3: srcval3 (block sizes for the angles)
  bool fc3_ = false;
  DualArray1D<RegionSize>::t_dev mbsz;
  bool fc3tab_ = false;            // evaluate from the coefficient table
  // fcorder3-1103 DIAGNOSTIC (env FC3_OVS): ghost-face OVERSHOOT beyond the range its own
  // stencil data allows (order 2: the 3 resample nodes; order 3: per source component,
  // weight sum x [min, max] of the faces read).  ovsv: max abs, max relative to the
  // stencil's max |face|, count relative > 1e-2, samples
  bool ovs_ = false;
  DvceArray1D<Real> ovsv;
  KOKKOS_INLINE_FUNCTION
  void ovs_note(const Real val, const Real lo, const Real hi, const Real sc) const {
    const Real o = fmax(0.0, fmax(val - hi, lo - val));
    Kokkos::atomic_max(&ovsv(0), o);
    const Real r = o/fmax(sc, 1.0e-300);
    Kokkos::atomic_max(&ovsv(1), r);
    if (r > 1.0e-2) Kokkos::atomic_add(&ovsv(2), 1.0);
    Kokkos::atomic_add(&ovsv(3), 1.0);
  }
  DvceArray1D<int> t_gbase, t_rq, t_cnt;
  DvceArray2D<int> t_code;
  DvceArray2D<Real> t_w;

  FCSeamPack(DvceFaceFld4D<Real> &b_, DvceFaceFld4D<Real> &cb_,
             DvceArray1D<MeshBufferDv> sb_, DvceArray1D<MeshBufferDv> rb_,
             DvceArray1D<PackSegFC> sg_, DvceArray1D<SubSegFC> su_,
             DvceArray1D<int> cs_, DvceArray1D<int> cb2_,
             DvceArray6D<Real> gt_, DvceArray6D<Real> st_, DvceArray1D<Real> pr_,
             int c0_, const RegionIndcs &ix_, Real sh_, bool nl_) :
             b(b_), cb(cb_), sbuf(sb_), rbuf(rb_), seg(sg_), sub(su_), cseg(cs_),
             cbeg(cb2_), cs_gtbl(gt_), cs_stbl(st_), pre(pr_), c0(c0_), cs_indcs(ix_),
             cs_shear_sgn(sh_), cs_nolim_(nl_) {}

  KOKKOS_FORCEINLINE_FUNCTION
  Real bx1(const SeamArgs &a, const int kk, const int jj, const int i) const {
    return a.cs_coar ? cb.x1f(a.m,kk,jj,i) : b.x1f(a.m,kk,jj,i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real bx2(const SeamArgs &a, const int kk, const int jj, const int i) const {
    return a.cs_coar ? cb.x2f(a.m,kk,jj,i) : b.x2f(a.m,kk,jj,i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real bx3(const SeamArgs &a, const int kk, const int jj, const int i) const {
    return a.cs_coar ? cb.x3f(a.m,kk,jj,i) : b.x3f(a.m,kk,jj,i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  static Real lag4(const Real t, const Real f0, const Real f1,
                   const Real f2, const Real f3) {
    return -(t-1.0)*(t-2.0)*(t-3.0)/6.0*f0
           + t*(t-2.0)*(t-3.0)/2.0*f1
           - t*(t-1.0)*(t-3.0)/2.0*f2
           + t*(t-1.0)*(t-2.0)/6.0*f3;
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real x3f_at_xiface(const SeamArgs &a, const int kf, const int jf, const int i) const {
    const int js_ = a.js_, je_ = a.je_;
    if (!a.wide2) {
      if (jf-1 < js_) { return 1.5*bx3(a,kf,js_,i) - 0.5*bx3(a,kf,js_+1,i); }
      if (jf   > je_) { return 1.5*bx3(a,kf,je_,i) - 0.5*bx3(a,kf,je_-1,i); }
      return 0.5*(bx3(a,kf,jf-1,i) + bx3(a,kf,jf,i));
    }
    int w = jf-2;
    w = (w < js_) ? js_ : ((w > je_-3) ? je_-3 : w);
    return lag4(static_cast<Real>(jf) - 0.5 - static_cast<Real>(w),
                bx3(a,kf,w,i), bx3(a,kf,w+1,i), bx3(a,kf,w+2,i), bx3(a,kf,w+3,i));
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real x2f_at_etaface(const SeamArgs &a, const int kf, const int jf, const int i) const {
    const int ks_ = a.ks_, ke_ = a.ke_;
    if (!a.wide3) {
      if (kf-1 < ks_) { return 1.5*bx2(a,ks_,jf,i) - 0.5*bx2(a,ks_+1,jf,i); }
      if (kf   > ke_) { return 1.5*bx2(a,ke_,jf,i) - 0.5*bx2(a,ke_-1,jf,i); }
      return 0.5*(bx2(a,kf-1,jf,i) + bx2(a,kf,jf,i));
    }
    int w = kf-2;
    w = (w < ks_) ? ks_ : ((w > ke_-3) ? ke_-3 : w);
    return lag4(static_cast<Real>(kf) - 0.5 - static_cast<Real>(w),
                bx2(a,w,jf,i), bx2(a,w+1,jf,i), bx2(a,w+2,jf,i), bx2(a,w+3,jf,i));
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real x3f_at_etacentre(const SeamArgs &a, const int kk, const int jj,
                        const int i) const {
    const int ks_ = a.ks_, ke_ = a.ke_;
    if (!a.wide3) {
      return 0.5*(x3f_at_xiface(a,kk,jj,i) + x3f_at_xiface(a,kk+1,jj,i));
    }
    int w = kk-1;
    w = (w < ks_) ? ks_ : ((w > ke_-2) ? ke_-2 : w);
    return lag4(static_cast<Real>(kk) + 0.5 - static_cast<Real>(w),
                x3f_at_xiface(a,w,jj,i),   x3f_at_xiface(a,w+1,jj,i),
                x3f_at_xiface(a,w+2,jj,i), x3f_at_xiface(a,w+3,jj,i));
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real x2f_at_xicentre(const SeamArgs &a, const int kk, const int jj,
                       const int i) const {
    const int js_ = a.js_, je_ = a.je_;
    if (!a.wide2) {
      return 0.5*(x2f_at_etaface(a,kk,jj,i) + x2f_at_etaface(a,kk,jj+1,i));
    }
    int w = jj-1;
    w = (w < js_) ? js_ : ((w > je_-2) ? je_-2 : w);
    return lag4(static_cast<Real>(jj) + 0.5 - static_cast<Real>(w),
                x2f_at_etaface(a,kk,w,i),   x2f_at_etaface(a,kk,w+1,i),
                x2f_at_etaface(a,kk,w+2,i), x2f_at_etaface(a,kk,w+3,i));
  }
  // seamexact-1103: the partner component at (kk,jj) -- x3f at an x2-face centre (x2 =
  // true) or x2f at an x3-face centre.  In the extended pass, at an along-seam ghost
  // position beyond an INTERIOR block end, the owning neighbour's value (q, see
  // CsSeamExtPartner); everywhere else this block's own interpolation, as before.
  KOKKOS_FORCEINLINE_FUNCTION
  Real q_partner(const SeamArgs &a, const int kk, const int jj, const int i,
                 const bool x2) const {
    // (radial ghosts i excluded: the same-panel exchange of q carries active i only, and
    // those seam ghosts are rewritten by the radial boundary condition anyway)
    if (useq_ && !a.cs_coar && i >= cs_indcs.is && i <= cs_indcs.ie) {
      const int m = a.m;
      if (a.cs_seam == 2) {          // along eta: kk is the along index
        const int ke1 = cs_indcs.ke + (x2 ? 0 : 1);
        if (kk < cs_indcs.ks && ext_jk(m,2) > 0) return q(m, x2 ? 0 : 2, kk, jj, i);
        if (kk > ke1 && ext_jk(m,3) > 0) {
          return x2 ? q(m,0,kk,jj,i) : q(m,3,kk-1,jj,i);
        }
      } else if (a.cs_seam == 3) {   // along xi: jj is the along index
        const int je1 = cs_indcs.je + (x2 ? 1 : 0);
        if (jj < cs_indcs.js && ext_jk(m,0) > 0) return q(m, x2 ? 0 : 2, kk, jj, i);
        if (jj > je1 && ext_jk(m,1) > 0) {
          return x2 ? q(m,1,kk,jj-1,i) : q(m,2,kk,jj,i);
        }
      }
    }
    return x2 ? x3f_at_etacentre(a,kk,jj,i) : x2f_at_xicentre(a,kk,jj,i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real srcval(const SeamArgs &a, const int kk, const int jj, const int i) const {
    const int m = a.m, vv = a.vv;
    if (!a.cs_xform) {
      if (vv == 0) return bx1(a,kk,jj,i)*a.signvar;
      if (vv == 1) return bx2(a,kk,jj,i)*a.signvar;
      return bx3(a,kk,jj,i)*a.signvar;
    }
    Real bxi, bet;
    if (vv == 1) {
      bxi = bx2(a,kk,jj,i);
      bet = q_partner(a,kk,jj,i,true);
    } else {
      bet = bx3(a,kk,jj,i);
      bxi = q_partner(a,kk,jj,i,false);
    }
    const int cs_seam = a.cs_seam;
    if (cs_seam != 0) {
      const int js_ = a.js_, ks_ = a.ks_, je_ = a.je_, ke_ = a.ke_;
      const bool along_eta = (cs_seam == 2);
      const int vwant = along_eta ? 2 : 1;
      const bool corr_et = (vwant == 2) && (vv == 2);
      const bool corr_xi = (vwant == 1) && (vv == 1);
      const int an0 = along_eta ? kk : jj;
      const int nn0 = along_eta ? jj : kk;
      const int alo = along_eta ? ks_ : js_;
      const int ahi = along_eta ? ke_ : je_;
      const int nlo = along_eta ? js_ : ks_;
      const int nhi = along_eta ? je_ : ke_;
      const int nn1 = (nn0-1 >= nlo) ? (nn0-1) : (nn0+1);
      const Real sgn_n = (nn0-1 >= nlo) ? 1.0 : -1.0;
      const bool ok_n = (nn1 >= nlo && nn1 <= nhi);
      const bool ok_a = (an0-1 >= alo && an0+1 <= ahi);
      if ((corr_et || corr_xi) && ok_n && ok_a) {
        const int sidx_ = along_eta ? 0 : 1;
        const Real sden = cs_stbl(m,a.cs_skind,sidx_,kk,jj,5);
        if (fabs(sden) > 1.0e-8) {
          const Real sig = cs_stbl(m,a.cs_skind,sidx_,kk,jj,6);
          auto fv = [&](const int aa, const int nn) {
            const int kq = along_eta ? aa : nn;
            const int jq = along_eta ? nn : aa;
            return corr_et ? bx3(a,kq,jq,i) : bx2(a,kq,jq,i);
          };
          const Real mx = sgn_n*((fv(an0+1,nn0) - fv(an0-1,nn0))
                               - (fv(an0+1,nn1) - fv(an0-1,nn1)));
          const Real dsh = cs_shear_sgn*sig*mx/24.0;
          if (corr_et) { bet += dsh; } else { bxi += dsh; }
        }
      }
    }
    Real oxi, oet;
    cubed_sphere::SeamFieldXform t;
    const int sl = a.cs_slotidx, gk = a.cs_gkind;
    for (int c=0; c<3; ++c) {
      t.s1[c] = cs_gtbl(m,sl,gk,kk,jj,c);
      t.s2[c] = cs_gtbl(m,sl,gk,kk,jj,3+c);
      t.d1[c] = cs_gtbl(m,sl,gk,kk,jj,6+c);
      t.d2[c] = cs_gtbl(m,sl,gk,kk,jj,9+c);
    }
    t.sns = cs_gtbl(m,sl,gk,kk,jj,12);
    t.csd = cs_gtbl(m,sl,gk,kk,jj,13);
    t.snd = cs_gtbl(m,sl,gk,kk,jj,14);
    cubed_sphere::ApplyFieldXform(t, bxi, bet, oxi, oet);
    return (a.v == 1) ? oxi : oet;
  }
  //--------------------------------------------------------------------------------------
  // fcorder3-1103 <mesh>/cs_seam_fc_order = 3: the destination ghost face fed by the
  // source sample (kk,jj,i) of component vv, as the THIRD-ORDER average over the
  // receiver's own ghost face (see cs_seam_fc3.hpp for the why).  All in the SOURCE chart:
  //   n = seam-normal angle, a = along-seam angle, u = offset from the seam in normal
  //   cells (the two charts share it exactly: the ghost layers ARE the source layers);
  //   the receiver's along angle at the sample's along index equals the source's grid
  //   angle there (the grids coincide on the seam; reversal is in the index map), and a
  //   receiver point at (u, a~) is the source point a = atan(tan a~ tan|n|) -- the
  //   along-seam resample's own map (BuildCubedSphereSeamTables).
  // The receiver's face is one of three kinds, fixed by vv's staggering:
  //   T1 vv staggered across the seam: the face at normal face u, along cell [a~0,a~1];
  //   T2 vv staggered along it: the face at along face a~, normal cell [u0,u0+1];
  //   T3 vv = 0 (x1f): the radial face over the cell [u0,u0+1] x [a~0,a~1].
  // 2-point Gauss per extended direction of B.n_dst at points rebuilt from the source's
  // face integrals (PtN: the across-staggered component, PtA: the along-staggered one,
  // PtR: x1f), weighted by the receiver face's arc / solid-angle density.
  // STENCILS (layout independence under cs_seam_blockend = extend): the along window is
  // centred on the target but clamped to within ng cells/faces of the SAMPLE's own along
  // index (a sample is never more than ng outside its block's active range, so with the
  // extended pass's bounds this never binds at an interior block end and every layout
  // picks the same window), then to the source bounds (binding only at a panel end, or
  // at a block end in the clamp mode).  Across the seam: the 4 faces / 3 cells next to it.
  struct Fc3Geo {
    int m, i, aeta, low, dir, nseam, ns, ne, as0, alo, ahi, clo, chi, flo, fhi, nxa;
    Real dn, da, a0, amin, amax;
  };
  KOKKOS_FORCEINLINE_FUNCTION
  Real fc3_aface(const Fc3Geo &g, const int f) const {
    return 0.25*M_PI*LeftEdgeX(f - g.as0, g.nxa, g.amin, g.amax);
  }
  // normal cell m (0..2 from the seam) -> block index
  KOKKOS_FORCEINLINE_FUNCTION
  int fc3_ncell(const Fc3Geo &g, const int mc) const {
    return g.low ? (g.ns + mc) : (g.ne - mc);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real fc3_bN(const Fc3Geo &g, const int aa, const int nn) const {
    return g.aeta ? b.x2f(g.m,aa,nn,g.i) : b.x3f(g.m,nn,aa,g.i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real fc3_bA(const Fc3Geo &g, const int aa, const int nn) const {
    return g.aeta ? b.x3f(g.m,aa,nn,g.i) : b.x2f(g.m,nn,aa,g.i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  Real fc3_bR(const Fc3Geo &g, const int aa, const int nn) const {
    return g.aeta ? b.x1f(g.m,aa,nn,g.i) : b.x1f(g.m,nn,aa,g.i);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  int fc3_cwin(const Fc3Geo &g, const Real xabs) const {   // first of 3 along cells
    int c0 = static_cast<int>(floor(xabs - 0.5)) - 1;
    c0 = (c0 < g.clo) ? g.clo : ((c0 > g.chi) ? g.chi : c0);
    c0 = (c0 < g.alo) ? g.alo : ((c0 > g.ahi - 2) ? g.ahi - 2 : c0);
    return c0;
  }
  KOKKOS_FORCEINLINE_FUNCTION
  int fc3_fwin(const Fc3Geo &g, const Real xabs) const {   // first of 4 along faces
    int w = static_cast<int>(floor(xabs + 0.5)) - 2;
    w = (w < g.flo) ? g.flo : ((w > g.fhi) ? g.fhi : w);
    w = (w < g.alo) ? g.alo : ((w > g.ahi - 2) ? g.ahi - 2 : w);
    return w;
  }
  // The point values are LINEAR in the source faces and every weight is geometry only,
  // so each reconstruction hands (component, along index, normal index, weight) terms to
  // an accumulator: Fc3Direct sums them on the fly, Fc3Tab records them (the per-sample
  // coefficient table of PackAndSendFCFlat, evaluated once per run).  comp: 0 = x1f,
  // 1 = the across-staggered component, 2 = the along-staggered one.
  // across-staggered component (normal faces 0..3 from the seam) at (u, along face
  // coordinate xabs, along angle tangent Aa), times sc
  template <class Acc>
  KOKKOS_FORCEINLINE_FUNCTION
  void fc3_PtN(const Fc3Geo &g, const Real u, const Real xabs, const Real Aa,
               const Real sc, Acc &acc) const {
    const int c0 = fc3_cwin(g, xabs);
    Real wq[3], wl[4];
    cs_fc3::QuadDensW(xabs - 0.5 - static_cast<Real>(c0 + 1), wq);
    cs_fc3::Lag4W(u, wl);
    Real tf[4];
    for (int l=0; l<4; ++l) tf[l] = tan(fc3_aface(g, c0 + l));
    for (int mc=0; mc<4; ++mc) {
      const int col = g.nseam + g.dir*mc;
      const Real C = tan(0.25*M_PI - static_cast<Real>(mc)*g.dn);
      const Real f = sc*wl[mc]/(cs_fc3::ArcDens(C, Aa)*g.da);
      for (int l=0; l<3; ++l) acc(g, 1, c0 + l, col, f*wq[l]*cs_fc3::Arc(C, tf[l], tf[l+1]));
    }
  }
  // along-staggered component (along faces w..w+3, normal cells 0..2) at (u, xabs,
  // tangent of |n| = Nn), times sc
  template <class Acc>
  KOKKOS_FORCEINLINE_FUNCTION
  void fc3_PtA(const Fc3Geo &g, const Real u, const Real xabs, const Real Nn,
               const Real sc, Acc &acc) const {
    const int w = fc3_fwin(g, xabs);
    Real wq[3], wl[4];
    cs_fc3::QuadDensW(u - 1.5, wq);
    cs_fc3::Lag4W(xabs - static_cast<Real>(w), wl);
    Real tn[4];
    for (int mc=0; mc<4; ++mc) tn[mc] = tan(0.25*M_PI - static_cast<Real>(mc)*g.dn);
    for (int rr=0; rr<4; ++rr) {
      const Real R = tan(fc3_aface(g, w + rr));
      const Real f = sc*wl[rr]/(cs_fc3::ArcDens(R, Nn)*g.dn);
      for (int mc=0; mc<3; ++mc) {
        acc(g, 2, w + rr, fc3_ncell(g, mc), f*wq[mc]*cs_fc3::Arc(R, tn[mc+1], tn[mc]));
      }
    }
  }
  // x1f (cell integrals, 3 x 3) at (u, xabs), point tangents (Nn, Aa), times sc
  template <class Acc>
  KOKKOS_FORCEINLINE_FUNCTION
  void fc3_PtR(const Fc3Geo &g, const Real u, const Real xabs, const Real Nn,
               const Real Aa, const Real sc, Acc &acc) const {
    const int c0 = fc3_cwin(g, xabs);
    Real wa[3], wn[3];
    cs_fc3::QuadDensW(xabs - 0.5 - static_cast<Real>(c0 + 1), wa);
    cs_fc3::QuadDensW(u - 1.5, wn);
    Real tf[4], tn[4];
    for (int l=0; l<4; ++l) tf[l] = tan(fc3_aface(g, c0 + l));
    for (int mc=0; mc<4; ++mc) tn[mc] = tan(0.25*M_PI - static_cast<Real>(mc)*g.dn);
    const Real f = sc/(cs_fc3::SolidDens(Nn, Aa)*g.dn*g.da);
    for (int l=0; l<3; ++l) {
      for (int mc=0; mc<3; ++mc) {
        acc(g, 0, c0 + l, fc3_ncell(g, mc),
            f*wa[l]*wn[mc]*cs_fc3::Solid(tn[mc+1], tn[mc], tf[l], tf[l+1]));
      }
    }
  }
  // the receiver face fed by sample (kk,jj) as terms of acc; false: not an order-3 sample
  // (block not on a panel edge, or the shared seam face) -- the caller uses srcval
  template <class Acc>
  KOKKOS_INLINE_FUNCTION
  bool srcval3_acc(const SeamArgs &a, const int psrc, const int pdst, const int kk,
                   const int jj, const int i, Acc &acc) const {
    const int m = a.m, vv = a.vv;
    Fc3Geo g;
    g.m = m; g.i = i;
    g.aeta = (a.cs_seam == 2) ? 1 : 0;           // along = k (eta), normal = j (xi)
    const RegionSize sz = mbsz(m);
    const Real nmin = g.aeta ? sz.x2min : sz.x3min, nmax = g.aeta ? sz.x2max : sz.x3max;
    g.amin = g.aeta ? sz.x3min : sz.x2min; g.amax = g.aeta ? sz.x3max : sz.x2max;
    const int nxn = g.aeta ? cs_indcs.nx2 : cs_indcs.nx3;
    g.nxa = g.aeta ? cs_indcs.nx3 : cs_indcs.nx2;
    g.ns = g.aeta ? cs_indcs.js : cs_indcs.ks;
    g.ne = g.aeta ? cs_indcs.je : cs_indcs.ke;
    g.as0 = g.aeta ? cs_indcs.ks : cs_indcs.js;
    g.alo = g.aeta ? a.ks_ : a.js_;
    g.ahi = g.aeta ? a.ke_ : a.je_;
    const int nn = g.aeta ? jj : kk, an = g.aeta ? kk : jj;
    // which panel edge is the seam: the one the sample sits next to
    const bool elo = (nmin < -1.0 + 1.0e-10), ehi = (nmax > 1.0 - 1.0e-10);
    if (!elo && !ehi) return false;
    g.low = (elo && (!ehi || (nn - g.ns <= g.ne + 1 - nn))) ? 1 : 0;
    g.dir = g.low ? 1 : -1;
    g.nseam = g.low ? g.ns : (g.ne + 1);
    g.dn = 0.25*M_PI*(nmax - nmin)/static_cast<Real>(nxn);
    g.da = 0.25*M_PI*(g.amax - g.amin)/static_cast<Real>(g.nxa);
    g.a0 = 0.25*M_PI*g.amin;
    // the receiver face kind
    const int kind = (vv == 0) ? 3 : (((g.aeta != 0) == (vv == 1)) ? 1 : 2);
    Real u0, u1;
    if (kind == 1) {
      u0 = static_cast<Real>(g.dir*(nn - g.nseam));
      if (u0 < 0.5) return false;    // the shared seam face is never a ghost
      u1 = u0;
    } else {
      u0 = static_cast<Real>(g.low ? (nn - g.ns) : (g.ne - nn));
      u1 = u0 + 1.0;
    }
    const int ngh = cs_indcs.ng;
    // along stencil limits relative to the sample (see the note above)
    if (kind == 2) {
      g.clo = an - ngh; g.chi = an + ngh - 3; g.flo = an - ngh; g.fhi = an + ngh - 3;
    } else {
      g.clo = an - ngh; g.chi = an + ngh - 2; g.flo = an - ngh; g.fhi = an + ngh - 2;
    }
    const Real at0 = fc3_aface(g, an);
    const Real at1 = (kind == 2) ? at0 : fc3_aface(g, an + 1);
    const int nga = (kind == 2) ? 1 : 2;
    const int ngn = (kind == 1) ? 1 : 2;
    for (int pass=0; pass<2; ++pass) {   // 0: the weight sum, 1: the terms
    Real den = 0.0;
    for (int gn=0; gn<ngn; ++gn) {
      const Real u = (kind == 1) ? u0 :
          0.5*(u0 + u1) + 0.5*(u1 - u0)*((gn == 0) ? -cs_fc3::kGauss : cs_fc3::kGauss);
      const Real nabs = 0.25*M_PI - u*g.dn;
      const Real Nn = tan(nabs);
      const Real X = 1.0/Nn;                          // receiver's normal tangent
      for (int ga=0; ga<nga; ++ga) {
        const Real at = (kind == 2) ? at0 : 0.5*(at0 + at1)
            + 0.5*(at1 - at0)*((ga == 0) ? -cs_fc3::kGauss : cs_fc3::kGauss);
        const Real Y = tan(at);
        const Real wgt = (kind == 3) ? cs_fc3::SolidDens(X, Y)
                       : ((kind == 1) ? cs_fc3::ArcDens(X, Y) : cs_fc3::ArcDens(Y, X));
        if (pass == 0) {den += wgt; continue;}
        const Real sc = wgt/acc.den;
        const Real Aa = Y*Nn;
        const Real asrc = atan(Aa);
        const Real xabs = static_cast<Real>(g.as0) + (asrc - g.a0)/g.da;
        if (kind == 3) {
          fc3_PtR(g, u, xabs, Nn, Aa, sc, acc);
        } else {
          // the receiver component is linear in the source pair: its two coefficients
          const Real nsg = g.low ? -nabs : nabs;
          const Real xi = g.aeta ? nsg : asrc, eta = g.aeta ? asrc : nsg;
          cubed_sphere::SeamFieldXform t;
          cubed_sphere::SeamFieldXformAt(psrc, pdst, xi, eta, t);
          Real ox1, oe1, ox2, oe2;
          cubed_sphere::ApplyFieldXform(t, 1.0, 0.0, ox1, oe1);
          cubed_sphere::ApplyFieldXform(t, 0.0, 1.0, ox2, oe2);
          const Real cxi = (a.v == 1) ? ox1 : oe1, cet = (a.v == 1) ? ox2 : oe2;
          fc3_PtN(g, u, xabs, Aa, sc*(g.aeta ? cxi : cet), acc);
          fc3_PtA(g, u, xabs, Nn, sc*(g.aeta ? cet : cxi), acc);
        }
      }
    }
    if (pass == 0) acc.den = den;
    }
    return true;
  }
  // the direct evaluation (reads the faces at radial index g.i)
  struct Fc3Direct {
    const FCSeamPack *p; Real sum = 0.0, den = 1.0;
    KOKKOS_INLINE_FUNCTION
    void operator()(const Fc3Geo &g, const int comp, const int aa, const int nn,
                    const Real w) {
      const Real v = (comp == 0) ? p->fc3_bR(g, aa, nn)
                   : ((comp == 1) ? p->fc3_bN(g, aa, nn) : p->fc3_bA(g, aa, nn));
      sum += w*v;
    }
  };
  // the recorded terms, merged by face (code: array 0/1/2 = x1f/x2f/x3f, k, j)
  struct Fc3Tab {
    int n = 0; Real den = 1.0;
    int code[kFc3Max]; Real w[kFc3Max];
    KOKKOS_INLINE_FUNCTION
    void operator()(const Fc3Geo &g, const int comp, const int aa, const int nn,
                    const Real wt) {
      const int arr = (comp == 0) ? 0 : (((comp == 1) == (g.aeta != 0)) ? 1 : 2);
      const int k = g.aeta ? aa : nn, j = g.aeta ? nn : aa;
      const int c = (arr << 26) | (k << 13) | j;
      for (int q=0; q<n; ++q) { if (code[q] == c) { w[q] += wt; return; } }
      if (n < kFc3Max) { code[n] = c; w[n] = wt; ++n; } else { n = kFc3Max + 1; }
    }
  };
  KOKKOS_INLINE_FUNCTION
  Real srcval3(const SeamArgs &a, const int psrc, const int pdst, const int kk,
               const int jj, const int i) const {
    Fc3Direct acc; acc.p = this;
    if (!srcval3_acc(a, psrc, pdst, kk, jj, i, acc)) return srcval(a,kk,jj,i);
    return acc.sum;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Real seamval_pre(const SeamArgs &a, const int kk, const int jj, const int i) const {
    const int cs_seam = a.cs_seam;
    const int m = a.m;
    const int sidx_ = (cs_seam == 2) ? 0 : 1;
    const int bs = static_cast<int>(cs_stbl(m,a.cs_skind,sidx_,kk,jj,0));
    const Real wm = cs_stbl(m,a.cs_skind,sidx_,kk,jj,1);
    const Real w0 = cs_stbl(m,a.cs_skind,sidx_,kk,jj,2);
    const Real wp = cs_stbl(m,a.cs_skind,sidx_,kk,jj,3);
    // r1mspeed-1102: the srcval values evaluated once by the FCSeamPre pass
    const int nn = (cs_seam == 2) ? jj : kk;
    const int64_t base = a.pvoff + (static_cast<int64_t>(nn - a.pnlo)*a.pnt
                                    + (bs - a.ptlo))*a.pni + (i - a.pil);
    const Real sv0 = pre(base);
    const Real sv1 = pre(base + a.pni);
    const Real sv2 = pre(base + 2*a.pni);
    const Real qq = wm*sv0 + w0*sv1 + wp*sv2;
    const Real dd1 = sv1 - sv0, dd2 = sv2 - sv1;
    const Real fpos = cs_stbl(m,a.cs_skind,sidx_,kk,jj,4);
    if (!cs_nolim_ && dd1*dd2 > 0.0 && fpos >= 0.0 && fpos <= 2.0) {
      const Real lo0 = fmin(sv0, fmin(sv1, sv2));
      const Real hi0 = fmax(sv0, fmax(sv1, sv2));
      return fmin(hi0, fmax(lo0, qq));
    }
    return qq;
  }

  KOKKOS_INLINE_FUNCTION
  void SetArgs(const PackSegFC &s, const SubSegFC &u, SeamArgs &a) const {
    a.m = s.m; a.v = s.v; a.vv = s.vv; a.signvar = s.signvar; a.cs_xform = s.cs_xform;
    a.cs_coar = s.cs_coar; a.cs_slotidx = s.cs_slotidx; a.cs_gkind = s.cs_gkind;
    a.cs_skind = s.cs_skind; a.cs_seam = s.cs_seam;
    a.js_ = s.cs_coar ? cs_indcs.cjs : cs_indcs.js;
    a.ks_ = s.cs_coar ? cs_indcs.cks : cs_indcs.ks;
    a.je_ = s.cs_coar ? cs_indcs.cje : cs_indcs.je;
    a.ke_ = s.cs_coar ? cs_indcs.cke : cs_indcs.ke;
    if (ext_ && !s.cs_coar) {
      a.js_ -= ext_jk(s.m,0); a.je_ += ext_jk(s.m,1);
      a.ks_ -= ext_jk(s.m,2); a.ke_ += ext_jk(s.m,3);
    }
    a.wide2 = (a.je_ - a.js_ >= 3);
    a.wide3 = (a.ke_ - a.ks_ >= 3);
    a.pvoff = u.pvoff; a.ptlo = u.ptlo; a.pnt = u.pnt; a.pnlo = u.pnlo;
    a.pil = u.pil; a.pni = u.pni;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int t) const {
    const int c = c0 + t/FC_PACK_CHUNK;
    const SubSegFC u = sub(cseg(c));
    const int idx = cbeg(c) + t%FC_PACK_CHUNK;
    if (idx >= u.nk*u.nj*u.ni) return;
    const PackSegFC s = seg(u.q);
    const int ni = s.ni;
    SeamArgs a;
    SetArgs(s, u, a);
    const int uji = u.nj*u.ni;
    int k = (idx)/uji;
    int j = (idx - k*uji)/u.ni;
    int i = (idx - k*uji - j*u.ni) + u.il;
    k += u.kl;
    j += u.jl;
    int kk = s.ak*k + s.bk;
    int jj = s.aj*j + s.bj;
    const int ix = s.ndat*s.v + i-s.il + ni*(s.sj*(j-s.jl) + s.sk*(k-s.kl));
    // r1mspeed-1102: never srcval here (FCSeamPre did it): small, stack-free kernel
    const Real val = (u.pnt == 0) ? pre(u.pvoff + idx) : seamval_pre(a,kk,jj,i);
    if (ovs_ && u.pnt != 0) {   // fcorder3-1103 diagnostic: the order-2 resample nodes
      const int sidx_ = (a.cs_seam == 2) ? 0 : 1;
      const int bs = static_cast<int>(cs_stbl(a.m,a.cs_skind,sidx_,kk,jj,0));
      const int nn = (a.cs_seam == 2) ? jj : kk;
      const int64_t base = a.pvoff + (static_cast<int64_t>(nn - a.pnlo)*a.pnt
                                      + (bs - a.ptlo))*a.pni + (i - a.pil);
      const Real s0 = pre(base), s1 = pre(base + a.pni), s2 = pre(base + 2*a.pni);
      ovs_note(val, fmin(s0, fmin(s1, s2)), fmax(s0, fmax(s1, s2)),
               fmax(fabs(s0), fmax(fabs(s1), fabs(s2))));
    }
    if (s.local) {
      rbuf(s.dn).vars(s.dm, ix) = val;
    } else {
      sbuf(s.n).vars(s.m, ix) = val;
    }
  }

  // seamexact-1103: the partner interpolants of every ACTIVE face of block m, with the
  // extended pass's (panel-clamped) bounds -- the values that pass would compute here
  KOKKOS_INLINE_FUNCTION
  void operator()(const FCSeamQ &, const int t) const {
    const int n1 = cs_indcs.nx1 + 2*cs_indcs.ng;
    const int n2 = cs_indcs.nx2 + 2*cs_indcs.ng;
    const int n3 = cs_indcs.nx3 + 2*cs_indcs.ng;
    const int m = t/(n3*n2*n1);
    const int k = (t - m*n3*n2*n1)/(n2*n1);
    const int j = (t - m*n3*n2*n1 - k*n2*n1)/n1;
    const int i = t - m*n3*n2*n1 - k*n2*n1 - j*n1;
    const int js = cs_indcs.js, je = cs_indcs.je, ks = cs_indcs.ks, ke = cs_indcs.ke;
    SeamArgs a;
    a.m = m; a.cs_coar = 0;
    a.js_ = js - ext_jk(m,0); a.je_ = je + ext_jk(m,1);
    a.ks_ = ks - ext_jk(m,2); a.ke_ = ke + ext_jk(m,3);
    a.wide2 = (a.je_ - a.js_ >= 3);
    a.wide3 = (a.ke_ - a.ks_ >= 3);
    const bool kin = (k >= ks && k <= ke), jin = (j >= js && j <= je);
    if (kin && jin) {
      q(m,0,k,j,i) = x3f_at_etacentre(a,k,j,i);     // face j
      q(m,1,k,j,i) = x3f_at_etacentre(a,k,j+1,i);   // face j+1
      q(m,2,k,j,i) = x2f_at_xicentre(a,k,j,i);      // face k
      q(m,3,k,j,i) = x2f_at_xicentre(a,k+1,j,i);    // face k+1
    }
  }

  // fcorder3-1103: the precompute pass of <mesh>/cs_seam_fc_order = 3 -- every seam
  // segment is sample-indexed (BuildPackSegmentsFC) and its cross-panel samples go
  // through srcval3.  A kernel of its own, so the order-2 pass keeps its registers.
  KOKKOS_INLINE_FUNCTION
  void operator()(const FCSeamPre3 &, const int t) const {
    const int c = c0 + t/FC_PACK_CHUNK;
    const SubSegFC u = sub(cseg(c));
    const int idx = cbeg(c) + t%FC_PACK_CHUNK;
    const PackSegFC s = seg(u.q);
    if (u.pnt != 0 || idx >= u.nk*u.nj*u.ni) return;
    SeamArgs a;
    SetArgs(s, u, a);
    const int uji = u.nj*u.ni;
    int k = idx/uji;
    int j = (idx - k*uji)/u.ni;
    const int i = (idx - k*uji - j*u.ni) + u.il;
    k += u.kl;
    j += u.jl;
    if (s.cs_seam != 0 && !s.cs_coar && (s.cs_xform || s.v == 0) &&
        s.cs_psrc != s.cs_pdst) {
      const int gr = fc3tab_ ? t_gbase(u.q) : -1;
      const int row = gr + (k - s.kl)*s.nj + (j - s.jl);
      const int nc = (gr >= 0) ? t_cnt(row) : -1;
      if (nc >= 0) {
        const int m = s.m;
        Real sum = 0.0;
        for (int q=0; q<nc; ++q) {
          const int c = t_code(row,q);
          const int ar = c >> 26, kq = (c >> 13) & 8191, jq = c & 8191;
          const Real v = (ar == 0) ? b.x1f(m,kq,jq,i)
                       : ((ar == 1) ? b.x2f(m,kq,jq,i) : b.x3f(m,kq,jq,i));
          sum += t_w(row,q)*v;
        }
        pre(u.pvoff + idx) = sum;
        if (ovs_) {   // fcorder3-1103 diagnostic: per-component weight sum x data range
          Real W[3] = {0.0, 0.0, 0.0}, vl[3] = {1e300, 1e300, 1e300};
          Real vh[3] = {-1e300, -1e300, -1e300}, sc = 0.0;
          for (int q=0; q<nc; ++q) {
            const int c = t_code(row,q);
            const int ar = c >> 26, kq = (c >> 13) & 8191, jq = c & 8191;
            const Real v = (ar == 0) ? b.x1f(m,kq,jq,i)
                         : ((ar == 1) ? b.x2f(m,kq,jq,i) : b.x3f(m,kq,jq,i));
            W[ar] += t_w(row,q); vl[ar] = fmin(vl[ar], v); vh[ar] = fmax(vh[ar], v);
            sc = fmax(sc, fabs(v));
          }
          Real lo = 0.0, hi = 0.0;
          for (int ar=0; ar<3; ++ar) {
            if (vh[ar] < vl[ar]) continue;
            lo += (W[ar] >= 0.0) ? W[ar]*vl[ar] : W[ar]*vh[ar];
            hi += (W[ar] >= 0.0) ? W[ar]*vh[ar] : W[ar]*vl[ar];
          }
          ovs_note(sum, lo, hi, sc);
        }
      } else {
        pre(u.pvoff + idx) = srcval3(a, s.cs_psrc, s.cs_pdst, s.ak*k + s.bk,
                                     s.aj*j + s.bj, i);
      }
    } else {
      pre(u.pvoff + idx) = srcval(a, s.ak*k + s.bk, s.aj*j + s.bj, i);
    }
  }

  // fcorder3-1103: one coefficient-table row (seam segment, sample k,j); cnt -1: the
  // sample is evaluated directly (not an order-3 sample, or > kFc3Max faces)
  KOKKOS_INLINE_FUNCTION
  void operator()(const FCSeamTab &, const int t) const {
    const int q = t_rq(t);
    const PackSegFC s = seg(q);
    SubSegFC u;
    u.q = q; u.il = s.il; u.ni = s.ni; u.jl = s.jl; u.nj = s.nj; u.kl = s.kl; u.nk = s.nk;
    u.pvoff = -1; u.ptlo = 0; u.pnt = 0; u.pnlo = 0; u.pnn = 0; u.pil = 0; u.pni = 0;
    SeamArgs a;
    SetArgs(s, u, a);
    const int loc = t - t_gbase(q);
    const int k = loc/s.nj + s.kl, j = loc%s.nj + s.jl;
    Fc3Tab acc;
    const bool ok = srcval3_acc(a, s.cs_psrc, s.cs_pdst, s.ak*k + s.bk, s.aj*j + s.bj,
                                s.il, acc);
    if (!ok || acc.n > kFc3Max) { t_cnt(t) = -1; return; }
    t_cnt(t) = acc.n;
    for (int n=0; n<acc.n; ++n) { t_code(t,n) = acc.code[n]; t_w(t,n) = acc.w[n]; }
  }

  // r1mspeed-1102: srcval over a SubSegFC's (normal, tangent, i) box, once per pack
  KOKKOS_INLINE_FUNCTION
  void operator()(const FCSeamPre &, const int t) const {
    const int c = c0 + t/FC_PACK_CHUNK;
    const SubSegFC u = sub(cseg(c));
    const int idx = cbeg(c) + t%FC_PACK_CHUNK;
    const PackSegFC s = seg(u.q);
    SeamArgs a;
    SetArgs(s, u, a);
    if (u.pnt == 0) {   // sample-indexed: the pack's own (k,j,i) -> (kk,jj,i) map
      if (idx >= u.nk*u.nj*u.ni) return;
      const int uji = u.nj*u.ni;
      int k = idx/uji;
      int j = (idx - k*uji)/u.ni;
      const int i = (idx - k*uji - j*u.ni) + u.il;
      k += u.kl;
      j += u.jl;
      pre(u.pvoff + idx) = srcval(a, s.ak*k + s.bk, s.aj*j + s.bj, i);
      return;
    }
    const int tni = u.pnt*u.pni;
    if (idx >= u.pnn*tni) return;
    const int nq = idx/tni;
    const int tq = (idx - nq*tni)/u.pni;
    const int i = (idx - nq*tni - tq*u.pni) + u.pil;
    const int nn = nq + u.pnlo, ta = tq + u.ptlo;
    const int kk = (s.cs_seam == 2) ? ta : nn;
    const int jj = (s.cs_seam == 2) ? nn : ta;
    pre(u.pvoff + idx) = srcval(a,kk,jj,i);
  }
};

}  // namespace

//----------------------------------------------------------------------------------------
//! \!fn void MeshBoundaryValuesFC::PackAndSendFCFlat()
//! \brief halofc-1102: the cubed-sphere face-field pack (see BuildPackSegmentsFC).

void MeshBoundaryValuesFC::PackAndSendFCFlat(DvceFaceFld4D<Real> &b,
                                             DvceFaceFld4D<Real> &cb) {
  if (!fc_seg_ready) BuildPackSegmentsFC();
  // same host-side switches as the team kernel (see the notes there)
  static const Real cs_shear_env = (std::getenv("CS_SHEAR") != nullptr)
                                 ? std::atof(std::getenv("CS_SHEAR")) : 0.0;
  const Real cs_shear_sgn = cs_shear_env;
  const bool cs_nolim_ = pmy_pack->pmesh->cs_seam_fc_nolimit;
  auto sbuf = SendBufDv();
  auto rbuf = RecvBufDv();
  // r1mspeed-1102: walk set 0 (all ng layers) or 1 (depth 1), see bvals.hpp
  // seamexact-1103: with cs_seam_blockend = extend and > 2 blocks per panel side the
  // extended pass and its partner interpolants (CsSeamExtPartner) read the along-seam
  // ghosts to depth ng, which a depth-1 pack (<mhd>/rkg_halo_d1) leaves stale: full depth
  // fcorder3-1103: order 3's along stencils read the along-seam ghosts to depth ng in the
  // extended pass, so it too needs the full-depth pack there
  // (only where a panel side is split between blocks: with one block per panel side no
  // stencil reaches an along-seam ghost)
  const bool fc3_ = (pmy_pack->pmesh->cs_seam_fc_order == 3);
  const auto &pmi_ = pmy_pack->pmesh->mesh_indcs;
  const auto &mbi_ = pmy_pack->pmesh->mb_indcs;
  const bool fc3split_ = fc3_ && (pmi_.nx2 > mbi_.nx2 || pmi_.nx3 > mbi_.nx3);
  const int d = (pack_depth == 1 &&
                 !(cs_blockend_ext_ && (cs_ext_q_ || fc3split_))) ? 1 : 0;
  last_pack_depth_ = pack_depth;   // layoutbug-1103
  const int np = fc_np[d], ns = fc_ns[d], nr = fc_npre[d];
  if (np > 0) {
    FCPlainPack f(b, cb, sbuf, rbuf, fc_seg, fc_sub[d], fc_cs[d], fc_cb[d], 0);
    Kokkos::parallel_for("SendBuffFCplain", Kokkos::Experimental::require(
        Kokkos::RangePolicy<>(DevExeSpace(), 0, np*FC_PACK_CHUNK),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
  }
  // layoutbug-1103: the second pass of <mesh>/cs_seam_blockend = extend
  auto stbl_ = cs_ext_pass_ ? cs_seam_stnc_x : cs_seam_stnc;
  // fcorder3-1103 DIAGNOSTIC env FC3_OVS=<every>: ghost-face overshoot, printed by rank 0
  static const char *ovse_ = std::getenv("FC3_OVS");
  static DvceArray1D<Real> *ovsv_ = nullptr;
  if (ovse_ != nullptr && ovsv_ == nullptr) ovsv_ = new DvceArray1D<Real>("fc3_ovs", 4);
  if (nr > 0) {
    FCSeamPack f(b, cb, sbuf, rbuf, fc_seg, fc_sub[d], fc_cs[d], fc_cb[d], cs_seam_geom,
                 stbl_, fc_preval, np + ns, pmy_pack->pmesh->mb_indcs,
                 cs_shear_sgn, cs_nolim_);
    f.ext_ = cs_ext_pass_; f.ext_jk = cs_ext_jk;
    f.useq_ = cs_ext_pass_ && cs_ext_q_; f.q = cs_ext_q;   // seamexact-1103
    if (ovse_ != nullptr) {f.ovs_ = true; f.ovsv = *ovsv_;}
    if (fc3_) {   // fcorder3-1103
      f.fc3_ = true; f.mbsz = pmy_pack->pmb->mb_size.d_view;
      // the coefficient table of this pass (built once; env FC3_DIRECT = the direct sum)
      static const bool fc3dir_ = (std::getenv("FC3_DIRECT") != nullptr);
      if (!fc3dir_) {
        Fc3Tables &T = (*fc3_tables_)[this];
        const int ip = cs_ext_pass_ ? 1 : 0;
        if (T.nrow == 0) {
          auto hs = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), fc_seg);
          const int nseg = static_cast<int>(hs.extent(0));
          Kokkos::realloc(T.gbase, nseg);
          auto hg = Kokkos::create_mirror_view(T.gbase);
          std::vector<int> rq;
          for (int q=0; q<nseg; ++q) {
            const PackSegFC &g = hs(q);
            hg(q) = -1;
            if (g.mode == 2 && g.cs_seam != 0 && !g.cs_coar && (g.cs_xform || g.v == 0) &&
                g.cs_psrc != g.cs_pdst) {
              hg(q) = static_cast<int>(rq.size());
              for (int r=0; r<g.nk*g.nj; ++r) rq.push_back(q);
            }
          }
          T.nrow = (rq.empty()) ? 1 : static_cast<int>(rq.size());
          if (rq.empty()) rq.push_back(0);
          Kokkos::realloc(T.rq, T.nrow);
          auto hr = Kokkos::create_mirror_view(T.rq);
          for (int r=0; r<T.nrow; ++r) hr(r) = rq[r];
          Kokkos::deep_copy(T.gbase, hg);
          Kokkos::deep_copy(T.rq, hr);
        }
        if (!T.ready[ip]) {
          Kokkos::realloc(T.cnt[ip], T.nrow);
          Kokkos::realloc(T.code[ip], T.nrow, kFc3Max);
          Kokkos::realloc(T.w[ip], T.nrow, kFc3Max);
          f.t_gbase = T.gbase; f.t_rq = T.rq; f.t_cnt = T.cnt[ip];
          f.t_code = T.code[ip]; f.t_w = T.w[ip];
          Kokkos::parallel_for("SendBuffFCtab3",
              Kokkos::RangePolicy<FCSeamTab>(DevExeSpace(), 0, T.nrow), f);
          T.ready[ip] = true;
        }
        f.fc3tab_ = true; f.t_gbase = T.gbase; f.t_rq = T.rq; f.t_cnt = T.cnt[ip];
        f.t_code = T.code[ip]; f.t_w = T.w[ip];
      }
      Kokkos::parallel_for("SendBuffFCsrc3", Kokkos::Experimental::require(
          Kokkos::RangePolicy<FCSeamPre3>(DevExeSpace(), 0, nr*FC_PACK_CHUNK),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
    } else {
    Kokkos::parallel_for("SendBuffFCsrc", Kokkos::Experimental::require(
        Kokkos::RangePolicy<FCSeamPre>(DevExeSpace(), 0, nr*FC_PACK_CHUNK),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
    }
  }
  if (ns > 0) {
    FCSeamPack f(b, cb, sbuf, rbuf, fc_seg, fc_sub[d], fc_cs[d], fc_cb[d], cs_seam_geom,
                 stbl_, fc_preval, np, pmy_pack->pmesh->mb_indcs, cs_shear_sgn,
                 cs_nolim_);
    f.ext_ = cs_ext_pass_; f.ext_jk = cs_ext_jk;
    f.useq_ = cs_ext_pass_ && cs_ext_q_; f.q = cs_ext_q;   // seamexact-1103
    if (ovse_ != nullptr) {f.ovs_ = true; f.ovsv = *ovsv_;}
    Kokkos::parallel_for("SendBuffFCseam", Kokkos::Experimental::require(
        Kokkos::RangePolicy<>(DevExeSpace(), 0, ns*FC_PACK_CHUNK),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
  }
  if (ovse_ != nullptr) {
    static int ocall = 0;
    if ((++ocall % std::max(1, std::atoi(ovse_))) == 0) {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), *ovsv_);
      Real v[4] = {h(0), h(1), h(2), h(3)};
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, v, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, v + 2, 2, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
      if (global_variable::my_rank == 0) {
        std::printf("### FC3_OVS order %d pack %d (cumulative): max abs %.4e max rel %.4e"
                    " frac rel>1e-2 %.3e of %.0f samples\n", pmy_pack->pmesh->cs_seam_fc_order,
                    ocall, v[0], v[1], v[2]/fmax(v[3], 1.0), v[3]);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \!fn void MeshBoundaryValuesFC::PackAndSendFC()
//! \brief Pack face-centered Mesh variables into boundary buffers and send to neighbors.
//!
//! As for cell-centered data, this routine packs ALL the buffers on ALL the faces, edges,
//! and corners simultaneously for all three components of face-fields on ALL the
//! MeshBlocks.
//!
//! Input array must be DvceFaceFld4D dimensioned (nmb, nx3, nx2, nx1)
//! DvceFaceFld4D of coarsened (restricted) fields also required with SMR/AMR

TaskStatus MeshBoundaryValuesFC::PackAndSendFC(DvceFaceFld4D<Real> &b,
                                               DvceFaceFld4D<Real> &cb) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;

  // CUBED-SPHERE SEAM GEOMETRY TABLE.  Built once, lazily, on the first call: it needs
  // the neighbour table, which is not yet set when this object is constructed, and it
  // can never go stale because a cubed sphere refuses adaptive refinement and refuses
  // any static level boundary that would lie on a seam.  See BuildCubedSphereSeamTables.
  if (pmy_pack->pmesh->use_cubed_sphere && !cs_seam_tbl_ready) {
    BuildCubedSphereSeamTables();
  }

  // halofc-1102: on a cubed sphere the pack is the flat, load-balanced pair of kernels
  // (bitwise the same buffers); every other mesh keeps the team kernel below.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    PackAndSendFCFlat(b, cb);
  } else {int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  const bool use_pole = pmy_pack->pmesh->use_polar_boundary;
  const bool cs_nolim_ = pmy_pack->pmesh->cs_seam_fc_nolimit;   // see mesh.cpp
  // Sign/enable switch for the seam SHEAR CORRECTION, while its orientation is being
  // pinned down empirically.  The derivation fixes the magnitude; which way the
  // receiver's face tilts relative to the source's depends on the seam orientation and
  // reversal, and the seam-jump gate settles it in one run.  0 disables it entirely.
  // Read the environment ONCE on the host, then take a LOCAL copy: a function-scope
  // `static` is a host variable, and a device lambda that references one fails to compile
  // on HIP/CUDA ("reference to __host__ variable in __host__ __device__ function").  The
  // local const is captured BY VALUE and works on every backend.  Same trap as f493c209.
  static const Real cs_shear_env = (std::getenv("CS_SHEAR") != nullptr)
                                 ? std::atof(std::getenv("CS_SHEAR")) : 0.0;
  const Real cs_shear_sgn = cs_shear_env;

  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto cs_gtbl = cs_seam_geom;
  auto cs_stbl = cs_seam_stnc;
  auto cs_sltb = cs_seam_slot;
  auto sbuf = SendBufDv();
  auto rbuf = RecvBufDv();

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(three field components), one team
  // each, as in PackAndSendCC.  This used to be 3*nmb teams, each walking the neighbours
  // serially with a team_barrier after every one ("to prevent race condition in
  // overlapping assignments", a comment that belongs to the UNPACK, where ghost faces are
  // shared between buffers).  The pack only READS b/cb and every (m,n,v) writes its own
  // slice ndat*v.. of its own buffer (sbuf[n] row m, or rbuf[dn] row dm), so the triples
  // are independent.  With 3*nmb = 72 teams on a 132-SM GH200 the serial walk made this
  // the most expensive kernel of a cubed-sphere MHD cycle (0.93 ms per call; 09-30).
  int nmnv = 3*nmb*nnghbr;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nmnv, Kokkos::AUTO);
  BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/(3*nnghbr);
    const int n = (tmember.league_rank()/3)%nnghbr;
    const int v = tmember.league_rank()%3;

    {
      // only load buffers when neighbor exists.  CUBED-SPHERE CUBE VERTEX: skipped on
      // both sides; FillPanelCornersFC overwrites this corner block. See bvals.hpp.
      if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
        // if neighbor is at coarser level, use cindices to pack buffer
        // Note indices can be different for each component of face-centered field.
        int il, iu, jl, ju, kl, ku, ndat;

        // WHICH COMPONENT'S INDEX RANGE DO WE READ THE SOURCE WITH?
        //
        // v is the RECEIVER's component: the packer writes into buffer slot v but reads
        // this block's array vv, and the two differ on a panel seam whose tangential axes
        // are interchanged (swap_ax) -- every equatorial-to-polar seam. Ranges are per
        // component and genuinely different: for a -x2 neighbour isame[1] (x2f) is ng
        // FACES in j by nx3 CELLS in k, while isame[2] (x3f) is ng CELLS in j by nx3+1
        // FACES in k. Reading array vv over isame[v] therefore took the wrong index set
        // AND the wrong NUMBER of elements, so the receiver unpacked a differently sized
        // region -- an O(1) error that did not converge under refinement. Read the source
        // with its OWN component's range; the sj/sk transposition below still lays it out
        // the way the receiver expects, and the counts match because nx2 == nx3 here.
        int vsrc = v;
        if (use_cs &&
            nghbr.d_view(m,n).panel != mbpanel.d_view(m) && v > 0) {
          PanelBoundaries pb0;
          pb0 = GetPanelBoundary(mbpanel.d_view(m),
                                                  nghbr.d_view(m,n).panel);
          if (pb0.swap_ax == 1) vsrc = (v == 1) ? 2 : 1;
        }

        if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
          il = sbuf[n].icoar[vsrc].bis;
          iu = sbuf[n].icoar[vsrc].bie;
          jl = sbuf[n].icoar[vsrc].bjs;
          ju = sbuf[n].icoar[vsrc].bje;
          kl = sbuf[n].icoar[vsrc].bks;
          ku = sbuf[n].icoar[vsrc].bke;
          ndat = sbuf[n].icoar_ndat;
        // if neighbor is at same level, use sindices to pack buffer
        } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
          il = sbuf[n].isame[vsrc].bis;
          iu = sbuf[n].isame[vsrc].bie;
          jl = sbuf[n].isame[vsrc].bjs;
          ju = sbuf[n].isame[vsrc].bje;
          kl = sbuf[n].isame[vsrc].bks;
          ku = sbuf[n].isame[vsrc].bke;
          ndat = sbuf[n].isame_ndat;
        // if neighbor is at finer level, use findices to pack buffer
        } else {
          il = sbuf[n].ifine[vsrc].bis;
          iu = sbuf[n].ifine[vsrc].bie;
          jl = sbuf[n].ifine[vsrc].bjs;
          ju = sbuf[n].ifine[vsrc].bje;
          kl = sbuf[n].ifine[vsrc].bks;
          ku = sbuf[n].ifine[vsrc].bke;
          ndat = sbuf[n].ifine_ndat;
        }
        const int ni = iu - il + 1;
        const int nj = ju - jl + 1;
        const int nk = ku - kl + 1;
        const int nkji = nk*nj*ni;
        const int nji  = nj*ni;

        // indices of recv'ing MB and buffer: assumes MB IDs are stored sequentially
        int dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
        int dn = nghbr.d_view(m,n).dest;

          const bool do_cs = use_cs &&
                             (nghbr.d_view(m,n).panel != mbpanel.d_view(m));
          const bool do_pole = use_pole &&
                               (nghbr.d_view(m,n).polar > 0);

          if (do_cs || do_pole) {
            int aj = 1, bj = 0;
            int ak = 1, bk = 0;
            int signvar = 1;
            int sj = 1, sk = nj;
            int vv = v;
            bool cs_xform = false;
            // A cross-panel neighbour at a COARSER level is served from the restricted
            // field cb, not b.  Before this existed the branch packed NOTHING, so the
            // buffer kept its zero-initialised contents and the coarse block read a
            // ghost field of exactly ZERO across the seam.  See the note in bvals_cc.cpp.
            const bool cs_coar = (nghbr.d_view(m,n).lev < mblev.d_view(m));
            int cs_srcpanel = 0, cs_dstpanel = 0;
            // where this (destination panel, staggering) pair sits in the seam tables
            int cs_slotidx = 0, cs_gkind = 0, cs_skind = 0;
            // 0 = no along-seam resample; 2 = x2-face seam (resample in k);
            // 3 = x3-face seam (resample in j). See the note on seamval below.
            int cs_seam = 0;

            if (do_cs) {
              const auto ngh = nghbr.d_view(m,n);
              const int my_panel = mbpanel.d_view(m);
              PanelBoundaries pb;
              pb = GetPanelBoundary(my_panel, ngh.panel);

              int map_vy = 1;
              int map_vz = 2;

              int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
              int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
              if (pb.swap_ax == 1) {
                map_vy = 2;
                map_vz = 1;
                sj = nk;
                sk = 1;
              }
              if (rev_a_preswap) {
                aj = -1;
                bj = jl + ju;
              }
              if (rev_b_preswap) {
                ak = -1;
                bk = kl + ku;
              }
              if (v == 1) vv = map_vy;
              if (v == 2) vv = map_vz;

              // x1 is RADIAL and rhat is common to both charts, so b.x1f crosses a seam
              // as a SCALAR: index map only. The two angular components do not. They are
              // projections on the panel's own FACE NORMALS, and the two charts carry
              // different normals at the same physical point -- differing by a shear that
              // is O(1) away from the seam midline, exactly as for the cell-centred
              // momentum (see cubed_sphere::TransformMomentum). A signed axis permutation
              // is exact for the INDICES and exact for the components only on the seam
              // midline, so it left an O(1) error in the halo field. Measured: it made a
              // uniform Cartesian field, an exact static state, blow up in one step.
              cs_xform = (v == 1) || (v == 2);
              cs_srcpanel = my_panel;
              cs_dstpanel = ngh.panel;
              cs_slotidx = cs_sltb.d_view(m, ngh.panel);
              cs_gkind = (cs_coar ? 2 : 0) + ((vv > 0) ? (vv - 1) : 0);
              cs_skind = (cs_coar ? 3 : 0) + vv;

              // Which buffer is this? A same-level FACE buffer is ng deep in its own
              // normal direction and spans the full active range in the other tangential
              // one; that is the only case the along-seam resample applies to. The edge
              // and corner buffers (ng x ng) stay a plain copy, as in bvals_cc.cpp.
              // WHICH BUFFER IS THIS, and therefore which way does the along-seam
              // resample run?  Read it off the SLOT INDEX, not off the buffer extents.
              // Slots are laid out (nghbr_index.hpp): 0-7 x1 faces, 8-15 x2 faces,
              // 16-23 x1x2 edges, 24-31 x3 faces, 32-39 x3x1 edges, 40-47 x2x3 edges,
              // 48-55 corners.  A buffer whose only tangential ghost direction is x2 has
              // its seam normal along x2 and resamples in x3 (cs_seam = 2), and vice
              // versa; one that is ghost in BOTH tangential directions -- the x2x3 edges
              // and the corners -- has no single along-seam axis and stays a plain copy.
              //
              // This USED TO BE INFERRED from the extents, as `nj == ng && nk > ng`.
              // That test is only valid at a SAME-LEVEL boundary.  At a coarse/fine
              // boundary a face-centred buffer carries one EXTRA layer in the direction
              // its component is staggered in, so the ghost-direction extent is ng+1 = 3
              // rather than ng, the equality fails, and cs_seam silently fell to 0 --
              // skipping the resample for exactly ONE of the two angular components on
              // every cross-level seam edge buffer.  Measured: the resampled source
              // sample then sat 1.82 CELLS away from the destination cell it fed (against
              // 0.0000 whenever the resample does run), which is O(h) and was the whole
              // of the first-order cross-level seam halo error.  The radial component is
              // not staggered in either tangential direction, kept its extent of ng, and
              // so kept the resample -- which is exactly why it was ten times less wrong.
              if (n >= 8 && n < 24) {
                cs_seam = 2;
              } else if (n >= 24 && n < 40) {
                cs_seam = 3;
              }
              // DOUBLY-GHOST BUFFERS (x2x3 edges 40-47, corners 48-55) are ghost in
              // BOTH tangential directions, so the slot alone does not name an along-seam
              // axis -- but the NEIGHBOUR TABLE does.  Only one of the two flanking faces
              // is a panel seam (if both were, this is a cube vertex whose exchange is
              // skipped as non-reciprocal), and that one fixes the seam normal exactly as
              // for a face buffer.
              if (n >= 40 && n < 56) {
                int sjc, skc;
                if (n < 48) {
                  const int cq = (n - 40)/2;
                  sjc = cq & 1;  skc = (cq >> 1) & 1;
                } else {
                  const int cq = n - 48;
                  sjc = (cq >> 1) & 1;  skc = (cq >> 2) & 1;
                }
                const int njid = sjc ? 12 : 8;
                const int nkid = skc ? 28 : 24;
                const bool sj_ = (nghbr.d_view(m,njid).gid >= 0 &&
                                  nghbr.d_view(m,njid).panel != my_panel);
                const bool sk_ = (nghbr.d_view(m,nkid).gid >= 0 &&
                                  nghbr.d_view(m,nkid).panel != my_panel);
                if (sj_ && !sk_) {
                  cs_seam = 2;
                } else if (sk_ && !sj_) {
                  cs_seam = 3;
                }
              }

              // THE RESAMPLE NEEDS 3 CELLS ALONG THE SEAM.  Its quadratic stencil is
              // clamped to [blo, bhi] = [start, start + extent - 3], which INVERTS when
              // the source's along-seam ACTIVE extent is under 3 -- the clamp then pins
              // the stencil one cell outside the active zone and it reads unfilled ghost.
              // That happens on the COARSE array as soon as cnx = nx/2 < 3, i.e. for a
              // MeshBlock under 6 cells wide, and MeshBlocks are allowed down to 4.
              // Measured on a refined 4x4 block: halo Linf 7.7e-01 against 2.4e-02 before
              // the resample was extended to these buffers -- the degeneracy predates
              // this session, but a01ace75 and e4be92c0 widened its blast radius by
              // resampling many more buffers, so it has to be guarded here.  Falling back
              // to the plain index copy costs the O(dx) seam offset the resample exists
              // to remove, which is bad but bounded; reading unfilled memory is not.
              const int seam_extent = cs_coar
                  ? ((cs_seam == 2) ? cs_indcs.cnx3 : cs_indcs.cnx2)
                  : ((cs_seam == 2) ? cs_indcs.nx3 : cs_indcs.nx2);
              if (cs_seam != 0 && seam_extent < 3) { cs_seam = 0; }
            } else if (do_pole) {
              aj = -1;
              bj = jl + ju;
              if (v == 1) signvar = -1;
              if (v == 2) signvar = -1;
            }

            // Value of destination component v, read from the source at (kk,jj,i).
            //
            // Producing either angular component needs BOTH source angular components at
            // the SAME point, and they are staggered differently: b.x2f sits on a xi face
            // and b.x3f on an eta face. The one the index map selects (vv) is taken where
            // it lives and the other is interpolated onto that location, which is second
            // order. The radial component is untouched.
            //
            // THAT INTERPOLATION MUST NOT REACH INTO THE SOURCE BLOCK'S OWN GHOSTS. The
            // secondary component is averaged across the ALONG-SEAM index, and the packed
            // range spans the whole seam, so at its two END faces one of the two cells is
            // the source block's ghost -- which on a panel seam is itself halo data, in
            // general stale and in a first exchange uninitialised. That produced an O(1)
            // error at exactly the two end faces of every seam and nowhere else
            // (measured: interior faces 3e-3..1.7e-2 at nx=32, end faces 0.38), and
            // since the halo check maxes over the seam that made all 24 seams look
            // O(1) and non-convergent.
            // At an end face the two nearest ACTIVE cells are extrapolated instead, which
            // keeps the same second-order accuracy without leaving the active zone.
            const int js_ = cs_coar ? cs_indcs.cjs : cs_indcs.js;
            const int ks_ = cs_coar ? cs_indcs.cks : cs_indcs.ks;
            const int je_ = cs_coar ? cs_indcs.cje : cs_indcs.je;
            const int ke_ = cs_coar ? cs_indcs.cke : cs_indcs.ke;
            // Every read below goes through these, so the whole seam transform works on
            // the restricted field cb unchanged when the neighbour is coarser.
            auto bx1 = [&](const int kk, const int jj, const int i) {
              return cs_coar ? cb.x1f(m,kk,jj,i) : b.x1f(m,kk,jj,i);
            };
            auto bx2 = [&](const int kk, const int jj, const int i) {
              return cs_coar ? cb.x2f(m,kk,jj,i) : b.x2f(m,kk,jj,i);
            };
            auto bx3 = [&](const int kk, const int jj, const int i) {
              return cs_coar ? cb.x3f(m,kk,jj,i) : b.x3f(m,kk,jj,i);
            };
            //
            // CO-LOCATING THE PARTNER COMPONENT MUST BE 4th ORDER ON *BOTH* AXES.
            //
            // The partner is moved onto the primary's location by two successive
            // interpolations on ORTHOGONAL axes (see `srcval`). A 2-point midpoint
            // average is second order in VALUE, and the resistive curl differences the
            // halo over h, so an O(h^2) ghost value is a FIRST ORDER current -- which is
            // exactly the measured seam order of 1.00. Two comparable O(h^2) terms on
            // orthogonal axes means raising ONE of them halves the amplitude and cannot
            // move the order; that was measured, and it is why the earlier one-axis
            // experiment failed. Both axes are raised here.
            //
            // PROVED BEFORE IT WAS WRITTEN, by an oracle A/B that replaced only the
            // partner with the analytic field at the primary's own (xi,eta): every
            // category that sat at 1.75-1.90 -- including the two SWAP ones that were
            // immune to every earlier switch -- went to 2.92-3.26, with the amplitude
            // down 9-78x, while x1f (which has no partner) stayed bitwise identical.
            // So the co-location is the WHOLE binding term in the seam halo.
            //
            // A CENTRED wide stencil is NOT ENOUGH, and the reason is structural: a
            // seam buffer packs the cells at the EDGE of the active range, so in the
            // seam-NORMAL direction a centred 4-point formula never fits and silently
            // degrades to the 2-point average over the WHOLE packed range -- which
            // measured 1.80/1.90, i.e. exactly the baseline, for a strictly more
            // expensive stencil. What is needed is a fourth-order interpolation whose
            // WINDOW SLIDES INWARD instead of degrading.
            //
            // `lag4` is the cubic through four consecutive samples starting at index w,
            // evaluated at offset t = target - w. The window is CLAMPED into the active
            // range, so it never reads the source block's own ghosts -- which on a panel
            // seam are themselves halo data, in general stale and uninitialised in a
            // first exchange -- and it stays fourth order right up to the edge, becoming
            // one-sided rather than low order. At the two ends of a seam the target lies
            // outside the active samples and the same window extrapolates, still cubic.
            auto lag4 = [](const Real t, const Real f0, const Real f1,
                                         const Real f2, const Real f3) {
              return -(t-1.0)*(t-2.0)*(t-3.0)/6.0*f0
                     + t*(t-2.0)*(t-3.0)/2.0*f1
                     - t*(t-1.0)*(t-3.0)/2.0*f2
                     + t*(t-1.0)*(t-2.0)/6.0*f3;
            };
            // Fewer than four samples on an axis leaves no cubic window; fall back to
            // the original second-order forms there. Same degeneracy the along-seam
            // resample guards against, and for the same reason.
            const bool wide2 = (je_ - js_ >= 3);
            const bool wide3 = (ke_ - ks_ >= 3);
            // b.x3f at (eta face kf, xi FACE jf), interpolated across xi. b.x3f lives at
            // xi CELL CENTRES, and xi face jf sits at cell-index jf-0.5.
            auto x3f_at_xiface = [&](const int kf, const int jf, const int i) {
              if (!wide2) {
                if (jf-1 < js_) { return 1.5*bx3(kf,js_,i) - 0.5*bx3(kf,js_+1,i); }
                if (jf   > je_) { return 1.5*bx3(kf,je_,i) - 0.5*bx3(kf,je_-1,i); }
                return 0.5*(bx3(kf,jf-1,i) + bx3(kf,jf,i));
              }
              int w = jf-2;
              w = (w < js_) ? js_ : ((w > je_-3) ? je_-3 : w);
              return lag4(static_cast<Real>(jf) - 0.5 - static_cast<Real>(w),
                          bx3(kf,w,i), bx3(kf,w+1,i), bx3(kf,w+2,i), bx3(kf,w+3,i));
            };
            // b.x2f at (eta FACE kf, xi face jf), interpolated across eta
            auto x2f_at_etaface = [&](const int kf, const int jf, const int i) {
              if (!wide3) {
                if (kf-1 < ks_) { return 1.5*bx2(ks_,jf,i) - 0.5*bx2(ks_+1,jf,i); }
                if (kf   > ke_) { return 1.5*bx2(ke_,jf,i) - 0.5*bx2(ke_-1,jf,i); }
                return 0.5*(bx2(kf-1,jf,i) + bx2(kf,jf,i));
              }
              int w = kf-2;
              w = (w < ks_) ? ks_ : ((w > ke_-3) ? ke_-3 : w);
              return lag4(static_cast<Real>(kf) - 0.5 - static_cast<Real>(w),
                          bx2(w,jf,i), bx2(w+1,jf,i), bx2(w+2,jf,i), bx2(w+3,jf,i));
            };
            // The OUTER interpolation, on the axis ORTHOGONAL to the one above: it takes
            // values living on the eta FACES (index range ks_..ke_+1) to the eta CENTRE
            // of cell kk, which sits at face-index kk+0.5. Same clamped cubic window.
            auto x3f_at_etacentre = [&](const int kk, const int jj, const int i) {
              if (!wide3) {
                return 0.5*(x3f_at_xiface(kk,jj,i) + x3f_at_xiface(kk+1,jj,i));
              }
              int w = kk-1;
              w = (w < ks_) ? ks_ : ((w > ke_-2) ? ke_-2 : w);
              return lag4(static_cast<Real>(kk) + 0.5 - static_cast<Real>(w),
                          x3f_at_xiface(w,jj,i),   x3f_at_xiface(w+1,jj,i),
                          x3f_at_xiface(w+2,jj,i), x3f_at_xiface(w+3,jj,i));
            };
            auto x2f_at_xicentre = [&](const int kk, const int jj, const int i) {
              if (!wide2) {
                return 0.5*(x2f_at_etaface(kk,jj,i) + x2f_at_etaface(kk,jj+1,i));
              }
              int w = jj-1;
              w = (w < js_) ? js_ : ((w > je_-2) ? je_-2 : w);
              return lag4(static_cast<Real>(jj) + 0.5 - static_cast<Real>(w),
                          x2f_at_etaface(kk,w,i),   x2f_at_etaface(kk,w+1,i),
                          x2f_at_etaface(kk,w+2,i), x2f_at_etaface(kk,w+3,i));
            };
            //
            // The (xi,eta) of a source sample, with the STAGGERING of component vv:
            // b.x2f sits on a xi face and b.x3f on an eta face, b.x1f on neither.  That
            // is what `cs_gkind` and `cs_skind` name, and everything the seam branch
            // used to compute from (xi,eta) per CELL -- the transform geometry, the
            // resample stencil and the shear term's sigma -- is read out of the tables
            // below at the sample's own (kk,jj).  The rows are indexed by the ABSOLUTE
            // (kk,jj), ghosts included: on a multilevel mesh a same-level send buffer is
            // widened by ng in its transverse directions and the packer reads source
            // cells outside the active range.  The fill evaluates exactly what this
            // lambda used to hand to LeftEdgeX/CellCenterX.
            auto srcval = [&](const int kk, const int jj, const int i) {
              if (!cs_xform) {
                if (vv == 0) return bx1(kk,jj,i)*signvar;
                if (vv == 1) return bx2(kk,jj,i)*signvar;
                return bx3(kk,jj,i)*signvar;
              }
              Real bxi, bet;
              if (vv == 1) {
                // primary on a XI face: (xi face jj, eta centre kk)
                bxi = bx2(kk,jj,i);
                bet = x3f_at_etacentre(kk,jj,i);
              } else {
                // primary on an ETA face: (xi centre jj, eta face kk)
                bet = bx3(kk,jj,i);
                bxi = x2f_at_xicentre(kk,jj,i);
              }
              // ---- SHEAR CORRECTION: the two charts' faces are not the same surface --
              //
              // The stored quantity is a FACE AVERAGE, and the receiver's ghost face is
              // not the source face that feeds it.  Take the non-swap +x / +y seam.  The
              // seam-NORMAL coordinate maps exactly -- xi_A + xi_B = pi/2 -- so the ghost
              // LAYERS coincide with the source's active layers and nothing is wrong in
              // that direction.  The faces are bounded in the other one, and there the
              // receiver's face is a surface of constant eta_A = atan(z/x) while the
              // source's is constant eta_B = atan(z/y).  Off the seam those are DIFFERENT
              // SURFACES, sheared relative to each other.
              //
              // Write the receiver's face in the source's chart: constant eta_A means
              // Y_B = X_B*tan(eta_A), so along the face eta_B drifts with xi_B at a rate
              //     sigma = d eta_B / d xi_B = sin(2 eta_B) / sin(2 xi_B),
              // which is ZERO at the seam midpoint (eta_B = 0) and grows monotonically to
              // 1 at the cube vertex.  Averaging f over a line tilted by sigma instead of
              // a flat one gives, with the linear term killed by centring,
              //  <f>_recv - <f>_src = sigma*(dxi^2/12)*d2f/dxi deta + O(sigma^2 h^2)
              //
              // MEASURED before this was written (cs_test's seam-jump gate, nlim=0, so
              // the halo machinery alone): the component whose face is the SHARED seam
              // surface reaches 3.13/3.13, while the other one sits at 1.77/1.89 and its
              // error grows 7.5x monotonically from the seam midpoint to the cube vertex
              // -- the sigma profile.  A halo that is second order in VALUE is a FIRST
              // order current, because the resistive curl differences it over h.
              //
              // Only the component whose face NORMAL is the ALONG-SEAM direction needs
              // this: the other tangential face IS the seam surface, shared by both
              // charts, and b.x1f's radial face is r = const in every chart.
              if (cs_seam != 0) {
                // WHICH AXIS IS ALONG THE SEAM, IN THE SOURCE'S FRAME.  cs_seam is read
                // off the DESTINATION's slot, and across a SWAP seam the neighbour
                // reaches this face through a face of its OTHER tangential axis, so the
                // source's normal and along-seam roles are exchanged.  A swap is exactly
                // vv != v (the component index the source is read with, set above).
                // Getting this wrong silently corrects the WRONG component -- which is
                // why the first version fired on only 4 of the 12 seams and moved the
                // ring EMF by 10%.
                // THE ANGLE AND INDEX ROLES COME FROM cs_seam, NOT FROM THE SWAP.  The
                // buffer's index map already absorbs the axis exchange, and the
                // along-seam resample -- validated to place its samples exactly -- reads
                // its normal and along-seam angles straight off cs_seam in the SOURCE's
                // (jj,kk) space.  Swapping the roles here instead made nrm the angle that
                // is NOT pinned at +-pi/4, so sin(2*nrm) went small and sigma blew up:
                // both signs then made the swap seams WORSE by the same amount, which is
                // an over-large term rather than a sign error.
                const bool along_eta = (cs_seam == 2);
                // The corrected component is the one whose face NORMAL is the along-seam
                // direction; the other tangential face IS the shared seam surface.
                // ONLY THE COMPONENT FLIPS.  The face that must be corrected is the one
                // whose normal is the along-seam direction -- but the transform maps
                // source components to destination components THROUGH the swap, so on a
                // swap seam it is the source's OTHER tangential face that lands on it.
                // MEASURED: the component rule does NOT change across a swap seam, and
                // neither does the sign.  Flipping it (the obvious guess, since the
                // transform maps components through the swap) measured WORSE on the swap
                // seams with either sign.  What had made the swap seams look untouched
                // was the stencil bug above, not the component choice.
                const int vwant = along_eta ? 2 : 1;
                const bool corr_et = (vwant == 2) && (vv == 2);
                const bool corr_xi = (vwant == 1) && (vv == 1);
                // THE SOURCE CELLS SIT AT THE EDGE OF THE ACTIVE RANGE -- that is what a
                // seam buffer packs -- so a CENTRED difference in the seam-NORMAL
                // direction reaches outside the data and never fires, silently and at any
                // amplitude.  One-sided stepping INWARD there, centred along the seam
                // where the full range is available.  One-sided is first order, but the
                // whole term is already O(h^2), so the residual is O(h^3).
                //
                // The NORMAL index is jj when the along-seam axis is eta, and kk when it
                // is xi -- it follows `along_eta`, not the slot.
                // THE STENCIL FOLLOWS THE AXIS ROLES, THE COMPONENT IS CHOSEN
                // SEPARATELY.  Keying the difference on which COMPONENT is corrected
                // instead silently indexes a j-neighbour as a k-neighbour the moment the
                // two disagree, which is exactly what the swap case does.
                //   along index  a = kk when along_eta, else jj
                //   normal index n = jj when along_eta, else kk
                const int an0 = along_eta ? kk : jj;          // along index, centre
                const int nn0 = along_eta ? jj : kk;          // normal index, at the edge
                const int alo = along_eta ? ks_ : js_;
                const int ahi = along_eta ? ke_ : je_;
                const int nlo = along_eta ? js_ : ks_;
                const int nhi = along_eta ? je_ : ke_;
                const int nn1 = (nn0-1 >= nlo) ? (nn0-1) : (nn0+1);
                const Real sgn_n = (nn0-1 >= nlo) ? 1.0 : -1.0;
                const bool ok_n = (nn1 >= nlo && nn1 <= nhi);
                const bool ok_a = (an0-1 >= alo && an0+1 <= ahi);
                if ((corr_et || corr_xi) && ok_n && ok_a) {
                  const int sidx_ = along_eta ? 0 : 1;
                  const Real sden = cs_stbl(m,cs_skind,sidx_,kk,jj,5);
                  if (fabs(sden) > 1.0e-8) {
                    const Real sig = cs_stbl(m,cs_skind,sidx_,kk,jj,6);
                    // f at (along index, normal index), whichever way round they are
                    auto fv = [&](const int aa, const int nn) {
                      const int kq = along_eta ? aa : nn;
                      const int jq = along_eta ? nn : aa;
                      return corr_et ? bx3(kq,jq,i) : bx2(kq,jq,i);
                    };
                    // d2f/dn da: centred ALONG the seam, one-sided ACROSS it.  sgn_n
                    // undoes the sign of the inward step so the derivative keeps its
                    // orientation whichever end of the block the source cells sit at.
                    const Real mx = sgn_n*((fv(an0+1,nn0) - fv(an0-1,nn0))
                                         - (fv(an0+1,nn1) - fv(an0-1,nn1)));
                    const Real dsh = cs_shear_sgn*sig*mx/24.0;
                    if (corr_et) { bet += dsh; } else { bxi += dsh; }
                  }
                }
              }
              // TransformFieldToDstNormals, with its geometry half read from the table
              // instead of recomputed: the two are the same call, so this is bitwise
              // what the inline version produced.
              Real oxi, oet;
              cubed_sphere::SeamFieldXform t;
              for (int c=0; c<3; ++c) {
                t.s1[c] = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,c);
                t.s2[c] = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,3+c);
                t.d1[c] = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,6+c);
                t.d2[c] = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,9+c);
              }
              t.sns = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,12);
              t.csd = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,13);
              t.snd = cs_gtbl(m,cs_slotidx,cs_gkind,kk,jj,14);
              cubed_sphere::ApplyFieldXform(t, bxi, bet, oxi, oet);
              return (v == 1) ? oxi : oet;
            };

            // ALONG-SEAM RESAMPLE, the face-centred twin of the one in bvals_cc.cpp.
            // Across a seam the two charts share the seam-NORMAL coordinate exactly but
            // not the seam-PARALLEL one: a source sample at seam-normal angle n and
            // seam-parallel angle a sits at seam-parallel angle atan(tan(a)/tan|n|) in
            // the destination chart. Inverting that, the value a ghost needs is the
            // source field at atan(tan(a)*tan|n|), always closer to the seam midline than
            // the sample's own angle, so the stencil never leaves the source block. The
            // offset is a fixed fraction of a cell -- resolution-INDEPENDENT in cell
            // units -- so the plain index copy leaves the halo field O(dx) wrong on every
            // component, b.x1f included, even though b.x1f crosses the seam as a scalar.
            // Quadratic (3-point) Lagrange, on the source's own index axis; the reversal
            // needs no special case because the map is odd in a and the index reversal
            // already carries the sign.
            //
            // The staggering matters here and it is not the same for the two angular
            // components: on a x3-face seam the along-seam axis is xi, which b.x2f
            // samples on FACES and b.x1f/b.x3f on cell CENTRES. The offset in INDEX units
            // is the same either way, so only the angle that enters the map changes,
            // which `angles` already knows.
            auto seamval = [&](const int kk, const int jj, const int i) {
              if (cs_seam == 0) return srcval(kk,jj,i);
              // The stencil is geometry: its index and its three weights come out of the
              // table, built with the SAME expressions this used to evaluate inline --
              // the seam-parallel map atan(tan(a)*tan|n|), the half-cell rounding, the
              // clamp to the SOURCE'S ACTIVE range (not the buffer's: for a face or
              // x1-edge buffer the two coincide, but a doubly-ghost buffer is only ng
              // deep along the seam and clamping to that would extrapolate from two
              // cells when the source block holds the data), and the quadratic weights.
              const int sidx_ = (cs_seam == 2) ? 0 : 1;
              const int bs = static_cast<int>(cs_stbl(m,cs_skind,sidx_,kk,jj,0));
              const Real wm = cs_stbl(m,cs_skind,sidx_,kk,jj,1);
              const Real w0 = cs_stbl(m,cs_skind,sidx_,kk,jj,2);
              const Real wp = cs_stbl(m,cs_skind,sidx_,kk,jj,3);
              const Real sv0 = (cs_seam == 2) ? srcval(bs,jj,i)   : srcval(kk,bs,i);
              const Real sv1 = (cs_seam == 2) ? srcval(bs+1,jj,i) : srcval(kk,bs+1,i);
              const Real sv2 = (cs_seam == 2) ? srcval(bs+2,jj,i) : srcval(kk,bs+2,i);
              // MONOTONICITY LIMIT, threshold-free.
              //
              // The rule is the standard one: monotone data must give a monotone
              // interpolant.  A SHOCK is monotone across the stencil, so clamping to the
              // stencil's own range stops the overshoot that puts a negative density in a
              // ghost cell.  A smooth EXTREMUM is NOT monotone, and is left alone --
              // clipping there is the classic way a limiter destroys accuracy (measured:
              // 1.4x on the smooth tangent-seam halo).  On smooth monotone data the
              // clamp is a numerical NO-OP, because the quadratic already lies inside the
              // range, so nothing is paid for it.
              //
              // It must also be skipped when the resample is EXTRAPOLATING -- `bs` is
              // clamped at the ends of the source range, so `pos` can fall outside
              // [bs, bs+2], where the correct value legitimately lies outside the node
              // range; clamping there cost 3.5x on the smooth seam halo.
              //
              // TWO EARLIER ATTEMPTS FAILED and are recorded so they are not retried: an
              // UNCONDITIONAL clamp (cost 3.5x, the extrapolation case), and a
              // second-difference roughness test gated on a RELATIVE span
              // (hi-lo) > 0.1*(|hi|+|lo|), which fires spuriously wherever the stencil
              // straddles ZERO -- it triggered ten million times on a smooth run.
              const Real qq = wm*sv0 + w0*sv1 + wp*sv2;
              const Real dd1 = sv1 - sv0, dd2 = sv2 - sv1;
              const Real fpos = cs_stbl(m,cs_skind,sidx_,kk,jj,4);
              if (!cs_nolim_ && dd1*dd2 > 0.0 && fpos >= 0.0 && fpos <= 2.0) {
                const Real lo0 = fmin(sv0, fmin(sv1, sv2));
                const Real hi0 = fmax(sv0, fmax(sv1, sv2));
                return fmin(hi0, fmax(lo0, qq));
              }
              return qq;
            };

          // copy field components directly into recv buffer if MeshBlocks on same rank
          // seamval() reads b0 or coarse_b0 according to cs_coar, so one expression
          // serves a neighbour at ANY level.
          if (nghbr.d_view(m,n).rank == my_rank) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
            [&](const int idx) {
              int k = (idx)/nji;
              int j = (idx - k*nji)/ni;
              int i = (idx - k*nji - j*ni) + il;
              k += kl;
              j += jl;
              int kk = ak*k + bk;
              int jj = aj*j + bj;
              rbuf[dn].vars(dm, ndat*v + i-il + ni*(sj*(j-jl) + sk*(k-kl)))
                = seamval(kk,jj,i);
            });

          // else copy field components into send buffer for MPI communication below
          } else {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
            [&](const int idx) {
              int k = (idx)/nji;
              int j = (idx - k*nji)/ni;
              int i = (idx - k*nji - j*ni) + il;
              k += kl;
              j += jl;
              int kk = ak*k + bk;
              int jj = aj*j + bj;
              sbuf[n].vars(m, ndat*v + i-il + ni*(sj*(j-jl) + sk*(k-kl)))
                = seamval(kk,jj,i);
            });
          }
          } else {   // normal boundary exchange
        // copy field components directly into recv buffer if MeshBlocks on same rank
        if (nghbr.d_view(m,n).rank == my_rank) {
          // if neighbor is at same or finer level, load data from b0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
            [&](const int idx) {
              int k = (idx)/nji;
              int j = (idx - k*nji)/ni;
              int i = (idx - k*nji - j*ni) + il;
              k += kl;
              j += jl;
              if (v==0) {
                rbuf[dn].vars(dm,i-il + ni*(j-jl + nj*(k-kl))) = b.x1f(m,k,j,i);
              } else if (v==1) {
                rbuf[dn].vars(dm,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = b.x2f(m,k,j,i);
              } else if (v==2) {
                rbuf[dn].vars(dm,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = b.x3f(m,k,j,i);
              }
            });
          // if neighbor is at coarser level, load data from coarse_b0
          } else {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
            [&](const int idx) {
              int k = (idx)/nji;
              int j = (idx - k*nji)/ni;
              int i = (idx - k*nji - j*ni) + il;
              k += kl;
              j += jl;
              if (v==0) {
                rbuf[dn].vars(dm,i-il + ni*(j-jl + nj*(k-kl))) = cb.x1f(m,k,j,i);
              } else if (v==1) {
                rbuf[dn].vars(dm,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = cb.x2f(m,k,j,i);
              } else if (v==2) {
                rbuf[dn].vars(dm,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = cb.x3f(m,k,j,i);
              }
            });
          }

        // else copy field components into send buffer for MPI communication below
        } else {
          // if neighbor is at same or finer level, load data from b0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
            [&](const int idx) {
              int k = (idx)/nji;
              int j = (idx - k*nji)/ni;
              int i = (idx - k*nji - j*ni) + il;
              k += kl;
              j += jl;
              if (v==0) {
                sbuf[n].vars(m,i-il + ni*(j-jl + nj*(k-kl))) = b.x1f(m,k,j,i);
              } else if (v==1) {
                sbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = b.x2f(m,k,j,i);
              } else if (v==2) {
                sbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = b.x3f(m,k,j,i);
              }
            });
          // if neighbor is at coarser level, load data from coarse_b0
          } else {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
            [&](const int idx) {
              int k = (idx)/nji;
              int j = (idx - k*nji)/ni;
              int i = (idx - k*nji - j*ni) + il;
              k += kl;
              j += jl;
              if (v==0) {
                sbuf[n].vars(m,i-il + ni*(j-jl + nj*(k-kl))) = cb.x1f(m,k,j,i);
              } else if (v==1) {
                sbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = cb.x2f(m,k,j,i);
              } else if (v==2) {
                sbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl))) = cb.x3f(m,k,j,i);
              }
            });
          }
        }
          } // end if-do-cs/do-pole block
      } // end if-neighbor-exists block
    }
  }); // end par_for_outer
  }

#if MPI_PARALLEL_ENABLED
  // Send boundary buffer to neighboring MeshBlocks using MPI
  // fence only before the first MPI_Isend: an exchange with no off-rank neighbour
  // (1 rank) has nothing to wait for -- the unpack kernel is on the same stream.
  bool fenced_ = false;
  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0 &&
          !(use_cs && IsCubeVertexCorner(nghbr.h_view, mbpanel.h_view, m, n))) {
        // index and rank of destination Neighbor
        int dn = nghbr.h_view(m,n).dest;
        int drank = nghbr.h_view(m,n).rank;
        if (drank != my_rank) {
          if (!fenced_) {
            Kokkos::fence();
            fenced_ = true;
          }
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int lid = nghbr.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);

          // get ptr to send buffer when neighbor is at coarser/same/fine level
          int data_size = 3;
          if ( nghbr.h_view(m,n).lev < pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= sendbuf[n].icoar_ndat;
          } else if ( nghbr.h_view(m,n).lev == pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= sendbuf[n].isame_ndat;
          } else {
            data_size *= sendbuf[n].ifine_ndat;
          }
          auto send_ptr = Kokkos::subview(sendbuf[n].vars, m, Kokkos::ALL);

          int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_vars, &(sendbuf[n].vars_req[m]));
          if (ierr != MPI_SUCCESS) {no_errors=false;}
        }
      }
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}


//----------------------------------------------------------------------------------------
// \!fn void MeshBoundaryValuesFC::FillPanelCornersFC()
// \brief Fill the ng x ng corner ghost blocks of a cubed-sphere panel corner.
//
// A panel corner is a CUBE VERTEX, where only THREE panels meet. There is no fourth
// block diagonally across it, so the generic corner buffer -- which assumes a rectangular
// index region of a single neighbour and applies the face seam's signed permutation to it
// -- reaches somewhere meaningless. Measured on the uniform-Cartesian-field gate: the
// corner ghosts were wrong by ~0.7 of |B| on EVERY component, the radial one included and
// flat in resolution, i.e. not a basis error but garbage. That is what drove the residual
// one-step velocity after the face halo was fixed: the max was in a corner cell on all
// six panels.
//
// The panel's own gnomonic map is perfectly well defined out there, though -- the ghost
// zone IS that map extended -- and the two face halos flanking the corner are now
// accurate. So each corner ghost is extrapolated quadratically from the face halo on one
// side, again from the one on the other side, and the two are averaged. Quadratic keeps
// the O(dx^3) local error the face halo already has, and the stencil only ever reads the
// two flanking face-halo strips, never the corner block being written.
//
// EACH EXTRAPOLANT GETS A SIGN-PRESERVING FALLBACK, the same rule as in
// FillPanelCornersCC and applied to every field component: if the three stencil nodes all
// share a strict sign and the extrapolant comes out with the opposite sign, or exactly
// zero, it is replaced by the node adjacent to the ghost (constant extrapolation);
// otherwise it is kept bit for bit. The raw quadratic weights are not positivity
// preserving -- with nghost = 3 the third ghost carries (10, -15, 6) -- which on the
// cell-centred side drove the red-giant density extrapolant NEGATIVE above the
// star/corona join. Because a smooth profile never trips the sign test, nothing changes
// there and the corner ghost keeps its order; the earlier [min,max] node clamp did not
// have that property. See the numbers quoted above FillPanelCornersCC.
//
// Called once at the end of RecvAndUnpackFC, when every face buffer is guaranteed
// unpacked. Same-panel corners are left entirely alone.
//
// THE RADIAL GHOST LAYERS ARE PART OF THIS CORNER BLOCK. The i loop below runs over the
// whole array, ghost zones included, not just the active radial range -- exactly as the
// cell-centred FillPanelCornersCC always has. With ONE MeshBlock spanning the radial
// direction the difference is invisible: x1 ends on a physical boundary at both ends, the
// corner slot has no neighbour at all, and the physical BC fills those cells afterwards.
// Split the radial direction and the (radial ghost) x (cube-vertex corner) block is a
// real ghost region that nothing else writes -- the generic 3D corner buffer that would
// have served it is skipped as non-reciprocal (see IsCubeVertexCorner in bvals.hpp) --
// so it kept whatever stale values it happened to hold, ~0.2 of |B| off the exact
// solution on the uniform-field gate against 6e-4 for the same block unsplit. The
// flanking strips the extrapolation reads there are the x1x2 and x3x1 EDGE halos, which
// are filled by their own buffers and were measured accurate. Ideal MHD never reaches
// this block, but the resistive curl does: it is why cubed-sphere resistivity did not
// reproduce its Linf under a plain radial split.

//----------------------------------------------------------------------------------------
//! \fn SignPreserveCS
//! \brief Sign-preserving fallback for a one-sided quadratic cube-vertex extrapolant.
//! q is the extrapolant and (n0,n1,n2) the three stencil nodes it was built from, with n0
//! the node ADJACENT to the ghost being filled. If n0, n1 and n2 all share a strict sign
//! and q has the opposite sign or is exactly zero, q is discarded for n0, i.e. constant
//! extrapolation; otherwise q is returned unchanged. See the note above the caller.

KOKKOS_INLINE_FUNCTION
static Real SignPreserveCS(const Real q, const Real n0, const Real n1, const Real n2) {
  if (((n0 > 0.0) && (n1 > 0.0) && (n2 > 0.0) && !(q > 0.0)) ||
      ((n0 < 0.0) && (n1 < 0.0) && (n2 < 0.0) && !(q < 0.0))) {
    return n0;
  }
  return q;
}

namespace cs_vertex_diag { extern Real fc_perturb, perturb_add; }   // cornerfix-1103

void MeshBoundaryValuesFC::FillPanelCornersFC(DvceFaceFld4D<Real> &b, bool coarse) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = coarse ? indcs.cjs : indcs.js;
  const int je = coarse ? indcs.cje : indcs.je;
  const int ks = coarse ? indcs.cks : indcs.ks;
  const int ke = coarse ? indcs.cke : indcs.ke;
  const int ng = indcs.ng;
  const int nmb = pmy_pack->nmb_thispack;
  const int nc1 = (coarse ? indcs.cnx1 : indcs.nx1) + 2*ng;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto b_ = b;
  const bool vfill_ = pmy_pack->pmesh->cs_vertex_fill;
  const bool coar_ = coarse;
  const int is_ = indcs.is, ie_ = indcs.ie;
  auto &mblev_ = pmy_pack->pmb->mb_lev;
  // cornerfix-1103 DIAGNOSTIC <mesh>/cs_vertex_fc_perturb / cs_vertex_perturb_add (0 unless
  // named): every corner-block face, sampled or extrapolated, -> b*(1 + fp) + fa
  const Real fp_ = cs_vertex_diag::fc_perturb, fa_ = cs_vertex_diag::perturb_add;
  const bool fpon_ = (fp_ != 0.0) || (fa_ != 0.0);

  par_for("cs_fill_corners", DevExeSpace(), 0,(nmb-1), 0,3, 0,nc1, 0,(ng-1), 0,(ng-1),
  KOKKOS_LAMBDA(const int m, const int c, const int i, const int gj, const int gk) {
    const int sj = (c & 1) ? 1 : -1;      // -1 = the -x2 side of the block
    const int sk = (c & 2) ? 1 : -1;      // -1 = the -x3 side
    // Buffer ids of the two FACE neighbours flanking this corner (see nghbr_index.hpp).
    const int nj_id = (sj < 0) ? 8 : 12;
    const int nk_id = (sk < 0) ? 24 : 28;
    const int mp = mbpanel.d_view(m);
    // ONLY A TRUE CUBE VERTEX.  The condition is that BOTH flanking faces are panel
    // seams -- exactly IsCubeVertexCorner's test, and exactly the case whose exchange is
    // skipped.  It used to fire when EITHER face was a seam, which with more than one
    // MeshBlock per panel also caught ordinary corners that DO have a real diagonal
    // neighbour, overwriting properly exchanged, seam-transformed data with a one-sided
    // extrapolation.
    bool seamj = false, seamk = false;
    if (nghbr.d_view(m,nj_id).gid >= 0 && nghbr.d_view(m,nj_id).panel != mp) seamj = true;
    if (nghbr.d_view(m,nk_id).gid >= 0 && nghbr.d_view(m,nk_id).panel != mp) seamk = true;
    if (!(seamj && seamk)) return;
    // THE CORNER FILL SUPERSEDES THIS.  When the flanking neighbours are at the SAME
    // level their widened buffers have already written this corner block with REAL data
    // sampled from the two panels that cover it (see buffs_fc.cpp), which is strictly
    // better than any extrapolation -- so do not overwrite it.  At a level boundary the
    // widening does not apply (only `isame` was widened) and the extrapolation below is
    // still the best available.
    //
    // WHICH BUFFER REACHES WHICH RADIAL LAYER.  Over the ACTIVE radial range the two
    // widened FACE buffers cover the block.  A RADIAL GHOST layer is triply ghost, so
    // only a buffer that is already ghost in x1 can reach it: the x1x2 and x3x1 EDGE
    // buffers, widened by the same rule.  Both must exist and be at this level, because
    // ownership splits the block between them along the diagonal and one alone would
    // leave half of it written by NOTHING -- the same trap as 397b4ad3, which put 0.2 of
    // |B| in the r x CORNER SEAM category.
    const bool vfill_here = (vfill_ && !coar_ &&
        nghbr.d_view(m,nj_id).lev == mblev_.d_view(m) &&
        nghbr.d_view(m,nk_id).lev == mblev_.d_view(m));
    const int s1 = (i < is_) ? -1 : 1;                 // radial side of a ghost layer
    const int e12 = 16 + (s1 + 1) + 2*(sj + 1);        // x1x2 edge (nghbr_index.hpp)
    const int e31 = 24 + (s1 + 9) + 2*(sk + 1);        // x3x1 edge
    const bool vfill_rg = vfill_here &&
        nghbr.d_view(m,e12).gid >= 0 && nghbr.d_view(m,e12).lev == mblev_.d_view(m) &&
        nghbr.d_view(m,e31).gid >= 0 && nghbr.d_view(m,e31).lev == mblev_.d_view(m);
    const bool vfill_cc = vfill_here && ((i >= is_ && i <= ie_) ? true : vfill_rg);
    const bool vfill_x1 = vfill_here && ((i >= is_ && i <= ie_ + 1) ? true : vfill_rg);

    // Quadratic Lagrange extrapolated d cells beyond an anchor, nodes 0,1,2 stepping
    // inward: w = ((d+1)(d+2)/2, -d(d+2), d(d+1)/2), which sums to 1.
    const Real dj = static_cast<Real>(gj + 1);
    const Real dk = static_cast<Real>(gk + 1);
    // QUADRATIC is the optimum of this family, measured: LINEAR (amplification 3 and 5
    // against quadratic's 7 and 17) loses to truncation -- corner 1.98e-03 -> 5.44e-03 --
    // and CUBIC (15 and 49) loses to amplification -- 1.98e-03 -> 2.46e-03.  Do not
    // re-propose either.
    const Real wj0 = 0.5*(dj+1.0)*(dj+2.0), wj1 = -dj*(dj+2.0), wj2 = 0.5*dj*(dj+1.0);
    const Real wk0 = 0.5*(dk+1.0)*(dk+2.0), wk1 = -dk*(dk+2.0), wk2 = 0.5*dk*(dk+1.0);
    const int stj = -sj, stk = -sk;       // step INWARD from the anchor
    // The two extrapolations do not reach equally far: the one along j travels dj cells
    // and the one along k travels dk.  A plain average weights them the same even when
    // one is a 1-cell reach and the other a 2-cell reach; weighting each by the OTHER's
    // distance favours the shorter, more accurate one.
    const Real wblend = dj/(dj + dk);

    // cell-index targets and anchors
    const int jtc = (sj < 0) ? (js-1-gj) : (je+1+gj);
    const int ktc = (sk < 0) ? (ks-1-gk) : (ke+1+gk);
    const int ajc = (sj < 0) ? js : je;
    const int akc = (sk < 0) ? ks : ke;
    // face-index targets and anchors (one more face than cells on the outer side)
    const int jtf = (sj < 0) ? (js-1-gj) : (je+2+gj);
    const int ktf = (sk < 0) ? (ks-1-gk) : (ke+2+gk);
    const int ajf = (sj < 0) ? js : (je+1);
    const int akf = (sk < 0) ? ks : (ke+1);

    // b.x1f: cell indices in both j and k, faces in i
    if (!vfill_x1) {
      const Real k0 = b_.x1f(m,akc,jtc,i);
      const Real k1 = b_.x1f(m,akc+stk,jtc,i);
      const Real k2 = b_.x1f(m,akc+2*stk,jtc,i);
      const Real j0 = b_.x1f(m,ktc,ajc,i);
      const Real j1 = b_.x1f(m,ktc,ajc+stj,i);
      const Real j2 = b_.x1f(m,ktc,ajc+2*stj,i);
      Real ek = wk0*k0 + wk1*k1 + wk2*k2;
      Real ej = wj0*j0 + wj1*j1 + wj2*j2;
      ek = SignPreserveCS(ek, k0, k1, k2);
      ej = SignPreserveCS(ej, j0, j1, j2);
      b_.x1f(m,ktc,jtc,i) = wblend*ek + (1.0-wblend)*ej;
    }
    if (fpon_) b_.x1f(m,ktc,jtc,i) = b_.x1f(m,ktc,jtc,i)*(1.0 + fp_) + fa_;
    if (i > (nc1-1)) return;
    // b.x2f: FACE index in j, cell index in k
    if (!vfill_cc) {
      const Real k0 = b_.x2f(m,akc,jtf,i);
      const Real k1 = b_.x2f(m,akc+stk,jtf,i);
      const Real k2 = b_.x2f(m,akc+2*stk,jtf,i);
      const Real j0 = b_.x2f(m,ktc,ajf,i);
      const Real j1 = b_.x2f(m,ktc,ajf+stj,i);
      const Real j2 = b_.x2f(m,ktc,ajf+2*stj,i);
      Real ek = wk0*k0 + wk1*k1 + wk2*k2;
      Real ej = wj0*j0 + wj1*j1 + wj2*j2;
      ek = SignPreserveCS(ek, k0, k1, k2);
      ej = SignPreserveCS(ej, j0, j1, j2);
      b_.x2f(m,ktc,jtf,i) = wblend*ek + (1.0-wblend)*ej;
    }
    if (fpon_) b_.x2f(m,ktc,jtf,i) = b_.x2f(m,ktc,jtf,i)*(1.0 + fp_) + fa_;
    // b.x3f: cell index in j, FACE index in k
    if (!vfill_cc) {
      const Real k0 = b_.x3f(m,akf,jtc,i);
      const Real k1 = b_.x3f(m,akf+stk,jtc,i);
      const Real k2 = b_.x3f(m,akf+2*stk,jtc,i);
      const Real j0 = b_.x3f(m,ktf,ajc,i);
      const Real j1 = b_.x3f(m,ktf,ajc+stj,i);
      const Real j2 = b_.x3f(m,ktf,ajc+2*stj,i);
      Real ek = wk0*k0 + wk1*k1 + wk2*k2;
      Real ej = wj0*j0 + wj1*j1 + wj2*j2;
      ek = SignPreserveCS(ek, k0, k1, k2);
      ej = SignPreserveCS(ej, j0, j1, j2);
      b_.x3f(m,ktf,jtc,i) = wblend*ek + (1.0-wblend)*ej;
    }
    if (fpon_) b_.x3f(m,ktf,jtc,i) = b_.x3f(m,ktf,jtc,i)*(1.0 + fp_) + fa_;
  });
  return;
}


//----------------------------------------------------------------------------------------
// \!fn void RecvBuffers()
// \brief Unpack boundary buffers

TaskStatus MeshBoundaryValuesFC::RecvAndUnpackFC(DvceFaceFld4D<Real> &b,
                                                 DvceFaceFld4D<Real> &cb) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto rbuf = RecvBufDv();
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  // For the cube-vertex corner fill's ownership test (see the note at its use below).
  const bool cs_vfill_ = pmy_pack->pmesh->cs_vertex_fill;
  auto &wi_ = pmy_pack->pmesh->mb_indcs;
  const int js_ = wi_.js, je_ = wi_.je, ks_ = wi_.ks, ke_ = wi_.ke;
  const int nx2_ = wi_.nx2, nx3_ = wi_.nx3;
  auto &mbsz = pmy_pack->pmb->mb_size;
#if MPI_PARALLEL_ENABLED
  //----- STEP 1: check that recv boundary buffer communications have all completed

  bool bflag = false;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0 &&
          !(use_cs && IsCubeVertexCorner(nghbr.h_view, mbpanel.h_view, m, n))) {
        if (nghbr.h_view(m,n).rank != global_variable::my_rank) {
          int test;
          int ierr = MPI_Test(&(recvbuf[n].vars_req[m]), &test, MPI_STATUS_IGNORE);
          if (ierr != MPI_SUCCESS) {no_errors=false;}
          if (!(static_cast<bool>(test))) {
            bflag = true;
          }
        }
      }
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // exit if recv boundary buffer communications have not completed
  if (bflag) {return TaskStatus::incomplete;}
#endif

  //----- STEP 2: buffers have all completed, so unpack 3-components of field

  // Outer loop over (# of MeshBlocks)*(three field components)*(radial chunks)
  // r1mspeed-1102: the radial index range is cut into chunks of FC_UNP_CHUNK, one team
  // each.  Buffers overlap only at the SAME (k,j,i), which stays in one team, so each
  // element still takes its writes in the buffer order n (bitwise); 72 teams were
  // ~half the device.
  constexpr int FC_UNP_CHUNK = 16;
  const int nic = (pmy_pack->pmesh->mb_indcs.nx1 + 2*pmy_pack->pmesh->mb_indcs.ng + 1
                   + FC_UNP_CHUNK - 1)/FC_UNP_CHUNK;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb*nic), Kokkos::AUTO);
  BvalsTeamFor("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/(3*nic);
    const int v = (tmember.league_rank()/nic)%3;
    const int ic0 = (tmember.league_rank()%nic)*FC_UNP_CHUNK;
    const int ic1 = ic0 + FC_UNP_CHUNK - 1;

    // scalar loop over neighbors to prevent race condition in overlapping assignments
    for (int n=0; n<nnghbr; ++n) {
      // only unpack buffers when neighbor exists (cube vertex skipped -- see bvals.hpp)
      if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
        // if neighbor is at coarser level, use cindices to unpack buffer
        int il, iu, jl, ju, kl, ku, ndat;
        if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
          il = rbuf[n].icoar[v].bis;
          iu = rbuf[n].icoar[v].bie;
          jl = rbuf[n].icoar[v].bjs;
          ju = rbuf[n].icoar[v].bje;
          kl = rbuf[n].icoar[v].bks;
          ku = rbuf[n].icoar[v].bke;
          ndat = rbuf[n].icoar_ndat;
        // if neighbor is at same level, use sindices to unpack buffer
        } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
          il = rbuf[n].isame[v].bis;
          iu = rbuf[n].isame[v].bie;
          jl = rbuf[n].isame[v].bjs;
          ju = rbuf[n].isame[v].bje;
          kl = rbuf[n].isame[v].bks;
          ku = rbuf[n].isame[v].bke;
          ndat = rbuf[n].isame_ndat;
        // if neighbor is at finer level, use findices to unpack buffer
        } else {
          il = rbuf[n].ifine[v].bis;
          iu = rbuf[n].ifine[v].bie;
          jl = rbuf[n].ifine[v].bjs;
          ju = rbuf[n].ifine[v].bje;
          kl = rbuf[n].ifine[v].bks;
          ku = rbuf[n].ifine[v].bke;
          ndat = rbuf[n].ifine_ndat;
        }
        const int ni = iu - il + 1;
        const int nj = ju - jl + 1;
        // r1mspeed-1102: this team's radial slice [cl, cu] of [il, iu]
        const int cl = (il > ic0) ? il : ic0;
        const int cu = (iu < ic1) ? iu : ic1;
        const int nci = (cu >= cl) ? (cu - cl + 1) : 0;
        const int nk = ku - kl + 1;
        const int nkji = nk*nj*nci;
        const int nji  = nj*nci;

        // if neighbor is at same or finer level, load data directly into b0
        if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
          [&](const int idx) {
            int k = (idx)/nji;
            int j = (idx - k*nji)/nci;
            int i = (idx - k*nji - j*nci) + cl;
            k += kl;
            j += jl;
            // CUBE-VERTEX CORNER FILL: OWNERSHIP.  The face buffers were widened by ng
            // along the seam (see buffs_fc.cpp) so the doubly-ghost corner block rides
            // the face exchange.  Both flanking face buffers reach that block, but each
            // panel only genuinely COVERS the half of it on its own side of the diagonal
            // -- so write a doubly-ghost cell only from the neighbour whose panel
            // actually contains that direction.  Decided here, from the RECEIVER's own
            // geometry, which needs nothing from the sender.  Cells outside the corner
            // block are unaffected, and a non-cubed-sphere run never enters this.
            // ONLY the slots that were actually widened (buffs_fc.cpp): the x2 faces
            // 8-15 and the x1x2 edges 16-23, whose seam normal is x2, and the x3 faces
            // 24-31 and the x3x1 edges 32-39, whose seam normal is x3.  The edges are the
            // ones that carry the fill into the corner's RADIAL GHOST layers.
            // This test MUST NOT reach the x2x3 edge (40-47) or corner (48-55) buffers:
            // their destinations are legitimately doubly-ghost, and gating them here
            // SUPPRESSES the real edge and corner exchange -- which is exactly the bug
            // that made a widening of ZERO cells still destroy the solution.
            const bool face_wide = (n >= 8 && n < 40);
            // STAGGERING.  b.x2f lives on x2 FACES and b.x3f on x3 faces, so for those
            // components je+1 / ke+1 is a REAL interior face, not a ghost.  Testing the
            // cell bounds alone misclassifies it and suppresses a legitimate write --
            // which is why a widening of ZERO cells still destroyed the solution.
            const int jghi = je_ + ((v == 1) ? 1 : 0);
            const int kghi = ke_ + ((v == 2) ? 1 : 0);
            if (cs_vfill_ && face_wide &&
                (j < js_ || j > jghi) && (k < ks_ || k > kghi)) {
              // ONLY a CUBE VERTEX needs the fill.  At an ordinary corner the diagonal
              // neighbour exists and its own buffer (40..55, unpacked LATER in this same
              // scalar n loop) is the correct filler; the panel-edge geometry used below
              // is not even valid there, because the block edge is not the panel edge.
              const int sj = (j < js_) ? 0 : 1;
              const int sk = (k < ks_) ? 0 : 1;
              if (!IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m,
                                      40 + 2*(sj + 2*sk))) return;
              // OWNERSHIP, exactly.  The corner ghost block is covered by the two
              // flanking panels, split by the seam BETWEEN THEM, and in this panel's own
              // gnomonic chart that seam is exactly the DIAGONAL.  Take P0 = +z, so
              // tan(xi) = x/z and tan(eta) = y/z; the +x/+y panel boundary is the plane
              // x = y, i.e. xi = eta.  So with a = s2*xi - pi/4 and b = s3*eta - pi/4
              // measuring how far the cell lies beyond each panel edge, a > b means the
              // direction belongs to the panel across the x2 face and a < b to the one
              // across the x3 face.  The other seven vertices follow by reflection, which
              // is what the s2, s3 signs are.  The test is COMPLEMENTARY by construction:
              // one of the two flanking buffers writes every cell of the block, and
              // neither writes a cell the other owns.  That matters because
              // FillPanelCornersFC stops extrapolating the whole block once the fill is
              // on -- a predicate with a gap would leave cells written by NOTHING.
              //
              // STAGGERING: b.x2f lives on x2 faces and b.x3f on x3 faces, so the
              // component being written decides whether its own index is a face or a
              // cell centre.  Half a cell is not negligible here: the two panels meet
              // AT this diagonal.
              const Real xi  = 0.25*M_PI*((v == 1)
                             ? LeftEdgeX(j-js_, nx2_, mbsz.d_view(m).x2min,
                                                      mbsz.d_view(m).x2max)
                             : CellCenterX(j-js_, nx2_, mbsz.d_view(m).x2min,
                                                        mbsz.d_view(m).x2max));
              const Real eta = 0.25*M_PI*((v == 2)
                             ? LeftEdgeX(k-ks_, nx3_, mbsz.d_view(m).x3min,
                                                      mbsz.d_view(m).x3max)
                             : CellCenterX(k-ks_, nx3_, mbsz.d_view(m).x3min,
                                                        mbsz.d_view(m).x3max));
              const Real a = (sj ? xi : -xi) - 0.25*M_PI;
              const Real b = (sk ? eta : -eta) - 0.25*M_PI;
              if ((n < 24) != (a >= b)) return;
            }
            if (v==0) {
              b.x1f(m,k,j,i) = rbuf[n].vars(m,i-il + ni*(j-jl + nj*(k-kl)));
            } else if (v==1) {
              b.x2f(m,k,j,i) = rbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl)));
            } else if (v==2) {
              b.x3f(m,k,j,i) = rbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl)));
            }
          });
        // if neighbor is at coarser level, load data into coarse_b0
        } else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
          [&](const int idx) {
            int k = (idx)/nji;
            int j = (idx - k*nji)/nci;
            int i = (idx - k*nji - j*nci) + cl;
            k += kl;
            j += jl;
            if (v==0) {
              cb.x1f(m,k,j,i) = rbuf[n].vars(m,i-il + ni*(j-jl + nj*(k-kl)));
            } else if (v==1) {
              cb.x2f(m,k,j,i) = rbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl)));
            } else if (v==2) {
              cb.x3f(m,k,j,i) = rbuf[n].vars(m,ndat*v + i-il + ni*(j-jl + nj*(k-kl)));
            }
          });
        }
        tmember.team_barrier();
      }  // end if-neighbor-exists block
    }
  });  // end par_for_outer

  // Every face buffer is unpacked by this point, which is what the corner fill needs.
  //
  // THE COARSE ARRAY NEEDS THE SAME FILL.  cb's cube-vertex corner is written by NOTHING:
  // its exchange is skipped as non-reciprocal exactly like the fine array's, but the fill
  // only ever ran on b.  ProlongateFC then interpolates that block into the fine ghosts.
  // A poison test -- setting cb's cube-vertex corner to 1e30 -- moves the prolongated
  // fine EDGE halo by 24% and the fine corner likewise, so it IS read; nothing blew up
  // only because the prolongation limiter clips it.  Measured on the refined iprob=11
  // halo, cb's corner goes 2.09e-01 (i.e. |B| -- garbage) to 1.86e-03, and the fine halo
  // follows: r x EDGE SEAM 2.3278e-03 -> 4.6508e-04, r x CORNER SEAM 1.8152e-02 ->
  // 2.0911e-03.  Both then land at the SAME-LEVEL halo's own accuracy, which is the
  // statement that a level boundary no longer degrades the seam halo.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    FillPanelCornersFC(b, false);
    if (pmy_pack->pmesh->multilevel) { FillPanelCornersFC(cb, true); }
  }

  // layoutbug-1103: <mesh>/cs_seam_blockend = extend, the second pass (see the note in
  // RecvAndUnpackCC): repeat the exchange with the panel-clamped stencil table, at the
  // depth of the pack it repeats.
  if (cs_blockend_ext_ && !cs_ext_pass_ && pmy_pack->pmesh->use_cubed_sphere) {
    if (cs_ext_q_) CsSeamExtPartner(b);   // seamexact-1103
    cs_ext_pass_ = true;
    const int pd0 = pack_depth;
    pack_depth = last_pack_depth_;
    (void) ClearSend();
    (void) InitRecv(3);
    while (PackAndSendFC(b, cb) != TaskStatus::complete) {}
    while (RecvAndUnpackFC(b, cb) != TaskStatus::complete) {}
    // seamexact-1103: wait on this pass's own sends here.  A caller that cleared its
    // sends before unpacking (rad_m1_coupling.cpp FluidRefresh) would otherwise leave
    // them in flight into its next exchange, which reuses the buffers and the requests.
    (void) ClearSend();
    pack_depth = pd0;
    cs_ext_pass_ = false;
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesFC::CsSeamExtPartner
//! \brief seamexact-1103: <mesh>/cs_seam_blockend = extend with > 2 MeshBlocks per panel
//! side.  The extended pass's partner interpolation (x3f at an x2-face centre, x2f at an
//! x3-face centre: a 4-point Lagrange in the along-seam direction of a 4-point Lagrange
//! across it) is centred only up to ng-1 cells into the sender's along-seam ghosts; at
//! the outermost ghost it would need face ng+2 beyond the block end and is one-sided,
//! where one block per panel interpolates.  The resample reaches that ghost wherever the
//! block end is off the panel midline (the seam map shifts the source position by up to
//! ~a cell there), so 4x4 and finer layouts differed from 1x1 at truncation level
//! (~1e-7 in B per exchange at nghost = 2; nghost = 3 hid it).  Here every block
//! evaluates the interpolants on its ACTIVE faces with the extended pass's own bounds,
//! and a same-panel ghost exchange (plain scalars) hands each block its neighbours'
//! values at the along-seam ghost positions; the extended pass reads them there
//! (FCSeamPack::q_partner): the very numbers one block per panel computes.  Called after
//! the first pass, so every face it reads is in place.

void MeshBoundaryValuesFC::CsSeamExtPartner(DvceFaceFld4D<Real> &b) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ncell = nmb*(id.nx3 + 2*id.ng)*(id.nx2 + 2*id.ng)*(id.nx1 + 2*id.ng);
  FCSeamPack f(b, b, SendBufDv(), RecvBufDv(), fc_seg, fc_sub[0], fc_cs[0], fc_cb[0],
               cs_seam_geom, cs_seam_stnc_x, fc_preval, 0, id, 0.0, false);
  f.ext_ = true; f.ext_jk = cs_ext_jk; f.q = cs_ext_q;
  Kokkos::parallel_for("CsSeamExtPartner",
      Kokkos::RangePolicy<FCSeamQ>(DevExeSpace(), 0, ncell), f);
  DvceArray5D<Real> q = cs_ext_q;
  while (pbq_->InitRecv(4) == TaskStatus::incomplete) {}
  while (pbq_->PackAndSendCC(q, q) == TaskStatus::incomplete) {}
  while (pbq_->RecvAndUnpackCC(q, q) == TaskStatus::incomplete) {}
  while (pbq_->ClearSend() == TaskStatus::incomplete) {}
}
