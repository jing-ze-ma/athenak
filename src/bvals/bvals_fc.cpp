//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for face-centered (FC)
//! Mesh variables.
//! Prolongation of FC variables occurs in ProlongateFC function called from task list

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <iomanip>    // std::setprecision

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"
#include "bvals/cs_seam_fc3.hpp"
#include "bvals/cs_seam_round.hpp"

//----------------------------------------------------------------------------------------
// BValFC constructor:

MeshBoundaryValuesFC::MeshBoundaryValuesFC(MeshBlockPack *pp, ParameterInput *pin) :
  MeshBoundaryValues(pp, pin, false) {
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesFC::BuildCubedSphereSeamTables
//! \brief Precomputes the (xi,eta)-dependent geometry and stencil tables of the
//! cubed-sphere seam packer, which do not depend on the radial index.

void MeshBoundaryValuesFC::BuildCubedSphereSeamTables() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  // TWO KIND AXES. The transform is needed only for the two ANGULAR components, the
  // resample for all three -- b.x1f crosses a seam as a scalar but is still resampled
  // along it -- so the geometry table carries 2 staggerings and the stencil table 3,
  // doubled when the mesh is multilevel and cross-level seam buffers exist.
  const int nkg = ml_ ? 4 : 2;
  const int nks = ml_ ? 6 : 3;
  // The tables are indexed by the ABSOLUTE (kk,jj) of the sample, GHOSTS INCLUDED. On a
  // multilevel mesh a same-level send buffer is widened by ng in its transverse
  // directions -- that is what feeds the neighbour's prolongation -- so the packer reads
  // source cells outside the active range, and indexing the table from the active corner
  // ran it off the end. One entry per face index of the block array covers everything.
  const int nk2 = indcs.nx2 + 2*indcs.ng + 1;
  const int nk3 = indcs.nx3 + 2*indcs.ng + 1;
  // <mesh>/cs_seam_fc_order = 3 needs one level, ng >= 2 and >= 4 cells
  // per block side (3 cells / 4 faces next to the seam, 3-4 along it)
  if (pmy_pack->pmesh->cs_seam_fc_order == 3 &&
      (ml_ || indcs.ng < 2 || indcs.nx2 < 4 || indcs.nx3 < 4)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mesh>/cs_seam_fc_order = 3 needs a single-level cubed sphere, "
              << "nghost >= 2 and MeshBlocks of >= 4 cells in x2 and x3" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Which destination panels does each block actually send to? At most four (the four
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
  // modify_host/sync_device, not the explicit-space templates: those fail the
  // DualView static_assert under hipcc (same trap as cdd7d2a5, and the same fix the
  // radiation tetrads already carry).
  cs_seam_slot.modify_host();
  cs_seam_slot.sync_device();

  Kokkos::realloc(cs_seam_geom, nmb, nslot, nkg, nk3, nk2, 15);
  Kokkos::realloc(cs_seam_stnc, nmb, nks, 2, nk3, nk2, 7);
  // the stencil of the second (extended) pass, <mesh>/cs_seam_blockend
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
    // with > 2 blocks per panel side a block end lies off the panel
    // midline, the along-seam resample shifts there and reaches the outermost along-seam
    // ghost, whose partner interpolation needs one face more than this block holds: those
    // interpolants come from the owning same-panel neighbour (CsSeamExtPartner).
    auto &pmi = pmy_pack->pmesh->mesh_indcs;
    cs_ext_q_ = (pmi.nx2 > 2*indcs.nx2) || (pmi.nx3 > 2*indcs.nx3);
    // order 3 keeps every stencil within ng of the sample (srcval3), so
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
    // reads, which is what fixes whether xi/eta are a face or a centre. vv = 0 is a real
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
    // orientations (cs_seam = 2 stored at sidx 0, cs_seam = 3 at sidx 1). Neither
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
      // offset of this block in the panel, for a layout-independent seam position
      const Real wp_ = (sidx == 0) ? (x3mx - x3mn) : (x2mx - x2mn);
      const Real xp_ = (sidx == 0) ? x3mn : x2mn;
      const int spn_ = static_cast<int>(floor((xp_ + 1.0)/wp_ + 0.5))
                       *((sidx == 0) ? nx3_ : nx2_);
      const Real pos = cs_seam::SeamPos(sc, spn_,
                                        (atan(tan(ang)*tan(fabs(nrm))) - ang)/dang);
      // TIE-SAFE stencil choice. At a panel END face (ang = +-pi/4, the
      // along-staggered component) the map shifts by exactly the normal offset, so at a
      // ghost cell centre pos is a HALF-INTEGER mathematically, and round(pos) went
      // either way with the last bit of pos -- which depends on the block's own
      // coordinates. That picked different 3-point stencils in different layouts (1x1 vs
      // 4x4: 4e-4 in x2f(ks-2,js) / x3f(ks,js-2) at 12 cube vertices). cs_seam::SeamRound
      // (cs_seam_round.hpp) breaks the tie the same way in every layout; any pos not
      // within 1e-9 of a tie is unchanged.
      int bs = cs_seam::SeamRound(pos) - 1;
      bs = (bs < blo) ? blo : ((bs > bhi) ? bhi : bs);
      const Real u = pos - static_cast<Real>(bs + 1);
      stnc_(m,kd,sidx,kk,jj,0) = static_cast<Real>(bs);
      stnc_(m,kd,sidx,kk,jj,1) = 0.5*u*(u - 1.0);
      stnc_(m,kd,sidx,kk,jj,2) = 1.0 - u*u;
      stnc_(m,kd,sidx,kk,jj,3) = 0.5*u*(u + 1.0);
      stnc_(m,kd,sidx,kk,jj,4) = pos - static_cast<Real>(bs);
      // SHEAR CORRECTION geometry. `along_eta` is (cs_seam == 2), i.e. sidx == 0, and
      // its (ang_,nrm_) are exactly the (ang,nrm) the resample already picked, so one
      // pair serves both. sigma is stored unconditionally; the packer keeps the
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
          if (ib_ > 0) {
            blx -= ngx_;
          }
          if (ib_ < nbp_ - 1) {
            bhx += ngx_;
          }
        }
        int bx = cs_seam::SeamRound(pos) - 1;   // tie-safe, as bs
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
//! \fn void MeshBoundaryValuesFC::BuildPackSegmentsFC
//! \brief Resolves the per-(m,n,v) setup of the cubed-sphere face-field pack once into
//! segment tables, so the pack kernel does not recompute it per call.

namespace {
constexpr int FC_PACK_CHUNK = 128;
constexpr int kFc3Max = 48;   // max distinct faces of one order-3 sample
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
          if (cs_seam != 0 && seam_extent < 3) {
            cs_seam = 0;
          }
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
  for (int c=0; c<nch; ++c) {
    h_cs(c) = cseg[c];
    h_cb(c) = cbeg[c];
  }
  Kokkos::deep_copy(fc_seg, h_seg);
  Kokkos::deep_copy(fc_chunk_seg, h_cs);
  Kokkos::deep_copy(fc_chunk_beg, h_cb);

  // the sample boxes, chunk lists and seam srcval boxes
  {
    const int js_ = cs_indcs.js, je_ = cs_indcs.je, ks_ = cs_indcs.ks, ke_ = cs_indcs.ke;
    // order 3 evaluates every seam sample itself (srcval3, no resample box)
    const bool fc3_ = (pmy_pack->pmesh->cs_seam_fc_order == 3);
    std::vector<SubSegFC> sub;
    std::vector<int> cs_p, cb_p, cs_s, cb_s, cs_r, cb_r;
    int64_t pv = 0;
    for (int q=0; q<nseg; ++q) {
      const PackSegFC &g = h_seg(q);
      SubSegFC u;
      u.q = q; u.il = g.il; u.ni = g.ni; u.jl = g.jl; u.nj = g.nj; u.kl = g.kl;
      u.nk = g.nk;
      u.pvoff = -1; u.ptlo = 0; u.pnt = 0; u.pnlo = 0; u.pnn = 0; u.pil = 0; u.pni = 0;
      if (g.mode == 2 && g.cs_seam != 0 && !g.cs_coar && !fc3_) {
        // the samples' resample reads srcval at (bs..bs+2) along the seam, bs within
        // [lo, hi-2] (BuildCubedSphereSeamTables), at the sample's own normal index
        u.pil = u.il; u.pni = u.ni;
        // with cs_seam_blockend = extend the box also covers ng cells
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
        for (int c0=0; c0<npts; c0+=FC_PACK_CHUNK) {
          cs_r.push_back(us);
          cb_r.push_back(c0);
        }
        pv += npts;
      } else if (g.mode == 2) {
        // every other seam / pole segment: srcval at each sample itself, indexed by the
        // sample (pnt = 0), so the sample pass never evaluates srcval
        u.pvoff = static_cast<int>(pv); u.pnt = 0;
        const int npts = u.nk*u.nj*u.ni;
        const int us = static_cast<int>(sub.size());
        for (int c0=0; c0<npts; c0+=FC_PACK_CHUNK) {
          cs_r.push_back(us);
          cb_r.push_back(c0);
        }
        pv += npts;
      }
      const int us = static_cast<int>(sub.size());
      const int nkji = u.nk*u.nj*u.ni;
      for (int c0=0; c0<nkji; c0+=FC_PACK_CHUNK) {
        if (g.mode != 2) {
          cs_p.push_back(us);
          cb_p.push_back(c0);
        } else {
          cs_s.push_back(us);
          cb_s.push_back(c0);
        }
      }
      sub.push_back(u);
    }
    fc_np = static_cast<int>(cs_p.size());
    fc_ns = static_cast<int>(cs_s.size());
    fc_npre = static_cast<int>(cs_r.size());
    // chunk lists: plain, then seam, then precompute
    std::vector<int> ca(cs_p), cb(cb_p);
    ca.insert(ca.end(), cs_s.begin(), cs_s.end());
    cb.insert(cb.end(), cb_s.begin(), cb_s.end());
    ca.insert(ca.end(), cs_r.begin(), cs_r.end());
    cb.insert(cb.end(), cb_r.begin(), cb_r.end());
    const int nsub = static_cast<int>(sub.size());
    Kokkos::realloc(fc_sub, (nsub > 0) ? nsub : 1);
    auto h_sub = Kokkos::create_mirror_view(fc_sub);
    for (int u=0; u<nsub; ++u) {
      h_sub(u) = sub[u];
    }
    Kokkos::deep_copy(fc_sub, h_sub);
    const int nc = static_cast<int>(ca.size());
    Kokkos::realloc(fc_cs, (nc > 0) ? nc : 1);
    Kokkos::realloc(fc_cb, (nc > 0) ? nc : 1);
    auto h_a = Kokkos::create_mirror_view(fc_cs);
    auto h_b = Kokkos::create_mirror_view(fc_cb);
    for (int c=0; c<nc; ++c) {
      h_a(c) = ca[c];
      h_b(c) = cb[c];
    }
    Kokkos::deep_copy(fc_cs, h_a);
    Kokkos::deep_copy(fc_cb, h_b);
    Kokkos::realloc(fc_preval, (pv > 0) ? pv : 1);
  }
  fc_seg_ready = true;
}

namespace {
//----------------------------------------------------------------------------------------
//! \struct FCPlainPack
//! \brief one thread per sample of a plain (same-panel) face-field buffer.

struct FCPlainPack {
  DvceFaceFld4D<Real> b, cb;
  DvceArray1D<MeshBufferDv> sbuf, rbuf;
  DvceArray1D<PackSegFC> seg;
  DvceArray1D<SubSegFC> sub;   // chunk -> sample box
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
//! \brief one thread per sample of a cross-panel (or polar) face-field
//! buffer. The member functions are the old kernel's lambdas (bx1..bx3, lag4, the four
//! partner interpolations, srcval, seamval) with every captured setup variable passed BY
//! VALUE in a SeamArgs; each expression is copied unchanged in form and order. See the
//! long comments at the old kernel (PackAndSendFC, non-cubed-sphere path) for the
//! physics of every step.

struct SeamArgs {
  int m, v, vv, signvar, cs_xform, cs_coar, cs_slotidx, cs_gkind, cs_skind, cs_seam;
  int js_, ks_, je_, ke_;
  bool wide2, wide3;
  int pvoff, ptlo, pnt, pnlo, pil, pni;   // (SubSegFC); pvoff < 0: none
};

struct FCSeamPre {};   // tag of the srcval precompute pass
struct FCSeamQ {};     // tag of the partner-interpolant pass
struct FCSeamPre3 {};  // tag of the order-3 precompute pass
struct FCSeamTab {};   // tag of the order-3 coefficient-table build

// the order-3 coefficient tables, one row per (seam segment, sample k,j),
// [0] with the first pass's bounds, [1] with the extended pass's (cs_seam_blockend =
// extend). Kept here, keyed by the boundary object, so bvals.hpp is untouched.
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
  DvceArray1D<Real> pre;   // srcval scratch (SubSegFC pv* boxes)
  int c0;
  RegionIndcs cs_indcs;
  bool cs_nolim_;
  // second pass of <mesh>/cs_seam_blockend = extend -- the source index
  // bounds are widened by (j lo, j hi, k lo, k hi) = ext_jk(m,:) at interior block ends
  bool ext_ = false;
  DvceArray2D<int> ext_jk;
  // the extended pass reads the partner interpolants at the along-seam
  // ghost positions of an interior block end from q (CsSeamExtPartner) instead of
  // evaluating them on this block's own faces (one-sided at the outermost ghost)
  bool useq_ = false;
  DvceArray5D<Real> q;
  // <mesh>/cs_seam_fc_order = 3: srcval3 (block sizes for the angles)
  bool fc3_ = false;
  DualArray1D<RegionSize>::t_dev mbsz;
  bool fc3tab_ = false;            // evaluate from the coefficient table
  DvceArray1D<int> t_gbase, t_rq, t_cnt;
  DvceArray2D<int> t_code;
  DvceArray2D<Real> t_w;

  FCSeamPack(DvceFaceFld4D<Real> &b_, DvceFaceFld4D<Real> &cb_,
             DvceArray1D<MeshBufferDv> sb_, DvceArray1D<MeshBufferDv> rb_,
             DvceArray1D<PackSegFC> sg_, DvceArray1D<SubSegFC> su_,
             DvceArray1D<int> cs_, DvceArray1D<int> cb2_,
             DvceArray6D<Real> gt_, DvceArray6D<Real> st_, DvceArray1D<Real> pr_,
             int c0_, const RegionIndcs &ix_, bool nl_) :
             b(b_), cb(cb_), sbuf(sb_), rbuf(rb_), seg(sg_), sub(su_), cseg(cs_),
             cbeg(cb2_), cs_gtbl(gt_), cs_stbl(st_), pre(pr_), c0(c0_), cs_indcs(ix_),
             cs_nolim_(nl_) {}

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
      if (jf-1 < js_) {
        return 1.5*bx3(a,kf,js_,i) - 0.5*bx3(a,kf,js_+1,i);
      }
      if (jf   > je_) {
        return 1.5*bx3(a,kf,je_,i) - 0.5*bx3(a,kf,je_-1,i);
      }
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
      if (kf-1 < ks_) {
        return 1.5*bx2(a,ks_,jf,i) - 0.5*bx2(a,ks_+1,jf,i);
      }
      if (kf   > ke_) {
        return 1.5*bx2(a,ke_,jf,i) - 0.5*bx2(a,ke_-1,jf,i);
      }
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
  // the partner component at (kk,jj) -- x3f at an x2-face centre (x2 =
  // true) or x2f at an x3-face centre. In the extended pass, at an along-seam ghost
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
  // <mesh>/cs_seam_fc_order = 3: third-order average of the destination ghost face fed
  // by source sample (kk,jj,i) of component vv (cs_seam_fc3.hpp), in the source chart
  struct Fc3Geo {
    int m, i, aeta, low, dir, nseam, ns, ne, as0, alo, ahi, clo, chi, flo, fhi, nxa,
        npan, ga0, off;
    Real dn, da, a0, amin, amax;
  };
  KOKKOS_FORCEINLINE_FUNCTION
  Real fc3_aface(const Fc3Geo &g, const int f) const {
    // from the GLOBAL along face index (f - as0 + ga0), so every MeshBlock
    // layout computes the very same bits (block-local LeftEdgeX differs in the last bit,
    // and the order-3 lever arms near a panel end amplified that to ~1e-13 in the ghosts)
    return 0.5*M_PI*static_cast<Real>(f - g.as0 + g.ga0)/static_cast<Real>(g.npan)
           - 0.25*M_PI;
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
    // xabs: the GLOBAL along face coordinate; off = as0 - ga0 turns global into local
    // rounded with the tie-safe seam helpers (, cs_seam_round.hpp)
    int c0 = cs_seam::SeamFloor(xabs - 0.5) - 1 + g.off;
    c0 = (c0 < g.clo) ? g.clo : ((c0 > g.chi) ? g.chi : c0);
    c0 = (c0 < g.alo) ? g.alo : ((c0 > g.ahi - 2) ? g.ahi - 2 : c0);
    return c0;
  }
  KOKKOS_FORCEINLINE_FUNCTION
  int fc3_fwin(const Fc3Geo &g, const Real xabs) const {   // first of 4 along faces
    int w = cs_seam::SeamRound(xabs) - 2 + g.off;
    w = (w < g.flo) ? g.flo : ((w > g.fhi) ? g.fhi : w);
    w = (w < g.alo) ? g.alo : ((w > g.ahi - 2) ? g.ahi - 2 : w);
    return w;
  }
  // The point values are LINEAR in the source faces and every weight is geometry only,
  // so each reconstruction hands (component, along index, normal index, weight) terms to
  // an accumulator: Fc3Direct sums them on the fly, Fc3Tab records them (the per-sample
  // coefficient table of PackAndSendFCFlat, evaluated once per run). comp: 0 = x1f,
  // 1 = the across-staggered component, 2 = the along-staggered one.
  // across-staggered component (normal faces 0..3 from the seam) at (u, along face
  // coordinate xabs, along angle tangent Aa), times sc
  template <class Acc>
  KOKKOS_FORCEINLINE_FUNCTION
  void fc3_PtN(const Fc3Geo &g, const Real u, const Real xabs, const Real Aa,
               const Real sc, Acc &acc) const {
    const int c0 = fc3_cwin(g, xabs);
    Real wq[3], wl[4];
    cs_fc3::QuadDensW(xabs - 0.5 - static_cast<Real>(c0 - g.off + 1), wq);
    cs_fc3::Lag4W(u, wl);
    Real tf[4];
    for (int l=0; l<4; ++l) tf[l] = tan(fc3_aface(g, c0 + l));
    for (int mc=0; mc<4; ++mc) {
      const int col = g.nseam + g.dir*mc;
      const Real C = tan(0.25*M_PI - static_cast<Real>(mc)*g.dn);
      const Real f = sc*wl[mc]/(cs_fc3::ArcDens(C, Aa)*g.da);
      for (int l=0; l<3; ++l) acc(g, 1, c0 + l, col,
                                  f*wq[l]*cs_fc3::Arc(C, tf[l], tf[l+1]));
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
    cs_fc3::Lag4W(xabs - static_cast<Real>(w - g.off), wl);
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
    cs_fc3::QuadDensW(xabs - 0.5 - static_cast<Real>(c0 - g.off + 1), wa);
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
    // the panel resolution and this block's global along offset; every
    // angle below from global indices (layout-independent bits)
    g.npan = static_cast<int>(floor(2.0*g.nxa/(g.amax - g.amin) + 0.5));
    g.ga0 = static_cast<int>(floor(0.5*(g.amin + 1.0)*g.npan + 0.5));
    {const int npn = static_cast<int>(floor(2.0*nxn/(nmax - nmin) + 0.5));
     g.dn = 0.5*M_PI/static_cast<Real>(npn);}
    g.da = 0.5*M_PI/static_cast<Real>(g.npan);
    g.off = g.as0 - g.ga0;
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
        if (pass == 0) {
          den += wgt;
          continue;
        }
        const Real sc = wgt/acc.den;
        const Real Aa = Y*Nn;
        const Real asrc = atan(Aa);
        // global along face coordinate: the same bits in every layout
        const Real xabs = (asrc + 0.25*M_PI)/g.da;
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
      for (int q=0; q<n; ++q) {
        if (code[q] == c) {
          w[q] += wt;
          return;
        }
      }
      if (n < kFc3Max) {
        code[n] = c;
        w[n] = wt;
        ++n;
      } else {
        n = kFc3Max + 1;
      }
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
    // the srcval values evaluated once by the FCSeamPre pass
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
    // never srcval here (FCSeamPre did it): small, stack-free kernel
    const Real val = (u.pnt == 0) ? pre(u.pvoff + idx) : seamval_pre(a,kk,jj,i);
    if (s.local) {
      rbuf(s.dn).vars(s.dm, ix) = val;
    } else {
      sbuf(s.n).vars(s.m, ix) = val;
    }
  }

  // the partner interpolants of every ACTIVE face of block m, with the
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

  // the precompute pass of <mesh>/cs_seam_fc_order = 3 -- every seam
  // segment is sample-indexed (BuildPackSegmentsFC) and its cross-panel samples go
  // through srcval3. A kernel of its own, so the order-2 pass keeps its registers.
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
      } else {
        pre(u.pvoff + idx) = srcval3(a, s.cs_psrc, s.cs_pdst, s.ak*k + s.bk,
                                     s.aj*j + s.bj, i);
      }
    } else {
      pre(u.pvoff + idx) = srcval(a, s.ak*k + s.bk, s.aj*j + s.bj, i);
    }
  }

  // one coefficient-table row (seam segment, sample k,j); cnt -1: the
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
    if (!ok || acc.n > kFc3Max) {
      t_cnt(t) = -1;
      return;
    }
    t_cnt(t) = acc.n;
    for (int n=0; n<acc.n; ++n) {
      t_code(t,n) = acc.code[n];
      t_w(t,n) = acc.w[n];
    }
  }

  // srcval over a SubSegFC's (normal, tangent, i) box, once per pack
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
//! \!fn void MeshBoundaryValuesFC::PackAndSendFCFlat
//! \brief the cubed-sphere face-field pack (see BuildPackSegmentsFC).

void MeshBoundaryValuesFC::PackAndSendFCFlat(DvceFaceFld4D<Real> &b,
                                             DvceFaceFld4D<Real> &cb) {
  if (!fc_seg_ready) BuildPackSegmentsFC();
  const bool cs_nolim_ = pmy_pack->pmesh->cs_seam_fc_nolimit;
  auto sbuf = SendBufDv();
  auto rbuf = RecvBufDv();
  const bool fc3_ = (pmy_pack->pmesh->cs_seam_fc_order == 3);
  const int np = fc_np, ns = fc_ns, nr = fc_npre;
  if (np > 0) {
    FCPlainPack f(b, cb, sbuf, rbuf, fc_seg, fc_sub, fc_cs, fc_cb, 0);
    Kokkos::parallel_for("SendBuffFCplain", Kokkos::Experimental::require(
        Kokkos::RangePolicy<>(DevExeSpace(), 0, np*FC_PACK_CHUNK),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
  }
  // the second pass of <mesh>/cs_seam_blockend = extend
  auto stbl_ = cs_ext_pass_ ? cs_seam_stnc_x : cs_seam_stnc;
  if (nr > 0) {
    FCSeamPack f(b, cb, sbuf, rbuf, fc_seg, fc_sub, fc_cs, fc_cb, cs_seam_geom,
                 stbl_, fc_preval, np + ns, pmy_pack->pmesh->mb_indcs, cs_nolim_);
    f.ext_ = cs_ext_pass_; f.ext_jk = cs_ext_jk;
    f.useq_ = cs_ext_pass_ && cs_ext_q_; f.q = cs_ext_q;
    if (fc3_) {
      f.fc3_ = true; f.mbsz = pmy_pack->pmb->mb_size.d_view;
      // the coefficient table of this pass (built once)
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
    FCSeamPack f(b, cb, sbuf, rbuf, fc_seg, fc_sub, fc_cs, fc_cb, cs_seam_geom,
                 stbl_, fc_preval, np, pmy_pack->pmesh->mb_indcs, cs_nolim_);
    f.ext_ = cs_ext_pass_; f.ext_jk = cs_ext_jk;
    f.useq_ = cs_ext_pass_ && cs_ext_q_; f.q = cs_ext_q;
    Kokkos::parallel_for("SendBuffFCseam", Kokkos::Experimental::require(
        Kokkos::RangePolicy<>(DevExeSpace(), 0, ns*FC_PACK_CHUNK),
        Kokkos::Experimental::WorkItemProperty::HintLightWeight), f);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesFC::PackAndSendFC
//! \brief Packs face-centred fields into boundary buffers and sends them to neighbours
//! (all buffers of all MeshBlocks at once), including the cubed-sphere seam transform
//! and resample.

TaskStatus MeshBoundaryValuesFC::PackAndSendFC(DvceFaceFld4D<Real> &b,
                                               DvceFaceFld4D<Real> &cb) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;

  // CUBED-SPHERE SEAM GEOMETRY TABLE. Built once, lazily, on the first call: it needs
  // the neighbour table, which is not yet set when this object is constructed, and it
  // can never go stale because a cubed sphere refuses adaptive refinement and refuses
  // any static level boundary that would lie on a seam. See BuildCubedSphereSeamTables.
  if (pmy_pack->pmesh->use_cubed_sphere && !cs_seam_tbl_ready) {
    BuildCubedSphereSeamTables();
  }

  // on a cubed sphere the pack is the flat, load-balanced pair of kernels
  // (bitwise the same buffers); every other mesh keeps the team kernel below.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    PackAndSendFCFlat(b, cb);
  } else {
  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  const bool use_pole = pmy_pack->pmesh->use_polar_boundary;
  const bool cs_nolim_ = pmy_pack->pmesh->cs_seam_fc_nolimit;   // see mesh.cpp

  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto cs_gtbl = cs_seam_geom;
  auto cs_stbl = cs_seam_stnc;
  auto cs_sltb = cs_seam_slot;
  auto sbuf = SendBufDv();
  auto rbuf = RecvBufDv();

  // outer loop over (# of MeshBlocks)*(# of buffers)*(three field components)
  int nmnv = 3*nmb*nnghbr;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nmnv, Kokkos::AUTO);
  BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/(3*nnghbr);
    const int n = (tmember.league_rank()/3)%nnghbr;
    const int v = tmember.league_rank()%3;

    {
      // only load buffers when neighbor exists. CUBED-SPHERE CUBE VERTEX: skipped on
      // both sides; FillPanelCornersFC overwrites this corner block. See bvals.hpp.
      if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
        // if neighbor is at coarser level, use cindices to pack buffer
        // Note indices can be different for each component of face-centered field.
        int il, iu, jl, ju, kl, ku, ndat;

        // v is the receiver's component; the source component vv differs across a swap
        // seam, and its index range is used to read the source
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
            // field cb, not b. Before this existed the branch packed NOTHING, so the
            // buffer kept its zero-initialised contents and the coarse block read a
            // ghost field of exactly ZERO across the seam. See the note in bvals_cc.cpp.
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

          // x1 is radial: b.x1f crosses a seam unchanged; the tangential pair is
          // transformed to the destination face normals
              cs_xform = (v == 1) || (v == 2);
              cs_srcpanel = my_panel;
              cs_dstpanel = ngh.panel;
              cs_slotidx = cs_sltb.d_view(m, ngh.panel);
              cs_gkind = (cs_coar ? 2 : 0) + ((vv > 0) ? (vv - 1) : 0);
              cs_skind = (cs_coar ? 3 : 0) + vv;

              // only same-level face buffers are resampled along the seam; edge and
              // corner buffers are plain copies
              if (n >= 8 && n < 24) {
                cs_seam = 2;
              } else if (n >= 24 && n < 40) {
                cs_seam = 3;
              }
              // DOUBLY-GHOST BUFFERS (x2x3 edges 40-47, corners 48-55) are ghost in
              // BOTH tangential directions, so the slot alone does not name an along-seam
              // axis -- but the NEIGHBOUR TABLE does. Only one of the two flanking faces
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

              // the resample needs 3 cells along the seam
              const int seam_extent = cs_coar
                  ? ((cs_seam == 2) ? cs_indcs.cnx3 : cs_indcs.cnx2)
                  : ((cs_seam == 2) ? cs_indcs.nx3 : cs_indcs.nx2);
              if (cs_seam != 0 && seam_extent < 3) {
                cs_seam = 0;
              }
            } else if (do_pole) {
              aj = -1;
              bj = jl + ju;
              if (v == 1) signvar = -1;
              if (v == 2) signvar = -1;
            }

            // value of destination component v from the source at (kk,jj,i): both source
            // angular components are needed at the same point, so the partner component
            // is interpolated onto the primary's staggered location
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
            // the partner component is moved to the primary's location by fourth-order
            // interpolation on both axes (lag4)
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
                if (jf-1 < js_) {
                  return 1.5*bx3(kf,js_,i) - 0.5*bx3(kf,js_+1,i);
                }
                if (jf   > je_) {
                  return 1.5*bx3(kf,je_,i) - 0.5*bx3(kf,je_-1,i);
                }
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
                if (kf-1 < ks_) {
                  return 1.5*bx2(ks_,jf,i) - 0.5*bx2(ks_+1,jf,i);
                }
                if (kf   > ke_) {
                  return 1.5*bx2(ke_,jf,i) - 0.5*bx2(ke_-1,jf,i);
                }
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
            // (xi,eta) of a source sample with the staggering of component vv (b.x2f on
            // a xi face, b.x3f on an eta face): cs_gkind / cs_skind
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

            // along-seam resample (face-centred version of the one in bvals_cc.cpp): the
            // ghost face value is interpolated along the seam from 3 source faces
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
                // monotonicity limit: clamp to the stencil range where the data are
                // monotone; smooth extrema are left unlimited
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
          // seamval reads b0 or coarse_b0 according to cs_coar, so one expression
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
          if (ierr != MPI_SUCCESS) {
            no_errors=false;
          }
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
//! \fn void MeshBoundaryValuesFC::FillPanelCornersFC
//! \brief Fills the ng x ng corner ghost faces at a cube vertex, where only three panels
//! meet and the generic corner buffer has no meaningful source: each component is
//! extrapolated from the two flanking face halos in the panel's own chart and averaged
//! (unless the cube-vertex corner fill already wrote it).

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
  const bool vfill_ = pmy_pack->pmesh->use_cubed_sphere;
  const bool coar_ = coarse;
  const int is_ = indcs.is, ie_ = indcs.ie;
  auto &mblev_ = pmy_pack->pmb->mb_lev;

  par_for("cs_fill_corners", DevExeSpace(), 0,(nmb-1), 0,3, 0,nc1, 0,(ng-1), 0,(ng-1),
  KOKKOS_LAMBDA(const int m, const int c, const int i, const int gj, const int gk) {
    const int sj = (c & 1) ? 1 : -1;      // -1 = the -x2 side of the block
    const int sk = (c & 2) ? 1 : -1;      // -1 = the -x3 side
    // Buffer ids of the two FACE neighbours flanking this corner (see nghbr_index.hpp).
    const int nj_id = (sj < 0) ? 8 : 12;
    const int nk_id = (sk < 0) ? 24 : 28;
    const int mp = mbpanel.d_view(m);
    // ONLY A TRUE CUBE VERTEX. The condition is that BOTH flanking faces are panel
    // seams -- exactly IsCubeVertexCorner's test, and exactly the case whose exchange is
    // skipped. It used to fire when EITHER face was a seam, which with more than one
    // MeshBlock per panel also caught ordinary corners that DO have a real diagonal
    // neighbour, overwriting properly exchanged, seam-transformed data with a one-sided
    // extrapolation.
    bool seamj = false, seamk = false;
    if (nghbr.d_view(m,nj_id).gid >= 0 && nghbr.d_view(m,nj_id).panel != mp) seamj = true;
    if (nghbr.d_view(m,nk_id).gid >= 0 && nghbr.d_view(m,nk_id).panel != mp) seamk = true;
    if (!(seamj && seamk)) return;
    // the sampled corner fill (same-level flanking neighbours) already wrote this block;
    // extrapolate only otherwise
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
    // and CUBIC (15 and 49) loses to amplification -- 1.98e-03 -> 2.46e-03. Do not
    // re-propose either.
    const Real wj0 = 0.5*(dj+1.0)*(dj+2.0), wj1 = -dj*(dj+2.0), wj2 = 0.5*dj*(dj+1.0);
    const Real wk0 = 0.5*(dk+1.0)*(dk+2.0), wk1 = -dk*(dk+2.0), wk2 = 0.5*dk*(dk+1.0);
    const int stj = -sj, stk = -sk;       // step INWARD from the anchor
    // The two extrapolations do not reach equally far: the one along j travels dj cells
    // and the one along k travels dk. A plain average weights them the same even when
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
  });
  return;
}


//----------------------------------------------------------------------------------------
// \!fn void RecvBuffers
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
  const bool cs_vfill_ = pmy_pack->pmesh->use_cubed_sphere;
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
          if (ierr != MPI_SUCCESS) {
            no_errors=false;
          }
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
  if (bflag) {
    return TaskStatus::incomplete;
  }
#endif

  //----- STEP 2: buffers have all completed, so unpack 3-components of field

  // Outer loop over (# of MeshBlocks)*(three field components)*(radial chunks)
  // the radial index range is cut into chunks of FC_UNP_CHUNK, one team
  // each. Buffers overlap only at the SAME (k,j,i), which stays in one team, so each
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
        // this team's radial slice [cl, cu] of [il, iu]
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
        // cube-vertex corner fill: the face buffers were widened by ng along the seam, so
        // both flanking buffers reach the corner block; each writes its own half
            const bool face_wide = (n >= 8 && n < 40);
            // STAGGERING. b.x2f lives on x2 FACES and b.x3f on x3 faces, so for those
            // components je+1 / ke+1 is a REAL interior face, not a ghost. Testing the
            // cell bounds alone misclassifies it and suppresses a legitimate write --
            // which is why a widening of ZERO cells still destroyed the solution.
            const int jghi = je_ + ((v == 1) ? 1 : 0);
            const int kghi = ke_ + ((v == 2) ? 1 : 0);
            if (cs_vfill_ && face_wide &&
                (j < js_ || j > jghi) && (k < ks_ || k > kghi)) {
              // ONLY a CUBE VERTEX needs the fill. At an ordinary corner the diagonal
              // neighbour exists and its own buffer (40..55, unpacked LATER in this same
              // scalar n loop) is the correct filler; the panel-edge geometry used below
              // is not even valid there, because the block edge is not the panel edge.
              const int sj = (j < js_) ? 0 : 1;
              const int sk = (k < ks_) ? 0 : 1;
              if (!IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m,
                                      40 + 2*(sj + 2*sk))) return;
            // ownership: in this panel's gnomonic chart the seam between the two flanking
            // panels is the diagonal of the corner block; each panel writes its own side
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

  // all face buffers are unpacked: fill the cube-vertex corners of the fine and (with
  // multilevel) the coarse array
  if (pmy_pack->pmesh->use_cubed_sphere) {
    FillPanelCornersFC(b, false);
    if (pmy_pack->pmesh->multilevel) {
      FillPanelCornersFC(cb, true);
    }
  }

  // <mesh>/cs_seam_blockend = extend, the second pass (see the note in
  // RecvAndUnpackCC): repeat the exchange with the panel-clamped stencil table.
  if (cs_blockend_ext_ && !cs_ext_pass_ && pmy_pack->pmesh->use_cubed_sphere) {
    if (cs_ext_q_) CsSeamExtPartner(b);
    cs_ext_pass_ = true;
    (void) ClearSend();
    (void) InitRecv(3);
    while (PackAndSendFC(b, cb) != TaskStatus::complete) {}
    while (RecvAndUnpackFC(b, cb) != TaskStatus::complete) {}
    // wait on this pass's own sends here. A caller that cleared its
    // sends before unpacking (rad_m1_coupling.cpp FluidRefresh) would otherwise leave
    // them in flight into its next exchange, which reuses the buffers and the requests.
    (void) ClearSend();
    cs_ext_pass_ = false;
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesFC::CsSeamExtPartner
//! \brief With <mesh>/cs_seam_blockend = extend and more than 2 MeshBlocks per panel
//! side, recomputes the partner-component interpolation of the extended seam pass where
//! its 4-point stencils reach beyond the sender's valid along-seam ghosts.

void MeshBoundaryValuesFC::CsSeamExtPartner(DvceFaceFld4D<Real> &b) {
  auto &id = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ncell = nmb*(id.nx3 + 2*id.ng)*(id.nx2 + 2*id.ng)*(id.nx1 + 2*id.ng);
  FCSeamPack f(b, b, SendBufDv(), RecvBufDv(), fc_seg, fc_sub, fc_cs, fc_cb,
               cs_seam_geom, cs_seam_stnc_x, fc_preval, 0, id, false);
  f.ext_ = true; f.ext_jk = cs_ext_jk; f.q = cs_ext_q;
  Kokkos::parallel_for("CsSeamExtPartner",
      Kokkos::RangePolicy<FCSeamQ>(DevExeSpace(), 0, ncell), f);
  DvceArray5D<Real> q = cs_ext_q;
  while (pbq_->InitRecv(4) == TaskStatus::incomplete) {}
  while (pbq_->PackAndSendCC(q, q) == TaskStatus::incomplete) {}
  while (pbq_->RecvAndUnpackCC(q, q) == TaskStatus::incomplete) {}
  while (pbq_->ClearSend() == TaskStatus::incomplete) {}
}
