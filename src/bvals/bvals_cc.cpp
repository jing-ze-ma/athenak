//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for cell-centered (CC)
//! Mesh variables.
//! Prolongation of CC variables  occurs in ProlongateCC() function called from task list

#include <cstdlib>
#include <iostream>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
// BValCC constructor:

MeshBoundaryValuesCC::MeshBoundaryValuesCC(MeshBlockPack *pp, ParameterInput *pin,
                                           bool z4c) :
  MeshBoundaryValues(pp, pin, z4c) {
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::PackAndSendCC()
//! \brief Pack cell-centered variables into boundary buffers and send to neighbors.
//!
//! This routine packs ALL the buffers on ALL the faces, edges, and corners simultaneously
//! for ALL the MeshBlocks. This reduces the number of kernel launches when there are a
//! large number of MeshBlocks per MPI rank. Buffer data are then sent (via MPI) or copied
//! directly for periodic or block boundaries.
//!
//! Input arrays must be 5D Kokkos View dimensioned (nmb, nvar, nx3, nx2, nx1)
//! 5D Kokkos View of coarsened (restricted) array data also required with SMR/AMR

TaskStatus MeshBoundaryValuesCC::PackAndSendCC(DvceArray5D<Real> &a,
                                               DvceArray5D<Real> &ca,
                                               int iwl, int iwu) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  int nvar = a.extent_int(1);  // TODO(@user): 2nd index from L of in array must be NVAR

  {int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  // skip the DIAGONAL transverse buffers entirely (see bvals.hpp)
  const bool skipd_ = skip_x2x3_diag;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  const bool use_pole = pmy_pack->pmesh->use_polar_boundary;
  // needed only on the cubed sphere, to give a source cell its (xi,eta)
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  int no_rs_cc = 0;
  { const char *e_ = std::getenv("CS_NORESAMP_CC");
    if (e_ != nullptr) { no_rs_cc = std::atoi(e_); } }
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;
  auto &is_z4c = is_z4c_;
  auto &multilevel = pmy_pack->pmesh->multilevel;
  // x1 index window (see the declaration): only buffers spanning the whole x1 extent
  const int wlo_ = iwl, whi_ = iwu;
  const int mbis_ = pmy_pack->pmesh->mb_indcs.is;
  const int mbie_ = pmy_pack->pmesh->mb_indcs.ie;
  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables)
  int nmnv = nmb*nnghbr*nvar;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nmnv, Kokkos::AUTO);
  Kokkos::parallel_for("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    // only load buffers when neighbor exists.  CUBED-SPHERE CUBE VERTEX: skipped on
    // both sides; FillPanelCornersCC overwrites exactly this corner block. See bvals.hpp.
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n)) &&
        !(skipd_ && IsX2X3DiagSlot(n))) {
      // if neighbor is at coarser level, use coar indices to pack buffer
      int il, iu, jl, ju, kl, ku;
      if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
        il = sbuf[n].icoar[0].bis;
        iu = sbuf[n].icoar[0].bie;
        jl = sbuf[n].icoar[0].bjs;
        ju = sbuf[n].icoar[0].bje;
        kl = sbuf[n].icoar[0].bks;
        ku = sbuf[n].icoar[0].bke;
      // if neighbor is at same level, use same indices to pack buffer
      } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
        il = sbuf[n].isame[0].bis;
        iu = sbuf[n].isame[0].bie;
        jl = sbuf[n].isame[0].bjs;
        ju = sbuf[n].isame[0].bje;
        kl = sbuf[n].isame[0].bks;
        ku = sbuf[n].isame[0].bke;
      // if neighbor is at finer level, use fine indices to pack buffer
      } else {
        il = sbuf[n].ifine[0].bis;
        iu = sbuf[n].ifine[0].bie;
        jl = sbuf[n].ifine[0].bjs;
        ju = sbuf[n].ifine[0].bje;
        kl = sbuf[n].ifine[0].bks;
        ku = sbuf[n].ifine[0].bke;
      }
      int ni = iu - il + 1;
      int nj = ju - jl + 1;
      int nk = ku - kl + 1;
      int nkj  = nk*nj;
      // the window narrows the LOOP only: ni and the index map stay as they were
      int ilp = il, iup = iu;
      // same-level x2/x3 halos only: coarse/fine ranges differ between the two sides
      // and the cubed-sphere corner fill reads the full strips
      if (wlo_ >= 0 && !use_cs && nghbr.d_view(m,n).lev == mblev.d_view(m) &&
          il <= mbis_ && iu >= mbie_) {
        ilp = (il > wlo_) ? il : wlo_;
        iup = (iu < whi_) ? iu : whi_;
      }

      // indices of recv'ing (destination) MB and buffer: MB IDs are stored sequentially
      // in MeshBlockPacks, so array index equals (target_id - first_id)
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
        // array ca, not a.  Before this existed the whole branch packed NOTHING, so the
        // buffer kept its zero-initialised contents and the coarse block read a ghost
        // state of exactly ZERO across the seam.  Invisible to any gate with v = 0
        // (zero IS the right answer there) and invisible to hydro (a dimensionally
        // split PLM sweep never reconstructs through a corner ghost), but the corner
        // EMF of the CT update reads DIAGONALLY, which is why it only showed in MHD.
        const bool cs_coar = (nghbr.d_view(m,n).lev < mblev.d_view(m));
        int cs_srcpanel = 0, cs_dstpanel = 0;
        // 0 = no along-seam resample; 2 = x2-face seam (resample in k);
        // 3 = x3-face seam (resample in j). See the note on seamval below.
        int cs_seam = 0;

        if (do_cs) {
          const auto ngh = nghbr.d_view(m,n);
          const int my_panel = mbpanel.d_view(m);
          PanelBoundaries pb;
          pb = GetPanelBoundary(my_panel, ngh.panel);

          // x1 is RADIAL on the cubed sphere and no seam crosses it, so i and IVX pass
          // through untouched. The panel-tangential pair is a = x2 (j, IVY) and
          // b = x3 (k, IVZ). See the axis note in mesh.hpp.
          // The INDEX map across a seam is a signed permutation and is handled below.
          // The VECTOR COMPONENTS are not: the two charts carry different tangent bases
          // at the same physical point, differing by a shear that is O(1) away from the
          // seam midline. They are transformed properly, per source cell, by
          // cubed_sphere::TransformMomentum -- see the note there. IVX is radial and
          // common to both charts, so it passes through.
          cs_xform = (v == IVY) || (v == IVZ);
          cs_dstpanel = ngh.panel;
          cs_srcpanel = my_panel;

          // pb.rev_a/rev_b are expressed in the DESTINATION panel's axes, so undo the
          // swap to get the reversal that applies to this block's own j and k.
          int rev_a_preswap = (pb.swap_ax == 1) ? pb.rev_b : pb.rev_a;
          int rev_b_preswap = (pb.swap_ax == 1) ? pb.rev_a : pb.rev_b;
          if (pb.swap_ax == 1) {
            // Transpose the two tangential axes IN THE BUFFER, so the receiver's generic
            // unpack lands them on its own j,k. The receiver's (nj,nk) are this block's
            // (nk,nj), which is why the i stride ni and the per-variable stride are the
            // same either way.
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

          // WHICH BUFFER IS THIS, and therefore which way does the along-seam resample
          // run?  Read it off the SLOT INDEX.  Slots are laid out (nghbr_index.hpp):
          // 0-7 x1 faces, 8-15 x2 faces, 16-23 x1x2 edges, 24-31 x3 faces, 32-39 x3x1
          // edges, 40-47 x2x3 edges, 48-55 corners.  A buffer whose only tangential
          // ghost direction is x2 has its seam normal along x2 and resamples in x3
          // (cs_seam = 2), and vice versa; one that is ghost in BOTH tangential
          // directions -- the x2x3 edges and the corners -- has no single along-seam
          // axis and stays a plain copy.
          //
          // This used to be inferred from the buffer EXTENTS, as `nj == ng && nk > ng`.
          // For cell-centred data that happens to give the same answer, but the identical
          // test in bvals_fc.cpp did NOT: a face-centred buffer at a coarse/fine boundary
          // carries one extra layer in the direction its component is staggered in, so
          // the equality failed and the resample was silently skipped there.  Both are
          // written the same way now so the two cannot drift apart again.
          if (n >= 8 && n < 24) {
            cs_seam = 2;
          } else if (n >= 24 && n < 40) {
            cs_seam = 3;
          } else if (n >= 40 && n < 56) {
            // x2x3 EDGE (40-47) and x1x2x3 CORNER (48-55) buffers.  Both are ghost in
            // BOTH tangential directions, so the slot alone does not name an along-seam
            // axis -- but the NEIGHBOUR TABLE does, and the treatment is identical.
            //
            // The CORNER half was added after the edges: x1 is RADIAL, so a corner ghost
            // is a radial ghost cell that is also tangentially ghost in x2 and x3, read
            // by multi-D stencils near a radial block boundary at a seam and by
            // prolongation.  bvals_fc.cpp already covered slots 40-55 in one test; this
            // brings bvals_cc.cpp to the same coverage.  Measured (static corner scan,
            // iprob 15, max |ghost - exact| over corner ghosts with exactly one seam
            // flank, 2 radial blocks): n = 32 per panel, 2x2 blocks 1.407e-4 -> 7.99e-8,
            // 4x4 blocks 7.414e-4 -> 5.08e-6; n = 64, 2x2 3.498e-5 -> 4.28e-9, 4x4
            // 3.538e-4 -> 5.64e-7 -- i.e. onto the x2x3 EDGE-seam bin to every digit.
            // At 1 x 1 block per panel every cross-panel corner is a CUBE VERTEX (both
            // flanks are seams), so cs_seam stays 0 there and 1 x 1 is bitwise unchanged.
            //
            // x2x3 EDGE buffers.  These used to be left as a plain copy, on the grounds
            // that a doubly-ghost buffer "has no single along-seam axis".  It does:
            // EXACTLY ONE of the two flanking faces is a panel seam, because if both were
            // this would be a CUBE VERTEX and the exchange is skipped altogether
            // (IsCubeVertexCorner).  The seam normal is that face's axis and the resample
            // runs along the other one, exactly as for the face buffer next to it.
            //
            // Leaving it a plain copy is the SECOND half of the 4x4-blocks-per-panel
            // defect (tests_seam4/README.md): with one block per panel every cross-panel
            // x2x3 edge IS a cube vertex, so nothing was ever wrong there, but as soon as
            // a panel is split these are ordinary diagonal ghosts with a real donor --
            // and they are what the metric cross term of the transverse operator reads,
            // which is why the EXPLICIT operator did not improve when only the face
            // halo was fixed.  Measured (static halo scan, n = 64 per panel, max
            // |ghost - exact| over the x2x3 edge ghosts that have a seam flank):
            // 2 x 2 blocks 3.50e-5 -> 5.6e-7, 4 x 4 blocks 3.54e-4 -> 5.6e-7.
            // Slot -> flanking tangential faces.  x2x3 edge: n = 40 + n1 + 2*(iy>0) +
            // 4*(iz>0).  Corner: n = 48 + (ix>0) + 2*(iy>0) + 4*(iz>0) -- bit 0 is the
            // RADIAL side and carries no tangential information.
            int fj, fk;
            if (n < 48) {
              const int q = n - 40;
              fj = (q/2) & 1;
              fk = (q/4) & 1;
            } else {
              const int q = n - 48;
              fj = (q >> 1) & 1;
              fk = (q >> 2) & 1;
            }
            const int nface2 = (fj == 0) ? 8 : 12;         // the flanking x2 face
            const int nface3 = (fk == 0) ? 24 : 28;        // the flanking x3 face
            const bool s2 = (nghbr.d_view(m,nface2).gid >= 0) &&
                            (nghbr.d_view(m,nface2).panel != my_panel);
            const bool s3 = (nghbr.d_view(m,nface3).gid >= 0) &&
                            (nghbr.d_view(m,nface3).panel != my_panel);
            if (s2 && !s3) {
              cs_seam = 2;
            } else if (s3 && !s2) {
              cs_seam = 3;
            }
          }
          // THE RESAMPLE NEEDS 3 CELLS ALONG THE SEAM -- see the note in bvals_fc.cpp,
          // where this degeneracy was measured.  The stencil bounds invert when the
          // source's along-seam ACTIVE extent is under 3, which happens on the COARSE
          // array as soon as cnx = nx/2 < 3, and MeshBlocks are allowed down to 4 cells.
          // No cell-centred gate resolves a difference here, but the code shape is the
          // same and reading unfilled ghost is not something to leave to luck.
          const int seam_extent = cs_coar
              ? ((cs_seam == 2) ? cs_indcs.cnx3 : cs_indcs.cnx2)
              : ((cs_seam == 2) ? cs_indcs.nx3 : cs_indcs.nx2);
          if (cs_seam != 0 && seam_extent < 3) { cs_seam = 0; }
          if (no_rs_cc) { cs_seam = 0; }
        } else if (do_pole) {
          aj = -1;
          bj = jl + ju;
          if (v == IVY) signvar = -1;
          if (v == IVZ) signvar = -1;
        }

        // THE SEAM GEOMETRY IS HOISTED OUT OF THE RADIAL LOOP.
        //
        // x1 is RADIAL on the cubed sphere and no seam crosses it, so (xi,eta) -- and
        // therefore EVERYTHING the transform and the resample compute from it -- is the
        // same for every cell in the i loop. It used to be recomputed per cell: four
        // PanelFrame switches, four tan/sqrt pairs, a PanelToCart and a CartToPanel with
        // two atan, plus the resample's own atan/tan, once for every radial cell of every
        // seam buffer. Profiled on a 300-cycle dhj MHD run, this pack cost 14.9x per call
        // on the cubed sphere against spherical polar, a quarter of the whole GPU time.
        //
        // Everything below is per (kk,jj) and none of it depends on the data, so the
        // hoist is BITWISE EXACT: each per-cell floating-point expression is unchanged
        // and still evaluated in the same order.  MEASURED on that run: this kernel went
        // 2267 ms -> 254 ms, 8.9x, for no change in any output byte.
        //
        // THE SAME HOIST IN bvals_fc.cpp MADE IT 6x SLOWER and is deliberately NOT there.
        // The face-centred pack needs THREE transforms per column (one per resample
        // stencil cell) rather than one, and its loop is flat over (k,j,i); turning it
        // into TeamThreadRange(k,j) x ThreadVectorRange(i) to get a place to hoist into
        // both serialises the radial direction (the TeamPolicy's vector length is 1) and
        // pushes the per-thread private segment from 908 to 2224 bytes.  Measured
        // 1432 ms -> 8596 ms.  Any retry needs an explicit vector length on the policy
        // and a way to carry the three transforms without a per-thread array.
        const int js_ = cs_coar ? cs_indcs.cjs : cs_indcs.js;
        const int ks_ = cs_coar ? cs_indcs.cks : cs_indcs.ks;
        const int nx2_ = cs_coar ? cs_indcs.cnx2 : cs_indcs.nx2;
        const int nx3_ = cs_coar ? cs_indcs.cnx3 : cs_indcs.nx3;
        const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
        const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;

        // Middle loop over k,j
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;
          int kk = ak*k + bk;
          int jj = aj*j + bj;

          // ALONG-SEAM RESAMPLE. Across a panel seam the two charts share the seam-normal
          // coordinate exactly but NOT the seam-parallel one. Writing the seam-normal
          // angle of a source cell as n and its seam-parallel angle as a, the physical
          // point of that cell sits at seam-parallel angle atan(tan(a)/tan|n|) in the
          // DESTINATION chart, not at a. So the plain index copy hands each ghost cell
          // the state of a point up to half a cell (layer 0) or ~1.5 cells (layer 1) away
          // along the seam -- an offset that does NOT shrink with resolution in cell
          // units, which makes the ghost value O(dx) wrong and the acceleration it drives
          // O(1)... i.e. the seam is only first-order accurate.
          //
          // Inverting that map, the value a ghost needs is the source field at
          // seam-parallel angle atan(tan(a)*tan|n|), which is always INSIDE the source
          // cell's own angle (|tan n| < 1), so the stencil never leaves the source block.
          // The same formula covers both orientations: a reversed seam flips the sign of
          // both a and the target, and the index reversal above already carries that.
          // Quadratic (3-point) Lagrange keeps the ghost error O(dx^3), which is what the
          // second-order flux difference needs.
          //
          // The stencil cells and their weights are a function of (kk,jj) alone. With no
          // resample there is a single "stencil cell", the source cell itself.
          int kst[3] = {kk, kk, kk};
          int jst[3] = {jj, jj, jj};
          int nst = 1;
          Real wm = 1.0, w0 = 0.0, wp = 0.0;
          // Is the stencil INTERPOLATING?  `b` below is clamped to the source's active
          // range, so `pos` can fall outside [b, b+2] and the quadratic then
          // EXTRAPOLATES -- see the note on the monotonicity limit in seamval().
          bool cs_interp = true;
          if (cs_seam != 0) {
            const Real xi  = 0.25*M_PI*CellCenterX(jj-js_, nx2_, x2mn, x2mx);
            const Real eta = 0.25*M_PI*CellCenterX(kk-ks_, nx3_, x3mn, x3mx);
            Real ang, nrm, dang;
            int sc, blo, bhi;
            if (cs_seam == 2) {
              ang = eta; nrm = xi;
              dang = 0.25*M_PI*(x3mx - x3mn)/static_cast<Real>(nx3_);
              sc = kk;
              blo = ks_; bhi = ks_ + nx3_ - 3;
            } else {
              ang = xi; nrm = eta;
              dang = 0.25*M_PI*(x2mx - x2mn)/static_cast<Real>(nx2_);
              sc = jj;
              blo = js_; bhi = js_ + nx2_ - 3;
            }
            // THE CLAMP IS TO THE SOURCE'S ACTIVE RANGE, NOT THE BUFFER'S.  For a face or
            // x1-edge buffer the two coincide (the along-seam extent of such a buffer IS
            // the active range), so this is bitwise for everything that resampled before;
            // an x2x3 edge buffer is only ng deep along the seam and clamping to that
            // would extrapolate from two cells while the source block holds the data.
            // The stencil reads the array directly, not the buffer, so any active cell is
            // available.  bvals_fc.cpp says the same thing about its own stencil table.
            const Real pos = sc + (atan(tan(ang)*tan(fabs(nrm))) - ang)/dang;
            int b = static_cast<int>(floor(pos + 0.5)) - 1;
            b = (b < blo) ? blo : ((b > bhi) ? bhi : b);
            const Real u = pos - static_cast<Real>(b + 1);
            cs_interp = (u >= -1.0) && (u <= 1.0);
            wm = 0.5*u*(u - 1.0);
            w0 = 1.0 - u*u;
            wp = 0.5*u*(u + 1.0);
            nst = 3;
            for (int s=0; s<3; ++s) {
              if (cs_seam == 2) {
                kst[s] = b + s;
              } else {
                jst[s] = b + s;
              }
            }
          }

          // The tangent-basis geometry of each stencil cell. Across a panel seam the two
          // charts carry different tangent bases at the same physical point, differing by
          // a shear that is O(1) away from the seam midline, so the two tangential
          // momenta must be re-expressed rather than permuted -- see cubed_sphere::
          // TransformMomentum. IVX is radial, common to both charts, and passes through,
          // as does every scalar; those take the plain copy path and need no geometry.
          cubed_sphere::SeamXform xf[3];
          if (cs_xform) {
            for (int s=0; s<nst; ++s) {
              const Real xis = 0.25*M_PI*CellCenterX(jst[s]-js_, nx2_, x2mn, x2mx);
              const Real etas = 0.25*M_PI*CellCenterX(kst[s]-ks_, nx3_, x3mn, x3mx);
              cubed_sphere::SeamXformAt(cs_srcpanel, cs_dstpanel, xis, etas, xf[s]);
            }
          }

          // Value of variable v at stencil cell s, radial index i. Both tangential
          // momenta of the source cell are needed to make either one, so this reads two
          // components and selects. Away from a panel seam it is exactly the old
          // `a(...)*signvar`. Only loads and the O(1) half of the transform remain here.
          auto sval = [&](const int s, const int i) {
            const int kq = kst[s], jq = jst[s];
            if (!cs_xform) {
              return (cs_coar ? ca(m,vv,kq,jq,i) : a(m,vv,kq,jq,i))*signvar;
            }
            Real m2o, m3o;
            const Real my_ = cs_coar ? ca(m,IVY,kq,jq,i) : a(m,IVY,kq,jq,i);
            const Real mz_ = cs_coar ? ca(m,IVZ,kq,jq,i) : a(m,IVZ,kq,jq,i);
            cubed_sphere::ApplyMomentumXform(xf[s], my_, mz_, m2o, m3o);
            return (v == IVY) ? m2o : m3o;
          };

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
          // It must also be skipped when the resample is EXTRAPOLATING -- `b` is
          // clamped at the ends of the source range, so `pos` can fall outside
          // [b, b+2], where the correct value legitimately lies outside the node
          // range; clamping there cost 3.5x on the smooth seam halo.
          //
          // TWO EARLIER ATTEMPTS FAILED and are recorded so they are not retried: an
          // UNCONDITIONAL clamp (cost 3.5x, the extrapolation case), and a
          // second-difference roughness test gated on a RELATIVE span
          // (hi-lo) > 0.1*(|hi|+|lo|), which fires spuriously wherever the stencil
          // straddles ZERO -- it triggered ten million times on a smooth run.
          //
          // CELL-CENTRED: clamp for every INTERPOLATING stencil.  The face-centred twin
          // guards this with "monotone stencil AND interpolating"; here the MONOTONE half
          // must be dropped, because a blast cap has a FLAT TOP, so a stencil astride its
          // edge is NOT monotone, and a guard meant for smooth extrema was then skipping
          // exactly the cells that needed clamping.  Non-monotone does not imply smooth.
          // Measured on iprob=12: with the monotone guard, a VERTEX-centred blast still
          // went entirely NaN at contrast 100; without it, it survives 1000.
          //
          // The INTERPOLATING half must stay.  It used to be dropped as well, and that is
          // HALF of the 4x4-MeshBlocks-per-panel seam defect (the other half is the x2x3
          // edge buffers above; see tests_seam4/README.md): `b` is clamped to the
          // source's ACTIVE range, which is one MeshBlock, not the panel.  The
          // along-seam resample pulls the donor position toward the PANEL CENTRE by up
          // to (layer + 1/2) cells, so at an
          // interior block end -- an end that exists only when a panel is split more
          // finely than 2 x 2, since at 2 x 2 the interior end sits exactly on the panel
          // midline where the pull is zero -- the stencil runs off the block and the
          // quadratic extrapolates.  Clamping an
          // EXTRAPOLANT to its node range replaces it by the nearest node, i.e. by a
          // FIRST-ORDER value, and that is what wrecked the halo.  Measured on the static
          // halo scan (cs_test iprob 15, <problem>/seam_halo_scan, n = 32 per panel,
          // max |ghost - exact|):
          //
          //   blocks/panel   clamped unconditionally   clamped only when interpolating
          //   1 x 1          4.103e-6                  4.103e-6   (unchanged, bitwise)
          //   2 x 2          4.671e-5                  4.103e-6   (= 1 x 1 exactly)
          //   4 x 4          6.950e-4                  5.077e-6
          //
          // The residual 1.2x at 4 x 4 is the extrapolation itself, which is O(dx^3) and
          // harmless; the 137x was the clamp.
          auto seamval = [&](const int i) {
            if (cs_seam == 0) return sval(0,i);
            const Real sv0 = sval(0,i);
            const Real sv1 = sval(1,i);
            const Real sv2 = sval(2,i);
            const Real qq = wm*sv0 + w0*sv1 + wp*sv2;
            if (!cs_interp) return qq;
            const Real lo0 = fmin(sv0, fmin(sv1, sv2));
            const Real hi0 = fmax(sv0, fmax(sv1, sv2));
            return fmin(hi0, fmax(lo0, qq));
          };

          // Inner (vector) loop over i
          // copy directly into recv buffer if MeshBlocks on same rank

          // seamval() reads u0 or coarse_u0 according to cs_coar, so one expression
          // serves a neighbour at ANY level.
          if (nghbr.d_view(m,n).rank == my_rank) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(i);
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*v);
              rbuf[dn].vars(dm, index) = val;
            });

          // else copy into send buffer for MPI communication below

          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(i);
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*v);
              sbuf[n].vars(m,index) = val;
            });
          }
        });

      } else {   // normal boundary exchange
      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
        int k = idx / nj;
        int j = (idx - k * nj) + jl;
        k += kl;

        // Inner (vector) loop over i
        // copy directly into recv buffer if MeshBlocks on same rank

        if (nghbr.d_view(m,n).rank == my_rank) {
          // if neighbor is at same or finer level, load data from u0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = ca(m,v,k,j,i);
            });
          }

        // else copy into send buffer for MPI communication below

        } else {
          // if neighbor is at same or finer level, load data from u0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = ca(m,v,k,j,i);
            });
          }
        }
      });
      }  // end if-do-cs/do-pole block
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  }); // end par_for_outer

  // The whole body of this kernel sits under (is_z4c && multilevel), both of which are
  // host-side constants, so when they do not hold the launch has nothing to do at all:
  // skip it.  It is otherwise a full-size empty team launch on EVERY halo exchange.
  if (is_z4c_ && ml_) {
  Kokkos::parallel_for("SendBuffZ4c", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    // only load buffers when neighbor exists.  CUBED-SPHERE CUBE VERTEX: skipped on
    // both sides; FillPanelCornersCC overwrites exactly this corner block. See bvals.hpp.
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n)) &&
        !(skipd_ && IsX2X3DiagSlot(n))) {
      int il, iu, jl, ju, kl, ku;
      // If neighbor is at same level and data is for Z4c module, append data from coarse
      // array for higher-order prolongation
      if ((nghbr.d_view(m,n).lev == mblev.d_view(m)) && (is_z4c) && (multilevel)) {
        il = sbuf[n].isame_z4c.bis;
        iu = sbuf[n].isame_z4c.bie;
        jl = sbuf[n].isame_z4c.bjs;
        ju = sbuf[n].isame_z4c.bje;
        kl = sbuf[n].isame_z4c.bks;
        ku = sbuf[n].isame_z4c.bke;
        int ni = iu - il + 1;
        int nj = ju - jl + 1;
        int nk = ku - kl + 1;
        int nkj  = nk*nj;
        int ndat = nvar*sbuf[n].isame_ndat; // size of same level data already in buff

        // indices of recv'ing (destination) MB and buffer: MB IDs are stored sequentially
        // in MeshBlockPacks, so array index equals (target_id - first_id)
        int dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
        int dn = nghbr.d_view(m,n).dest;

        // Middle loop over k,j
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;

          // Inner (vector) loop over i
          // copy directly into recv buffer if MeshBlocks on same rank
          if (nghbr.d_view(m,n).rank == my_rank) {
            // load data from coarse_u0
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              rbuf[dn].vars(dm,ndat+ (i-il + ni*(j-jl + nj*(k-kl + nk*v))))=ca(m,v,k,j,i);
            });

          // else copy into send buffer for MPI communication below
          } else {
            // load data from coarse_u0
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              sbuf[n].vars(m,ndat+ (i-il + ni*(j-jl + nj*(k-kl + nk*v))) )=ca(m,v,k,j,i);
            });
          }
        });
      }
    } // end if-neighbor-exists block
    tmember.team_barrier();
  }); // end par_for_outer
  }
  }

#if MPI_PARALLEL_ENABLED
  // Send boundary buffer to neighboring MeshBlocks using MPI
  Kokkos::fence();
  auto &is_z4c = is_z4c_;
  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  // skip the DIAGONAL transverse buffers entirely (see bvals.hpp)
  const bool skipd_ = skip_x2x3_diag;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0 &&
          !(use_cs && IsCubeVertexCorner(nghbr.h_view, mbpanel.h_view, m, n)) &&
          !(skipd_ && IsX2X3DiagSlot(n))) {
        // index and rank of destination Neighbor
        int dn = nghbr.h_view(m,n).dest;
        int drank = nghbr.h_view(m,n).rank;
        if (drank != my_rank) {
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int lid = nghbr.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);

          // get ptr to send buffer when neighbor is at coarser/same/fine level
          int data_size = nvar;
          if ( nghbr.h_view(m,n).lev < pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= sendbuf[n].icoar_ndat;
          } else if ( nghbr.h_view(m,n).lev == pmy_pack->pmb->mb_lev.h_view(m) ) {
            if (is_z4c) {
              data_size *= sendbuf[n].isame_z4c_ndat;
            } else {
              data_size *= sendbuf[n].isame_ndat;
            }
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
// \!fn void RecvBuffers()
// \brief Unpack boundary buffers

TaskStatus MeshBoundaryValuesCC::RecvAndUnpackCC(DvceArray5D<Real> &a,
                                                 DvceArray5D<Real> &ca,
                                                 int iwl, int iwu) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
  auto &is_z4c = is_z4c_;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  // skip the DIAGONAL transverse buffers entirely (see bvals.hpp)
  const bool skipd_ = skip_x2x3_diag;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  auto &multilevel = pmy_pack->pmesh->multilevel;
  // x1 index window: see PackAndSendCC and the declaration in bvals.hpp
  const int wlo_ = iwl, whi_ = iwu;
  const int mbis_ = pmy_pack->pmesh->mb_indcs.is;
  const int mbie_ = pmy_pack->pmesh->mb_indcs.ie;
#if MPI_PARALLEL_ENABLED
  //----- STEP 1: check that recv boundary buffer communications have all completed

  bool bflag = false;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0 &&
          !(use_cs && IsCubeVertexCorner(nghbr.h_view, mbpanel.h_view, m, n)) &&
          !(skipd_ && IsX2X3DiagSlot(n))) {
        if (nghbr.h_view(m,n).rank != global_variable::my_rank) {
          int test;
          int ierr = MPI_Test(&(rbuf[n].vars_req[m]), &test, MPI_STATUS_IGNORE);
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

  //----- STEP 2: buffers have all completed, so unpack

  int nvar = a.extent_int(1);  // TODO(@user): 2nd index from L of in array must be NVAR

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (nmb*nnghbr*nvar), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    // only unpack buffers when neighbor exists (cube vertex skipped -- see bvals.hpp)
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n)) &&
        !(skipd_ && IsX2X3DiagSlot(n))) {
      int il, iu, jl, ju, kl, ku;
      // if neighbor is at coarser level, use coar indices to unpack buffer
      if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
        il = rbuf[n].icoar[0].bis;
        iu = rbuf[n].icoar[0].bie;
        jl = rbuf[n].icoar[0].bjs;
        ju = rbuf[n].icoar[0].bje;
        kl = rbuf[n].icoar[0].bks;
        ku = rbuf[n].icoar[0].bke;
      // if neighbor is at same level, use same indices to unpack buffer
      } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
        il = rbuf[n].isame[0].bis;
        iu = rbuf[n].isame[0].bie;
        jl = rbuf[n].isame[0].bjs;
        ju = rbuf[n].isame[0].bje;
        kl = rbuf[n].isame[0].bks;
        ku = rbuf[n].isame[0].bke;
      // if neighbor is at finer level, use fine indices to unpack buffer
      } else {
        il = rbuf[n].ifine[0].bis;
        iu = rbuf[n].ifine[0].bie;
        jl = rbuf[n].ifine[0].bjs;
        ju = rbuf[n].ifine[0].bje;
        kl = rbuf[n].ifine[0].bks;
        ku = rbuf[n].ifine[0].bke;
      }
      int ni = iu - il + 1;
      int nj = ju - jl + 1;
      int nk = ku - kl + 1;
      int nkj  = nk*nj;
      // the window narrows the LOOP only: ni and the index map stay as they were
      int ilp = il, iup = iu;
      // same-level x2/x3 halos only: coarse/fine ranges differ between the two sides
      // and the cubed-sphere corner fill reads the full strips
      if (wlo_ >= 0 && !use_cs && nghbr.d_view(m,n).lev == mblev.d_view(m) &&
          il <= mbis_ && iu >= mbie_) {
        ilp = (il > wlo_) ? il : wlo_;
        iup = (iu < whi_) ? iu : whi_;
      }

      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
        int k = idx / nj;
        int j = (idx - k * nj) + jl;
        k += kl;

        // if neighbor is at same or finer level, load data directly into u0
        if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
          [&](const int i) {
            a(m,v,k,j,i) = rbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) );
          });

        // if neighbor is at coarser level, load data into coarse_u0
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
          [&](const int i) {
            ca(m,v,k,j,i) = rbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) );
          });
        }
      });
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables).  Whole body is
  // under the host-side (is_z4c && multilevel): skip the launch when it cannot fire.
  if (is_z4c_ && ml_) {
  Kokkos::parallel_for("RecvBuffZ4c", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);
    // only unpack buffers when neighbor exists (cube vertex skipped -- see bvals.hpp)
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n)) &&
        !(skipd_ && IsX2X3DiagSlot(n))) {
      int il, iu, jl, ju, kl, ku;
      // If neighbor is at same level and data is for Z4c module, unpack data from coarse
      // array for higher-order prolongation
      if ((nghbr.d_view(m,n).lev == mblev.d_view(m)) && (is_z4c) && (multilevel)) {
        il = rbuf[n].isame_z4c.bis;
        iu = rbuf[n].isame_z4c.bie;
        jl = rbuf[n].isame_z4c.bjs;
        ju = rbuf[n].isame_z4c.bje;
        kl = rbuf[n].isame_z4c.bks;
        ku = rbuf[n].isame_z4c.bke;
        int ni = iu - il + 1;
        int nj = ju - jl + 1;
        int nk = ku - kl + 1;
        int nkj  = nk*nj;
        int ndat = nvar*rbuf[n].isame_ndat; // size of same level data packed in buff

        // Middle loop over k,j
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;

          // load data into coarse_u0
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
          [&](const int i) {
            ca(m,v,k,j,i) = rbuf[n].vars(m,ndat + (i-il + ni*(j-jl + nj*(k-kl + nk*v))) );
          });
        });
      }
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer
  }

  // Every face buffer is unpacked by this point, which is what the corner fill needs.
  // The COARSE array needs it too -- see the note in bvals_fc.cpp's RecvAndUnpackFC.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    FillPanelCornersCC(a, false);
    if (pmy_pack->pmesh->multilevel) { FillPanelCornersCC(ca, true); }
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
// \!fn void MeshBoundaryValuesCC::FillPanelCornersCC()
// \brief Fill the ng x ng corner ghost cells of a cubed-sphere panel corner.
//
// The cell-centred twin of MeshBoundaryValuesFC::FillPanelCornersFC, and needed for the
// same reason: a panel corner is a CUBE VERTEX where only THREE panels meet, so the
// generic corner buffer -- a rectangular index region of one diagonal neighbour under the
// face seam's signed permutation -- reaches somewhere meaningless. The value copied there
// is wrong by O(1) in EVERY variable, density included, so it is not a basis error that
// TransformMomentum could repair.
//
// For hydro this was measured not to matter: a dimensionally split PLM+HLLC sweep never
// reconstructs through the ng x ng diagonal block. MHD does. The edge-centred EMF at a
// panel-corner edge is built from fluxes on the faces meeting there, whose reconstruction
// does read the diagonal cells, so O(1) garbage in these ghosts becomes an O(1) EMF on a
// corner edge and the CT update then drives the corner field away in a single step. That
// is what made the uniform-field gate -- an exact static state -- lose 1.3% of its
// magnetic energy on cycle 1 and blow up on cycle 2.
//
// The panel's own gnomonic map is well defined out there and the two flanking face halos
// are accurate, so each corner ghost is extrapolated quadratically from one, again from
// the other, and the two are averaged. The stencil reads only the flanking strips, never
// the corner block being written. Same-panel corners are left entirely alone.
//
// EACH EXTRAPOLANT GETS A SIGN-PRESERVING FALLBACK. The raw quadratic weights are not
// positivity preserving -- with nghost = 3 the third ghost carries (10, -15, 6) -- and in
// the red-giant run the density extrapolant went NEGATIVE in the layer above the
// star-corona join, firing the density and energy floors ~1e6 times per output interval
// in these ghost cells. The rule applied here is the weakest one that rules that out: if
// the three stencil nodes all share a STRICT sign and the extrapolant comes out with the
// opposite sign, or exactly zero, the extrapolant is replaced by the node ADJACENT to the
// ghost -- it falls back to constant extrapolation there. In every other case it is kept
// bit for bit. A positive density or internal energy above the red-giant star/corona join
// therefore stays positive, while a smooth profile, where the quadratic never crosses
// zero, is untouched, so the angular diffusion operator that reads this block through the
// conduction cross term keeps its order. An earlier version instead CLAMPED every
// extrapolant to the [min,max] of its three nodes; that also truncates smooth MONOTONE
// data, which is what that operator is built from, and it pushed the nx = 32 Linf
// residual of tst/test_suite/rad/test_rad_cs_raddiff_cpu.py from 0.0363 to 0.189 against
// a 0.1 gate. The sign fallback leaves that test where it was.

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

void MeshBoundaryValuesCC::FillPanelCornersCC(DvceArray5D<Real> &a, bool coarse) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = coarse ? indcs.cjs : indcs.js;
  const int je = coarse ? indcs.cje : indcs.je;
  const int ks = coarse ? indcs.cks : indcs.ks;
  const int ke = coarse ? indcs.cke : indcs.ke;
  const int ng = indcs.ng;
  const int nmb = pmy_pack->nmb_thispack;
  const int nvar = a.extent_int(1);
  const int n1 = a.extent_int(4);
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  const bool poison_ = pmy_pack->pmesh->cs_corner_poison;
  auto a_ = a;

  // par_for takes at most five ranges, so the ng x ng corner block is flattened into g
  par_for("cs_fill_corners_cc", DevExeSpace(), 0,(nmb-1), 0,3, 0,(nvar-1), 0,(n1-1),
          0,(ng*ng-1),
  KOKKOS_LAMBDA(const int m, const int c, const int v, const int i, const int g) {
    const int gj = g/ng;
    const int gk = g - gj*ng;
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

    // Quadratic Lagrange extrapolated d cells beyond an anchor, nodes 0,1,2 stepping
    // inward: w = ((d+1)(d+2)/2, -d(d+2), d(d+1)/2), which sums to 1.
    const Real dj = static_cast<Real>(gj + 1);
    const Real dk = static_cast<Real>(gk + 1);
    const Real wj0 = 0.5*(dj+1.0)*(dj+2.0), wj1 = -dj*(dj+2.0), wj2 = 0.5*dj*(dj+1.0);
    const Real wk0 = 0.5*(dk+1.0)*(dk+2.0), wk1 = -dk*(dk+2.0), wk2 = 0.5*dk*(dk+1.0);
    const int stj = -sj, stk = -sk;       // step INWARD from the anchor
    const int jt = (sj < 0) ? (js-1-gj) : (je+1+gj);
    const int kt = (sk < 0) ? (ks-1-gk) : (ke+1+gk);
    const int aj = (sj < 0) ? js : je;
    const int ak = (sk < 0) ? ks : ke;

    const Real k0 = a_(m,v,ak,jt,i);
    const Real k1 = a_(m,v,ak+stk,jt,i);
    const Real k2 = a_(m,v,ak+2*stk,jt,i);
    const Real j0 = a_(m,v,kt,aj,i);
    const Real j1 = a_(m,v,kt,aj+stj,i);
    const Real j2 = a_(m,v,kt,aj+2*stj,i);
    Real ek = wk0*k0 + wk1*k1 + wk2*k2;
    Real ej = wj0*j0 + wj1*j1 + wj2*j2;
    // SIGN-PRESERVING FALLBACK -- see the note above the function.  It only fires when
    // the three nodes share a strict sign and the extrapolant flips it (or zeroes it),
    // and then drops to the node next to the ghost; smooth data is kept exactly.
    ek = SignPreserveCS(ek, k0, k1, k2);
    ej = SignPreserveCS(ej, j0, j1, j2);
    a_(m,v,kt,jt,i) = 0.5*(ek + ej);
    if (poison_) { a_(m,v,kt,jt,i) = 1.0e30; }   // DEBUG: see <mesh>/cs_corner_poison
  });
  return;
}

