//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file bvals_cc_pack.inc
//! \brief d6fast4-1112: the body of the PackAndSendCC team kernel (pack_, moved
//! verbatim from bvals_cc.cpp).  With BVCC_PACK_KSEAM / BVCC_PACK_KWB defined (CUDA)
//! kseam_ and kwb_ are compile-time constants inside, so nvcc drops the dead seam /
//! wb_seam branches as the std::true_type / false_type arguments do on HIP; otherwise
//! the generic (HIP) signature.
//! Needs BVCC_PACK_NAME defined.
#if defined(BVCC_PACK_KSEAM)
  auto BVCC_PACK_NAME = KOKKOS_LAMBDA(const TeamMember_t &tmember) {
    constexpr bool kseam_ = BVCC_PACK_KSEAM;
    constexpr bool kwb_ = BVCC_PACK_KWB;
#else
  auto BVCC_PACK_NAME = KOKKOS_LAMBDA(const TeamMember_t &tmember, const auto kseam_,
                             const auto kwb_) {
#endif
    const int lr_ = tmember.league_rank()/nsp_;
    const int sp_ = tmember.league_rank() - lr_*nsp_;
    const int m = (lr_)/(nnghbr*nvar);
    const int n = (lr_ - m*(nnghbr*nvar))/nvar;
    const int vl = (lr_ - m*(nnghbr*nvar) - n*nvar);
    const int v = vlst_.v[vl];

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
      // speed2-1031 bvsplit: this team's chunk of the (k, j) range
      const int q0_ = static_cast<int>((static_cast<int64_t>(nkj)*sp_)/nsp_);
      const int q1_ = static_cast<int>((static_cast<int64_t>(nkj)*(sp_ + 1))/nsp_);
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

      if (static_cast<bool>(kseam_) && (do_cs || do_pole)) {
        int aj = 1, bj = 0;
        int ak = 1, bk = 0;
        int signvar = 1;
        int sj = 1, sk = nj;
        int vv = v;
        bool cs_xform = false;
        // cs-s4-1017: the FACE-NORMAL pair transform of cs_perm_pairs (see below) and
        // the half-range quad of cs_hrq0
        bool cs_fxform = false;
        bool cs_hrquad = false;
        // THE COMPONENT ROLE of slot v: is it a member of a tangential vector pair, and
        // which slots hold the pair's x2 (va_) and x3 (vb_) members.  Without a table
        // (SetVectorPairs) this is the historical rule, slots (IVY, IVZ) of any array.
        bool vec_ = (v == IVY) || (v == IVZ);
        int va_ = IVY, vb_ = IVZ;
        if (vtab_) {
          const int r_ = vrole(v);
          vec_ = (r_ >= 0);
          if (vec_) {
            if ((r_ & 1) == 0) {
              va_ = v; vb_ = r_ >> 1;
            } else {
              va_ = r_ >> 1; vb_ = v;
            }
          }
        }
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
          cs_xform = vec_;
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
          // cs_perm_pairs (default false): the tangential pairs of this object are
          // FACE-NORMAL components (v.nhat_xi, v.nhat_eta: the implicit M1 work array,
          // stage CS2), not covariant ones.  cs-s4-1017: they are re-expressed on the
          // DESTINATION panel's face normals at the source cell's own point
          // (cubed_sphere::TransformFieldToDstNormals: normals -> contravariant ->
          // Cartesian -> destination normals), the face-centred field's transform.
          // Until cs-s4-1017 they took the signed AXIS permutation of the index map
          // instead: exact for the SEAM-NORMAL component (the seam plane is common to
          // both charts), but the ALONG-SEAM one came out in the source panel's mirrored
          // basis, O(1) wrong where the vector has a seam-normal part (rad_m1
          // cs_vec_check: 0.95-1.07 of |v| at N = 8..64, not converging).  No cs
          // operator reads the along-seam component of a seam ghost today, so the change
          // is round-off in what runs (the seam-normal component, now through the
          // arithmetic).  These
          // objects run without the along-seam resample (cs_noresample), so the
          // transform is taken at the single source cell (cs_seam = 0 below).
          if (cperm_ && vec_) {
            cs_xform = false;
            cs_fxform = true;
          }
          // cs-s4-1017: the half-range quad (cs_hrq0, see bvals.hpp)
          if (hrq_ >= 0 && v >= hrq_ && v < hrq_ + 4) {
            cs_xform = false;
            cs_fxform = false;
            cs_hrquad = true;
            const int ad = (v - hrq_)/2, pm = (v - hrq_) % 2;
            const bool sw = (pb.swap_ax == 1);
            const int as = sw ? (1 - ad) : ad;
            const int rva = sw ? rev_b_preswap : rev_a_preswap;
            const int rvb = sw ? rev_a_preswap : rev_b_preswap;
            const bool rv = ((ad == 0) ? rva : rvb) != 0;
            vv = hrq_ + 2*as + (rv ? (1 - pm) : pm);
            signvar = rv ? -1 : 1;
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
          if (no_rs_cc || nors_) { cs_seam = 0; }
          if (cs_fxform || cs_hrquad) { cs_seam = 0; }
        } else if (do_pole) {
          aj = -1;
          bj = jl + ju;
          if (vec_) signvar = -1;
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
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, q0_, q1_), [&](const int idx) {
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
          Real wpos = 0.0;    // wball-1104 wb_seam: the target's along-seam index position
          int winw = 0;       //   and the inward seam-normal step in the source chart
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
            // layoutbug-1103 (<mesh>/cs_seam_blockend = extend, second pass): at an
            // INTERIOR block end along the seam the stencil continues into this block's
            // own along-seam ghosts (same-panel copies of the next block's active cells,
            // fresh after the first pass), up to ng cells, so the stencil is the one a
            // single block per panel would use: clamped to the PANEL, not the block.
            if (ext_ && !cs_coar) {
              const Real w_ = (cs_seam == 2) ? (x3mx - x3mn) : (x2mx - x2mn);
              const Real x0_ = (cs_seam == 2) ? x3mn : x2mn;
              const int nbp_ = static_cast<int>(floor(2.0/w_ + 0.5));
              const int ib_ = static_cast<int>(floor((x0_ + 1.0)/w_ + 0.5));
              if (ib_ > 0) {blo -= ngx_;}
              if (ib_ < nbp_ - 1) {bhi += ngx_;}
            }
            // THE CLAMP IS TO THE SOURCE'S ACTIVE RANGE, NOT THE BUFFER'S.  For a face or
            // x1-edge buffer the two coincide (the along-seam extent of such a buffer IS
            // the active range), so this is bitwise for everything that resampled before;
            // an x2x3 edge buffer is only ng deep along the seam and clamping to that
            // would extrapolate from two cells while the source block holds the data.
            // The stencil reads the array directly, not the buffer, so any active cell is
            // available.  bvals_fc.cpp says the same thing about its own stencil table.
            // layoutfix-1106: cs_seam_pos_global, the offset of this block in the panel
            int spn_ = 0;
            if (posg_) {
              const Real wp_ = (cs_seam == 2) ? (x3mx - x3mn) : (x2mx - x2mn);
              const Real xp_ = (cs_seam == 2) ? x3mn : x2mn;
              spn_ = static_cast<int>(floor((xp_ + 1.0)/wp_ + 0.5))
                     *((cs_seam == 2) ? nx3_ : nx2_);
            }
            const Real pos = cs_seam::SeamPos(sc, spn_,
                                              (atan(tan(ang)*tan(fabs(nrm))) - ang)/dang);
            wpos = pos;
            winw = (nrm > 0.0) ? -1 : 1;
            int b = cs_seam::SeamRound(pos) - 1;   // fc2vtx-1103: tie-safe
            b = (b < blo) ? blo : ((b > bhi) ? bhi : b);
            const Real u = pos - static_cast<Real>(b + 1);
            cs_interp = (u >= -1.0) && (u <= 1.0);
            wm = 0.5*u*(u - 1.0);
            w0 = 1.0 - u*u;
            wp = 0.5*u*(u + 1.0);
            if (rs_lin && cs_interp) {
              // <mesh>/cs_seam_resample = linear: the two bracketing donors, weights >= 0
              wm = (u < 0.0) ? -u : 0.0;
              wp = (u > 0.0) ? u : 0.0;
              w0 = 1.0 - wm - wp;
            }
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
          // cs-s4-1017: the face-normal pair (cs_perm_pairs), one source cell (nst = 1)
          cubed_sphere::SeamFieldXform fxf;
          if (cs_fxform) {
            const Real xis = 0.25*M_PI*CellCenterX(jst[0]-js_, nx2_, x2mn, x2mx);
            const Real etas = 0.25*M_PI*CellCenterX(kst[0]-ks_, nx3_, x3mn, x3mx);
            cubed_sphere::SeamFieldXformAt(cs_srcpanel, cs_dstpanel, xis, etas, fxf);
          }

          // Value of variable v at stencil cell s, radial index i. Both tangential
          // momenta of the source cell are needed to make either one, so this reads two
          // components and selects. Away from a panel seam it is exactly the old
          // `a(...)*signvar`. Only loads and the O(1) half of the transform remain here.
          auto sval = [&](const int s, const int i) {
            const int kq = kst[s], jq = jst[s];
            if (cs_fxform) {
              Real f2o, f3o;
              const Real fy_ = cs_coar ? ca(m,va_,kq,jq,i) : a(m,va_,kq,jq,i);
              const Real fz_ = cs_coar ? ca(m,vb_,kq,jq,i) : a(m,vb_,kq,jq,i);
              cubed_sphere::ApplyFieldXform(fxf, fy_, fz_, f2o, f3o);
              return (v == va_) ? f2o : f3o;
            }
            if (!cs_xform) {
              return (cs_coar ? ca(m,vv,kq,jq,i) : a(m,vv,kq,jq,i))*signvar;
            }
            Real m2o, m3o;
            const Real my_ = cs_coar ? ca(m,va_,kq,jq,i) : a(m,va_,kq,jq,i);
            const Real mz_ = cs_coar ? ca(m,vb_,kq,jq,i) : a(m,vb_,kq,jq,i);
            cubed_sphere::ApplyMomentumXform(xf[s], my_, mz_, m2o, m3o);
            return (v == va_) ? m2o : m3o;
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
          // <mesh>/cs_seam_rho_guard (see mesh.cpp): the plain copy of the source cell,
          // with its own basis transform, when the resampled DENSITY and the source
          // cell's density differ by more than the factor rg_.  Off (rg_ = 0): bitwise.
          cubed_sphere::SeamXform xfc;
          const bool rgon = (rg_ > 0.0) && (cs_seam != 0) && (v >= 0);
          if (rgon && cs_xform) {
            const Real xic = 0.25*M_PI*CellCenterX(jj-js_, nx2_, x2mn, x2mx);
            const Real etac = 0.25*M_PI*CellCenterX(kk-ks_, nx3_, x3mn, x3mx);
            cubed_sphere::SeamXformAt(cs_srcpanel, cs_dstpanel, xic, etac, xfc);
          }
          auto rgbad = [&](const int i) {
            const Real d0 = cs_coar ? ca(m,0,kst[0],jst[0],i) : a(m,0,kst[0],jst[0],i);
            const Real d1 = cs_coar ? ca(m,0,kst[1],jst[1],i) : a(m,0,kst[1],jst[1],i);
            const Real d2 = cs_coar ? ca(m,0,kst[2],jst[2],i) : a(m,0,kst[2],jst[2],i);
            Real dq = wm*d0 + w0*d1 + wp*d2;
            if (cs_interp) {
              dq = fmin(fmax(d0, fmax(d1, d2)), fmax(fmin(d0, fmin(d1, d2)), dq));
            }
            const Real dc = cs_coar ? ca(m,0,kk,jj,i) : a(m,0,kk,jj,i);
            return (dq > rg_*dc) || (dc > rg_*dq);
          };
          auto plainval = [&](const int i) {
            if (!cs_xform) {
              return (cs_coar ? ca(m,vv,kk,jj,i) : a(m,vv,kk,jj,i))*signvar;
            }
            Real m2o, m3o;
            const Real my_ = cs_coar ? ca(m,va_,kk,jj,i) : a(m,va_,kk,jj,i);
            const Real mz_ = cs_coar ? ca(m,vb_,kk,jj,i) : a(m,vb_,kk,jj,i);
            cubed_sphere::ApplyMomentumXform(xfc, my_, mz_, m2o, m3o);
            return (v == va_) ? m2o : m3o;
          };
          // wball-1104 wb_seam: the state of stencil cell s walked to the ghost's potential
          // (density for v = 0, total energy for v = wbsien_); see bvals.hpp
          const bool wbsv = static_cast<bool>(kwb_) && wbs_ && (cs_seam != 0) &&
                            !cs_coar && (v == 0 || v == wbsien_);
          Real lw5[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
          int c5 = 0;
          if (wbsv) {
            c5 = static_cast<int>(floor(wpos + 0.5));
            const Real x = wpos - static_cast<Real>(c5);
            for (int q = 0; q < 5; ++q) {
              Real l = 1.0;
              for (int r = 0; r < 5; ++r) {
                if (r != q) l *= (x - static_cast<Real>(r - 2))/static_cast<Real>(q - r);
              }
              lw5[q] = l;
            }
          }
          auto wbval = [&](const int s, const int i) {
            const int kq = kst[s], jq = jst[s];
            Real phig = 0.0;
            const int nA = (cs_seam == 2) ? wbsphi_.extent_int(1) : wbsphi_.extent_int(2);
            for (int q = 0; q < 5; ++q) {
              int ix = c5 - 2 + q;
              ix = (ix < 0) ? 0 : ((ix > nA - 1) ? nA - 1 : ix);
              phig += lw5[q]*((cs_seam == 2) ? wbsphi_(m,ix,jj,i) : wbsphi_(m,kk,ix,i));
            }
            auto prim = [&](const int kc, const int jc, Real &d, Real &e, Real &ke,
                            Real &ph) {
              d = a(m,0,kc,jc,i);
              const Real m1 = a(m,1,kc,jc,i), m2 = a(m,2,kc,jc,i), m3 = a(m,3,kc,jc,i);
              const Real xs = 0.25*M_PI*CellCenterX(jc-js_, nx2_, x2mn, x2mx);
              const Real es = 0.25*M_PI*CellCenterX(kc-ks_, nx3_, x3mn, x3mx);
              const Real cc = -sin(xs)*sin(es);
              ke = 0.5*(m1*m1 + (m2*m2 - 2.0*cc*m2*m3 + m3*m3)/(1.0 - cc*cc))/d;
              ph = wbsphi_(m,kc,jc,i);
              e = a(m,wbsien_,kc,jc,i) - ke - (wbsetg_ ? d*ph : 0.0);
              if (wbsmhd_) {
                const Real b1 = wbsbcc_(m,0,kc,jc,i), b2 = wbsbcc_(m,1,kc,jc,i);
                const Real b3 = wbsbcc_(m,2,kc,jc,i);
                ke += 0.5*(b1*b1 + b2*b2 + b3*b3);   // carried as part of E like KE
                e -= 0.5*(b1*b1 + b2*b2 + b3*b3);
              }
            };
            Real d, e, ke, ph;
            prim(kq, jq, d, e, ke, ph);
            Real dd = d, ee = e, tt = -1.0;
            if (wbsopt_ == 3) {
              const int ki = (cs_seam == 3) ? kq + winw : kq;
              const int ji = (cs_seam == 2) ? jq + winw : jq;
              Real di, ei, kei, phi_;
              prim(ki, ji, di, ei, kei, phi_);
              const Real ts = wbseos_.Temperature(d, e), ti = wbseos_.Temperature(di, ei);
              const Real dph = ph - phi_;
              const Real al = (fabs(dph) > 1.0e-12*(fabs(ph) + 1.0e-300)) ? (ts - ti)/dph : 0.0;
              WBAdvance(wbseos_, 3, d, e, phig - ph, dd, ee, tt, ts, al, ts);
            } else {
              WBAdvance(wbseos_, 1, d, e, phig - ph, dd, ee, tt);
            }
            if (v == 0) return dd;
            return ee + ke + (wbsetg_ ? dd*phig : 0.0);
          };
          auto seamval = [&](const int i) {
            if (cs_seam == 0) return sval(0,i);
            if (wbsv) {
              const Real sv0 = wbval(0,i), sv1 = wbval(1,i), sv2 = wbval(2,i);
              const Real qq = wm*sv0 + w0*sv1 + wp*sv2;
              if (!cs_interp) return qq;
              return fmin(fmax(sv0, fmax(sv1, sv2)), fmax(fmin(sv0, fmin(sv1, sv2)), qq));
            }
            if (rgon && rgbad(i)) return plainval(i);
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
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*vl);
              rbuf[dn].vars(dm, index) = val;
            });

          // else copy into send buffer for MPI communication below

          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(i);
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*vl);
              sbuf[n].vars(m,index) = val;
            });
          }
        });

      } else {   // normal boundary exchange
      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, q0_, q1_), [&](const int idx) {
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
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = ca(m,v,k,j,i);
            });
          }

        // else copy into send buffer for MPI communication below

        } else {
          // if neighbor is at same or finer level, load data from u0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) ) = ca(m,v,k,j,i);
            });
          }
        }
      });
      }  // end if-do-cs/do-pole block
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  };
