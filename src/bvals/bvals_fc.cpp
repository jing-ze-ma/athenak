//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_fc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for face-centered (FC)
//! Mesh variables.
//! Prolongation of FC variables  occurs in ProlongateFC() function called from task list

#include <cstdlib>
#include <iostream>
#include <utility>
#include <iomanip>    // std::setprecision()

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
// BValFC constructor:

MeshBoundaryValuesFC::MeshBoundaryValuesFC(MeshBlockPack *pp, ParameterInput *pin) :
  MeshBoundaryValues(pp, pin, false) {
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

  {int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  const bool use_pole = pmy_pack->pmesh->use_polar_boundary;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(three field components)
  int nmnv = 3*nmb;
  Kokkos::TeamPolicy<> policy(DevExeSpace(), nmnv, Kokkos::AUTO);
  Kokkos::parallel_for("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/3;
    const int v = tmember.league_rank()%3;

    // scalar loop over neighbors to prevent race condition in overlapping assignments
    for (int n=0; n<nnghbr; ++n) {
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
            const int nx2_ = cs_coar ? cs_indcs.cnx2 : cs_indcs.nx2;
            const int nx3_ = cs_coar ? cs_indcs.cnx3 : cs_indcs.nx3;
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
            // b.x3f at (eta face kf, xi FACE jf), interpolated across xi
            auto x3f_at_xiface = [&](const int kf, const int jf, const int i) {
              if (jf-1 < js_) {
                return 1.5*bx3(kf,js_,i) - 0.5*bx3(kf,js_+1,i);
              } else if (jf > je_) {
                return 1.5*bx3(kf,je_,i) - 0.5*bx3(kf,je_-1,i);
              }
              return 0.5*(bx3(kf,jf-1,i) + bx3(kf,jf,i));
            };
            // b.x2f at (eta FACE kf, xi face jf), interpolated across eta
            auto x2f_at_etaface = [&](const int kf, const int jf, const int i) {
              if (kf-1 < ks_) {
                return 1.5*bx2(ks_,jf,i) - 0.5*bx2(ks_+1,jf,i);
              } else if (kf > ke_) {
                return 1.5*bx2(ke_,jf,i) - 0.5*bx2(ke_-1,jf,i);
              }
              return 0.5*(bx2(kf-1,jf,i) + bx2(kf,jf,i));
            };
            //
            // The (xi,eta) of a source sample, with the STAGGERING of component vv:
            // b.x2f sits on a xi face and b.x3f on an eta face, b.x1f on neither.
            const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
            const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
            auto angles = [&](const int kk, const int jj, Real &xi, Real &eta) {
              xi  = 0.25*M_PI*((vv == 1) ? LeftEdgeX(jj-js_, nx2_, x2mn, x2mx)
                                         : CellCenterX(jj-js_, nx2_, x2mn, x2mx));
              eta = 0.25*M_PI*((vv == 2) ? LeftEdgeX(kk-ks_, nx3_, x3mn, x3mx)
                                         : CellCenterX(kk-ks_, nx3_, x3mn, x3mx));
            };
            auto srcval = [&](const int kk, const int jj, const int i) {
              if (!cs_xform) {
                if (vv == 0) return bx1(kk,jj,i)*signvar;
                if (vv == 1) return bx2(kk,jj,i)*signvar;
                return bx3(kk,jj,i)*signvar;
              }
              Real xi, eta, bxi, bet;
              angles(kk, jj, xi, eta);
              if (vv == 1) {
                // primary on a XI face: (xi face jj, eta centre kk)
                bxi = bx2(kk,jj,i);
                bet = 0.5*(x3f_at_xiface(kk,jj,i) + x3f_at_xiface(kk+1,jj,i));
              } else {
                // primary on an ETA face: (xi centre jj, eta face kk)
                bet = bx3(kk,jj,i);
                bxi = 0.5*(x2f_at_etaface(kk,jj,i) + x2f_at_etaface(kk,jj+1,i));
              }
              Real oxi, oet;
              cubed_sphere::TransformFieldToDstNormals(cs_srcpanel, cs_dstpanel, xi, eta,
                                                       bxi, bet, oxi, oet);
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
              Real xi, eta;
              angles(kk, jj, xi, eta);
              Real ang, nrm, dang;
              int sc, blo, bhi;
              if (cs_seam == 2) {
                ang = eta; nrm = xi;
                dang = 0.25*M_PI*(x3mx - x3mn)/static_cast<Real>(nx3_);
                // Bounds from the SOURCE'S ACTIVE range, not the buffer's.  For a face
                // or x1-edge buffer the two coincide, so this is a no-op there; for a
                // doubly-ghost buffer the along-seam direction is only ng deep in the
                // BUFFER and clamping to that would extrapolate from two cells when the
                // source block holds the data.
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
              const Real wm = 0.5*u*(u - 1.0);
              const Real w0 = 1.0 - u*u;
              const Real wp = 0.5*u*(u + 1.0);
              if (cs_seam == 2) {
                return wm*srcval(bs,jj,i) + w0*srcval(bs+1,jj,i) + wp*srcval(bs+2,jj,i);
              }
              return wm*srcval(kk,bs,i) + w0*srcval(kk,bs+1,i) + wp*srcval(kk,bs+2,i);
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
      tmember.team_barrier();
    }
  }); // end par_for_outer
  }

#if MPI_PARALLEL_ENABLED
  // Send boundary buffer to neighboring MeshBlocks using MPI
  Kokkos::fence();
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
    {
      const Real ek = wk0*b_.x1f(m,akc,jtc,i) + wk1*b_.x1f(m,akc+stk,jtc,i)
                    + wk2*b_.x1f(m,akc+2*stk,jtc,i);
      const Real ej = wj0*b_.x1f(m,ktc,ajc,i) + wj1*b_.x1f(m,ktc,ajc+stj,i)
                    + wj2*b_.x1f(m,ktc,ajc+2*stj,i);
      b_.x1f(m,ktc,jtc,i) = wblend*ek + (1.0-wblend)*ej;
    }
    if (i > (nc1-1)) return;
    // b.x2f: FACE index in j, cell index in k
    {
      const Real ek = wk0*b_.x2f(m,akc,jtf,i) + wk1*b_.x2f(m,akc+stk,jtf,i)
                    + wk2*b_.x2f(m,akc+2*stk,jtf,i);
      const Real ej = wj0*b_.x2f(m,ktc,ajf,i) + wj1*b_.x2f(m,ktc,ajf+stj,i)
                    + wj2*b_.x2f(m,ktc,ajf+2*stj,i);
      b_.x2f(m,ktc,jtf,i) = wblend*ek + (1.0-wblend)*ej;
    }
    // b.x3f: cell index in j, FACE index in k
    {
      const Real ek = wk0*b_.x3f(m,akf,jtc,i) + wk1*b_.x3f(m,akf+stk,jtc,i)
                    + wk2*b_.x3f(m,akf+2*stk,jtc,i);
      const Real ej = wj0*b_.x3f(m,ktf,ajc,i) + wj1*b_.x3f(m,ktf,ajc+stj,i)
                    + wj2*b_.x3f(m,ktf,ajc+2*stj,i);
      b_.x3f(m,ktf,jtc,i) = wblend*ek + (1.0-wblend)*ej;
    }
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
  auto &rbuf = recvbuf;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  auto &mblev = pmy_pack->pmb->mb_lev;
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

  //----- STEP 2: buffers have all completed, so unpack 3-components of field

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(three field components)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/3;
    const int v = tmember.league_rank()%3;

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
        const int nk = ku - kl + 1;
        const int nkji = nk*nj*ni;
        const int nji  = nj*ni;

        // if neighbor is at same or finer level, load data directly into b0
        if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkji),
          [&](const int idx) {
            int k = (idx)/nji;
            int j = (idx - k*nji)/ni;
            int i = (idx - k*nji - j*ni) + il;
            k += kl;
            j += jl;
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
            int j = (idx - k*nji)/ni;
            int i = (idx - k*nji - j*ni) + il;
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

  return TaskStatus::complete;
}
