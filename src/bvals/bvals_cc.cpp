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
                                               DvceArray5D<Real> &ca) {
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
  const bool ml_ = pmy_pack->pmesh->multilevel;
  const bool use_pole = pmy_pack->pmesh->use_polar_boundary;
  // needed only on the cubed sphere, to give a source cell its (xi,eta)
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &cs_indcs = pmy_pack->pmesh->mb_indcs;
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;
  auto &is_z4c = is_z4c_;
  auto &multilevel = pmy_pack->pmesh->multilevel;
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
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
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

          // Which buffer is this? A same-level FACE buffer is ng cells deep in its own
          // normal direction and spans the full active range in the other tangential
          // one; that is the only case the along-seam resample applies to. The edge and
          // corner buffers (ng x ng) are left as a plain copy -- a dimensionally split
          // PLM+HLLC sweep never reconstructs through them.
          const int ngh_ = cs_indcs.ng;
          if (nj == ngh_ && nk > ngh_) {
            cs_seam = 2;
          } else if (nk == ngh_ && nj > ngh_) {
            cs_seam = 3;
          }
        } else if (do_pole) {
          aj = -1;
          bj = jl + ju;
          if (v == IVY) signvar = -1;
          if (v == IVZ) signvar = -1;
        }

        // Value of variable v at the source cell, with the seam transform applied. Both
        // tangential momenta of the source cell are needed to produce either one, so
        // this reads two components and selects; IVX and every scalar take the plain
        // copy path. Away from a panel seam it is exactly the old `a(...)*signvar`.
        auto srcval = [&](const int kk, const int jj, const int i) {
          if (!cs_xform) {
            return (cs_coar ? ca(m,vv,kk,jj,i) : a(m,vv,kk,jj,i))*signvar;
          }
          const int js_ = cs_coar ? cs_indcs.cjs : cs_indcs.js;
          const int ks_ = cs_coar ? cs_indcs.cks : cs_indcs.ks;
          const int nx2_ = cs_coar ? cs_indcs.cnx2 : cs_indcs.nx2;
          const int nx3_ = cs_coar ? cs_indcs.cnx3 : cs_indcs.nx3;
          const Real xi = 0.25*M_PI*CellCenterX(jj-js_, nx2_,
                            mbsize.d_view(m).x2min, mbsize.d_view(m).x2max);
          const Real eta = 0.25*M_PI*CellCenterX(kk-ks_, nx3_,
                            mbsize.d_view(m).x3min, mbsize.d_view(m).x3max);
          Real m2o, m3o;
          const Real my_ = cs_coar ? ca(m,IVY,kk,jj,i) : a(m,IVY,kk,jj,i);
          const Real mz_ = cs_coar ? ca(m,IVZ,kk,jj,i) : a(m,IVZ,kk,jj,i);
          cubed_sphere::TransformMomentum(cs_srcpanel, cs_dstpanel, xi, eta,
                                          my_, mz_, m2o, m3o);
          return (v == IVY) ? m2o : m3o;
        };

        // ALONG-SEAM RESAMPLE. Across a panel seam the two charts share the seam-normal
        // coordinate exactly but NOT the seam-parallel one. Writing the seam-normal
        // angle of a source cell as n and its seam-parallel angle as a, the physical
        // point of that cell sits at seam-parallel angle atan(tan(a)/tan|n|) in the
        // DESTINATION chart, not at a. So the plain index copy hands each ghost cell the
        // state of a point up to half a cell (layer 0) or ~1.5 cells (layer 1) away
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
        auto seamval = [&](const int kk, const int jj, const int i) {
          if (cs_seam == 0) return srcval(kk,jj,i);
          const int js_ = cs_coar ? cs_indcs.cjs : cs_indcs.js;
          const int ks_ = cs_coar ? cs_indcs.cks : cs_indcs.ks;
          const int nx2_ = cs_coar ? cs_indcs.cnx2 : cs_indcs.nx2;
          const int nx3_ = cs_coar ? cs_indcs.cnx3 : cs_indcs.nx3;
          const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
          const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
          const Real xi  = 0.25*M_PI*CellCenterX(jj-js_, nx2_, x2mn, x2mx);
          const Real eta = 0.25*M_PI*CellCenterX(kk-ks_, nx3_, x3mn, x3mx);
          Real ang, nrm, dang;
          int sc, blo, bhi;
          if (cs_seam == 2) {
            ang = eta; nrm = xi;
            dang = 0.25*M_PI*(x3mx - x3mn)/static_cast<Real>(nx3_);
            sc = kk; blo = kl; bhi = ku - 2;
          } else {
            ang = xi; nrm = eta;
            dang = 0.25*M_PI*(x2mx - x2mn)/static_cast<Real>(nx2_);
            sc = jj; blo = jl; bhi = ju - 2;
          }
          const Real pos = sc + (atan(tan(ang)*tan(fabs(nrm))) - ang)/dang;
          int b = static_cast<int>(floor(pos + 0.5)) - 1;
          b = (b < blo) ? blo : ((b > bhi) ? bhi : b);
          const Real u = pos - static_cast<Real>(b + 1);
          const Real wm = 0.5*u*(u - 1.0);
          const Real w0 = 1.0 - u*u;
          const Real wp = 0.5*u*(u + 1.0);
          if (cs_seam == 2) {
            return wm*srcval(b,jj,i) + w0*srcval(b+1,jj,i) + wp*srcval(b+2,jj,i);
          }
          return wm*srcval(kk,b,i) + w0*srcval(kk,b+1,i) + wp*srcval(kk,b+2,i);
        };

        // Middle loop over k,j
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;
          int kk = ak*k + bk;
          int jj = aj*j + bj;

          // Inner (vector) loop over i
          // copy directly into recv buffer if MeshBlocks on same rank

          // seamval() reads u0 or coarse_u0 according to cs_coar, so one expression
          // serves a neighbour at ANY level.
          if (nghbr.d_view(m,n).rank == my_rank) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(kk,jj,i);
              int index = i-il + ni*(sj*(j-jl) + sk*(k-kl) + nk*nj*v);
              rbuf[dn].vars(dm, index) = val;
            });

          // else copy into send buffer for MPI communication below

          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              Real val = seamval(kk,jj,i);
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
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              rbuf[dn].vars(dm, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = ca(m,v,k,j,i);
            });
          }

        // else copy into send buffer for MPI communication below

        } else {
          // if neighbor is at same or finer level, load data from u0
          if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              sbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) ) = a(m,v,k,j,i);
            });
          // if neighbor is at coarser level, load data from coarse_u0
          } else {
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
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

  Kokkos::parallel_for("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    // only load buffers when neighbor exists.  CUBED-SPHERE CUBE VERTEX: skipped on
    // both sides; FillPanelCornersCC overwrites exactly this corner block. See bvals.hpp.
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
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
                                                 DvceArray5D<Real> &ca) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
  auto &is_z4c = is_z4c_;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  const bool use_cs = pmy_pack->pmesh->use_cubed_sphere;
  const bool ml_ = pmy_pack->pmesh->multilevel;
  auto &multilevel = pmy_pack->pmesh->multilevel;
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
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
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

      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
        int k = idx / nj;
        int j = (idx - k * nj) + jl;
        k += kl;

        // if neighbor is at same or finer level, load data directly into u0
        if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
          [&](const int i) {
            a(m,v,k,j,i) = rbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) );
          });

        // if neighbor is at coarser level, load data into coarse_u0
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
          [&](const int i) {
            ca(m,v,k,j,i) = rbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*v))) );
          });
        }
      });
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables)
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);
    // only unpack buffers when neighbor exists (cube vertex skipped -- see bvals.hpp)
    if (nghbr.d_view(m,n).gid >= 0 &&
        !(use_cs && IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m, n))) {
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

  // Every face buffer is unpacked by this point, which is what the corner fill needs.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    FillPanelCornersCC(a);
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

void MeshBoundaryValuesCC::FillPanelCornersCC(DvceArray5D<Real> &a) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int ng = indcs.ng;
  const int nmb = pmy_pack->nmb_thispack;
  const int nvar = a.extent_int(1);
  const int n1 = a.extent_int(4);
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mblev = pmy_pack->pmb->mb_lev;
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
    bool seam = false;
    if (nghbr.d_view(m,nj_id).gid >= 0 && nghbr.d_view(m,nj_id).panel != mp) seam = true;
    if (nghbr.d_view(m,nk_id).gid >= 0 && nghbr.d_view(m,nk_id).panel != mp) seam = true;
    if (!seam) return;

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

    const Real ek = wk0*a_(m,v,ak,jt,i) + wk1*a_(m,v,ak+stk,jt,i)
                  + wk2*a_(m,v,ak+2*stk,jt,i);
    const Real ej = wj0*a_(m,v,kt,aj,i) + wj1*a_(m,v,kt,aj+stj,i)
                  + wj2*a_(m,v,kt,aj+2*stj,i);
    a_(m,v,kt,jt,i) = 0.5*(ek + ej);
  });
  return;
}

