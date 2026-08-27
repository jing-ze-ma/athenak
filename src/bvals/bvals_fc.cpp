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
      // only load buffers when neighbor exists
      if (nghbr.d_view(m,n).gid >= 0) {
        // if neighbor is at coarser level, use cindices to pack buffer
        // Note indices can be different for each component of face-centered field.
        int il, iu, jl, ju, kl, ku, ndat;

        // WHICH COMPONENT'S INDEX RANGE DO WE READ THE SOURCE WITH?
        //
        // v is the RECEIVER's component: the packer writes into buffer slot v but reads
        // this block's array vv, and the two differ on a panel seam whose tangential axes
        // are interchanged (swap_ax) -- every equatorial-to-polar seam. The ranges are per
        // component and genuinely different: for a -x2 neighbour isame[1] (x2f) is ng
        // FACES in j by nx3 CELLS in k, while isame[2] (x3f) is ng CELLS in j by nx3+1
        // FACES in k. Reading array vv over isame[v] therefore took the wrong index set
        // AND the wrong NUMBER of elements, so the receiver unpacked a differently sized
        // region -- an O(1) error that did not converge under refinement. Read the source
        // with its OWN component's range; the sj/sk transposition below still lays it out
        // the way the receiver expects, and the counts match because nx2 == nx3 here.
        int vsrc = v;
        if (pmy_pack->pmesh->use_cubed_sphere &&
            nghbr.d_view(m,n).panel != mbpanel.d_view(m) && v > 0) {
          PanelBoundaries pb0;
          pb0 = pmy_pack->pmesh->GetPanelBoundary(mbpanel.d_view(m),
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
          
          const bool do_cs = pmy_pack->pmesh->use_cubed_sphere && (nghbr.d_view(m,n).panel != mbpanel.d_view(m));
          const bool do_pole = pmy_pack->pmesh->use_polar_boundary && (nghbr.d_view(m,n).polar > 0);
            
          if (do_cs || do_pole) {

            int aj = 1, bj = 0;
            int ak = 1, bk = 0;
            int signvar = 1;
            int sj = 1, sk = nj;
            int vv = v;
            bool cs_xform = false;
            int cs_srcpanel = 0, cs_dstpanel = 0;

            if (do_cs) {
              const auto ngh = nghbr.d_view(m,n);
              const int my_panel = mbpanel.d_view(m);
              PanelBoundaries pb;
              pb = pmy_pack->pmesh->GetPanelBoundary(my_panel, ngh.panel);

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
            // it lives and the other is averaged onto that location from its four
            // neighbours, which is second order. The radial component is untouched.
            const int js_ = cs_indcs.js, ks_ = cs_indcs.ks;
            auto seamval = [&](const int kk, const int jj, const int i) {
              if (!cs_xform) {
                if (vv == 0) return b.x1f(m,kk,jj,i)*signvar;
                if (vv == 1) return b.x2f(m,kk,jj,i)*signvar;
                return b.x3f(m,kk,jj,i)*signvar;
              }
              const Real x2mn = mbsize.d_view(m).x2min, x2mx = mbsize.d_view(m).x2max;
              const Real x3mn = mbsize.d_view(m).x3min, x3mx = mbsize.d_view(m).x3max;
              Real xi, eta, bxi, bet;
              if (vv == 1) {
                // primary on a XI face: (xi face jj, eta centre kk)
                xi  = 0.25*M_PI*LeftEdgeX(jj-js_, cs_indcs.nx2, x2mn, x2mx);
                eta = 0.25*M_PI*CellCenterX(kk-ks_, cs_indcs.nx3, x3mn, x3mx);
                bxi = b.x2f(m,kk,jj,i);
                bet = 0.25*(b.x3f(m,kk,jj-1,i) + b.x3f(m,kk,jj,i)
                          + b.x3f(m,kk+1,jj-1,i) + b.x3f(m,kk+1,jj,i));
              } else {
                // primary on an ETA face: (xi centre jj, eta face kk)
                xi  = 0.25*M_PI*CellCenterX(jj-js_, cs_indcs.nx2, x2mn, x2mx);
                eta = 0.25*M_PI*LeftEdgeX(kk-ks_, cs_indcs.nx3, x3mn, x3mx);
                bet = b.x3f(m,kk,jj,i);
                bxi = 0.25*(b.x2f(m,kk-1,jj,i) + b.x2f(m,kk-1,jj+1,i)
                          + b.x2f(m,kk,jj,i) + b.x2f(m,kk,jj+1,i));
              }
              Real oxi, oet;
              cubed_sphere::TransformFieldToDstNormals(cs_srcpanel, cs_dstpanel, xi, eta,
                                                       bxi, bet, oxi, oet);
              return (v == 1) ? oxi : oet;
            };

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
                int kk = ak*k + bk;
                int jj = aj*j + bj;
                rbuf[dn].vars(dm, ndat*v + i-il + ni*(sj*(j-jl) + sk*(k-kl)))
                  = seamval(kk,jj,i);
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
                int kk = ak*k + bk;
                int jj = aj*j + bj;
                sbuf[n].vars(m, ndat*v + i-il + ni*(sj*(j-jl) + sk*(k-kl)))
                  = seamval(kk,jj,i);
              });
            }
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
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0) {  // neighbor exists and not a physical boundary
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
// \!fn void RecvBuffers()
// \brief Unpack boundary buffers

TaskStatus MeshBoundaryValuesFC::RecvAndUnpackFC(DvceFaceFld4D<Real> &b,
                                                 DvceFaceFld4D<Real> &cb) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
#if MPI_PARALLEL_ENABLED
  //----- STEP 1: check that recv boundary buffer communications have all completed

  bool bflag = false;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (nghbr.h_view(m,n).gid >= 0) { // ID != -1, so not a physical boundary
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

  auto &mblev = pmy_pack->pmb->mb_lev;
  // Outer loop over (# of MeshBlocks)*(# of buffers)*(three field components)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/3;
    const int v = tmember.league_rank()%3;

    // scalar loop over neighbors to prevent race condition in overlapping assignments
    for (int n=0; n<nnghbr; ++n) {
      // only unpack buffers when neighbor exists
      if (nghbr.d_view(m,n).gid >= 0) {
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

  return TaskStatus::complete;
}
