//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for cell-centered (CC)
//! Mesh variables.
//! Prolongation of CC variables occurs in ProlongateCC function called from task list

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals.hpp"
#include "bvals/cs_seam_round.hpp"

//----------------------------------------------------------------------------------------
// BValCC constructor:

MeshBoundaryValuesCC::MeshBoundaryValuesCC(MeshBlockPack *pp, ParameterInput *pin,
                                           bool z4c) :
  MeshBoundaryValues(pp, pin, z4c) {
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::SetVectorPairs
//! \brief Replace the default (IVY, IVZ) tangential-pair rule of this exchange by an
//! explicit table: each {a,b} is a pair (a = x2 member, b = x3 member), every other slot
//! of the nvar exchanged is a scalar (no polar sign flip, no seam basis transform).

void MeshBoundaryValuesCC::SetVectorPairs(const int nvar,
                                          const std::vector<std::pair<int,int>> &pairs) {
  std::vector<int> r(nvar, -1);
  for (const auto &p : pairs) {
    const int a = p.first, b = p.second;
    if (a < 0 || b < 0 || a >= nvar || b >= nvar || a == b || r[a] != -1 || r[b] != -1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "SetVectorPairs: bad pair (" << a << "," << b
                << ") for nvar = " << nvar << std::endl;
      std::exit(EXIT_FAILURE);
    }
    r[a] = 2*b;
    r[b] = 2*a + 1;
  }
  Kokkos::realloc(vrole_, std::max(nvar, 1));
  auto h = Kokkos::create_mirror_view(vrole_);
  for (int v=0; v<nvar; ++v) {
    h(v) = r[v];
  }
  Kokkos::deep_copy(vrole_, h);
  vrole_set_ = true;
}

//----------------------------------------------------------------------------------------
//! \fn static int BvCcSplit
//! \brief Number of teams each (block, buffer, variable) of the CC pack/unpack is split
//! into over its (k,j) range (env BVALS_CC_SPLIT, else the object's cc_split).
static int BvCcSplit(const int def) {
  static int s = -1;
  if (s < 0) {
    const char *e = std::getenv("BVALS_CC_SPLIT");
    s = (e != nullptr) ? std::max(std::atoi(e), 1) : 0;
  }
  return (s > 0) ? s : std::max(def, 1);
}
//----------------------------------------------------------------------------------------
//! \fn static int BvCcVl
//! \brief Vector length of the CC pack/unpack (env BVALS_CC_VL; default 32 on CUDA, 16
//! otherwise).
static int BvCcVl() {
  static int v = -1;
  if (v < 0) {
    const char *e = std::getenv("BVALS_CC_VL");
#if defined(KOKKOS_ENABLE_CUDA)
    v = (e != nullptr) ? std::atoi(e) : 32;
#else
    v = (e != nullptr) ? std::atoi(e) : 16;
#endif
    if (v < 1) {
      v = 1;
    }
  }
  return v;
}

// force-inline the seam-pack stencil lambdas on CUDA (avoids local-memory spills)
#if defined(KOKKOS_ENABLE_CUDA)
#define BV_CUDA_INLINE __attribute__((always_inline))
#else
#define BV_CUDA_INLINE
#endif

//----------------------------------------------------------------------------------------
//! \fn static BvCcTeamPolicy BvCcPolicy
//! \brief Team policy of the CC pack/unpack (env BVALS_CC_TS = team size, 0 = AUTO).
static BvCcTeamPolicy BvCcPolicy(const int nleague) {
  static int t = -1;
  if (t < 0) {
    const char *e = std::getenv("BVALS_CC_TS");
    t = (e != nullptr) ? std::atoi(e) : 0;
    if (t < 0) {
      t = 0;
    }
  }
  if (t > 0) {
    return BvCcTeamPolicy(DevExeSpace(), nleague, t, BvCcVl());
  }
  return BvCcTeamPolicy(DevExeSpace(), nleague, Kokkos::AUTO, BvCcVl());
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::PackAndSendCC
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
                                               int iwl, int iwu, int vmask) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  // the variables moved (default vmask = 0: all of a, in order); see bvals.hpp
  if (vmask != 0 && (pmy_pack->pmesh->multilevel || a.extent_int(1) > 32)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "PackAndSendCC: variable list needs a uniform grid and "
              << "<= 32 variables" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const BvVarList vlst_ = MakeBvVarList(a.extent_int(1), vmask);
  int nvar = (vmask != 0) ? vlst_.n : a.extent_int(1);

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
  const Real rg_ = cs_rho_guard;
  const bool rs_lin = cs_lin_resample;
  const bool nors_ = cs_noresample;
  const bool ext_ = cs_ext_pass_;   // second pass, see cs_seam_blockend
  const int ngx_ = pmy_pack->pmesh->mb_indcs.ng;
  const bool cperm_ = cs_perm_pairs;
  const int hrq_ = cs_hrq0;
  auto sbuf = SendBufDv();
  auto rbuf = RecvBufDv();
  auto &is_z4c = is_z4c_;
  auto &multilevel = pmy_pack->pmesh->multilevel;
  // the component-role table (SetVectorPairs); unset = the (IVY, IVZ) rule
  const bool vtab_ = vrole_set_;
  auto vrole = vrole_;
  // x1 index window (see the declaration): only buffers spanning the whole x1 extent
  const int wlo_ = iwl, whi_ = iwu;
  const int mbis_ = pmy_pack->pmesh->mb_indcs.is;
  const int mbie_ = pmy_pack->pmesh->mb_indcs.ie;
  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables)
  int nmnv = nmb*nnghbr*nvar;
  const int vmask_ = vmask;   // (BvVarOf, CUDA)
  const int nsp_ = BvCcSplit(cc_split);
  BvCcTeamPolicy policy = BvCcPolicy(nmnv*nsp_);
  Kokkos::TeamPolicy<> policy1(DevExeSpace(), nmnv, Kokkos::AUTO);   // (z4c kernel)
  // the seam/pole branch needs per-thread scratch; the body is compiled twice (with and
  // without it, selected at launch) so ordinary packs do not pay the scratch attach
#if defined(KOKKOS_ENABLE_CUDA)
  // nvcc has no generic extended lambda, and a bool flag did not compile the dead
  // branches out: two kernels with kseam_ a compile-time constant (bvals_cc_pack.inc),
  // the same arithmetic
#define BVCC_PACK_NAME pack_t_
#define BVCC_PACK_KSEAM true
#include "bvals_cc_pack.inc"
#undef BVCC_PACK_KSEAM
#undef BVCC_PACK_NAME
#define BVCC_PACK_NAME pack_f_
#define BVCC_PACK_KSEAM false
#include "bvals_cc_pack.inc"  // NOLINT(build/include)
#undef BVCC_PACK_KSEAM
#undef BVCC_PACK_NAME
#else
#define BVCC_PACK_NAME pack_
#include "bvals_cc_pack.inc"  // NOLINT(build/include)
#undef BVCC_PACK_NAME
#endif
  if (use_cs || use_pole) {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
#if defined(KOKKOS_ENABLE_CUDA)
      pack_t_(tmember);
#else
      pack_(tmember, std::true_type());
#endif
    });
  } else {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
#if defined(KOKKOS_ENABLE_CUDA)
      pack_f_(tmember);
#else
      pack_(tmember, std::false_type());
#endif
    });
  }  // end par_for_outer

  // the kernel body only acts for z4c with multilevel: skip the launch otherwise
  if (is_z4c_ && ml_) {
  BvalsTeamFor("SendBuffZ4c", policy1, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int vl = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);
    const int v = vlst_.v[vl];

    // only load buffers when neighbor exists. CUBED-SPHERE CUBE VERTEX: skipped on
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
              rbuf[dn].vars(dm,ndat+ (i-il + ni*(j-jl + nj*(k-kl + nk*vl)))) =
                ca(m,v,k,j,i);
            });

          // else copy into send buffer for MPI communication below
          } else {
            // load data from coarse_u0
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,il,iu+1),
            [&](const int i) {
              sbuf[n].vars(m,ndat+ (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) )=ca(m,v,k,j,i);
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
  // fence only before the first MPI_Isend: an exchange with no off-rank neighbour
  // (1 rank) has nothing to wait for -- the unpack kernel is on the same stream.
  bool fenced_ = false;
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
          if (!fenced_) {
            Kokkos::fence();
            fenced_ = true;
          }
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int lid = nghbr.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);

          // get ptr to send buffer when neighbor is at coarser/same/fine level
          int data_size = nvar;
          int per_var = 0;
          if ( nghbr.h_view(m,n).lev < pmy_pack->pmb->mb_lev.h_view(m) ) {
            per_var = sendbuf[n].icoar_ndat;
          } else if ( nghbr.h_view(m,n).lev == pmy_pack->pmb->mb_lev.h_view(m) ) {
            if (is_z4c) {
              per_var = sendbuf[n].isame_z4c_ndat;
            } else {
              per_var = sendbuf[n].isame_ndat;
            }
          } else {
            per_var = sendbuf[n].ifine_ndat;
          }
          data_size *= per_var;
          auto send_ptr = Kokkos::subview(sendbuf[n].vars, m, Kokkos::ALL);

          // variable list: the selected variables are the first nvar chunks
          int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL,
                               drank, tag, comm_vars, &(sendbuf[n].vars_req[m]));
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
// \!fn void RecvBuffers
// \brief Unpack boundary buffers

TaskStatus MeshBoundaryValuesCC::RecvAndUnpackCC(DvceArray5D<Real> &a,
                                                 DvceArray5D<Real> &ca,
                                                 int iwl, int iwu, int vmask) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto rbuf = RecvBufDv();
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
  // For the cell-centred cube-vertex corner fill's ownership test (see its use below).
  const bool cs_vfill_cc_ = pmy_pack->pmesh->cs_vertex_fill_cc;
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
          !(use_cs && IsCubeVertexCorner(nghbr.h_view, mbpanel.h_view, m, n)) &&
          !(skipd_ && IsX2X3DiagSlot(n))) {
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

  //----- STEP 2: buffers have all completed, so unpack

  // the variables moved (default vmask = 0: all of a, in order); see bvals.hpp
  const BvVarList vlst_ = MakeBvVarList(a.extent_int(1), vmask);
  int nvar = (vmask != 0) ? vlst_.n : a.extent_int(1);

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables)
  const int nsp_ = BvCcSplit(cc_split);
  BvCcTeamPolicy policy = BvCcPolicy((nmb*nnghbr*nvar)*nsp_);
  Kokkos::TeamPolicy<> policy1(DevExeSpace(), (nmb*nnghbr*nvar), Kokkos::AUTO);   // (z4c)
  BvalsTeamFor("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int lr_ = tmember.league_rank()/nsp_;
    const int sp_ = tmember.league_rank() - lr_*nsp_;
    const int m = (lr_)/(nnghbr*nvar);
    const int n = (lr_ - m*(nnghbr*nvar))/nvar;
    const int vl = (lr_ - m*(nnghbr*nvar) - n*nvar);
    const int v = vlst_.v[vl];

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
      // bvsplit: this team's chunk of the (k, j) range
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

      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, q0_, q1_),
                           [&](const int idx) {
        int k = idx / nj;
        int j = (idx - k * nj) + jl;
        k += kl;

        // cube-vertex corner fill (<mesh>/cs_vertex_fill_cc): the widened x2/x3 face and
        // edge buffers both reach the corner block; each panel writes the half on its own
        // side of the diagonal
        if (cs_vfill_cc_ && (n >= 8) && (n < 40) &&
            (nghbr.d_view(m,n).lev == mblev.d_view(m)) &&
            ((j < js_) || (j > je_)) && ((k < ks_) || (k > ke_))) {
          // only a cube vertex needs the fill (an ordinary corner has a diagonal
          // neighbour whose buffer is unpacked later)
          const int sj = (j < js_) ? 0 : 1;
          const int sk = (k < ks_) ? 0 : 1;
          if (!IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m,
                                  40 + 2*(sj + 2*sk))) return;
          // ownership: in this panel's gnomonic chart the seam between the two flanking
          // panels is the diagonal of the corner block; each panel writes its own side
          const Real xi  = 0.25*M_PI*CellCenterX(j-js_, nx2_, mbsz.d_view(m).x2min,
                                                              mbsz.d_view(m).x2max);
          const Real eta = 0.25*M_PI*CellCenterX(k-ks_, nx3_, mbsz.d_view(m).x3min,
                                                              mbsz.d_view(m).x3max);
          const Real aa = (sj ? xi : -xi) - 0.25*M_PI;
          const Real bb = (sk ? eta : -eta) - 0.25*M_PI;
          if ((n < 24) != (aa >= bb)) return;
        }

        // if neighbor is at same or finer level, load data directly into u0
        if (nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
          [&](const int i) {
            a(m,v,k,j,i) = rbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) );
          });

        // if neighbor is at coarser level, load data into coarse_u0
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(tmember,ilp,iup+1),
          [&](const int i) {
            ca(m,v,k,j,i) = rbuf[n].vars(m, (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) );
          });
        }
      });
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables). Whole body is
  // under the host-side (is_z4c && multilevel): skip the launch when it cannot fire.
  if (is_z4c_ && ml_) {
  BvalsTeamFor("RecvBuffZ4c", policy1, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int vl = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);
    const int v = vlst_.v[vl];
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
            ca(m,v,k,j,i) =
              rbuf[n].vars(m,ndat + (i-il + ni*(j-jl + nj*(k-kl + nk*vl))) );
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
    FillPanelCornersCC(a, false, vmask);
    if (pmy_pack->pmesh->multilevel) {
      FillPanelCornersCC(ca, true, vmask);
    }
  }

  // <mesh>/cs_seam_blockend = extend: a second pass re-packs the seam halo with the
  // ghosts of the sending block filled, so at an internal block end the along-seam
  // resample interpolates instead of extrapolating (layout-independent seam halo)
  if (cs_blockend_ext_ && !cs_ext_pass_ && !cs_noresample &&
      pmy_pack->pmesh->use_cubed_sphere) {
    cs_ext_pass_ = true;
    (void) ClearSend();
    (void) InitRecv(a.extent_int(1), vmask);
    while (PackAndSendCC(a, ca, iwl, iwu, vmask) != TaskStatus::complete) {}
    while (RecvAndUnpackCC(a, ca, iwl, iwu, vmask) != TaskStatus::complete) {}
    // wait on this pass's own sends (callers may reuse the buffers and requests)
    (void) ClearSend();
    cs_ext_pass_ = false;
  }

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::FillPanelCornersCC
//! \brief Fills the ng x ng corner ghost cells at a cube vertex, where only three panels
//! meet and the generic corner buffer has no meaningful source. Each corner ghost is
//! extrapolated quadratically from each of the two flanking face halos (in the panel's
//! own gnomonic chart) and the two are averaged, with a sign-preserving fallback.
//! Same-panel corners are not touched.

//----------------------------------------------------------------------------------------
//! \fn Real SignPreserveCS
//! \brief Sign-preserving fallback for a one-sided quadratic cube-vertex extrapolant q
//! built from the nodes (n0,n1,n2), n0 adjacent to the ghost: if the nodes share a strict
//! sign and q has the opposite sign or is zero, n0 is returned (constant extrapolation);
//! otherwise q.

namespace cs_vertex_diag { extern Real cc_perturb; }   // mesh.cpp

KOKKOS_INLINE_FUNCTION
static Real SignPreserveCS(const Real q, const Real n0, const Real n1, const Real n2) {
  if (((n0 > 0.0) && (n1 > 0.0) && (n2 > 0.0) && !(q > 0.0)) ||
      ((n0 < 0.0) && (n1 < 0.0) && (n2 < 0.0) && !(q < 0.0))) {
    return n0;
  }
  return q;
}

void MeshBoundaryValuesCC::FillPanelCornersCC(DvceArray5D<Real> &a, bool coarse,
                                              int vmask) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int js = coarse ? indcs.cjs : indcs.js;
  const int je = coarse ? indcs.cje : indcs.je;
  const int ks = coarse ? indcs.cks : indcs.ks;
  const int ke = coarse ? indcs.cke : indcs.ke;
  const int ng = indcs.ng;
  const int nmb = pmy_pack->nmb_thispack;
  // the variables filled (default vmask = 0: all of a); see bvals.hpp
  const BvVarList vlst_ = MakeBvVarList(a.extent_int(1), vmask);
  const int nvar = (vmask != 0) ? vlst_.n : a.extent_int(1);
  const int n1 = a.extent_int(4);
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbpanel = pmy_pack->pmb->mb_panel;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto a_ = a;
  const bool vfill_ = pmy_pack->pmesh->cs_vertex_fill_cc;
  const bool coar_ = coarse;
  const Real vpert_ = cs_vertex_diag::cc_perturb;   // (mesh.cpp)
  const bool vpon_ = (vpert_ != 0.0);
  const int is_ = indcs.is, ie_ = indcs.ie;

  // par_for takes at most five ranges, so the ng x ng corner block is flattened into g
  par_for("cs_fill_corners_cc", DevExeSpace(), 0,(nmb-1), 0,3, 0,(nvar-1), 0,(n1-1),
          0,(ng*ng-1),
  KOKKOS_LAMBDA(const int m, const int c, const int vl, const int i, const int g) {
    const int v = vlst_.v[vl];
    const int gj = g/ng;
    const int gk = g - gj*ng;
    const int sj = (c & 1) ? 1 : -1;      // -1 = the -x2 side of the block
    const int sk = (c & 2) ? 1 : -1;      // -1 = the -x3 side
    // Buffer ids of the two FACE neighbours flanking this corner (see nghbr_index.hpp).
    const int nj_id = (sj < 0) ? 8 : 12;
    const int nk_id = (sk < 0) ? 24 : 28;
    const int mp = mbpanel.d_view(m);
    // only a true cube vertex: both flanking faces are panel seams
    bool seamj = false, seamk = false;
    if (nghbr.d_view(m,nj_id).gid >= 0 && nghbr.d_view(m,nj_id).panel != mp) seamj = true;
    if (nghbr.d_view(m,nk_id).gid >= 0 && nghbr.d_view(m,nk_id).panel != mp) seamk = true;
    if (!(seamj && seamk)) return;
    // the sampled fill (same-level flanking neighbours with widened buffers) already
    // wrote this block with real data; extrapolate only otherwise
    const bool vfill_here = (vfill_ && !coar_ &&
        nghbr.d_view(m,nj_id).lev == mblev.d_view(m) &&
        nghbr.d_view(m,nk_id).lev == mblev.d_view(m));
    if (vfill_here) {
      const int s1 = (i < is_) ? -1 : 1;              // radial side of a ghost layer
      const int e12 = 16 + (s1 + 1) + 2*(sj + 1);     // x1x2 edge (nghbr_index.hpp)
      const int e31 = 24 + (s1 + 9) + 2*(sk + 1);     // x3x1 edge
      const bool vfill_rg =
          nghbr.d_view(m,e12).gid >= 0 && nghbr.d_view(m,e12).lev == mblev.d_view(m) &&
          nghbr.d_view(m,e31).gid >= 0 && nghbr.d_view(m,e31).lev == mblev.d_view(m);
      if ((i >= is_ && i <= ie_) || vfill_rg) {
        if (vpon_) {   // diagnostic: the sampled block is perturbed too
          const int jt = (sj < 0) ? (js-1-gj) : (je+1+gj);
          const int kt = (sk < 0) ? (ks-1-gk) : (ke+1+gk);
          a_(m,v,kt,jt,i) = a_(m,v,kt,jt,i)*(1.0 + vpert_);
        }
        return;
      }
    }
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
    // sign-preserving fallback (see SignPreserveCS)
    ek = SignPreserveCS(ek, k0, k1, k2);
    ej = SignPreserveCS(ej, j0, j1, j2);
    a_(m,v,kt,jt,i) = 0.5*(ek + ej)*(1.0 + vpert_);   // vpert_ = 0: bitwise
  });
  return;
}

