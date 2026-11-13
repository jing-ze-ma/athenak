//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//! \brief functions to pack/send and recv/unpack boundary values for cell-centered (CC)
//! Mesh variables.
//! Prolongation of CC variables  occurs in ProlongateCC() function called from task list

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
#include "utils/wb_background.hpp"

//----------------------------------------------------------------------------------------
// BValCC constructor:

MeshBoundaryValuesCC::MeshBoundaryValuesCC(MeshBlockPack *pp, ParameterInput *pin,
                                           bool z4c) :
  MeshBoundaryValues(pp, pin, z4c) {
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::SetVectorPairs()
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
  for (int v=0; v<nvar; ++v) {h(v) = r[v];}
  Kokkos::deep_copy(vrole_, h);
  vrole_set_ = true;
}

//----------------------------------------------------------------------------------------
// speed2-1031 (bvsplit): the team decomposition of the CC pack / unpack.  They ran one
// team per (block, buffer, variable) with vector length 1: each thread walked the whole
// radial i range of its (k, j) column serially, and only the few teams of the populated
// buffers worked (HG cs, 3 blocks/rank: ~600 us per 1-variable pack on 168 teams).  Now
// each team is split into S teams over chunks of its (k, j) range and the i loop runs on
// VL vector lanes.  Every buffer / ghost value is the same expression of the same inputs
// (per-element copies / transforms): bitwise.  Env BVALS_CC_SPLIT (default 1 since
// defaults-1105, USER 11-05: 3.6 % faster than 16 on the 2-GPU dhj restart, bitwise; was 16)
// and BVALS_CC_VL (default 32 on CUDA, 16 elsewhere since fast2nd-1112); both = 1 is the
// old launch.  The env
// variable still overrides (BVALS_CC_SPLIT=16 = the 10-31 split).
static int BvCcSplit() {
  static int s = -1;
  if (s < 0) {
    const char *e = std::getenv("BVALS_CC_SPLIT");
    s = (e != nullptr) ? std::atoi(e) : 1;
    if (s < 1) {s = 1;}
  }
  return s;
}
// fast2nd-1112: BVALS_CC_VL is read on every backend; HIP default 16 (was fixed at 1).
// Small HG cs deck, P0 keys, 40 cycles, 1 MI300A (apudev 13066772, bitwise identical
// output): VL 1 27.5 / 28.1 s, VL 16 20.4 s, VL 64 20.6 / 22.0 s
// (SPLIT 4 + VL 64: 27.7 s)
static int BvCcVl() {
  static int v = -1;
  if (v < 0) {
    const char *e = std::getenv("BVALS_CC_VL");
#if defined(KOKKOS_ENABLE_CUDA)
    v = (e != nullptr) ? std::atoi(e) : 32;
#else
    v = (e != nullptr) ? std::atoi(e) : 16;
#endif
    if (v < 1) {v = 1;}
  }
  return v;
}

#if defined(KOKKOS_ENABLE_HIP)
// d6fast6-1113 (HIP only): the per-(k,j) seam transforms of the cubed-sphere pack can be
// kept in a per-THREAD slot of team scratch (LDS) instead of registers / the private
// segment.  The seam pack holds three SeamXform (stencil cells, an array indexed by the
// runtime stencil loop: 416 B of scratch per lane), the rho-guard SeamXform and a
// SeamFieldXform across the radial loop (228 VGPRs, occupancy 2 waves/SIMD on gfx942;
// rocprofv3, HG small cs deck: SendBuff 0.74 s, RecvBuff 0.19 s for the same grids).
// The transforms come from the same SeamXformAt / SeamFieldXformAt calls, are stored and
// re-read exactly (doubles through memory) and applied by the same ApplyMomentumXform /
// ApplyFieldXform: bitwise.  The loads are volatile so they are not hoisted back into
// registers.  Slot: xf[0..2] (3 x 14 Reals), xfc (14), fxf (15).
// env BVALS_CC_LDS (default 1); 0 = the register version (the launch before d6fast6).
constexpr int kBvXfN = 14, kBvXsSlot = 72;
KOKKOS_INLINE_FUNCTION
void BvXfPut(volatile Real *p, const cubed_sphere::SeamXform &t) {
  for (int c=0; c<3; ++c) {
    p[c] = t.s1[c]; p[3+c] = t.s2[c]; p[6+c] = t.d1[c]; p[9+c] = t.d2[c];
  }
  p[12] = t.cs; p[13] = t.det;
}
KOKKOS_INLINE_FUNCTION
void BvXfGet(const volatile Real *p, cubed_sphere::SeamXform &t) {
  for (int c=0; c<3; ++c) {
    t.s1[c] = p[c]; t.s2[c] = p[3+c]; t.d1[c] = p[6+c]; t.d2[c] = p[9+c];
  }
  t.cs = p[12]; t.det = p[13];
}
KOKKOS_INLINE_FUNCTION
void BvFxPut(volatile Real *p, const cubed_sphere::SeamFieldXform &t) {
  for (int c=0; c<3; ++c) {
    p[c] = t.s1[c]; p[3+c] = t.s2[c]; p[6+c] = t.d1[c]; p[9+c] = t.d2[c];
  }
  p[12] = t.sns; p[13] = t.csd; p[14] = t.snd;
}
KOKKOS_INLINE_FUNCTION
void BvFxGet(const volatile Real *p, cubed_sphere::SeamFieldXform &t) {
  for (int c=0; c<3; ++c) {
    t.s1[c] = p[c]; t.s2[c] = p[3+c]; t.d1[c] = p[6+c]; t.d2[c] = p[9+c];
  }
  t.sns = p[12]; t.csd = p[13]; t.snd = p[14];
}
static int BvCcLds() {
  static int l = -1;
  if (l < 0) {
    const char *e = std::getenv("BVALS_CC_LDS");
    l = (e != nullptr) ? std::atoi(e) : 1;
  }
  return l;
}
// the team size of the LDS seam pack (it sizes the scratch): env BVALS_CC_TS, default
// 256 / vector length = what Kokkos::AUTO gives there (rocprofv3: 16 x 16 workgroups).
// Which thread packs which (k, j) row does not change any value: bitwise.
static int BvCcTs() {
  static int t = -1;
  if (t < 0) {
    const char *e = std::getenv("BVALS_CC_TS");
    t = (e != nullptr) ? std::atoi(e) : 0;
    if (t <= 0) {t = 256/BvCcVl();}
    if (t < 1) {t = 1;}
  }
  return t;
}
#endif

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
  int no_rs_cc = 0;
  { const char *e_ = std::getenv("CS_NORESAMP_CC");
    if (e_ != nullptr) { no_rs_cc = std::atoi(e_); } }
  const Real rg_ = cs_rho_guard;
  const bool rs_lin = cs_lin_resample;
  const bool nors_ = cs_noresample;
  // wball-1104 wb_seam (see bvals.hpp)
  const bool wbs_ = wb_seam;
  auto wbsphi_ = wbs_phi;
  const EOS_Data wbseos_ = wbs_eos;
  const int wbsopt_ = wbs_opt, wbsien_ = wbs_ien;
  const bool wbsetg_ = wbs_etg;
  auto wbsbcc_ = wbs_bcc;
  const bool wbsmhd_ = (wbs_bcc.extent_int(0) > 0);
  const bool ext_ = cs_ext_pass_;   // layoutbug-1103: second pass, see cs_seam_blockend
  const bool posg_ = pmy_pack->pmesh->cs_seam_pos_global;   // layoutfix-1106
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
  const int nsp_ = BvCcSplit();
  BvCcTeamPolicy policy(DevExeSpace(), nmnv*nsp_, Kokkos::AUTO, BvCcVl());
  Kokkos::TeamPolicy<> policy1(DevExeSpace(), nmnv, Kokkos::AUTO);   // (z4c kernel)
  // The seam/pole branch below keeps per-thread arrays (stencil indices, three SeamXform)
  // that give the kernel a PRIVATE SEGMENT (scratch).  A dispatch that needs scratch
  // right after ones that do not stalls the GPU queue ~140 us while the runtime
  // (re)attaches scratch -- measured before every pack of the implicit M1 halo exchange
  // (tests_m1/runs_3k_gpu3d/README_HALO.md).  So the body is compiled twice: with the
  // branch (cubed sphere / polar) and without it (everything else), where do_cs and
  // do_pole are false anyway, i.e. the same arithmetic; still one launch per call.
  // nvcc rejects generic (auto) extended lambdas: under CUDA the flag is a plain bool
  // (and the call sites pass true/false), the same arithmetic.
  // d6fast3-1112: a second compile-time flag kwb_ for the wb_seam branch (WBAdvance with
  // the EOS): with wb_seam off (every cs deck but the dhj wb runs) it is dead code, but
  // it took the seam pack from 816 to 4816 bytes of scratch per lane (HIP, 128 VGPRs).
  // The seam variant without it is launched when wb_seam is off: wbsv is false there
  // anyway, so the arithmetic is the same (bitwise).
#if defined(KOKKOS_ENABLE_CUDA)
  // d6fast4-1112: nvcc has no generic extended lambda, and the bool flags did not compile
  // the dead branches out (the CUDA pack had 255 registers + 3 kB stack): three kernels
  // with kseam_ / kwb_ compile-time constants (bvals_cc_pack.inc), the same arithmetic
#define BVCC_PACK_NAME pack_tt_
#define BVCC_PACK_KSEAM true
#define BVCC_PACK_KWB true
#include "bvals_cc_pack.inc"
#undef BVCC_PACK_KWB
#undef BVCC_PACK_NAME
#define BVCC_PACK_NAME pack_tf_
#define BVCC_PACK_KWB false
#include "bvals_cc_pack.inc"
#undef BVCC_PACK_KWB
#undef BVCC_PACK_KSEAM
#undef BVCC_PACK_NAME
#define BVCC_PACK_NAME pack_ff_
#define BVCC_PACK_KSEAM false
#define BVCC_PACK_KWB false
#include "bvals_cc_pack.inc"
#undef BVCC_PACK_KWB
#undef BVCC_PACK_KSEAM
#undef BVCC_PACK_NAME
#else
#define BVCC_PACK_NAME pack_
#include "bvals_cc_pack.inc"
#undef BVCC_PACK_NAME
#endif
#if defined(KOKKOS_ENABLE_HIP)
  // d6fast6-1113: the seam variants with the transform slots (explicit team size)
  if ((use_cs || use_pole) && BvCcLds() != 0) {
    BvCcTeamPolicy policy_s(DevExeSpace(), nmnv*nsp_, BvCcTs(), BvCcVl());
    policy_s.set_scratch_size(0, Kokkos::PerTeam(sizeof(Real)*kBvXsSlot*BvCcTs()));
    if (wbs_) {
      BvalsTeamFor("SendBuff", policy_s, KOKKOS_LAMBDA(TeamMember_t tmember) {
        pack_(tmember, std::true_type(), std::true_type(), std::true_type());
      });
    } else {
      BvalsTeamFor("SendBuff", policy_s, KOKKOS_LAMBDA(TeamMember_t tmember) {
        pack_(tmember, std::true_type(), std::false_type(), std::true_type());
      });
    }
  } else if ((use_cs || use_pole) && wbs_) {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
      pack_(tmember, std::true_type(), std::true_type(), std::false_type());
    });
  } else if (use_cs || use_pole) {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
      pack_(tmember, std::true_type(), std::false_type(), std::false_type());
    });
  } else {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
      pack_(tmember, std::false_type(), std::false_type(), std::false_type());
    });
  }  // end par_for_outer
#else
  if ((use_cs || use_pole) && wbs_) {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
#if defined(KOKKOS_ENABLE_CUDA)
      pack_tt_(tmember);
#else
      pack_(tmember, std::true_type(), std::true_type());
#endif
    });
  } else if (use_cs || use_pole) {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
#if defined(KOKKOS_ENABLE_CUDA)
      pack_tf_(tmember);
#else
      pack_(tmember, std::true_type(), std::false_type());
#endif
    });
  } else {
    BvalsTeamFor("SendBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
#if defined(KOKKOS_ENABLE_CUDA)
      pack_ff_(tmember);
#else
      pack_(tmember, std::false_type(), std::false_type());
#endif
    });
  }  // end par_for_outer
#endif

  // The whole body of this kernel sits under (is_z4c && multilevel), both of which are
  // host-side constants, so when they do not hold the launch has nothing to do at all:
  // skip it.  It is otherwise a full-size empty team launch on EVERY halo exchange.
  if (is_z4c_ && ml_) {
  BvalsTeamFor("SendBuffZ4c", policy1, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int vl = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);
    const int v = vlst_.v[vl];

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

  // the variables moved (default vmask = 0: all of a, in order); see bvals.hpp
  const BvVarList vlst_ = MakeBvVarList(a.extent_int(1), vmask);
  int nvar = (vmask != 0) ? vlst_.n : a.extent_int(1);

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables)
  const int nsp_ = BvCcSplit();
  BvCcTeamPolicy policy(DevExeSpace(), (nmb*nnghbr*nvar)*nsp_, Kokkos::AUTO, BvCcVl());
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

      // Middle loop over k,j
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, q0_, q1_), [&](const int idx) {
        int k = idx / nj;
        int j = (idx - k * nj) + jl;
        k += kl;

        // CUBE-VERTEX CORNER FILL: OWNERSHIP.  With `<mesh>/cs_vertex_fill_cc` the x2/x3
        // face and x1x2/x3x1 edge buffers were widened by ng along the seam (see
        // WidenForCubeVertex in buffs_cc.cpp) so the doubly-ghost corner block rides the
        // face exchange.  Both flanking buffers reach that block, but each panel only
        // genuinely COVERS the half of it on its own side of the diagonal, so write a
        // doubly-ghost cell only from the neighbour whose panel contains that direction.
        // Decided here, from the RECEIVER's own geometry, so it needs nothing from the
        // sender and is identical on every rank.  ONLY the slots that were widened:
        // 8-23 (x2 faces and x1x2 edges) and 24-39 (x3 faces and x3x1 edges).  The test
        // MUST NOT reach the x2x3 edge (40-47) or corner (48-55) buffers, whose
        // destinations are legitimately doubly-ghost -- gating them here would suppress
        // the real edge and corner exchange.  Cells outside the corner block, and any
        // grid that is not a cubed sphere, never enter this.
        if (cs_vfill_cc_ && (n >= 8) && (n < 40) &&
            (nghbr.d_view(m,n).lev == mblev.d_view(m)) &&
            ((j < js_) || (j > je_)) && ((k < ks_) || (k > ke_))) {
          // ONLY a CUBE VERTEX needs the fill.  At an ordinary block corner the diagonal
          // neighbour exists and its own buffer (40..55, unpacked LATER in this same
          // scalar n loop) is the correct filler, and the panel-edge geometry below is
          // not even valid there.
          const int sj = (j < js_) ? 0 : 1;
          const int sk = (k < ks_) ? 0 : 1;
          if (!IsCubeVertexCorner(nghbr.d_view, mbpanel.d_view, m,
                                  40 + 2*(sj + 2*sk))) return;
          // OWNERSHIP, exactly.  The corner ghost block is covered by the two flanking
          // panels, split by the seam BETWEEN THEM, and in this panel's own gnomonic
          // chart that seam is exactly the DIAGONAL.  Take P = +z: tan(xi) = x/z and
          // tan(eta) = y/z, and the +x/+y panel boundary is the plane x = y, i.e.
          // xi = eta.  With a = s2*xi - pi/4 and b = s3*eta - pi/4 measuring how far the
          // cell lies beyond each panel edge, a > b means the direction belongs to the
          // panel across the x2 face and a < b to the one across the x3 face; the other
          // seven vertices follow by reflection, which is what the signs are.  The test
          // is COMPLEMENTARY by construction, so one of the two buffers writes every
          // cell of the block and neither writes a cell the other owns -- which matters
          // because FillPanelCornersCC stops extrapolating the block once the fill is on.
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

  // Outer loop over (# of MeshBlocks)*(# of buffers)*(# of variables).  Whole body is
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
    if (pmy_pack->pmesh->multilevel) { FillPanelCornersCC(ca, true, vmask); }
  }

  // layoutbug-1103: <mesh>/cs_seam_blockend = extend.  The along-seam resample of the
  // pack above was clamped to the SENDING block's active range, so at an interior block
  // end on a seam (more than one MeshBlock per panel side) the ghost was EXTRAPOLATED
  // from that block alone while a single block per panel interpolates: a MeshBlock-
  // layout dependence of the seam halo at truncation level.  Every ghost is now in place,
  // in particular each sender's along-seam ghosts next to such an end (same-panel copies),
  // so the whole exchange is repeated once with the stencil allowed into them (blocking;
  // same buffers, same variables and window).  Non-seam buffers and every stencil that
  // did not reach a block end come out bitwise as in the first pass.
  if (cs_blockend_ext_ && !cs_ext_pass_ && !cs_noresample &&
      pmy_pack->pmesh->use_cubed_sphere) {
    cs_ext_pass_ = true;
    (void) ClearSend();
    (void) InitRecv(a.extent_int(1), vmask);
    while (PackAndSendCC(a, ca, iwl, iwu, vmask) != TaskStatus::complete) {}
    while (RecvAndUnpackCC(a, ca, iwl, iwu, vmask) != TaskStatus::complete) {}
    // seamexact-1103: wait on this pass's own sends here.  A caller that cleared its
    // sends before unpacking (rad_m1_coupling.cpp FluidRefresh) would otherwise leave
    // them in flight into its next exchange, which reuses the buffers and the requests.
    (void) ClearSend();
    cs_ext_pass_ = false;
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

namespace cs_vertex_diag { extern Real cc_perturb, perturb_add; }   // vertexacc-1103 (mesh.cpp)

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
  const Real vpert_ = cs_vertex_diag::cc_perturb;   // vertexacc-1103 (mesh.cpp)
  const Real vadd_ = cs_vertex_diag::perturb_add;   // cornerfix-1103 (mesh.cpp)
  const bool vpon_ = (vpert_ != 0.0) || (vadd_ != 0.0);
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
    // THE SAMPLED FILL SUPERSEDES THIS, exactly as in FillPanelCornersFC.  When the two
    // flanking neighbours are at the SAME level their widened buffers have already
    // written this corner block with REAL data from the two panels that cover it, which
    // is strictly better than any extrapolation, so do not overwrite it.  At a level
    // boundary only `isame` was widened, so the extrapolation below is still the best
    // available, and the same holds on the COARSE array, which is never widened.
    //
    // WHICH BUFFER REACHES WHICH RADIAL LAYER.  Over the ACTIVE radial range the two
    // widened FACE buffers cover the block.  A RADIAL GHOST layer is triply ghost, so
    // only a buffer already ghost in x1 can reach it: the x1x2 and x3x1 EDGE buffers,
    // widened by the same rule.  BOTH must exist at this level, because ownership splits
    // the block along the diagonal and one alone would leave half of it written by
    // nothing.
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
        if (vpon_) {   // cornerfix-1103 diagnostic: the sampled block is perturbed too
          const int jt = (sj < 0) ? (js-1-gj) : (je+1+gj);
          const int kt = (sk < 0) ? (ks-1-gk) : (ke+1+gk);
          a_(m,v,kt,jt,i) = a_(m,v,kt,jt,i)*(1.0 + vpert_) + vadd_;
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
    // SIGN-PRESERVING FALLBACK -- see the note above the function.  It only fires when
    // the three nodes share a strict sign and the extrapolant flips it (or zeroes it),
    // and then drops to the node next to the ghost; smooth data is kept exactly.
    ek = SignPreserveCS(ek, k0, k1, k2);
    ej = SignPreserveCS(ej, j0, j1, j2);
    a_(m,v,kt,jt,i) = 0.5*(ek + ej)*(1.0 + vpert_);   // vpert_ = 0: bitwise
    if (vadd_ != 0.0) a_(m,v,kt,jt,i) += vadd_;
  });
  return;
}

