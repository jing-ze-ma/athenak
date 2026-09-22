//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_tau.cpp
//! \brief <rad_m1>/closure = tau: the Eddington tensor of the multi-D implicit solve from
//! the COLUMN OPTICAL DEPTH, never from the cell's own flux.
//!
//!   tau(x1) = int_{x1}^{x1max} rho (kappa_F + kappa_s) dx1'   (from the domain top)
//!   chi     = (tau + q_inf) / (3 (tau + q(tau))),  q(tau) = 0.710446 - 0.133054
//!             exp(-3.4488 tau)  (exact grey plane-parallel K/J, Hopf-function fit)
//!   n       = -grad tau / |grad tau|,  D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b
//!
//! It is the DIAGNOSTIC <rad_m1>/dbg_tensor = tau of milestone 3i made a real option with
//! MeshBlock stacks along x1 over MPI ranks.  On ONE MeshBlock along x1 the arithmetic is
//! the diagnostic's operation for operation (gate: bitwise identical).
//!
//! PARTITION.  Each MeshBlock integrates its own part of every column (its interior
//! cells is..ie, for its interior columns AND the x2/x3 face-ghost columns) and adds the
//! column integrals S_q of the blocks q above it in the same x1 stack, summed from the
//! top down in a fixed order.  The S_q are exchanged directly between the members of a
//! stack (one message of (nx3+2ng)(nx2+2ng) reals per ordered pair of members on
//! different ranks, per hydro step; same-rank members are read in place on the device).
//!
//! Why not the gather of implicit_partition = gather (ImplicitGatherSolve): that moves
//! 4 nx1 reals per column per block to a root and back each Picard pass because a
//! tridiagonal solve is a global recurrence.  The optical depth is a plain prefix sum,
//! and a prefix sum decomposes exactly into local partial sums plus an offset: one
//! scalar per column per block is all that has to move, once per step, and every block
//! can form its own offset without a root and a scatter.  Only the stack TOPOLOGY idea
//! (members = same (lx2,lx3,level), ordered by lx1) is reused; it is rebuilt here so that
//! this file does not depend on the gather's internal tables.
//!
//! GHOST COLUMNS.  -grad tau across an x2/x3 block face needs tau in the ghost column.
//! It is computed HERE from the ghost values of M1_IW_KT, which the transverse halo
//! (ImplicitTransverseHalo, run before the Picard loop) fills with the very numbers the
//! neighbour holds, and from the column integrals of the ghost columns of the blocks
//! above, which on a uniform mesh are the neighbour stack's own integrals over the same
//! x1 range.  The ghost tau is therefore bit-identical to the tau the neighbour computes
//! for that column itself, without another exchange.  A PHYSICAL (non-periodic) x2/x3
//! face uses a one-sided difference.
//!
//! ROUND-OFF.  Across x2/x3 decompositions the result is bitwise invariant (above).
//! Along x1 the association of the sum changes (block partial sums + offset instead of
//! one running sum), so 1 vs N blocks along x1 agree to round-off only.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

namespace {
void TauFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_m1>/closure = tau: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
// MPI tags of the column-integral exchange: base + the SENDER's local MeshBlock id
// (< 2^14, see NUM_BITS_LID), which keeps every tag <= 32767, the MPI minimum of
// MPI_TAG_UB.  Every message of this exchange is completed before TauClosureBuild
// returns, so nothing of it is in flight together with the gather solve's messages.
constexpr int TAU_TAG_BASE = 16384;
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::TauClosureInit
//! \brief check the mesh, build the x1 stack topology, allocate.  Called on the first
//! TauClosureBuild (by every rank: all ranks run the implicit solve together).

void RadiationM1::TauClosureInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto &mindcs = pm->mesh_indcs;
  if (!trans_on) {
    TauFatal("needs <rad_m1>/transport = implicit on a MULTI-D mesh (the 1-D branch "
             "of the solve never reads the tensor)");
  }
  if (pm->multilevel) {
    TauFatal("SMR/AMR is not supported (uniform mesh only)");
  }
  if ((mindcs.nx1 % indcs.nx1) != 0) {
    TauFatal("<mesh>/nx1 must be an exact multiple of <meshblock>/nx1");
  }
  if (ibc_x1min == M1_IBC_PERIODIC) {
    TauFatal("tau is measured from the top of the domain (x1max); a PERIODIC x1 "
             "boundary has no top");
  }
  tau_nblk = mindcs.nx1/indcs.nx1;
  const int nblk = tau_nblk;
  const int nmb = pmy_pack->nmb_thispack;
  const int g0 = pmy_pack->gids;
  const int me = global_variable::my_rank;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;

  // the stack of every local block: members = same (lx2, lx3, level), position = lx1
  Kokkos::realloc(tau_pos, nmb);
  Kokkos::realloc(tau_loc, nmb*nblk);
  tau_mrank.assign(nmb*nblk, -1);
  tau_mlid.assign(nmb*nblk, -1);
  tau_any_mpi = false;
  for (int m=0; m<nmb; ++m) {
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    tau_pos.h_view(m) = static_cast<int>(lm.lx1);
    for (int q=0; q<nblk; ++q) {tau_loc.h_view(m*nblk+q) = -1;}
    for (int g=0; g<pm->nmb_total; ++g) {
      LogicalLocation &lg = pm->lloc_eachmb[g];
      if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level &&
          lg.lx1 >= 0 && lg.lx1 < nblk) {
        int q = static_cast<int>(lg.lx1);
        int rk = pm->rank_eachmb[g];
        tau_mrank[m*nblk+q] = rk;
        tau_mlid[m*nblk+q] = g - pm->gids_eachrank[rk];
        if (rk == me) {
          tau_loc.h_view(m*nblk+q) = g - g0;
        } else if (q != tau_pos.h_view(m)) {
          tau_any_mpi = true;
        }
      }
    }
    for (int q=0; q<nblk; ++q) {
      if (tau_mrank[m*nblk+q] < 0) {
        TauFatal("the x1 stack of a MeshBlock is incomplete (a non-uniform mesh?)");
      }
    }
  }
  tau_pos.modify_host();
  tau_pos.sync_device();
  tau_loc.modify_host();
  tau_loc.sync_device();

  Kokkos::realloc(tau_col, nmb, n3, n2, n1);
  Kokkos::deep_copy(tau_col, 0.0);
  Kokkos::realloc(tau_ten, nmb, 4, n3, n2, n1);
  Kokkos::deep_copy(tau_ten, 0.0);
  if (nblk > 1) {
    Kokkos::realloc(tau_sum, nmb, n3, n2);
    Kokkos::deep_copy(tau_sum, 0.0);
    Kokkos::realloc(tau_rcv, nmb, nblk, n3, n2);
    Kokkos::deep_copy(tau_rcv, 0.0);
    if (tau_any_mpi) {
      Kokkos::realloc(tau_sum_h, nmb, n3, n2);
      Kokkos::realloc(tau_rcv_h, nmb, nblk, n3, n2);
      Kokkos::deep_copy(tau_rcv_h, 0.0);
    }
  }
  tau_ready = true;
  if (me == 0) {
    std::cout << "<rad_m1> closure = tau: column optical depth from x1max, " << nblk
              << " MeshBlock(s) per x1 stack" << (tau_any_mpi ? " (across ranks)" : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::TauClosureBuild
//! \brief (chi, n) of every active cell into tau_ten, from M1_IW_KT of the work array.

void RadiationM1::TauClosureBuild() {
  if (!tau_ready) {TauClosureInit();}
  Kokkos::fence();
  Kokkos::Timer timer;

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const bool md2 = (indcs.nx2 > 1), md3 = (indcs.nx3 > 1);
  // the columns whose tau is needed: the interior and the x2/x3 FACE ghosts
  const int jl = md2 ? (js - 1) : js, ju = md2 ? (je + 1) : je;
  const int kl = md3 ? (ks - 1) : ks, ku = md3 ? (ke + 1) : ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int nblk = tau_nblk;
  const bool thrd = trans_x3;
  auto iw_ = iw;
  auto tc_ = tau_col;
  auto tt_ = tau_ten;
  auto pos_ = tau_pos;
  auto loc_ = tau_loc;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;

  // (1) the column integral of this block, and the exchange within the stack
  if (nblk > 1) {
    auto ts_ = tau_sum;
    par_for("m1_tau_sum", DevExeSpace(), 0, nmb1, kl, ku, jl, ju,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real s = 0.0;
      for (int i=is; i<=ie; ++i) {
        s += iw_(m,M1_IW_KT,k,j,i)*dx1;
      }
      ts_(m,k,j) = s;
    });
#if MPI_PARALLEL_ENABLED
    if (tau_any_mpi) {
      const int nmb = pmy_pack->nmb_thispack;
      const int me = global_variable::my_rank;
      const int ncol = static_cast<int>(tau_sum.extent(1)*tau_sum.extent(2));
      Kokkos::deep_copy(tau_sum_h, tau_sum);
      std::vector<MPI_Request> req;
      for (int m=0; m<nmb; ++m) {
        int p = tau_pos.h_view(m);
        // receive from every member ABOVE on another rank
        for (int q=p+1; q<nblk; ++q) {
          int rk = tau_mrank[m*nblk+q];
          if (rk == me) continue;
          req.push_back(MPI_REQUEST_NULL);
          MPI_Irecv(&tau_rcv_h(m,q,0,0), ncol, MPI_ATHENA_REAL, rk,
                    TAU_TAG_BASE + tau_mlid[m*nblk+q], MPI_COMM_WORLD, &req.back());
        }
      }
      for (int m=0; m<nmb; ++m) {
        int p = tau_pos.h_view(m);
        // send to every member BELOW on another rank.  Two receivers on one rank get
        // two messages with the same (source, tag) but the SAME payload, so the order
        // in which they are matched does not matter.
        for (int q=0; q<p; ++q) {
          int rk = tau_mrank[m*nblk+q];
          if (rk == me) continue;
          req.push_back(MPI_REQUEST_NULL);
          MPI_Isend(&tau_sum_h(m,0,0), ncol, MPI_ATHENA_REAL, rk, TAU_TAG_BASE + m,
                    MPI_COMM_WORLD, &req.back());
        }
      }
      MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
      Kokkos::deep_copy(tau_rcv, tau_rcv_h);
    }
#endif
  }

  // (2) tau at the cell centres of the interior and the face-ghost columns: the running
  // sum from the cell to the top of the block (the dbg_tensor = tau arithmetic), plus the
  // integrals of the blocks above, summed from the top down
  {
    auto ts_ = tau_sum;
    auto tr_ = tau_rcv;
    par_for("m1_tau_col", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx1 = mbsize.d_view(m).dx1;
      Real t = 0.5*iw_(m,M1_IW_KT,k,j,i)*dx1;
      for (int ii = i+1; ii <= ie; ++ii) {
        t += iw_(m,M1_IW_KT,k,j,ii)*dx1;
      }
      int p = pos_.d_view(m);
      if (p < nblk - 1) {
        Real off = 0.0;
        for (int q=nblk-1; q>p; --q) {
          int ml = loc_.d_view(m*nblk+q);
          Real sq = (ml >= 0) ? ts_(ml,k,j) : tr_(m,q,k,j);
          off = (q == nblk-1) ? sq : (off + sq);
        }
        t += off;
      }
      tc_(m,k,j,i) = t;
    });
  }

  // (3) the tensor of the active cells
  par_for("m1_tau_ten", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real tc = tc_(m,k,j,i);
    // exact grey plane-parallel K/J = (tau + q_inf)/(3 (tau + q(tau))), q = Hopf
    Real qh = 0.710446 - 0.133054*exp(-3.4488*tc);
    Real chi = (tc + 0.710446)/(3.0*(tc + qh));
    Real g1 = -iw_(m,M1_IW_KT,k,j,i);
    Real g2 = 0.0, g3 = 0.0;
    if (md2) {
      Real dx2 = mbsize.d_view(m).dx2;
      BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x2);
      BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x2);
      bool plo = (j == js) && (blo != BoundaryFlag::block) &&
                 (blo != BoundaryFlag::periodic);
      bool phi = (j == je) && (bhi != BoundaryFlag::block) &&
                 (bhi != BoundaryFlag::periodic);
      Real tq = tc_(m,k,j+1,i), tm = tc_(m,k,j-1,i);
      if (plo && phi) {
        g2 = 0.0;
      } else if (plo) {
        g2 = (tq - tc)/dx2;
      } else if (phi) {
        g2 = (tc - tm)/dx2;
      } else {
        g2 = (tq - tm)/(2.0*dx2);
      }
    }
    if (md3) {
      Real dx3 = mbsize.d_view(m).dx3;
      BoundaryFlag blo = mbbcs.d_view(m,BoundaryFace::inner_x3);
      BoundaryFlag bhi = mbbcs.d_view(m,BoundaryFace::outer_x3);
      bool plo = (k == ks) && (blo != BoundaryFlag::block) &&
                 (blo != BoundaryFlag::periodic);
      bool phi = (k == ke) && (bhi != BoundaryFlag::block) &&
                 (bhi != BoundaryFlag::periodic);
      Real tq = tc_(m,k+1,j,i), tm = tc_(m,k-1,j,i);
      if (plo && phi) {
        g3 = 0.0;
      } else if (plo) {
        g3 = (tq - tc)/dx3;
      } else if (phi) {
        g3 = (tc - tm)/dx3;
      } else {
        g3 = (tq - tm)/(2.0*dx3);
      }
    }
    Real gg = thrd ? (g1*g1 + g2*g2 + g3*g3) : (g1*g1 + g2*g2);
    Real ign = 1.0/fmax(sqrt(gg), 1.0e-300);
    tt_(m,0,k,j,i) = chi;
    tt_(m,1,k,j,i) = -g1*ign;
    tt_(m,2,k,j,i) = -g2*ign;
    tt_(m,3,k,j,i) = thrd ? (-g3*ign) : 0.0;
  });

  Kokkos::fence();
  tau_time += timer.seconds();
  tau_ncall += 1.0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::TauClosureReport
//! \brief one cost line at the end of the run (rank 0)

void RadiationM1::TauClosureReport() {
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_m1> closure = tau: " << tau_ncall << " tensor builds, "
            << tau_time << " s (rank 0, fenced), "
            << ((tau_ncall > 0.0) ? (1.0e3*tau_time/tau_ncall) : 0.0)
            << " ms per build" << std::endl;
}

} // namespace radm1
