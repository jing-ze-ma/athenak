//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_tau.cpp
//! \brief <rad_vet>/closure = tau: the Eddington tensor of the multi-D implicit solve
//! from the COLUMN OPTICAL DEPTH, never from the cell's own flux.
//!
//!   tau(x1) = int_{x1}^{x1max} rho (kappa_F + kappa_s) dx1'   (from the domain top)
//!   chi     = (tau + q_inf) / (3 (tau + q(tau))),  q(tau) = 0.710446 - 0.133054
//!             exp(-3.4488 tau)  (exact grey plane-parallel K/J, Hopf-function fit)
//!   n       = -grad tau / |grad tau|,  D_ab = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b
//!
//! It is the DIAGNOSTIC <rad_vet>/dbg_tensor = tau of milestone 3i made a real option
//! with MeshBlock stacks along x1 over MPI ranks.  On ONE MeshBlock along x1 the
//! arithmetic is the diagnostic's operation for operation (gate: bitwise identical).
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
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {
void TauFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_vet>/closure = tau: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
// MPI tags of the column-integral exchange: base + the SENDER's local MeshBlock id
// (< 2^14, see NUM_BITS_LID), which keeps every tag <= 32767, the MPI minimum of
// MPI_TAG_UB.  Every message of this exchange is completed before TauClosureBuild
// returns, so nothing of it is in flight together with the gather solve's messages.
constexpr int TAU_TAG_BASE = 16384;
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::TauClosureInit
//! \brief check the mesh, build the x1 stack topology, allocate.  Called on the first
//! TauClosureBuild (by every rank: all ranks run the implicit solve together).

void RadiationVet::TauClosureInit() {
  if (vet_col) {VetColInit(); return;}   // closure = vet_col (rad_vet_vetcol.cpp)
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto &mindcs = pm->mesh_indcs;
  if (!trans_on) {
    TauFatal("needs <rad_vet>/transport = implicit on a MULTI-D mesh (the 1-D branch "
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
    std::cout << "<rad_vet> closure = tau: column optical depth from x1max, " << nblk
              << " MeshBlock(s) per x1 stack" << (tau_any_mpi ? " (across ranks)" : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::TauClosureReport
//! \brief one cost line at the end of the run (rank 0)

void RadiationVet::TauClosureReport() {
  if (vet_col) {VetColReport(); return;}
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_vet> closure = tau: " << tau_ncall << " tensor builds, "
            << tau_time << " s (rank 0, fenced), "
            << ((tau_ncall > 0.0) ? (1.0e3*tau_time/tau_ncall) : 0.0)
            << " ms per build" << std::endl;
}

} // namespace radvet
