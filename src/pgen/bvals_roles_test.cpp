//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file bvals_roles_test.cpp
//! \brief Unit test of the COMPONENT-ROLE TABLE of MeshBoundaryValuesCC
//! (MeshBoundaryValuesCC::SetVectorPairs), on a cubed-sphere or a spherical-polar mesh.
//!
//! A 6-component cell-centred array is exchanged three ways, from identical data:
//!  A: table {{4,5}}  -- slots 0-3 scalars, (4,5) the tangential pair
//!  B: default rule   -- the same fields PERMUTED so the pair sits in (IVY,IVZ) = (2,3),
//!                       which is exactly how the hydro's momentum/velocity is exchanged,
//!                       and the four scalars in slots 0,1,4,5
//!  C: default rule on A's layout (the trap: slots 2,3 of A treated as a vector)
//!  D: table {{2,3}} on B's layout (must equal the default rule bitwise)
//! Checks: A[v] == B[perm(v)] in every cell (scalars pass as the hydro's scalars do, the
//! pair rotates/flips as the hydro's velocity does); D == B bitwise; the exchange filled
//! ghosts; and C differs from A in slots 2-5 (the table is not vacuous on this mesh).
//! The result is printed as one line "BVALS_ROLES ... PASS|FAIL".  Run with nlim = 0.

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "bvals/bvals.hpp"
#include "pgen/pgen.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

void Exchange(MeshBoundaryValuesCC *pb, DvceArray5D<Real> &a, DvceArray5D<Real> &ca,
              const int nvar) {
  while (pb->InitRecv(nvar) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(a, ca) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(a, ca) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
}

}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  if (restart) return;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "bvals_roles_test requires a <hydro> block" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // a quiet hydro state, so the (unused) evolution is harmless
  {
    auto &u0 = pmbp->phydro->u0;
    const Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
    Kokkos::deep_copy(u0, 0.0);
    auto u0_ = u0;
    par_for("bvr_init", DevExeSpace(), 0, pmbp->nmb_thispack-1, 0, u0.extent_int(2)-1,
            0, u0.extent_int(3)-1, 0, u0.extent_int(4)-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      u0_(m,IDN,k,j,i) = 1.0;
      u0_(m,IEN,k,j,i) = 1.0/gm1;
    });
  }

  const int nvar = 6;
  const int nmb = pmbp->nmb_thispack;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const Real sentinel = -7.25;

  // six fields: four scalars s0..s3 and the tangential pair (va, vb), all generic
  // (no symmetry), all functions of (gid, k, j, i)
  auto field = [](const int f, const int g, const int k, const int j, const int i) {
    const Real x = 0.37*i + 0.71*j + 1.13*k + 0.53*g;
    switch (f) {
      case 0: return 1.0 + 0.2*std::sin(x);          // s0
      case 1: return 0.5*std::cos(1.7*x + 0.3);      // s1
      case 2: return 2.0 + std::sin(0.6*x - 1.1);    // s2 (the "F1" of a halo)
      case 3: return -0.8 + 0.3*std::cos(2.3*x);     // s3 (the "KT" of a halo)
      case 4: return std::sin(1.9*x + 0.7);          // va (x2 member)
      default: return 0.6*std::cos(0.9*x - 0.2);     // vb (x3 member)
    }
  };
  // slot of field f in layout A (pair in 4,5) and in layout B (pair in 2,3)
  const int slotA[6] = {0, 1, 2, 3, 4, 5};
  const int slotB[6] = {0, 1, 4, 5, 2, 3};

  auto fill = [&](DvceArray5D<Real> &arr, const int *slot) {
    auto h = Kokkos::create_mirror_view(arr);
    auto &gid = pmbp->pmb->mb_gid;
    for (int m=0; m<nmb; ++m) {
      const int g = gid.h_view(m);
      for (int f=0; f<nvar; ++f) {
        for (int k=0; k<n3; ++k) for (int j=0; j<n2; ++j) for (int i=0; i<n1; ++i) {
          const bool act = (i >= is && i <= ie && j >= js && j <= je &&
                            k >= ks && k <= ke);
          h(m,slot[f],k,j,i) = act ? field(f,g,k,j,i) : sentinel;
        }
      }
    }
    Kokkos::deep_copy(arr, h);
  };

  DvceArray5D<Real> a("bvr_a", nmb, nvar, n3, n2, n1), b("bvr_b", nmb, nvar, n3, n2, n1);
  DvceArray5D<Real> c("bvr_c", nmb, nvar, n3, n2, n1), d("bvr_d", nmb, nvar, n3, n2, n1);
  DvceArray5D<Real> dum("bvr_dum", nmb, nvar, 1, 1, 1);
  fill(a, slotA); fill(c, slotA); fill(b, slotB); fill(d, slotB);

  auto *pba = new MeshBoundaryValuesCC(pmbp, pin, false);
  auto *pbb = new MeshBoundaryValuesCC(pmbp, pin, false);
  auto *pbc = new MeshBoundaryValuesCC(pmbp, pin, false);
  auto *pbd = new MeshBoundaryValuesCC(pmbp, pin, false);
  pba->InitializeBuffers(nvar); pbb->InitializeBuffers(nvar);
  pbc->InitializeBuffers(nvar); pbd->InitializeBuffers(nvar);
  pba->SetVectorPairs(nvar, {{4, 5}});
  pbd->SetVectorPairs(nvar, {{2, 3}});
  Exchange(pba, a, dum, nvar);
  Exchange(pbb, b, dum, nvar);
  Exchange(pbc, c, dum, nvar);
  Exchange(pbd, d, dum, nvar);

  auto ha = Kokkos::create_mirror_view_and_copy(HostMemSpace(), a);
  auto hb = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b);
  auto hc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), c);
  auto hd = Kokkos::create_mirror_view_and_copy(HostMemSpace(), d);
  // counters: [0] ghosts filled, [1] A vs B(perm) mismatches, [2] D vs B mismatches,
  // [3] ghosts where C != A in the scalar slots 2,3, [4] where C != A in the pair 4,5
  int64_t cnt[5] = {0, 0, 0, 0, 0};
  for (int m=0; m<nmb; ++m) {
    for (int k=0; k<n3; ++k) for (int j=0; j<n2; ++j) for (int i=0; i<n1; ++i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      for (int f=0; f<nvar; ++f) {
        const Real va = ha(m,slotA[f],k,j,i), vb = hb(m,slotB[f],k,j,i);
        if (!act && f == 0 && va != sentinel) {cnt[0]++;}
        if (va != vb) {cnt[1]++;}
        if (hd(m,f,k,j,i) != hb(m,f,k,j,i)) {cnt[2]++;}
      }
      if (!act) {
        if (hc(m,2,k,j,i) != ha(m,2,k,j,i) || hc(m,3,k,j,i) != ha(m,3,k,j,i)) cnt[3]++;
        if (hc(m,4,k,j,i) != ha(m,4,k,j,i) || hc(m,5,k,j,i) != ha(m,5,k,j,i)) cnt[4]++;
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, cnt, 5, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
#endif
  const bool pass = (cnt[0] > 0 && cnt[1] == 0 && cnt[2] == 0 && cnt[3] > 0 &&
                     cnt[4] > 0);
  if (global_variable::my_rank == 0) {
    std::cout << "BVALS_ROLES geometry="
              << (pmy_mesh_->use_cubed_sphere ? "cs" :
                  (pmy_mesh_->use_spherical_polar ? "sp" : "cart"))
              << " ghosts_filled=" << cnt[0]
              << " table_vs_hydro_rule_mismatch=" << cnt[1]
              << " pair23_table_vs_default_mismatch=" << cnt[2]
              << " default_rule_moves_scalars_2_3=" << cnt[3]
              << " default_rule_misses_pair_4_5=" << cnt[4]
              << (pass ? " PASS" : " FAIL") << std::endl;
  }
  delete pba; delete pbb; delete pbc; delete pbd;
  if (!pass) {std::exit(EXIT_FAILURE);}
}
