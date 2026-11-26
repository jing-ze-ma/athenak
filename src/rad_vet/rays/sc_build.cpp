//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_build.cpp
//! \brief curved-mesh rays: one build (source, tables, sweep, moments).
//!
//!  called from: rad_vet_closure.cpp, rad_vet_closure_ring.cpp
//!  calls into: rays/sc_cs_setup.cpp, rays/sc_directions.cpp, rays/sc_setup.cpp,
//!    rays/sc_source.cpp, rays/sc_sph_halo.cpp, rays/sc_sph_sweep.cpp, rays/sc_tables.cpp

#include <iostream>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_halo.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSource
//! \brief source of a gd build: ln chi, ln S of every cell (VetRaySource) into the band
//! index space of vgd_cs with its lateral band filled exactly (one multi-shell
//! exchange); with several MeshBlocks along x1 also the x1 ghost shells

void RadiationVet::ScSource() {
  VetRaySource();
  // from the first shell up (gd sweep), or every shell for the clean scheme, whose rays
  // and closure read the whole column
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int og = vgd_w - indcs.ng, ok = vgd_w - indcs.ks, x1g = (vgd_nb1 > 1) ? 1 : 0;
  const int cut = 0;
  const int ilo = indcs.is + cut - x1g, ie = indcs.ie + x1g;
  auto cs_ = vray_cs;
  auto cw_ = vgd_cs;
  par_for("vet_vgd_csw", DevExeSpace(), 0, pmy_pack->nmb_thispack - 1, indcs.ks,
          indcs.ke, indcs.js, indcs.je, ilo, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    cw_(m,0,k+ok,j+og,i) = cs_(m,0,k,j,i);
    cw_(m,1,k+ok,j+og,i) = cs_(m,1,k,j,i);
  });
  vgd_cur = DevExeSpace();
  if (vgd_csg) {   // the cs band copies
    ScCsSource();
    return;
  }
  ScHaloCs(vgd_cs, ilo, ie);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScBuild
//! \brief the gd build: source and first shell (VetRaySource),
//! sweep(s) with the lagged lateral inflow, moments; timed per part (fenced)

void RadiationVet::ScBuild() {
  Kokkos::Timer tm;
  ScSource();
  Kokkos::fence();
  vgd_tsrc += tm.seconds();
  ScTabPoll();      // the prefetched tables, uploaded when ready
  // vet_gd_rotate_every: a new z-angle at the first build of a block of N cycles; the
  // stored intensities belong to the old directions, so the lagged inflow re-converges
  const Real al = RotAngle(pmy_pack->pmesh->ncycle);
  // no vet_rotate_avg or cross-fade here: the clean scheme sweeps rotating sets on its
  // ring (RingMoments); without ring they are refused (InitArrays)
  if (al != vgd_alpha) {
    Kokkos::Timer tt;
    ScTables(al);
    vgd_ttab += tt.seconds();
  }
  ScTabPrefetch(pmy_pack->pmesh->ncycle);
  ScSweepPost();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSweepPost
//! \brief sweep part of a gd build with the current tables, and ScPost

void RadiationVet::ScSweepPost() {
  Kokkos::Timer tm;
  // the per-shell halo makes one sweep exact
  vgd_cur = DevExeSpace();
  vgd_scut = vray_icut;
  tm.reset();
  ScSweep();
  Kokkos::fence();
  vgd_tswp += tm.seconds();
  ScPost();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScPost
//! \brief after the sweep(s) of a build (also ring's): the GPU halo's error
//! flag, the clamp check over all ranks

void RadiationVet::ScPost() {
  Kokkos::fence();
  ScIpcCheck();
  vray_ncall += 1.0;
  // persistent pinned buffer
  if (vgd_pin_cnt.extent(0) != vray_cnt.extent(0)) {
    vgd_pin_cnt = Kokkos::View<Real*, Kokkos::SharedHostPinnedSpace>("vet_vgd_pin_cnt",
                                                                     vray_cnt.extent(0));
  }
  Kokkos::deep_copy(vgd_pin_cnt, vray_cnt);
  vray_nclamp = vgd_pin_cnt(0);
  // all ranks (sweep exactness: a read beyond the band is clamped and counted)
  vgd_nclamp_all = vray_nclamp;
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &vgd_nclamp_all, 1, MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  // the band is capped at the block size: a read beyond it is clamped (inexact).
  // Fatal unless vet_gd_allow_clamp (warning); on the cubed sphere every clamped read
  // counts (the depth there is a bound with margin, not an exact reach)
  if ((vgd_capped || vgd_csg) && vgd_nclamp_all > vgd_nclamp_seen) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_gd: " << (vgd_nclamp_all - vgd_nclamp_seen)
                << " lateral reads were CLAMPED to the band (all ranks) at cycle "
                << pmy_pack->pmesh->ncycle << ": the MeshBlocks are smaller than the "
                << "near-tangent reach, the sweep is not exact" << std::endl;
    }
    if (!vgd_allow_clamp) {
      ScFatal("clamped lateral reads (MeshBlocks too small for the reach); use larger "
              "MeshBlocks or <rad_vet>/vet_gd_allow_clamp = true");
    }
    vgd_nclamp_seen = vgd_nclamp_all;
  }
}

}  // namespace radvet
