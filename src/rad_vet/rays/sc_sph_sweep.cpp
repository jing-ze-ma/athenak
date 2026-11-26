//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_sweep.cpp
//! \brief <rad_vet>/vet_gd: the Eddington tensor on the spherical-polar grid from a
//! short-characteristics formal solution along globally fixed (Cartesian) directions.
//!
//! The grey equation n . grad I = chi (S - I) is swept shell by shell for a fixed set of
//! straight-ray directions (HEALPix or product set; weights integrate 1, n, nn exactly).
//! Two passes, one launch per shell, one thread per (block, cell, direction):
//!   inward (n . r_hat < 0), top-down: upwind point on shell i+1 (vacuum at the top);
//!   outward (n . r_hat >= 0), bottom-up: shell i-1, or beyond the tangent point on the
//!     same shell (inward value), or the inner face (diffusion intensity).
//! Upwind I, ln chi, ln S are bilinear in (theta, phi) index space; S is linear in tau.
//! Periodic theta/phi walls: ghost intensities are re-indexed by the local-frame map
//! n -> M n (nearest set direction). Block boundaries: lagged inflow (block-Jacobi).
//! Output: the intensities vgd_i (ClosureFromSc, RayFaceSums).
//!
//!  called from: rad_vet_closure_ring.cpp, rays/sc_build.cpp
//!  calls into: rad_vet_params_closure.cpp, rays/sc_bandring.cpp, rays/sc_cart_setup.cpp,
//!    rays/sc_cs_sweep.cpp, rays/sc_setup.cpp, rays/sc_sph_globalband.cpp,
//!    rays/sc_sph_halo.cpp, rays/sc_sph_halo_exact.cpp, rays/sc_sph_halo_fused.cpp,
//!    rays/sc_sph_walls.cpp, rays/sc_tables.cpp

#include <math.h>

#include <algorithm>
#include <array>
#include <type_traits>
#include <vector>

#include "athena.hpp"
// FMA contraction only within one expression in this TU
#if defined(__clang__)
#pragma clang fp contract(on)
#endif
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rays/sc_ray.hpp"
#include "rad_vet/rays/sc_halo.hpp"
#include "rad_vet/rays/sc_cs.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ScSweep
//! \brief one sweep (inward then outward) of every direction

void RadiationVet::ScSweep() {
  if (vgd_csg) {   // the cubed-sphere sweep (rays/sc_cs_sweep.cpp)
    ScCsSweep();
    return;
  }
  if (vgd_slab == 1) {
    ScSweepSlab();
    return;
  }
  ScSweepSp();
}

#define VET_VGD_SWEEP_NAME ScSweepSp
#define VET_VGD_SWEEP_SLAB 0
#include "rad_vet/rays/sc_sph_sweep_shell.inc"
#undef VET_VGD_SWEEP_NAME
#undef VET_VGD_SWEEP_SLAB
#define VET_VGD_SWEEP_NAME ScSweepSlab
#define VET_VGD_SWEEP_SLAB 1
#include "rad_vet/rays/sc_sph_sweep_shell.inc"
#undef VET_VGD_SWEEP_NAME
#undef VET_VGD_SWEEP_SLAB

}  // namespace radvet
