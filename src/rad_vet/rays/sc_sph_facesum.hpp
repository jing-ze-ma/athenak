#ifndef RAD_VET_RAYS_SC_SPH_FACESUM_HPP_
#define RAD_VET_RAYS_SC_SPH_FACESUM_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_sph_facesum.hpp
//! \brief spherical polar: the closure_consistency r-face sums made in the sweep.  The
//! slots of a ray's scalars in vgd_rgs (shell kernel ScSweepSpFs, then ScSpFaceSum).
//!
//!  called from: rays/sc_sph_sweep.cpp, rays/sc_sph_facesum.cpp

// build-time A/B switch: 1 the in-sweep face sums on spherical polar, 0 RayFaceSums'
// upwind rule
#ifndef VET_SP_FS
#define VET_SP_FS 1
#endif

namespace radvet {

// a ray's scalars: the intensity it carries through its r face (xv; -1 not swept this
// pass) and the one where it leaves through the inner or top face (xe, -1 none)
enum ScSpFsSlot {SPFS_XV, SPFS_XE, VGD_SPFS_NS};

}  // namespace radvet

#endif  // RAD_VET_RAYS_SC_SPH_FACESUM_HPP_
