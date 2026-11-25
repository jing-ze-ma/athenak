#ifndef UTILS_FOFC_FLAGS_HPP_
#define UTILS_FOFC_FLAGS_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file fofc_flags.hpp
//! \brief the owner's FOFC flags in the ghosts (<hydro|mhd>/fofc_flag_exchange; the edge
//! ghosts for the corner EMFs).  The
//! trial state that FOFC tests is the update of a cell, also of a ghost cell from the
//! block's own copy of its neighbour: the two differ by the block's positions and areas
//! (round-off) and, where a panel seam meets a block edge, by the seam states.  A trial
//! that sits on a floor (floor_from_ekin: e = e(rho, tfloor)) can then be flagged on one
//! side of a face only: the two sides carry different fluxes and mass and energy are not
//! conserved.  Ghosts across a physical boundary or a cubed-sphere panel seam (resampled
//! ghosts) keep their own flags.

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "bvals/bvals.hpp"
#include "utils/recipe.hpp"

namespace fofc_flags {

//! true where a ghost's trial can differ from its owner's: more than one block (or per
//! panel on the cubed sphere), or a periodic wrap of one block on spherical polar (a
//! Cartesian block that wraps holds the same numbers on both sides).  Not with mesh
//! refinement or a shearing box (no exchange of the flags is defined there).
inline bool Needed(Mesh *pm) {
  bool per = false, shear = false;
  for (int n = 0; n < 6; ++n) {
    per = per || (pm->mesh_bcs[n] == BoundaryFlag::periodic);
    shear = shear || (pm->mesh_bcs[n] == BoundaryFlag::shear_periodic);
  }
  if (pm->multilevel || shear) return false;
  if (pm->use_cubed_sphere) return pm->nmb_total > 6;
  return pm->nmb_total > 1 || (pm->use_spherical_polar && per);
}

//! the boundary object of the flags of <blk> (default true where Needed), else nullptr
inline MeshBoundaryValuesCC *Init(MeshBlockPack *pack, ParameterInput *pin,
                                  const std::string &blk) {
  if (!Needed(pack->pmesh) || pack->pcoord->is_dynamical_relativistic) return nullptr;
  recipe::OldStored(pin, blk, "fofc_flag_exchange", "false", "true");
  if (!pin->GetOrAddBoolean(blk, "fofc_flag_exchange",
          recipe::Fix(pin, blk, "fofc_flag_exchange", "true", "false"))) return nullptr;
  auto *pbv = new MeshBoundaryValuesCC(pack, pin, false);
  pbv->SetVectorPairs(1, {});   // (the flag is a scalar)
  pbv->InitializeBuffers(1);
  return pbv;
}

//! flags of the ghosts of fofc (m,k,j,i) from their owner; fx: scratch, pbv: a
//! boundary object of one scalar
inline void Exchange(MeshBlockPack *pack, MeshBoundaryValuesCC *pbv,
                     DvceArray5D<Real> &fx, DvceArray4D<bool> &fofc) {
  auto &indcs = pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke, ng = indcs.ng;
  const int dj = pack->pmesh->multi_d, dk = pack->pmesh->three_d;
  const int nmb1 = pack->nmb_thispack - 1;
  if (fx.extent_int(0) != nmb1 + 1) {
    Kokkos::realloc(fx, nmb1 + 1, 1, indcs.nx3 + 2*ng*dk, indcs.nx2 + 2*ng*dj,
                    indcs.nx1 + 2*ng);
  }
  auto fx_ = fx;
  auto fofc_ = fofc;
  auto &size = pack->pmb->mb_size;
  const bool cs = pack->pmesh->use_cubed_sphere;
  par_for("FOFC-flag", DevExeSpace(), 0, nmb1, ks-dk, ke+dk, js-dj, je+dj, is-1, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    fx_(m,0,k,j,i) = fofc_(m,k,j,i) ? 1.0 : 0.0;
  });
  DvceArray5D<Real> dum("fofc_xc", 1, 1, 1, 1, 1);   // (no refinement)
  while (pbv->InitRecv(1) == TaskStatus::incomplete) {}
  while (pbv->PackAndSendCC(fx_, dum) == TaskStatus::incomplete) {}
  while (pbv->RecvAndUnpackCC(fx_, dum) == TaskStatus::incomplete) {}
  while (pbv->ClearSend() == TaskStatus::incomplete) {}
  while (pbv->ClearRecv() == TaskStatus::incomplete) {}
  par_for("FOFC-owner", DevExeSpace(), 0, nmb1, ks-dk, ke+dk, js-dj, je+dj, is-1, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const auto &sz = size.d_view(m);
    const bool seam = cs && ((j < js && sz.x2min <= -1.0 + 1.0e-12) ||
                             (j > je && sz.x2max >= 1.0 - 1.0e-12) ||
                             (k < ks && sz.x3min <= -1.0 + 1.0e-12) ||
                             (k > ke && sz.x3max >= 1.0 - 1.0e-12));
    if (!seam) fofc_(m,k,j,i) = (fx_(m,0,k,j,i) != 0.0);
  });
}

}  // namespace fofc_flags
#endif  // UTILS_FOFC_FLAGS_HPP_
