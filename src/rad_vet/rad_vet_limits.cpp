//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_limits.cpp
//! \brief the closure-limit pass of the M1 moments (ApplyClosureLimits): E >= e_floor
//! and |F| <= c E; scheme = clean keeps its solution.

#include <math.h>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "rad_vet/vet_fluid.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::ApplyClosureLimits
//! \brief enforce E >= e_floor and |F| <= c E on the evolved state, over the whole
//! array (active cells and ghosts alike) so that the reconstruction below only ever
//! sees admissible states.

TaskStatus RadiationVet::ApplyClosureLimits(Driver *pdrive, int stage) {
  // implicit_timers: the M1 stage chain starts here
  tmr_on = (tmr_c0 > 0) && (pmy_pack->pmesh->ncycle >= tmr_c0);
  if (tmr_on) {tmr_cnt[7] += 1.0;}
  TmrMark(-1);
  // scheme = clean keeps its solution: a clip would change the energy and momentum totals
  if (clean) {
    TmrMark(0);
    return TaskStatus::complete;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto u0_ = u0;
  Real cl = c_light;
  Real efl = e_floor;

  // STAGE CS2: on the cubed sphere (f2, f3) are covariant: the metric norm
  const bool csl = pmy_pack->pmesh->use_cubed_sphere;
  auto ccl = pmy_pack->pcoord->cos_cell;
  // F4: clip counters over the active cells (output only)
  if (acl_cnt_d.extent_int(0) != 2) {
    Kokkos::realloc(acl_cnt_d, 2);
    Kokkos::deep_copy(acl_cnt_d, 0.0);
  }
  auto ac_ = acl_cnt_d;
  const int ais = indcs.is, aie = indcs.ie, ajs = indcs.js, aje = indcs.je;
  const int aks = indcs.ks, ake = indcs.ke;
  par_for("vet_limits", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = u0_(m,M1_E,k,j,i);
    Real f1 = u0_(m,M1_F1,k,j,i);
    Real f2 = u0_(m,M1_F2,k,j,i);
    Real f3 = u0_(m,M1_F3,k,j,i);
    const Real ein = e, f1in = f1, f2in = f2, f3in = f3;
    const bool lim = csl ? M1ApplyLimitsCs(cl, efl, ccl(m,k,j), e, f1, f2, f3)
                         : M1ApplyLimits(cl, efl, e, f1, f2, f3);
    if (lim && i >= ais && i <= aie && j >= ajs && j <= aje && k >= aks && k <= ake) {
      if (ein < efl) {Kokkos::atomic_add(&ac_(0), 1.0);}
      if (f1 != f1in || f2 != f2in || f3 != f3in) {Kokkos::atomic_add(&ac_(1), 1.0);}
    }
    if (lim) {
      u0_(m,M1_E,k,j,i) = e;
      u0_(m,M1_F1,k,j,i) = f1;
      u0_(m,M1_F2,k,j,i) = f2;
      u0_(m,M1_F3,k,j,i) = f3;
    }
  });
  TmrMark(0);
  return TaskStatus::complete;
}

} // namespace radvet
