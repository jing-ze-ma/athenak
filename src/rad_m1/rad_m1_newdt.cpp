//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_newdt.cpp
//! \brief radiation substep timestep, dt_rad = cfl_rad*min(dx)/chat (design sect. 5).
//!
//! NOTE the CFL number is the module's own <rad_m1>/cfl_rad and is already contained in
//! dtnew; Mesh::NewTimeStep therefore uses dtnew directly and does NOT multiply it by
//! <time>/cfl_number again.

#include <float.h>

#include <algorithm>
#include <limits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::NewTimeStep

TaskStatus RadiationM1::NewTimeStep(Driver *pdriver, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  Real dt1 = std::numeric_limits<float>::max();
  Real dt2 = std::numeric_limits<float>::max();
  Real dt3 = std::numeric_limits<float>::max();

  auto &size = pmy_pack->pmb->mb_size;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;

  Kokkos::parallel_reduce("RadM1Nudt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
    int m = (idx)/nkji;
    min_dt1 = fmin((size.d_view(m).dx1), min_dt1);
    min_dt2 = fmin((size.d_view(m).dx2), min_dt2);
    min_dt3 = fmin((size.d_view(m).dx3), min_dt3);
  }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2), Kokkos::Min<Real>(dt3));

  Real dxmin = dt1;
  if (pmy_pack->pmesh->multi_d) {dxmin = std::min(dxmin, dt2);}
  if (pmy_pack->pmesh->three_d) {dxmin = std::min(dxmin, dt3);}
  if (sph_geom) {
    // stage S1, spherical polar: mb_size.dx2/dx3 are ANGLES and dx1 is the unstretched
    // spacing, so take the smallest physical cell width dr, r dtheta, r sin(theta) dphi
    // of Coordinates over the active cells (reduced over ranks by the driver)
    auto cx1 = pmy_pack->pcoord->dx1;
    auto cx2 = pmy_pack->pcoord->dx2;
    auto cx3 = pmy_pack->pcoord->dx3;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nji = nx2*nx1;
    const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
    Real smin = std::numeric_limits<float>::max();
    Kokkos::parallel_reduce("RadM1NudtSph",
                            Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &lmin) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/nx1;
      int i = r - j*nx1;
      k += ks; j += js; i += is;
      Real d = cx1(m,k,j,i);
      if (md) {d = fmin(d, cx2(m,k,j,i));}
      if (td) {d = fmin(d, cx3(m,k,j,i));}
      lmin = fmin(d, lmin);
    }, Kokkos::Min<Real>(smin));
    dxmin = smin;
  }

  // ---- MILESTONE 3a hook: the implicit scheme has no light-speed CFL.  It takes the
  // step it is given; <rad_m1>/implicit_cfl > 0 asks for a radiation CFL c dt/dx, which
  // is how the radiation-only tests choose their step.
  if (transport >= M1_TRANSPORT_IMPLICIT_X1) {
    dtnew = (impl_cfl > 0.0) ? (impl_cfl*dxmin/c_light)
                             : static_cast<Real>(std::numeric_limits<float>::max());
    return TaskStatus::complete;
  }
  // ---- end of the 3a hook

  dtnew = cfl_rad*dxmin/chat;
  return TaskStatus::complete;
}

} // namespace radm1
