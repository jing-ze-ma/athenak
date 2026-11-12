//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_consaudit.cpp
//! \brief consaudit-1112 <rad_m1>/dbg_cons_audit (output only): the energy, L_z and
//! radial-momentum change of gas + radiation per segment of the M1 chain (rad_m1.hpp).
//! Reads the state only: a run with the key on is bitwise the run without it.

#include <cmath>
#include <cstdio>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CAuditMark
//! \brief consaudit-1112: called at the top of EAuditMark(pt).  Per cell
//!   q = (U, M, L_z,gas, L_z,rad, P_r,gas, P_r,rad) x V,
//! U = IEN + (c/chat) E, L_z = (r_c x p).z with p the Cartesian vector of the radial +
//! contravariant tangential components (cs: covariant m, F raised with cos_cell; sp:
//! orthonormal), radiation momentum F/(c chat).  The change since the previous mark goes to
//! segment pt (6 = a discarded attempt).  Mark 5 prints CAUDIT and CAUDITB (rank 0) and
//! restarts the sums.  Collective (every rank calls it at the same points).

void RadiationM1::CAuditMark(int pt) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int ni = indcs.nx1, nj = indcs.nx2, nk = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int n3 = u0.extent_int(2), n2 = u0.extent_int(3), n1 = u0.extent_int(4);
  if (caud_q.extent_int(0) != nmb) {
    Kokkos::realloc(caud_q, nmb, 6, n3, n2, n1);
    caud_init = false;
  }
  if (caud_bin.extent_int(0) != 144) {Kokkos::realloc(caud_bin, 144);}
  FluidRef fl = FluidRef::Get(pmy_pack);
  if (!fl.on) return;
  auto uh = fl.u0;
  auto ur = u0;
  auto cq = caud_q;
  auto cb = caud_bin;
  auto vol = pmy_pack->pcoord->volume;
  auto x1v = pmy_pack->pcoord->x1v;
  auto x2v = pmy_pack->pcoord->x2v;
  auto x3v = pmy_pack->pcoord->x3v;
  const bool csm = pmy_pack->pmesh->use_cubed_sphere;
  auto ccell = pmy_pack->pcoord->cos_cell;
  auto &size = pmy_pack->pmb->mb_size;
  auto &mbpan = pmy_pack->pmb->mb_panel;
  const Real ccr = c_light/chat, icc = 1.0/(c_light*chat);
  const bool init = caud_init;
  // bin group: 0 the solve, 1 the fix segments, -1 none
  const int grp = (pt == 3) ? 0 : ((pt == 1 || pt == 2 || pt == 4 || pt == 5) ? 1 : -1);
  if (!init) {Kokkos::deep_copy(cb, 0.0);}   // later resets: after each print
  array_sum::GlobalSum s6;
  Kokkos::parallel_reduce("m1_caudit", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                          nmb*nk*nj*ni),
  KOKKOS_LAMBDA(const int n, array_sum::GlobalSum &sum) {
    const int m = n/(nk*nj*ni);
    const int k = (n/(nj*ni))%nk + ks, j = (n/ni)%nj + js, i = n%ni + is;
    Real rh[3], ta[3], tb[3];
    Real c = 0.0;
    if (csm) {
      const Real xi = 0.25*M_PI*CellCenterX(j-js, nj, size.d_view(m).x2min,
                                             size.d_view(m).x2max);
      const Real eta = 0.25*M_PI*CellCenterX(k-ks, nk, size.d_view(m).x3min,
                                              size.d_view(m).x3max);
      cubed_sphere::PanelToCart(mbpan.d_view(m), xi, eta, rh);
      cubed_sphere::PanelTangents(mbpan.d_view(m), xi, eta, ta, tb);
      c = ccell(m,k,j);
    } else {
      const Real th = x2v(m,j), ph = x3v(m,k);
      rh[0] = sin(th)*cos(ph); rh[1] = sin(th)*sin(ph); rh[2] = cos(th);
      ta[0] = cos(th)*cos(ph); ta[1] = cos(th)*sin(ph); ta[2] = -sin(th);
      tb[0] = -sin(ph); tb[1] = cos(ph); tb[2] = 0.0;
    }
    const Real r = x1v(m,i), v = vol(m,k,j,i), s2 = 1.0 - c*c;
    auto lz = [&](const Real a1, const Real a2c, const Real a3c) {
      const Real a2 = (a2c - c*a3c)/s2, a3 = (a3c - c*a2c)/s2;   // contravariant
      const Real p0 = a1*rh[0] + a2*ta[0] + a3*tb[0];
      const Real p1 = a1*rh[1] + a2*ta[1] + a3*tb[1];
      return r*(rh[0]*p1 - rh[1]*p0);
    };
    Real q[6];
    q[0] = (uh(m,IEN,k,j,i) + ccr*ur(m,M1_E,k,j,i))*v;
    q[1] = uh(m,IDN,k,j,i)*v;
    q[2] = lz(uh(m,IM1,k,j,i), uh(m,IM2,k,j,i), uh(m,IM3,k,j,i))*v;
    q[3] = lz(ur(m,M1_F1,k,j,i), ur(m,M1_F2,k,j,i), ur(m,M1_F3,k,j,i))*v*icc;
    q[4] = uh(m,IM1,k,j,i)*v;
    q[5] = ur(m,M1_F1,k,j,i)*v*icc;
    array_sum::GlobalSum h;
    for (int r2=0; r2<NREDUCTION_VARIABLES; ++r2) h.the_array[r2] = 0.0;
    if (init) {
      Real d[6];
      for (int a=0; a<6; ++a) {
        d[a] = q[a] - cq(m,a,k,j,i);
        h.the_array[a] = d[a];
      }
      if (grp >= 0) {
        const int rb = (8*(i - is))/ni;
        const Real az = fabs(rh[2]);
        const int lb = (az < 0.5) ? 0 : ((az < 0.8660254037844386) ? 1 : 2);
        const int b = rb*3 + lb;
        Kokkos::atomic_add(&cb(grp*72 + b), d[0]);
        Kokkos::atomic_add(&cb(grp*72 + 24 + b), d[2] + d[3]);
        Kokkos::atomic_add(&cb(grp*72 + 48 + b), d[4] + d[5]);
      }
    }
    for (int a=0; a<6; ++a) {cq(m,a,k,j,i) = q[a];}
    sum += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(s6));
  const bool was = caud_init;
  caud_init = true;
  if (!was) return;
  Real loc[6];
  for (int a=0; a<6; ++a) {loc[a] = s6.the_array[a];}
#if MPI_PARALLEL_ENABLED
  Real glb[6];
  MPI_Allreduce(loc, glb, 6, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  for (int a=0; a<6; ++a) {loc[a] = glb[a];}
#endif
  const int sg = (pt >= 0 && pt <= 6) ? pt : 0;
  for (int a=0; a<6; ++a) {caud_seg[sg][a] += loc[a];}
  if (pt != 5) return;
  auto hb = Kokkos::create_mirror_view_and_copy(HostMemSpace(), cb);
  Real bl[144], bg[144];
  for (int b=0; b<144; ++b) {bl[b] = hb(b); bg[b] = bl[b];}
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(bl, bg, 144, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
  Kokkos::deep_copy(cb, 0.0);
  if (global_variable::my_rank == 0) {
    std::printf("CAUDIT cycle=%d t2s=%d dt=%.10e", pmy_pack->pmesh->ncycle, t2_solve,
                dt_sub);
    for (int s=0; s<7; ++s) {
      std::printf(" s%d=%.10e,%.10e,%.10e,%.10e,%.10e,%.10e", s, caud_seg[s][0],
                  caud_seg[s][1], caud_seg[s][2], caud_seg[s][3], caud_seg[s][4],
                  caud_seg[s][5]);
    }
    std::printf("\nCAUDITB cycle=%d", pmy_pack->pmesh->ncycle);
    for (int b=0; b<144; ++b) {std::printf(" %.6e", bg[b]);}
    std::printf("\n");
  }
  for (int s=0; s<7; ++s) {for (int a=0; a<6; ++a) {caud_seg[s][a] = 0.0;}}
}

}  // namespace radm1
