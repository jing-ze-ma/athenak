//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_source.cpp
//! \brief curved-mesh rays: set-up, and extinction and source (ln chi, ln S) of every
//! cell.
//!
//!  called from: rad_vet_closure.cpp, rad_vet_tasks.cpp, rays/sc_build.cpp
//!  calls into: rad_vet_closure_rays.cpp, rad_vet_imex.cpp, rays/sc_setup.cpp

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_parfor.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {

void VrayFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_vet> rays: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRayInit
//! \brief checks, buffers, exchange objects and the direction set (ScInit)

void RadiationVet::VetRayInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  if (!trans_on) {VrayFatal("needs <rad_vet>/transport = implicit on a MULTI-D mesh");}
  if (pm->multilevel) {VrayFatal("SMR/AMR is not supported");}
  if (ibc_x1min == VET_IBC_PERIODIC) {VrayFatal("a periodic x1 boundary has no top");}
  // a stretched theta and the full sphere take the node table
  if (!pm->use_cubed_sphere) {
    vgd_pol = pm->use_polar_boundary;
    vgd_sth = pm->use_grid_stretch_theta || vgd_pol;
  }
  // the SC upwind interpolation is bilinear in lateral INDEX space: exact only for
  // uniform lateral spacing (radial stretch is fine).  Refuse a theta stretch, the
  // Athena++ ratio keys, and any non-uniform built xx2f / xx3f (cs panel coordinates
  // and a Cartesian slab are uniform by construction).
  if (vgd_slab != 1) {
    if (!vray_rat_keys.empty()) {
      VrayFatal("needs a laterally UNIFORM mesh; this input sets" + vray_rat_keys +
                " (Athena++ ratio keys: AthenaK does not read them, its theta/phi grid "
                "stays uniform -- remove them, or set 1)");
    }
    auto &ix = pm->mb_indcs;
    auto f2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx2f);
    auto f3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->xx3f);
    // per direction: min and max face spacing over every active face of every block
    double dmin[2] = {1.0e300, 1.0e300}, dmax[2] = {-1.0e300, -1.0e300};
    for (int m = 0; m < pmy_pack->nmb_thispack; ++m) {
      for (int j = ix.js; j <= ix.je; ++j) {
        const double d = f2(m,j+1) - f2(m,j);
        dmin[0] = std::min(dmin[0], d); dmax[0] = std::max(dmax[0], d);
      }
      for (int k = ix.ks; k <= ix.ke; ++k) {
        const double d = f3(m,k+1) - f3(m,k);
        dmin[1] = std::min(dmin[1], d); dmax[1] = std::max(dmax[1], d);
      }
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, dmin, 2, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, dmax, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    for (int a = (vgd_sth ? 1 : 0); a < 2; ++a) {   // theta via the table
      const double rtol = (sizeof(Real) == 4) ? 1.0e-3 : 1.0e-9;  // LeftEdgeX round-off
      if (!(dmin[a] > 0.0) || dmax[a] - dmin[a] > rtol*dmax[a]) {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s cell widths range %.17g .. %.17g",
                      (a == 0) ? "x2 (theta)" : "x3 (phi)", dmin[a], dmax[a]);
        VrayFatal(std::string("needs a laterally UNIFORM mesh; the built grid has ") +
                  buf + " (the SC upwind interpolation is bilinear in (theta, phi) index "
                  "space; only radial stretches are supported)");
      }
    }
  }
  if (indcs.nx1 < 3) {VrayFatal("needs at least 3 cells along x1");}
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int c2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int c3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  Kokkos::realloc(vray_cs, nmb, 2, c3, c2, c1);
  Kokkos::deep_copy(vray_cs, 0.0);
  Kokkos::realloc(vray_cs_c, nmb, 2, 1, 1, 1);
  Kokkos::realloc(vray_cnt, 1);
  Kokkos::deep_copy(vray_cnt, 0.0);
  pbval_vs = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
  pbval_vs->InitializeBuffers(2);
  pbval_vs->SetVectorPairs(2, {});
  // follow <mesh>/cs_seam_blockend as pbval_u does
  pbval_vs->cs_blockend_ext_ = pbval_u->cs_blockend_ext_;
  vray_ready = true;
  ScInit();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRayExchange
//! \brief cell-centred ghost exchange of a (m, n, k, j, i) array

void RadiationVet::VetRayExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                                  MeshBoundaryValuesCC *pb) {
  const int nv = static_cast<int>(a.extent(1));
  while (pb->InitRecv(nv) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(a, ac) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(a, ac) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRaySourceEnd
//! \brief ghosts of vray_cs, then the ray cut (after VetRaySource's source kernel)

void RadiationVet::VetRaySourceEnd() {
  VetRayExchange(vray_cs, vray_cs_c, pbval_vs);
  RayCut();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRaySource
//! \brief source of a vet_gd build: ln chi, ln S of every active cell, their ghosts and
//! the ray cut

void RadiationVet::VetRaySource() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto opac_ = opac;
  auto cs_ = vray_cs;

  const bool thermal = fl_on && coupling && !opac_zero;
  // vet_source_noesrc: the physical E^n = u0
  const bool noes = VsrcNoes();   // rad_vet.hpp
  auto ur_ = u0;
  const bool vsf = vscat;   // vet_scatter
  const Real ar = arad, efl = e_floor;
  par_for_lb("vet_ray_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real chi = fmax(iw_(m,VET_IW_KT,k,j,i), 1.0e-300);
    Real e = fmax(noes ? ur_(m,VET_E,k,j,i) : iw_(m,VET_IW_EN,k,j,i), efl);
    Real s = e;
    if (thermal) {
      Real tg = iw_(m,VET_IW_TP,k,j,i);
      Real t2 = tg*tg;
      Real eth = fmin(opac_(m,VET_OP_P,k,j,i)/chi, 1.0);
      if (vsf) {
        const Real rkp = fmin(opac_(m,VET_OP_P,k,j,i), chi);   // (S between J and B)
        s = VetScatterMix(rkp, chi, ar*t2*t2, e);
      } else {
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
    }
    cs_(m,0,k,j,i) = log(chi);
    cs_(m,1,k,j,i) = log(fmax(s, 1.0e-300));
  });
  VetRaySourceEnd();
}

}  // namespace radvet
