//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetlat.cpp
//! \brief the vet_col_lat interface that vet_gd rides on (rad_vet_vetgd.cpp): the
//! start-up checks and buffers (VetLatInit), ln chi and ln S of every cell with their
//! ghosts and the ray cut (VetLatSweep), and the largest lateral D_r,lat
//! (VetLatOdMax).  The lateral tensor itself is vet_gd's (tau_ten slots M1_TT_LAT0..).
//!
//! LIMITS (fatal at start-up): the spherical-polar wedge (sph_geom) or the cubed sphere
//! through vet_gd, vet_col_axis = radial, a laterally uniform mesh.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet_parfor.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_vetseg.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {
void VlatFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_vet>/vet_col_lat: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatInit
//! \brief checks, buffers, exchange objects (called at the end of VetColInit)

void RadiationVet::VetLatInit() {
  Mesh *pm = pmy_pack->pmesh;
  // (CS_VETGD.md stage S3): vet_gd runs on the cubed sphere (its own
  // sweep, rad_vet_vetgd_cs.cpp); the vet_col_lat SC correction itself does not
  if (pm->use_cubed_sphere && !vgd_on) {
    VlatFatal("the cubed sphere: only through vet_gd; the vet_col_lat SC "
              "correction is wedge-only");
  }
  if (!sph_geom) {
    VlatFatal("needs the spherical-polar wedge (Cartesian: closure = vet_sc)");
  }
  // S1: rad_vet runs on the full sphere (polar x2 faces) with vet_col only.
  // The SC reads of vet_col_lat / vet_gd wrap theta periodically at block level, their
  // band formula has 1/sin(theta_min) and no pole partner map: stage S2 (latitude cap).
  // vet_gd on the spherical-polar mesh takes a stretched theta and the full
  // sphere through its node table (rad_vet_vetgd_sph2.inc); vet_col_lat keeps both
  // refusals
  if (vgd_on && !pm->use_cubed_sphere) {
    vgd_pol = pm->use_polar_boundary;
    vgd_sth = pm->use_grid_stretch_theta || vgd_pol;
  }
  if (pm->use_polar_boundary && !vgd_sth) {
    VlatFatal(std::string(vgd_on ? "vet_gd" : "vet_col_lat (lateral SC)") +
              " on the FULL sphere (mesh/use_polar_boundary) is not "
              "implemented: use closure = vet_col without vet_gd / vet_col_lat");
  }
  // (USER 10-17): the SC upwind interpolation of vet_col_lat / vet_gd is
  // bilinear in (theta, phi) INDEX space, which is exact only for UNIFORM lateral
  // spacing (a radial stretch is fine: the sweep uses the physical radii).  The only
  // lateral stretch of this code is mesh/use_grid_stretch_theta (rad_vet.cpp refuses it
  // for every M1 run today; refused here again so lifting that guard cannot silently
  // re-enable it); the Athena++ keys x2rat / x3rat are not read by AthenaK; and the
  // built face arrays xx2f / xx3f are checked so that any future lateral map is caught.
  // on the cubed sphere xx2f / xx3f are the panel coordinates on [-1, 1]
  // (equiangular: xi = pi/4 x2, coordinates.cpp), uniform by construction, and the cs
  // sweep's reads are bilinear in each panel's equiangular INDEX space: the same
  // requirement, so the same check applies unchanged (it cannot fire on a cs mesh).
  {
    const std::string who = vgd_on ? "vet_gd" : "vet_col_lat (lateral SC)";
    if (pm->use_grid_stretch_theta && !vgd_sth) {
      VlatFatal(who + " needs a theta-UNIFORM mesh: mesh/use_grid_stretch_theta = true "
                "(the SC upwind interpolation is bilinear in (theta, phi) index space; "
                "only radial stretches are supported)");
    }
    if (!vlat_rat_keys.empty()) {
      VlatFatal(who + " needs a laterally UNIFORM mesh; this input sets" + vlat_rat_keys +
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
        VlatFatal(who + " needs a laterally UNIFORM mesh; the built grid has " + buf +
                  " (the SC upwind interpolation is bilinear in (theta, phi) index "
                  "space; only radial stretches are supported)");
      }
    }
  }
  if (vcol_axis_flux) {VlatFatal("needs vet_col_axis = radial");}
  auto &indcs = pm->mb_indcs;
  if (indcs.nx1 < 3) {VlatFatal("needs at least 3 cells along x1");}
  const int nmb = pmy_pack->nmb_thispack;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int c2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int c3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  Kokkos::realloc(vlat_cs, nmb, 2, c3, c2, c1);
  Kokkos::deep_copy(vlat_cs, 0.0);
  Kokkos::realloc(vlat_cs_c, nmb, 2, 1, 1, 1);
  Kokkos::realloc(vlat_cnt, 1);
  Kokkos::deep_copy(vlat_cnt, 0.0);
  pbval_vs = new MeshBoundaryValuesCC(pmy_pack, nullptr, false);
  pbval_vs->InitializeBuffers(2);
  pbval_vs->SetVectorPairs(2, {});
  // the resampled lateral-VET exchanges follow <mesh>/cs_seam_blockend
  // as the moments' pbval_u does
  pbval_vs->cs_blockend_ext_ = pbval_u->cs_blockend_ext_;
  vlat_ready = true;
  if (vgd_on) {VetGdInit(); return;}
  // without vet_gd no lateral sweep runs (the LAT slots stay 0)
  if (global_variable::my_rank == 0) {
    std::cout << "### NOTE: <rad_vet>/vet_col_lat = true has no effect without vet_gd"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatOdMax
//! \brief vlat_odmax = max over ranks of |D_rt|, |D_rp| (tau_ten LAT slots 1, 2) over the
//! active cells: 0 keeps the operator form untouched (vet_col_lat and vet_gd)

void RadiationVet::VetLatOdMax() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto tt_ = tau_ten;
  // + the tangential slots (vet_gd); on the cubed sphere D_r,lat is the Cartesian
  // vector in LAT1..LAT3, + Q (LAT4, LAT5)
  const bool c3 = cs_geom;
  const bool tn = vgd_on && vgd_tan && !c3, cst = vgd_on && vgd_tan && c3;
  // cs_vertex_mono (rad_vet.hpp): D_r,lat (LAT1..3) = 0 in the cells within
  // cs_vertex_mono_n + 1 of a cube vertex, set after every build and before the ghost
  // exchange (the callers), so the operator, its stencil, the lagged assembly and the
  // time2 slot extrapolation all see the same zero, and every face of a cell within
  // cs_vertex_mono_n (incl. the seam faces' foot values) carries no D_r,lat
  if (c3 && cs_vtx_mono) {
    auto cse = cs_seam.d_view;
    const int vr = cs_vtx_n + 1;
    par_for("vet_vlat_csvtx", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (M1CsVtxDist(k, j, ks, ke, js, je, cse(m,0) != 0, cse(m,1) != 0, cse(m,2) != 0,
                      cse(m,3) != 0) <= vr) {
        tt_(m,M1_TT_LAT0+1,k,j,i) = 0.0;
        tt_(m,M1_TT_LAT0+2,k,j,i) = 0.0;
        tt_(m,M1_TT_LAT0+3,k,j,i) = 0.0;
        if (cst) {
          tt_(m,M1_TT_LAT0+4,k,j,i) = 0.0;
          tt_(m,M1_TT_LAT0+5,k,j,i) = 0.0;
        }
      }
    });
  }
  {
    Real omx = 0.0;
    const int nk = ke - ks + 1, nj = je - js + 1, ni = ie - is + 1;
    Kokkos::parallel_reduce("vet_vlat_odmax",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nk*nj*ni),
      KOKKOS_LAMBDA(const int idx, Real &mx) {
        int t = idx/ni;
        const int i = is + (idx - t*ni);
        const int j = js + (t % nj);
        t /= nj;
        const int k = ks + (t % nk);
        const int m = t/nk;
        Real v = fmax(fabs(tt_(m,M1_TT_LAT0+1,k,j,i)),
                      fabs(tt_(m,M1_TT_LAT0+2,k,j,i)));
        if (tn) {
          v = fmax(v, fmax(fabs(tt_(m,M1_TT_LAT0+3,k,j,i)),
                           fabs(tt_(m,M1_TT_LAT0+4,k,j,i))));
        }
        if (c3) {v = fmax(v, fabs(tt_(m,M1_TT_LAT0+3,k,j,i)));}
        if (cst) {
          v = fmax(v, fmax(fabs(tt_(m,M1_TT_LAT0+4,k,j,i)),
                           fabs(tt_(m,M1_TT_LAT0+5,k,j,i))));
        }
        mx = (v > mx) ? v : mx;
      }, Kokkos::Max<Real>(omx));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &omx, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
    vlat_odmax = omx;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatExchange
//! \brief the ordinary cell-centred ghost exchange of a (m, n, k, j, i) array

void RadiationVet::VetLatExchange(DvceArray5D<Real> &a, DvceArray5D<Real> &ac,
                                  MeshBoundaryValuesCC *pb) {
  const int nv = static_cast<int>(a.extent(1));
  while (pb->InitRecv(nv) == TaskStatus::incomplete) {}
  while (pb->PackAndSendCC(a, ac) == TaskStatus::incomplete) {}
  while (pb->RecvAndUnpackCC(a, ac) == TaskStatus::incomplete) {}
  while (pb->ClearSend() == TaskStatus::incomplete) {}
  while (pb->ClearRecv() == TaskStatus::incomplete) {}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatCsFinish
//! \brief after the source kernel of VetLatSweep: the ghosts of vlat_cs, then the
//! ray cut

void RadiationVet::VetLatCsFinish() {
  VetLatExchange(vlat_cs, vlat_cs_c, pbval_vs);
  CleanRayCut();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatSweep
//! \brief the source of a vet_gd build: ln chi, ln S of every active cell (vet_col's
//! extinction and source), their ghosts and the ray cut (VetLatCsFinish)

void RadiationVet::VetLatSweep() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto iw_ = iw;
  auto opac_ = opac;
  auto cs_ = vlat_cs;
  auto cx1f = pmy_pack->pcoord->xx1f;
  auto &mbsize = pmy_pack->pmb->mb_size;

  const bool thermal = fl_on && coupling && !opac_zero;
  // vet_source_noesrc (F1): the physical E^n = u0 (see VetColBuild)
  const bool noes = VsrcNoes();   // rad_vet.hpp
  auto ur_ = u0;
  const int vsf = vscat ? vscat_form : 0;   // vet_scatter
  const Real ar = arad, efl = e_floor;
  const bool srcep = vlat_src_ep;   // vet_gd_rebuild_every: the Picard iterate's E
  const bool cvx = cln_srccvx;   // clean_source_convex (the relaxed source: removed)
  par_for_lb("vet_vlat_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real chi = fmax(iw_(m,M1_IW_KT,k,j,i), 1.0e-300);
    Real e = fmax((noes && !srcep) ? ur_(m,M1_E,k,j,i)
                                   : iw_(m,srcep ? M1_IW_EP : M1_IW_EN,k,j,i), efl);
    Real s = e;
    if (thermal) {
      Real tg = iw_(m,M1_IW_TP,k,j,i);
      Real t2 = tg*tg;
      Real eth = fmin(opac_(m,M1_OP_P,k,j,i)/chi, 1.0);
      const Real rks = (vsf != 0) ? opac_(m,M1_OP_S,k,j,i) : 0.0;
      if (vsf != 0) {
        const Real rkp = cvx ? fmin(opac_(m,M1_OP_P,k,j,i), chi) : opac_(m,M1_OP_P,k,j,i);
        s = VetScatterMix(vsf, rkp, chi, rks, ar*t2*t2, e);
      } else {
        s = eth*ar*t2*t2 + (1.0 - eth)*e;
      }
    }
    cs_(m,0,k,j,i) = log(chi);
    cs_(m,1,k,j,i) = log(fmax(s, 1.0e-300));
  });
  VetLatCsFinish();
}

} // namespace radvet
