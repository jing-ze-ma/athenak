//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetcol.cpp
//! \brief Set-up read by the vet_gd closure (VetColInit: checks, outer Marshak q
//! vcol_q, spherical Marshak coefficients) and the end-of-run cost report; then the checks,
//! the buffers and the source (ln chi, ln S) of the vet_gd rays.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/vet_fluid.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_vetseg.hpp"
#include "rad_vet/rad_vet_parfor.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {
void VcolFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_vet>/closure = vet_col: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetColInit
//! \brief the checks, the spherical Marshak coefficients and vcol_q

void RadiationVet::VetColInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  if (!trans_on) {VcolFatal("needs <rad_vet>/transport = implicit on a MULTI-D mesh");}
  if (pm->multilevel) {VcolFatal("SMR/AMR is not supported");}
  if (ibc_x1min == VET_IBC_PERIODIC) {VcolFatal("a periodic x1 boundary has no top");}
  if (indcs.nx1 < 2) {VcolFatal("needs at least 2 cells along x1");}
  vcol_sph = sph_geom;
  vcol_sqf = vcol_sqf && vcol_sph && vcol_sq && impl_mface_lin;
  const int ie = indcs.ie;
  if (vcol_sph) {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    const double rtop = x1f_h(0, ie + 1);
    VetSphMarshakCoef(x1v_h(0, ie), x1v_h(0, ie - 1), rtop, vcol_jca, vcol_jcb);
    vcol_jcap = 2.0*SQR(x1v_h(0, ie)/rtop);
  }
  const int nmb = pmy_pack->nmb_thispack;
  int c2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int c3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  if (vcol_sq) {
    Kokkos::realloc(vcol_q, nmb, c3, c2);
    Kokkos::deep_copy(vcol_q, marshak_q);
    if (ibc_x1max != VET_IBC_MARSHAK && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_col_surface_q: the outer x1 boundary is not marshak; "
                << "the q is built but unused" << std::endl;
    }
  }
  tau_ready = true;
  VetLatInit();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetColReport
//! \brief one cost line at the end of the run (rank 0)

void RadiationVet::VetColReport() {
  if (global_variable::my_rank != 0 || !vgd_on) return;
  std::cout << "<rad_vet> vet_gd: " << vlat_ncall << " builds checked; lateral "
            << "reads clamped to the ghost band (rank 0, all sweeps): " << vlat_nclamp
            << std::endl;
  std::cout << "<rad_vet> vet_gd: rank 0 fenced seconds: source " << vgd_tsrc
            << ", sweeps " << vgd_tswp << " (of which per-shell halo+walls "
            << vgd_thalo << ", of which MPI wait " << vgd_tmpi
            << "), tables (rotation) " << vgd_ttab
            << "; in-solve rebuilds " << vgd_nrb << " ("
            << vgd_trb << " s)" << std::endl;
  std::cout << "<rad_vet> vet_gd halo (rank 0): exchanges " << vgd_nexch
            << ", bytes sent "
            << vgd_nbyte << " (" << ((vgd_nexch > 0.0) ? vgd_nbyte/vgd_nexch : 0.0)
            << " per exchange), post " << vgd_tpost << " s, partners "
            << vgd_prk.size() << "; band-clamped lateral reads (ALL ranks, all "
            << "sweeps) "
            << vgd_nclamp_all << std::endl;
  if (vgd_csg || vgd_gbc) {   // cs or sp global band
    std::cout << "<rad_vet> vet_gd " << (vgd_csg ? "cs" : "global")
              << " band (rank 0): " << ((vgd_cscomp || vgd_gbc) ? "compact" : "dense")
              << ", intensity values sent "
              << vgd_csvsent << " of dense " << vgd_csvdense << " ("
              << ((vgd_csvdense > 0.0) ? vgd_csvsent/vgd_csvdense : 0.0)
              << "); list builds " << vgd_csnbuild << " (" << vgd_cstbuild << " s)"
              << std::endl;
  }
  if (vgd_csg && vgd_cscomp) {   // where the list builds go
    std::cout << "<rad_vet> cs list builds (rank 0 s): mask kernel "
              << vgd_cstb_k << ", download+check " << vgd_cstb_d << ", mark exchange "
              << vgd_cstb_m << ", host lists " << vgd_cstb_h << ", upload "
              << vgd_cstb_u << "; mask stride " << vgd_csper << " (padded)"
              << std::endl;
  }
  {
    // ragged band: accesses outside a shell's slab (must be 0)
    Real no = 0.0;
    if (vgd_i.oos.extent(0) > 0) {
      auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i.oos);
      no = h(0);
    }
    std::cout << "<rad_vet> vet_gd ragged band (rank 0): out-of-slab accesses "
              << no << std::endl;
    if (vgd_ring) {
      std::cout << "<rad_vet> vet_gd_band_ring (rank 0): sweeps " << vgd_rg_nk
                << ", dense band saves " << vgd_rg_ndense << ", save store "
                << vgd_rg_store.extent(0)*sizeof(Real)/1.0e9 << " GB" << std::endl;
    }
    if (vgd_xff) {
      double oos2 = 0.0;
      if (vgd_i_xb.oos.extent(0) > 0) {
        auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), vgd_i_xb.oos);
        oos2 = h(0);
      }
      std::cout << "<rad_vet> vet_gd_rotate_xfade_fuse (rank 0): fused builds "
                << vgd_xff_n << ", two-sweep fallbacks " << vgd_xff_nfb
                << ", merged exchanges " << vgd_xff_nmsg << ", out-of-slab accesses "
                << "(other set) " << oos2 << std::endl;
    }
  }
  if (vgd_hcomp > 0) {
    std::cout << "<rad_vet> (rank 0) vet_gd halo lists "
              << vgd_hl_nmade << " made, " << vgd_hl_bytes/1048576.0 << " MB of "
              << vgd_hl_cap/1048576.0 << " (key " << vgd_hl_mb
              << ", GPU: <= free/2); exchanges by list " << vgd_hl_nuse << ", dense "
              << vgd_hl_ndense << std::endl;
  }
  if (vgd_hx_on) {
    std::cout << "<rad_vet> (rank 0) vet_gd_halo_exact lists made "
              << vgd_hx_nmade << ", fallbacks " << vgd_hx_nfall
              << " (marked reads outside the masks' region " << vgd_hx_nout
              << "), read ghosts with an unsent source " << vgd_hx_nun
              << "; received values per list: exact/mask "
              << ((vgd_hx_vmask > 0.0) ? vgd_hx_vexact/vgd_hx_vmask : 0.0)
              << std::endl;
  }
  std::cout << "<rad_vet> (rank 0) exact-list batches " << vgd_hxb_ngroup
            << " (group " << vgd_hxb_g << ", " << vgd_hxb_t << " s, list reallocs "
            << vgd_hxb_nrealloc << "); table prefetch 1"
            << " swapped in " << vgd_tab_nswap << " of the rotations" << std::endl;
}

//----------------------------------------------------------------------------------------
// lateral VET (vet_col_lat) and the vet_gd source

namespace {
void VlatFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<rad_vet>/vet_gd: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatInit
//! \brief checks, buffers, exchange objects (called at the end of VetColInit)

void RadiationVet::VetLatInit() {
  Mesh *pm = pmy_pack->pmesh;
  if (!sph_geom) {
    VlatFatal("needs the spherical-polar wedge (Cartesian: closure = vet_sc)");
  }
  // a stretched theta and the full sphere take the node table
  if (!pm->use_cubed_sphere) {
    vgd_pol = pm->use_polar_boundary;
    vgd_sth = pm->use_grid_stretch_theta || vgd_pol;
  }
  // the SC upwind interpolation is bilinear in lateral INDEX space: exact only for
  // uniform lateral spacing (radial stretch is fine).  Refuse a theta stretch, the
  // Athena++ ratio keys, and any non-uniform built xx2f / xx3f (cs panel coordinates
  // are uniform by construction).
  {
    if (!vlat_rat_keys.empty()) {
      VlatFatal("needs a laterally UNIFORM mesh; this input sets" + vlat_rat_keys +
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
        VlatFatal(std::string("needs a laterally UNIFORM mesh; the built grid has ") + buf +
                  " (the SC upwind interpolation is bilinear in (theta, phi) index "
                  "space; only radial stretches are supported)");
      }
    }
  }
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
  // follow <mesh>/cs_seam_blockend as pbval_u does
  pbval_vs->cs_blockend_ext_ = pbval_u->cs_blockend_ext_;
  VetGdInit();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatExchange
//! \brief cell-centred ghost exchange of a (m, n, k, j, i) array

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
//! \brief ghosts of vlat_cs, then the ray cut (after VetLatSweep's source kernel)

void RadiationVet::VetLatCsFinish() {
  VetLatExchange(vlat_cs, vlat_cs_c, pbval_vs);
  CleanRayCut();
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetLatSweep
//! \brief source of a vet_gd build: ln chi, ln S of every active cell, their ghosts and
//! the ray cut

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
  // vet_source_noesrc: the physical E^n = u0
  const bool noes = VsrcNoes();   // rad_vet.hpp
  auto ur_ = u0;
  const int vsf = vscat ? vscat_form : 0;   // vet_scatter
  const Real ar = arad, efl = e_floor;
  const bool cvx = cln_srccvx;   // clean_source_convex: kappa_P capped at chi
  par_for_lb("vet_vlat_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real chi = fmax(iw_(m,VET_IW_KT,k,j,i), 1.0e-300);
    Real e = fmax(noes ? ur_(m,VET_E,k,j,i) : iw_(m,VET_IW_EN,k,j,i), efl);
    Real s = e;
    if (thermal) {
      Real tg = iw_(m,VET_IW_TP,k,j,i);
      Real t2 = tg*tg;
      Real eth = fmin(opac_(m,VET_OP_P,k,j,i)/chi, 1.0);
      const Real rks = (vsf != 0) ? opac_(m,VET_OP_S,k,j,i) : 0.0;
      if (vsf != 0) {
        const Real rkp = cvx ? fmin(opac_(m,VET_OP_P,k,j,i), chi) :
            opac_(m,VET_OP_P,k,j,i);
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
