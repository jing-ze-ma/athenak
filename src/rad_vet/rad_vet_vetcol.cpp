//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_vetcol.cpp
//! \brief <rad_vet>/closure = vet_col under scheme = clean: the set-up the vet_gd closure
//! reads (VetColInit: the checks, tau_ten, the outer Marshak q vcol_q, the spherical
//! Marshak coefficients) and the end-of-run report.  The per-column 1-D formal solution
//! itself (the old production closure) is not in this build.

#include <algorithm>
#include <cmath>
#include <cstdio>
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
//! \brief the checks, the spherical Marshak coefficients, vcol_q and tau_ten

void RadiationVet::VetColInit() {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto &mindcs = pm->mesh_indcs;
  if (!trans_on) {VcolFatal("needs <rad_vet>/transport = implicit on a MULTI-D mesh");}
  if (pm->multilevel) {VcolFatal("SMR/AMR is not supported");}
  // (the clean closure from vet_gd sets up but never builds the columns)
  if (mindcs.nx1 != indcs.nx1 && !vgd_on) {
    VcolFatal("needs ONE MeshBlock along x1 (meshblock/nx1 = mesh/nx1)");
  }
  if (ibc_x1min == M1_IBC_PERIODIC) {VcolFatal("a periodic x1 boundary has no top");}
  // STAGE CS2 (m1-cs-implicit): the cubed sphere with the RADIAL axis.  x1 is r on every
  // panel and one block spans the whole column, so the column is the radial line of a
  // panel cell and the p-ray tables of block 0 serve every block, as on sp.  The tensor
  // is uniaxial about r_hat; its tangential part (1-f_K)/2 is an isotropic pressure,
  // which the cs transverse rows read as the scalar D_tt (M1DDiag, CS1).  The flux axis
  // needs the cell F direction on the skewed panel basis: refused.
  if (pm->use_cubed_sphere && vcol_axis_flux) {
    VcolFatal("vet_col_axis = flux is not supported on the cubed sphere");
  }
  if (indcs.nx1 < 2) {VcolFatal("needs at least 2 cells along x1");}
  if (vcol_nc < 1 || vcol_np < 1 || vcol_nmu < 1 || vcol_every < 1) {
    VcolFatal("vet_col_ncore, vet_col_nsub, vet_col_nmu, vet_col_every must be >= 1");
  }
  vcol_sph = sph_geom;
  // vet_col_reflect_top: only on sp with a reflecting outer x1 of the M1 solve
  vcol_rtop = vcol_rtop && vcol_sph && (ibc_x1max == M1_IBC_REFLECT);
  vcol_sqf = vcol_sqf && vcol_sph && vcol_sq && impl_mface_lin;
  const int ie = indcs.ie;
  // (the per-column formal solution's quadrature tables went with VetColBuild: under
  // scheme = clean the closure is vet_gd's, which reads tau_ten, vcol_q and the
  // spherical Marshak coefficients set here)
  if (vcol_sph) {
    auto x1v_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->x1v);
    auto x1f_h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                     pmy_pack->pcoord->xx1f);
    const double rtop = x1f_h(0, ie + 1);
    M1SphMarshakCoef(x1v_h(0, ie), x1v_h(0, ie - 1), rtop, vcol_jca, vcol_jcb);
    vcol_jcap = 2.0*SQR(x1v_h(0, ie)/rtop);
  }
  const int nmb = pmy_pack->nmb_thispack;
  int c1 = indcs.nx1 + 2*(indcs.ng);
  int c2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int c3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  if (vcol_sq) {
    Kokkos::realloc(vcol_q, nmb, c3, c2);
    Kokkos::deep_copy(vcol_q, marshak_q);
    if (ibc_x1max != M1_IBC_MARSHAK && global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_col_surface_q: the outer x1 boundary is not marshak; "
                << "the q is built but unused" << std::endl;
    }
  }
  // vet_col_lat (rad_vet_vetlat.cpp) and the cubed sphere: M1_TT_NLAT more slots, the
  // lateral correction D_ab (cs: in the basis r, xi, eta).  Elsewhere 4 slots.
  Kokkos::realloc(tau_ten, nmb, (vlat_on || pm->use_cubed_sphere) ?
                  (M1_TT_LAT0 + M1_TT_NLAT) : 4, c3, c2, c1);
  Kokkos::deep_copy(tau_ten, 0.0);
  tau_ready = true;
  if (vlat_on) {VetLatInit();}
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetColReport
//! \brief one cost line at the end of the run (rank 0)

void RadiationVet::VetColReport() {
  if (global_variable::my_rank != 0) return;
  std::cout << "<rad_vet> closure = vet_col: " << vcol_ncall << " tensor builds ("
            << vcol_nskip << " skipped), " << vcol_time << " s (rank 0, fenced), "
            << ((vcol_ncall > 0.0) ? (1.0e3*vcol_time/vcol_ncall) : 0.0)
            << " ms per build" << std::endl;
  if (vlat_on) {
    std::cout << "<rad_vet> vet_col_lat: " << vlat_ncall << " sweeps, " << vlat_time
              << " s (rank 0, fenced), "
              << ((vlat_ncall > 0.0) ? (1.0e3*vlat_time/vlat_ncall) : 0.0)
              << " ms per sweep; lateral reads clamped to the ghost band (rank 0, "
              << "all sweeps): " << vlat_nclamp << std::endl;
    if (vgd_on) {
      std::cout << "<rad_vet> vet_gd: rank 0 fenced seconds: source " << vgd_tsrc
                << ", sweeps " << vgd_tswp << " (of which per-shell halo+walls "
                << vgd_thalo << ", of which MPI wait " << vgd_tmpi
                << "), tables (rotation) " << vgd_ttab
                << ", moments " << vgd_tmom << "; in-solve rebuilds " << vgd_nrb << " ("
                << vgd_trb << " s)" << std::endl;
      std::cout << "<rad_vet> vet_gd halo (rank 0): exchanges " << vgd_nexch
                << ", bytes sent "
                << vgd_nbyte << " (" << ((vgd_nexch > 0.0) ? vgd_nbyte/vgd_nexch : 0.0)
                << " per exchange), post " << vgd_tpost << " s, partners "
                << vgd_prk.size() << "; band-clamped lateral reads (ALL ranks, all "
                << "sweeps) "
                << vgd_nclamp_all << std::endl;
      if (vgd_csg || vgd_gbc) {   // (+ the sp global band)
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
  }
}

} // namespace radvet
