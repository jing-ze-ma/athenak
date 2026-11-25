//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_sph.cpp
//! \brief Option checks for implicit VET on curvilinear meshes (spherical-polar wedge or
//! full sphere, cubed sphere) and the cubed-sphere seam set-up (CubedS1Init).
//!
//! The geometry lives in the implicit kernels: the E row takes dt A_f/V_i per face and
//! each face equation the centre-to-centre distance dxface (dr, r dtheta, r sin dphi).
//! F is stored in physical orthonormal components (r, theta, phi), as the momentum.

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"

namespace radvet {
//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SphericalS1Check
//! \brief fatal on any option the sp / cs implicit operator does not support; sets
//! impl_offdiag


void RadiationVet::SphericalS1Check(ParameterInput *pin) {
  std::string why;
  auto *pm = pmy_pack->pmesh;
  // vet_sc is refused
  if (vet_sc) {why += " closure = vet_sc;";}
  if (pm->multilevel) {why += " SMR/AMR;";}
  if (pm->mesh_indcs.nx1 != pm->mb_indcs.nx1) {
    why += " more than one MeshBlock along x1 (meshblock/nx1 must equal mesh/nx1);";
  }
  {
  // berthon | blend on the radial faces (ImplicitInit refuses them on cs)
    if (impl_flux != VET_IFLUX_CENTRAL && impl_flux != VET_IFLUX_BERTHON &&
        impl_flux != VET_IFLUX_BLEND) {
      why += " implicit_flux = ap_hll;";
    }
    // plm_dc corrects the berthon part: needs berthon | blend
    if (impl_recon != VET_IRECON_DC && impl_flux == VET_IFLUX_CENTRAL) {
      why += " implicit_recon != dc;";
    }
    // lpr allowed on the wedge, lp not; cs refuses both (skewed rows ignore thx2/thx3)
    if (impl_tlim == VET_TLIM_LP || (impl_tlim != VET_TLIM_NONE && cs_geom)) {
      why += " implicit_trans_limit != none (lpr only, wedge only);";
    }
    if (pin->DoesParameterExist("rad_vet","implicit_offdiag")) {
      std::string s = pin->GetString("rad_vet","implicit_offdiag");
      // `operator` (19-point stencil, od cache, uniform-dx VetOffDiv) is not converted
      bool ok = (s.compare("auto") == 0 || s.compare("none") == 0 ||
                 s.compare("lagged") == 0);
      if (!ok) {
        why += " implicit_offdiag = " + s + ";";
      }
    }
  }
  // cubed sphere: no lagged off-diagonal
  if (cs_geom) {
    std::string wcs;
    if (pin->DoesParameterExist("rad_vet","implicit_offdiag") &&
        pin->GetString("rad_vet","implicit_offdiag").compare("lagged") == 0) {
      wcs += " implicit_offdiag = lagged;";
    }
    // seam geometry needs the equiangular panel [-1, 1]^2 with nx2 = nx3, and the
    // one-sided tangential derivatives at a seam need >= 3 cells per block
    if (pm->mesh_indcs.nx2 != pm->mesh_indcs.nx3 ||
        pm->mesh_size.x2min != -1.0 || pm->mesh_size.x2max != 1.0 ||
        pm->mesh_size.x3min != -1.0 || pm->mesh_size.x3max != 1.0) {
      wcs += " a panel other than x2, x3 in [-1, 1] with mesh nx2 = nx3;";
    }
    if (pm->mb_indcs.nx2 < 3 || pm->mb_indcs.nx3 < 3) {
      wcs += " meshblock nx2 or nx3 < 3;";
    }
    if (!wcs.empty()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_vet> on the cubed sphere supports only "
        << "transport = implicit, closure = vet_col | none, "
        << "implicit_offdiag none, time_scheme = be | hesdirk2, implicit_solver "
        << "= bicgstab, implicit_precond = rbgs_fwd | mg (and the spherical-polar "
        << "list below); this input has:"
        << wcs << why << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // full sphere: a polar face has zero area and must be reflecting (F_theta = 0);
  // lagged off-diagonal terms and trans_limit are refused there
  if (!cs_geom && pm->use_polar_boundary) {
    const bool p2lo = (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar);
    const bool p2hi = (pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar);
    if ((p2lo && ibc_x2min != VET_IBC_REFLECT) || (p2hi && ibc_x2max != VET_IBC_REFLECT))
        {
      why += " implicit_bc_x2min|max != auto | reflect at a polar face;";
    }
    if (pin->DoesParameterExist("rad_vet","implicit_offdiag") &&
        pin->GetString("rad_vet","implicit_offdiag").compare("lagged") == 0) {
      why += " implicit_offdiag = lagged on the full sphere;";
    }
    if (impl_tlim != VET_TLIM_NONE) {
      why += " implicit_trans_limit != none on the full sphere;";
    }
  }
  if (!why.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> on a spherical-polar mesh supports "
      << "only transport = implicit, closure = vet_col | none, "
      << "implicit_offdiag = auto | none | lagged, "
      << "time_scheme = be | hesdirk2, one MeshBlock "
      << "along x1, no SMR, implicit_flux = central | berthon | blend, "
      << "implicit_recon = dc, implicit_trans_limit = none; this input has:" << why
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  {
    std::string s = pin->GetOrAddString("rad_vet","implicit_offdiag","auto");
    impl_offdiag = (s.compare("lagged") == 0) ? VET_OD_LAGGED : VET_OD_NONE;
    od_now = impl_offdiag;
  }
  if (global_variable::my_rank == 0 && cs_geom) {
    std::cout << "<rad_vet>: cubed sphere (covariant cell F and gas coupling; skewed "
              << "transverse faces with the lagged cross term, mirror-pair seam faces "
              << "without the resample; implicit with areas/volumes/face distances from "
              << "Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  ", stretched radial grid)" : ")") << std::endl;
  } else if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: spherical-polar "
              << (pm->use_polar_boundary ? "FULL SPHERE [polar x2 faces reflecting at "
                                           "zero area] " : "wedge ")
              << (pm->use_grid_stretch_theta ? "[stretched theta] " : "")
              << "("
              << "implicit with areas/volumes/face distances from Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  " (stretched radial grid)" : "") << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::CubedS1Init
//! \brief cubed sphere, after ImplicitInit: open the panel seams of the implicit
//! operator, tabulate the seam geometry (VetCsSeamGeom, incl. ghost rows), effective
//! transverse distances/areas and the face and cell-centre normals

void RadiationVet::CubedS1Init() {
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ng = indcs.ng;
  const int n2 = indcs.nx2 + 2*ng, n3 = indcs.nx3 + 2*ng;
  const int npan = pm->mesh_indcs.nx2;   // cells per panel edge (= nx3, checked)
  auto &bc = pmy_pack->pmb->mb_bcs;
  auto &msz = pmy_pack->pmb->mb_size;
  vetbcs = DualArray2D<BoundaryFlag>("vetbcs", nmb, 6);
  cs_seam = DualArray2D<int>("vetcsseam", nmb, 4);
  csg2 = DvceArray4D<Real>("vetcsg2", nmb, 4, n3, 2);
  csg3 = DvceArray4D<Real>("vetcsg3", nmb, 4, n2, 2);
  auto g2h = Kokkos::create_mirror_view(csg2);
  auto g3h = Kokkos::create_mirror_view(csg3);
  Kokkos::deep_copy(g2h, 0.0);
  Kokkos::deep_copy(g3h, 0.0);
  for (int m = 0; m < nmb; ++m) {
    for (int f = 0; f < 6; ++f) {
      BoundaryFlag b = bc.h_view(m,f);
      vetbcs.h_view(m,f) = (b == BoundaryFlag::panel) ? BoundaryFlag::block : b;
    }
    for (int f = 0; f < 4; ++f) {
      cs_seam.h_view(m,f) = (bc.h_view(m,2+f) == BoundaryFlag::panel) ? 1 : 0;
    }
    // the first active cell's index along the panel edge
    const int j0 = static_cast<int>(std::lround(0.5*(msz.h_view(m).x2min + 1.0)*npan));
    const int k0 = static_cast<int>(std::lround(0.5*(msz.h_view(m).x3min + 1.0)*npan));
    for (int sd = 0; sd < 2; ++sd) {
      for (int k = 0; k < n3; ++k) {
        const int q = k - indcs.ks + k0;
        if (q < 0 || q >= npan) continue;
        Real a, b, c, d;
        VetCsSeamGeom(q, npan, a, b, c, d);
        g2h(m,0,k,sd) = a; g2h(m,1,k,sd) = b; g2h(m,2,k,sd) = c; g2h(m,3,k,sd) = d;
      }
      for (int j = 0; j < n2; ++j) {
        const int q = j - indcs.js + j0;
        if (q < 0 || q >= npan) continue;
        Real a, b, c, d;
        VetCsSeamGeom(q, npan, a, b, c, d);
        g3h(m,0,j,sd) = a; g3h(m,1,j,sd) = b; g3h(m,2,j,sd) = c; g3h(m,3,j,sd) = d;
      }
    }
  }
  vetbcs.modify_host();
  vetbcs.sync_device();
  cs_seam.modify_host();
  cs_seam.sync_device();
  Kokkos::deep_copy(csg2, g2h);
  Kokkos::deep_copy(csg3, g3h);
  // effective transverse distances and areas (rad_vet.hpp)
  {
    auto *pc = pmy_pack->pcoord;
    const int n1 = indcs.nx1 + 2*ng;
    cs_dxf_eff = DvceFaceFld4D<Real>("vetcsdxf", nmb, n3, n2, n1);
    cs_area_eff = DvceFaceFld4D<Real>("vetcsarea", nmb, n3, n2, n1);
    Kokkos::deep_copy(cs_dxf_eff.x1f, pc->dxface.x1f);
    Kokkos::deep_copy(cs_dxf_eff.x2f, pc->dxface.x2f);
    Kokkos::deep_copy(cs_dxf_eff.x3f, pc->dxface.x3f);
    Kokkos::deep_copy(cs_area_eff.x1f, pc->area.x1f);
    Kokkos::deep_copy(cs_area_eff.x2f, pc->area.x2f);
    Kokkos::deep_copy(cs_area_eff.x3f, pc->area.x3f);
    auto dx2 = cs_dxf_eff.x2f, dx3 = cs_dxf_eff.x3f;
    auto a2 = cs_area_eff.x2f, a3 = cs_area_eff.x3f;
    auto sn2 = pc->sin_face_xi, sn3 = pc->sin_face_eta;
    auto x1v = pc->x1v, x1f = pc->xx1f;
    auto g2 = csg2, g3 = csg3;
    auto cse = cs_seam.d_view;
    const int js = indcs.js, je = indcs.je, ks = indcs.ks, ke = indcs.ke;
    const int is = indcs.is, ie = indcs.ie;
    par_for("vet_cs_eff2", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool sm = (j == js && cse(m,0) != 0) || (j == je+1 && cse(m,1) != 0);
      const int sd = (j == js) ? 0 : 1;
      if (sm) {
        dx2(m,k,j,i) = x1v(m,i)*g2(m,1,k,sd);
        a2(m,k,j,i) = 0.5*(x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i))*g2(m,0,k,sd);
      } else {
        dx2(m,k,j,i) = dx2(m,k,j,i)*sn2(m,k,j);
      }
    });
    par_for("vet_cs_eff3", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool sm = (k == ks && cse(m,2) != 0) || (k == ke+1 && cse(m,3) != 0);
      const int sd = (k == ks) ? 0 : 1;
      if (sm) {
        dx3(m,k,j,i) = x1v(m,i)*g3(m,1,j,sd);
        a3(m,k,j,i) = 0.5*(x1f(m,i+1)*x1f(m,i+1) - x1f(m,i)*x1f(m,i))*g3(m,0,j,sd);
      } else {
        dx3(m,k,j,i) = dx3(m,k,j,i)*sn3(m,k,j);
      }
    });
  }
  // planar face normals nhat_xi = (a - tan(xi_f) q)/sqrt(1 + tan^2), likewise eta; the
  // face plane holds the origin, so the normal is constant along the column
  {
    Kokkos::realloc(cs_fn2, nmb, n2 + 1, 3);
    Kokkos::realloc(cs_fn3, nmb, n3 + 1, 3);
    auto f2h = Kokkos::create_mirror_view(cs_fn2);
    auto f3h = Kokkos::create_mirror_view(cs_fn3);
    auto &pan = pmy_pack->pmb->mb_panel;
    for (int m = 0; m < nmb; ++m) {
      Real fa[3], fb[3], fq[3];
      cubed_sphere::PanelFrame(pan.h_view(m), fa, fb, fq);
      for (int j = 0; j <= n2; ++j) {
        const Real x = std::tan(0.25*M_PI*LeftEdgeX(j - indcs.js, indcs.nx2,
                                                   msz.h_view(m).x2min,
                                                   msz.h_view(m).x2max));
        const Real nn = 1.0/std::sqrt(1.0 + x*x);
        for (int c = 0; c < 3; ++c) {f2h(m,j,c) = (fa[c] - x*fq[c])*nn;}
      }
      for (int k = 0; k <= n3; ++k) {
        const Real y = std::tan(0.25*M_PI*LeftEdgeX(k - indcs.ks, indcs.nx3,
                                                   msz.h_view(m).x3min,
                                                   msz.h_view(m).x3max));
        const Real nn = 1.0/std::sqrt(1.0 + y*y);
        for (int c = 0; c < 3; ++c) {f3h(m,k,c) = (fb[c] - y*fq[c])*nn;}
      }
    }
    Kokkos::deep_copy(cs_fn2, f2h);
    Kokkos::deep_copy(cs_fn3, f3h);
    // cell-centre normals (same formula at xi_c, eta_c)
    Kokkos::realloc(cs_cn2, nmb, n2, 3);
    Kokkos::realloc(cs_cn3, nmb, n3, 3);
    auto c2h = Kokkos::create_mirror_view(cs_cn2);
    auto c3h = Kokkos::create_mirror_view(cs_cn3);
    for (int m = 0; m < nmb; ++m) {
      Real fa[3], fb[3], fq[3];
      cubed_sphere::PanelFrame(pan.h_view(m), fa, fb, fq);
      for (int j = 0; j < n2; ++j) {
        const Real x = std::tan(0.25*M_PI*CellCenterX(j - indcs.js, indcs.nx2,
                                                     msz.h_view(m).x2min,
                                                     msz.h_view(m).x2max));
        const Real nn = 1.0/std::sqrt(1.0 + x*x);
        for (int c = 0; c < 3; ++c) {c2h(m,j,c) = (fa[c] - x*fq[c])*nn;}
      }
      for (int k = 0; k < n3; ++k) {
        const Real y = std::tan(0.25*M_PI*CellCenterX(k - indcs.ks, indcs.nx3,
                                                     msz.h_view(m).x3min,
                                                     msz.h_view(m).x3max));
        const Real nn = 1.0/std::sqrt(1.0 + y*y);
        for (int c = 0; c < 3; ++c) {c3h(m,k,c) = (fb[c] - y*fq[c])*nn;}
      }
    }
    Kokkos::deep_copy(cs_cn2, c2h);
    Kokkos::deep_copy(cs_cn3, c3h);
  }
  // cs_seam_avg (read only when named): seam average of the stored face state
  cs_seam_avg_on = true;
  if (pin_cs_ != nullptr && pin_cs_->DoesParameterExist("rad_vet","cs_seam_avg")) {
    cs_seam_avg_on = pin_cs_->GetBoolean("rad_vet","cs_seam_avg");
  }
}

} // namespace radvet
