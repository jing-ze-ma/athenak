//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_sph.cpp
//! \brief STAGE S1 of docs/dev/rad_m1_curvilinear_design.md (branch m1-curv-design):
//! implicit M1 on a spherical-polar WEDGE (tests_m1/runs_5a_sp_s1).
//!
//! What S1 supports, and nothing else (every other sp configuration is a startup fatal,
//! and the cubed sphere stays refused in the constructor):
//!   * mesh/use_spherical_polar without the polar boundary (the theta range must stay
//!     clear of the poles; the theta faces are periodic or reflecting like any other
//!     boundary), the radial stretches use_grid_stretch_r / _r_poly allowed;
//!   * <rad_m1>/transport = implicit, closure = eddington, time_scheme = be (m1-sph2:
//!     also hesdirk2, and implicit_vimp; tests_m1/runs_5h_sph2);
//!   * one MeshBlock along x1 per column (implicit_partition = none), no SMR/AMR;
//!   * implicit_halo_direct, implicit_halo_mpi, implicit_halo_overlap and
//!     implicit_halo_ovl_faces as on the Cartesian mesh (without a pole every ghost is
//!     the plain copy; halo_mpi and the overlap allowed since m1-sphhalo,
//!     tests_m1/runs_5i_sphhalo);
//!   * implicit_flux = central, implicit_recon = dc, implicit_trans_limit = none,
//!     no dbg_tensor; implicit_offdiag is set to none (the Eddington
//!     tensor has no off-diagonal part, so the terms it drops are identically zero).
//!
//! The geometry itself lives in the implicit kernels (rad_m1_implicit.cpp, `sph`
//! branches): the E row takes dt A_f / V_i for each face (the Coordinates areas and
//! volumes, radial stretch included) instead of dt/dx_d, and each face-flux equation
//! takes the centre-to-centre distance dxface (dr, r dtheta, r sin(theta) dphi) instead
//! of dx_d.  F is stored in PHYSICAL orthonormal components (r, theta, phi), as the hydro
//! momentum, so the radiation force is handed to the hydro unchanged.
//!
//! STAGE S2 (tests_m1/runs_5b_sp_s2) adds closure = m1 | minerbo | kershaw on the same
//! wedge (sph_q).  With P = p I + q n n the radial face takes the INTEGRATING FACTOR of
//! design sect. 1.4 (M1SphDrr: D_rr -> (1-chi)/2 + (3 chi-1)/2 n_r^2 (r_c/r_f)^2 per cell
//! and face, which keeps the M-matrix and makes r^2 E = const the exact discrete free-
//! streaming state), and every face equation takes the LAGGED rest of (div P)_d at the
//! two cell centres (M1SphCurv, rad_m1_implicit.hpp): the curvature of the diagonal
//! components always (-q(1-n_r^2)/r radially, cot(P_tt-P_pp)/r on theta faces), the
//! off-diagonal terms with their spherical factors under implicit_offdiag = lagged
//! (auto = none on sp, see SphericalS1Check; `operator` stays refused).
//! Cartesian M1OffDiv values are overwritten wherever they are formed, and the od
//! cache and the 19-point stencil (operator only) never run on this mesh.
//!
//! m1-sp-order2 (tests_m1/runs_5o_sporder2): the interior operator on the wedge is
//! already SECOND order in space (centred face-normal F0 on the faces, A_f/V_i
//! divergence, exact integrating factor, face-mean curvature; implicit_enthalpy = plm
//! for the enthalpy flux, uniform and stretched r).  implicit_recon / implicit_flux /
//! implicit_trans_limit stay refused: under transport = implicit the face flux is the
//! central form on every mesh (ImplicitInit), so plm_dc has nothing to act on.  What
//! was first order is the Marshak end face (c q E of the end CELL) -> option
//! <rad_m1>/implicit_marshak_face = linear (sp only), and the vet_col formal solution
//! (trapezoid in mu, S and chi linear in the path, cell E/F at the core) -> option
//! <rad_m1>/vet_col_order2 = true.
//!
//! m1-sp-order2b (tests_m1/runs_5q_sporder2b): on the wedge both are the DEFAULT
//! (implicit_marshak_face = linear with the fixed-tensor closures eddington, vet_sc,
//! tau and vet_col; cell with the lagged m1 closures; vet_col_order2 = true), with
//! vet_col_fk_min = 0 (no 1/3 clamp of f_K), vet_col_reflect_top = true (the mirrored
//! top under a reflecting outer x1), vet_col_surface_face = true (q from the face J),
//! time2_vet_col = predict (the vet_col tensor of the hesdirk2 stages at t^{n+1}, not
//! U^n), time2_vstage = true (implicit_vimp: the stage's own velocity in the cell flux
//! and the gas work in the E row) and implicit_precond = mg_gc (implicit_mg_levels = 1)
//! where rbgs_fwd was the default (mg_gc reverted to rbgs_fwd by m1-wedge, 09-26:
//! -8.7 % on the rad-hydro wedge with moving gas).  A restart whose file lacks a key
//! keeps its old value (the resolved values are echoed); Cartesian meshes are unchanged.

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
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

void RadiationM1::SphericalS1Check(ParameterInput *pin) {
  std::string why;
  auto *pm = pmy_pack->pmesh;
  if (transport != M1_TRANSPORT_IMPLICIT) {why += " transport != implicit;";}
  // STAGE S2: the uniaxial chi(f) closures (m1, minerbo, kershaw) as well; vet_sc and
  // the tau closure (plane-parallel column depth on index space) stay refused
  // STAGE S5: closure = vet_col (the per-column spherical formal solution) is allowed
  if (vet_sc || (tau_closure && !vet_col)) {why += " closure = vet_sc | tau;";}
  if (pm->multilevel) {why += " SMR/AMR;";}
  if (pm->mesh_indcs.nx1 != pm->mb_indcs.nx1) {
    why += " more than one MeshBlock along x1 (meshblock/nx1 must equal mesh/nx1);";
  }
  if (transport == M1_TRANSPORT_IMPLICIT) {
    // sp-blend-1008: berthon | blend on the radial faces (the wedge only; ImplicitInit
    // refuses them on the cubed sphere and the tau-only weight everywhere multi-D)
    if (impl_flux != M1_IFLUX_CENTRAL && impl_flux != M1_IFLUX_BERTHON &&
        impl_flux != M1_IFLUX_BLEND) {
      why += " implicit_flux = ap_hll;";
    }
    // rad-beam-1008: the plm deferred correction of the berthon part (plm_dc) on the
    // radial faces, with berthon | blend only (central has no upwind part to correct)
    if (impl_recon != M1_IRECON_DC && impl_flux == M1_IFLUX_CENTRAL) {
      why += " implicit_recon != dc;";
    }
    // m1lpr-1019: lpr (the realizability bound) is allowed on the wedge, lp is not; the
    // cubed sphere refuses both (its skewed rows do not read thx2/thx3)
    if (impl_tlim == M1_TLIM_LP || (impl_tlim != M1_TLIM_NONE && cs_geom)) {
      why += " implicit_trans_limit != none (lpr only, wedge only);";
    }
    // m1-sph2 (tests_m1/runs_5h_sph2): time_scheme = hesdirk2 and implicit_vimp are
    // allowed on the wedge (the stage solves' old vector carries no geometry; the vimp
    // rows take the sp areas, volumes and face distances in ImplicitVimpBuild)
    if (dbg_tensor != 0) {why += " dbg_tensor != none;";}
    if (pin->DoesParameterExist("rad_m1","implicit_offdiag")) {
      std::string s = pin->GetString("rad_m1","implicit_offdiag");
      // S2: `lagged` for the non-Eddington closures; `operator` (the 19-point stencil,
      // the od cache and M1OffDiv on the uniform mb_size) is not converted
      bool ok = (s.compare("auto") == 0 || s.compare("none") == 0 ||
                 (!eddington && s.compare("lagged") == 0));
      if (!ok) {
        why += " implicit_offdiag = " + s + ";";
      }
    }
  }
  // STAGE CS0 (m1-cs-implicit, rt_design_1003/CS_IMPLICIT_VET.md sect. 4): on the cubed
  // sphere only closure = eddington (isotropic tangential pressure, so no transverse
  // curvature and no Christoffel face sum), time_scheme = be, the fused bicgstab with
  // implicit_precond = rbgs_fwd (mg has no seam-aware coarse grid), no implicit_vimp,
  // no implicit_halo_mpi; the transverse faces are the plain two-point form (stage CS1
  // adds the skew 1/sin factor and the lagged cross term)
  if (cs_geom) {
    std::string wcs;
    // STAGE CS2: closure = vet_col with the radial axis as well (uniaxial about r_hat:
    // no transverse curvature, the radial integrating factor of S2 as on sp)
    if ((!eddington && !vet_col) || vet_sc || (tau_closure && !vet_col)) {
      wcs += " closure != eddington | vet_col;";
    }
    if (vet_col && vcol_axis_flux) {wcs += " vet_col_axis = flux;";}
    if (pin->DoesParameterExist("rad_m1","implicit_offdiag") &&
        pin->GetString("rad_m1","implicit_offdiag").compare("lagged") == 0) {
      wcs += " implicit_offdiag = lagged;";
    }
    // STAGE CS3 (C9): time_scheme = hesdirk2 as well.  The stage solves' old vector and
    // slopes are cell- and face-local (rad_m1_time2.cpp), so they carry no geometry; the
    // seam face state is averaged at the start of each solve (CubedSeamFaceAverage).
    // implicit_vimp too (its cs form: ImplicitVimpBuild, the face-normal work row)
    if (impl_solver != M1_ISOLV_BICGSTAB) {wcs += " implicit_solver != bicgstab;";}
    // cs-port-1016 (CS_VETGD.md S2): implicit_precond = mg as well.  mg is algebraic and
    // block-local (x1 line solves, one red-black transverse sweep, Galerkin 2x2
    // aggregates inside a MeshBlock, rad_m1_precond.cpp): it reads no geometry, and its
    // only communication, the z0 halo of the fine residual, is ImplicitKrylovHalo, i.e.
    // the cs seam (no-resample, mirror-pair) exchange of the Krylov vector itself.
    // mg_gc / mg_gf stay refused (band / Fourier coarse spaces of a periodic x2/x3 mesh;
    // the cs analogue is mg_ylm, stage 3 list).  The cs DEFAULT stays rbgs_fwd.
    if (impl_prec != 2 && impl_prec != 3) {
      wcs += " implicit_precond != rbgs_fwd | mg;";
    }
    if (impl_prec_fb >= 0) {   // precond-fb-1023: rbgs / line are refused on cs
      wcs += " implicit_precond_fallback != none;";
    }
    // cs-s4-1017 (USER 10-18: production sp keys must run on cs): implicit_halo_mpi
    // (and its overlap) is a communication layout, not a method; on cs it falls back to
    // the ordinary seam-aware exchange (the same numbers), as implicit_halo_direct does
    if (impl_halo_mpi) {
      impl_halo_mpi = false;
      impl_halo_ovl = false;
      impl_ovl_faces = false;
      impl_halo_ipc = false;
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> cubed sphere: implicit_halo_mpi (+ overlap) falls back to "
                  << "the ordinary seam-aware exchange (the seam halo_mpi port is open)"
                  << std::endl;
      }
    }
    // CS1: the canonical seam geometry needs the equiangular panel of -1..1 with the same
    // cell count along both tangential axes, and the one-sided tangential derivatives
    // next to a seam three cells of the block
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
        << std::endl << "<rad_m1> on the cubed sphere (stages CS0-CS3) supports only "
        << "transport = implicit, closure = eddington | vet_col (radial axis), "
        << "implicit_offdiag none, time_scheme = be | hesdirk2, implicit_solver "
        << "= bicgstab, implicit_precond = rbgs_fwd | mg (and the spherical-polar "
        << "list below); this input has:"
        << wcs << why << std::endl
        << "See /viper/ptmp2/jinma/rt_design_1003/CS_IMPLICIT_VET.md sect. 4."
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  // fullsph-1027 S1: the FULL sphere (polar x2 faces) and/or a theta stretch.  A polar
  // face is a physical face of every implicit kernel (flag != block, periodic) with area 0,
  // so it must be REFLECTING (F_theta = 0 there; auto resolves to reflect on a non-periodic
  // face): marshak / outflow / flux would put a boundary flux on the axis.  Out of S1
  // scope on that mesh: implicit_offdiag = lagged (one-sided off-diagonal differences at
  // the axis untested), implicit_trans_limit = lpr; vet_col_lat / vet_gd (S2) are refused
  // in VetLatInit, vet_sc / tau above.
  if (!cs_geom && pm->use_polar_boundary) {
    const bool p2lo = (pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar);
    const bool p2hi = (pm->mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar);
    if ((p2lo && ibc_x2min != M1_IBC_REFLECT) || (p2hi && ibc_x2max != M1_IBC_REFLECT)) {
      why += " implicit_bc_x2min|max != auto | reflect at a polar face;";
    }
    if (pin->DoesParameterExist("rad_m1","implicit_offdiag") &&
        pin->GetString("rad_m1","implicit_offdiag").compare("lagged") == 0) {
      why += " implicit_offdiag = lagged on the full sphere;";
    }
    if (impl_tlim != M1_TLIM_NONE) {why += " implicit_trans_limit != none on the full sphere;";}
  }
  if (!why.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> on a spherical-polar mesh (stages S1, S2) supports "
      << "only transport = implicit, closure = eddington | m1 | minerbo | kershaw | "
      << "vet_col, "
      << "implicit_offdiag = auto | none (| lagged for m1/minerbo/kershaw), "
      << "time_scheme = be | hesdirk2, one MeshBlock "
      << "along x1, no SMR, implicit_flux = central | berthon | blend, "
      << "implicit_recon = dc, implicit_trans_limit = none, no "
      << "dbg_tensor; this input has:" << why << std::endl
      << "See tests_m1/runs_5a_sp_s1/README.md and tests_m1/runs_5b_sp_s2/README.md."
      << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // the Eddington tensor is diagonal: the off-diagonal terms are identically zero, so
  // dropping them changes no number and keeps the uniform-dx M1OffDiv off this mesh
  impl_offdiag = M1_OD_NONE;
  od_now = M1_OD_NONE;
  // STAGE S2: a chi(f) closure.  The curvature of the diagonal components (M1SphCurv)
  // is always on; the off-diagonal terms are dropped (auto = none, the M-matrix form)
  // or, with implicit_offdiag = lagged, the sp M1SphCurv form lagged into the face
  // equations.  Lagged off-diagonal terms have no fixed point in an optically thin
  // region at c dt >> dx (tests_m1/runs_5b_sp_s2: the free-streaming source diverges
  // with them on the wedge, and so does the Cartesian lagged form with reflecting x2
  // walls), which is why `none` is the default here.
  sph_q = !eddington;
  if (sph_q) {
    std::string s = pin->GetOrAddString("rad_m1","implicit_offdiag","auto");
    impl_offdiag = (s.compare("lagged") == 0) ? M1_OD_LAGGED : M1_OD_NONE;
    od_now = impl_offdiag;
  }
  if (global_variable::my_rank == 0 && cs_geom) {
    std::cout << "<rad_m1>: cubed sphere (stage CS2: " << (vet_col ? "vet_col, radial "
              "axis" : "eddington") << "; covariant cell F and gas coupling; skewed "
              << "transverse faces with the lagged cross term, mirror-pair seam faces "
              << "without the resample; implicit with areas/volumes/face distances from "
              << "Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  ", stretched radial grid)" : ")") << std::endl;
  } else if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: spherical-polar "
              << (pm->use_polar_boundary ? "FULL SPHERE [polar x2 faces reflecting at "
                                           "zero area] " : "wedge ")
              << (pm->use_grid_stretch_theta ? "[stretched theta] " : "")
              << "("
              << (sph_q ? "stage S2: chi(f) closure, radial integrating factor, "
                          "lagged curvature, offdiag " : "stage S1: ")
              << (sph_q ? ((impl_offdiag == M1_OD_NONE) ? "none; " : "lagged; ") : "")
              << "implicit with areas/volumes/face distances from Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  " (stretched radial grid)" : "") << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CubedS1Init
//! \brief STAGE CS1 (cubed sphere, after ImplicitInit): open the panel seams of the
//! implicit operator (m1bcs), record which block faces are seams, tabulate the canonical
//! seam geometry (M1CsSeamGeom) of every x2/x3 seam face incl. the ghost rows along the
//! seam, and switch the implicit scratch exchanges to the no-resample (mirror) halo.

void RadiationM1::CubedS1Init() {
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ng = indcs.ng;
  const int n2 = indcs.nx2 + 2*ng, n3 = indcs.nx3 + 2*ng;
  const int npan = pm->mesh_indcs.nx2;   // cells per panel edge (= nx3, checked)
  auto &bc = pmy_pack->pmb->mb_bcs;
  auto &msz = pmy_pack->pmb->mb_size;
  m1bcs = DualArray2D<BoundaryFlag>("m1bcs", nmb, 6);
  cs_seam = DualArray2D<int>("m1csseam", nmb, 4);
  csg2 = DvceArray4D<Real>("m1csg2", nmb, 4, n3, 2);
  csg3 = DvceArray4D<Real>("m1csg3", nmb, 4, n2, 2);
  auto g2h = Kokkos::create_mirror_view(csg2);
  auto g3h = Kokkos::create_mirror_view(csg3);
  Kokkos::deep_copy(g2h, 0.0);
  Kokkos::deep_copy(g3h, 0.0);
  for (int m = 0; m < nmb; ++m) {
    for (int f = 0; f < 6; ++f) {
      BoundaryFlag b = bc.h_view(m,f);
      m1bcs.h_view(m,f) = (b == BoundaryFlag::panel) ? BoundaryFlag::block : b;
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
        M1CsSeamGeom(q, npan, a, b, c, d);
        g2h(m,0,k,sd) = a; g2h(m,1,k,sd) = b; g2h(m,2,k,sd) = c; g2h(m,3,k,sd) = d;
      }
      for (int j = 0; j < n2; ++j) {
        const int q = j - indcs.js + j0;
        if (q < 0 || q >= npan) continue;
        Real a, b, c, d;
        M1CsSeamGeom(q, npan, a, b, c, d);
        g3h(m,0,j,sd) = a; g3h(m,1,j,sd) = b; g3h(m,2,j,sd) = c; g3h(m,3,j,sd) = d;
      }
    }
  }
  m1bcs.modify_host();
  m1bcs.sync_device();
  cs_seam.modify_host();
  cs_seam.sync_device();
  Kokkos::deep_copy(csg2, g2h);
  Kokkos::deep_copy(csg3, g3h);
  // STAGE CS3: the effective transverse distances and areas (rad_m1.hpp)
  {
    auto *pc = pmy_pack->pcoord;
    const int n1 = indcs.nx1 + 2*ng;
    cs_dxf_eff = DvceFaceFld4D<Real>("m1csdxf", nmb, n3, n2, n1);
    cs_area_eff = DvceFaceFld4D<Real>("m1csarea", nmb, n3, n2, n1);
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
    par_for("m1_cs_eff2", DevExeSpace(), 0, nmb-1, ks, ke, js, je+1, is, ie,
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
    par_for("m1_cs_eff3", DevExeSpace(), 0, nmb-1, ks, ke+1, js, je, is, ie,
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
  // cs-s4-1017: the planar face normals, nhat_xi = (a - tan(xi_f) n)/sqrt(1 + tan^2) and
  // nhat_eta = (b - tan(eta_f) n)/sqrt(1 + tan^2) (gradient directions of xi, eta; the
  // face plane xi = xi_f holds the origin, so its normal is the same along the column)
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
    // the cell-centre normals (the same formula at xi_c, eta_c)
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
  // DEBUG <rad_m1>/cs_seam_avg (read only when named, default true): the C5 seam average
  // of the stored face state
  cs_seam_avg_on = true;
  if (pin_cs_ != nullptr && pin_cs_->DoesParameterExist("rad_m1","cs_seam_avg")) {
    cs_seam_avg_on = pin_cs_->GetBoolean("rad_m1","cs_seam_avg");
  }
  if (pbval_th != nullptr) {
    pbval_th->cs_noresample = true;
    // STAGE CS2: (N2,N3), (A2,A3), (V2,V3) of the work array are face-normal components
    // (cs-s4-1017: re-expressed on the destination's face normals at a seam,
    // bvals_cc.cpp)
    pbval_th->cs_perm_pairs = true;
  }
  if (pbval_tq != nullptr) {pbval_tq->cs_noresample = true;}
  if (pbval_kr != nullptr) {pbval_kr->cs_noresample = true;}
  // STAGE CS3: implicit_vimp's exchange; (DV2, DV3) are face-normal components
  if (impl_vimp && pbval_vm != nullptr) {
    pbval_vm->cs_noresample = true;
    pbval_vm->cs_perm_pairs = true;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CubedSeamFaceAverage
//! \brief STAGE CS1 (C5): make the stored transverse face state f0x2/f0x3 single valued
//! at the panel seams, each panel keeping the mean of the two outward values
//! (MeshBoundaryValuesCC::*FluxSeamCC, the hydro's seam-flux exchange, on pbval_kr).
//! Called at the start of each implicit step, before f0x* is copied to the old state.

void RadiationM1::CubedSeamFaceAverage() {
  if (!cs_geom || pbval_kr == nullptr || !cs_seam_avg_on) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = indcs.nx2 + 2*indcs.ng;
  const int n3 = indcs.nx3 + 2*indcs.ng;
  DvceFaceFld5D<Real> fl("m1_csfa", nmb, 1, n3, n2, n1);
  auto f2 = f0x2;
  auto f3 = f0x3;
  auto x2 = fl.x2f;
  auto x3 = fl.x3f;
  par_for("m1_csfa_in2", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x2(m,0,k,j,i) = f2(m,k,j,i);
  });
  par_for("m1_csfa_in3", DevExeSpace(), 0, nmb-1, 0, n3, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    x3(m,0,k,j,i) = f3(m,k,j,i);
  });
  pbval_kr->InitFluxSeamRecv(1);
  pbval_kr->PackAndSendFluxSeamCC(fl);
  while (pbval_kr->RecvAndUnpackFluxSeamCC(fl) != TaskStatus::complete) {}
  pbval_kr->ClearFluxSend();
  pbval_kr->ClearFluxRecv();
  // diagnostic: the largest change the average makes, relative to max |F0| (cs_seam_dmax)
  Real dmx = 0.0, fmx = 0.0;
  Kokkos::parallel_reduce("m1_csfa_d2", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
      DevExeSpace(), {0, 0, 0, 0}, {nmb, n3, n2 + 1, n1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &ld, Real &lf) {
    ld = fmax(ld, fabs(f2(m,k,j,i) - x2(m,0,k,j,i)));
    lf = fmax(lf, fabs(f2(m,k,j,i)));
  }, Kokkos::Max<Real>(dmx), Kokkos::Max<Real>(fmx));
  cs_seam_dmax = fmax(cs_seam_dmax, dmx/fmax(fmx, 1.0e-300));
  par_for("m1_csfa_out2", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f2(m,k,j,i) = x2(m,0,k,j,i);
  });
  par_for("m1_csfa_out3", DevExeSpace(), 0, nmb-1, 0, n3, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    f3(m,k,j,i) = x3(m,0,k,j,i);
  });
  if (cs_seam_avg_n == 0 && global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: CS1 seam average of f0x2/f0x3 active" << std::endl;
  }
  ++cs_seam_avg_n;
}


//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::CubedVecCheck
//! \brief cs-s4-1017 (TEST ONLY, cubed sphere, <rad_m1>/cs_vec_check = true read only
//! when named): the seam/vertex halo of the TANGENTIAL VECTOR PAIRS against an analytic
//! field.  v = Omega x rhat (a solid-body rotation about a generic axis, cs_vec_axis
//! "ax ay az", default (0.3,-0.5,0.8)) and the scalar s = 1 + 0.2 rhat.(0.6,0.0,0.8) are
//! written into the ACTIVE cells, the ghosts are poisoned, and the run's own exchanges fill them:
//!  (1) u0 through pbval_u (F2, F3 = covariant components on the unit panel tangents
//!      e_xi, e_eta, TransformMomentum + the along-seam quadratic resample): every ghost
//!      is compared with the exact s and v . e_a at the GHOST's own position in this
//!      block's chart;
//!  (2) the implicit work array through pbval_th (V2, V3 = FACE-NORMAL components v .
//!      nhat_xi, v . nhat_eta, no resample; cs-s4-1017: re-expressed on the destination's
//!      face normals, before it the signed axis permutation of CS2): a seam ghost is
//!      compared with (A) the exact face-normal components at the ghost position
//!      and (B) at the no-resample donor's position (the mirror image of this block's own
//!      cell at the same depth across the seam plane), separately for the SEAM-NORMAL
//!      component (the only one the cs operator reads across a seam) and the ALONG-SEAM
//!      one.
//! Per class (same panel, seam depth 1..ng, cube vertex) the largest error over the
//! mesh is printed (CSVEC lines).  Everything is restored afterwards (u0, iw): no effect
//! on the run.  Called once after pbval_u is built (rad_m1.cpp).

void RadiationM1::CubedVecCheck(ParameterInput *pin) {
  if (!cs_geom || !pin->DoesParameterExist("rad_m1","cs_vec_check") ||
      !pin->GetBoolean("rad_m1","cs_vec_check")) {return;}
  auto *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  const int nmb = pmy_pack->nmb_thispack;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nx2 = indcs.nx2, nx3 = indcs.nx3;
  const int n1 = indcs.nx1 + 2*ng, n2 = nx2 + 2*ng, n3 = nx3 + 2*ng;
  const int npan = pm->mesh_indcs.nx2;
  Real om[3] = {0.3, -0.5, 0.8};
  if (pin->DoesParameterExist("rad_m1","cs_vec_axis")) {
    std::string s = pin->GetString("rad_m1","cs_vec_axis");
    std::sscanf(s.c_str(), "%lf %lf %lf", &om[0], &om[1], &om[2]);
  }
  {
    const Real on = std::sqrt(om[0]*om[0] + om[1]*om[1] + om[2]*om[2]);
    for (int c = 0; c < 3; ++c) {om[c] /= on;}
  }
  const Real o0 = om[0], o1 = om[1], o2 = om[2];
  auto &msz = pmy_pack->pmb->mb_size;
  auto &mpan = pmy_pack->pmb->mb_panel;
  // panel-global index of the first active cell of each block
  DualArray2D<int> g0("m1_csvec_g0", nmb, 2);
  for (int m = 0; m < nmb; ++m) {
    g0.h_view(m,0) = static_cast<int>(std::lround(0.5*(msz.h_view(m).x2min + 1.0)*npan));
    g0.h_view(m,1) = static_cast<int>(std::lround(0.5*(msz.h_view(m).x3min + 1.0)*npan));
  }
  g0.modify_host();
  g0.sync_device();
  auto g0d = g0.d_view;
  auto mszd = msz.d_view;
  auto mpd = mpan.d_view;
  // the analytic field at chart point (p, xi, eta): s, (v.e1, v.e2), (v.n1, v.n2)
  auto field = KOKKOS_LAMBDA(const int p, const Real xi, const Real eta, Real &s,
                             Real &c1, Real &c2, Real &f1, Real &f2) {
    Real q[3], e1[3], e2[3];
    cubed_sphere::PanelToCart(p, xi, eta, q);
    cubed_sphere::PanelTangents(p, xi, eta, e1, e2);
    const Real v0 = o1*q[2] - o2*q[1], v1 = o2*q[0] - o0*q[2], v2 = o0*q[1] - o1*q[0];
    s = 1.0 + 0.2*(0.6*q[0] + 0.8*q[2]);
    c1 = v0*e1[0] + v1*e1[1] + v2*e1[2];
    c2 = v0*e2[0] + v1*e2[1] + v2*e2[2];
    const Real cs = e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2];
    const Real sn = sqrt(1.0 - cs*cs);
    // nhat_xi = (e1 - cs e2)/sn, nhat_eta = (e2 - cs e1)/sn
    f1 = (c1 - cs*c2)/sn;
    f2 = (c2 - cs*c1)/sn;
  };
  // chart coordinates of cell (m,k,j) of this block (ghosts: its own extended chart)
  auto chart = KOKKOS_LAMBDA(const int m, const int k, const int j, Real &xi, Real &eta) {
    xi = 0.25*M_PI*CellCenterX(j - js, nx2, mszd(m).x2min, mszd(m).x2max);
    eta = 0.25*M_PI*CellCenterX(k - ks, nx3, mszd(m).x3min, mszd(m).x3max);
  };
  // class of a ghost (k,j) with i active: -1 active, 0 same panel, d = 1..ng seam depth,
  // ng+1 vertex (out of the panel along both tangential axes); sd: 2 / 3 the seam axis
  auto gclass = KOKKOS_LAMBDA(const int m, const int k, const int j, int &sd) -> int {
    sd = 0;
    if (j >= js && j <= je && k >= ks && k <= ke) {return -1;}
    const int jg = j - js + g0d(m,0), kg = k - ks + g0d(m,1);
    const int dj = (jg < 0) ? -jg : ((jg >= npan) ? (jg - npan + 1) : 0);
    const int dk = (kg < 0) ? -kg : ((kg >= npan) ? (kg - npan + 1) : 0);
    if (dj > 0 && dk > 0) {return ng + 1;}
    if (dj > 0) {sd = 2; return dj;}
    if (dk > 0) {sd = 3; return dk;}
    return 0;
  };
  const int ncl = ng + 2;
  std::vector<Real> res(8*ncl, 0.0);
  // ---------------- (1) u0 through pbval_u
  {
    DvceArray5D<Real> ubak("m1_csvec_u", u0.extent(0), u0.extent(1), u0.extent(2),
                           u0.extent(3), u0.extent(4));
    Kokkos::deep_copy(ubak, u0);
    auto u = u0;
    par_for("m1_csvec_f", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      if (!act) {
        for (int v = 0; v < M1_NVAR; ++v) {u(m,v,k,j,i) = 1.0e20;}
        return;
      }
      Real xi, eta, s, c1, c2, f1, f2;
      chart(m, k, j, xi, eta);
      field(mpd(m), xi, eta, s, c1, c2, f1, f2);
      u(m,M1_E,k,j,i) = s;
      u(m,M1_F1,k,j,i) = 0.5*s;
      u(m,M1_F2,k,j,i) = c1;
      u(m,M1_F3,k,j,i) = c2;
    });
    while (pbval_u->InitRecv(M1_NVAR) == TaskStatus::incomplete) {}
    while (pbval_u->PackAndSendCC(u0, coarse_u0) == TaskStatus::incomplete) {}
    while (pbval_u->RecvAndUnpackCC(u0, coarse_u0) == TaskStatus::incomplete) {}
    while (pbval_u->ClearSend() == TaskStatus::incomplete) {}
    while (pbval_u->ClearRecv() == TaskStatus::incomplete) {}
    for (int c = 0; c < ncl; ++c) {
      Real es = 0.0, ev = 0.0, ecnt = 0.0, ep = 0.0;
      const int cc = c;
      Kokkos::parallel_reduce("m1_csvec_r1", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExeSpace(), {0, 0, 0, is}, {nmb, n3, n2, ie+1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a1,
                    Real &a2, Real &a3, Real &a4) {
        int sd;
        if (gclass(m, k, j, sd) != cc) return;
        Real xi, eta, s, c1, c2, f1, f2;
        chart(m, k, j, xi, eta);
        field(mpd(m), xi, eta, s, c1, c2, f1, f2);
        a1 = fmax(a1, fabs(u(m,M1_E,k,j,i) - s));
        a2 = fmax(a2, sqrt((u(m,M1_F2,k,j,i) - c1)*(u(m,M1_F2,k,j,i) - c1)
                           + (u(m,M1_F3,k,j,i) - c2)*(u(m,M1_F3,k,j,i) - c2)));
        a3 += 1.0;
        // the plain signed-permutation of the donor (no basis transform): the error a
        // copy of the source components would make, for scale
        a4 = fmax(a4, fabs(u(m,M1_F1,k,j,i) - 0.5*s));
      }, Kokkos::Max<Real>(es), Kokkos::Max<Real>(ev), Kokkos::Sum<Real>(ecnt),
         Kokkos::Max<Real>(ep));
      res[0*ncl + c] = es; res[1*ncl + c] = ev; res[2*ncl + c] = ecnt;
      res[3*ncl + c] = ep;
    }
    Kokkos::deep_copy(u0, ubak);
  }
  // ---------------- (2) the work-array face-normal pairs through pbval_th
  const bool wk = (pbval_th != nullptr) && (iw.extent_int(1) > M1_IW_V3);
  if (wk) {
    DvceArray5D<Real> wbak("m1_csvec_w", iw.extent(0), iw.extent(1), iw.extent(2),
                           iw.extent(3), iw.extent(4));
    Kokkos::deep_copy(wbak, iw);
    auto w = iw;
    par_for("m1_csvec_fw", DevExeSpace(), 0, nmb-1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const bool act = (i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
      Real xi, eta, s, c1, c2, f1, f2;
      chart(m, k, j, xi, eta);
      field(mpd(m), xi, eta, s, c1, c2, f1, f2);
      for (int n = 0; n < M1_NHALO_T; ++n) {
        w(m,M1HaloCompT(n),k,j,i) = act ? s : 1.0e20;
      }
      w(m,M1_IW_N1,k,j,i) = act ? 1.0 : 1.0e20;   // a unit n = rhat (the renorm no-op)
      w(m,M1_IW_N2,k,j,i) = act ? 0.0 : 1.0e20;
      w(m,M1_IW_N3,k,j,i) = act ? 0.0 : 1.0e20;
      w(m,M1_IW_V2,k,j,i) = act ? f1 : 1.0e20;
      w(m,M1_IW_V3,k,j,i) = act ? f2 : 1.0e20;
    });
    ImplicitHaloExchange(M1_NHALO_T, -1);
    for (int c = 0; c < ncl; ++c) {
      Real an = 0.0, at = 0.0, bn = 0.0, bt = 0.0;
      const int cc = c;
      Kokkos::parallel_reduce("m1_csvec_r2", Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          DevExeSpace(), {0, 0, 0, is}, {nmb, n3, n2, ie+1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &a1,
                    Real &a2, Real &a3, Real &a4) {
        int sd;
        if (gclass(m, k, j, sd) != cc) return;
        const int p = mpd(m);
        Real xi, eta, s, c1, c2, f1, f2;
        chart(m, k, j, xi, eta);
        field(p, xi, eta, s, c1, c2, f1, f2);
        const Real g2 = w(m,M1_IW_V2,k,j,i), g3 = w(m,M1_IW_V3,k,j,i);
        // (A) at the ghost position; the seam-normal component is that of the seam axis
        const Real gn = (sd == 3) ? g3 : g2, gt = (sd == 3) ? g2 : g3;
        const Real fn = (sd == 3) ? f2 : f1, ft = (sd == 3) ? f1 : f2;
        a1 = fmax(a1, fabs(gn - fn));
        a2 = fmax(a2, fabs(gt - ft));
        // (B) at the donor = the mirror of this block's own cell at the same depth
        if (sd == 2 || sd == 3) {
          const int jg = j - js + g0d(m,0), kg = k - ks + g0d(m,1);
          int jm = j, km = k;
          if (sd == 2) {jm = (jg < 0) ? (j + 2*(-jg) - 1) : (j - 2*(jg - npan + 1) + 1);}
          if (sd == 3) {km = (kg < 0) ? (k + 2*(-kg) - 1) : (k - 2*(kg - npan + 1) + 1);}
          Real xm, em;
          chart(m, km, jm, xm, em);
          Real q[3], fa[3], fb[3], fq[3];
          cubed_sphere::PanelToCart(p, xm, em, q);
          cubed_sphere::PanelFrame(p, fa, fb, fq);
          // the seam plane's normal u = (n -+ a) or (n -+ b) (unnormalised, |u|^2 = 2)
          Real uu[3];
          const bool hi = (sd == 2) ? (jg >= npan) : (kg >= npan);
          for (int d = 0; d < 3; ++d) {
            const Real ax = (sd == 2) ? fa[d] : fb[d];
            uu[d] = fq[d] - (hi ? ax : -ax);
          }
          const Real qu = q[0]*uu[0] + q[1]*uu[1] + q[2]*uu[2];
          for (int d = 0; d < 3; ++d) {q[d] -= qu*uu[d];}
          Real xd, ed;
          cubed_sphere::CartToPanel(p, q, xd, ed);
          Real sb, b1, b2, h1, h2;
          field(p, xd, ed, sb, b1, b2, h1, h2);
          const Real hn = (sd == 3) ? h2 : h1, ht = (sd == 3) ? h1 : h2;
          a3 = fmax(a3, fabs(gn - hn));
          a4 = fmax(a4, fabs(gt - ht));
        }
      }, Kokkos::Max<Real>(an), Kokkos::Max<Real>(at), Kokkos::Max<Real>(bn),
         Kokkos::Max<Real>(bt));
      res[4*ncl + c] = an; res[5*ncl + c] = at; res[6*ncl + c] = bn; res[7*ncl + c] = bt;
    }
    Kokkos::deep_copy(iw, wbak);
  }
#if MPI_PARALLEL_ENABLED
  {
    std::vector<Real> g(res.size());
    MPI_Allreduce(res.data(), g.data(), static_cast<int>(res.size()), MPI_ATHENA_REAL,
                  MPI_MAX, MPI_COMM_WORLD);
    std::vector<Real> cn(ncl), cg(ncl);
    for (int c = 0; c < ncl; ++c) {cn[c] = res[2*ncl + c];}
    MPI_Allreduce(cn.data(), cg.data(), ncl, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    res = g;
    for (int c = 0; c < ncl; ++c) {res[2*ncl + c] = cg[c];}
  }
#endif
  if (global_variable::my_rank == 0) {
    std::printf("CSVEC N=%d ng=%d axis (%.3f %.3f %.3f): |v| <= 1, s in [0.8,1.2]\n",
                npan, ng, o0, o1, o2);
    for (int c = 0; c < ncl; ++c) {
      char nm[32];
      if (c == 0) {std::snprintf(nm, sizeof(nm), "same-panel");
      } else if (c <= ng) {std::snprintf(nm, sizeof(nm), "seam depth %d", c);
      } else {std::snprintf(nm, sizeof(nm), "cube vertex");}
      std::printf("CSVEC %-13s ghosts %8.0f | u0 (covariant, resample): |ds| %.3e |dv| "
                  "%.3e |dF1| %.3e", nm, res[2*ncl + c], res[0*ncl + c],
                  res[1*ncl + c], res[3*ncl + c]);
      if (wk) {
        std::printf(" | iw (face-normal pairs): ghost-pos n %.3e t %.3e; donor-pos n %.3e"
                    " t %.3e", res[4*ncl + c], res[5*ncl + c], res[6*ncl + c],
                    res[7*ncl + c]);
      }
      std::printf("\n");
    }
  }
}

} // namespace radm1
