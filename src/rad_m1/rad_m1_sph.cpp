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
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
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
    if (impl_recon != M1_IRECON_DC) {why += " implicit_recon != dc;";}
    if (impl_tlim != M1_TLIM_NONE) {why += " implicit_trans_limit != none;";}
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
    if (impl_prec != 2) {wcs += " implicit_precond != rbgs_fwd;";}
    if (impl_halo_mpi) {wcs += " implicit_halo_mpi = true;";}
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
        << "= bicgstab, implicit_precond = rbgs_fwd, no "
        << "implicit_halo_mpi (and the spherical-polar list below); this input has:"
        << wcs << why << std::endl
        << "See /viper/ptmp2/jinma/rt_design_1003/CS_IMPLICIT_VET.md sect. 4."
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
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
    std::cout << "<rad_m1>: spherical-polar wedge ("
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
  // DEBUG <rad_m1>/cs_seam_avg (read only when named, default true): the C5 seam average
  // of the stored face state
  cs_seam_avg_on = true;
  if (pin_cs_ != nullptr && pin_cs_->DoesParameterExist("rad_m1","cs_seam_avg")) {
    cs_seam_avg_on = pin_cs_->GetBoolean("rad_m1","cs_seam_avg");
  }
  if (pbval_th != nullptr) {
    pbval_th->cs_noresample = true;
    // STAGE CS2: (N2,N3), (A2,A3), (V2,V3) of the work array are face-normal components
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

} // namespace radm1
