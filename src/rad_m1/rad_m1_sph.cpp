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
//!   * <rad_m1>/transport = implicit, closure = eddington, time_scheme = be;
//!   * one MeshBlock along x1 per column (implicit_partition = none), no SMR/AMR;
//!   * implicit_halo_mpi = false (the direct same-rank halo is allowed: without a pole
//!     it is the plain copy of the Cartesian mesh);
//!   * implicit_flux = central, implicit_recon = dc, implicit_trans_limit = none,
//!     implicit_vimp off, no dbg_tensor; implicit_offdiag is set to none (the Eddington
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

#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

namespace radm1 {

void RadiationM1::SphericalS1Check(ParameterInput *pin) {
  std::string why;
  auto *pm = pmy_pack->pmesh;
  if (transport != M1_TRANSPORT_IMPLICIT) {why += " transport != implicit;";}
  // STAGE S2: the uniaxial chi(f) closures (m1, minerbo, kershaw) as well; vet_sc and
  // the tau closure (plane-parallel column depth on index space) stay refused
  if (vet_sc || tau_closure) {why += " closure = vet_sc | tau;";}
  if (pm->multilevel) {why += " SMR/AMR;";}
  if (pm->mesh_indcs.nx1 != pm->mb_indcs.nx1) {
    why += " more than one MeshBlock along x1 (meshblock/nx1 must equal mesh/nx1);";
  }
  if (transport == M1_TRANSPORT_IMPLICIT) {
    if (impl_halo_mpi) {why += " implicit_halo_mpi = true;";}
    if (impl_flux != M1_IFLUX_CENTRAL) {why += " implicit_flux != central;";}
    if (impl_recon != M1_IRECON_DC) {why += " implicit_recon != dc;";}
    if (impl_tlim != M1_TLIM_NONE) {why += " implicit_trans_limit != none;";}
    if (impl_vimp) {why += " implicit_vimp = true;";}
    if (dbg_tensor != 0) {why += " dbg_tensor != none;";}
    if (time_scheme != M1_TIME_BE) {why += " time_scheme != be;";}
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
  if (!why.empty()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> on a spherical-polar mesh (stages S1, S2) supports "
      << "only transport = implicit, closure = eddington | m1 | minerbo | kershaw, "
      << "implicit_offdiag = auto | none (| lagged for m1/minerbo/kershaw), "
      << "time_scheme = be, one MeshBlock "
      << "along x1, no SMR, implicit_halo_mpi = false, implicit_flux = central, "
      << "implicit_recon = dc, implicit_trans_limit = none, no implicit_vimp, no "
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
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: spherical-polar wedge ("
              << (sph_q ? "stage S2: chi(f) closure, radial integrating factor, "
                          "lagged curvature, offdiag " : "stage S1: ")
              << (sph_q ? ((impl_offdiag == M1_OD_NONE) ? "none; " : "lagged; ") : "")
              << "implicit with areas/volumes/face distances from Coordinates"
              << (pm->use_grid_stretch_r || pm->use_grid_stretch_r_poly ?
                  " (stretched radial grid)" : "") << std::endl;
  }
}

} // namespace radm1
