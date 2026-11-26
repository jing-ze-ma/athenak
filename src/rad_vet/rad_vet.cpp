//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet.cpp
//! \brief RadiationVet constructor: reads the <rad_vet> keys, checks them and allocates
//! the module's arrays.
//!
//!  called from: rad_vet_rhs.cpp, rad_vet_steady.cpp, rad_vet_update.cpp,
//!    driver/driver.cpp, pgen/box_convection.cpp, pgen/ry_per_accretor.cpp,
//!    pgen/star_vet.cpp, pgen/tests/rad_vet_rhse.cpp, pgen/tests/rad_vet_tests2.cpp
//!  calls into: rad_vet_geometry.cpp, rad_vet_ledger.cpp, rad_vet_params.cpp,
//!    rad_vet_params_closure.cpp, rad_vet_params_implicit.cpp, rad_vet_params_matter.cpp,
//!    rays/sc_cart_setup.cpp

#include <float.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "utils/recipe.hpp"

namespace radvet {

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::RadiationVet
//! \brief reads <rad_vet>, checks the configuration and allocates the arrays.
//! A scheme key a restart file lacks takes the fresh default, with a NOTE (recipe::Fix).

RadiationVet::RadiationVet(MeshBlockPack *ppack, ParameterInput *pin) :
    u0("vet_u0",1,1,1,1,1),
    u1("vet_u1",1,1,1,1,1),
    coarse_u0("vet_coarse_u0",1,1,1,1,1),
    uflx("vet_uflx",1,1,1,1,1),
    pmy_pack(ppack) {
  // (1) parameters
  c_light = pin->GetReal("rad_vet","c_light");
  if (!(c_light > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/c_light must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat = c_light;   // (no reduced speed of light: chat_over_c was removed)
  // transport: implicit only
  {std::string tr = pin->GetOrAddString("rad_vet","transport","explicit");
  if (tr.compare("implicit") != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet>/transport = '" << tr << "': only transport = implicit "
      << "(the clean scheme)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }
  // (2) keys and checks, in this order (it is the order of the parameter dump)
  ReadSchemeKeys(pin);    // rad_vet_params.cpp
  CheckGeometry(pin);     // rad_vet_geometry.cpp
  ReadMatterKeys(pin);    // rad_vet_params_matter.cpp
  ReadClosureKeys(pin);   // rad_vet_params_closure.cpp

  // the VET face states read two cells beyond the face
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  if (indcs.ng < 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> requires at least 2 ghost zones, but <mesh>/nghost="
      << indcs.ng << std::endl;
    std::exit(EXIT_FAILURE);
  }

  dtnew = (FLT_MAX);
  dt_sub = 0.0;
  nsub = 1;

  // allocate arrays
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(u0, nmb, VET_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(u1, nmb, VET_NVAR, ncells3, ncells2, ncells1);
  // face fluxes of the explicit stage tasks (never in the implicit task list): dummies
  Kokkos::realloc(uflx.x1f, nmb, VET_NVAR, 1, 1, 1);
  Kokkos::realloc(uflx.x2f, nmb, VET_NVAR, 1, 1, 1);
  Kokkos::realloc(uflx.x3f, nmb, VET_NVAR, 1, 1, 1);
  // rho*kappa in the ghost cells too; vet_scatter: a 4th component, rho kappa_e
  Kokkos::realloc(opac, nmb, vscat ? (VET_NOPAC + 1) : VET_NOPAC, ncells3, ncells2,
                  ncells1);
  Kokkos::deep_copy(opac, 0.0);
  // force_reference: a valid View to capture; the pgen replaces it (SetForceReference)
  if (force_ref != VET_FREF_NONE) {
    Kokkos::realloc(arad_ref, nmb, ncells3, ncells2, ncells1);
    Kokkos::deep_copy(arad_ref, 0.0);
  }
  // MHD: |B|^2/2 of the current (b0) and U^n (b1) fields, every cell incl. ghosts
  if (fl_mhd) {
    Kokkos::realloc(emag0, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(emag1, nmb, ncells3, ncells2, ncells1);
  } else {
    Kokkos::realloc(emag0, 1, 1, 1, 1);
    Kokkos::realloc(emag1, 1, 1, 1, 1);
  }
  Kokkos::realloc(cnt, VET_NCNT);
  for (int n=0; n<VET_NCNT; ++n) {cnt.h_view(n) = 0.0;}
  cnt.modify_host();
  cnt.sync_device();

  if (ppack->pmesh->multilevel) {
    int nccells1 = indcs.cnx1 + 2*(indcs.ng);
    int nccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int nccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, VET_NVAR, nccells3, nccells2, nccells1);
  }

  // the implicit solver's parameters, checks and arrays
  ImplicitInit(pin);
  // clean sc rays with a reflecting top: the mirror reads the last build
  cln_scmem = vgd_on && sph_geom && ibc_x1max == VET_IBC_REFLECT;
  // the ray cut (RayCut): sc rays on 3-D sp, one MeshBlock along x1, no reflecting
  // top
  const Mesh *pm = ppack->pmesh;
  const bool ok = vgd_on && sph_geom &&
                  pm->three_d && !pm->multilevel && pm->mesh_indcs.nx1 == indcs.nx1 &&
                  !cln_scmem;
  cln_rcut_on = ok;
  if (!ok && global_variable::my_rank == 0) {
    std::cout << "<rad_vet> clean_ray_cut: not applicable here (staggered sc rays on "
              << "3-D spherical polar or the cubed sphere, one MeshBlock along x1 only)"
              << std::endl;
  }
  if (cln_rcut_on) {
    cln_rcut_wait = pin->GetOrAddInteger("rad_vet","ray_cut_wait",10);
    recipe::Fix(pin, "rad_vet", "ray_cut_r", "0.003", "0.01");
    recipe::OldStored(pin, "rad_vet", "ray_cut_r", "0.01", "0.003");
    cln_rcut_r = pin->GetOrAddReal("rad_vet","ray_cut_r",0.003);
    // the cut and its count in a restart file (RayCut)
    cln_rcut = pin->GetOrAddInteger("rad_vet","ray_cut_now",-1);
    cln_rcut_n = pin->GetOrAddInteger("rad_vet","ray_cut_count",0);
    pin_rcut_ = pin;
  }
  if (sph_geom) {SphericalCheck(pin);}
  if (cs_geom) {
    CubedInit();
  }

  // (5) boundary buffers
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  // <mesh>/cs_seam_blockend = extend also on E, F (bvals_cc.cpp's along-seam resample)
  pbval_u->cs_blockend_ext_ = CsSeamBlockendExtend(ppack, pin);
  pbval_u->InitializeBuffers(VET_NVAR);
  // (F2, F3) is the tangential pair; E and F1 (radial) are scalars
  pbval_u->SetVectorPairs(VET_NVAR, {{VET_F2, VET_F3}});

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet>: c=" << c_light
              << " closure=" << ((vet_sc || vgd_on) ? "sc" : "analytic")
              << std::endl;
    std::cout << "         kappa_p=" << kappa_p << " kappa_e=" << kappa_e
              << " kappa_f=" << kappa_f << " kappa_s=" << kappa_s
              << " coupling=" << (coupling ? "true" : "false")
              << " gas_feedback=" << (gas_feedback ? "true" : "false")
              << " arad=" << arad << std::endl;
    std::cout << "         opacity=" << (opacity_type == VET_OPAC_TABLE ? "table" : "")
              << " force_reference="
              << (force_ref == VET_FREF_WB_ARAD ? "wb_arad" :
                  (force_ref == VET_FREF_CLEAN ? "clean" :
                  (force_ref == VET_FREF_WB_MD ? "wb_arad_md" : "none")))
              << " temp_unit_kelvin=" << otab.tunit
              << " rho_unit_cgs=" << otab.dunit << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetOpacityTables
//! \brief hand the module the two opacity tables (opacity = table)

void RadiationVet::SetOpacityTables(const DvceArray2D<Real> &kr,
                                    const DvceArray2D<Real> &kp,
                                    const DvceArray1D<Real> &lT,
                                    const DvceArray1D<Real> &lD,
                                    const int nT, const int nD) {
  otab.kr = kr;
  otab.kp = kp;
  otab.lT = lT;
  otab.lD = lD;
  otab.nT = nT;
  otab.nD = nD;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetForceReference
//! \brief the reference acceleration of force_reference = wb_arad; resolves
//! force_reference_work = auto

void RadiationVet::SetForceReference(const DvceArray4D<Real> &a) {
  arad_ref = a;
  // force_reference_work = auto: split only if the pgen gives the reference work
  if (fref_wsplit_auto) {
    fref_wsplit = fref_wsplit_ok;
    fref_wsplit_auto = false;
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> force_reference_work = auto -> "
                << (fref_wsplit ? "split" : "full") << std::endl;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetX1maxBathColumns
//! \brief per-column incident bath E_bath(m,k,j) of the outer x1 Marshak face; also
//! allocates ibadv_x1max (lab-frame advective energy flux through that face)

void RadiationVet::SetX1maxBathColumns(const DvceArray3D<Real> &eb) {
  iebath_hi_col = eb;
  iebath_hi_col_on = true;
  Kokkos::realloc(ibadv_x1max, eb.extent_int(0), eb.extent_int(1), eb.extent_int(2));
  Kokkos::deep_copy(ibadv_x1max, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::SetX2FaceFlux
//! \brief the imposed x2-face fluxes (m,k,i) of implicit_bc_x2min|max = flux

void RadiationVet::SetX2FaceFlux(const DvceArray3D<Real> &flo,
                                 const DvceArray3D<Real> &fhi) {
  ifx2_lo = flo;
  ifx2_hi = fhi;
  ifx2_on = true;
}

bool RadiationVet::FrefWaccOn() const {
  // also the backward-Euler steps of a hesdirk2 run (first step, a
  // redo; not a helper solve)
  const bool t2be = (time_scheme == VET_TIME_HESDIRK2) &&
                    (!t2_ok || FrefWaccStage()) && !t2_inprobe;
  return fref_wsplit && ((time_scheme == VET_TIME_BE) || t2be);
}

bool RadiationVet::FrefWaccStage() const {
  return t2_ok;
}

Real RadiationVet::FrefWaccGam0(Real bdt, Real dt_mesh) const {
  if (FrefWaccStage()) return 0.0;   // each stage solve pays its hydro stage's work
  for (int s = 0; s < fref_hnst; ++s) {
    if (std::abs(bdt - fref_hbeta[s]*dt_mesh) <= 1.0e-12*fref_hbeta[s]*dt_mesh) {
      return fref_hgam0[s];
    }
  }
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
            << "<rad_vet>/force_reference_work = split with time_scheme = be: no hydro "
            << "stage has beta dt = " << bdt << " (dt = " << dt_mesh << "; integrators "
            << "rk1, rk2, rk3 supported)"
            << std::endl;
  std::exit(EXIT_FAILURE);
  return 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::FrefWaccOn, FrefWaccGam0, FrefWacc
//! \brief force_reference_work = split: the accumulator of the reference work the gas
//! receives in the hydro stages (fref_wacc)

DvceArray4D<Real> RadiationVet::FrefWacc() {
  const int n0 = arad_ref.extent_int(0), n1 = arad_ref.extent_int(1);
  const int n2 = arad_ref.extent_int(2), n3 = arad_ref.extent_int(3);
  if (fref_wacc.extent_int(0) != n0 || fref_wacc.extent_int(1) != n1 ||
      fref_wacc.extent_int(2) != n2 || fref_wacc.extent_int(3) != n3) {
    Kokkos::realloc(fref_wacc, n0, n1, n2, n3);
    Kokkos::deep_copy(fref_wacc, 0.0);
  }
  return fref_wacc;
}

//----------------------------------------------------------------------------------------
//! \fn RadiationVet::~RadiationVet
//! \brief waits for pending table work, frees IPC buffers, prints the reports

RadiationVet::~RadiationVet() {
  if (vgd_tab_state == 1) {vgd_tabf.wait();}   // the table prefetch
  if (vgd_tab_state == 3) {vgd_tex.fence();}
#if defined(KOKKOS_ENABLE_CUDA)
  if (vgd_ipc_state == 1) {     // vet_gd_halo_ipc
    (void)cudaDeviceSynchronize();
    for (void *q : vgd_ipc_open) {if (q != nullptr) {(void)cudaIpcCloseMemHandle(q);}}
    (void)cudaFree(vgd_ipc_rbuf);
    (void)cudaFree(vgd_ipc_sig);
    (void)cudaFree(vgd_ipc_err);
  }
#endif
#if defined(KOKKOS_ENABLE_HIP)
  if (vgd_ipc_state == 1) {     // vet_gd_halo_ipc (HIP)
    (void)hipDeviceSynchronize();
    for (void *q : vgd_ipc_open) {if (q != nullptr) {(void)hipIpcCloseMemHandle(q);}}
    (void)hipFree(vgd_ipc_rbuf);
    (void)hipFree(vgd_ipc_sig);
    (void)hipFree(vgd_ipc_err);
  }
#endif
  VetRayReport();
  ReportCounters();
  ImplicitReport();
  Time2Report();      // time_scheme = hesdirk2
  if (vet_sc) {ScCartReport();}
  ScCartFree();
  delete pbval_u;
  if (pbval_d3 != nullptr) {delete pbval_d3;}
  if (pbval_sg != nullptr) {delete pbval_sg;}
  if (pbval_vs != nullptr) {delete pbval_vs;}
  if (pbval_hr != nullptr) {delete pbval_hr;}
  if (pbval_gd != nullptr) {delete pbval_gd;}
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::SetSubsteps
//! \brief the stage dt of a mesh step: one solve per step (no sub-cycling)

int RadiationVet::SetSubsteps(Real dt_mesh) {
  nsub = 1;
  dt_sub = dt_mesh;
  return nsub;
}

}  // namespace radvet
