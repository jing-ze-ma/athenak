//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1.cpp
//! \brief constructor, parameter reading and allocation for the RadiationM1 class

#include <float.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_m1/rad_m1.hpp"

namespace radm1 {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

RadiationM1::RadiationM1(MeshBlockPack *ppack, ParameterInput *pin) :
    u0("m1_u0",1,1,1,1,1),
    u1("m1_u1",1,1,1,1,1),
    coarse_u0("m1_coarse_u0",1,1,1,1,1),
    uflx("m1_uflx",1,1,1,1,1),
    pmy_pack(ppack) {
  // (1) parameters
  c_light = pin->GetReal("rad_m1","c_light");
  if (!(c_light > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/c_light must be positive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat_over_c = pin->GetOrAddReal("rad_m1","chat_over_c",1.0);
  if (!(chat_over_c > 0.0) || chat_over_c > 1.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/chat_over_c must be in (0,1]" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  chat = chat_over_c*c_light;
  cfl_rad = pin->GetOrAddReal("rad_m1","cfl_rad",0.4);
  e_floor = pin->GetOrAddReal("rad_m1","e_floor",(FLT_MIN));
  subcycle = pin->GetOrAddBoolean("rad_m1","subcycle",true);

  // thick-limit flux.  Design sect. 3 makes ap_hll the eventual default; milestone 1a
  // implements the plain HLL flux only, so anything else is refused rather than
  // silently run as "none".
  thick_flux_str = pin->GetOrAddString("rad_m1","thick_flux","none");
  if (thick_flux_str.compare("none") != 0) {
    if (thick_flux_str.compare("ap_hll") == 0 ||
        thick_flux_str.compare("scaled") == 0) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/thick_flux = '" << thick_flux_str
        << "' is not implemented yet (milestone 1a implements 'none' only)" << std::endl;
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/thick_flux = '" << thick_flux_str
        << "' is not a valid choice (none | ap_hll | scaled)" << std::endl;
    }
    std::exit(EXIT_FAILURE);
  }

  // closure
  {std::string cl = pin->GetOrAddString("rad_m1","closure","m1");
  if (cl.compare("m1") == 0) {
    eddington = false;
  } else if (cl.compare("eddington") == 0) {
    eddington = true;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/closure = '" << cl << "' not implemented "
      << "(m1 | eddington)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }

  // reconstruction (PLM only for now: the design reconstructs (E, f_i) with PLM)
  {std::string xorder = pin->GetOrAddString("rad_m1","reconstruct","plm");
  if (xorder.compare("plm") == 0) {
    recon_method = ReconstructionMethod::plm;
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/reconstruct = '" << xorder << "' not implemented "
      << "(plm only in milestone 1a)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  }

  // PLM of (E, f_i) reads two cells beyond the face, so two ghost zones are the minimum
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  if (indcs.ng < 2) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> requires at least 2 ghost zones, but <mesh>/nghost="
      << indcs.ng << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2) who owns the mesh timestep.  With sub-cycling on, the radiation CFL is NOT
  // folded into Mesh::NewTimeStep -- unless nothing else sets a timestep at all, in
  // which case the mesh step would be unbounded.
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (!subcycle) || (!other_sets_dt);

  // (3) PD-ARS explicit weights.  With S == 0 the tableau (design sect. 5) is
  //   U1   = U^n + dt L(U^n)
  //   U^n1 = U^n + dt/2 [L(U^n) + L(U1)]
  // i.e. Heun / SSP-RK2 written in the AthenaK gam0/gam1/beta form.
  gam0[0] = 1.0;  gam1[0] = 0.0;  beta[0] = 1.0;
  gam0[1] = 0.5;  gam1[1] = 0.5;  beta[1] = 0.5;

  dtnew = (FLT_MAX);
  dt_sub = 0.0;
  nsub = 1;

  // (4) allocate arrays
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(u0, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(u1, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(uflx.x1f, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(uflx.x2f, nmb, M1_NVAR, ncells3, ncells2, ncells1);
  Kokkos::realloc(uflx.x3f, nmb, M1_NVAR, ncells3, ncells2, ncells1);

  if (ppack->pmesh->multilevel) {
    int nccells1 = indcs.cnx1 + 2*(indcs.ng);
    int nccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int nccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, M1_NVAR, nccells3, nccells2, nccells1);
  }

  // (5) boundary buffers
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  pbval_u->InitializeBuffers(M1_NVAR);

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: c=" << c_light << " chat/c=" << chat_over_c
              << " cfl_rad=" << cfl_rad << " thick_flux=" << thick_flux_str
              << " closure=" << (eddington ? "eddington" : "m1")
              << " subcycle=" << (subcycle ? "true" : "false") << std::endl;
  }
}

//----------------------------------------------------------------------------------------
// destructor

RadiationM1::~RadiationM1() {
  delete pbval_u;
}

//----------------------------------------------------------------------------------------
//! \fn int RadiationM1::SetSubsteps
//! \brief Set the substep dt and the number of substeps for a mesh step dt_mesh.
//! N_sub = ceil(dt_mesh/dt_rad) with dt_rad = cfl_rad*min(dx)/chat (= dtnew), and every
//! substep uses dt_mesh/N_sub.  When this module owns the mesh timestep there is nothing
//! to sub-cycle and N_sub = 1.

int RadiationM1::SetSubsteps(Real dt_mesh) {
  if (sets_mesh_dt || !subcycle) {
    nsub = 1;
  } else {
    nsub = 1;
    if (dtnew > 0.0 && dt_mesh > dtnew) {
      nsub = static_cast<int>(std::ceil(dt_mesh/dtnew));
      if (nsub < 1) {nsub = 1;}
    }
  }
  dt_sub = dt_mesh/static_cast<Real>(nsub);
  return nsub;
}

} // namespace radm1
