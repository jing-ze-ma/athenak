//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro.cpp
//! \brief implementation of Hydro class constructor and assorted other functions

#include <iostream>
#include <string>
#include <algorithm>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"
#include "shearing_box/shearing_box.hpp"
#include "shearing_box/orbital_advection.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"

namespace hydro {
//----------------------------------------------------------------------------------------
// constructor, initializes data structures and parameters

Hydro::Hydro(MeshBlockPack *ppack, ParameterInput *pin) :
    pmy_pack(ppack),
    u0("cons",1,1,1,1,1),
    w0("prim",1,1,1,1,1),
    coarse_u0("ccons",1,1,1,1,1),
    coarse_w0("cprim",1,1,1,1,1),
    u1("cons1",1,1,1,1,1),
    uflx("uflx",1,1,1,1,1),
    phi0("phi_fc",1,1,1,1),
    phicc0("phi_cc",1,1,1,1),
    u0wb("conswb",1,1,1,1,1),
    w0wb("primwb",1,1,1,1,1),
    w0facewb("primfwb",1,1,1,1,1),
    utest("utest",1,1,1,1,1),
    fofc("fofc",1,1,1,1) {
  // Total number of MeshBlocks on this rank to be used in array dimensioning
  int nmb = std::max((ppack->nmb_thispack), (ppack->pmesh->nmb_maxperrank));

  // (1) construct EOS object (no default)
  std::string eqn_of_state = pin->GetString("hydro","eos");
  // ideal gas EOS
  if (eqn_of_state.compare("ideal") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic) {
      peos = new IdealSRHydro(ppack, pin);
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      peos = new IdealGRHydro(ppack, pin);
    } else {
      peos = new IdealHydro(ppack, pin);
    }
    nhydro = 5;
  // general EOS (non-relativistic only)
  } else if (eqn_of_state.compare("general") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "<hydro>/eos = general cannot be used with SR/GR" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    peos = new GeneralHydro(ppack, pin);
    nhydro = 5;
  // isothermal EOS
  } else if (eqn_of_state.compare("isothermal") == 0) {
    if (pmy_pack->pcoord->is_special_relativistic ||
        pmy_pack->pcoord->is_general_relativistic) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "<hydro>/eos = isothermal cannot be used with SR/GR" << std::endl;
      std::exit(EXIT_FAILURE);
    } else {
      peos = new IsothermalHydro(ppack, pin);
      nhydro = 4;
    }
  // EOS string not recognized
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<hydro>/eos = '" << eqn_of_state << "' not implemented" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // (2) Initialize scalars, diffusion, source terms
  nscalars = pin->GetOrAddInteger("hydro","nscalars",0);

  // Viscosity (if requested in input file)
  if (pin->DoesParameterExist("hydro","isotropic_viscosity")) {
    pvisc = new Viscosity("hydro", ppack, pin);
  } else {
    pvisc = nullptr;
  }

  // Thermal conduction (if requested in input file)
  if (pin->DoesParameterExist("hydro","isotropic_conduction")) {
    if (peos->eos_data.is_ideal) {
      pcond = new Conduction("hydro", ppack, pin);
    } else {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Thermal conduction in hydro requires ideal gas EOS" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    pcond = nullptr;
  }

  // Source terms (if needed)
  if (pin->DoesBlockExist("hydro_srcterms")) {
    psrc = new SourceTerms("hydro_srcterms", ppack, pin);
  }

  // (3) read time-evolution option [already error checked in driver constructor]
  // Then initialize memory and algorithms for reconstruction and Riemann solvers
  std::string evolution_t = pin->GetString("time","evolution");

  // allocate memory for conserved and primitive variables
  // With AMR, maximum size of Views are limited by total device memory through an input
  // parameter, which in turn limits max number of MBs that can be created.
  {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(u0, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
    Kokkos::realloc(w0, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
    // derived thermodynamic variables are only needed for an expensive (general) EOS
    if (peos->eos_data.IsGeneral()) {
      Kokkos::realloc(wder, nmb, NDERIVED, ncells3, ncells2, ncells1);
      Kokkos::realloc(wtemp, nmb, ncells3, ncells2, ncells1);
    }
  }

  // allocate memory for conserved variables on coarse mesh
  if (ppack->pmesh->multilevel) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int n_ccells1 = indcs.cnx1 + 2*(indcs.ng);
    int n_ccells2 = (indcs.cnx2 > 1)? (indcs.cnx2 + 2*(indcs.ng)) : 1;
    int n_ccells3 = (indcs.cnx3 > 1)? (indcs.cnx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(coarse_u0, nmb, (nhydro+nscalars), n_ccells3, n_ccells2, n_ccells1);
    Kokkos::realloc(coarse_w0, nmb, (nhydro+nscalars), n_ccells3, n_ccells2, n_ccells1);
  }

  // allocate boundary buffers for conserved (cell-centered) variables
  pbval_u = new MeshBoundaryValuesCC(ppack, pin, false);
  pbval_u->InitializeBuffers((nhydro+nscalars));

  // Orbital advection and shearing box BCs (if requested in input file)
  if (pin->DoesBlockExist("shearing_box")) {
    porb_u = new OrbitalAdvectionCC(ppack, pin, (nhydro+nscalars));
    psbox_u = new ShearingBoxCC(ppack, pin, (nhydro+nscalars));
  } else {
    porb_u = nullptr;
    psbox_u = nullptr;
  }
        
  // determine if etotgrav and local well-balanced scheme is enabled
  use_etotgrav = pin->GetOrAddBoolean("hydro","etotgrav",false);
  use_wellbalance_dynamic = pin->GetOrAddBoolean("hydro","wellbalance_dynamic",false);
  use_wb_x1 = pin->GetOrAddBoolean("hydro","wb_x1",false);
  use_wb_x2 = pin->GetOrAddBoolean("hydro","wb_x2",false);
  use_wb_x3 = pin->GetOrAddBoolean("hydro","wb_x3",false);
  use_wb_rho = pin->GetOrAddBoolean("hydro","wb_rho",false);
  // allocate array of flags used with etotgrav
  if (use_etotgrav || use_wellbalance_dynamic) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(phicc0, nmb, ncells3, ncells2, ncells1);
    Kokkos::realloc(phi0.x1f, nmb, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(phi0.x2f, nmb, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(phi0.x3f, nmb, ncells3+1, ncells2, ncells1);
  }
  if (use_wellbalance_dynamic) {
    // select well-balanced scheme assumption (no default).  Test for compatibility of options
    std::string wb_opt = pin->GetString("hydro","wb_option");
    if (wb_opt.compare("isothermal") == 0) {
      wb_option = WBOption::isothermal;
    } else if (wb_opt.compare("isodensity") == 0) {
      wb_option = WBOption::isodensity;
    } else if (wb_opt.compare("isentropic") == 0) {
      wb_option = WBOption::isentropic;
    } else if (wb_opt.compare("adaptive") == 0) {
      wb_option = WBOption::adaptive;
    // Error for anything else
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro> wb_option = '" << wb_opt
                << "' not implemented for hydro" << std::endl;
      std::exit(EXIT_FAILURE);
    }
      // select well-balanced scheme direction
      if (!use_wb_x1 && !use_wb_x2 && !use_wb_x3) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro> wb_direction not set!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }
        
  // determine if static wellbalance is enabled
  use_wellbalance_static = pin->GetOrAddBoolean("hydro","wellbalance_static",false);
  use_wellbalance_static_reconst_perturb = pin->GetOrAddBoolean("hydro","wellbalance_static_reconst",false);
  // allocate array of flags used with wellbalance
  if (use_wellbalance_static) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(u0wb, nmb, nhydro, ncells3, ncells2, ncells1);
    Kokkos::realloc(w0wb, nmb, nhydro, ncells3, ncells2, ncells1);
    Kokkos::realloc(w0facewb.x1f, nmb, nhydro, ncells3, ncells2, ncells1+1);
    Kokkos::realloc(w0facewb.x2f, nmb, nhydro, ncells3, ncells2+1, ncells1);
    Kokkos::realloc(w0facewb.x3f, nmb, nhydro, ncells3+1, ncells2, ncells1);
  }

  // for time-evolving problems, continue to construct methods, allocate arrays
  if (evolution_t.compare("stationary") != 0) {
    // determine if FOFC is enabled
    use_fofc = pin->GetOrAddBoolean("hydro","fofc",false);

    // select reconstruction method (default PLM)
    std::string xorder = pin->GetOrAddString("hydro","reconstruct","plm");
    if (xorder.compare("dc") == 0) {
      recon_method = ReconstructionMethod::dc;
    } else if (xorder.compare("plm") == 0) {
      recon_method = ReconstructionMethod::plm;
      // check that nghost > 2 with PLM+FOFC
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      if (use_fofc && indcs.ng < 3) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "FOFC and " << xorder << " reconstruction requires at "
          << "least 3 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
        std::exit(EXIT_FAILURE);
      }
    } else if (xorder.compare("ppm4") == 0 ||
               xorder.compare("ppmx") == 0 ||
               xorder.compare("wenoz") == 0) {
      // check that nghost > 2
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      if (indcs.ng < 3) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << xorder << " reconstruction requires at least 3 ghost zones, "
          << "but <mesh>/nghost=" << indcs.ng << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // check that nghost > 3 with PPM4(or PPMX or WENOZ)+FOFC
      if (use_fofc && indcs.ng < 4) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "FOFC and " << xorder << " reconstruction requires at "
          << "least 4 ghost zones, but <mesh>/nghost=" << indcs.ng << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (xorder.compare("ppm4") == 0) {
        recon_method = ReconstructionMethod::ppm4;
      } else if (xorder.compare("ppmx") == 0) {
        recon_method = ReconstructionMethod::ppmx;
      } else if (xorder.compare("wenoz") == 0) {
        recon_method = ReconstructionMethod::wenoz;
      }
    } else {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "<hydro> reconstruct = '" << xorder << "' not implemented"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // select Riemann solver (no default).  Test for compatibility of options
    std::string rsolver = pin->GetString("hydro","rsolver");
    // Special relativistic dynamic solvers
    if (pmy_pack->pcoord->is_special_relativistic) {
      if (evolution_t.compare("dynamic") == 0) {
        if (rsolver.compare("llf") == 0) {
          rsolver_method = Hydro_RSolver::llf_sr;
        } else if (rsolver.compare("hlle") == 0) {
          rsolver_method = Hydro_RSolver::hlle_sr;
        } else if (rsolver.compare("hllc") == 0) {
          rsolver_method = Hydro_RSolver::hllc_sr;
        // Error for anything else
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro> rsolver = '" << rsolver
                    << "' not implemented for SR dynamics" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "kinematic dynamics not implemented for SR" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // General relativistic dynamic solvers
    } else if (pmy_pack->pcoord->is_general_relativistic) {
      if (evolution_t.compare("dynamic") == 0) {
        if (rsolver.compare("llf") == 0) {
          rsolver_method = Hydro_RSolver::llf_gr;
        } else if (rsolver.compare("hlle") == 0) {
          rsolver_method = Hydro_RSolver::hlle_gr;
        // Error for anything else
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro> rsolver = '" << rsolver
                    << "' not implemented for GR dynamics" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "kinematic dynamics not implemented for GR" <<std::endl;
        std::exit(EXIT_FAILURE);
      }

    // Non-relativistic dynamic solvers
    } else if (evolution_t.compare("dynamic") == 0) {
      // LLF solver
      if (rsolver.compare("llf") == 0) {
        rsolver_method = Hydro_RSolver::llf;
      // HLLE solver
      } else if (rsolver.compare("hlle") == 0) {
        rsolver_method = Hydro_RSolver::hlle;
      // HLLC solver
      } else if (rsolver.compare("hllc") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::hllc;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = hllc cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // LHLLC solver
      } else if (rsolver.compare("lhllc") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::lhllc;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = lhllc cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // HLLCLM solver
      } else if (rsolver.compare("hllclm") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::hllclm;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = hllclm cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // AUSMPUP solver
      } else if (rsolver.compare("ausmpup") == 0) {
        if (peos->eos_data.is_ideal) {
          rsolver_method = Hydro_RSolver::ausmpup;
        } else {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = ausmpup cannot be used with "
                    << "isothermal EOS" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      // Roe solver
      } else if (rsolver.compare("roe") == 0) {
        // The Roe solver linearizes the flux Jacobian about a Roe-averaged state whose
        // construction assumes an ideal gas throughout (see RoeFluxAdb, which is written
        // entirely in terms of gm1). Substituting a local Gamma_1 would NOT give a valid
        // Roe average: the linearized matrix would no longer satisfy the jump condition
        // F_R - F_L = A (U_R - U_L), which is the defining property of the scheme, and
        // shock speeds would be silently wrong. A correct general-EOS Roe solver needs
        // Vinokur-Montagne averaging, which is not implemented. Refuse instead.
        if (peos->eos_data.IsGeneral()) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "<hydro>/rsolver = roe is not supported with "
                    << "<hydro>/eos = general: the Roe average is only defined for an "
                    << "ideal gas. Use llf, hlle, hllc, hllclm, lhllc or ausmpup."
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        rsolver_method = Hydro_RSolver::roe;
      // Error for anything else
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro> rsolver = '" << rsolver << "' not implemented"
                  << " for dynamic problems" << std::endl;
        std::exit(EXIT_FAILURE);
      }

    // Non-relativistic kinematic solvers
    } else {
      // Advect solver
      if (rsolver.compare("advect") == 0) {
        rsolver_method = Hydro_RSolver::advect;
      } else {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<hydro> rsolver = '" << rsolver << "' not implemented"
                  << " for kinematic problems" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    // Final memory allocations
    {
      // allocate second registers, fluxes
      auto &indcs = pmy_pack->pmesh->mb_indcs;
      int ncells1 = indcs.nx1 + 2*(indcs.ng);
      int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
      int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
      Kokkos::realloc(u1,       nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(uflx.x1f, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(uflx.x2f, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);
      Kokkos::realloc(uflx.x3f, nmb, (nhydro+nscalars), ncells3, ncells2, ncells1);

      // allocate array of flags used with FOFC
      if (use_fofc) {
        Kokkos::realloc(fofc,  nmb, ncells3, ncells2, ncells1);
        Kokkos::realloc(utest, nmb, nhydro, ncells3, ncells2, ncells1);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// destructor

Hydro::~Hydro() {
  delete peos;
  delete pbval_u;
  if (pvisc != nullptr) {delete pvisc;}
  if (pcond != nullptr) {delete pcond;}
  if (psrc != nullptr) {delete psrc;}
}

} // namespace hydro
