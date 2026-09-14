//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.cpp
//! \brief Implements functions for Conduction class. This includes isotropic thermal
//! conduction, in which heat flux is proportional to negative local temperature gradient.
//! Conduction may be added to Hydro and/or MHD independently.

#include <float.h>
#include <cstdio>
#include <algorithm>
#include <limits>
#include <string>
#include <iostream> // cout

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "conduction.hpp"
#include "bvals/bvals.hpp"
#include "utils/rosseland.hpp"
#include "utils/eint_from_cons.hpp"

// <problem>/nan_report: the first x1 radiative-conduction face whose flux comes out
// non-finite, and the inputs that made it.  Pointers, not Views: a file-scope View
// outlives Kokkos::finalize.  Nothing is allocated or run unless the switch is on.
namespace {
DvceArray1D<int> *cnd_nanrep_cnt = nullptr;
DvceArray1D<Real> *cnd_nanrep_rec = nullptr;
int cnd_nanrep_lines = 0;
const int cnd_nanrep_maxlines = 400;
}  // namespace
#include "coordinates/coordinates.hpp"
#include "units/units.hpp"

// VanLeer Limiter which takes 2 slopes
KOKKOS_INLINE_FUNCTION
Real VLL2State(const Real a, const Real b) {
  if (a*b > 0) {
    return 2.0*a*b/(a+b);
  } else {
    return 0.0;
  }
}

// VanLeer Limiter which takes 4 slopes
KOKKOS_INLINE_FUNCTION
Real VLL4State(const Real a, const Real b, const Real c, const Real d) {
  return VLL2State(VLL2State(a,b), VLL2State(c,d));
}

//----------------------------------------------------------------------------------------
//! \fn Real TempDepKappa()
//! \brief Temperature-dependent conductivity given by Parker (1953) and Spitzer (1962)

KOKKOS_INLINE_FUNCTION
Real TempDepKappa(Real temp, Real limit) {
  if (temp < 6.5e4) {
    return 2.5e3 * pow(temp, 0.5);
  } else {
    return fmin(6.0e-7*pow(temp, 2.5), limit);
  }
}

//----------------------------------------------------------------------------------------
//! \brief Conduction constructor
// Note first argument passes string ("hydro" or "mhd") denoting in wihch class this
// object is being constructed, and therefore which <block> in the input file from which
// the parameters are read.
// Note that the coefficient of thermal conduction, kappa, corresponds to conductivity,
// not diffusivity. This is different from the coefficient used in Athena++.

Conduction::Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin) :
    pmy_pack(pp), my_block(block) {
  // Read parameters for isotropic thermal conduction (if any)
  if (pin->DoesParameterExist(block,"isotropic_conduction")) {
    iso_cond_type = pin->GetString(block,"isotropic_conduction");
    // Check for valid type
    if ((iso_cond_type.compare("constant") != 0) &&
        (iso_cond_type.compare("spitzer") != 0) &&
        (iso_cond_type.compare("spitzer_limited") != 0) &&
        (iso_cond_type.compare("radiative") != 0)) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Invalid choice for isotropic thermal conduction type" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // constant conductivity
    if (iso_cond_type.compare("constant") == 0) {
      kappa_iso = pin->GetReal(block,"kappa_iso");
    }
    kappa_iso_limit = pin->GetOrAddReal(block,"kappa_iso_limit",
                      static_cast<Real>(std::numeric_limits<float>::max()));
    if (iso_cond_type.compare("radiative") == 0) {
      if (pp->punit == nullptr) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "radiative conduction needs a <units> block"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      rad_met = pin->GetOrAddReal(block,"rad_met",0.0);
      // pressure cut in bar; the flux through the wall in erg/cm^2/s
      rad_pcut = pin->GetOrAddReal(block,"rad_pcut_bar",0.0)*1.0e6
                 /pp->punit->pressure_cgs();
      // negative: the problem generator sets it (deep_hot_jupiter_rt: sigma T_int^4)
      rad_flux_inner = pin->GetOrAddReal(block,"rad_flux_inner",0.0);
      if (rad_flux_inner > 0.0) {
        rad_flux_inner /= (pp->punit->pressure_cgs()*pp->punit->velocity_cgs());
      }
      rad_kappa_fac = pin->GetOrAddReal(block,"rad_kappa_fac",1.0);
      nan_report = pin->GetOrAddBoolean("problem","nan_report",false);
      rad_flux_limit = pin->GetOrAddBoolean(block,"rad_flux_limit",true);
      // the free-streaming ceiling of that limiter: c a T^4 = 4 sigma T^4, or the
      // historical sigma T^4 with rad_flim_legacy.  See conduction.hpp.
      rad_flim_legacy = pin->GetOrAddBoolean(block,"rad_flim_legacy",false);
      rad_flim_fac = rad_flim_legacy ? 1.0 : 4.0;
      rad_tau_lo = pin->GetOrAddReal(block,"rad_tau_lo",0.0);
      rad_tau_hi = pin->GetOrAddReal(block,"rad_tau_hi",0.0);
      rad_tau_mode = (rad_tau_hi > 0.0);
      // ceiling on the temperature entering kappa_rad; 0 = off (bitwise inert)
      rad_tmax = pin->GetOrAddReal(block,"rad_tmax_kappa",0.0);
      if (rad_tmax > 0.0 && global_variable::my_rank == 0) {
        std::cout << "Conduction: radiative kappa temperature ceiling rad_tmax_kappa = "
                  << rad_tmax << " K" << std::endl;
      }
      rad_cs_exact = pin->GetOrAddBoolean(block,"rad_cs_exact",true);
      rad_kappa_rmax = pin->GetOrAddReal(block,"rad_kappa_rmax",0.0);
      rad_kappa_above = pin->GetOrAddReal(block,"rad_kappa_above",0.0);
      rad_gate_rho = pin->GetOrAddReal(block,"rad_gate_rho",0.0);
      rad_gate_dex = pin->GetOrAddReal(block,"rad_gate_dex",0.5);
      if (rad_gate_rho > 0.0 && rad_kappa_rmax > 0.0) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_gate_rho and rad_kappa_rmax are two forms of the "
                  << "same switch: set one or the other, not both" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rad_gate_rho > 0.0 && !(rad_gate_dex > 0.0)) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_gate_dex must be positive" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      rad_blend_radial = pin->GetOrAddBoolean(block,"rad_blend_radial",true);
      // rad_blend_use_2s: the ramp faces carry w*F_2s, not w*(-K dT/dz).  See
      // conduction.hpp for the whole argument.
      rad_blend_use_2s = pin->GetOrAddInteger(block,"rad_blend_use_2s",0);
      if (rad_blend_use_2s < 0 || rad_blend_use_2s > 2) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_blend_use_2s is 0 (off), 1 (prescribed flux) or "
                  << "2 (defect correction)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      rad_implicit_x1 = pin->GetOrAddBoolean(block,"rad_implicit_x1",false);
      rad_cap_ang = pin->GetOrAddReal(block,"rad_cap_ang",0.0);
      rad_implicit_ang = pin->GetOrAddBoolean(block,"rad_implicit_ang",false);
      // rad_sts_all: ONE RKL1 loop over x1, x2 and x3.  It IS the transverse treatment,
      // so it implies rad_implicit_ang -- everything rad_implicit_ang allocates, checks,
      // drops from the timestep and builds in BuildAngularCoeffs is needed here too, and
      // the extra x1 faces are the only difference.  Set before any of the checks below,
      // so that the Cartesian / uniform-grid / SMR / nghost guards of rad_implicit_ang
      // cover it without being repeated.
      rad_sts_all = pin->GetOrAddBoolean(block,"rad_sts_all",false);
      if (rad_sts_all) rad_implicit_ang = true;
      // the STIFFNESS SPLIT and the once-per-cycle application of the RKL1 operator,
      // and the round-off margin of its substage count.  All three are switches on the
      // super-time-stepped operator and do nothing at all when it is off; see
      // conduction.hpp for what each one does and why.
      rad_sts_split = pin->GetOrAddBoolean(block,"rad_sts_split",false);
      rad_sts_split_x = pin->GetOrAddReal(block,"rad_sts_split_x",0.5);
      rad_sts_once = pin->GetOrAddBoolean(block,"rad_sts_once",false);
      rad_sts_margin = pin->GetOrAddReal(block,"rad_sts_margin",0.10);
      rad_ang_maxit = pin->GetOrAddInteger(block,"rad_ang_maxit",200);
      rad_ang_verbose = pin->GetOrAddBoolean(block,"rad_ang_verbose",false);
      // DIAGNOSTIC ONLY: the T-linearisation audit of ImplicitRadialUpdate
      rad_x1_verbose = pin->GetOrAddBoolean(block,"rad_x1_verbose",false);
      rad_x1_every = pin->GetOrAddInteger(block,"rad_x1_every",1);
      {
        std::string ksrc = pin->GetOrAddString(block,"rad_kappa_src","freedman");
        if (ksrc.compare("table") == 0) {
          rad_kappa_tab = true;
        } else if (ksrc.compare("table_rho") == 0) {
          rad_kappa_tab = true;
          rad_kappa_rho = true;
        } else if (ksrc.compare("freedman") != 0) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_kappa_src must be freedman, table or table_rho"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (rad_tau_mode) {
        if (!(rad_tau_lo > 0.0 && rad_tau_lo < rad_tau_hi)) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "need 0 < rad_tau_lo < rad_tau_hi" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        auto &indcs = pp->pmesh->mb_indcs;
        const int nmb = pp->nmb_thispack;
        const int ncells1 = indcs.nx1 + 2*indcs.ng;
        const int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
        const int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
        Kokkos::realloc(rad_w, nmb, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(rad_tauf, nmb, ncells3, ncells2, ncells1+1);
        if (rad_blend_use_2s) {
          if (!rad_blend_radial) {
            std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                      << std::endl << "rad_blend_use_2s replaces the RADIAL blend flux "
                      << "and needs rad_blend_radial = true" << std::endl;
            std::exit(EXIT_FAILURE);
          }
          Kokkos::realloc(rad_f2s, nmb, ncells3, ncells2, ncells1+1);
          Kokkos::deep_copy(rad_f2s, 0.0);
        }
      } else if (rad_blend_use_2s) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_blend_use_2s is a property of the tau blend and "
                  << "needs rad_tau_hi > 0" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rad_implicit_x1 && rad_blend_use_2s > 0) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_blend_use_2s puts part of the ramp's radial flux "
                  << "back on an EXPLICIT, lagged footing, and rad_implicit_x1 is on "
                  << "because the explicit radial radiative dt on this column is orders "
                  << "below the step. Measured on the He-star 1-D arms: dt collapses at "
                  << "cycle 2-3 in both modes. See conduction.hpp, rad_blend_use_2s."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rad_implicit_x1) {
        // hydro and MHD both: ImplicitRadialUpdate subtracts the magnetic energy from
        // the conserved state when the block is <mhd> (see MagEnergyCC), so the frozen
        // internal energy it linearises about is the internal energy in either system.
        // the solve is column-local, so the whole radial extent must be in one MeshBlock
        if (pp->pmesh->mb_indcs.nx1 != pp->pmesh->mesh_indcs.nx1) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_implicit_x1 needs the whole radial extent in "
                    << "one MeshBlock: <meshblock>/nx1 = " << pp->pmesh->mb_indcs.nx1
                    << " but <mesh>/nx1 = " << pp->pmesh->mesh_indcs.nx1 << std::endl;
          std::exit(EXIT_FAILURE);
        }
        auto &indcs = pp->pmesh->mb_indcs;
        const int nmb = pp->nmb_thispack;
        const int ncells1 = indcs.nx1 + 2*indcs.ng;
        const int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
        const int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
        Kokkos::realloc(imp_wrk, nmb, nimpw, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(imp_flag, 1);
        Kokkos::realloc(imp_rec, 8);
        Kokkos::realloc(imp_x1dg, 16);
      }
      if (rad_cap_ang > 0.0 && rad_implicit_ang) {
        // one treatment of the transverse operator or the other: the cap throttles the
        // explicit flux, the implicit solve removes it from the fluxes altogether, and
        // running both would cap an operator that is no longer there
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_cap_ang and rad_implicit_ang are two treatments "
                  << "of the same operator: set one or the other, not both" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rad_implicit_ang) {
        // v1 is CARTESIAN: the cubed-sphere face-normal derivative carries a metric
        // cross term that is not part of the 5-point stencil the solver inverts, and
        // the spherical-polar pole rows need their own treatment.  Both are step 2.
        if (pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_implicit_ang is Cartesian-only in this version"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // ... and uniform-grid only: the increment is exchanged through its own
        // cell-centred boundary object, which would have to prolongate/restrict it at a
        // level boundary for the operator to stay conservative there
        if (pp->pmesh->multilevel) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_implicit_ang does not support SMR/AMR"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (!pp->pmesh->multi_d) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_implicit_ang needs a 2D or 3D mesh: there is no "
                    << "transverse operator in 1D" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (rad_sts_all) {
        // the RKL1 stencil owns every interior x1 face, so a prescribed-flux face would
        // have to be cut out of IT as well as out of the tridiagonal solve; not done
        if (rad_blend_use_2s) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_blend_use_2s is not implemented for the "
                    << "super-time-stepped radial operator: use rad_implicit_x1 or the "
                    << "explicit x1 path" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // ONE operator over all three directions: the tridiagonal radial solve is not
        // part of it and running both would apply the radial operator twice
        if (rad_implicit_x1) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_sts_all and rad_implicit_x1 are two treatments "
                    << "of the radial direction: set one or the other, not both"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // the x1 faces of a column are all in the stencil and the PHYSICAL x1 faces are
        // left to the explicit path, exactly as ImplicitRadialUpdate leaves them, so the
        // whole x1 extent has to be in one MeshBlock: an interior x1 face that fell on a
        // block boundary would be closed by both operators and carry nothing at all
        if (pp->pmesh->mb_indcs.nx1 != pp->pmesh->mesh_indcs.nx1) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_sts_all needs the whole x1 extent in one "
                    << "MeshBlock: <meshblock>/nx1 = " << pp->pmesh->mb_indcs.nx1
                    << " but <mesh>/nx1 = " << pp->pmesh->mesh_indcs.nx1 << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (rad_sts_split || rad_sts_once) {
        // both are switches ON the super-time-stepped operator: there is nothing to
        // split off and nothing to defer when that operator is not running
        if (!rad_implicit_ang) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_sts_split/rad_sts_once need rad_implicit_ang "
                    << "or rad_sts_all" << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (rad_sts_split && !(rad_sts_split_x > 0.0 && rad_sts_split_x <= 1.0)) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_sts_split_x is the explicit row-sum budget x_i "
                  << "and must be in (0,1]" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (!(rad_sts_margin >= 0.0)) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_sts_margin must be >= 0" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rad_cap_ang > 0.0) {
        // the angular cap only makes sense once the radial direction is unconditionally
        // stable: with x1 still explicit the run dies on dt1 long before dt2/dt3 matter,
        // and dropping dt2/dt3 while dt1 still applies would be a silent half-measure
        if (!rad_implicit_x1) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_cap_ang > 0 requires rad_implicit_x1 = true"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (rad_cap_ang > 0.0 || rad_implicit_ang) {
        auto &indcs = pp->pmesh->mb_indcs;
        if (indcs.ng < 2) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_cap_ang needs nghost >= 2: the stiffness of a "
                    << "ghost cell reads its own outer face" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        const int nmb = pp->nmb_thispack;
        const int ncells1 = indcs.nx1 + 2*indcs.ng;
        const int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
        const int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
        Kokkos::realloc(cap_x, nmb, ncells3, ncells2, ncells1);
        if (rad_sts_all) Kokkos::realloc(cap_c1, nmb, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(cap_c2, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(cap_c3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::realloc(cap_cnt, 2);
        Kokkos::realloc(cap_rec, 6);
        // the stiffness split: one explicit-fraction array per direction in the stencil,
        // and one flag per MeshBlock.  Nothing is allocated when the switch is off.
        if (rad_sts_split) {
          if (rad_sts_all) Kokkos::realloc(cap_f1, nmb, ncells3, ncells2, ncells1+1);
          Kokkos::realloc(cap_f2, nmb, ncells3, ncells2+1, ncells1);
          Kokkos::realloc(cap_f3, nmb, ncells3+1, ncells2, ncells1);
          Kokkos::realloc(sts_blk, nmb);
          sts_blk_used = true;
        }
        if (rad_implicit_ang) {
          Kokkos::realloc(tr_st, nmb, ntrs, ncells3, ncells2, ncells1);
          Kokkos::realloc(tr_ya, nmb, 1, ncells3, ncells2, ncells1);
          Kokkos::realloc(tr_yb, nmb, 1, ncells3, ncells2, ncells1);
          Kokkos::realloc(tr_yc, nmb, 1, ncells3, ncells2, ncells1);
          // the coarse register every MeshBoundaryValuesCC call takes.  Never used
          // (SMR/AMR is refused above) but it has to exist and be the right shape.
          const int cc1 = indcs.cnx1 + 2*(indcs.ng);
          const int cc2 = (indcs.cnx2 > 1) ? (indcs.cnx2 + 2*(indcs.ng)) : 1;
          const int cc3 = (indcs.cnx3 > 1) ? (indcs.cnx3 + 2*(indcs.ng)) : 1;
          Kokkos::realloc(tr_ycoar, nmb, 1, cc3, cc2, cc1);
          // its own MeshBoundaryValues object, hence its own MPI_Comm_dup'd
          // communicator, so the substage traffic cannot collide with u0 or b0
          pbval_tr = new MeshBoundaryValuesCC(pp, pin, false);
          pbval_tr->InitializeBuffers(1);
        }
      }
    } else if (pin->GetOrAddBoolean(block,"rad_implicit_x1",false) ||
               pin->GetOrAddBoolean(block,"rad_implicit_ang",false) ||
               pin->GetOrAddBoolean(block,"rad_sts_all",false)) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                << std::endl << "rad_implicit_x1/rad_implicit_ang/rad_sts_all need "
                << "isotropic_conduction = radiative" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  Kokkos::realloc(dt_diag, ndtdiag);
}

//----------------------------------------------------------------------------------------
//! \brief Conduction destructor

Conduction::~Conduction() {
  // the COST of the super-time-stepped operator over the whole run, in one line: the
  // substage count is the number of stencil sweeps and halo exchanges it took, and it is
  // what rad_sts_split, rad_sts_once and rad_sts_margin are there to reduce.  Only with
  // rad_ang_verbose, and only from rank 0.
  if (sts_ncall > 0 && rad_ang_verbose && global_variable::my_rank == 0) {
    std::cout << "### rad_sts totals: " << sts_ncall << " calls, "
              << sts_nsub_tot << " substages, mean "
              << (static_cast<double>(sts_nsub_tot)/static_cast<double>(sts_ncall))
              << " per call" << std::endl;
  }
  if (pbval_tr != nullptr) delete pbval_tr;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::EnableDiag()
//! \brief Allocate the per-cycle diagnostic array.  NOT done in the constructor: the
//! diagnostic is a problem-generator option, and the pgen runs after AddPhysics has
//! built this object.  With AMR the pack size can change and this array would be stale,
//! but the runs this exists for (cubed-sphere deep hot Jupiter) have no AMR.

void Conduction::EnableDiag(int nmb, int n3, int n2, int n1) {
  diag = true;
  Kokkos::realloc(cond_diag, nmb, 6, n3, n2, n1);
  Kokkos::deep_copy(cond_diag, 0.0);
}

//----------------------------------------------------------------------------------------
//! \fn void AddHeatFluxes()
//! \brief Wrapper function that adds heat fluxes for different types of thermal
//! conduction to face-centered fluxes of conserved variables

void Conduction::AddHeatFluxes(const DvceArray5D<Real> &w0, const EOS_Data &eos,
    DvceFaceFld5D<Real> &flx) {
  if (iso_cond_type.compare("constant") == 0) {
    AddIsotropicHeatFluxConstCond(w0, eos, flx);
  } else if ((iso_cond_type.compare("spitzer") == 0) ||
             (iso_cond_type.compare("spitzer_limited") == 0)) {
    AddIsotropicHeatFluxSpitzerCond(w0, eos, flx);
  } else if (iso_cond_type.compare("radiative") == 0) {
    if (rad_tau_mode) BuildRadWeights(w0, eos);
    AddIsotropicHeatFluxRadiative(w0, eos, flx);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::BuildRadWeights
//! \brief the Rosseland optical depth of every x1 face from the top of the domain,
//! tau_R = sum kappa_R rho dr over the cells above, and the blend weight w(tau_R).
//! One serial sweep per column, parallel over columns; the top face has tau = 0 (the
//! two-stream RT starts from the same point).

void Conduction::BuildRadWeights(const DvceArray5D<Real> &w0, const EOS_Data &eos) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, ng = indcs.ng;
  const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &dx1_ = pmy_pack->pcoord->dx1;
  const Real gm1 = eos.gamma-1.0;
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real met = rad_met, kfac = rad_kappa_fac, lo = rad_tau_lo, hi = rad_tau_hi;
  auto &wf = rad_w;
  auto &tf = rad_tauf;
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;   // table axis is log rho, not log p
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  // the radiatively inert region above rad_kappa_rmax: tau accumulates at rad_kappa_above
  // there, so it does not grow through the corona and the blend weight stays 0
  const Real krmax = rad_kappa_rmax, kabove = rad_kappa_above;
  // the DENSITY gate (rad_gate_rho): the same inert medium, selected by what the gas is
  // rather than by where it is.  kappa_eff = G kappa_table + (1 - G) rad_kappa_above.
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  auto &x1v_t = pmy_pack->pcoord->x1v;
  par_for("radtau", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real tau = 0.0;
    tf(m,k,j,ie+1) = 0.0;
    wf(m,k,j,ie+1) = RadBlendWeight(0.0, lo, hi);
    for (int i=ie; i>=is-ng; --i) {
      const Real t = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      const Real p = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      const Real rho = w0(m,IDN,k,j,i)*dens_unit;
      const Real dr = (curv ? dx1_(m,k,j,i) : size.d_view(m).dx1)*len_unit;
      Real kr = (krmax > 0.0 && x1v_t(m,i) > krmax)
          ? kabove
          : (ktab
             ? RosselandTable(krt, krlT, krlP, krnT, krnP, t*temp_unit,
                              krho ? rho : p*pres_unit)
             : RosselandFreedman2014(t*temp_unit, p*pres_unit, met));
      if (gaterho > 0.0) {
        const Real g = RadGate(rho, gaterho, gatedex);
        kr = g*kr + (1.0 - g)*kabove;
      }
      tau += kfac*kr*rho*dr;
      tf(m,k,j,i) = tau;
      wf(m,k,j,i) = RadBlendWeight(tau, lo, hi);
    }
  });
  rad_w_built = true;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiativeKappa
//! \brief the radiative conductivity 16 sigma T^3/(3 kappa_R rho) in cgs, T in K, p in
//! dyn/cm^2, rho in g/cm^3

KOKKOS_INLINE_FUNCTION
Real RadiativeKappa(const Real tk, const Real pcgs, const Real rhocgs, const Real met,
                    const Real kfac) {
  const Real sigma_sb = 5.670374419e-5;
  const Real kr = kfac*RosselandFreedman2014(tk, pcgs, met);
  return 16.0*sigma_sb*tk*tk*tk/(3.0*kr*rhocgs);
}

//! \brief the same from a supplied kappa_R [cm^2/g]
KOKKOS_INLINE_FUNCTION
Real RadiativeKappaKR(const Real tk, const Real rhocgs, const Real kfac, const Real kr) {
  const Real sigma_sb = 5.670374419e-5;
  return 16.0*sigma_sb*tk*tk*tk/(3.0*kfac*kr*rhocgs);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadFaceKappa
//! \brief the radiative conductivity of a face in cgs from its temperature, pressure and
//! density, taking kappa_R either from the tabulated lookup (on a log p or a log rho
//! axis) or from the Freedman fit.  This is the ONE place the face conductivity is
//! defined: the explicit face flux (AddIsotropicHeatFluxRadiative) and the implicit
//! radial solve (ImplicitRadialUpdate) both call it, so the two operators cannot drift
//! apart.  `tmax` is <hydro>/rad_tmax_kappa: the temperature entering the conductivity
//! (its T^3 and the kappa_R lookup) is capped there, so the ceiling reaches the explicit
//! faces, the angular cap and the implicit radial solve through this one function.  The
//! gradient and the free-streaming limiter keep the true temperature.

KOKKOS_INLINE_FUNCTION
Real RadFaceKappa(const Real tk, const Real pcgs, const Real rhocgs, const bool ktab,
                  const DvceArray2D<Real> &krt, const DvceArray1D<Real> &krlT,
                  const DvceArray1D<Real> &krlP, const int krnT, const int krnP,
                  const bool krho, const Real met, const Real kfac,
                  const Real tmax) {
  const Real tka = KappaTemp(tk, tmax);
  return ktab
      ? RadiativeKappaKR(tka, rhocgs, kfac,
                         RosselandTable(krt, krlT, krlP, krnT, krnP, tka,
                                        krho ? rhocgs : pcgs))
      : RadiativeKappa(tka, pcgs, rhocgs, met, kfac);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadFaceKCode
//! \brief the face CONDUCTIVITY in CODE units: the cgs conductivity `kap` that
//! RadFaceKappa returns, multiplied by the flux limiter evaluated on the frozen
//! face-normal gradient `gradn` (code units, T per length) and converted to code units.
//! It is the coefficient of (T_j - T_i)/dl in the face flux, i.e. what the explicit face
//! flux -K grad T contributes per unit gradient.
//!
//! THE ONE DEFINITION used by every operator that needs K itself rather than the flux:
//! BuildAngularCoeffs forms the x1 (rad_sts_all), x2 and x3 face coefficients of the
//! super-time-stepped operator with it, so the three directions of that operator cannot
//! disagree about what K is, and none of them can drift from the explicit face flux.
//! (Conduction::ImplicitRadialUpdate keeps its own inline copy of the same three lines:
//! its expression multiplies the tau-blend weight in at a different place in the product,
//! so folding it in here would move its last bits, and that solve is bitwise frozen.)

KOKKOS_INLINE_FUNCTION
Real RadFaceKCode(const Real kap, const Real tk, const Real gradn, const bool limit,
                  const Real temp_unit, const Real len_unit, const Real eflx_unit,
                  const Real ffac) {
  Real lf = 1.0;
  if (limit) {
    const Real f = -kap*gradn*temp_unit/len_unit;
    // ffac = Conduction::rad_flim_fac: 4 (c a T^4, the true free-streaming flux) or 1
    // (sigma T^4, the pre-2026-09-14 behaviour).  See rad_flim_legacy.
    const Real ffree = ffac*5.670374419e-5*tk*tk*tk*tk;
    lf = (ffree > 0.0) ? 1.0/sqrt(1.0 + SQR(f/ffree)) : 0.0;
  }
  return kap*lf*temp_unit/len_unit/eflx_unit;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::BuildAngularCoeffs
//! \brief the FROZEN coefficients of the TRANSVERSE (x2/x3) radiative operator: the face
//! coefficient C_f = A_f K_f/(dl_f sin alpha) on every x2 and x3 face, the per-cell
//! explicit stiffness x_i = beta_dt alpha_i (sum over its 4 transverse faces of C_f)/V_i
//! = dt/dt_cond,transverse, and -- when the transverse operator is solved implicitly
//! (rad_implicit_ang) -- the frozen temperature and 1/(rho c_v) of every cell.  K_f is
//! the SAME RadFaceKappa x flux limiter x tau-blend weight the x2/x3 flux kernels use
//! (divided by sin(alpha) on the cubed sphere), so the capped, the implicit and the
//! explicit operator cannot disagree about what K is.
//!
//! Called from AddIsotropicHeatFluxRadiative, i.e. inside the stage and BEFORE the RK
//! update, so w0 is the state every consumer linearises about and its ghost cells are
//! the ones the last exchange filled.  Everything is built one cell into the x2/x3
//! ghosts, because the first interior face needs the cell on its other side.

void Conduction::BuildAngularCoeffs(const DvceArray5D<Real> &w0, const EOS_Data &eos,
                                    const Real beta_dt) {
  const Real capang = rad_cap_ang;
  const bool impang = rad_implicit_ang;
  if (!(capang > 0.0) && !impang) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  const Real gm1 = eos.gamma-1.0;
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real eflx_unit = pres_unit*pmy_pack->punit->velocity_cgs();
  const Real met = rad_met, kfac = rad_kappa_fac;
  const Real tmax = rad_tmax;
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  auto &wf = rad_w;
  const bool limit = rad_flux_limit;
  const Real ffac = rad_flim_fac;      // 4 sigma T^4, or 1 with rad_flim_legacy
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  const Real krmax = rad_kappa_rmax;
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  // the face conductivity in code units, IDENTICAL in form and order to the face_kcode
  // of AddIsotropicHeatFluxRadiative (which is the coefficient of (T_j - T_i)/dl in the
  // face flux it adds)
  auto face_kcode = [=] (const Real tl, const Real tr, const Real pl, const Real pr,
                         const Real dl_, const Real dr_, const Real gradn) {
    const Real pf = 0.5*(pl + pr);
    if (pf < pcut) return 0.0;
    const Real tk = 0.5*(tl + tr)*temp_unit;
    const Real rhof = 0.5*(dl_ + dr_)*dens_unit;
    const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                  krnT, krnP, krho, met, kfac, tmax);
    return RadFaceKCode(kap, tk, gradn, limit, temp_unit, len_unit, eflx_unit, ffac);
  };
  const bool cs = pmy_pack->pmesh->use_cubed_sphere && rad_cs_exact;
  auto &sinc_ = pmy_pack->pcoord->sin_cell;
  auto &cosc_ = pmy_pack->pcoord->cos_cell;
  auto tcell = [=] (const int m, const int k, const int j, const int i) {
    return gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
  };
  const Real capbdt = beta_dt;
  // rad_sts_all: the x1 faces join the stencil, and the radial tau-blend weight comes
  // with them (rad_blend_radial), exactly as it does in the explicit x1 face flux
  const bool sts1 = rad_sts_all;
  const bool blend_r = rad_blend_radial;
  auto capx = cap_x;
  auto capc1 = cap_c1;
  auto capc2 = cap_c2;
  auto capc3 = cap_c3;
  auto capcnt = cap_cnt;
  auto caprec = cap_rec;
  auto &vol_ = pmy_pack->pcoord->volume;
  auto &area2_ = pmy_pack->pcoord->area.x2f;
  auto &area3_ = pmy_pack->pcoord->area.x3f;
  auto eos_ = eos;
  if (capang > 0.0) Kokkos::deep_copy(capcnt, 0);
  par_for("radcapc2", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+2, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    capc2(m,k,j,i) = 0.0;
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = tcell(m,k,j-1,i), tr = tcell(m,k,j,i);
    const Real pl = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx2_(m,k,j-1,i) + dx2_(m,k,j,i)) : size.d_view(m).dx2;
    Real wt = taumode ? 0.25*(wf(m,k,j-1,i) + wf(m,k,j-1,i+1)
                              + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k,j-1,i) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
    Real gradn = (tr - tl)/dl;
    Real sn = 1.0;
    if (cs && three_d) {
      const Real c = 0.5*(cosc_(m,k,j-1) + cosc_(m,k,j));
      sn = 0.5*(sinc_(m,k,j-1) + sinc_(m,k,j));
      const Real ge = 0.5*((tcell(m,k+1,j-1,i) - tcell(m,k-1,j-1,i))
                            /(0.5*dx3_(m,k-1,j-1,i) + dx3_(m,k,j-1,i)
                              + 0.5*dx3_(m,k+1,j-1,i))
                         + (tcell(m,k+1,j,i) - tcell(m,k-1,j,i))
                            /(0.5*dx3_(m,k-1,j,i) + dx3_(m,k,j,i)
                              + 0.5*dx3_(m,k+1,j,i)));
      gradn = (gradn - c*ge)/sn;
    }
    // the coefficient of (T_j - T_i)/dl in the face flux is K/sin(alpha): the metric
    // normalisation of grad(xi) is part of the stiffness, exactly as it is part of the
    // dt2 the cap replaces (SQR(d2)*s2/k2 in NewTimeStep)
    const Real kc = wt*face_kcode(tl, tr, pl, pr, w0(m,IDN,k,j-1,i),
                                  w0(m,IDN,k,j,i), gradn);
    const Real af = curv ? area2_(m,k,j,i) : 1.0/size.d_view(m).dx2;
    capc2(m,k,j,i) = kc*af/(dl*sn);
  });

  if (three_d) {
    par_for("radcapc3", DevExeSpace(), 0, nmb1, ks-1, ke+2, js-1, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      capc3(m,k,j,i) = 0.0;
      if (krmax > 0.0 && x1v_(m,i) > krmax) return;
      const Real tl = tcell(m,k-1,j,i), tr = tcell(m,k,j,i);
      const Real pl = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
      const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      const Real dl = curv ? 0.5*(dx3_(m,k-1,j,i) + dx3_(m,k,j,i))
                           : size.d_view(m).dx3;
      Real wt = taumode ? 0.25*(wf(m,k-1,j,i) + wf(m,k-1,j,i+1)
                                + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
      if (gaterho > 0.0) {
        wt *= RadGate(0.5*(w0(m,IDN,k-1,j,i) + w0(m,IDN,k,j,i))*dens_unit,
                      gaterho, gatedex);
      }
      Real gradn = (tr - tl)/dl;
      Real sn = 1.0;
      if (cs) {
        const Real c = 0.5*(cosc_(m,k-1,j) + cosc_(m,k,j));
        sn = 0.5*(sinc_(m,k-1,j) + sinc_(m,k,j));
        const Real gx = 0.5*((tcell(m,k-1,j+1,i) - tcell(m,k-1,j-1,i))
                              /(0.5*dx2_(m,k-1,j-1,i) + dx2_(m,k-1,j,i)
                                + 0.5*dx2_(m,k-1,j+1,i))
                           + (tcell(m,k,j+1,i) - tcell(m,k,j-1,i))
                              /(0.5*dx2_(m,k,j-1,i) + dx2_(m,k,j,i)
                                + 0.5*dx2_(m,k,j+1,i)));
        gradn = (gradn - c*gx)/sn;
      }
      const Real kc = wt*face_kcode(tl, tr, pl, pr, w0(m,IDN,k-1,j,i),
                                    w0(m,IDN,k,j,i), gradn);
      const Real af = curv ? area3_(m,k,j,i) : 1.0/size.d_view(m).dx3;
      capc3(m,k,j,i) = kc*af/(dl*sn);
    });
  }

  // rad_sts_all: the x1 face coefficient, from the SAME face_kcode (and therefore the
  // same RadFaceKappa, the same flux limiter and the same tau-blend weight) as the x2/x3
  // faces above and as the explicit x1 face flux.  Cartesian only -- rad_sts_all refuses
  // curvilinear meshes -- so A_f = 1/dx1 and V_i = 1, the form the flux-divergence of
  // Hydro::RKUpdate applies on a Cartesian grid.
  if (sts1) {
    par_for("radstsc1", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      capc1(m,k,j,i) = 0.0;
      // the two PHYSICAL x1 faces are NOT in the stencil.  They keep the explicit
      // treatment of AddIsotropicHeatFluxRadiative -- the imposed internal flux
      // rad_flux_inner at the bottom, the ghost-based gradient at the top -- which is
      // exactly what ImplicitRadialUpdate leaves them, and it is what makes the interior
      // fluxes telescope so that the operator moves energy without creating it.
      if (i == is || i == ie+1) return;
      if (krmax > 0.0 && x1v_(m,i) > krmax) return;
      const Real tl = tcell(m,k,j,i-1), tr = tcell(m,k,j,i);
      const Real pl = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
      const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      const Real dl = size.d_view(m).dx1;
      Real wt = (taumode && blend_r) ? wf(m,k,j,i) : 1.0;
      if (gaterho > 0.0) {
        wt *= RadGate(0.5*(w0(m,IDN,k,j,i-1) + w0(m,IDN,k,j,i))*dens_unit,
                      gaterho, gatedex);
      }
      const Real kc = wt*face_kcode(tl, tr, pl, pr, w0(m,IDN,k,j,i-1),
                                    w0(m,IDN,k,j,i), (tr - tl)/dl);
      const Real af = 1.0/size.d_view(m).dx1;
      capc1(m,k,j,i) = kc*af/dl;
    });
  }

  par_for("radcapx", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    Real rcv = d/gm1;
    if (gen) rcv = d*eos_.SpecificHeatCv(d, w0(m,IEN,k,j,i), wtemp_(m,k,j,i));
    Real sum = capc2(m,k,j,i) + capc2(m,k,j+1,i);
    if (three_d) sum += capc3(m,k,j,i) + capc3(m,k+1,j,i);
    const Real vi = curv ? vol_(m,k,j,i) : 1.0;
    capx(m,k,j,i) = (rcv > 0.0 && vi > 0.0 && sum > 0.0) ? capbdt*sum/(rcv*vi) : 0.0;
  });

  if (capang > 0.0) {
    // the diagnostic: how stiff the angular operator actually is, and where.  Active
    // cells only -- a ghost is somebody else's cell and would be double counted.
    {
      const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
      const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;
      Kokkos::ValLocScalar<Real, int> xloc;
      int nover = 0;
      Kokkos::parallel_reduce("radcapdiag",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
      KOKKOS_LAMBDA(const int &idx, Kokkos::ValLocScalar<Real, int> &xres, int &nov) {
        const int m = idx/nkji_;
        const int k = (idx - m*nkji_)/nji_ + ks;
        const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
        const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
        const Real xv = capx(m,k,j,i);
        if (xv > capang) ++nov;
        if (xv > xres.val) { xres.val = xv; xres.loc = idx; }
      }, Kokkos::MaxLoc<Real, int>(xloc), nover);
      const int xl = xloc.loc;
      const Real xv = xloc.val;
      par_for("radcaploc", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
        caprec(0) = xv;
        if (xl >= 0) {
          const int m = xl/nkji_;
          const int k = (xl - m*nkji_)/nji_ + ks;
          const int j = (xl - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
          const int i = (xl - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
          caprec(1) = static_cast<Real>(m);
          caprec(2) = static_cast<Real>(k);
          caprec(3) = static_cast<Real>(j);
          caprec(4) = static_cast<Real>(i);
          caprec(5) = x1v_(m,i);
        } else {
          caprec(1) = -1.0; caprec(2) = -1.0; caprec(3) = -1.0;
          caprec(4) = -1.0; caprec(5) = -1.0;
        }
      });
      cap_diag_x = xv;
      cap_diag_over = nover;
    }
  }
  // rad_implicit_ang: the implicit transverse solve linearises T about this same frozen
  // w0 state, so it needs T*_i and alpha_i = 1/(rho_i c_v,i) on exactly the range the
  // stiffness above covers.  Nothing here is allocated or run when the flag is off.
  if (impang) {
    auto trst = tr_st;
    const int it_ = TRST, ia_ = TRSA;
    par_for("radtrst", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = w0(m,IDN,k,j,i);
      Real rcv = d/gm1;
      if (gen) rcv = d*eos_.SpecificHeatCv(d, w0(m,IEN,k,j,i), wtemp_(m,k,j,i));
      trst(m,it_,k,j,i) = tcell(m,k,j,i);
      trst(m,ia_,k,j,i) = (rcv > 0.0 && isfinite(rcv)) ? 1.0/rcv : 0.0;
    });
  }

  // ------------------------------------------------------------------------------------
  // rad_sts_split: THE STIFFNESS SPLIT.  Everything above this point is untouched, so
  // the operator is bitwise what it was when the switch is off.
  //
  // THE BUDGET.  NewTimeStep limits the EXPLICIT operator by
  //     dt <= cfl fac dx_d^2 rho c_v/K_d   per direction d,   fac = 1/(2 ndim),
  // and on a Cartesian mesh C_f = K/dx^2 and V_i = 1, so a direction at that limit
  // contributes 2 C dt alpha = cfl/ndim to the row sum
  //     x_i = beta_dt alpha_i (sum_f C_f)/V_i
  // and the whole row is bounded by cfl.  The row-sum form is the direction-independent
  // statement of the same limit, and it is the one that can be applied per face: give the
  // explicit part the budget x_i <= xsplit, share it evenly over the nf = 2 ndim faces of
  // the row, and hold each face to the smaller of the two budgets it sees,
  //     C_max,f = xsplit min(V_i/alpha_i, V_j/alpha_j)/(nf beta_dt),
  //     C_exp,f = min(C_f, C_max,f),   C_sts,f = C_f - C_exp,f >= 0.
  // Both cells then satisfy sum_f C_exp,f <= nf (budget/nf) = budget by construction --
  // for ANY dt, which is why the dt limiter needs no change -- and the expression is
  // symmetric in the two cells, so they (and the two MeshBlocks at a block boundary) form
  // bitwise the same number and the explicit part is exactly conservative.
  //
  // WHAT IS STORED.  cap_f* holds the explicit FRACTION C_exp,f/C_f, which is what the
  // face-flux kernels multiply their flux by (the flux they form is C_f (T_j - T_i) in
  // disguise, see AddIsotropicHeatFluxRadiative), and cap_c* is overwritten with C_sts so
  // that the RKL1 loop needs no change at all.  A face the RKL1 loop treats as CLOSED
  // gets fraction 0: the split redistributes the operator, it does not open faces.
  if (rad_sts_split) {
    const Real xsplit = rad_sts_split_x;
    // THE STEP THE BUDGET IS FOR.  Normally the stage's own beta_dt: the explicit part
    // is added to THIS stage's fluxes and has to be stable over it.  With rad_sts_once
    // the RKL1 loop runs once over the FULL dt, and the division has to be the same in
    // every stage or the two parts would not add up to the whole operator -- the
    // explicit part would be sized for one step and the remainder subtracted for
    // another.  Budgeting the full dt in that case is both consistent and stricter than
    // each stage needs, so the explicit part stays stable a fortiori.
    const Real spbdt = rad_sts_once ? pmy_pack->pmesh->dt : beta_dt;
    const int ndim = (three_d ? 2 : 1) + (sts1 ? 1 : 0) + 1;   // x2 [+x3] [+x1]
    const Real nf = 2.0*static_cast<Real>(ndim);
    auto capf1 = cap_f1;
    auto capf2 = cap_f2;
    auto capf3 = cap_f3;
    auto trst = tr_st;
    auto blk = sts_blk;
    const int ia_ = TRSA;
    auto &mb_bcs = pmy_pack->pmb->mb_bcs;
    Kokkos::deep_copy(blk, 0);
    // V_i/alpha_i = V_i rho_i c_v,i, the cell's heat capacity; a cell the linearisation
    // dropped (alpha = 0) carries no flux either way
    auto cap_i = [=] (const int m, const int k, const int j, const int i) {
      const Real ai = trst(m,ia_,k,j,i);
      const Real vi = curv ? vol_(m,k,j,i) : 1.0;
      return (ai > 0.0 && vi > 0.0) ? vi/ai : 0.0;
    };
    // the explicit fraction of a face, and the C_sts it leaves behind
    auto split_f = [=] (const Real cf, const Real hi, const Real hj, const bool open) {
      if (!open || !(cf > 0.0) || !(spbdt > 0.0)) return 0.0;
      const Real hmin = fmin(hi, hj);
      if (!(hmin > 0.0)) return 0.0;
      const Real cmax = xsplit*hmin/(nf*spbdt);
      return (cf > cmax) ? cmax/cf : 1.0;
    };
    par_for("radstssp2", DevExeSpace(), 0, nmb1, ks-1, ke+1, js, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      bool op = true;
      if (j == js) {
        const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::inner_x2);
        op = (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
      } else if (j == je+1) {
        const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::outer_x2);
        op = (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
      }
      const Real cf = capc2(m,k,j,i);
      const Real fr = split_f(cf, cap_i(m,k,j-1,i), cap_i(m,k,j,i), op);
      capf2(m,k,j,i) = fr;
      capc2(m,k,j,i) = cf*(1.0 - fr);
      if (capc2(m,k,j,i) > 0.0) Kokkos::atomic_fetch_max(&blk(m), 1);
    });
    if (three_d) {
      par_for("radstssp3", DevExeSpace(), 0, nmb1, ks, ke+1, js-1, je+1, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        bool op = true;
        if (k == ks) {
          const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::inner_x3);
          op = (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
        } else if (k == ke+1) {
          const BoundaryFlag f = mb_bcs.d_view(m,BoundaryFace::outer_x3);
          op = (f == BoundaryFlag::block || f == BoundaryFlag::periodic);
        }
        const Real cf = capc3(m,k,j,i);
        const Real fr = split_f(cf, cap_i(m,k-1,j,i), cap_i(m,k,j,i), op);
        capf3(m,k,j,i) = fr;
        capc3(m,k,j,i) = cf*(1.0 - fr);
        if (capc3(m,k,j,i) > 0.0) Kokkos::atomic_fetch_max(&blk(m), 1);
      });
    }
    if (sts1) {
      par_for("radstssp1", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie+1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        // the two PHYSICAL x1 faces are not in the stencil and keep the full explicit
        // treatment AddIsotropicHeatFluxRadiative gives them; fraction 1 would double
        // them, fraction 0 is what the flux kernel is told to ignore
        if (i == is || i == ie+1) { capf1(m,k,j,i) = 0.0; return; }
        const Real cf = capc1(m,k,j,i);
        const Real fr = split_f(cf, cap_i(m,k,j,i-1), cap_i(m,k,j,i), true);
        capf1(m,k,j,i) = fr;
        capc1(m,k,j,i) = cf*(1.0 - fr);
        if (capc1(m,k,j,i) > 0.0) Kokkos::atomic_fetch_max(&blk(m), 1);
      });
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsotropicHeatFluxRadiative()
//! \brief the diffusion approximation to radiative transport in the optically thick
//! layers: heat flux -kappa_rad grad T on every face at or below the pressure cut, with
//! a smooth saturation at the free-streaming limit sigma T^4, and the imposed internal
//! flux
//! through the inner x1 wall.  On a curvilinear grid the radial gradient uses the
//! centroid spacing; the tangential gradients are negligible in an atmosphere and use
//! the cell widths.

void Conduction::AddIsotropicHeatFluxRadiative(const DvceArray5D<Real> &w0,
                                               const EOS_Data &eos,
                                               DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  Real gm1 = eos.gamma-1.0;
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real eflx_unit = pres_unit*pmy_pack->punit->velocity_cgs();   // erg/cm^2/s
  const Real met = rad_met, kfac = rad_kappa_fac, fin = rad_flux_inner;
  const Real tmax = rad_tmax;   // temperature ceiling in kappa_rad only (0 = off)
  // pressure cut, or the tau blend: with the blend every face is masked by its weight
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  auto &wf = rad_w;
  const bool limit = rad_flux_limit;
  const Real ffac = rad_flim_fac;      // 4 sigma T^4, or 1 with rad_flim_legacy
  const Real sigma_sb = 5.670374419e-5;

  // the heat flux across one face in CODE units, from the two adjacent cell states and
  // the centroid distance dl (code units); zero above the pressure cut
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;   // table axis is log rho, not log p
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  // no radiative diffusion above rad_kappa_rmax: the corona is transparent by fiat, and
  // an opacity of zero would divide by zero in RadiativeKappaKR, so the face is skipped
  // outright rather than handed kappa = 0
  const Real krmax = rad_kappa_rmax;
  // the DENSITY gate: kappa follows the gas, not the radius.  For a FACE the gate reads
  // the same face-averaged density the flux does, so the two operators switch off over
  // exactly the same faces.  Folded into the tau-blend weight, which is where every
  // conduction face flux is already scaled.
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  // gradn is the FACE-NORMAL temperature derivative in code units (T per length); the
  // caller forms it, which is where the grid enters
  auto face_flux = [=] (const Real tl, const Real tr, const Real pl, const Real pr,
                        const Real dl_, const Real dr_, const Real gradn) {
    const Real pf = 0.5*(pl + pr);
    if (pf < pcut) return 0.0;
    const Real tk = 0.5*(tl + tr)*temp_unit;
    const Real rhof = 0.5*(dl_ + dr_)*dens_unit;
    // the CAPPED temperature (rad_tmax_kappa) enters kappa_rad only; the gradient
    // below and the free-streaming limit use the true face temperature tk
    const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                  krnT, krnP, krho, met, kfac, tmax);
    Real f = -kap*gradn*temp_unit/len_unit;     // erg/cm^2/s, positive outward
    if (limit) {
      // saturate smoothly at the free-streaming flux sigma T^4: 0.3 % at F = 0.08 sigma
      // T^4,
      // where the diffusion approximation is still exact, and never above sigma T^4
      const Real ffree = ffac*sigma_sb*tk*tk*tk*tk;
      f /= sqrt(1.0 + SQR(f/ffree));
    }
    return f/eflx_unit;
  };

  // (The same coefficient AS A CONDUCTIVITY -- what the capped, implicit and
  // super-time-stepped operators need instead of the flux -- is RadFaceKCode above;
  // BuildAngularCoeffs forms every face coefficient with it.)

  // CUBED SPHERE: the xi and eta coordinate lines meet at an angle alpha (cos_cell,
  // sin_cell), so the face-normal derivative on a xi-face is
  //   dT/dn = (dT/dl_xi - cos(alpha) dT/dl_eta) / sin(alpha)
  // with dl the arc lengths -- the second term is the metric cross term, the 1/sin the
  // normalisation of grad(xi). The transverse derivative is centred across the two
  // cells the face separates. Radial faces are orthogonal to both and need nothing.
  const bool cs = pmy_pack->pmesh->use_cubed_sphere && rad_cs_exact;
  auto &sinc_ = pmy_pack->pcoord->sin_cell;
  auto &cosc_ = pmy_pack->pcoord->cos_cell;
  auto tcell = [=] (const int m, const int k, const int j, const int i) {
    return gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
  };

  // <problem>/nan_report device counter/record; untouched when the switch is off
  const bool nanrep_c = nan_report;
  if (nanrep_c && cnd_nanrep_cnt == nullptr) {
    cnd_nanrep_cnt = new DvceArray1D<int>("cnd_nanrep_cnt", 1);
    cnd_nanrep_rec = new DvceArray1D<Real>("cnd_nanrep_rec", 16);
  }
  auto cndcnt = nanrep_c ? *cnd_nanrep_cnt : DvceArray1D<int>("d", 1);
  auto cndrec = nanrep_c ? *cnd_nanrep_rec : DvceArray1D<Real>("d", 1);
  if (nanrep_c) {
    Kokkos::deep_copy(cndcnt, 0);
    Kokkos::deep_copy(cndrec, 0.0);
  }
  // per-cycle diagnostic (problem/diag_gid): locals only, never `this`
  const bool diag_ = diag;
  auto cdg = cond_diag;

  auto &flx1 = flx.x1f;
  // rad_implicit_x1: the INTERIOR x1 faces are handled by ImplicitRadialUpdate after the
  // RK update, so nothing is added to them here.  The two boundary faces stay explicit:
  // the inner one is the imposed wall flux (or the ghost-based gradient), the outer one
  // the ghost-based gradient, and neither is part of the tridiagonal system.
  // rad_sts_all does exactly the same, for the same reason: its RKL1 loop owns the
  // interior x1 faces and leaves the two physical ones here.
  const bool impx1 = rad_implicit_x1 || rad_sts_all;
  // rad_blend_use_2s: faces inside the tau ramp carry the TWO-STREAM's own flux, scaled
  // by w, instead of w*(-K dT/dz).  See conduction.hpp.  Applied here for both radial
  // paths, since a prescribed-flux face is outside the tridiagonal system either way.
  const bool use2s = (rad_blend_use_2s > 0) && taumode && blend_r && rad_f2s_ready;
  // mode 1 replaces the face flux outright; mode 2 adds only the defect and leaves the
  // face in the implicit system.  Without rad_implicit_x1 the two coincide.
  const bool pres2s = use2s && (rad_blend_use_2s == 1 || !(rad_implicit_x1||rad_sts_all));
  auto f2s_ = use2s ? rad_f2s : DvceArray4D<Real>("radf2sdummy", 1, 1, 1, 1);
  // rad_sts_split: the interior x1 faces of the rad_sts_all stencil DO carry a flux
  // here -- the explicit part C_exp of the split, as the fraction cap_f1 of the full
  // face flux (see BuildAngularCoeffs).  The two physical x1 faces are untouched: they
  // were never in the stencil and are already explicit at full strength.
  const bool splt1 = rad_sts_split && rad_sts_all;
  auto capf1_ = splt1 ? cap_f1 : DvceArray4D<Real>("radsp1dummy", 1, 1, 1, 1);
  // ...which means the frozen coefficients have to exist BEFORE this kernel, not after
  // it as they do on every other path (the call below is skipped when this one runs).
  // rad_sts_all implies rad_implicit_ang, which is refused on a 1D mesh, so the x2/x3
  // stencils of BuildAngularCoeffs are always in range here.
  if (splt1) BuildAngularCoeffs(w0, eos, stage_beta_dt);
  par_for("radcond1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // the imposed internal flux through the inner wall replaces the gradient there
    if (i == is && fin != 0.0 &&
        (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user ||
         mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect)) {
      flx1(m,IEN,k,j,i) += fin;
      if (diag_) cdg(m,0,k,j,i) = fin;
      return;
    }
    // rad_blend_use_2s: the ramp face is a PRESCRIBED-FLUX face.  Both cells sharing it
    // see this one number, so the exchange is conservative to round-off; and because the
    // two-stream's own share is (1 - w) of the SAME number, the two sum to F_2s exactly
    // and the handover term -d/dz[(1 - w)(F_2s - F_diff)] is identically zero.
    // MODE 2 marks the face for the defect correction below; it stays in the implicit
    // system, so the early return for impx1 must NOT be taken on it.
    const bool ramp2s = use2s && i > is && i < ie+1 &&
                        wf(m,k,j,i) > 0.0 && wf(m,k,j,i) < 1.0;
    if (pres2s && ramp2s) {
      const Real f2 = wf(m,k,j,i)*f2s_(m,k,j,i);
      flx1(m,IEN,k,j,i) += f2;
      if (diag_) cdg(m,0,k,j,i) = f2;
      return;
    }
    Real fsp = 1.0;
    if (impx1 && i > is && i < ie+1 && !ramp2s) {
      if (!splt1) return;
      fsp = capf1_(m,k,j,i);
      if (!(fsp > 0.0)) return;
    }
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? (x1v_(m,i) - x1v_(m,i-1)) : size.d_view(m).dx1;
    Real wt = (taumode && blend_r) ? wf(m,k,j,i) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k,j,i-1) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
    Real fcnd = fsp*wt*face_flux(tl, tr, pl, pr, w0(m,IDN,k,j,i-1),
                                 w0(m,IDN,k,j,i), (tr - tl)/dl);
    // MODE 2: replace the explicit contribution by the DEFECT w (F_2s - F_diff*).  The
    // implicit solve supplies w F_diff(T_new) across this same face with the same weight
    // wt, so the two sum to w F_2s at the frozen state.  Both cells see this one number.
    if (ramp2s) fcnd = wt*f2s_(m,k,j,i) - fcnd;
    flx1(m,IEN,k,j,i) += fcnd;
    if (diag_) cdg(m,0,k,j,i) = fcnd;
    // --- <problem>/nan_report: record the first face whose conduction flux, or the
    // total flux it lands in, is not finite, with the inputs that produced it.
    if (nanrep_c) {
      if (!isfinite(fcnd) || !isfinite(flx1(m,IEN,k,j,i))) {
        if (Kokkos::atomic_fetch_add(&cndcnt(0), 1) == 0) {
          const Real tk = 0.5*(tl + tr)*temp_unit;
          const Real rhof = 0.5*(w0(m,IDN,k,j,i-1) + w0(m,IDN,k,j,i))*dens_unit;
          const Real pf = 0.5*(pl + pr);
          const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                        krnT, krnP, krho, met, kfac, tmax);
          cndrec(0) = static_cast<Real>(m);
          cndrec(1) = static_cast<Real>(k);
          cndrec(2) = static_cast<Real>(j);
          cndrec(3) = static_cast<Real>(i);
          cndrec(4) = tl*temp_unit;
          cndrec(5) = tr*temp_unit;
          cndrec(6) = kap;
          cndrec(7) = wt;
          cndrec(8) = (tr - tl)/dl*temp_unit/len_unit;
          cndrec(9) = dl*len_unit;
          cndrec(10) = fcnd;
          cndrec(11) = flx1(m,IEN,k,j,i);
          cndrec(12) = w0(m,IDN,k,j,i-1);
          cndrec(13) = w0(m,IDN,k,j,i);
          cndrec(14) = pl;
          cndrec(15) = pr;
        }
      }
    }
  });
  if (nanrep_c && cnd_nanrep_lines < cnd_nanrep_maxlines) {
    auto hcc = Kokkos::create_mirror_view(cndcnt);
    Kokkos::deep_copy(hcc, cndcnt);
    if (hcc(0) > 0) {
      auto hcr = Kokkos::create_mirror_view(cndrec);
      Kokkos::deep_copy(hcr, cndrec);
      ++cnd_nanrep_lines;
      const int mb = static_cast<int>(hcr(0));
      std::cout << "### nan_report [conduction_x1_face] rank "
                << global_variable::my_rank << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time << ": " << hcc(0)
                << " bad face(s); first (m,k,j,i) = (" << mb << ","
                << static_cast<int>(hcr(1)) << "," << static_cast<int>(hcr(2)) << ","
                << static_cast<int>(hcr(3)) << ") gid = " << (pmy_pack->gids + mb)
                << " T_l = " << hcr(4) << " T_r = " << hcr(5) << " kappa = " << hcr(6)
                << " w_tau = " << hcr(7) << " dT/dx = " << hcr(8)
                << " dl = " << hcr(9) << " F_cond = " << hcr(10)
                << " flx1(IEN) = " << hcr(11) << " d_l = " << hcr(12)
                << " d_r = " << hcr(13) << " p_l = " << hcr(14)
                << " p_r = " << hcr(15) << std::endl;
    }
  }
  if (!multi_d) return;

  // ------------------------------------------------------------------------------------
  // rad_cap_ang: the CONSERVATIVE per-face cap on the explicit transverse operator.
  // BuildAngularCoeffs evaluates the frozen face coefficient A_f K_f/(dl_f sin alpha) on
  // the x2 and x3 faces and the per-cell stiffness x_i = dt/dt_cond,transverse; the flux
  // kernels below then multiply each face by min(1, cap/max(x_i,x_j)) -- one number per
  // face, so nothing is created or destroyed, only moved more slowly than an explicit
  // step could resolve.
  const Real capang = rad_cap_ang;
  auto capx = cap_x;
  auto capcnt = cap_cnt;
  auto caprec = cap_rec;
  if ((capang > 0.0 || rad_implicit_ang) && !splt1) {
    BuildAngularCoeffs(w0, eos, stage_beta_dt);
  }
  // rad_implicit_ang: the transverse fluxes are NOT added here at all.  The operator is
  // applied after the RK update, by Conduction::ImplicitTransverseUpdate (see
  // conduction_transverse.cpp), and the x2/x3 conduction timestep goes with it.
  // rad_sts_split is the exception: the part of each face the current step CAN carry
  // explicitly is added here, inside the stage, and only the remainder C_sts is left to
  // the super-time-stepped operator.  cap_f2/cap_f3 hold that fraction.
  const bool splt = rad_sts_split;
  if (rad_implicit_ang && !splt) return;
  auto capf2_ = splt ? cap_f2 : DvceArray4D<Real>("radsp2dummy", 1, 1, 1, 1);
  auto capf3_ = splt ? cap_f3 : DvceArray4D<Real>("radsp3dummy", 1, 1, 1, 1);

  auto &flx2 = flx.x2f;
  par_for("radcond2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx2_(m,k,j-1,i) + dx2_(m,k,j,i)) : size.d_view(m).dx2;
    Real wt = taumode ? 0.25*(wf(m,k,j-1,i) + wf(m,k,j-1,i+1)
                              + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k,j-1,i) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
    Real gradn = (tr - tl)/dl;
    if (cs && three_d) {
      const Real c = 0.5*(cosc_(m,k,j-1) + cosc_(m,k,j));
      const Real sn = 0.5*(sinc_(m,k,j-1) + sinc_(m,k,j));
      const Real ge = 0.5*((tcell(m,k+1,j-1,i) - tcell(m,k-1,j-1,i))
                            /(0.5*dx3_(m,k-1,j-1,i) + dx3_(m,k,j-1,i)
                              + 0.5*dx3_(m,k+1,j-1,i))
                         + (tcell(m,k+1,j,i) - tcell(m,k-1,j,i))
                            /(0.5*dx3_(m,k-1,j,i) + dx3_(m,k,j,i) + 0.5*dx3_(m,k+1,j,i)));
      gradn = (gradn - c*ge)/sn;
    }
    // rad_cap_ang: one scale factor per face, applied to both of its cells
    Real wtc = wt;
    if (capang > 0.0) {
      const Real xm = fmax(capx(m,k,j-1,i), capx(m,k,j,i));
      if (xm > capang) {
        wtc = wt*(capang/xm);
        Kokkos::atomic_fetch_add(&capcnt(1), 1);
      }
    }
    if (splt) wtc = wt*capf2_(m,k,j,i);
    const Real fadd = wtc*face_flux(tl, tr, pl, pr, w0(m,IDN,k,j-1,i),
                                    w0(m,IDN,k,j,i), gradn);
    flx2(m,IEN,k,j,i) += fadd;
    if (diag_) cdg(m,1,k,j,i) = fadd;
  });
  if (!three_d) return;

  auto &flx3 = flx.x3f;
  par_for("radcond3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx3_(m,k-1,j,i) + dx3_(m,k,j,i)) : size.d_view(m).dx3;
    Real wt = taumode ? 0.25*(wf(m,k-1,j,i) + wf(m,k-1,j,i+1)
                              + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k-1,j,i) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
    Real gradn = (tr - tl)/dl;
    if (cs) {
      const Real c = 0.5*(cosc_(m,k-1,j) + cosc_(m,k,j));
      const Real sn = 0.5*(sinc_(m,k-1,j) + sinc_(m,k,j));
      const Real gx = 0.5*((tcell(m,k-1,j+1,i) - tcell(m,k-1,j-1,i))
                            /(0.5*dx2_(m,k-1,j-1,i) + dx2_(m,k-1,j,i)
                              + 0.5*dx2_(m,k-1,j+1,i))
                         + (tcell(m,k,j+1,i) - tcell(m,k,j-1,i))
                            /(0.5*dx2_(m,k,j-1,i) + dx2_(m,k,j,i) + 0.5*dx2_(m,k,j+1,i)));
      gradn = (gradn - c*gx)/sn;
    }
    Real wtc = wt;
    if (capang > 0.0) {
      const Real xm = fmax(capx(m,k-1,j,i), capx(m,k,j,i));
      if (xm > capang) {
        wtc = wt*(capang/xm);
        Kokkos::atomic_fetch_add(&capcnt(1), 1);
      }
    }
    if (splt) wtc = wt*capf3_(m,k,j,i);
    const Real fadd = wtc*face_flux(tl, tr, pl, pr, w0(m,IDN,k-1,j,i),
                                    w0(m,IDN,k,j,i), gradn);
    flx3(m,IEN,k,j,i) += fadd;
    if (diag_) cdg(m,2,k,j,i) = fadd;
  });

  // rad_cap_ang report: how stiff the angular operator was and how much of it had to be
  // slowed down.  Same style and line budget as the rad_implicit_x1 report.
  if (capang > 0.0 && cap_lines < 400) {
    auto hcc = Kokkos::create_mirror_view(capcnt);
    Kokkos::deep_copy(hcc, capcnt);
    const bool active = (cap_diag_over > 0 || hcc(1) > 0);
    // the first 20 calls print from EVERY rank, so the startup lines carry the global
    // maximum of x_i and not just rank 0's share; after that only a rank that actually
    // had to cap something says anything
    if (active || cap_lines < 20) {
      ++cap_lines;
      auto hcr = Kokkos::create_mirror_view(caprec);
      Kokkos::deep_copy(hcr, caprec);
      const int mb = static_cast<int>(hcr(1));
      std::cout << "### rad_cap_ang rank " << global_variable::my_rank
                << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time
                << ": max x_i = " << hcr(0) << " (cap " << capang << ") at (m,k,j,i) = ("
                << mb << "," << static_cast<int>(hcr(2)) << ","
                << static_cast<int>(hcr(3)) << "," << static_cast<int>(hcr(4))
                << ") gid = " << (pmy_pack->gids + ((mb >= 0) ? mb : 0))
                << " r = " << hcr(5) << " | cells over cap = " << cap_diag_over
                << ", faces capped = " << hcc(1) << std::endl;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::ImplicitRadialUpdate
//! \brief BACKWARD-EULER solve of the radial flux-limited radiative diffusion, one
//! tridiagonal system per (m,k,j) column.  Runs after the RK update and the source
//! terms, and before the ghost exchange, so the values it writes are the ones that are
//! communicated.
//!
//! WHY.  The explicit radial operator is limited by dt ~ dr^2 rho c_v/kappa_rad with
//! kappa_rad = 16 sigma T^3/(3 kappa_R rho).  In an evacuated cell just above the
//! photosphere (rho 60x below its shell median, tau_cell ~ 1) that limit is ~1e-4 s
//! while the run takes 30 s steps, and since dt is set from the PREVIOUS state the
//! overshoot is unbounded: the cell heats past its neighbours, kappa_rad grows as T^3,
//! and T reaches 1e10 K in a single step.  Implicit in the ENERGY removes the
//! constraint entirely: the operator below is a symmetric M-matrix, so it is
//! unconditionally stable and cannot overshoot the neighbour equilibrium.
//!
//! WHAT IS SOLVED.  With e the internal energy density, T linearised about the frozen
//! post-update state,  T_i = T*_i + (e_i - e*_i)/(rho_i c_v,i),  and the frozen face
//! conductivity K_f (the SAME RadFaceKappa, tau-blend weight and flux limiter the
//! explicit face flux uses, the limiter evaluated on the frozen gradient),
//!     F_f = -K_f (T_i - T_{i-1})/dl_f
//!     e_i = e*_i + beta_dt (A_f F_f - A_{f+1} F_{f+1})/V_i
//! with exactly the areas and volumes Hydro::RKUpdate divides by, so the update is the
//! same flux-form divergence.  The two boundary faces (is and ie+1) are NOT in the
//! system -- they are added explicitly in AddIsotropicHeatFluxRadiative -- so the
//! interior fluxes telescope and sum_i V_i (e_i - e*_i) = 0 per column to round-off.
//! Only u0(IEN) is written; w0 is rebuilt by ConToPrim later in the same stage.
//!
//! THE UNKNOWN IS dT, NOT de.  The system above written on the energy increment is NOT
//! diagonally dominant: row i carries the cell's own alpha_i = 1/(rho c_v)_i on the
//! diagonal but the NEIGHBOURS' alpha on the off-diagonals, so dominance needs
//! alpha_i (C_l + C_r) >= C_l alpha_{i-1} + C_r alpha_{i+1}, which fails wherever
//! 1/(rho c_v) rises steeply across a face (the top of a convection zone, a steep
//! density drop) with beta_dt C/V >> 1.  The solve then undershoots and drives e_int
//! negative.  Solving instead for dT_i, with the heat capacity ON THE DIAGONAL,
//!     (rho c_v V)_i/beta_dt dT_i - sum_f C_f (dT_j - dT_i) = sum_f C_f (T*_j - T*_i),
//!     de_i = (rho c_v)_i dT_i   (the same frozen c_v),
//! is the SAME linearised problem -- it is the old row i multiplied through by
//! (rho c_v)_i -- but it is symmetric, an M-matrix, and diagonally dominant for ANY
//! alpha contrast, so dT cannot reach a new extremum and de cannot empty a cell.
//! Results change only by round-off wherever the old solve converged.
//!
//! HYDRO AND MHD.  In MHD u0(IEN) carries the magnetic energy too, so the frozen
//! internal energy subtracts 0.5|bcc0|^2 as well (MagEnergyCC); nothing else in the
//! routine changes, since what is solved for is the INCREMENT and what is written is
//! u0(IEN) += x.  The field is frozen over the step exactly as T*, c_v and K_f are.

void Conduction::ImplicitRadialUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                                      const Real beta_dt, const bool rt_on,
                                      const int rt_pass) {
  if (!rad_implicit_x1) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;
  const bool curvg = pmy_pack->pmesh->use_cubed_sphere
                     || pmy_pack->pmesh->use_spherical_polar;
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  auto &cosc_ = pmy_pack->pcoord->cos_cell;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &vol_ = pmy_pack->pcoord->volume;
  auto &area1_ = pmy_pack->pcoord->area.x1f;
  auto eos_ = eos;
  const bool gen = eos.IsGeneral();
  // hydro or MHD: the cached temperature, the gravitational potential and, in MHD, the
  // cell-centred field whose energy has to come out of u0(IEN) as well
  const bool ismhd = (my_block.compare("mhd") == 0);
  auto &wtemp_ = ismhd ? pmy_pack->pmhd->wtemp : pmy_pack->phydro->wtemp;
  auto &phicc_ = ismhd ? pmy_pack->pmhd->phicc0 : pmy_pack->phydro->phicc0;
  const bool etg = ismhd ? pmy_pack->pmhd->use_etotgrav
                         : pmy_pack->phydro->use_etotgrav;
  // bcc0 is the cell-centred form of the CURRENT b0 and, on the cubed sphere, is already
  // in the orthonormal frame -- see MagEnergyCC.  MHD::ImplicitConduction runs from the
  // stage before MHD::CT, exactly as the hydro one runs before the ghost exchange, so
  // b0 has not moved since the ConToPrim that filled it.  A zero-size dummy in hydro,
  // captured by the kernel and never read.
  DvceArray5D<Real> bcc_("imp_bcc_dummy", 1, 1, 1, 1, 1);
  if (ismhd) bcc_ = pmy_pack->pmhd->bcc0;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real eflx_unit = pres_unit*pmy_pack->punit->velocity_cgs();
  const Real met = rad_met, kfac = rad_kappa_fac;
  // the rad_tmax_kappa ceiling reaches the implicit radial solve too: it forms its
  // face conductivities with the same RadFaceKappa as the explicit operator
  const Real tmax = rad_tmax;   // temperature ceiling in kappa_rad only (0 = off)
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  const bool limit = rad_flux_limit;
  const Real ffac = rad_flim_fac;      // 4 sigma T^4, or 1 with rad_flim_legacy
  const Real krmax = rad_kappa_rmax;
  // the density gate, on the same face-averaged density the explicit x1 face uses
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  const Real sigma_sb = 5.670374419e-5;
  auto &wf = rad_w;
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  // rad_blend_use_2s: the ramp faces are prescribed-flux and leave this system
  // mode 1 ONLY: mode 2 keeps the ramp faces in the system and corrects them explicitly
  const bool use2s_ = (rad_blend_use_2s == 1) && taumode && blend_r && rad_f2s_ready;
  auto wrk = imp_wrk;
  auto iflag = imp_flag;
  auto irec = imp_rec;
  // ---- DIAGNOSTIC ONLY: the T-linearisation audit (rad_x1_verbose) -----------------
  const bool x1dbg_ = rad_x1_verbose && (rad_x1_every > 0) &&
                      (pmy_pack->pmesh->ncycle % rad_x1_every == 0);
  auto x1dg_ = imp_x1dg;
  if (x1dbg_) Kokkos::deep_copy(x1dg_, 0.0);
  // ---- the merged two-stream column solve (<problem>/rt_implicit_column) ------------
  // rt_on is passed true ONLY by two_stream_rt, which has just filled rt_col_res /
  // rt_col_jac / rt_col_dbdt for the CURRENT state.  Everything below is behind it, so
  // with the switch off not one expression of the original solve changes.
  const bool rtc_ = rt_on && rt_col_active && rt_col_alloc;
  const int rtpass_ = rt_pass;
  auto rtres_ = rtc_ ? rt_col_res : DvceArray4D<Real>("rtc_res_d", 1, 1, 1, 1);
  auto rtjac_ = rtc_ ? rt_col_jac : DvceArray5D<Real>("rtc_jac_d", 1, 1, 1, 1, 1);
  auto rtdbt_ = rtc_ ? rt_col_dbdt : DvceArray4D<Real>("rtc_dbt_d", 1, 1, 1, 1);
  auto rttn_ = rtc_ ? rt_col_tn : DvceArray4D<Real>("rtc_tn_d", 1, 1, 1, 1);
  auto rtdx_ = rtc_ ? rt_col_dtex : DvceArray4D<Real>("rtc_dx_d", 1, 1, 1, 1);
  const Real rtdtmax_ = rt_col_dtmax;
  const bool rtdbg_ = rtc_ && rt_col_verbose;
  auto rtdg_ = rtc_ ? rt_col_diag : DvceArray1D<Real>("rtc_dg_d", 8);
  if (rtc_) Kokkos::deep_copy(rtdg_, 0.0);
  // slot indices as plain locals: a static constexpr member would capture `this`
  const int e_ = IMPE, t_ = IMPT, al_ = IMPA, pr_ = IMPP;
  const int c_ = IMPC, cp_ = IMPCP, dp_ = IMPDP;
  Kokkos::deep_copy(iflag, 0);

  // THE COEFFICIENTS ARE FORMED IN THEIR OWN KERNELS, one thread per CELL.
  //
  // The solve itself is a column sweep -- one thread per (m,k,j) -- and there are only
  // ~1.2e4 columns on a rank, which is not enough work to fill a GPU. That did not matter
  // while the thread only walked a tridiagonal recurrence, but the first two sweeps are
  // the expensive ones: per cell an EintFromCons, an EOS Temperature (an ITERATIVE
  // inversion for the general EOS), a SpecificHeatCv and a Pressure, and per face a
  // RadFaceKappa table lookup. Profiled, "radimpx1" was 21.6% of the GPU time of a dhj
  // MHD run. Those two sweeps are perfectly parallel over cells and faces -- the state
  // sweep reads only u0 and the frozen caches, and the face sweep reads only what the
  // state sweep wrote at i-1 and i -- so they become two par_for kernels over the full
  // (m,k,j,i) range and the recurrence keeps the column kernel to itself. Every
  // expression is unchanged and in the same order, so the update is bitwise identical.
  // MEASURED on that run: "radimpx1" went 3179 ms -> 368 ms, 8.6x, and 21.6% -> 2.2% of
  // the GPU time, for no change in any output byte.
  //
  // ---- the FROZEN state of every cell: internal energy, temperature, pressure and
  // 1/(rho c_v).  The internal energy is extracted exactly as ConToPrim extracts it:
  // the gravitational term (etotgrav) and, on the cubed sphere, the kinetic energy
  // formed with the non-orthogonal metric (GnomonicEquiangleRaiseVel).  A cell whose
  // internal energy is not positive is marked with a negative 1/(rho c_v) and is
  // dropped from the system rather than handed to the EOS inversion.
  par_for("radimpx1_coef", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = u0(m,IDN,k,j,i);
    // the one shared extraction (utils/eint_from_cons.hpp): identical arithmetic in
    // identical order to what stood here, so this refactor is bit-for-bit a no-op
    const Real ei = EintFromCons(u0, m, k, j, i, cs_ ? cosc_(m,k,j) : 0.0, cs_, etg,
                                 etg ? phicc_(m,k,j,i) : 0.0,
                                 ismhd ? MagEnergyCC(bcc_,m,k,j,i) : 0.0);
    wrk(m,e_,k,j,i) = ei;
    wrk(m,t_,k,j,i) = 0.0;
    wrk(m,pr_,k,j,i) = 0.0;
    wrk(m,al_,k,j,i) = -1.0;
    if (!(ei > 0.0) || !(d > 0.0) || !isfinite(ei)) return;
    const Real tt = eos_.Temperature(d, ei, gen ? wtemp_(m,k,j,i) : -1.0);
    if (!(tt > 0.0) || !isfinite(tt)) return;
    const Real cv = eos_.SpecificHeatCv(d, ei, tt);
    if (!(cv > 0.0) || !isfinite(cv)) return;
    wrk(m,t_,k,j,i) = tt;
    wrk(m,pr_,k,j,i) = eos_.Pressure(d, ei, tt);
    wrk(m,al_,k,j,i) = 1.0/(d*cv);
  });

  // T^n, the temperature the OUTER iteration started from.  Passes 2..k must keep the
  // heat-capacity term anchored on it, or each pass would take another FULL backward-
  // Euler conduction step instead of correcting the one already taken.
  if (rtc_ && rtpass_ == 0) {
    auto wrk_tn = wrk;
    const int t_tn = t_;
    par_for("radimpx1_tn", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      rttn_(m,k,j,i) = wrk_tn(m,t_tn,k,j,i);
    });
  }

  // ---- the frozen face coefficient A_f K_f/dl_f.  The two boundary faces are
  // outside the system: they were added explicitly with the ghost states.
  par_for("radimpx1_face", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (i == is || i == ie+1) {
      wrk(m,c_,k,j,i) = 0.0;
      return;
    }
    // rad_blend_use_2s: a face inside the tau ramp carries a PRESCRIBED flux, added
    // explicitly by the face-flux kernel.  Zero conductance takes it out of the
    // tridiagonal coupling entirely, which is what "prescribed" means for this solve; the
    // rows on either side then see only their remaining faces and stay diagonally
    // dominant.  The mask must be the SAME test the flux kernel used, and it is: rad_w is
    // built once per stage, before both.
    if (use2s_) {
      const Real wv = wf(m,k,j,i);
      if (wv > 0.0 && wv < 1.0) {
        wrk(m,c_,k,j,i) = 0.0;
        return;
      }
    }
    const Real dx1c = size.d_view(m).dx1;
    Real ca = 0.0;
    const Real all = wrk(m,al_,k,j,i-1), alr = wrk(m,al_,k,j,i);
    if (all > 0.0 && alr > 0.0 && !(krmax > 0.0 && x1v_(m,i) > krmax)) {
      const Real tl = wrk(m,t_,k,j,i-1), tr = wrk(m,t_,k,j,i);
      const Real pf = 0.5*(wrk(m,pr_,k,j,i-1) + wrk(m,pr_,k,j,i));
      if (!(pf < pcut)) {
        const Real dl = curvg ? (x1v_(m,i) - x1v_(m,i-1)) : dx1c;
        const Real tk = 0.5*(tl + tr)*temp_unit;
        const Real rhof = 0.5*(u0(m,IDN,k,j,i-1) + u0(m,IDN,k,j,i))*dens_unit;
        const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                      krnT, krnP, krho, met, kfac, tmax);
        // the flux limiter, evaluated on the FROZEN gradient and then held fixed:
        // F = -kap g/sqrt(1 + (kap g/F_free)^2) linearises to a diffusion coefficient
        // kap/sqrt(1 + s^2) at fixed s, which is what keeps the system linear
        Real lf = 1.0;
        if (limit) {
          const Real ffree = ffac*sigma_sb*tk*tk*tk*tk;
          if (ffree > 0.0) {
            const Real fu = -kap*((tr - tl)/dl)*temp_unit/len_unit;
            lf = 1.0/sqrt(1.0 + SQR(fu/ffree));
          } else {
            lf = 0.0;
          }
        }
        Real wt = (taumode && blend_r) ? wf(m,k,j,i) : 1.0;
        if (gaterho > 0.0) wt *= RadGate(rhof, gaterho, gatedex);
        const Real af = curvg ? area1_(m,k,j,i) : 1.0;
        ca = wt*kap*lf*temp_unit/len_unit/eflx_unit*af/dl;
        if (!isfinite(ca) || ca < 0.0) ca = 0.0;
      }
    }
    wrk(m,c_,k,j,i) = ca;
  });

  const int nkj = (ke - ks + 1)*(je - js + 1);
  const int nj = (je - js + 1);
  Real maxviol = 0.0;
  int nfail = 0;
  int nclip = 0;
  Kokkos::parallel_reduce("radimpx1",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkj),
  KOKKOS_LAMBDA(const int &idx, Real &mviol, int &nbad, int &nclp) {
    const int m = idx/nkj;
    const int k = (idx - m*nkj)/nj + ks;
    const int j = (idx - m*nkj - (k - ks)*nj) + js;
    const Real dx1c = size.d_view(m).dx1;

    // ---- forward sweep of the Thomas algorithm on the TEMPERATURE increment
    // y_i = T_i - T*_i.  Row i, with D_i = (rho c_v V)_i/beta_dt = V_i/(beta_dt alpha_i):
    //   -A_i C_i y_{i-1} + (D_i + A_i C_i + A_{i+1} C_{i+1}) y_i
    //                                                      - A_{i+1} C_{i+1} y_{i+1}
    // = A_i C_i (T*_{i-1} - T*_i) + A_{i+1} C_{i+1} (T*_{i+1} - T*_i)
    // which is the old row on the energy increment multiplied through by (rho c_v)_i:
    // the same linearised problem, but symmetric, diagonally dominant for any alpha
    // contrast, and an M-matrix.  The energy increment comes back as de_i = y_i/alpha_i.
    // A cell dropped by the state sweep (alpha <= 0) has both of its face conductances
    // zero, so the row degenerates to y_i = 0 for any positive diagonal.
    for (int i=is; i<=ie; ++i) {
      const Real vi = curvg ? vol_(m,k,j,i) : dx1c;
      const Real cl = wrk(m,c_,k,j,i), cr = wrk(m,c_,k,j,i+1);
      const Real ac = wrk(m,al_,k,j,i);
      const Real dg = (ac > 0.0) ? vi/(beta_dt*ac) : 1.0;
      Real aa = -cl;
      Real bb = dg + cl + cr;
      Real cc = -cr;
      const Real tc = wrk(m,t_,k,j,i);
      const Real tm = (i > is) ? wrk(m,t_,k,j,i-1) : 0.0;
      const Real tp = (i < ie) ? wrk(m,t_,k,j,i+1) : 0.0;
      Real rhs = cr*(tp - tc) + cl*(tm - tc);
      // ---- the two-stream's nearest-neighbour linearisation, folded in ------------
      // R_i is the FULL explicit source at this state (everything >= 2 cells away, the
      // stellar beam and the tau-blend handover included), so the fixed point of the
      // outer iteration is the exact backward-Euler balance whatever the Jacobian gets
      // wrong.  J_ij = dR_i/dB_j dB_j/dT_j, converted from Kelvin to code temperature
      // and volume-integrated so it lives in the same units as the face conductances.
      // J_ii <= 0 (the cell's own emission) raises the diagonal and J_{i,i+-1} >= 0 (it
      // absorbs what its neighbour emits) lowers the off-diagonals, so the row stays an
      // M-matrix; the audit below counts any row where that fails.
      if (rtc_ && ac > 0.0) {
        const Real jm = (i > is) ? vi*rtjac_(m,0,k,j,i)*rtdbt_(m,k,j,i-1)*temp_unit : 0.0;
        const Real j0 = vi*rtjac_(m,1,k,j,i)*rtdbt_(m,k,j,i)*temp_unit;
        const Real jp = (i < ie) ? vi*rtjac_(m,2,k,j,i)*rtdbt_(m,k,j,i+1)*temp_unit : 0.0;
        if (isfinite(jm) && isfinite(j0) && isfinite(jp)) {
          aa -= jm;
          bb -= j0;
          cc -= jp;
        }
        rhs += vi*rtres_(m,k,j,i);
        // A neighbour left OUT of the system -- a thin cell that relaxed itself, or the
        // ghost beyond the first/last row -- has dbdt = 0 above, so its exchange term
        // dropped out of the matrix.  Put back the part of it that is already KNOWN: the
        // change in that cell's Planck function, dB, which its own relaxation has just
        // applied.  R_i was formed before that change, so without this the exchange
        // across a thick/thin interface would be counted with the stale neighbour.
        if (i > is) rhs += vi*rtjac_(m,0,k,j,i)*rtdx_(m,k,j,i-1);
        if (i < ie) rhs += vi*rtjac_(m,2,k,j,i)*rtdx_(m,k,j,i+1);
        if (rtpass_ > 0) rhs -= dg*(tc - rttn_(m,k,j,i));
        const Real offs = fabs(aa) + fabs(cc);
        if (bb > 0.0) {
          const Real rat = offs/bb;
          if (rat > 1.0) {
            Kokkos::atomic_fetch_add(&rtdg_(0), 1.0);
            Kokkos::atomic_max(&rtdg_(1), rat);
          }
        } else {
          Kokkos::atomic_fetch_add(&rtdg_(0), 1.0);
        }
      }
      if (rtdbg_ && m == 0 && k == ks && j == js) {
        const Real j0d = vi*rtjac_(m,1,k,j,i)*rtdbt_(m,k,j,i)*temp_unit;
        Kokkos::printf("### rtcol_row i=%d V=%.5e dg=%.6e cl=%.4e cr=%.4e VJii=%.6e "
                       "aa=%.6e bb=%.6e cc=%.6e rhs=%.6e VR=%.6e T*=%.6e tu=%.4e\n",
                       i, vi, dg, cl, cr, j0d, aa, bb, cc, rhs, vi*rtres_(m,k,j,i),
                       tc, temp_unit);
      }
      if (i == is) {
        wrk(m,cp_,k,j,i) = cc/bb;
        wrk(m,dp_,k,j,i) = rhs/bb;
      } else {
        const Real den = bb - aa*wrk(m,cp_,k,j,i-1);
        wrk(m,cp_,k,j,i) = cc/den;
        wrk(m,dp_,k,j,i) = (rhs - aa*wrk(m,dp_,k,j,i-1))/den;
      }
    }

    // ---- back substitution, write-back, and the per-column conservation residual.
    // Any positivity clipping is CONSERVATIVE: the energy a clipped cell is not allowed
    // to give up is taken from the neighbour it is most strongly coupled to, so
    // sum_i V_i de_i is preserved to round-off.  The sweep runs downwards, so a debt
    // owed to i-1 is carried in `pend` (a volume-integrated energy) and paid on the next
    // iteration, while a debt owed to i+1 is applied to u0 directly -- one thread owns
    // the whole column, so there is no race.
    Real xnext = 0.0, csum = 0.0, cabs = 0.0, pend = 0.0, eprev = 0.0;
    // the merged solve is NOT source-free: the column exchanges energy with the two
    // boundaries through the radiation field, so sum_i V_i de_i must come out equal to
    // beta_dt times the LINEARISED two-stream source summed over the column, not to
    // zero.  rtexp accumulates exactly that; it is identically zero with the switch off,
    // so the residual reported below is the same number it always was.
    Real rtexp = 0.0;
    for (int i=ie; i>=is; --i) {
      const Real vi = curvg ? vol_(m,k,j,i) : dx1c;
      const Real ac = wrk(m,al_,k,j,i);
      Real y = wrk(m,dp_,k,j,i) - wrk(m,cp_,k,j,i)*xnext;
      bool bad = !isfinite(y);
      if (bad) y = 0.0;
      // ---- rt_impl_dtmax: bound the linearised Newton step -----------------------
      // B ~ T^4 is CONVEX, so a row that needs heating is handed a tangent whose root
      // lies past the true one, and in an optically thin cell (E -> 0) that overshoot is
      // unbounded -- the same divergence the per-cell relaxation avoids by relaxing
      // toward the true fixed point.  rt_impl_tau_min keeps those cells out of the
      // system; this is the belt for the ones that are in it.
      if (rtc_ && rtdtmax_ > 0.0) {
        const Real tcb = wrk(m,t_,k,j,i);
        if (tcb > 0.0) {
          const Real ycap = rtdtmax_*tcb;
          if (fabs(y) > ycap) {
            Kokkos::atomic_fetch_add(&rtdg_(4), 1.0);
            Kokkos::atomic_max(&rtdg_(5), fabs(y)/tcb);
            y = (y > 0.0) ? ycap : -ycap;
          }
        }
      }
      if (rtdbg_ && m == 0 && k == ks && j == js) {
        Kokkos::printf("### rtcol_sol i=%d dT=%.6e T*=%.6e dT/T=%.4e de=%.6e\n",
                       i, y, wrk(m,t_,k,j,i),
                       (wrk(m,t_,k,j,i) > 0.0) ? y/wrk(m,t_,k,j,i) : 0.0,
                       (wrk(m,al_,k,j,i) > 0.0) ? y/wrk(m,al_,k,j,i) : 0.0);
      }
      // DIAGNOSTIC ONLY: park the solved increment where the audit below can read it.
      // dp_ at this i has already been consumed by the line above and is never read
      // again, so this is a dead slot from here on.
      if (x1dbg_) wrk(m,dp_,k,j,i) = y;
      const Real yprev = xnext;             // y_{i+1}, already solved
      xnext = y;
      if (rtc_ && ac > 0.0) {
        // R_i + J_ii y_i + J_{i,i+1} y_{i+1}, and cell i+1's absorption of THIS cell's
        // emission, V_{i+1} J_{i+1,i} y_i -- the only term of row i+1 still outstanding
        rtexp += vi*rtres_(m,k,j,i);
        rtexp += vi*rtjac_(m,1,k,j,i)*rtdbt_(m,k,j,i)*temp_unit*y;
        // the same known dB of a neighbour left out of the system
        if (i > is) rtexp += vi*rtjac_(m,0,k,j,i)*rtdx_(m,k,j,i-1);
        if (i < ie) rtexp += vi*rtjac_(m,2,k,j,i)*rtdx_(m,k,j,i+1);
        if (i < ie) {
          const Real vp = curvg ? vol_(m,k,j,i+1) : dx1c;
          rtexp += vi*rtjac_(m,2,k,j,i)*rtdbt_(m,k,j,i+1)*temp_unit*yprev;
          rtexp += vp*rtjac_(m,0,k,j,i+1)*rtdbt_(m,k,j,i)*temp_unit*y;
        }
      }
      Real x = (ac > 0.0) ? y/ac : 0.0;
      if (!isfinite(x)) {
        x = 0.0;
        bad = true;
      }
      x += pend/vi;
      pend = 0.0;
      const Real es = wrk(m,e_,k,j,i);
      if (x != 0.0 && !((es + x) > 0.0)) {
        // with the M-matrix on T this should be unreachable from the solve itself; it
        // can still be reached by a debt handed down from i+1.  Clip to a positive
        // sliver and move the difference onto the stiffest neighbouring face.
        const Real cl = wrk(m,c_,k,j,i), cr = wrk(m,c_,k,j,i+1);
        const Real xn = -(1.0 - 1.0e-10)*es;
        const Real amt = (xn - x)*vi;   // energy kept here, owed by a neighbour
        // the stiffest neighbour is asked first and the other one second.  Downwards
        // the debt is safe unconditionally -- cell i-1 has not been tested yet, so if
        // it cannot afford it either it clips in turn and passes the rest on -- while
        // upwards it has to fit in what cell i+1 has left, because that cell is done.
        bool paid = false;
        for (int p = 0; p < 2 && !paid; ++p) {
          if ((cr >= cl) == (p == 0)) {
            const Real vp = curvg ? vol_(m,k,j,i+1) : dx1c;
            const Real take = amt/vp;
            if (es > 0.0 && cr > 0.0 && i < ie && take < (1.0 - 1.0e-10)*eprev) {
              u0(m,IEN,k,j,i+1) -= take;
              eprev -= take;
              csum -= amt;
              cabs += fabs(amt);
              paid = true;
            }
          } else {
            if (es > 0.0 && cl > 0.0 && i > is) {
              pend = -amt;
              paid = true;
            }
          }
        }
        if (paid) {
          x = xn;
          ++nclp;
        } else {
          // nowhere to put it: an isolated cell, a cell with no positive energy left,
          // or a neighbour that cannot afford the debt.  This is the only remaining
          // non-conservative path and it is counted as a fallback.
          x = 0.0;
          bad = true;
        }
      }
      if (bad) {
        ++nbad;
        if (Kokkos::atomic_fetch_add(&iflag(0), 1) == 0) {
          irec(0) = static_cast<Real>(m);
          irec(1) = static_cast<Real>(k);
          irec(2) = static_cast<Real>(j);
          irec(3) = static_cast<Real>(i);
          irec(4) = x1v_(m,i);
          irec(5) = wrk(m,t_,k,j,i)*temp_unit;
          irec(6) = u0(m,IDN,k,j,i)*dens_unit;
          irec(7) = y;
        }
      }
      if (rtc_ && rtres_(m,k,j,i) != 0.0) {
        // the L-infinity energy increment the tridiagonal itself applies to a row that
        // carries a two-stream source: the THICK half of the per-pass contraction
        Kokkos::atomic_max(&rtdg_(6), fabs(x));
        const Real tcc = wrk(m,t_,k,j,i);
        if (tcc > 0.0) Kokkos::atomic_max(&rtdg_(7), fabs(y)/tcc);
      }
      u0(m,IEN,k,j,i) += x;
      eprev = es + x;
      csum += vi*x;
      cabs += fabs(vi*x);
    }
    if (rtc_) {
      rtexp *= beta_dt;
      Kokkos::atomic_fetch_add(&rtdg_(2), csum);
      Kokkos::atomic_fetch_add(&rtdg_(3), rtexp);
    }
    if (cabs > 0.0) mviol = fmax(mviol, fabs(csum - rtexp)/cabs);
  }, Kokkos::Max<Real>(maxviol), nfail, nclip);

  // ==================================================================================
  // DIAGNOSTIC ONLY (rad_x1_verbose).  How wrong is the T-linearisation?
  //
  // The solve froze the face conductance C_f = A_f K_f/dl at the OLD state and solved a
  // linear system for dT.  Three errors are measured, per interior face, over every
  // column on the rank:
  //   (1) the K-nonlinearity: recompute C_f with the SAME kappa table / limiter / blend
  //       at T + dT (and p scaled by T_new/T*, rho frozen) and compare the flux
  //       C_f(T+dT) * dT_grad with the flux the linear solve actually applied,
  //       C_f(T) * dT_grad;
  //   (2) the u-form check: radiative diffusion is EXACTLY linear in u = T^4 at frozen
  //       opacity, F = -(ac/3 kappa rho) du/dz = -(K/(4T^3)) du/dz.  Compare
  //       (K_f/(4 T_f^3)) (u_i - u_j) with the T-form K_f (T_i - T_j) the solve used.
  //   (3) the column-integrated |dF| of (1), normalised by the flux through the lowest
  //       interior face.
  // plus max |dT/T| and where it sits (i, cell tau measured down from the top, T[K]).
  // Two passes: pass 0 takes the maxima, pass 1 records where each maximum sits.
  // ==================================================================================
  if (x1dbg_) {
    for (int pass = 0; pass < 2; ++pass) {
      const int pss = pass;
      Kokkos::parallel_for("radimpx1_dbg",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkj),
      KOKKOS_LAMBDA(const int &idx) {
        const int m = idx/nkj;
        const int k = (idx - m*nkj)/nj + ks;
        const int j = (idx - m*nkj - (k - ks)*nj) + js;
        const Real dx1c = size.d_view(m).dx1;
        Real tau = 0.0;
        Real colabs = 0.0, fbot = 0.0;
        for (int i=ie; i>=is; --i) {
          const Real tc = wrk(m,t_,k,j,i);
          const Real ac = wrk(m,al_,k,j,i);
          const Real dTc = wrk(m,dp_,k,j,i);
          // cell optical depth accumulated from the top of the column
          if (tc > 0.0 && ac > 0.0) {
            const Real rho = u0(m,IDN,k,j,i)*dens_unit;
            const Real kc = RadFaceKappa(tc*temp_unit, wrk(m,pr_,k,j,i)*pres_unit, rho,
                                         ktab, krt, krlT, krlP, krnT, krnP, krho,
                                         met, kfac, tmax);
            const Real dlc = (curvg && i < ie) ? (x1v_(m,i+1) - x1v_(m,i)) : dx1c;
            // RadFaceKappa returns the radiative CONDUCTIVITY 16 sigma T^3/(3 kfac
            // kappa_R rho), so the true opacity x density is 16 sigma T^3/(3 kfac kc)
            const Real tk3 = tc*temp_unit;
            if (kc > 0.0) {
              tau += 16.0*5.670374419e-5*tk3*tk3*tk3/(3.0*kfac*kc)*dlc*len_unit;
            }
          }
          if (tc > 0.0 && ac > 0.0 && isfinite(dTc)) {
            const Real r = fabs(dTc/tc);
            if (pss == 0) {
              Kokkos::atomic_max(&x1dg_(0), r);
              Kokkos::atomic_add(&x1dg_(11), r);
              Kokkos::atomic_add(&x1dg_(12), 1.0);
            } else if (r == x1dg_(0) && r > 0.0) {
              x1dg_(1) = static_cast<Real>(i);
              x1dg_(2) = tau;
              x1dg_(3) = tc*temp_unit;
            }
          }
          // ---- the face between cell i and cell i+1 (face index i+1) --------------
          if (i < ie) {
            const int f = i + 1;
            const Real cf = wrk(m,c_,k,j,f);
            const Real tl = wrk(m,t_,k,j,f-1), tr = wrk(m,t_,k,j,f);
            const Real dl_l = wrk(m,dp_,k,j,f-1), dl_r = wrk(m,dp_,k,j,f);
            if (cf > 0.0 && tl > 0.0 && tr > 0.0 && isfinite(dl_l) && isfinite(dl_r)) {
              const Real tln = tl + dl_l, trn = tr + dl_r;
              const Real flin = cf*(tln - trn);
              if (tln > 0.0 && trn > 0.0) {
                // (1) recompute the conductance at the NEW temperature, same recipe
                const Real dl = curvg ? (x1v_(m,f) - x1v_(m,f-1)) : dx1c;
                const Real tkn = 0.5*(tln + trn)*temp_unit;
                const Real rhof = 0.5*(u0(m,IDN,k,j,f-1) + u0(m,IDN,k,j,f))*dens_unit;
                // pressure carried along with the temperature at frozen density
                const Real pfn = 0.5*(wrk(m,pr_,k,j,f-1)*(tln/tl)
                                    + wrk(m,pr_,k,j,f)*(trn/tr));
                const Real kapn = RadFaceKappa(tkn, pfn*pres_unit, rhof, ktab, krt,
                                               krlT, krlP, krnT, krnP, krho, met,
                                               kfac, tmax);
                Real lfn = 1.0;
                if (limit) {
                  const Real ffree = ffac*sigma_sb*tkn*tkn*tkn*tkn;
                  if (ffree > 0.0) {
                    const Real fu = -kapn*((trn - tln)/dl)*temp_unit/len_unit;
                    lfn = 1.0/sqrt(1.0 + SQR(fu/ffree));
                  } else {
                    lfn = 0.0;
                  }
                }
                Real wtn = (taumode && blend_r) ? wf(m,k,j,f) : 1.0;
                if (gaterho > 0.0) wtn *= RadGate(rhof, gaterho, gatedex);
                const Real afn = curvg ? area1_(m,k,j,f) : 1.0;
                Real cfn = wtn*kapn*lfn*temp_unit/len_unit/eflx_unit*afn/dl;
                if (!isfinite(cfn) || cfn < 0.0) cfn = 0.0;
                const Real fnl = cfn*(tln - trn);
                const Real den = fmax(fabs(flin), 1.0e-300);
                const Real e1 = fabs(fnl - flin)/den;
                // (2) the exact u-form flux at the SAME frozen opacity
                const Real tf3 = 0.5*(tl + tr); // T_f of the frozen state
                const Real fu_ = (tf3 > 0.0)
                    ? cf/(4.0*tf3*tf3*tf3)*(tln*tln*tln*tln - trn*trn*trn*trn) : flin;
                const Real e2 = fabs(fu_ - flin)/den;
                colabs += fabs(fnl - flin);
                if (f == is + 1) fbot = fabs(flin);
                if (pss == 0) {
                  Kokkos::atomic_max(&x1dg_(4), e1);
                  Kokkos::atomic_max(&x1dg_(7), e2);
                  Kokkos::atomic_max(&x1dg_(15), fabs(flin));
                  Kokkos::atomic_add(&x1dg_(13), e1);
                  Kokkos::atomic_add(&x1dg_(14), 1.0);
                } else {
                  if (e1 == x1dg_(4) && e1 > 0.0) {
                    x1dg_(5) = static_cast<Real>(f); x1dg_(6) = tau;
                  }
                  if (e2 == x1dg_(7) && e2 > 0.0) {
                    x1dg_(8) = static_cast<Real>(f); x1dg_(9) = tau;
                  }
                }
              }
            }
          }
        }
        if (pss == 0 && fbot > 0.0) Kokkos::atomic_max(&x1dg_(10), colabs/fbot);
      });
    }
    if (global_variable::my_rank == 0 && x1dbg_lines < 4000) {
      ++x1dbg_lines;
      auto hx = Kokkos::create_mirror_view(imp_x1dg);
      Kokkos::deep_copy(hx, imp_x1dg);
      const Real nc = (hx(12) > 0.0) ? hx(12) : 1.0;
      const Real nf = (hx(14) > 0.0) ? hx(14) : 1.0;
      std::cout << "### rad_x1_lin cycle " << pmy_pack->pmesh->ncycle
                << " t= " << pmy_pack->pmesh->time
                << " dt= " << beta_dt
                << " maxdToT= " << hx(0) << " @i= " << static_cast<int>(hx(1))
                << " tau= " << hx(2) << " T= " << hx(3)
                << " meandToT= " << hx(11)/nc
                << " | maxKerr= " << hx(4) << " @i= " << static_cast<int>(hx(5))
                << " tau= " << hx(6) << " meanKerr= " << hx(13)/nf
                << " | maxUerr= " << hx(7) << " @i= " << static_cast<int>(hx(8))
                << " tau= " << hx(9)
                << " | colint= " << hx(10) << " Fmax= " << hx(15) << std::endl;
    }
  }

  // <problem>/nan_report: the conservation residual of the tridiagonal solve, and any
  // cell that had to fall back.  Nothing is printed unless the switch is on, except a
  // fallback, which is always worth a line.
  if (imp_lines < 400) {
    const bool bad = (nfail > 0);
    const bool tell = nan_report && (imp_lines < 20 || maxviol > 1.0e-10 || nclip > 0);
    if (bad || (tell && global_variable::my_rank == 0)) {
      ++imp_lines;
      std::cout << "### rad_implicit_x1 rank " << global_variable::my_rank
                << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time
                << ": max |sum V de|/sum V|de| = " << maxviol
                << ", fallback cells = " << nfail
                << ", conservative clips = " << nclip << std::endl;
      if (bad) {
        auto hr = Kokkos::create_mirror_view(irec);
        Kokkos::deep_copy(hr, irec);
        const int mb = static_cast<int>(hr(0));
        std::cout << "    first fallback (m,k,j,i) = (" << mb << ","
                  << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
                  << static_cast<int>(hr(3)) << ") gid = " << (pmy_pack->gids + mb)
                  << " r = " << hr(4) << " T* = " << hr(5) << " rho = " << hr(6)
                  << " de = " << hr(7) << std::endl;
      }
    }
  }
  // ---- the M-matrix audit of the merged rows, and the column energy budget ----------
  if (rtc_) {
    auto hg = Kokkos::create_mirror_view(rt_col_diag);
    Kokkos::deep_copy(hg, rt_col_diag);
    const bool viol = (hg(0) > 0.0);
    if (viol || (rt_col_verbose && rt_col_lines < 60) ||
        (nan_report && rt_col_lines < 20)) {
      if (global_variable::my_rank == 0 && rt_col_lines < 200) {
        ++rt_col_lines;
        const Real den = fabs(hg(3)) + fabs(hg(2));
        std::cout << "### rt_implicit_column rank " << global_variable::my_rank
                  << " cycle " << pmy_pack->pmesh->ncycle
                  << " pass " << rtpass_
                  << ": M-matrix violations = " << static_cast<int64_t>(hg(0))
                  << ", worst |offdiag|/diag = " << hg(1)
                  << ", sum V de = " << hg(2) << ", expected = " << hg(3)
                  << ", rel = " << ((den > 0.0) ? fabs(hg(2) - hg(3))/den : 0.0)
                  << ", dT caps = " << static_cast<int64_t>(hg(4))
                  << ", worst |dT|/T = " << hg(5)
                  << ", max|de_thick| = " << hg(6)
                  << ", max|dT/T|_thick = " << hg(7)
                  << std::endl;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! n void Conduction::EnableRTColumn
//! rief allocate the four arrays the merged two-stream column solve exchanges with
//! two_stream_rt.  Called once, from the RT header, the first time the switch is seen.

void Conduction::EnableRTColumn() {
  if (rt_col_alloc) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng) : 1;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng) : 1;
  const int nmb = pmy_pack->nmb_thispack;
  Kokkos::realloc(rt_col_res, nmb, n3, n2, n1);
  Kokkos::realloc(rt_col_jac, nmb, 3, n3, n2, n1);
  Kokkos::realloc(rt_col_dbdt, nmb, n3, n2, n1);
  Kokkos::realloc(rt_col_tn, nmb, n3, n2, n1);
  Kokkos::realloc(rt_col_dtex, nmb, n3, n2, n1);
  Kokkos::realloc(rt_col_diag, 8);
  Kokkos::deep_copy(rt_col_res, 0.0);
  Kokkos::deep_copy(rt_col_jac, 0.0);
  Kokkos::deep_copy(rt_col_dbdt, 0.0);
  Kokkos::deep_copy(rt_col_tn, 0.0);
  Kokkos::deep_copy(rt_col_dtex, 0.0);
  Kokkos::deep_copy(rt_col_diag, 0.0);
  rt_col_active = true;
  rt_col_alloc = true;
  if (!rad_implicit_x1) {
    std::cout << "### FATAL ERROR in Conduction::EnableRTColumn: "
              << "<problem>/rt_implicit_column needs rad_implicit_x1 = true."
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsotropicHeatFluxConstCond()
//! \brief Adds isotropic heat flux computed using constant conductivity to face-centered
//! fluxes of conserved variables

void Conduction::AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w0,
    const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  Real gm1 = eos.gamma-1.0;
  // General EOS: temperature and pressure were evaluated once per cell in ConsToPrim.
  // Reading the cached values matters here -- T(d,e) is a root find, and the stencils
  // below touch up to ~30 neighbouring cells per face.
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  Real &kappa_ = kappa_iso;

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("conduct1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtempdx;
    if (gen) {
      dtempdx = (wtemp_(m,k,j,i) - wtemp_(m,k,j,i-1))
                / size.d_view(m).dx1;
    } else {
      dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                 - w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1))
                * gm1 / size.d_view(m).dx1;
    }
    flx1(m,IEN,k,j,i) -= kappa_ * dtempdx;
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("conduct2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtempdx;
    if (gen) {
      dtempdx = (wtemp_(m,k,j,i) - wtemp_(m,k,j-1,i))
                / size.d_view(m).dx2;
    } else {
      dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                 - w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i))
                * gm1 / size.d_view(m).dx2;
    }
    flx2(m,IEN,k,j,i) -= kappa_ * dtempdx;
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("conduct3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtempdx;
    if (gen) {
      dtempdx = (wtemp_(m,k,j,i) - wtemp_(m,k-1,j,i))
                / size.d_view(m).dx3;
    } else {
      dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                 - w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i))
                * gm1 / size.d_view(m).dx3;
    }
    flx3(m,IEN,k,j,i) -= kappa_ * dtempdx;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TempDependentHeatFlux()
//! \brief Adds heat flux to face-centered fluxes of conserved variables with
//! temperature-dependent conductivity

void Conduction::AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w0,
    const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
/*
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool &sat_hflux_ = sat_hflux;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  Real gm1 = eos.gamma-1.0;
  // General EOS: temperature and pressure were evaluated once per cell in ConsToPrim.
  // Reading the cached values matters here -- T(d,e) is a root find, and the stencils
  // below touch up to ~30 neighbouring cells per face.
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  Real kappaceil = kappa_ceiling;
  Real temp_unit = pmy_pack->punit->temperature_cgs();
  Real kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                    pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("conduct1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
    Real temp_r = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    Real pres_l = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
    Real pres_r = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx1 = (temp_r-temp_l)/size.d_view(m).dx1;
    Real hflx = kappaf*dtempdx1;
    // Saturation of thermal conduction by harmonic mean
    if (sat_hflux_) {
      Real dtempdx2 = 0.0, dtempdx3 = 0.0;
      if (multi_d) {
        temp_ll = (gen ? wtemp_(m,k,j-1,i-1)
                       : w0(m,IEN,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1)*gm1);
        temp_lr = (gen ? wtemp_(m,k,j+1,i-1)
                       : w0(m,IEN,k,j+1,i-1)/w0(m,IDN,k,j+1,i-1)*gm1);
        temp_rl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
        temp_rr = (gen ? wtemp_(m,k,j+1,i) : w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)*gm1);
        dtempdx2 = VanLeerLimiter4State(temp_rr-temp_r,temp_r-temp_rl,
                                        temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx2;
      }
      if (three_d) {
        temp_ll = (gen ? wtemp_(m,k-1,j,i-1)
                       : w0(m,IEN,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1)*gm1);
        temp_lr = (gen ? wtemp_(m,k+1,j,i-1)
                       : w0(m,IEN,k+1,j,i-1)/w0(m,IDN,k+1,j,i-1)*gm1);
        temp_rl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
        temp_rr = (gen ? wtemp_(m,k+1,j,i) : w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)*gm1);
        dtempdx3 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                              temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx3;
      }
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx1(m,IEN,k,j,i) -= hflx;
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("conduct2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = 0.0, temp_r = 0.0, pres_l = 0.0, pres_r = 0.0;
    Real temp_ll = 0.0, temp_lr = 0.0, temp_rl = 0.0, temp_rr = 0.0;
    temp_l = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
    temp_r = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    pres_l = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
    pres_r = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx2 = (temp_r-temp_l)/size.d_view(m).dx2;
    Real hflx = kappaf*dtempdx2;
    // Saturation of thermal conduction
    if (sat_hflux_) {
      Real dtempdx1 = 0.0, dtempdx3 = 0.0;
      temp_ll = (gen ? wtemp_(m,k,j-1,i-1) : w0(m,IEN,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1)*gm1);
      temp_lr = (gen ? wtemp_(m,k,j-1,i+1) : w0(m,IEN,k,j-1,i+1)/w0(m,IDN,k,j-1,i+1)*gm1);
      temp_rl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
      temp_rr = (gen ? wtemp_(m,k,j,i+1) : w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)*gm1);
      dtempdx1 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx1;
      if (three_d) {
        temp_ll = (gen ? wtemp_(m,k-1,j-1,i)
                       : w0(m,IEN,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i)*gm1);
        temp_lr = (gen ? wtemp_(m,k+1,j-1,i)
                       : w0(m,IEN,k+1,j-1,i)/w0(m,IDN,k+1,j-1,i)*gm1);
        temp_rl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
        temp_rr = (gen ? wtemp_(m,k+1,j,i) : w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)*gm1);
        dtempdx3 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                              temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx3;
      }
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx2(m,IEN,k,j,i) -= hflx;
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("conduct3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = 0.0, temp_r = 0.0, pres_l = 0.0, pres_r = 0.0;
    Real temp_ll = 0.0, temp_lr = 0.0, temp_rl = 0.0, temp_rr = 0.0;
    temp_l = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
    temp_r = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    pres_l = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
    pres_r = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx3 = (temp_r-temp_l)/size.d_view(m).dx3;
    Real hflx = kappaf*dtempdx3;
    // Saturation of thermal conduction
    if (sat_hflux_) {
      Real dtempdx1 = 0.0, dtempdx2 = 0.0;
      temp_ll = (gen ? wtemp_(m,k-1,j,i-1) : w0(m,IEN,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1)*gm1);
      temp_lr = (gen ? wtemp_(m,k-1,j,i+1) : w0(m,IEN,k-1,j,i+1)/w0(m,IDN,k-1,j,i+1)*gm1);
      temp_rl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
      temp_rr = (gen ? wtemp_(m,k,j,i+1) : w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)*gm1);
      dtempdx1 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx1;
      temp_ll = (gen ? wtemp_(m,k-1,j-1,i) : w0(m,IEN,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i)*gm1);
      temp_lr = (gen ? wtemp_(m,k-1,j+1,i) : w0(m,IEN,k-1,j+1,i)/w0(m,IDN,k-1,j+1,i)*gm1);
      temp_rl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
      temp_rr = (gen ? wtemp_(m,k,j+1,i) : w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)*gm1);
      dtempdx2 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx2;
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx3(m,IEN,k,j,i) -= hflx;
  });

*/
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::NewTimeStep()
//! \brief Compute new time step for thermal conduction.

void Conduction::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  // at initialisation the blend weights do not exist yet; without them this timestep is
  // unconstrained and the first cycle runs the operator far past its stability limit
  if (rad_tau_mode && !rad_w_built) {
    BuildRadWeights(w0, eos_data);
  }
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
//  if (sat_hflux == true) {
//    dtnew = static_cast<Real>(std::numeric_limits<float>::max());
//    return;
//  }

  // set flag for Spitzer conductivity
  bool spitzer = false;
  if ((iso_cond_type.compare("spitzer") == 0) ||
      (iso_cond_type.compare("spitzer_limited") == 0)) {
    spitzer = true;
  }
  Real limit_ = kappa_iso_limit;
  Real temp_unit=0.0, kappa_unit=0.0, pres_unit=1.0, dens_unit=1.0;
  Real vel_unit = 1.0;
  const bool radiative = (iso_cond_type.compare("radiative") == 0);
  const Real met = rad_met, kfac = rad_kappa_fac;
  const Real tmax = rad_tmax;   // temperature ceiling in kappa_rad only (0 = off)
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  const bool limit = rad_flux_limit;
  const Real ffac = rad_flim_fac;      // 4 sigma T^4, or 1 with rad_flim_legacy
  // the radial operator is unconditionally stable when it is solved implicitly, so it
  // carries no timestep constraint; x2/x3 are still explicit and still do
  // rad_sts_all removes it in the same way, by putting the x1 faces in the RKL1 loop
  const bool impx1 = rad_implicit_x1 || rad_sts_all;
  // ...and the transverse operator carries no constraint either when rad_cap_ang caps
  // every face, or when rad_implicit_ang solves it implicitly, so dt2 and dt3 go with it
  const bool capa = (rad_cap_ang > 0.0) || rad_implicit_ang;
  auto &wf = rad_w;
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;   // table axis is log rho, not log p
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;

  if (spitzer || radiative) {
    temp_unit = pmy_pack->punit->temperature_cgs();
    kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                 pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();
    pres_unit = pmy_pack->punit->pressure_cgs();
    dens_unit = pmy_pack->punit->density_cgs();
    vel_unit = pmy_pack->punit->velocity_cgs();
  }
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &dx1_ = pmy_pack->pcoord->dx1;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  const bool cs = pmy_pack->pmesh->use_cubed_sphere && rad_cs_exact;
  auto &sinc_ = pmy_pack->pcoord->sin_cell;
  // the inert region carries no flux, so it must not carry a constraint either
  const Real krmax = rad_kappa_rmax;
  // ...and neither does gas the density gate has made inert: the face fluxes it feeds
  // are all scaled by G, so the cell's constraint scales with G too
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  auto &x1v_n = pmy_pack->pcoord->x1v;

  // capture variables for kernel
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  auto &w0_ = w0;
  auto &multi_d = pmy_pack->pmesh->multi_d;
  auto &three_d = pmy_pack->pmesh->three_d;
  auto &size = pmy_pack->pmb->mb_size;
  Real gm1 = eos_data.gamma-1.0;
  // General EOS: temperature and pressure were evaluated once per cell in ConsToPrim.
  // Reading the cached values matters here -- T(d,e) is a root find, and the stencils
  // below touch up to ~30 neighbouring cells per face.
  const bool gen = eos_data.IsGeneral();
  auto eos_ = eos_data;   // by-value copy, capturable in the device lambda
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  Real kappa0 = kappa_iso;

  // find smallest timestep for thermal conduction in each cell
  // Note loop over all cells needed even for constant conductivity
  // MinLOC, not Min: when this timestep collapses the only question that matters is
  // WHICH cell did it, and reconstructing that afterwards from a dump means redoing the
  // opacity lookup, the tau blend and the flux limiter outside the code.  The location
  // rides along for free.
  // per-cycle diagnostic (problem/diag_gid): locals only, never `this`
  const bool diag_ = diag;
  auto cdg = cond_diag;
  const Real dt_huge = static_cast<Real>(std::numeric_limits<float>::max());

  Kokkos::ValLocScalar<Real, int> mloc;
  Kokkos::parallel_reduce("cond_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Kokkos::ValLocScalar<Real, int> &mres) {
    Real min_dt = mres.val;
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    // the diagnostic slots for cells that take one of the early returns below stay at
    // zero, which is what "this cell puts no constraint on the timestep" means
    if (diag_) {
      cdg(m,3,k,j,i) = 0.0;
      cdg(m,4,k,j,i) = 0.0;
      cdg(m,5,k,j,i) = 0.0;
    }

    Real kappa_ = kappa0;
    Real wmax = 1.0;
    if (spitzer) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      kappa_ = TempDepKappa(temp*temp_unit, limit_)/kappa_unit;
    } else if (radiative) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      Real pres = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      if (pres < pcut) return;   // no flux above the cut: no constraint
      if (krmax > 0.0 && x1v_n(m,i) > krmax) return;     // radiatively inert corona
      const Real tkap = KappaTemp(temp*temp_unit, tmax);
      kappa_ = (ktab
          ? RadiativeKappaKR(tkap, w0_(m,IDN,k,j,i)*dens_unit, kfac,
                             RosselandTable(krt, krlT, krlP, krnT, krnP, tkap,
                                            krho ? w0_(m,IDN,k,j,i)*dens_unit
                                                 : pres*pres_unit))
          : RadiativeKappa(tkap, pres*pres_unit, w0_(m,IDN,k,j,i)*dens_unit,
                           met, kfac))/kappa_unit;
      if (gaterho > 0.0) {
        kappa_ *= RadGate(w0_(m,IDN,k,j,i)*dens_unit, gaterho, gatedex);
      }
      // the blend: the face flux is w*kappa*grad T, so the explicit limit is on w*kappa,
      // and a cell whose faces carry no weight carries no constraint
      if (taumode) {
        wmax = fmax(wf(m,k,j,i), wf(m,k,j,i+1));
        if (wmax == 0.0 && blend_r) return;
      }
    }

    // the heat diffusion time is dx^2 rho c_v / kappa. For an ideal gas c_v = 1/(gamma-1)
    // which is what the rho/gm1 below amounts to; a general EOS has c_v(d,e), and it can
    // be an order of magnitude larger inside an ionization zone, so this is not a
    // cosmetic substitution -- it directly sets the conduction-limited timestep.
    Real rcv = w0_(m,IDN,k,j,i)/gm1;
    if (gen) {
      rcv = w0_(m,IDN,k,j,i)*eos_.SpecificHeatCv(w0_(m,IDN,k,j,i), w0_(m,IEN,k,j,i));
    }

    // THE FLUX-LIMITED DIFFUSIVITY.  With the limiter the face flux is
    //   F = kappa g / sqrt(1 + (kappa g / F_free)^2),   g = |grad T|,
    // and what an explicit step has to resolve is its linearisation
    //   dF/dg = kappa (1 + s^2)^(-3/2),   s = kappa g / F_free,
    // not kappa itself, which is unbounded where kappa_R rho -> 0: in the optically
    // thin atmosphere of a star kappa is ~1e8 x the limited value and the timestep
    // collapsed to microseconds.  Central differences of the cell temperature give g
    // per direction (the ghosts are filled when this runs).
    Real ffree = 0.0;
    if (radiative && limit) {
      const Real tk = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1)
                      *temp_unit;
      ffree = ffac*5.670374419e-5*tk*tk*tk*tk/(pres_unit*vel_unit);  // code units
    }
    auto tc = [&] (const int kk, const int jj, const int ii) {
      return (gen ? wtemp_(m,kk,jj,ii) : w0(m,IEN,kk,jj,ii)/w0(m,IDN,kk,jj,ii)*gm1);
    };
    auto keff = [&] (const Real dl, const Real tm, const Real tp) {
      if (!(radiative && limit) || !(ffree > 0.0)) return kappa_;
      const Real s = kappa_*fabs(tp - tm)/(2.0*dl)/ffree;
      return kappa_/((1.0 + s*s)*sqrt(1.0 + s*s));
    };

    // the cell's OWN dt candidate, kept separate from the running reduction minimum
    Real dtc = dt_huge;
    const Real d1 = (curv && radiative) ? dx1_(m,k,j,i) : size.d_view(m).dx1;
    if (!impx1) {
      const Real w1 = (taumode && blend_r) ? wmax : 1.0;
      const Real k1 = w1*keff(d1, tc(k,j,i-1), tc(k,j,i+1));
      if (k1 > 0.0) {
        min_dt = fmin(min_dt, SQR(d1)/k1*rcv);
        dtc = fmin(dtc, SQR(d1)/k1*rcv);
      }
      if (diag_) {
        cdg(m,3,k,j,i) = kappa_;
        cdg(m,4,k,j,i) = keff(d1, tc(k,j,i-1), tc(k,j,i+1));
      }
    }
    // on a curvilinear grid size.dx2/dx3 are ANGLES; the physical widths are pcoord's
    // cubed sphere: the exact operator's angular diffusivity is kappa/sin^2(alpha)
    const Real s2 = (cs && radiative) ? SQR(sinc_(m,k,j)) : 1.0;
    const Real wa = taumode ? wmax : 1.0;
    if (multi_d && !capa) {
      const Real d2 = (curv && radiative) ? dx2_(m,k,j,i) : size.d_view(m).dx2;
      const Real k2 = wa*keff(d2, tc(k,j-1,i), tc(k,j+1,i));
      if (k2 > 0.0) {
        min_dt = fmin(min_dt, SQR(d2)*s2/k2*rcv);
        dtc = fmin(dtc, SQR(d2)*s2/k2*rcv);
      }
    }
    if (three_d && !capa) {
      const Real d3 = (curv && radiative) ? dx3_(m,k,j,i) : size.d_view(m).dx3;
      const Real k3 = wa*keff(d3, tc(k-1,j,i), tc(k+1,j,i));
      if (k3 > 0.0) {
        min_dt = fmin(min_dt, SQR(d3)*s2/k3*rcv);
        dtc = fmin(dtc, SQR(d3)*s2/k3*rcv);
      }
    }
    if (diag_) cdg(m,5,k,j,i) = dtc;
    if (min_dt < mres.val) { mres.val = min_dt; mres.loc = idx; }
  }, Kokkos::MinLoc<Real, int>(mloc));
  dtnew = mloc.val*fac;
  // decode the winning cell, for the collapse report in Mesh::NewTimeStep
  if (mloc.loc >= 0 && mloc.loc < nmkji) {
    dtnew_m = (mloc.loc)/nkji;
    dtnew_k = (mloc.loc - dtnew_m*nkji)/nji + ks;
    dtnew_j = (mloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji)/nx1 + js;
    dtnew_i = (mloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji - (dtnew_j-js)*nx1) + is;
  } else {
    dtnew_m = dtnew_k = dtnew_j = dtnew_i = -1;
  }

  // THE STATE OF THAT CELL.  Where the cell is does not say why it is slow: the blend
  // weight, the flux limiter's s and the heat capacity all set the number and none of
  // them appears in any output.  Recompute them for the one winning cell, and only when
  // dtnew has just collapsed (or on the first call), so a healthy run pays nothing.
  dt_diag_valid = false;
  if (radiative && dtnew_m >= 0 && (dtnew_prev < 0.0 || dtnew < 0.25*dtnew_prev)) {
    auto dd = dt_diag;
    const int dm = dtnew_m, dk = dtnew_k, dj = dtnew_j, di = dtnew_i;
    auto &x1v_ = pmy_pack->pcoord->x1v;
    auto &tfd = rad_tauf;
    par_for("cond_dtdiag", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      const Real temp = (gen ? wtemp_(dm,dk,dj,di)
                             : w0_(dm,IEN,dk,dj,di)/w0_(dm,IDN,dk,dj,di)*gm1);
      const Real pres = (gen ? wder_(dm,IDPR,dk,dj,di) : w0_(dm,IEN,dk,dj,di)*gm1);
      const Real dens = w0_(dm,IDN,dk,dj,di);
      const Real tkap = KappaTemp(temp*temp_unit, tmax);
      const Real kr = (ktab ? RosselandTable(krt, krlT, krlP, krnT, krnP, tkap,
                                             krho ? dens*dens_unit : pres*pres_unit)
                            : -1.0);
      Real kappa_ = (ktab
          ? RadiativeKappaKR(tkap, dens*dens_unit, kfac, kr)
          : RadiativeKappa(tkap, pres*pres_unit, dens*dens_unit, met, kfac))
          /kappa_unit;
      if (gaterho > 0.0) kappa_ *= RadGate(dens*dens_unit, gaterho, gatedex);
      Real rcv = dens/gm1;
      if (gen) {
        rcv = dens*eos_.SpecificHeatCv(dens, w0_(dm,IEN,dk,dj,di));
      }
      Real ffree = 0.0;
      if (limit) {
        const Real tk = temp*temp_unit;
        ffree = ffac*5.670374419e-5*tk*tk*tk*tk/(pres_unit*vel_unit);
      }
      auto tc = [&] (const int kk, const int jj, const int ii) {
        return (gen ? wtemp_(dm,kk,jj,ii)
                    : w0_(dm,IEN,kk,jj,ii)/w0_(dm,IDN,kk,jj,ii)*gm1);
      };
      auto sof = [&] (const Real dl, const Real tm, const Real tp) {
        return (limit && ffree > 0.0) ? kappa_*fabs(tp - tm)/(2.0*dl)/ffree : 0.0;
      };
      auto keff = [&] (const Real s) {
        return (limit && ffree > 0.0) ? kappa_/((1.0 + s*s)*sqrt(1.0 + s*s)) : kappa_;
      };
      const Real wmax = taumode ? fmax(wf(dm,dk,dj,di), wf(dm,dk,dj,di+1)) : 1.0;
      const Real d1 = curv ? dx1_(dm,dk,dj,di) : size.d_view(dm).dx1;
      const Real d2 = curv ? dx2_(dm,dk,dj,di) : size.d_view(dm).dx2;
      const Real d3 = curv ? dx3_(dm,dk,dj,di) : size.d_view(dm).dx3;
      const Real s2 = cs ? SQR(sinc_(dm,dk,dj)) : 1.0;
      // GUARD the neighbours the mesh does not have: in 2-D ncells3 == 1 and in 1-D
      // ncells2 == 1, so tc(dk-1,..) / tc(..,dj-1,..) read off the end of w0.
      const Real s1v = sof(d1, tc(dk,dj,di-1), tc(dk,dj,di+1));
      const Real s2v = multi_d ? sof(d2, tc(dk,dj-1,di), tc(dk,dj+1,di)) : 0.0;
      const Real s3v = three_d ? sof(d3, tc(dk-1,dj,di), tc(dk+1,dj,di)) : 0.0;
      const Real w1 = (taumode && blend_r) ? wmax : 1.0;
      const Real wa = taumode ? wmax : 1.0;
      // x1v is only allocated on a curvilinear mesh; on Cartesian it is a 1x1 dummy
      dd.d_view(0)  = curv ? x1v_(dm,di) : -1.0;
      dd.d_view(1)  = dens*dens_unit;
      dd.d_view(2)  = temp*temp_unit;
      dd.d_view(3)  = pres*pres_unit;
      dd.d_view(4)  = kr;
      dd.d_view(5)  = kappa_;
      dd.d_view(6)  = rcv;
      dd.d_view(7)  = taumode ? wf(dm,dk,dj,di)   : 1.0;
      dd.d_view(8)  = taumode ? wf(dm,dk,dj,di+1) : 1.0;
      dd.d_view(9)  = taumode ? tfd(dm,dk,dj,di)  : -1.0;
      dd.d_view(10) = taumode ? tfd(dm,dk,dj,di+1): -1.0;
      dd.d_view(11) = s1v;
      dd.d_view(12) = (impx1 || !(w1*keff(s1v) > 0.0))
                      ? -1.0 : SQR(d1)/(w1*keff(s1v))*rcv*fac;
      dd.d_view(13) = (capa || !(multi_d && wa*keff(s2v) > 0.0))
                      ? -1.0 : SQR(d2)*s2/(wa*keff(s2v))*rcv*fac;
      dd.d_view(14) = (capa || !(three_d && wa*keff(s3v) > 0.0))
                      ? -1.0 : SQR(d3)*s2/(wa*keff(s3v))*rcv*fac;
      dd.d_view(15) = ffree;
    });
    dt_diag.modify_device();
    dt_diag.sync_host();
    dt_diag_valid = true;
  }
  dtnew_prev = dtnew;

  return;
}
