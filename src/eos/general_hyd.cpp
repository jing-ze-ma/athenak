//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_hyd.cpp
//! \brief derived class that implements a general EOS in nonrelativistic hydro
//!
//! Unlike an ideal gas, a general EOS is expensive to evaluate (the temperature inversion
//! is a root find) and is NOT scale free: ionization, dissociation and radiation pressure
//! all depend on absolute density and temperature. This class therefore
//!   (1) requires the run to declare its physical scale via the <units> block, and
//!   (2) evaluates pressure and Gamma_1 once per cell in ConsToPrim, storing them in the
//!       derived-variable array so that the Riemann solvers never call the EOS.

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "hydro/hydro.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_hyd.hpp"
#include "eos/ideal_c2p_hyd.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls EOS base class constructor

GeneralHydro::GeneralHydro(MeshBlockPack *pp, ParameterInput *pin) :
    EquationOfState("hydro", pp, pin) {
  eos_data.is_ideal = true;   // a general EOS still evolves an energy equation
  eos_data.eos_type = EOSType::general;
  eos_data.gamma = pin->GetReal("hydro","gamma");
  eos_data.iso_cs = 0.0;

  // A general EOS depends on ABSOLUTE density and temperature, so the run must declare
  // its physical scale. Without a <units> block the code units are arbitrary and the
  // ionization/dissociation/radiation terms would be silently meaningless.
  if (pp->punit == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<hydro>/eos = general requires a <units> block: a general EOS is not "
              << "scale free and needs the physical density and temperature scale."
              << std::endl << "Running in cgs (length_cgs = mass_cgs = time_cgs = 1) is "
              << "the recommended setup." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  eos_data.dens_cgs = pp->punit->density_cgs();
  eos_data.pres_cgs = pp->punit->pressure_cgs();
  // NOTE: Units::temperature_cgs() folds in the constant mean molecular weight <units>/mu
  // The general EOS carries composition itself (mu varies with the ionization state), so
  // the reference mu is divided back out here and the temperature scale is defined with
  // mu_ref = 1, i.e. T_cgs = v_code^2*m_u/k_B.
  eos_data.temp_cgs = pp->punit->temperature_cgs()/pp->punit->mu();

  // Build the interpolation table, unless the run asked for the gamma-law mode. The
  // gamma law is the default because it makes the general code path reproduce the ideal
  // one exactly, which is how the general interface is regression tested; a run that
  // wants real thermodynamics selects general_eos = table.
  BuildGeneralEOS("hydro", pin);
}

//----------------------------------------------------------------------------------------
//! \fn void ConsToPrim()
//! \brief Converts conserved into primitive variables, and evaluates the derived
//! thermodynamic variables (pressure, Gamma_1) that are reconstructed to interfaces.
//! Operates over range of cells given in argument list. Number of times floors used
//! stored into event counters.

void GeneralHydro::ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                              const bool only_testfloors,
                              const int il, const int iu, const int jl, const int ju,
                              const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;
  auto &eos = eos_data;
  auto &fofc_ = pmy_pack->phydro->fofc;
  auto &wder_ = pmy_pack->phydro->wder;
  auto &wtemp_ = pmy_pack->phydro->wtemp;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0;
  Kokkos::parallel_reduce("hyd_c2p_gen",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    j += jl;
    k += kl;

    // load single state conserved variables
    HydCons1D u;
    u.d  = cons(m,IDN,k,j,i);
    u.mx = cons(m,IM1,k,j,i);
    u.my = cons(m,IM2,k,j,i);
    u.mz = cons(m,IM3,k,j,i);
    u.e  = cons(m,IEN,k,j,i);

    // call c2p function, which also returns the derived thermodynamic quantities
    // (inline function in general_c2p_hyd.hpp file)
    HydPrim1D w;
    Real pgas, g1;
    Real temp;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false;
    // the cached temperature in this cell warm starts the T(d,e) root find
    SingleC2P_GeneralHyd(u, eos, w, wtemp_(m,k,j,i), temp, pgas, g1,
                         dfloor_used, efloor_used, tfloor_used);

    // set FOFC flag and quit loop if this function called only to check floors
    if (only_testfloors) {
      if (dfloor_used || efloor_used || tfloor_used) {
        fofc_(m,k,j,i) = true;
        sumd++;  // use dfloor as counter for when either is true
      }
    } else {
      // update counter, reset conserved if floor was hit
      if (dfloor_used) {
        cons(m,IDN,k,j,i) = u.d;
        sumd++;
      }
      if (efloor_used) {
        cons(m,IEN,k,j,i) = u.e;
        sume++;
      }
      if (tfloor_used) {
        cons(m,IEN,k,j,i) = u.e;
        sumt++;
      }
      // store primitive state in 3D array
      prim(m,IDN,k,j,i) = w.d;
      prim(m,IVX,k,j,i) = w.vx;
      prim(m,IVY,k,j,i) = w.vy;
      prim(m,IVZ,k,j,i) = w.vz;
      prim(m,IEN,k,j,i) = w.e;
      // store derived thermodynamic variables for reconstruction
      wder_(m,IDPR,k,j,i) = pgas;
      wder_(m,IDG1,k,j,i) = g1;
      // cache the temperature solved for above; it is both the value other modules read
      // and the warm start for this cell's inversion at the next stage
      wtemp_(m,k,j,i) = temp;
      // convert scalars (if any)
      for (int n=nhyd; n<(nhyd+nscal); ++n) {
        // apply scalar floor
        if (cons(m,n,k,j,i) < 0.0) {
          cons(m,n,k,j,i) = 0.0;
        }
        prim(m,n,k,j,i) = cons(m,n,k,j,i)/u.d;
      }
    }
  }, Kokkos::Sum<int>(nfloord_), Kokkos::Sum<int>(nfloore_), Kokkos::Sum<int>(nfloort_));

  // store appropriate counters
  if (only_testfloors) {
    pmy_pack->pmesh->ecounter.nfofc += nfloord_;
  } else {
    pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
    pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
    pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void PrimToCons()
//! \brief Converts primitive into conserved variables. Operates over range of cells given
//! in argument list.  Floors never needed.
//! NOTE: for non-relativistic hydrodynamics this transformation is EOS independent
//! (E = e + E_kin involves no thermodynamics), so the ideal-gas kernel is reused.

void GeneralHydro::PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                              const int il, const int iu, const int jl, const int ju,
                              const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  par_for("hyd_p2c_gen", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // load single state primitive variables
    HydPrim1D w;
    w.d  = prim(m,IDN,k,j,i);
    w.vx = prim(m,IVX,k,j,i);
    w.vy = prim(m,IVY,k,j,i);
    w.vz = prim(m,IVZ,k,j,i);
    w.e  = prim(m,IEN,k,j,i);

    // call p2c function
    HydCons1D u;
    SingleP2C_IdealHyd(w, u);

    // store conserved state in 3D array
    cons(m,IDN,k,j,i) = u.d;
    cons(m,IM1,k,j,i) = u.mx;
    cons(m,IM2,k,j,i) = u.my;
    cons(m,IM3,k,j,i) = u.mz;
    cons(m,IEN,k,j,i) = u.e;

    // convert scalars (if any)
    for (int n=nhyd; n<(nhyd+nscal); ++n) {
      cons(m,n,k,j,i) = u.d*prim(m,n,k,j,i);
    }
  });

  return;
}
