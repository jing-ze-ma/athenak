#ifndef DIFFUSION_CONDUCTION_HPP_
#define DIFFUSION_CONDUCTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.hpp
//! \brief Contains data and functions that implement various formulations for conduction.
//  Currently only isotropic conduction implemented

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"

//----------------------------------------------------------------------------------------
//! \class Conduction
//! \brief data and functions that implement thermal conduction in Hydro and MHD

class Conduction {
 public:
  Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~Conduction();

  // data
  Real dtnew;
  std::string iso_cond_type; // "constant", "spitzer", "spitzer_limited", "radiative"
  // radiative conductivity kappa_rad = 16 sigma T^3/(3 kappa_R rho), kappa_R the
  // Freedman+2014 Rosseland mean (utils/rosseland.hpp), applied on faces whose pressure
  // is at or above rad_pcut (the layers the two-stream RT does not reach), flux-limited
  // to
  // sigma T^4, with an imposed heat flux through the inner x1 wall (the planet's internal
  // flux). rad_kappa_fac scales kappa_R (tests).
  Real rad_met = 0.0;          // [M/H] for the Rosseland fit
  Real rad_pcut = 0.0;         // code units; faces with p < rad_pcut get no flux
  Real rad_flux_inner = 0.0;   // code units; heat flux through the inner x1 wall
  Real rad_kappa_fac = 1.0;
  bool rad_flux_limit = true;
  Real kappa_iso;            // isotropic thermal conductivity
  Real kappa_iso_limit;      // limit to isotropic thermal conductivity

  // functions
  void AddHeatFluxes(const DvceArray5D<Real> &w, const EOS_Data &eos,
                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxRadiative(const DvceArray5D<Real> &w, const EOS_Data &eos,
    DvceFaceFld5D<Real> &flx);
  void AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                       DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);

 private:
  MeshBlockPack* pmy_pack;
  // "hydro" or "mhd": identifies which physics module owns this Conduction object, so
  // that the general EOS derived-variable arrays (cached temperature and pressure) of the
  // right module can be found.
  std::string my_block;
};
#endif // DIFFUSION_CONDUCTION_HPP_
