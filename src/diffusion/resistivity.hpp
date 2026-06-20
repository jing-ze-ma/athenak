#ifndef DIFFUSION_RESISTIVITY_HPP_
#define DIFFUSION_RESISTIVITY_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity.hpp
//  \brief Contains data and functions that implement various non-ideal MHD (resistive)
//  processes, such as Ohmic diffusion. TODO(@user): add ambipolar diffusion, Hall effect

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/meshblock.hpp"

//----------------------------------------------------------------------------------------
//! \class Resistivity
//  \brief data and functions that implement various resistive physics

class Resistivity {
 public:
  Resistivity(MeshBlockPack *pp, ParameterInput *pin);
  ~Resistivity();

  // data
  Real dtnew;
  std::string iso_resist_type;  // "constant" or "
  Real eta_ohm_const;
  DvceArray4D<Real> eta_b; // total resistivity of non-ideal MHD
  Real min_xe;
    
  DvceEdgeFld4D<Real> efld_resist;   // edge-centered electric fields due to non-ideal effects (E_{resistive} = \eta J)

  // functions to add resistive E-Field and energy flux
  void AddResistiveEMFs(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddResistiveFluxes(const DvceFaceFld4D<Real> &b0, const DvceArray5D<Real> &bc, DvceFaceFld5D<Real> &flx);
  void AddEMFConstantResist(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddEMFGeneralResist(const DvceFaceFld4D<Real> &b0, DvceEdgeFld4D<Real> &efld);
  void AddEMFDirect(const DvceEdgeFld4D<Real> &efld_resist, DvceEdgeFld4D<Real> &efld);
  void AddFluxConstantGridResist(const DvceFaceFld4D<Real> &b, DvceFaceFld5D<Real> &flx);
  void AddFluxGeneralResist(const DvceFaceFld4D<Real> &b, const DvceArray5D<Real> &bc, DvceFaceFld5D<Real> &flx);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void NewTimeStepConstantGridResist(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void NewTimeStepGeneralResist(const DvceArray5D<Real> &w, const EOS_Data &eos_data);
  void ClearResistiveEMFs(DvceEdgeFld4D<Real> &efld);
  void SetResistivity(const DvceArray5D<Real> &w, const Real &gamma, const Real &Rgas, DvceArray4D<Real> &eta_b, const int il, const int iu, const int jl, const int ju, const int kl, const int ku);

  KOKKOS_INLINE_FUNCTION
  void ResistivityPerna(const Real &nn, const Real &T, Real &eta) {
      // Perna+2010 <- Balbus & Hawley 2000
      Real ak = 1.0e-7;
      Real td25 = sqrt(sqrt(T/1.0e3));
      Real xe = 6.47e-13/(1.15e-11)*sqrt(ak/1.0e-7)*td25*td25*td25*sqrt(2.4e15/nn)*exp(-25188.0/T);
      xe = (xe < min_xe) ? min_xe : xe;
      eta = 230.0*sqrt(T)/xe;
  };

 private:
  MeshBlockPack* pmy_pack;
};

#endif // DIFFUSION_RESISTIVITY_HPP_
