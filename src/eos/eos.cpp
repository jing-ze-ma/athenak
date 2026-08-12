//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos.cpp
//! \brief implements constructor and some fns for EquationOfState abstract base class

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
// EquationOfState constructor

EquationOfState::EquationOfState(std::string bk, MeshBlockPack* pp, ParameterInput *pin) :
    pmy_pack(pp) {
  eos_data.dfloor = pin->GetOrAddReal(bk,"dfloor",(FLT_MIN));
  eos_data.pfloor = pin->GetOrAddReal(bk,"pfloor",(FLT_MIN));
  eos_data.tfloor = pin->GetOrAddReal(bk,"tfloor",(FLT_MIN));
  eos_data.sfloor = pin->GetOrAddReal(bk,"sfloor",(FLT_MIN));
}

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::BuildGeneralEOS()
//! \brief select what the general EOS interface evaluates, and build the table if asked.
//!
//! Two modes, chosen by <block>/general_eos:
//!   gamma  (default) -- the interface evaluates the same ideal gamma law as the ideal
//!                       EOS. This is not a physics option so much as a control: it makes
//!                       the general code path -- reconstruction of p and Gamma_1, the
//!                       general Riemann solvers, the general well-balanced background --
//!                       exercisable against an answer that is known exactly, which is
//!                       what the general-vs-ideal regression tests assert.
//!   table            -- the tabulated analytic EOS: H2 dissociation, H and He
//!                       ionization, and optionally radiation pressure.
//!
//! Called only from the GeneralHydro/GeneralMHD constructors, and only after dens_cgs,
//! pres_cgs and temp_cgs have been set, because the table is built in physical units.

void EquationOfState::BuildGeneralEOS(std::string block, ParameterInput *pin) {
  std::string mode = pin->GetOrAddString(block, "general_eos", "gamma");
  if (mode.compare("gamma") == 0) {
    eos_data.tbl.active = false;
  } else if (mode.compare("table") == 0) {
    BuildEOSTable(eos_data.tbl, pin, block, eos_data.dens_cgs, eos_data.pres_cgs,
                  eos_data.temp_cgs, eos_data.pfloor);
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/general_eos = '" << mode << "' not recognized; "
              << "must be 'gamma' or 'table'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ConsToPrim()
//! \brief No-Op versions of hydro and MHD conservative to primitive functions.
//! Required because each derived class overrides only one.

void EquationOfState::ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                                 const bool only_testfloors,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}

void EquationOfState::ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                                 DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                                 const bool only_testfloors,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}

//----------------------------------------------------------------------------------------
//! \fn void PrimToCon()
//! \brief No-Op versions of hydro and MHD primitive to conservative functions.
//! Required because each derived class overrides only one.

void EquationOfState::PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}
void EquationOfState::PrimToCons(const DvceArray5D<Real> &prim,
                                 const DvceArray5D<Real> &bcc, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
}
