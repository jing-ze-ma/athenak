//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_mhd.cpp
//! \brief derived class that implements a general EOS in nonrelativistic MHD
//!
//! See general_hyd.cpp for the design rationale: a general EOS is expensive and is not
//! scale free, so it requires a <units> block and evaluates pressure and Gamma_1 once per
//! cell in ConsToPrim, storing them for reconstruction to interfaces.

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "mhd/mhd.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"
#include "eos/general_c2p_mhd.hpp"
#include "eos/ideal_c2p_mhd.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls EOS base class constructor

GeneralMHD::GeneralMHD(MeshBlockPack *pp, ParameterInput *pin) :
    EquationOfState("mhd", pp, pin) {
  eos_data.is_ideal = true;   // a general EOS still evolves an energy equation
  eos_data.eos_type = EOSType::general;
  eos_data.gamma = pin->GetReal("mhd","gamma");
  eos_data.iso_cs = 0.0;

  // A general EOS depends on ABSOLUTE density and temperature, so the run must declare
  // its physical scale (see general_hyd.cpp).
  if (pp->punit == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<mhd>/eos = general requires a <units> block: a general EOS is not "
              << "scale free and needs the physical density and temperature scale."
              << std::endl << "Running in cgs (length_cgs = mass_cgs = time_cgs = 1) is "
              << "the recommended setup." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  eos_data.dens_cgs = pp->punit->density_cgs();
  eos_data.pres_cgs = pp->punit->pressure_cgs();
  // reference mu divided back out; see general_hyd.cpp
  eos_data.temp_cgs = pp->punit->temperature_cgs()/pp->punit->mu();

  // gamma-law or tabulated; see EquationOfState::BuildGeneralEOS()
  BuildGeneralEOS("mhd", pin);
}

//----------------------------------------------------------------------------------------
//! \!fn void ConsToPrim()
//! \brief Converts conserved into primitive variables, and evaluates the derived
//! thermodynamic variables (pressure, Gamma_1) that are reconstructed to interfaces.
//! Operates over range of cells given in argument list.

void GeneralMHD::ConsToPrim(DvceArray5D<Real> &cons, const DvceFaceFld4D<Real> &b,
                            DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                            const bool only_testfloors,
                            const int il, const int iu, const int jl, const int ju,
                            const int kl, const int ku) {
  int &nmhd  = pmy_pack->pmhd->nmhd;
  int &nscal = pmy_pack->pmhd->nscalars;
  int &nmb = pmy_pack->nmb_thispack;
  auto &eos = eos_data;
  auto &fofc_ = pmy_pack->pmhd->fofc;
  auto &wder_ = pmy_pack->pmhd->wder;
  auto &wtemp_ = pmy_pack->pmhd->wtemp;

  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;
  auto &x2v_ = pmy_pack->pcoord->x2v;
  auto &x2f_ = pmy_pack->pcoord->xx2f;
  auto &x3v_ = pmy_pack->pcoord->x3v;
  auto &x3f_ = pmy_pack->pcoord->xx3f;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0;
  Kokkos::parallel_reduce("mhd_c2p_gen",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    j += jl;
    k += kl;

    // load single state conserved variables
    MHDCons1D u;
    u.d  = cons(m,IDN,k,j,i);
    u.mx = cons(m,IM1,k,j,i);
    u.my = cons(m,IM2,k,j,i);
    u.mz = cons(m,IM3,k,j,i);
    u.e  = cons(m,IEN,k,j,i);

    // load cell-centered fields into conserved state
    // use input CC fields if only testing floors with FOFC
    if (only_testfloors) {
      u.bx = bcc(m,IBX,k,j,i);
      u.by = bcc(m,IBY,k,j,i);
      u.bz = bcc(m,IBZ,k,j,i);
    // else use simple linear average of face-centered fields
    } else {
      if (use_spherical_polar) {
        Real lw, rw;
        lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
        rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
        u.bx = lw*b.x1f(m,k,j,i) + rw*b.x1f(m,k,j,i+1);
        lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
        rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
        u.by = lw*b.x2f(m,k,j,i) + rw*b.x2f(m,k,j+1,i);
        lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
        rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
        u.bz = lw*b.x3f(m,k,j,i) + rw*b.x3f(m,k+1,j,i);
      } else {
        u.bx = 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
        u.by = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
        u.bz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
      }
    }

    // call c2p function, which also returns the derived thermodynamic quantities
    // (inline function in general_c2p_mhd.hpp file)
    HydPrim1D w;
    Real pgas, g1;
    Real temp;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false;
    // the cached temperature in this cell warm starts the T(d,e) root find
    SingleC2P_GeneralMHD(u, eos, w, wtemp_(m,k,j,i), temp, pgas, g1,
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
      // store cell-centered fields in 3D array
      bcc(m,IBX,k,j,i) = u.bx;
      bcc(m,IBY,k,j,i) = u.by;
      bcc(m,IBZ,k,j,i) = u.bz;
      // store derived thermodynamic variables for reconstruction
      wder_(m,IDPR,k,j,i) = pgas;
      wder_(m,IDG1,k,j,i) = g1;
      // cache the temperature solved for above; it is both the value other modules read
      // and the warm start for this cell's inversion at the next stage
      wtemp_(m,k,j,i) = temp;
      // convert scalars (if any), always stored at end of cons and prim arrays.
      for (int n=nmhd; n<(nmhd+nscal); ++n) {
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
//! \!fn void PrimToCons()
//! \brief Converts primitive into conserved variables. Operates over range of cells given
//! in argument list. Does not change cell- or face-centered magnetic fields.
//! NOTE: for non-relativistic MHD this transformation is EOS independent, so the
//! ideal-gas kernel is reused.

void GeneralMHD::PrimToCons(const DvceArray5D<Real> &prim, const DvceArray5D<Real> &bcc,
                            DvceArray5D<Real> &cons, const int il, const int iu,
                            const int jl, const int ju, const int kl, const int ku) {
  int &nmhd  = pmy_pack->pmhd->nmhd;
  int &nscal = pmy_pack->pmhd->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  par_for("mhd_p2c_gen", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    // load single state primitive variables
    MHDPrim1D w;
    w.d  = prim(m,IDN,k,j,i);
    w.vx = prim(m,IVX,k,j,i);
    w.vy = prim(m,IVY,k,j,i);
    w.vz = prim(m,IVZ,k,j,i);
    w.e  = prim(m,IEN,k,j,i);

    // load cell-centered fields into primitive state
    w.bx = bcc(m,IBX,k,j,i);
    w.by = bcc(m,IBY,k,j,i);
    w.bz = bcc(m,IBZ,k,j,i);

    // call p2c function
    HydCons1D u;
    SingleP2C_IdealMHD(w, u);

    // store conserved state in 3D array
    cons(m,IDN,k,j,i) = u.d;
    cons(m,IM1,k,j,i) = u.mx;
    cons(m,IM2,k,j,i) = u.my;
    cons(m,IM3,k,j,i) = u.mz;
    cons(m,IEN,k,j,i) = u.e;

    // convert scalars (if any), always stored at end of cons and prim arrays.
    for (int n=nmhd; n<(nmhd+nscal); ++n) {
      cons(m,n,k,j,i) = u.d*prim(m,n,k,j,i);
    }
  });

  return;
}
