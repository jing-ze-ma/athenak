//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_hyd_floors.cpp
//! \brief GeneralHydro::ConsToPrimFloors(), the inversion used when a floor switch of
//! EOS_Data (dfloor_keep_velocity, vceil, efloor_from_ekin, efloor_as_tfloor) is on.  It
//! is compiled apart from general_hyd.cpp so that the default kernel keeps its arithmetic
//! bitwise (EOS_Data::floors_legacy); changes meant for the default path go into both.

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
//! \fn void ConsToPrimFloors()
//! \brief Converts conserved into primitive variables with the floor switches active;
//! same contract as GeneralHydro::ConsToPrim()

void GeneralHydro::ConsToPrimFloors(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                                    const bool only_testfloors,
                                    const int il, const int iu, const int jl,
                                    const int ju, const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;
  auto &eos = eos_data;
  auto &fofc_ = pmy_pack->phydro->fofc;
  auto &wder_ = pmy_pack->phydro->wder;
  auto &wtemp_ = pmy_pack->phydro->wtemp;
  const bool keepv_defer_ = eos.dfloor_keep_velocity && eos.defer_cons_floors;
  auto dfl_fv_ = pmy_pack->phydro->dfl_fv;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0, nceilv_=0, ntclamp_=0, ntset_=0;
  Real efloor_de_=0.0;
  Real vceil_de_=0.0;
  Real dfloor_ke_=0.0;
  const bool kvtally = eos.dfloor_keep_velocity && !eos.defer_cons_floors;
  Kokkos::parallel_reduce("hyd_c2p_gen",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt, int &sumv,
                int &sumc, int &sumts, Real &sumde, Real &sumvde, Real &sumdke) {
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

    // call c2p function, which also returns the derived quantities (general_c2p_hyd.hpp)
    HydPrim1D w;
    Real pgas, g1;
    Real temp;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, mom_scaled=false;
    bool vceil_used=false, vceil_test=false, tclamp_used=false, tset_used=false;
    Real efloor_de=0.0;
    Real vceil_de=0.0;
    Real dfloor_fv=1.0;
    // floor-test pass (FOFC): a non-finite or non-positive trial state has no
    // temperature; flag the cell and skip the solve
    if (only_testfloors) {
      if (!Kokkos::isfinite(u.d) || !Kokkos::isfinite(u.mx) ||
          !Kokkos::isfinite(u.my) || !Kokkos::isfinite(u.mz) ||
          !Kokkos::isfinite(u.e) || (u.d <= 0.0) || (u.e <= 0.0)) {
        fofc_(m,k,j,i) = true;
        sumd++;
        return;
      }
    }
    // kinetic energy removed by dfloor_keep_velocity (output only)
    const Real ke_in = (kvtally && u.d > 0.0)
        ? 0.5*(SQR(u.mx) + SQR(u.my) + SQR(u.mz))/u.d : 0.0;
    // the cached temperature in this cell warm starts the T(d,e) root find
    SingleC2P_GeneralHyd(u, eos, w, wtemp_(m,k,j,i), temp, pgas, g1,
                         dfloor_used, efloor_used, tfloor_used, efloor_de, mom_scaled,
                         dfloor_fv, vceil_used, vceil_test, tclamp_used, tset_used,
                         vceil_de);
    // dfloor_keep_velocity on the cubed sphere: GnomonicEquiangleRaiseVel applies the
    // metric part of the correction with fv per cell (not on the floor-test pass, which
    // is followed by no RaiseVel)
    if (keepv_defer_ && !only_testfloors) dfl_fv_(m,k,j,i) = dfloor_fv;
    sumde += efloor_de;
    // the energy the velocity ceiling clipped: see EventCounters::vceil_de
    sumvde += vceil_de;
    if (kvtally && !only_testfloors && dfloor_fv < 1.0) {
      sumdke += (1.0 - dfloor_fv*dfloor_fv*dfloor_fv)*ke_in;
    }

    // set FOFC flag and quit loop if this function called only to check floors
    if (only_testfloors) {
      // the velocity ceiling counts as a floor event for FOFC (see ideal_hyd_floors.cpp)
      if (dfloor_used || efloor_used || tfloor_used || vceil_test) {
        fofc_(m,k,j,i) = true;
        sumd++;  // use dfloor as counter for when either is true
      }
    } else {
      // update counter, reset conserved if floor was hit
      if (dfloor_used) {
        // <hydro>/floor_tracer: the floor's mass joins the floor-mass scalar
        if (eos.ifltr >= 0) cons(m,eos.ifltr,k,j,i) += u.d - cons(m,IDN,k,j,i);
        cons(m,IDN,k,j,i) = u.d;
        sumd++;
      }
      if (efloor_used) {
        cons(m,IEN,k,j,i) = u.e;
        sume++;
      }
      // efloor_from_ekin (and dfloor_keep_velocity) rescaled the momentum: write it back
      if (mom_scaled) {
        cons(m,IM1,k,j,i) = u.mx;
        cons(m,IM2,k,j,i) = u.my;
        cons(m,IM3,k,j,i) = u.mz;
      }
      // dfloor_keep_velocity off the cubed sphere: write back the energy the c2p removed
      // with the kinetic energy
      if ((dfloor_fv < 1.0 || vceil_used) && !eos.defer_cons_floors) {
        cons(m,IEN,k,j,i) = u.e;
      }
      if (vceil_used) {
        sumv++;
      }
      if (tfloor_used) {
        cons(m,IEN,k,j,i) = u.e;
        sumt++;
      }
      // the inversion was clamped to the table's temperature range
      if (tclamp_used) {
        sumc++;
      }
      // efloor_as_tfloor rebuilt e from the temperature: the conserved energy follows
      // (off the cubed sphere; there GnomonicEquiangleRaiseVel writes it)
      if (tset_used) {
        if (!eos.defer_cons_floors) cons(m,IEN,k,j,i) = u.e;
        sumts++;
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
      // cache the temperature for other modules and as next stage's warm start
      wtemp_(m,k,j,i) = temp;
      // convert scalars (if any)
      for (int n=nhyd; n<(nhyd+nscal); ++n) {
        // apply scalar floor
        if (cons(m,n,k,j,i) < 0.0) {
          cons(m,n,k,j,i) = 0.0;
        }
        // the floor mass cannot exceed the cell mass
        if (n == eos.ifltr && cons(m,n,k,j,i) > u.d) {
          cons(m,n,k,j,i) = u.d;
        }
        prim(m,n,k,j,i) = cons(m,n,k,j,i)/u.d;
      }
    }
  }, Kokkos::Sum<int>(nfloord_), Kokkos::Sum<int>(nfloore_), Kokkos::Sum<int>(nfloort_),
     Kokkos::Sum<int>(nceilv_), Kokkos::Sum<int>(ntclamp_), Kokkos::Sum<int>(ntset_),
     Kokkos::Sum<Real>(efloor_de_), Kokkos::Sum<Real>(vceil_de_),
     Kokkos::Sum<Real>(dfloor_ke_));

  // store appropriate counters
  if (only_testfloors) {
    pmy_pack->pmesh->ecounter.nfofc += nfloord_;
  } else {
    pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
    pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
    pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;
    pmy_pack->pmesh->ecounter.neos_vceil  += nceilv_;
    pmy_pack->pmesh->ecounter.neos_tclamp += ntclamp_;
    pmy_pack->pmesh->ecounter.neos_tset += ntset_;
    pmy_pack->pmesh->ecounter.efloor_de += efloor_de_;
    pmy_pack->pmesh->ecounter.vceil_de += vceil_de_;
    pmy_pack->pmesh->ecounter.dfloor_ke += dfloor_ke_;
  }

  return;
}
