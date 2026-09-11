//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_hyd_floors.cpp
//! \brief GeneralHydro::ConsToPrimFloors(), the conserved-to-primitive inversion used
//! when any of the floor switches of EOS_Data (dfloor_keep_velocity, vceil,
//! eos_floor_consistent, efloor_from_ekin) is enabled.
//!
//! WHY IT IS A FILE OF ITS OWN.  The default path is GeneralHydro::ConsToPrim() in
//! general_hyd.cpp, which is kept exactly as it was before these switches existed so
//! that a build carrying them reproduces an older run bit for bit.  That only works
//! while the two kernels are compiled APART: with both in one translation unit the
//! device compiler contracts the arithmetic of the untouched one differently, and the
//! default answer moves by a unit in the last place -- which on the cubed sphere grows
//! into a visibly different run within a few thousand cycles.  Keep the two files
//! separate, and mirror into general_hyd.cpp any change that is meant to apply by
//! DEFAULT.

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
//! \brief Converts conserved into primitive variables with the floor switches active.
//! Same contract as GeneralHydro::ConsToPrim(); see the note above.

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

  int nfloord_=0, nfloore_=0, nfloort_=0, nceilv_=0;
  Real efloor_de_=0.0;
  Kokkos::parallel_reduce("hyd_c2p_gen",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt, int &sumv,
                Real &sumde) {
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
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, mom_scaled=false;
    bool vceil_used=false;
    Real efloor_de=0.0;
    Real dfloor_fv=1.0;
    // the cached temperature in this cell warm starts the T(d,e) root find
    SingleC2P_GeneralHyd(u, eos, w, wtemp_(m,k,j,i), temp, pgas, g1,
                         dfloor_used, efloor_used, tfloor_used, efloor_de, mom_scaled,
                         dfloor_fv, vceil_used);
    // <hydro>/dfloor_keep_velocity on the cubed sphere: the metric-correct part of the
    // correction is applied by GnomonicEquiangleRaiseVel, which needs fv per cell.  Not
    // on the floor-TEST pass (FOFC): that one is handed scratch conserved data and is
    // followed by no RaiseVel, so the fv it computed must not reach one.
    if (keepv_defer_ && !only_testfloors) dfl_fv_(m,k,j,i) = dfloor_fv;
    sumde += efloor_de;

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
      // <hydro>/efloor_from_ekin (and the density floor under dfloor_keep_velocity)
      // paid by rescaling the momentum, so the conserved momentum has to follow the
      // primitive velocity down
      if (mom_scaled) {
        cons(m,IM1,k,j,i) = u.mx;
        cons(m,IM2,k,j,i) = u.my;
        cons(m,IM3,k,j,i) = u.mz;
      }
      // dfloor_keep_velocity OFF the cubed sphere takes the removed kinetic energy out
      // of u.e inside the c2p; without this write the energy write above happens only
      // when a pressure or temperature floor also fired, and the cell would keep the
      // energy it no longer has any kinetic energy for.
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
  }, Kokkos::Sum<int>(nfloord_), Kokkos::Sum<int>(nfloore_), Kokkos::Sum<int>(nfloort_),
     Kokkos::Sum<int>(nceilv_), Kokkos::Sum<Real>(efloor_de_));

  // store appropriate counters
  if (only_testfloors) {
    pmy_pack->pmesh->ecounter.nfofc += nfloord_;
  } else {
    pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
    pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
    pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;
    pmy_pack->pmesh->ecounter.neos_vceil  += nceilv_;
    pmy_pack->pmesh->ecounter.efloor_de += efloor_de_;
  }

  return;
}
