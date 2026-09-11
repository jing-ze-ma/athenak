//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ideal_hyd.cpp
//! \brief derived class that implements ideal gas EOS in nonrelativistic hydro

#include "athena.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"
#include "eos/ideal_c2p_hyd.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls EOS base class constructor

IdealHydro::IdealHydro(MeshBlockPack *pp, ParameterInput *pin) :
    EquationOfState("hydro", pp, pin) {
  eos_data.is_ideal = true;
  eos_data.gamma = pin->GetReal("hydro","gamma");
  eos_data.iso_cs = 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn void ConsToPrim()
//! \brief Converts conserved into primitive variables. Operates over range of cells given
//! in argument list. Number of times floors used stored into event counters.

void IdealHydro::ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                            const bool only_testfloors,
                            const int il, const int iu, const int jl, const int ju,
                            const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;
  auto &eos = eos_data;
  const bool keepv_defer_ = eos.dfloor_keep_velocity && eos.defer_cons_floors;
  auto dfl_fv_ = pmy_pack->phydro->dfl_fv;
  auto &fofc_ = pmy_pack->phydro->fofc;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  int nfloord_=0, nfloore_=0, nfloort_=0, nceilv_=0;
  Kokkos::parallel_reduce("hyd_c2p",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, int &sumd, int &sume, int &sumt, int &sumv) {
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

    // call c2p function
    // (inline function in ideal_c2p_hyd.hpp file)
    HydPrim1D w;
    bool dfloor_used=false, efloor_used=false, tfloor_used=false, vceil_used=false;
    Real dfloor_fv=1.0;
    SingleC2P_IdealHyd(u, eos, w, dfloor_used, efloor_used, tfloor_used, dfloor_fv,
                       vceil_used);
    // The floor-TEST pass (FOFC) must leave no trace: it is handed scratch conserved
    // data and is followed by no GnomonicEquiangleRaiseVel, so neither the momentum
    // rescale nor the fv it would hand on may be written from here.
    if (!only_testfloors) {
      if ((eos.dfloor_keep_velocity && dfloor_fv < 1.0) || vceil_used) {
        cons(m,IM1,k,j,i) = u.mx;
        cons(m,IM2,k,j,i) = u.my;
        cons(m,IM3,k,j,i) = u.mz;
        if (!eos.defer_cons_floors) cons(m,IEN,k,j,i) = u.e;
      }
      if (keepv_defer_) dfl_fv_(m,k,j,i) = dfloor_fv;
    }

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
      if (vceil_used) {
        sumv++;
      }
      // store primitive state in 3D array
      prim(m,IDN,k,j,i) = w.d;
      prim(m,IVX,k,j,i) = w.vx;
      prim(m,IVY,k,j,i) = w.vy;
      prim(m,IVZ,k,j,i) = w.vz;
      prim(m,IEN,k,j,i) = w.e;
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
     Kokkos::Sum<int>(nceilv_));

  // store appropriate counters
  if (only_testfloors) {
    pmy_pack->pmesh->ecounter.nfofc += nfloord_;
  } else {
    pmy_pack->pmesh->ecounter.neos_dfloor += nfloord_;
    pmy_pack->pmesh->ecounter.neos_efloor += nfloore_;
    pmy_pack->pmesh->ecounter.neos_tfloor += nfloort_;
    pmy_pack->pmesh->ecounter.neos_vceil  += nceilv_;
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void PrimToCons()
//! \brief Converts primitive into conserved variables. Operates over range of cells given
//! in argument list.  Floors never needed.

void IdealHydro::PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                            const int il, const int iu, const int jl, const int ju,
                            const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  par_for("hyd_p2c", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
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
