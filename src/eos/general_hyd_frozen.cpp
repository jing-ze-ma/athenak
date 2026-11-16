//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_hyd_frozen.cpp
//! \brief GeneralHydro::ConsToPrimFrozen(), the conversion of the first ConToPrim after a
//! restart.  The table inversion of T(d, e) is not idempotent (a converged answer moves
//! by a few ulp when inverted again), while the straight run carries wtemp and wder into
//! the cycle it resumes.  This kernel therefore does only the algebraic half (density,
//! velocities, internal energy) and reads, never writes, wtemp and wder; floors are not
//! re-applied (the restored state was floored before it was written).
//! Hydro::c2p_freeze_derived is cleared after the call. Compiled apart from
//! general_hyd.cpp so the default kernel keeps its arithmetic.

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "hydro/hydro.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GeneralHydro::ConsToPrimFrozen()
//! \brief the algebraic half of ConsToPrim() with the thermodynamic cache frozen

void GeneralHydro::ConsToPrimFrozen(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                                    const int il, const int iu, const int jl,
                                    const int ju, const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  const int ni   = (iu - il + 1);
  const int nji  = (ju - jl + 1)*ni;
  const int nkji = (ku - kl + 1)*nji;
  const int nmkji = nmb*nkji;

  Kokkos::parallel_for("hyd_c2p_frozen",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx) {
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/ni;
    int i = (idx - m*nkji - k*nji - j*ni) + il;
    j += jl;
    k += kl;

    const Real ud = cons(m,IDN,k,j,i);
    const Real di = 1.0/ud;
    const Real mx = cons(m,IM1,k,j,i);
    const Real my = cons(m,IM2,k,j,i);
    const Real mz = cons(m,IM3,k,j,i);
    // the same expressions, in the same order, as SingleC2P_GeneralHydLegacy
    prim(m,IDN,k,j,i) = ud;
    prim(m,IVX,k,j,i) = di*mx;
    prim(m,IVY,k,j,i) = di*my;
    prim(m,IVZ,k,j,i) = di*mz;
    prim(m,IEN,k,j,i) = cons(m,IEN,k,j,i) - 0.5*di*(SQR(mx) + SQR(my) + SQR(mz));
    for (int n=nhyd; n<(nhyd+nscal); ++n) {
      prim(m,n,k,j,i) = cons(m,n,k,j,i)/ud;
    }
  });

  return;
}
