//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_hyd_frozen.cpp
//! \brief GeneralHydro::ConsToPrimFrozen(), the conserved-to-primitive conversion used
//! ONCE, by the first ConToPrim of a RESTARTED run.
//!
//! WHY IT EXISTS.  A restart is meant to be a bitwise continuation, and for a general
//! EOS it was not.  The conversion has two halves: the algebraic one (density,
//! velocities, internal energy), which is a pure function of the conserved variables the
//! restart file restores exactly, and the THERMODYNAMIC one -- solve e(d,T) for T, then
//! take p and Gamma_1 at that T -- which is an iterative table inversion warm started
//! from the cached temperature.  That inversion is NOT idempotent: handed back its own
//! converged answer it takes another step and lands a few ULP away (measured: 24% of
//! cells move, by up to 5e-14 relative).  The straight run does not do that conversion
//! at the start of the cycle it is resuming -- it carries wtemp and wder from the last
//! stage of the previous cycle -- while Driver::Initialize does, so the two runs enter
//! the first cycle with different temperatures and different reconstructed pressures.
//! On the He box that shows up in the first post-restart history row at 1e-11 in the
//! transverse momentum, and it grows.
//!
//! WHAT IT DOES.  Exactly the algebraic half, and nothing else: wtemp and wder are read,
//! not written, so the state the restart file restored survives into the first Fluxes
//! call untouched.  Floors are not re-applied: the conserved variables in the file were
//! written AFTER the conversion that floored them, so every floor test would be false
//! anyway, and re-running them could only move the answer.  The next ConToPrim, at the
//! end of the first stage, is the ordinary one again (Hydro::c2p_freeze_derived is
//! cleared after this call).
//!
//! WHY IT IS A FILE OF ITS OWN.  See the note at the top of general_hyd_floors.cpp: with
//! two kernels in one translation unit the device compiler contracts the arithmetic of
//! the untouched one differently, and the DEFAULT answer moves by a unit in the last
//! place.  general_hyd.cpp is not touched by this feature at all.

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "hydro/hydro.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GeneralHydro::ConsToPrimFrozen()
//! \brief the algebraic half of ConsToPrim(), with the thermodynamic cache frozen.

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
