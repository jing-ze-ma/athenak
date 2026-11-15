//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file general_mhd_frozen.cpp
//! \brief GeneralMHD::ConsToPrimFrozen(), the conversion of the first ConToPrim after a
//! restart.  The table inversion of T(d, e) is not idempotent (a converged answer moves
//! by a few ulp when inverted again), while the straight run carries wtemp and wder into
//! the cycle it resumes.  This kernel therefore does only the algebraic half (density,
//! velocities, cell-centred field, internal energy) and reads, never writes, wtemp and
//! wder; floors are not re-applied (the restored state was floored before it was
//! written).  MHD::c2p_freeze_derived is cleared after the call. Compiled apart from
//! general_mhd.cpp so the default kernel keeps its arithmetic.

#include <float.h>

#include <iostream>
#include <string>

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "mhd/mhd.hpp"
#include "units/units.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void GeneralMHD::ConsToPrimFrozen()
//! \brief the algebraic half of ConsToPrim() with the thermodynamic cache frozen

void GeneralMHD::ConsToPrimFrozen(DvceArray5D<Real> &cons,
                                  const DvceFaceFld4D<Real> &b,
                                  DvceArray5D<Real> &prim, DvceArray5D<Real> &bcc,
                                  const int il, const int iu, const int jl,
                                  const int ju, const int kl, const int ku) {
  int &nmhd  = pmy_pack->pmhd->nmhd;
  int &nscal = pmy_pack->pmhd->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  // the same branch, with the same weights, as GeneralMHD::ConsToPrim: on the cubed
  // sphere and on spherical polar the cell centre is not the midpoint of its two faces
  const bool cs_str_x1_ = pmy_pack->pmesh->use_cubed_sphere;
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

  Kokkos::parallel_for("mhd_c2p_frozen",
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

    // cell-centred field from the face field, exactly as ConsToPrim forms it
    Real bx, by, bz;
    if (use_spherical_polar) {
      Real lw, rw;
      lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
      rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
      bx = lw*b.x1f(m,k,j,i) + rw*b.x1f(m,k,j,i+1);
      lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
      rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
      by = lw*b.x2f(m,k,j,i) + rw*b.x2f(m,k,j+1,i);
      lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
      rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
      bz = lw*b.x3f(m,k,j,i) + rw*b.x3f(m,k+1,j,i);
    } else {
      bx = cs_str_x1_ ?
          CellCenteredRadialFld(b.x1f(m,k,j,i), b.x1f(m,k,j,i+1),
                                x1f_(m,i), x1f_(m,i+1), x1v_(m,i)) :
          0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
      by = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
      bz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
    }

    // the same expressions, in the same order, as SingleC2P_GeneralMHD
    const Real e_k = 0.5*di*(SQR(mx) + SQR(my) + SQR(mz));
    const Real e_m = 0.5*(SQR(bx) + SQR(by) + SQR(bz));
    prim(m,IDN,k,j,i) = ud;
    prim(m,IVX,k,j,i) = di*mx;
    prim(m,IVY,k,j,i) = di*my;
    prim(m,IVZ,k,j,i) = di*mz;
    prim(m,IEN,k,j,i) = (cons(m,IEN,k,j,i) - e_k - e_m);
    bcc(m,IBX,k,j,i) = bx;
    bcc(m,IBY,k,j,i) = by;
    bcc(m,IBZ,k,j,i) = bz;
    for (int n=nmhd; n<(nmhd+nscal); ++n) {
      prim(m,n,k,j,i) = cons(m,n,k,j,i)/ud;
    }
  });

  return;
}
