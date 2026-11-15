//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_rkgdt.cpp
//! \brief the diffusive dt of the RKG super-step from the forward-Euler limit
//! 2/lambda_max of the resistive operator (<mhd>/rkg_dt_safety), instead of the
//! isotropic cfl_number x (1/6) min_d dx_d^2/eta, which is small in thin radial cells.
//! The bound, per cell and maximised over the mesh:
//!   lambda_c = 4 eta_r (1/dx1)(1/dc_- + 1/dc_+)/2 + 4 eta_h (1/dx2^2 + 1/dx3^2)/s^2
//! dc_-/+ are the centroid spacings to the radial neighbours; eta_r is the maximum eta
//! over the cells of the edges carrying the radial derivative, eta_h the maximum over the
//! 3x3x3 block; s is the sine of the cubed-sphere cell angle, and on spherical polar dx3
//! is taken at the theta face nearer the pole.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "bvals/bvals.hpp"
#include "pgen/pgen.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/resistivity.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn Real Resistivity::LambdaBoundMax
//! \brief max over the mesh (all ranks) of the per-cell bound lambda_c above

Real Resistivity::LambdaBoundMax() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool cs = pmy_pack->pmesh->use_cubed_sphere;
  auto dx1_ = pmy_pack->pcoord->dx1;
  auto dx2_ = pmy_pack->pcoord->dx2;
  auto dx3_ = pmy_pack->pcoord->dx3;
  auto x1v_ = pmy_pack->pcoord->x1v;
  auto scell = pmy_pack->pcoord->sin_cell;
  auto x2f_ = pmy_pack->pcoord->xx2f;
  auto x3f_ = pmy_pack->pcoord->xx3f;
  const bool sp = pmy_pack->pmesh->use_spherical_polar;
  auto eta_ = eta_b;
  // eta_O + eta_A per cell (eta_A = min(ad_coef |B|^2, max_eta_ad)): the ambipolar
  // operator is the Ohmic one with J projected perpendicular to B, so its spectrum is
  // bounded by the isotropic one with eta_A added
  const bool ad_on = (ad_type != 0);
  auto coef_ = ad_coef;
  auto bcc_ = pmy_pack->pmhd->bcc0;
  const Real cap_ad = max_eta_ad;

  Real lmax = 0.0;
  Kokkos::parallel_reduce("res_rkgdt_bound",
                          Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, Real &lm) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    const int k0 = three_d ? k-1 : k, k1 = three_d ? k+1 : k;
    const int j0 = multi_d ? j-1 : j, j1 = multi_d ? j+1 : j;
    // eta_r: the cells of the x2/x3 edges that carry the radial derivative (one
    // transverse offset at a time); eta_h: the whole 3x3x3 block
    Real er = 0.0, eh = 0.0;
    for (int kk = k0; kk <= k1; ++kk) {
      for (int jj = j0; jj <= j1; ++jj) {
        const bool diag = (kk != k) && (jj != j);
        for (int ii = i-1; ii <= i+1; ++ii) {
          Real e = eta_(m,kk,jj,ii);
          if (ad_on) {
            const Real b2 = SQR(bcc_(m,IBX,kk,jj,ii)) + SQR(bcc_(m,IBY,kk,jj,ii))
                          + SQR(bcc_(m,IBZ,kk,jj,ii));
            e += fmin(coef_(m,kk,jj,ii)*b2, cap_ad);
          }
          eh = fmax(eh, e);
          if (!diag) er = fmax(er, e);
        }
      }
    }
    const Real dcm = x1v_(m,i) - x1v_(m,i-1);
    const Real dcp = x1v_(m,i+1) - x1v_(m,i);
    Real lam = 2.0*(1.0/dcm + 1.0/dcp)/dx1_(m,k,j,i);
    Real h = 0.0;
    if (multi_d) h += 4.0/SQR(dx2_(m,k,j,i));
    if (three_d) {
      Real d3 = dx3_(m,k,j,i);
      if (sp) {
        // the phi length at the theta face nearer the pole; a cell on the pole takes
        // sin(dtheta/2)
        const Real s0 = fabs(sin(x2f_(m,j))), s1 = fabs(sin(x2f_(m,j+1)));
        const Real sh = fabs(sin(0.5*(x2f_(m,j+1) - x2f_(m,j))));
        const Real se = fmax(fmin(s0, s1), sh);
        d3 = x1v_(m,i)*(x3f_(m,k+1) - x3f_(m,k))*se;
      }
      h += 4.0/SQR(d3);
    }
    if (cs) h /= SQR(scell(m,k,j));
    lam = er*lam + eh*h;
    lm = fmax(lm, lam);
  }, Kokkos::Max<Real>(lmax));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &lmax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  return lmax;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::NewTimeStepRKGFE
//! \brief dtnew = rkg_dt_safety x 2/lambda_bound (mesh.cpp takes it as dt_diff without
//! the CFL factor)

void Resistivity::NewTimeStepRKGFE() {
  const Real lam = LambdaBoundMax();
  dtnew = (lam > 0.0) ? rkg_dt_safety*2.0/lam
                      : static_cast<Real>(std::numeric_limits<float>::max());
}
