//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity_rkgdt.cpp
//! \brief the diffusive dt of the RKG super-step from the forward-Euler limit 2/lambda_max
//! of the resistive operator (<mhd>/rkg_dt_safety), instead of the isotropic
//! cfl_number x (1/6) min_d dx_d^2/eta, which is ~10x small in thin radial cells.
//!
//! The bound, per cell and maximised over the mesh:
//!   lambda_c = 4 eta_r (1/dx1)(1/dc_- + 1/dc_+)/2 + 4 eta_h (1/dx2^2 + 1/dx3^2)/s^2
//! dc_-/+ are the centroid spacings to the radial neighbours; eta_r is the maximum edge
//! eta over the cells carrying the radial derivative (3 radial x 5 transverse), eta_h the
//! maximum over the 3x3x3 block (the diagonal ghosts, unfilled at a cube vertex, enter
//! eta_h only); s is the sine of the cs cell angle, and on sp dx3 is taken at the theta
//! face nearer the pole.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

#include "athena.hpp"
#include "globals.hpp"
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
//! \brief max over the mesh (all ranks) of the per-cell bound lambda_c above; optionally
//! the radius and eta of the cell where it is attained.

Real Resistivity::LambdaBoundMax(Real *rloc, Real *etaloc) {
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
  // ambipolar-1101 merge: eta_O + eta_A per cell (eta_A = min(ad_coef |B|^2, max_eta_ad),
  // the coefficient AddEMFAmbipolar uses).  The AD operator eta_A (J - (J.B)B/B^2) is the
  // Ohmic one with J projected perpendicular to B, so its frozen-coefficient spectrum is
  // bounded by the isotropic one with eta_A added.
  const bool ad_on = (ad_type != 0);
  auto coef_ = ad_coef;
  auto bcc_ = pmy_pack->pmhd->bcc0;
  const Real cap_ad = max_eta_ad;

  using ML = Kokkos::MaxLoc<Real, int>;
  ML::value_type res;
  Kokkos::parallel_reduce("res_rkgdt_bound", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int idx, ML::value_type &v) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    const int k0 = three_d ? k-1 : k, k1 = three_d ? k+1 : k;
    const int j0 = multi_d ? j-1 : j, j1 = multi_d ? j+1 : j;
    // eta_r: the cells of the x2/x3 edges that carry the radial derivative (one transverse
    // offset at a time: the edge eta is the mean of the 4 cells around the edge); eta_h:
    // the whole 3x3x3 block, for the horizontal term (the radial x1 edges at a cube vertex
    // also average the vertex ghost, which no neighbour fills)
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
          if (!diag) {er = fmax(er, e);}
        }
      }
    }
    const Real dcm = x1v_(m,i) - x1v_(m,i-1);
    const Real dcp = x1v_(m,i+1) - x1v_(m,i);
    Real lam = 2.0*(1.0/dcm + 1.0/dcp)/dx1_(m,k,j,i);
    Real h = 0.0;
    if (multi_d) {h += 4.0/SQR(dx2_(m,k,j,i));}
    if (three_d) {
      Real d3 = dx3_(m,k,j,i);
      if (sp) {
        // the phi length at the theta face nearer the pole, not at the centroid; a cell
        // on the pole takes sin(dtheta/2) (sp_test iprob 11: the pole cells' mode is 1.43x
        // the centroid value 4 eta/dx3^2, this bound 1.36x above it)
        const Real s0 = fabs(sin(x2f_(m,j))), s1 = fabs(sin(x2f_(m,j+1)));
        const Real sh = fabs(sin(0.5*(x2f_(m,j+1) - x2f_(m,j))));
        const Real se = fmax(fmin(s0, s1), sh);
        d3 = x1v_(m,i)*(x3f_(m,k+1) - x3f_(m,k))*se;
      }
      h += 4.0/SQR(d3);
    }
    if (cs) {h /= SQR(scell(m,k,j));}
    lam = er*lam + eh*h;
    if (lam > v.val) {v.val = lam; v.loc = idx;}
  }, ML(res));

  Real lmax = res.val;
  Real r = 0.0, e = 0.0;
  if (rloc != nullptr || etaloc != nullptr) {
    if (nmkji > 0) {
      const int idx = res.loc;
      int m = idx/nkji;
      int k = (idx - m*nkji)/nji;
      int j = (idx - m*nkji - k*nji)/nx1;
      int i = (idx - m*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                    Kokkos::subview(x1v_, m, i));
      auto he = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                    Kokkos::subview(eta_, m, k, j, i));
      r = hr();
      e = he();
    }
  }
#if MPI_PARALLEL_ENABLED
  struct {double v; int r;} in, out;
  in.v = static_cast<double>(lmax);
  in.r = global_variable::my_rank;
  MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
  lmax = static_cast<Real>(out.v);
  if (rloc != nullptr || etaloc != nullptr) {
    double re[2] = {static_cast<double>(r), static_cast<double>(e)};
    MPI_Bcast(re, 2, MPI_DOUBLE, out.r, MPI_COMM_WORLD);
    r = re[0];
    e = re[1];
  }
#endif
  if (rloc != nullptr) {*rloc = r;}
  if (etaloc != nullptr) {*etaloc = e;}
  return lmax;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::NewTimeStepRKGFE
//! \brief dtnew = rkg_dt_safety x 2/lambda_bound (local max; mesh.cpp takes it as dt_diff
//! without the CFL factor, and the MPI min over ranks as before).

void Resistivity::NewTimeStepRKGFE() {
  const Real lam = LambdaBoundMax();
  dtnew = (lam > 0.0) ? rkg_dt_safety*2.0/lam
                      : static_cast<Real>(std::numeric_limits<float>::max());
}
