//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file fofc_etotgrav.cpp
//! \brief the device kernel of the FOFC etotgrav trial correction (see fofc_etotgrav.hpp)

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/fofc_etotgrav.hpp"

//! utest(IEN) -= div(F_rho*Phi_f) + utest(IDN)*Phi_c over [il,iu]x[jl,ju]x[kl,ku]
void FofcTrialRemoveGrav(MeshBlockPack *pp, DvceArray5D<Real> utest,
                         const DvceFaceFld5D<Real> &flx, const DvceArray4D<Real> &phicc,
                         const DvceFaceFld4D<Real> &phif, const Real beta_dt,
                         const int il, const int iu, const int jl, const int ju,
                         const int kl, const int ku) {
  const bool curv = pp->pmesh->use_cubed_sphere || pp->pmesh->use_spherical_polar;
  const bool multi_d = pp->pmesh->multi_d;
  const bool three_d = pp->pmesh->three_d;
  auto vol = pp->pcoord->volume;
  auto area = pp->pcoord->area;
  auto &size = pp->pmb->mb_size;
  auto flx_ = flx;
  auto phif_ = phif;
  auto phicc_ = phicc;
  par_for("FOFC-etotgrav", DevExeSpace(), 0, pp->nmb_thispack-1, kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real divg = FofcGravDiv(m, k, j, i, curv, multi_d, three_d, beta_dt, flx_,
                                  phif_, vol, area, size.d_view(m).dx1,
                                  size.d_view(m).dx2, size.d_view(m).dx3);
    utest(m,IEN,k,j,i) -= divg + utest(m,IDN,k,j,i)*phicc_(m,k,j,i);
  });
}

