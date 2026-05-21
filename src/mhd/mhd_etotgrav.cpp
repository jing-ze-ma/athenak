//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_etotgrav.cpp
//! \brief Implements functions for including time-independent gravitational potential energy in the energy conservation.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "mhd.hpp"

namespace mhd {
//----------------------------------------------------------------------------------------
//! \fn void MHD::AddGravFlux
//! \brief Implements the gravity term in the energy flux.

void MHD::AddGravFlux(const DvceFaceFld4D<Real> &phi0, DvceFaceFld5D<Real> &flx) {
     auto &indcs = pmy_pack->pmesh->mb_indcs;
     int is = indcs.is, ie = indcs.ie;
     int js = indcs.js, je = indcs.je;
     int ks = indcs.ks, ke = indcs.ke;
     int nmb1 = pmy_pack->nmb_thispack - 1;
     auto size = pmy_pack->pmb->mb_size;

     //--------------------------------------------------------------------------------------
     // fluxes in x1-direction

     auto &flx1 = flx.x1f;
     auto &phi1 = phi0.x1f;

     par_for("etotgrav1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
       flx1(m,IEN,k,j,i) += flx1(m,IDN,k,j,i) * phi1(m,k,j,i);
     });
     if (pmy_pack->pmesh->one_d) {return;}

     //--------------------------------------------------------------------------------------
     // fluxes in x2-direction

     auto &flx2 = flx.x2f;
     auto &phi2 = phi0.x2f;

     par_for("etotgrav2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        flx2(m,IEN,k,j,i) += flx2(m,IDN,k,j,i) * phi2(m,k,j,i);
     });
     if (pmy_pack->pmesh->two_d) {return;}

     //--------------------------------------------------------------------------------------
     // fluxes in x3-direction

     auto &flx3 = flx.x3f;
     auto &phi3 = phi0.x3f;

     par_for("etotgrav3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        flx3(m,IEN,k,j,i) += flx3(m,IDN,k,j,i) * phi3(m,k,j,i);
     });

     return;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::AddGravEtot
//! \brief Adds the gravitational potential energy in the total energy.

void MHD::AddGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku) {
    int nmb1 = pmy_pack->nmb_thispack - 1;

    par_for("etotgravadd", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      cons(m,IEN,k,j,i) += cons(m,IDN,k,j,i)*phicc0(m,k,j,i);
    });

     return;
}
    
//----------------------------------------------------------------------------------------
//! \fn void MHD::RemoveGravEtot
//! \brief Removes the gravitational potential energy from the total energy.

void MHD::RemoveGravEtot(const DvceArray4D<Real> &phicc0, DvceArray5D<Real> &cons, const int il, const int iu, const int jl, const int ju, const int kl, const int ku) {
    int nmb1 = pmy_pack->nmb_thispack - 1;
        
    par_for("etotgravrem", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, il, iu,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      cons(m,IEN,k,j,i) -= cons(m,IDN,k,j,i)*phicc0(m,k,j,i);
    });

     return;
}
        
} // namespace mhd
