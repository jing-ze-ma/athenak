//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_wellbalance.cpp
//! \brief Implements functions for deviation-based well-balanced scheme.

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
//! \fn void MHD::RemoveWbFlux
//! \brief Removes the background state flux from the total flux.

void MHD::RemoveWbFlux(const DvceFaceFld5D<Real> &w0facewb, DvceFaceFld5D<Real> &flx) {
     auto &indcs = pmy_pack->pmesh->mb_indcs;
     int is = indcs.is, ie = indcs.ie;
     int js = indcs.js, je = indcs.je;
     int ks = indcs.ks, ke = indcs.ke;
     int nmb1 = pmy_pack->nmb_thispack - 1;
     auto size = pmy_pack->pmb->mb_size;
    
     Real gamma = peos->eos_data.gamma;
     Real gm1 = gamma - 1.0;

     //--------------------------------------------------------------------------------------
     // fluxes in x1-direction

     auto &flx1 = flx.x1f;
     auto &w0facewb1 = w0facewb.x1f;

     par_for("wbremflux1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
       flx1(m,IM1,k,j,i) -= w0facewb1(m,IEN,k,j,i) * gm1;
     });
     if (pmy_pack->pmesh->one_d) {return;}

     //--------------------------------------------------------------------------------------
     // fluxes in x2-direction

     auto &flx2 = flx.x2f;
     auto &w0facewb2 = w0facewb.x2f;

     par_for("wbremflux2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
         flx2(m,IM2,k,j,i) -= w0facewb2(m,IEN,k,j,i) * gm1;
     });
     if (pmy_pack->pmesh->two_d) {return;}

     //--------------------------------------------------------------------------------------
     // fluxes in x3-direction

     auto &flx3 = flx.x3f;
     auto &w0facewb3 = w0facewb.x3f;

     par_for("wbremflux3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
         flx3(m,IM3,k,j,i) -= w0facewb3(m,IEN,k,j,i) * gm1;
     });

     return;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::AddWbVar
//! \brief Adds the background variables onto perturbed variables.

void MHD::AddWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int nmb1 = pmy_pack->nmb_thispack - 1;
    int &ng = indcs.ng;
    int n1m1 = indcs.nx1 + 2*ng - 1;
    int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    int nvar = nmhd;

    par_for("wbaddvar", DevExeSpace(), 0, nmb1, 0, nvar-1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
      var(m,n,k,j,i) += varwb(m,n,k,j,i);
    });

     return;
}
    
//----------------------------------------------------------------------------------------
//! \fn void MHD::RemoveWbVar
//! \brief Removes the background variables onto perturbed variables.

void MHD::RemoveWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int nmb1 = pmy_pack->nmb_thispack - 1;
    int &ng = indcs.ng;
    int n1m1 = indcs.nx1 + 2*ng - 1;
    int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    int nvar = nmhd;

    par_for("wbaddvar", DevExeSpace(), 0, nmb1, 0, nvar-1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
      var(m,n,k,j,i) -= varwb(m,n,k,j,i);
    });

     return;
}
        
} // namespace mhd
