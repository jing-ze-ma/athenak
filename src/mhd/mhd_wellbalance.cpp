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
//! \fn void MHD::SetWbBackgroundPressure
//! \brief evaluates the gas pressure of the static well-balanced background, once.
//!
//! The background (w0wb, w0facewb) is handed over by the problem generator and does not
//! evolve, so its pressure never changes either. Evaluating it here, at initialization,
//! keeps the EOS out of the deviation reconstruction, the flux correction and the
//! coordinate source terms, each of which would otherwise ask for the background pressure
//! in every cell of every stage -- of order twenty calls per cell per stage, every one of
//! them a root find once the EOS stops being a gamma law.
//!
//! Filled over the FULL arrays including ghost zones: the deviation reconstruction reads
//! the background at i-1 and i+1 across the whole tile, not just the active zone.

void MHD::SetWbBackgroundPressure() {
  if (!use_wellbalance_static) {return;}
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  auto eos = peos->eos_data;

  auto &wc = w0wb;
  auto &pc = pwb;
  par_for("wbsetpres", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    pc(m,k,j,i) = eos.Pressure(wc(m,IDN,k,j,i), wc(m,IEN,k,j,i));
  });

  auto &wf1 = w0facewb.x1f;
  auto &pf1 = pfacewb.x1f;
  par_for("wbsetpresf1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    pf1(m,k,j,i) = eos.Pressure(wf1(m,IDN,k,j,i), wf1(m,IEN,k,j,i));
  });

  auto &wf2 = w0facewb.x2f;
  auto &pf2 = pfacewb.x2f;
  par_for("wbsetpresf2", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1+1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    pf2(m,k,j,i) = eos.Pressure(wf2(m,IDN,k,j,i), wf2(m,IEN,k,j,i));
  });

  auto &wf3 = w0facewb.x3f;
  auto &pf3 = pfacewb.x3f;
  par_for("wbsetpresf3", DevExeSpace(), 0, nmb1, 0, n3m1+1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    pf3(m,k,j,i) = eos.Pressure(wf3(m,IDN,k,j,i), wf3(m,IEN,k,j,i));
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MHD::RemoveWbFlux
//! \brief Removes the background state flux from the total flux.

void MHD::RemoveWbFlux(const DvceFaceFld4D<Real> &pfacewb, DvceFaceFld5D<Real> &flx) {
     auto &indcs = pmy_pack->pmesh->mb_indcs;
     int is = indcs.is, ie = indcs.ie;
     int js = indcs.js, je = indcs.je;
     int ks = indcs.ks, ke = indcs.ke;
     int nmb1 = pmy_pack->nmb_thispack - 1;
     auto size = pmy_pack->pmb->mb_size;
    
     // The background pressure is NOT (gamma-1)*e_bg -- the background primitives carry
     // their own density -- but it is also not asked of the EOS here: the background is
     // static, so SetWbBackgroundPressure() evaluated it once at startup.

     //--------------------------------------------------------------------------------------
     // fluxes in x1-direction

     auto &flx1 = flx.x1f;
     auto &pfacewb1 = pfacewb.x1f;

     par_for("wbremflux1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
       flx1(m,IM1,k,j,i) -= pfacewb1(m,k,j,i);
     });
     if (pmy_pack->pmesh->one_d) {return;}

     //--------------------------------------------------------------------------------------
     // fluxes in x2-direction

     auto &flx2 = flx.x2f;
     auto &pfacewb2 = pfacewb.x2f;

     par_for("wbremflux2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
         flx2(m,IM2,k,j,i) -= pfacewb2(m,k,j,i);
     });
     if (pmy_pack->pmesh->two_d) {return;}

     //--------------------------------------------------------------------------------------
     // fluxes in x3-direction

     auto &flx3 = flx.x3f;
     auto &pfacewb3 = pfacewb.x3f;

     par_for("wbremflux3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
     KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
         flx3(m,IM3,k,j,i) -= pfacewb3(m,k,j,i);
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
