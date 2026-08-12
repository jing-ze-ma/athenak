//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_wellbalance.cpp
//! \brief Implements functions for deviation-based well-balanced scheme.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro.hpp"

namespace hydro {
//----------------------------------------------------------------------------------------
//! \fn void Hydro::SetWbBackgroundPressure
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

void Hydro::SetWbBackgroundPressure() {
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
//! \fn void Hydro::RemoveWbFlux
//! \brief Removes the background state flux from the total flux.

void Hydro::RemoveWbFlux(const DvceFaceFld4D<Real> &pfacewb, DvceFaceFld5D<Real> &flx) {
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
//! \fn void Hydro::AddWbVar
//! \brief Adds the background variables onto perturbed variables.

void Hydro::AddWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int nmb1 = pmy_pack->nmb_thispack - 1;
    int &ng = indcs.ng;
    int n1m1 = indcs.nx1 + 2*ng - 1;
    int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    int nvar = nhydro;

    par_for("wbaddvar", DevExeSpace(), 0, nmb1, 0, nvar-1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
      var(m,n,k,j,i) += varwb(m,n,k,j,i);
    });

     return;
}
    
//----------------------------------------------------------------------------------------
//! \fn void Hydro::RemoveWbVar
//! \brief Removes the background variables onto perturbed variables.

void Hydro::RemoveWbVar(const DvceArray5D<Real> &varwb, DvceArray5D<Real> &var) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    int nmb1 = pmy_pack->nmb_thispack - 1;
    int &ng = indcs.ng;
    int n1m1 = indcs.nx1 + 2*ng - 1;
    int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
    int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    int nvar = nhydro;

    par_for("wbaddvar", DevExeSpace(), 0, nmb1, 0, nvar-1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int n, int k, int j, int i) {
      var(m,n,k,j,i) -= varwb(m,n,k,j,i);
    });

     return;
}

////----------------------------------------------------------------------------------------
////! \fn AddWbPrimFace()
////! \brief Adds background face-centered variables onto ql(i+1) and qr(i).
//
//KOKKOS_INLINE_FUNCTION
//void Hydro::AddWbPrimFace(const Real &qlwb_ip1, const Real &qrwb_i,
//         Real &ql_ip1, Real &qr_i) {
//  ql_ip1 += qlwb_ip1;
//  qr_i   += qrwb_i;
//  return;
//}
//
////----------------------------------------------------------------------------------------
////! \fn AddWbPrimFaceX1()
////! \brief Adds background states onto face-centered primitive variables in x1-direction.
////! This function should be called over [is-1,ie+1] to get BOTH L/R states over [is,ie]
//
//KOKKOS_INLINE_FUNCTION
//void Hydro::AddWbPrimFaceX1(TeamMember_t const &member, const int m, const int k, const int j,
//     const int il, const int iu, const DvceArray5D<Real> &q,
//     ScrArray2D<Real> &ql, ScrArray2D<Real> &qr) {
//  int nvar = q.extent_int(1);
//  for (int n=0; n<nvar; ++n) {
//    par_for_inner(member, il, iu, [&](const int i) {
//      AddWbPrimFace(q(m,n,k,j,i+1), q(m,n,k,j,i), ql(n,i+1), qr(n,i));
//    });
//  }
//  return;
//}
//
////----------------------------------------------------------------------------------------
////! \fn AddWbPrimFaceX2()
////! \brief Adds background states onto face-centered primitive variables in x2-direction.
////! This function should be called over [js-1,je+1] to get BOTH L/R states over [js,je]
//
//KOKKOS_INLINE_FUNCTION
//void Hydro::AddWbPrimFaceX2(TeamMember_t const &member, const int m, const int k, const int j,
//     const int il, const int iu, const DvceArray5D<Real> &q,
//     ScrArray2D<Real> &ql_jp1, ScrArray2D<Real> &qr_j) {
//  int nvar = q.extent_int(1);
//  for (int n=0; n<nvar; ++n) {
//    par_for_inner(member, il, iu, [&](const int i) {
//      AddWbPrimFace(q(m,n,k,j+1,i), q(m,n,k,j,i), ql_jp1(n,i), qr_j(n,i));
//    });
//  }
//  return;
//}
//
////----------------------------------------------------------------------------------------
////! \fn AddWbPrimFaceX3()
////! \brief Adds background states onto face-centered primitive variables in x3-direction.
////! This function should be called over [ks-1,ke+1] to get BOTH L/R states over [ks,ke]
//
//KOKKOS_INLINE_FUNCTION
//void Hydro::AddWbPrimFaceX3(TeamMember_t const &member, const int m, const int k, const int j,
//     const int il, const int iu, const DvceArray5D<Real> &q,
//     ScrArray2D<Real> &ql_kp1, ScrArray2D<Real> &qr_k) {
//  int nvar = q.extent_int(1);
//  for (int n=0; n<nvar; ++n) {
//    par_for_inner(member, il, iu, [&](const int i) {
//      AddWbPrimFace(q(m,n,k+1,j,i), q(m,n,k,j,i), ql_kp1(n,i), qr_k(n,i));
//    });
//  }
//  return;
//}
        
} // namespace hydro
