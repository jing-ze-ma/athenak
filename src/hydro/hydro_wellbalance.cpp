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
        

//----------------------------------------------------------------------------------------
//! \fn void Hydro::BuildWBCache
//! \brief walks the local hydrostatic background once per cell and stage and stores the
//! five stencil states of its density, energy and pressure channels in wbq0.
//!
//! Before this, the walk was repeated for the density channel, the energy channel, the
//! pressure channel and once more in the problem generator's source term -- four full
//! stencils per cell per stage, each of them, under a tabulated EOS, a dozen table
//! evaluations.  The reconstruction and the source now read this cache instead.  Built
//! for i in [is-1, ie+1] (the x1 sweep reconstructs those cells) over the rows the sweep
//! covers (passed in); the
//! stencil reaches one cell further, which is inside the ghost zone.  ConsToPrim's cached
//! temperature seeds every inversion, so the isothermal, polytropic and isodensity
//! backgrounds do no root find at all.  Bitwise identical to the uncached path: the same
//! calls in the same order, just once.

void Hydro::BuildWBCache(const int jl, const int ju, const int kl, const int ku) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto eos = peos->eos_data;
  const WBOption wbo = wb_option;
  const bool gen = eos.IsGeneral();
  auto &w0_ = w0;
  auto &phicc = phicc0;
  auto &phi = phi0.x1f;
  auto &wt = wtemp;
  auto &c = wbq0;
  par_for("wbcache", DevExeSpace(), 0, nmb1, kl, ku, jl, ju, is-1, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (gen) {
      // ONE walk gives all three channels
      WBState s_im1, s_imh, s_i, s_iph, s_ip1;
      WBBackgroundStencil(eos, wbo,
              w0_(m,IDN,k,j,i-1), w0_(m,IDN,k,j,i), w0_(m,IDN,k,j,i+1),
              w0_(m,IEN,k,j,i-1), w0_(m,IEN,k,j,i), w0_(m,IEN,k,j,i+1),
              phicc(m,k,j,i-1), phi(m,k,j,i), phicc(m,k,j,i), phi(m,k,j,i+1),
              phicc(m,k,j,i+1), s_im1, s_imh, s_i, s_iph, s_ip1,
              WBT(wt,m,k,j,i-1), WBT(wt,m,k,j,i), WBT(wt,m,k,j,i+1));
      c(m,0,k,j,i) = s_im1.d; c(m,1,k,j,i) = s_imh.d; c(m,2,k,j,i) = s_i.d;
      c(m,3,k,j,i) = s_iph.d; c(m,4,k,j,i) = s_ip1.d;
      c(m,5,k,j,i) = s_im1.e; c(m,6,k,j,i) = s_imh.e; c(m,7,k,j,i) = s_i.e;
      c(m,8,k,j,i) = s_iph.e; c(m,9,k,j,i) = s_ip1.e;
      c(m,10,k,j,i) = s_im1.p; c(m,11,k,j,i) = s_imh.p; c(m,12,k,j,i) = s_i.p;
      c(m,13,k,j,i) = s_iph.p; c(m,14,k,j,i) = s_ip1.p;
    } else {
      // ideal gas: the closed forms, density and energy; pressure is (gamma-1) e
      for (int var=0; var<2; ++var) {
        Real a, b, cc, d, e;
        getWBq0(eos, wbo, var,
                w0_(m,IDN,k,j,i-1), w0_(m,IDN,k,j,i), w0_(m,IDN,k,j,i+1),
                w0_(m,IEN,k,j,i-1), w0_(m,IEN,k,j,i), w0_(m,IEN,k,j,i+1),
                phicc(m,k,j,i-1), phi(m,k,j,i), phicc(m,k,j,i), phi(m,k,j,i+1),
                phicc(m,k,j,i+1), a, b, cc, d, e);
        c(m,5*var,k,j,i) = a; c(m,5*var+1,k,j,i) = b; c(m,5*var+2,k,j,i) = cc;
        c(m,5*var+3,k,j,i) = d; c(m,5*var+4,k,j,i) = e;
      }
      const Real gm1 = eos.gamma - 1.0;
      for (int q=0; q<5; ++q) c(m,10+q,k,j,i) = gm1*c(m,5+q,k,j,i);
    }
  });
  return;
}

} // namespace hydro
