#ifndef HYDRO_FOFC_ETOTGRAV_HPP_
#define HYDRO_FOFC_ETOTGRAV_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file fofc_etotgrav.hpp
//! \brief the FOFC trial state under <hydro|mhd>/etotgrav: turn its energy into
//! e + KE (+ emag) before the floor test.  Shared by Hydro::FOFC and MHD::FOFC.
//!
//! With etotgrav the conserved energy is E = e + KE (+ emag) + rho*Phi_c and the real
//! update adds the gravitational energy flux F_rho*Phi_f to every face (AddGravFlux,
//! applied AFTER FOFC).  The trial state utest = gam0*u0 + gam1*u1 - div(F_Riemann) has
//! neither: its energy is
//!     E_t = (e + KE)_pre + rho_pre*Phi_c - div(F_E),
//! so handing it to ConsToPrim reads rho_pre*Phi_c as internal energy (with Phi >= 0,
//! red_giant/box, an energy or temperature floor hit is never detected), and removing
//! only rho_t*Phi_c leaves Phi_c*div(F_rho) -- the potential energy of the mass the
//! stage moves, which where rho*Phi >> e is far larger than e and flags healthy cells.
//! The consistent trial energy is the one the real update produces from these fluxes:
//!     E_t - div(F_rho*Phi_f) - rho_t*Phi_c
//!       = (e + KE)_pre - div(F_E) - [div(F_rho*Phi_f) - Phi_c*div(F_rho)],
//! i.e. the Riemann update plus the gravitational work, with the same face areas and
//! volumes (curvilinear) or widths (Cartesian) as the trial divergence.  utest is
//! scratch; nothing is re-added.  With etotgrav off this is never called.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"

KOKKOS_INLINE_FUNCTION
Real FofcGravDiv(const int m, const int k, const int j, const int i, const bool curv,
                 const bool multi_d, const bool three_d, const Real beta_dt,
                 const DvceFaceFld5D<Real> &flx, const DvceFaceFld4D<Real> &phi,
                 const DvceArray4D<Real> &vol, const DvceFaceFld4D<Real> &area,
                 const Real dx1, const Real dx2, const Real dx3) {
  const int n = IDN;
  Real div;
  if (curv) {
    const Real dtodv = beta_dt/vol(m,k,j,i);
    div = dtodv*(flx.x1f(m,n,k,j,i+1)*phi.x1f(m,k,j,i+1)*area.x1f(m,k,j,i+1) -
                 flx.x1f(m,n,k,j,i  )*phi.x1f(m,k,j,i  )*area.x1f(m,k,j,i  ));
    if (multi_d) {
      div += dtodv*(flx.x2f(m,n,k,j+1,i)*phi.x2f(m,k,j+1,i)*area.x2f(m,k,j+1,i) -
                    flx.x2f(m,n,k,j  ,i)*phi.x2f(m,k,j  ,i)*area.x2f(m,k,j  ,i));
    }
    if (three_d) {
      div += dtodv*(flx.x3f(m,n,k+1,j,i)*phi.x3f(m,k+1,j,i)*area.x3f(m,k+1,j,i) -
                    flx.x3f(m,n,k  ,j,i)*phi.x3f(m,k  ,j,i)*area.x3f(m,k  ,j,i));
    }
  } else {
    div = (beta_dt/dx1)*(flx.x1f(m,n,k,j,i+1)*phi.x1f(m,k,j,i+1) -
                         flx.x1f(m,n,k,j,i  )*phi.x1f(m,k,j,i  ));
    if (multi_d) {
      div += (beta_dt/dx2)*(flx.x2f(m,n,k,j+1,i)*phi.x2f(m,k,j+1,i) -
                            flx.x2f(m,n,k,j  ,i)*phi.x2f(m,k,j  ,i));
    }
    if (three_d) {
      div += (beta_dt/dx3)*(flx.x3f(m,n,k+1,j,i)*phi.x3f(m,k+1,j,i) -
                            flx.x3f(m,n,k  ,j,i)*phi.x3f(m,k  ,j,i));
    }
  }
  return div;
}

//! utest(IEN) -= div(F_rho*Phi_f) + utest(IDN)*Phi_c over [il,iu]x[jl,ju]x[kl,ku]
//! Defined once in fofc_etotgrav.cpp: an inline header function holding a KOKKOS_LAMBDA
//! compiled into two TUs (hydro_fofc, mhd_fofc) breaks under nvcc (one linker copy, the
//! other caller jumps through an uninitialised lambda wrapper; see c2p_track).
void FofcTrialRemoveGrav(MeshBlockPack *pp, DvceArray5D<Real> utest,
                         const DvceFaceFld5D<Real> &flx, const DvceArray4D<Real> &phicc,
                         const DvceFaceFld4D<Real> &phif, const Real beta_dt,
                         const int il, const int iu, const int jl, const int ju,
                         const int kl, const int ku);

#endif  // HYDRO_FOFC_ETOTGRAV_HPP_
