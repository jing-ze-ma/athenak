#ifndef RAD_VET_VET_FLUID_HPP_
#define RAD_VET_VET_FLUID_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file vet_fluid.hpp
//! \brief FluidRef: one accessor for the fluid M1 couples to (<hydro> or <mhd>).
//! FluidRef::Get(pmbp) returns shallow View handles plus flags and EOS. Under MHD the
//! conserved energy carries |B|^2/2, removed with the C2P cell field (M1EmagCell).

#include <utility>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

namespace radvet {

struct FluidRef {
  hydro::Hydro *ph = nullptr;
  mhd::MHD *pm = nullptr;
  bool on = false;             // a fluid exists
  bool mhd = false;            // it is <mhd>
  DvceArray5D<Real> u0, u1, w0;
  DvceArray4D<Real> wtemp;
  DvceFaceFld5D<Real> *uflx = nullptr;   // face fields have no default ctor: pointers
  DvceArray4D<Real> phicc0;
  DvceFaceFld4D<Real> *phi0 = nullptr;
  DvceArray4D<Real> phicc_wb, phi_wb_x1f;
  DvceArray5D<Real> wbq0;
  EOS_Data eos;
  bool use_etotgrav = false;
  bool use_wellbalance_dynamic = false;
  bool use_wb_x1 = false;
  bool use_phi_wb = false;
  WBOption wb_option;
  Conduction *pcond = nullptr;
  Viscosity *pvisc = nullptr;
  EquationOfState *peos = nullptr;

  static FluidRef Get(MeshBlockPack *pp) {
    FluidRef f;
    if (pp->pmhd != nullptr) {
      f.pm = pp->pmhd;
      f.Fill(pp->pmhd);
      f.mhd = true;
    } else if (pp->phydro != nullptr) {
      f.ph = pp->phydro;
      f.Fill(pp->phydro);
    }
    return f;
  }
  template <class F>
  void Fill(F *p) {
    on = true;
    u0 = p->u0; u1 = p->u1; w0 = p->w0; wtemp = p->wtemp; uflx = &(p->uflx);
    phicc0 = p->phicc0; phi0 = &(p->phi0); phicc_wb = p->phicc_wb;
    phi_wb_x1f = p->phi_wb_x1f; wbq0 = p->wbq0;
    peos = p->peos; eos = p->peos->eos_data;
    use_etotgrav = p->use_etotgrav;
    use_wellbalance_dynamic = p->use_wellbalance_dynamic;
    use_wb_x1 = p->use_wb_x1;
    use_phi_wb = p->use_phi_wb;
    wb_option = p->wb_option;
    pcond = p->pcond;
    pvisc = p->pvisc;
  }
  void EnableWBEffectivePotential() {
    if (pm != nullptr) {
      pm->EnableWBEffectivePotential();
      phicc_wb = pm->phicc_wb; phi_wb_x1f = pm->phi_wb_x1f; use_phi_wb = true;
    } else if (ph != nullptr) {
      ph->EnableWBEffectivePotential();
      phicc_wb = ph->phicc_wb; phi_wb_x1f = ph->phi_wb_x1f; use_phi_wb = true;
    }
  }
};

//----------------------------------------------------------------------------------------
//! \fn Real M1EmagCell
//! \brief |B|^2/2 of the cell field, averaged from faces exactly as the MHD ConsToPrim
//! does (sph: x_v-weighted; cs: CellCenteredRadialFld in x1; else plain means).
//! A hydro caller may pass any 4-D View as never-read face arrays.

KOKKOS_INLINE_FUNCTION
Real M1EmagCell(const DvceArray4D<Real> &b1f, const DvceArray4D<Real> &b2f,
                const DvceArray4D<Real> &b3f, const bool sph, const bool csr,
                const DvceArray2D<Real> &x1v, const DvceArray2D<Real> &x1f,
                const DvceArray2D<Real> &x2v, const DvceArray2D<Real> &x2f,
                const DvceArray2D<Real> &x3v, const DvceArray2D<Real> &x3f,
                const int m, const int k, const int j, const int i) {
  Real bx, by, bz;
  if (sph) {
    Real lw, rw;
    lw = (x1f(m,i+1)-x1v(m,i))/(x1f(m,i+1)-x1f(m,i));
    rw = (x1v(m,i)-x1f(m,i))/(x1f(m,i+1)-x1f(m,i));
    bx = lw*b1f(m,k,j,i) + rw*b1f(m,k,j,i+1);
    lw = (x2f(m,j+1)-x2v(m,j))/(x2f(m,j+1)-x2f(m,j));
    rw = (x2v(m,j)-x2f(m,j))/(x2f(m,j+1)-x2f(m,j));
    by = lw*b2f(m,k,j,i) + rw*b2f(m,k,j+1,i);
    lw = (x3f(m,k+1)-x3v(m,k))/(x3f(m,k+1)-x3f(m,k));
    rw = (x3v(m,k)-x3f(m,k))/(x3f(m,k+1)-x3f(m,k));
    bz = lw*b3f(m,k,j,i) + rw*b3f(m,k+1,j,i);
  } else {
    bx = csr ? CellCenteredRadialFld(b1f(m,k,j,i), b1f(m,k,j,i+1),
                                     x1f(m,i), x1f(m,i+1), x1v(m,i))
             : 0.5*(b1f(m,k,j,i) + b1f(m,k,j,i+1));
    by = 0.5*(b2f(m,k,j,i) + b2f(m,k,j+1,i));
    bz = 0.5*(b3f(m,k,j,i) + b3f(m,k+1,j,i));
  }
  return 0.5*(bx*bx + by*by + bz*bz);
}

} // namespace radvet

#endif // RAD_VET_VET_FLUID_HPP_
