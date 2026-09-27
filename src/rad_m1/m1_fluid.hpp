#ifndef RAD_M1_M1_FLUID_HPP_
#define RAD_M1_M1_FLUID_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file m1_fluid.hpp
//! \brief m1-mhd (docs/dev/m1_mhd_0927.md): ONE accessor for the fluid the grey M1
//! module couples to -- <hydro> or <mhd> -- so that every M1 kernel and the two M1
//! problem generators (box_convection, rad_m1_wedge) run one code path.
//!
//! FluidRef is a host-side value built by FluidRef::Get(pmbp) each time it is needed:
//! shallow copies of the fluid's Views (a Kokkos View assignment copies the handle, not
//! the data) plus its flags and EOS.  Kernels copy the members they need into locals
//! exactly as before (`auto uh = fl.u0;`), so a hydro run executes the same arithmetic
//! on the same memory (bitwise).  Under MHD the conserved energy also carries |B|^2/2:
//! the M1 kernels subtract it with the cell field the MHD C2P builds (M1EmagCell).

#include <utility>

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

namespace radm1 {

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
  }
  // swap the OWNER's u0 and u1 (the handles; used to evaluate a U^n-state quantity with
  // a routine that reads u0).  The FluidRef copies are swapped too.
  void SwapU01() {
    if (pm != nullptr) {
      std::swap(pm->u0, pm->u1);
    } else if (ph != nullptr) {
      std::swap(ph->u0, ph->u1);
    }
    std::swap(u0, u1);
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
//! \brief |B|^2/2 of the cell field built from the faces EXACTLY as the MHD C2P builds it
//! (general_mhd.cpp / ideal_mhd.cpp ConsToPrim): spherical polar = the x_v-weighted mean
//! in all three directions; cubed sphere = CellCenteredRadialFld in x1 and plain means in
//! x2, x3; otherwise plain face means.  The temperature M1 forms is then the C2P's.

KOKKOS_INLINE_FUNCTION
Real M1EmagCell(const DvceFaceFld4D<Real> &b, const bool sph, const bool csr,
                const DvceArray2D<Real> &x1v, const DvceArray2D<Real> &x1f,
                const DvceArray2D<Real> &x2v, const DvceArray2D<Real> &x2f,
                const DvceArray2D<Real> &x3v, const DvceArray2D<Real> &x3f,
                const int m, const int k, const int j, const int i) {
  Real bx, by, bz;
  if (sph) {
    Real lw, rw;
    lw = (x1f(m,i+1)-x1v(m,i))/(x1f(m,i+1)-x1f(m,i));
    rw = (x1v(m,i)-x1f(m,i))/(x1f(m,i+1)-x1f(m,i));
    bx = lw*b.x1f(m,k,j,i) + rw*b.x1f(m,k,j,i+1);
    lw = (x2f(m,j+1)-x2v(m,j))/(x2f(m,j+1)-x2f(m,j));
    rw = (x2v(m,j)-x2f(m,j))/(x2f(m,j+1)-x2f(m,j));
    by = lw*b.x2f(m,k,j,i) + rw*b.x2f(m,k,j+1,i);
    lw = (x3f(m,k+1)-x3v(m,k))/(x3f(m,k+1)-x3f(m,k));
    rw = (x3v(m,k)-x3f(m,k))/(x3f(m,k+1)-x3f(m,k));
    bz = lw*b.x3f(m,k,j,i) + rw*b.x3f(m,k+1,j,i);
  } else {
    bx = csr ? CellCenteredRadialFld(b.x1f(m,k,j,i), b.x1f(m,k,j,i+1),
                                     x1f(m,i), x1f(m,i+1), x1v(m,i))
             : 0.5*(b.x1f(m,k,j,i) + b.x1f(m,k,j,i+1));
    by = 0.5*(b.x2f(m,k,j,i) + b.x2f(m,k,j+1,i));
    bz = 0.5*(b.x3f(m,k,j,i) + b.x3f(m,k+1,j,i));
  }
  return 0.5*(bx*bx + by*by + bz*bz);
}

} // namespace radm1

#endif // RAD_M1_M1_FLUID_HPP_
