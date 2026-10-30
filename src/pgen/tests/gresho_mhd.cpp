//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gresho_mhd.cpp
//! \brief Problem generator for a magnetized Gresho vortex at low Mach number.
//!
//! The Gresho vortex is a time-independent solution of the Euler equations: the
//! centrifugal force of a rotating ring of gas is balanced by a radial pressure
//! gradient.  Adding a UNIFORM field leaves the equilibrium untouched (a uniform field
//! exerts no force), so the exact solution is stationary in MHD as well, and any decay
//! of the kinetic energy is numerical dissipation.  Lowering the Mach number (i.e.
//! raising the background pressure at fixed peak rotation speed) makes this a standard
//! test of how a Riemann solver behaves in the low Mach number limit: the dissipation of
//! a standard Godunov solver grows like 1/M, and the vortex is destroyed in a turnover.
//!
//! The azimuthal velocity and pressure are (Gresho & Chan 1990; Liska & Wendroff 2003)
//!   u_phi = 5 r,                 p = p0 + 12.5 r^2                      (r < 0.2)
//!   u_phi = 2 - 5 r,             p = p0 + 12.5 r^2 + 4 - 20 r + 4 ln(5r)(0.2 < r < 0.4)
//!   u_phi = 0,                   p = p0 - 2 + 4 ln 2                    (r > 0.4)
//! with rho = d0, so that max(u_phi) = 1 at r = 0.2 and the core rotates rigidly with
//! angular frequency 5 (one turnover in t = 2 pi/5 = 1.2566).  The background pressure
//! follows from the requested Mach number, p0 = d0/(gamma M^2), and the uniform field
//! along x1 from the requested plasma beta, B0 = sqrt(2 p0/beta).
//!
//! REFERENCES:
//! - P. M. Gresho & S. T. Chan, "On the theory of semi-implicit projection methods...",
//!   Int. J. Numer. Methods Fluids, 11, 621 (1990)
//! - T. Minoshima & T. Miyoshi, "A low-dissipation HLLD approximate Riemann solver for a
//!   very wide range of Mach numbers", JCP, 446, 110639 (2021)

// C++ headers
#include <math.h>
#include <cstdint>
#include <iostream>   // endl

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real GreshoVphi(const Real r)
//! \brief azimuthal velocity of the Gresho vortex

KOKKOS_INLINE_FUNCTION
Real GreshoVphi(const Real r) {
  if (r < 0.2) {
    return 5.0*r;
  } else if (r < 0.4) {
    return 2.0 - 5.0*r;
  }
  return 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real GreshoPres(const Real r, const Real p0)
//! \brief pressure of the Gresho vortex, in hydrostatic balance with GreshoVphi

KOKKOS_INLINE_FUNCTION
Real GreshoPres(const Real r, const Real p0) {
  if (r < 0.2) {
    return p0 + 12.5*r*r;
  } else if (r < 0.4) {
    return p0 + 12.5*r*r + 4.0 - 20.0*r + 4.0*std::log(5.0*r);
  }
  return p0 - 2.0 + 4.0*std::log(2.0);
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::GreshoMHD(ParameterInput *pin, const bool restart)
//! \brief Problem generator for the magnetized Gresho vortex.  Assumes the domain is
//! [-0.5,0.5] x [-0.5,0.5], so that the vortex is centered on the origin.

void ProblemGenerator::GreshoMHD(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // lowmach-noise-1030: also runs in hydro (<hydro> block, no field)
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (!is_mhd && pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Gresho vortex test needs an <mhd> or a <hydro> block" << std::endl;
    exit(EXIT_FAILURE);
  }

  EOS_Data &eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  Real d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  Real mach = pin->GetOrAddReal("problem", "mach", 0.01);
  Real beta = pin->GetOrAddReal("problem", "beta", 100.0);
  // peak rotation speed is 1, so M = 1/c_s fixes the background pressure
  Real p0 = d0/(eos.gamma*mach*mach);
  Real b0mag = is_mhd ? std::sqrt(2.0*p0/beta) : 0.0;
  // lowmach-noise-1030 (TEST keys, read only when set): uniform field at angle b_angle
  // (degrees) from x1 in the x1-x2 plane (default 0 = along x1, as before), and a
  // grid-scale velocity seed vpert_grid*(-1)^(i+j) (in units of the peak speed 1) on
  // both velocity components, or with vpert_rand > 0 a white-noise seed of that rms
  // (cell-hashed, reproducible)
  Real bang = 0.0, vgrid = 0.0, vrand = 0.0;
  if (pin->DoesParameterExist("problem", "b_angle")) {
    bang = pin->GetReal("problem", "b_angle")*M_PI/180.0;
  }
  if (pin->DoesParameterExist("problem", "vpert_grid")) {
    vgrid = pin->GetReal("problem", "vpert_grid");
  }
  if (pin->DoesParameterExist("problem", "vpert_rand")) {
    vrand = pin->GetReal("problem", "vpert_rand");
  }
  Real b1 = b0mag*std::cos(bang), b2 = b0mag*std::sin(bang), b3 = 0.0;
  // ... and with b_incl (degrees from x3) set, a field out of the x1-x2 plane, e.g.
  // b_incl = 0: B along x3, perpendicular to every face normal of the 2-D mesh
  if (pin->DoesParameterExist("problem", "b_incl")) {
    const Real binc = pin->GetReal("problem", "b_incl")*M_PI/180.0;
    b1 = b0mag*std::sin(binc)*std::cos(bang);
    b2 = b0mag*std::sin(binc)*std::sin(bang);
    b3 = b0mag*std::cos(binc);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &size = pmbp->pmb->mb_size;

  par_for("pgen_gresho", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real rad = std::sqrt(x1v*x1v + x2v*x2v);
    Real vphi = GreshoVphi(rad);
    Real vx = 0.0, vy = 0.0;
    if (rad > 0.0) {
      vx = -vphi*x2v/rad;
      vy =  vphi*x1v/rad;
    }
    if (vgrid != 0.0) {
      // global parity from the cell centre (same on every MeshBlock layout)
      const int gi = static_cast<int>(std::floor((x1v + 1000.0)*nx1/(x1max - x1min)));
      const int gj = static_cast<int>(std::floor((x2v + 1000.0)*nx2/(x2max - x2min)));
      const Real sgn = ((gi + gj) % 2 == 0) ? 1.0 : -1.0;
      vx += vgrid*sgn;
      vy += vgrid*sgn;
    }
    if (vrand != 0.0) {
      // cell-hashed uniform deviates in [-1,1), rms 1/sqrt(3)
      const uint64_t gi = static_cast<uint64_t>((x1v + 1000.0)*nx1/(x1max - x1min));
      const uint64_t gj = static_cast<uint64_t>((x2v + 1000.0)*nx2/(x2max - x2min));
      uint64_t h = (gi*0x9E3779B97F4A7C15ULL)
                   ^ ((gj + 0x632BE59BD9B4E019ULL)*0xBF58476D1CE4E5B9ULL);
      h ^= h >> 31; h *= 0x94D049BB133111EBULL; h ^= h >> 29;
      const Real r1 = static_cast<Real>(h >> 11)*(1.0/9007199254740992.0);
      h *= 0xD6E8FEB86659FD93ULL; h ^= h >> 32;
      const Real r2 = static_cast<Real>(h >> 11)*(1.0/9007199254740992.0);
      vx += vrand*std::sqrt(3.0)*(2.0*r1 - 1.0);
      vy += vrand*std::sqrt(3.0)*(2.0*r2 - 1.0);
    }

    u0(m,IDN,k,j,i) = d0;
    u0(m,IM1,k,j,i) = d0*vx;
    u0(m,IM2,k,j,i) = d0*vy;
    u0(m,IM3,k,j,i) = 0.0;
    // total energy, with the uniform field contributing 0.5*B0^2 everywhere
    u0(m,IEN,k,j,i) = GreshoPres(rad,p0)/gm1 + 0.5*d0*(vx*vx + vy*vy)
                      + 0.5*b0mag*b0mag;
  });
  if (!is_mhd) return;

  // uniform field at angle b_angle from x1 in the x1-x2 plane (divergence-free)
  auto &b0 = pmbp->pmhd->b0;
  par_for("pgen_gresho_b", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m,k,j,i) = b1;
    b0.x2f(m,k,j,i) = b2;
    b0.x3f(m,k,j,i) = b3;
    if (i==ie) {
      b0.x1f(m,k,j,i+1) = b1;
    }
    if (j==je) {
      b0.x2f(m,k,j+1,i) = b2;
    }
    if (k==ke) {
      b0.x3f(m,k+1,j,i) = b3;
    }
  });

  return;
}
