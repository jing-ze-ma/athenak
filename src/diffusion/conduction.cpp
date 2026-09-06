//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.cpp
//! \brief Implements functions for Conduction class. This includes isotropic thermal
//! conduction, in which heat flux is proportional to negative local temperature gradient.
//! Conduction may be added to Hydro and/or MHD independently.

#include <float.h>
#include <algorithm>
#include <limits>
#include <string>
#include <iostream> // cout

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "eos/eos.hpp"
#include "conduction.hpp"
#include "utils/rosseland.hpp"
#include "coordinates/coordinates.hpp"
#include "units/units.hpp"

// VanLeer Limiter which takes 2 slopes
KOKKOS_INLINE_FUNCTION
Real VLL2State(const Real a, const Real b) {
  if (a*b > 0) {
    return 2.0*a*b/(a+b);
  } else {
    return 0.0;
  }
}

// VanLeer Limiter which takes 4 slopes
KOKKOS_INLINE_FUNCTION
Real VLL4State(const Real a, const Real b, const Real c, const Real d) {
  return VLL2State(VLL2State(a,b), VLL2State(c,d));
}

//----------------------------------------------------------------------------------------
//! \fn Real TempDepKappa()
//! \brief Temperature-dependent conductivity given by Parker (1953) and Spitzer (1962)

KOKKOS_INLINE_FUNCTION
Real TempDepKappa(Real temp, Real limit) {
  if (temp < 6.5e4) {
    return 2.5e3 * pow(temp, 0.5);
  } else {
    return fmin(6.0e-7*pow(temp, 2.5), limit);
  }
}

//----------------------------------------------------------------------------------------
//! \brief Conduction constructor
// Note first argument passes string ("hydro" or "mhd") denoting in wihch class this
// object is being constructed, and therefore which <block> in the input file from which
// the parameters are read.
// Note that the coefficient of thermal conduction, kappa, corresponds to conductivity,
// not diffusivity. This is different from the coefficient used in Athena++.

Conduction::Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin) :
    pmy_pack(pp), my_block(block) {
  // Read parameters for isotropic thermal conduction (if any)
  if (pin->DoesParameterExist(block,"isotropic_conduction")) {
    iso_cond_type = pin->GetString(block,"isotropic_conduction");
    // Check for valid type
    if ((iso_cond_type.compare("constant") != 0) &&
        (iso_cond_type.compare("spitzer") != 0) &&
        (iso_cond_type.compare("spitzer_limited") != 0) &&
        (iso_cond_type.compare("radiative") != 0)) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Invalid choice for isotropic thermal conduction type" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // constant conductivity
    if (iso_cond_type.compare("constant") == 0) {
      kappa_iso = pin->GetReal(block,"kappa_iso");
    }
    kappa_iso_limit = pin->GetOrAddReal(block,"kappa_iso_limit",
                      static_cast<Real>(std::numeric_limits<float>::max()));
    if (iso_cond_type.compare("radiative") == 0) {
      if (pp->punit == nullptr) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "radiative conduction needs a <units> block"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      rad_met = pin->GetOrAddReal(block,"rad_met",0.0);
      // pressure cut in bar; the flux through the wall in erg/cm^2/s
      rad_pcut = pin->GetOrAddReal(block,"rad_pcut_bar",0.0)*1.0e6
                 /pp->punit->pressure_cgs();
      rad_flux_inner = pin->GetOrAddReal(block,"rad_flux_inner",0.0)
                       /(pp->punit->pressure_cgs()*pp->punit->velocity_cgs());
      rad_kappa_fac = pin->GetOrAddReal(block,"rad_kappa_fac",1.0);
      rad_flux_limit = pin->GetOrAddBoolean(block,"rad_flux_limit",true);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief Conduction destructor

Conduction::~Conduction() {
}

//----------------------------------------------------------------------------------------
//! \fn void AddHeatFluxes()
//! \brief Wrapper function that adds heat fluxes for different types of thermal
//! conduction to face-centered fluxes of conserved variables

void Conduction::AddHeatFluxes(const DvceArray5D<Real> &w0, const EOS_Data &eos,
    DvceFaceFld5D<Real> &flx) {
  if (iso_cond_type.compare("constant") == 0) {
    AddIsotropicHeatFluxConstCond(w0, eos, flx);
  } else if ((iso_cond_type.compare("spitzer") == 0) ||
             (iso_cond_type.compare("spitzer_limited") == 0)) {
    AddIsotropicHeatFluxSpitzerCond(w0, eos, flx);
  } else if (iso_cond_type.compare("radiative") == 0) {
    AddIsotropicHeatFluxRadiative(w0, eos, flx);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiativeKappa
//! \brief the radiative conductivity 16 sigma T^3/(3 kappa_R rho) in cgs, T in K, p in
//! dyn/cm^2, rho in g/cm^3

KOKKOS_INLINE_FUNCTION
Real RadiativeKappa(const Real tk, const Real pcgs, const Real rhocgs, const Real met,
                    const Real kfac) {
  const Real sigma_sb = 5.670374419e-5;
  const Real kr = kfac*RosselandFreedman2014(tk, pcgs, met);
  return 16.0*sigma_sb*tk*tk*tk/(3.0*kr*rhocgs);
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsotropicHeatFluxRadiative()
//! \brief the diffusion approximation to radiative transport in the optically thick
//! layers: heat flux -kappa_rad grad T on every face at or below the pressure cut, with
//! a smooth saturation at the free-streaming limit sigma T^4, and the imposed internal
//! flux
//! through the inner x1 wall.  On a curvilinear grid the radial gradient uses the
//! centroid spacing; the tangential gradients are negligible in an atmosphere and use
//! the cell widths.

void Conduction::AddIsotropicHeatFluxRadiative(const DvceArray5D<Real> &w0,
                                               const EOS_Data &eos,
                                               DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  auto &mb_bcs = pmy_pack->pmb->mb_bcs;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  Real gm1 = eos.gamma-1.0;
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real eflx_unit = pres_unit*pmy_pack->punit->velocity_cgs();   // erg/cm^2/s
  const Real met = rad_met, pcut = rad_pcut, kfac = rad_kappa_fac, fin = rad_flux_inner;
  const bool limit = rad_flux_limit;
  const Real sigma_sb = 5.670374419e-5;

  // the heat flux across one face in CODE units, from the two adjacent cell states and
  // the centroid distance dl (code units); zero above the pressure cut
  auto face_flux = [=] (const Real tl, const Real tr, const Real pl, const Real pr,
                        const Real dl_, const Real dr_, const Real dl) {
    const Real pf = 0.5*(pl + pr);
    if (pf < pcut) return 0.0;
    const Real tk = 0.5*(tl + tr)*temp_unit;
    const Real kap = RadiativeKappa(tk, pf*pres_unit, 0.5*(dl_ + dr_)*dens_unit, met,
                                    kfac);
    Real f = -kap*(tr - tl)*temp_unit/(dl*len_unit);     // erg/cm^2/s, positive outward
    if (limit) {
      // saturate smoothly at the free-streaming flux sigma T^4: 0.3 % at F = 0.08 sigma
      // T^4,
      // where the diffusion approximation is still exact, and never above sigma T^4
      const Real ffree = sigma_sb*tk*tk*tk*tk;
      f /= sqrt(1.0 + SQR(f/ffree));
    }
    return f/eflx_unit;
  };

  auto &flx1 = flx.x1f;
  par_for("radcond1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // the imposed internal flux through the inner wall replaces the gradient there
    if (i == is && fin != 0.0 &&
        (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user ||
         mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect)) {
      flx1(m,IEN,k,j,i) += fin;
      return;
    }
    const Real tl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? (x1v_(m,i) - x1v_(m,i-1)) : size.d_view(m).dx1;
    flx1(m,IEN,k,j,i) += face_flux(tl, tr, pl, pr, w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i),
                                   dl);
  });
  if (!multi_d) return;

  auto &flx2 = flx.x2f;
  par_for("radcond2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx2_(m,k,j-1,i) + dx2_(m,k,j,i)) : size.d_view(m).dx2;
    flx2(m,IEN,k,j,i) += face_flux(tl, tr, pl, pr, w0(m,IDN,k,j-1,i), w0(m,IDN,k,j,i),
                                   dl);
  });
  if (!three_d) return;

  auto &flx3 = flx.x3f;
  par_for("radcond3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real tl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx3_(m,k-1,j,i) + dx3_(m,k,j,i)) : size.d_view(m).dx3;
    flx3(m,IEN,k,j,i) += face_flux(tl, tr, pl, pr, w0(m,IDN,k-1,j,i), w0(m,IDN,k,j,i),
                                   dl);
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddIsotropicHeatFluxConstCond()
//! \brief Adds isotropic heat flux computed using constant conductivity to face-centered
//! fluxes of conserved variables

void Conduction::AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w0,
    const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  Real gm1 = eos.gamma-1.0;
  // General EOS: temperature and pressure were evaluated once per cell in ConsToPrim.
  // Reading the cached values matters here -- T(d,e) is a root find, and the stencils
  // below touch up to ~30 neighbouring cells per face.
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  Real &kappa_ = kappa_iso;

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("conduct1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtempdx;
    if (gen) {
      dtempdx = (wtemp_(m,k,j,i) - wtemp_(m,k,j,i-1))
                / size.d_view(m).dx1;
    } else {
      dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                 - w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1))
                * gm1 / size.d_view(m).dx1;
    }
    flx1(m,IEN,k,j,i) -= kappa_ * dtempdx;
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("conduct2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtempdx;
    if (gen) {
      dtempdx = (wtemp_(m,k,j,i) - wtemp_(m,k,j-1,i))
                / size.d_view(m).dx2;
    } else {
      dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                 - w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i))
                * gm1 / size.d_view(m).dx2;
    }
    flx2(m,IEN,k,j,i) -= kappa_ * dtempdx;
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("conduct3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real dtempdx;
    if (gen) {
      dtempdx = (wtemp_(m,k,j,i) - wtemp_(m,k-1,j,i))
                / size.d_view(m).dx3;
    } else {
      dtempdx = (w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)
                 - w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i))
                * gm1 / size.d_view(m).dx3;
    }
    flx3(m,IEN,k,j,i) -= kappa_ * dtempdx;
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TempDependentHeatFlux()
//! \brief Adds heat flux to face-centered fluxes of conserved variables with
//! temperature-dependent conductivity

void Conduction::AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w0,
    const EOS_Data &eos, DvceFaceFld5D<Real> &flx) {
/*
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool &sat_hflux_ = sat_hflux;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  Real gm1 = eos.gamma-1.0;
  // General EOS: temperature and pressure were evaluated once per cell in ConsToPrim.
  // Reading the cached values matters here -- T(d,e) is a root find, and the stencils
  // below touch up to ~30 neighbouring cells per face.
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  Real kappaceil = kappa_ceiling;
  Real temp_unit = pmy_pack->punit->temperature_cgs();
  Real kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                    pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();

  // fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("conduct1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
    Real temp_r = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    Real pres_l = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
    Real pres_r = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx1 = (temp_r-temp_l)/size.d_view(m).dx1;
    Real hflx = kappaf*dtempdx1;
    // Saturation of thermal conduction by harmonic mean
    if (sat_hflux_) {
      Real dtempdx2 = 0.0, dtempdx3 = 0.0;
      if (multi_d) {
        temp_ll = (gen ? wtemp_(m,k,j-1,i-1)
                       : w0(m,IEN,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1)*gm1);
        temp_lr = (gen ? wtemp_(m,k,j+1,i-1)
                       : w0(m,IEN,k,j+1,i-1)/w0(m,IDN,k,j+1,i-1)*gm1);
        temp_rl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
        temp_rr = (gen ? wtemp_(m,k,j+1,i) : w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)*gm1);
        dtempdx2 = VanLeerLimiter4State(temp_rr-temp_r,temp_r-temp_rl,
                                        temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx2;
      }
      if (three_d) {
        temp_ll = (gen ? wtemp_(m,k-1,j,i-1)
                       : w0(m,IEN,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1)*gm1);
        temp_lr = (gen ? wtemp_(m,k+1,j,i-1)
                       : w0(m,IEN,k+1,j,i-1)/w0(m,IDN,k+1,j,i-1)*gm1);
        temp_rl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
        temp_rr = (gen ? wtemp_(m,k+1,j,i) : w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)*gm1);
        dtempdx3 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                              temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx3;
      }
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx1(m,IEN,k,j,i) -= hflx;
  });
  if (pmy_pack->pmesh->one_d) {return;}

  // fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("conduct2",DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = 0.0, temp_r = 0.0, pres_l = 0.0, pres_r = 0.0;
    Real temp_ll = 0.0, temp_lr = 0.0, temp_rl = 0.0, temp_rr = 0.0;
    temp_l = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
    temp_r = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    pres_l = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
    pres_r = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx2 = (temp_r-temp_l)/size.d_view(m).dx2;
    Real hflx = kappaf*dtempdx2;
    // Saturation of thermal conduction
    if (sat_hflux_) {
      Real dtempdx1 = 0.0, dtempdx3 = 0.0;
      temp_ll = (gen ? wtemp_(m,k,j-1,i-1) : w0(m,IEN,k,j-1,i-1)/w0(m,IDN,k,j-1,i-1)*gm1);
      temp_lr = (gen ? wtemp_(m,k,j-1,i+1) : w0(m,IEN,k,j-1,i+1)/w0(m,IDN,k,j-1,i+1)*gm1);
      temp_rl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
      temp_rr = (gen ? wtemp_(m,k,j,i+1) : w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)*gm1);
      dtempdx1 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx1;
      if (three_d) {
        temp_ll = (gen ? wtemp_(m,k-1,j-1,i)
                       : w0(m,IEN,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i)*gm1);
        temp_lr = (gen ? wtemp_(m,k+1,j-1,i)
                       : w0(m,IEN,k+1,j-1,i)/w0(m,IDN,k+1,j-1,i)*gm1);
        temp_rl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
        temp_rr = (gen ? wtemp_(m,k+1,j,i) : w0(m,IEN,k+1,j,i)/w0(m,IDN,k+1,j,i)*gm1);
        dtempdx3 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                              temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx3;
      }
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx2(m,IEN,k,j,i) -= hflx;
  });
  if (pmy_pack->pmesh->two_d) {return;}

  // fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("conduct3",DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // Add heat fluxes into fluxes of conserved variables: energy
    Real temp_l = 0.0, temp_r = 0.0, pres_l = 0.0, pres_r = 0.0;
    Real temp_ll = 0.0, temp_lr = 0.0, temp_rl = 0.0, temp_rr = 0.0;
    temp_l = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
    temp_r = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    pres_l = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
    pres_r = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    Real kappaf = 0.5*(TempDepKappa(temp_unit*temp_l,kappaceil)+
                  TempDepKappa(temp_unit*temp_r,kappaceil))/kappa_unit;
    Real dtempdx3 = (temp_r-temp_l)/size.d_view(m).dx3;
    Real hflx = kappaf*dtempdx3;
    // Saturation of thermal conduction
    if (sat_hflux_) {
      Real dtempdx1 = 0.0, dtempdx2 = 0.0;
      temp_ll = (gen ? wtemp_(m,k-1,j,i-1) : w0(m,IEN,k-1,j,i-1)/w0(m,IDN,k-1,j,i-1)*gm1);
      temp_lr = (gen ? wtemp_(m,k-1,j,i+1) : w0(m,IEN,k-1,j,i+1)/w0(m,IDN,k-1,j,i+1)*gm1);
      temp_rl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
      temp_rr = (gen ? wtemp_(m,k,j,i+1) : w0(m,IEN,k,j,i+1)/w0(m,IDN,k,j,i+1)*gm1);
      dtempdx1 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx1;
      temp_ll = (gen ? wtemp_(m,k-1,j-1,i) : w0(m,IEN,k-1,j-1,i)/w0(m,IDN,k-1,j-1,i)*gm1);
      temp_lr = (gen ? wtemp_(m,k-1,j+1,i) : w0(m,IEN,k-1,j+1,i)/w0(m,IDN,k-1,j+1,i)*gm1);
      temp_rl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
      temp_rr = (gen ? wtemp_(m,k,j+1,i) : w0(m,IEN,k,j+1,i)/w0(m,IDN,k,j+1,i)*gm1);
      dtempdx2 = VL4Limiter(temp_rr-temp_r,temp_r-temp_rl,
                            temp_lr-temp_l,temp_l-temp_ll)/size.d_view(m).dx2;
      Real tempgrad = sqrt(SQR(dtempdx1)+SQR(dtempdx2)+SQR(dtempdx3));
      Real pres_cs = 0.5*(pres_l*sqrt(temp_l)+pres_r*sqrt(temp_r));
      Real sat_fac = 1.0/(1.0+kappaf*tempgrad/(1.5*pres_cs));
      hflx *= sat_fac;
    }
    flx3(m,IEN,k,j,i) -= hflx;
  });

*/
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::NewTimeStep()
//! \brief Compute new time step for thermal conduction.

void Conduction::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
//  if (sat_hflux == true) {
//    dtnew = static_cast<Real>(std::numeric_limits<float>::max());
//    return;
//  }

  // set flag for Spitzer conductivity
  bool spitzer = false;
  if ((iso_cond_type.compare("spitzer") == 0) ||
      (iso_cond_type.compare("spitzer_limited") == 0)) {
    spitzer = true;
  }
  Real limit_ = kappa_iso_limit;
  Real temp_unit=0.0, kappa_unit=0.0, pres_unit=1.0, dens_unit=1.0;
  const bool radiative = (iso_cond_type.compare("radiative") == 0);
  const Real met = rad_met, pcut = rad_pcut, kfac = rad_kappa_fac;

  if (spitzer || radiative) {
    temp_unit = pmy_pack->punit->temperature_cgs();
    kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                 pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();
    pres_unit = pmy_pack->punit->pressure_cgs();
    dens_unit = pmy_pack->punit->density_cgs();
  }
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &dx1_ = pmy_pack->pcoord->dx1;

  // capture variables for kernel
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
  auto &w0_ = w0;
  auto &multi_d = pmy_pack->pmesh->multi_d;
  auto &three_d = pmy_pack->pmesh->three_d;
  auto &size = pmy_pack->pmb->mb_size;
  Real gm1 = eos_data.gamma-1.0;
  // General EOS: temperature and pressure were evaluated once per cell in ConsToPrim.
  // Reading the cached values matters here -- T(d,e) is a root find, and the stencils
  // below touch up to ~30 neighbouring cells per face.
  const bool gen = eos_data.IsGeneral();
  auto eos_ = eos_data;   // by-value copy, capturable in the device lambda
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  Real kappa0 = kappa_iso;

  // find smallest timestep for thermal conduction in each cell
  // Note loop over all cells needed even for constant conductivity
  Kokkos::parallel_reduce("cond_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt) {
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    Real kappa_ = kappa0;
    if (spitzer) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      kappa_ = TempDepKappa(temp*temp_unit, limit_)/kappa_unit;
    } else if (radiative) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      Real pres = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      if (pres < pcut) return;   // no flux above the cut: no constraint
      kappa_ = RadiativeKappa(temp*temp_unit, pres*pres_unit, w0_(m,IDN,k,j,i)*dens_unit,
                              met, kfac)/kappa_unit;
    }

    // the heat diffusion time is dx^2 rho c_v / kappa. For an ideal gas c_v = 1/(gamma-1)
    // which is what the rho/gm1 below amounts to; a general EOS has c_v(d,e), and it can
    // be an order of magnitude larger inside an ionization zone, so this is not a
    // cosmetic substitution -- it directly sets the conduction-limited timestep.
    Real rcv = w0_(m,IDN,k,j,i)/gm1;
    if (gen) {
      rcv = w0_(m,IDN,k,j,i)*eos_.SpecificHeatCv(w0_(m,IDN,k,j,i), w0_(m,IEN,k,j,i));
    }

    const Real d1 = (curv && radiative) ? dx1_(m,k,j,i) : size.d_view(m).dx1;
    min_dt = fmin(min_dt, SQR(d1)/kappa_*rcv);
    if (multi_d) {
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx2)/kappa_*rcv);
    }
    if (three_d) {
      min_dt = fmin(min_dt, SQR(size.d_view(m).dx3)/kappa_*rcv);
    }
  }, Kokkos::Min<Real>(dtnew));
  dtnew *= fac;

  return;
}
