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
#include "globals.hpp"
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
      // negative: the problem generator sets it (deep_hot_jupiter_rt: sigma T_int^4)
      rad_flux_inner = pin->GetOrAddReal(block,"rad_flux_inner",0.0);
      if (rad_flux_inner > 0.0) {
        rad_flux_inner /= (pp->punit->pressure_cgs()*pp->punit->velocity_cgs());
      }
      rad_kappa_fac = pin->GetOrAddReal(block,"rad_kappa_fac",1.0);
      rad_flux_limit = pin->GetOrAddBoolean(block,"rad_flux_limit",true);
      rad_tau_lo = pin->GetOrAddReal(block,"rad_tau_lo",0.0);
      rad_tau_hi = pin->GetOrAddReal(block,"rad_tau_hi",0.0);
      rad_tau_mode = (rad_tau_hi > 0.0);
      // ceiling on the temperature entering kappa_rad; 0 = off (bitwise inert)
      rad_tmax = pin->GetOrAddReal(block,"rad_tmax_kappa",0.0);
      if (rad_tmax > 0.0 && global_variable::my_rank == 0) {
        std::cout << "Conduction: radiative kappa temperature ceiling rad_tmax_kappa = "
                  << rad_tmax << " K" << std::endl;
      }
      rad_cs_exact = pin->GetOrAddBoolean(block,"rad_cs_exact",true);
      rad_blend_radial = pin->GetOrAddBoolean(block,"rad_blend_radial",true);
      {
        std::string ksrc = pin->GetOrAddString(block,"rad_kappa_src","freedman");
        if (ksrc.compare("table") == 0) {
          rad_kappa_tab = true;
        } else if (ksrc.compare("table_rho") == 0) {
          rad_kappa_tab = true;
          rad_kappa_rho = true;
        } else if (ksrc.compare("freedman") != 0) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_kappa_src must be freedman, table or table_rho"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }
      if (rad_tau_mode) {
        if (!(rad_tau_lo > 0.0 && rad_tau_lo < rad_tau_hi)) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "need 0 < rad_tau_lo < rad_tau_hi" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        auto &indcs = pp->pmesh->mb_indcs;
        const int nmb = pp->nmb_thispack;
        const int ncells1 = indcs.nx1 + 2*indcs.ng;
        const int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
        const int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
        Kokkos::realloc(rad_w, nmb, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(rad_tauf, nmb, ncells3, ncells2, ncells1+1);
      }
    }
  }
  Kokkos::realloc(dt_diag, ndtdiag);
}

//----------------------------------------------------------------------------------------
//! \brief Conduction destructor

Conduction::~Conduction() {
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::EnableDiag()
//! \brief Allocate the per-cycle diagnostic array.  NOT done in the constructor: the
//! diagnostic is a problem-generator option, and the pgen runs after AddPhysics has
//! built this object.  With AMR the pack size can change and this array would be stale,
//! but the runs this exists for (cubed-sphere deep hot Jupiter) have no AMR.

void Conduction::EnableDiag(int nmb, int n3, int n2, int n1) {
  diag = true;
  Kokkos::realloc(cond_diag, nmb, 6, n3, n2, n1);
  Kokkos::deep_copy(cond_diag, 0.0);
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
    if (rad_tau_mode) BuildRadWeights(w0, eos);
    AddIsotropicHeatFluxRadiative(w0, eos, flx);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::BuildRadWeights
//! \brief the Rosseland optical depth of every x1 face from the top of the domain,
//! tau_R = sum kappa_R rho dr over the cells above, and the blend weight w(tau_R).
//! One serial sweep per column, parallel over columns; the top face has tau = 0 (the
//! two-stream RT starts from the same point).

void Conduction::BuildRadWeights(const DvceArray5D<Real> &w0, const EOS_Data &eos) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, ng = indcs.ng;
  const int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &dx1_ = pmy_pack->pcoord->dx1;
  const Real gm1 = eos.gamma-1.0;
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wtemp
                                                : pmy_pack->phydro->wtemp;
  auto &wder_ = (my_block.compare("mhd") == 0) ? pmy_pack->pmhd->wder
                                               : pmy_pack->phydro->wder;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real met = rad_met, kfac = rad_kappa_fac, lo = rad_tau_lo, hi = rad_tau_hi;
  auto &wf = rad_w;
  auto &tf = rad_tauf;
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;   // table axis is log rho, not log p
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  par_for("radtau", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    Real tau = 0.0;
    tf(m,k,j,ie+1) = 0.0;
    wf(m,k,j,ie+1) = RadBlendWeight(0.0, lo, hi);
    for (int i=ie; i>=is-ng; --i) {
      const Real t = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      const Real p = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      const Real rho = w0(m,IDN,k,j,i)*dens_unit;
      const Real dr = (curv ? dx1_(m,k,j,i) : size.d_view(m).dx1)*len_unit;
      const Real kr = ktab
          ? RosselandTable(krt, krlT, krlP, krnT, krnP, t*temp_unit,
                           krho ? rho : p*pres_unit)
          : RosselandFreedman2014(t*temp_unit, p*pres_unit, met);
      tau += kfac*kr*rho*dr;
      tf(m,k,j,i) = tau;
      wf(m,k,j,i) = RadBlendWeight(tau, lo, hi);
    }
  });
  rad_w_built = true;
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

//! \brief the same from a supplied kappa_R [cm^2/g]
KOKKOS_INLINE_FUNCTION
Real RadiativeKappaKR(const Real tk, const Real rhocgs, const Real kfac, const Real kr) {
  const Real sigma_sb = 5.670374419e-5;
  return 16.0*sigma_sb*tk*tk*tk/(3.0*kfac*kr*rhocgs);
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
  const Real met = rad_met, kfac = rad_kappa_fac, fin = rad_flux_inner;
  const Real tmax = rad_tmax;   // temperature ceiling in kappa_rad only (0 = off)
  // pressure cut, or the tau blend: with the blend every face is masked by its weight
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  auto &wf = rad_w;
  const bool limit = rad_flux_limit;
  const Real sigma_sb = 5.670374419e-5;

  // the heat flux across one face in CODE units, from the two adjacent cell states and
  // the centroid distance dl (code units); zero above the pressure cut
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;   // table axis is log rho, not log p
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  // gradn is the FACE-NORMAL temperature derivative in code units (T per length); the
  // caller forms it, which is where the grid enters
  auto face_flux = [=] (const Real tl, const Real tr, const Real pl, const Real pr,
                        const Real dl_, const Real dr_, const Real gradn) {
    const Real pf = 0.5*(pl + pr);
    if (pf < pcut) return 0.0;
    const Real tk = 0.5*(tl + tr)*temp_unit;
    const Real rhof = 0.5*(dl_ + dr_)*dens_unit;
    // the CAPPED temperature enters kappa_rad only; the gradient below and the
    // free-streaming limit use the true face temperature tk
    const Real tkap = KappaTemp(tk, tmax);
    const Real kap = ktab
        ? RadiativeKappaKR(tkap, rhof, kfac,
                           RosselandTable(krt, krlT, krlP, krnT, krnP, tkap,
                                          krho ? rhof : pf*pres_unit))
        : RadiativeKappa(tkap, pf*pres_unit, rhof, met, kfac);
    Real f = -kap*gradn*temp_unit/len_unit;     // erg/cm^2/s, positive outward
    if (limit) {
      // saturate smoothly at the free-streaming flux sigma T^4: 0.3 % at F = 0.08 sigma
      // T^4,
      // where the diffusion approximation is still exact, and never above sigma T^4
      const Real ffree = sigma_sb*tk*tk*tk*tk;
      f /= sqrt(1.0 + SQR(f/ffree));
    }
    return f/eflx_unit;
  };

  // CUBED SPHERE: the xi and eta coordinate lines meet at an angle alpha (cos_cell,
  // sin_cell), so the face-normal derivative on a xi-face is
  //   dT/dn = (dT/dl_xi - cos(alpha) dT/dl_eta) / sin(alpha)
  // with dl the arc lengths -- the second term is the metric cross term, the 1/sin the
  // normalisation of grad(xi). The transverse derivative is centred across the two
  // cells the face separates. Radial faces are orthogonal to both and need nothing.
  const bool cs = pmy_pack->pmesh->use_cubed_sphere && rad_cs_exact;
  auto &sinc_ = pmy_pack->pcoord->sin_cell;
  auto &cosc_ = pmy_pack->pcoord->cos_cell;
  auto tcell = [=] (const int m, const int k, const int j, const int i) {
    return gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1;
  };

  // per-cycle diagnostic (problem/diag_gid): locals only, never `this`
  const bool diag_ = diag;
  auto cdg = cond_diag;

  auto &flx1 = flx.x1f;
  par_for("radcond1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // the imposed internal flux through the inner wall replaces the gradient there
    if (i == is && fin != 0.0 &&
        (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user ||
         mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect)) {
      flx1(m,IEN,k,j,i) += fin;
      if (diag_) cdg(m,0,k,j,i) = fin;
      return;
    }
    const Real tl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? (x1v_(m,i) - x1v_(m,i-1)) : size.d_view(m).dx1;
    const Real wt = (taumode && blend_r) ? wf(m,k,j,i) : 1.0;
    const Real fadd = wt*face_flux(tl, tr, pl, pr, w0(m,IDN,k,j,i-1),
                                   w0(m,IDN,k,j,i), (tr - tl)/dl);
    flx1(m,IEN,k,j,i) += fadd;
    if (diag_) cdg(m,0,k,j,i) = fadd;
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
    const Real wt = taumode ? 0.25*(wf(m,k,j-1,i) + wf(m,k,j-1,i+1)
                                    + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    Real gradn = (tr - tl)/dl;
    if (cs && three_d) {
      const Real c = 0.5*(cosc_(m,k,j-1) + cosc_(m,k,j));
      const Real sn = 0.5*(sinc_(m,k,j-1) + sinc_(m,k,j));
      const Real ge = 0.5*((tcell(m,k+1,j-1,i) - tcell(m,k-1,j-1,i))
                            /(0.5*dx3_(m,k-1,j-1,i) + dx3_(m,k,j-1,i)
                              + 0.5*dx3_(m,k+1,j-1,i))
                         + (tcell(m,k+1,j,i) - tcell(m,k-1,j,i))
                            /(0.5*dx3_(m,k-1,j,i) + dx3_(m,k,j,i) + 0.5*dx3_(m,k+1,j,i)));
      gradn = (gradn - c*ge)/sn;
    }
    const Real fadd = wt*face_flux(tl, tr, pl, pr, w0(m,IDN,k,j-1,i),
                                   w0(m,IDN,k,j,i), gradn);
    flx2(m,IEN,k,j,i) += fadd;
    if (diag_) cdg(m,1,k,j,i) = fadd;
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
    const Real wt = taumode ? 0.25*(wf(m,k-1,j,i) + wf(m,k-1,j,i+1)
                                    + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    Real gradn = (tr - tl)/dl;
    if (cs) {
      const Real c = 0.5*(cosc_(m,k-1,j) + cosc_(m,k,j));
      const Real sn = 0.5*(sinc_(m,k-1,j) + sinc_(m,k,j));
      const Real gx = 0.5*((tcell(m,k-1,j+1,i) - tcell(m,k-1,j-1,i))
                            /(0.5*dx2_(m,k-1,j-1,i) + dx2_(m,k-1,j,i)
                              + 0.5*dx2_(m,k-1,j+1,i))
                         + (tcell(m,k,j+1,i) - tcell(m,k,j-1,i))
                            /(0.5*dx2_(m,k,j-1,i) + dx2_(m,k,j,i) + 0.5*dx2_(m,k,j+1,i)));
      gradn = (gradn - c*gx)/sn;
    }
    const Real fadd = wt*face_flux(tl, tr, pl, pr, w0(m,IDN,k-1,j,i),
                                   w0(m,IDN,k,j,i), gradn);
    flx3(m,IEN,k,j,i) += fadd;
    if (diag_) cdg(m,2,k,j,i) = fadd;
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
  // at initialisation the blend weights do not exist yet; without them this timestep is
  // unconstrained and the first cycle runs the operator far past its stability limit
  if (rad_tau_mode && !rad_w_built) {
    BuildRadWeights(w0, eos_data);
  }
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
  Real vel_unit = 1.0;
  const bool radiative = (iso_cond_type.compare("radiative") == 0);
  const Real met = rad_met, kfac = rad_kappa_fac;
  const Real tmax = rad_tmax;   // temperature ceiling in kappa_rad only (0 = off)
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  const bool limit = rad_flux_limit;
  auto &wf = rad_w;
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;   // table axis is log rho, not log p
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;

  if (spitzer || radiative) {
    temp_unit = pmy_pack->punit->temperature_cgs();
    kappa_unit = pmy_pack->punit->pressure_cgs()*pmy_pack->punit->velocity_cgs()*
                 pmy_pack->punit->length_cgs()/pmy_pack->punit->temperature_cgs();
    pres_unit = pmy_pack->punit->pressure_cgs();
    dens_unit = pmy_pack->punit->density_cgs();
    vel_unit = pmy_pack->punit->velocity_cgs();
  }
  const bool curv = pmy_pack->pmesh->use_spherical_polar
                    || pmy_pack->pmesh->use_cubed_sphere;
  auto &dx1_ = pmy_pack->pcoord->dx1;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  const bool cs = pmy_pack->pmesh->use_cubed_sphere && rad_cs_exact;
  auto &sinc_ = pmy_pack->pcoord->sin_cell;

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
  // MinLOC, not Min: when this timestep collapses the only question that matters is
  // WHICH cell did it, and reconstructing that afterwards from a dump means redoing the
  // opacity lookup, the tau blend and the flux limiter outside the code.  The location
  // rides along for free.
  // per-cycle diagnostic (problem/diag_gid): locals only, never `this`
  const bool diag_ = diag;
  auto cdg = cond_diag;
  const Real dt_huge = static_cast<Real>(std::numeric_limits<float>::max());

  Kokkos::ValLocScalar<Real, int> mloc;
  Kokkos::parallel_reduce("cond_newdt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Kokkos::ValLocScalar<Real, int> &mres) {
    Real min_dt = mres.val;
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;

    // the diagnostic slots for cells that take one of the early returns below stay at
    // zero, which is what "this cell puts no constraint on the timestep" means
    if (diag_) {
      cdg(m,3,k,j,i) = 0.0;
      cdg(m,4,k,j,i) = 0.0;
      cdg(m,5,k,j,i) = 0.0;
    }

    Real kappa_ = kappa0;
    Real wmax = 1.0;
    if (spitzer) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      kappa_ = TempDepKappa(temp*temp_unit, limit_)/kappa_unit;
    } else if (radiative) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      Real pres = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      if (pres < pcut) return;   // no flux above the cut: no constraint
      const Real tkap = KappaTemp(temp*temp_unit, tmax);
      kappa_ = (ktab
          ? RadiativeKappaKR(tkap, w0_(m,IDN,k,j,i)*dens_unit, kfac,
                             RosselandTable(krt, krlT, krlP, krnT, krnP, tkap,
                                            krho ? w0_(m,IDN,k,j,i)*dens_unit
                                                 : pres*pres_unit))
          : RadiativeKappa(tkap, pres*pres_unit, w0_(m,IDN,k,j,i)*dens_unit,
                           met, kfac))/kappa_unit;
      // the blend: the face flux is w*kappa*grad T, so the explicit limit is on w*kappa,
      // and a cell whose faces carry no weight carries no constraint
      if (taumode) {
        wmax = fmax(wf(m,k,j,i), wf(m,k,j,i+1));
        if (wmax == 0.0 && blend_r) return;
      }
    }

    // the heat diffusion time is dx^2 rho c_v / kappa. For an ideal gas c_v = 1/(gamma-1)
    // which is what the rho/gm1 below amounts to; a general EOS has c_v(d,e), and it can
    // be an order of magnitude larger inside an ionization zone, so this is not a
    // cosmetic substitution -- it directly sets the conduction-limited timestep.
    Real rcv = w0_(m,IDN,k,j,i)/gm1;
    if (gen) {
      rcv = w0_(m,IDN,k,j,i)*eos_.SpecificHeatCv(w0_(m,IDN,k,j,i), w0_(m,IEN,k,j,i));
    }

    // THE FLUX-LIMITED DIFFUSIVITY.  With the limiter the face flux is
    //   F = kappa g / sqrt(1 + (kappa g / F_free)^2),   g = |grad T|,
    // and what an explicit step has to resolve is its linearisation
    //   dF/dg = kappa (1 + s^2)^(-3/2),   s = kappa g / F_free,
    // not kappa itself, which is unbounded where kappa_R rho -> 0: in the optically
    // thin atmosphere of a star kappa is ~1e8 x the limited value and the timestep
    // collapsed to microseconds.  Central differences of the cell temperature give g
    // per direction (the ghosts are filled when this runs).
    Real ffree = 0.0;
    if (radiative && limit) {
      const Real tk = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1)
                      *temp_unit;
      ffree = 5.670374419e-5*tk*tk*tk*tk/(pres_unit*vel_unit);   // code units
    }
    auto tc = [&] (const int kk, const int jj, const int ii) {
      return (gen ? wtemp_(m,kk,jj,ii) : w0(m,IEN,kk,jj,ii)/w0(m,IDN,kk,jj,ii)*gm1);
    };
    auto keff = [&] (const Real dl, const Real tm, const Real tp) {
      if (!(radiative && limit) || !(ffree > 0.0)) return kappa_;
      const Real s = kappa_*fabs(tp - tm)/(2.0*dl)/ffree;
      return kappa_/((1.0 + s*s)*sqrt(1.0 + s*s));
    };

    // the cell's OWN dt candidate, kept separate from the running reduction minimum
    Real dtc = dt_huge;
    const Real d1 = (curv && radiative) ? dx1_(m,k,j,i) : size.d_view(m).dx1;
    {
      const Real w1 = (taumode && blend_r) ? wmax : 1.0;
      const Real k1 = w1*keff(d1, tc(k,j,i-1), tc(k,j,i+1));
      if (k1 > 0.0) {
        min_dt = fmin(min_dt, SQR(d1)/k1*rcv);
        dtc = fmin(dtc, SQR(d1)/k1*rcv);
      }
      if (diag_) {
        cdg(m,3,k,j,i) = kappa_;
        cdg(m,4,k,j,i) = keff(d1, tc(k,j,i-1), tc(k,j,i+1));
      }
    }
    // on a curvilinear grid size.dx2/dx3 are ANGLES; the physical widths are pcoord's
    // cubed sphere: the exact operator's angular diffusivity is kappa/sin^2(alpha)
    const Real s2 = (cs && radiative) ? SQR(sinc_(m,k,j)) : 1.0;
    const Real wa = taumode ? wmax : 1.0;
    if (multi_d) {
      const Real d2 = (curv && radiative) ? dx2_(m,k,j,i) : size.d_view(m).dx2;
      const Real k2 = wa*keff(d2, tc(k,j-1,i), tc(k,j+1,i));
      if (k2 > 0.0) {
        min_dt = fmin(min_dt, SQR(d2)*s2/k2*rcv);
        dtc = fmin(dtc, SQR(d2)*s2/k2*rcv);
      }
    }
    if (three_d) {
      const Real d3 = (curv && radiative) ? dx3_(m,k,j,i) : size.d_view(m).dx3;
      const Real k3 = wa*keff(d3, tc(k-1,j,i), tc(k+1,j,i));
      if (k3 > 0.0) {
        min_dt = fmin(min_dt, SQR(d3)*s2/k3*rcv);
        dtc = fmin(dtc, SQR(d3)*s2/k3*rcv);
      }
    }
    if (diag_) cdg(m,5,k,j,i) = dtc;
    if (min_dt < mres.val) { mres.val = min_dt; mres.loc = idx; }
  }, Kokkos::MinLoc<Real, int>(mloc));
  dtnew = mloc.val*fac;
  // decode the winning cell, for the collapse report in Mesh::NewTimeStep
  if (mloc.loc >= 0 && mloc.loc < nmkji) {
    dtnew_m = (mloc.loc)/nkji;
    dtnew_k = (mloc.loc - dtnew_m*nkji)/nji + ks;
    dtnew_j = (mloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji)/nx1 + js;
    dtnew_i = (mloc.loc - dtnew_m*nkji - (dtnew_k-ks)*nji - (dtnew_j-js)*nx1) + is;
  } else {
    dtnew_m = dtnew_k = dtnew_j = dtnew_i = -1;
  }

  // THE STATE OF THAT CELL.  Where the cell is does not say why it is slow: the blend
  // weight, the flux limiter's s and the heat capacity all set the number and none of
  // them appears in any output.  Recompute them for the one winning cell, and only when
  // dtnew has just collapsed (or on the first call), so a healthy run pays nothing.
  dt_diag_valid = false;
  if (radiative && dtnew_m >= 0 && (dtnew_prev < 0.0 || dtnew < 0.25*dtnew_prev)) {
    auto dd = dt_diag;
    const int dm = dtnew_m, dk = dtnew_k, dj = dtnew_j, di = dtnew_i;
    auto &x1v_ = pmy_pack->pcoord->x1v;
    auto &tfd = rad_tauf;
    par_for("cond_dtdiag", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
      const Real temp = (gen ? wtemp_(dm,dk,dj,di)
                             : w0_(dm,IEN,dk,dj,di)/w0_(dm,IDN,dk,dj,di)*gm1);
      const Real pres = (gen ? wder_(dm,IDPR,dk,dj,di) : w0_(dm,IEN,dk,dj,di)*gm1);
      const Real dens = w0_(dm,IDN,dk,dj,di);
      const Real tkap = KappaTemp(temp*temp_unit, tmax);
      const Real kr = (ktab ? RosselandTable(krt, krlT, krlP, krnT, krnP, tkap,
                                             krho ? dens*dens_unit : pres*pres_unit)
                            : -1.0);
      const Real kappa_ = (ktab
          ? RadiativeKappaKR(tkap, dens*dens_unit, kfac, kr)
          : RadiativeKappa(tkap, pres*pres_unit, dens*dens_unit, met, kfac))
          /kappa_unit;
      Real rcv = dens/gm1;
      if (gen) {
        rcv = dens*eos_.SpecificHeatCv(dens, w0_(dm,IEN,dk,dj,di));
      }
      Real ffree = 0.0;
      if (limit) {
        const Real tk = temp*temp_unit;
        ffree = 5.670374419e-5*tk*tk*tk*tk/(pres_unit*vel_unit);
      }
      auto tc = [&] (const int kk, const int jj, const int ii) {
        return (gen ? wtemp_(dm,kk,jj,ii)
                    : w0_(dm,IEN,kk,jj,ii)/w0_(dm,IDN,kk,jj,ii)*gm1);
      };
      auto sof = [&] (const Real dl, const Real tm, const Real tp) {
        return (limit && ffree > 0.0) ? kappa_*fabs(tp - tm)/(2.0*dl)/ffree : 0.0;
      };
      auto keff = [&] (const Real s) {
        return (limit && ffree > 0.0) ? kappa_/((1.0 + s*s)*sqrt(1.0 + s*s)) : kappa_;
      };
      const Real wmax = taumode ? fmax(wf(dm,dk,dj,di), wf(dm,dk,dj,di+1)) : 1.0;
      const Real d1 = curv ? dx1_(dm,dk,dj,di) : size.d_view(dm).dx1;
      const Real d2 = curv ? dx2_(dm,dk,dj,di) : size.d_view(dm).dx2;
      const Real d3 = curv ? dx3_(dm,dk,dj,di) : size.d_view(dm).dx3;
      const Real s2 = cs ? SQR(sinc_(dm,dk,dj)) : 1.0;
      const Real s1v = sof(d1, tc(dk,dj,di-1), tc(dk,dj,di+1));
      const Real s2v = sof(d2, tc(dk,dj-1,di), tc(dk,dj+1,di));
      const Real s3v = sof(d3, tc(dk-1,dj,di), tc(dk+1,dj,di));
      const Real w1 = (taumode && blend_r) ? wmax : 1.0;
      const Real wa = taumode ? wmax : 1.0;
      dd.d_view(0)  = x1v_(dm,di);
      dd.d_view(1)  = dens*dens_unit;
      dd.d_view(2)  = temp*temp_unit;
      dd.d_view(3)  = pres*pres_unit;
      dd.d_view(4)  = kr;
      dd.d_view(5)  = kappa_;
      dd.d_view(6)  = rcv;
      dd.d_view(7)  = taumode ? wf(dm,dk,dj,di)   : 1.0;
      dd.d_view(8)  = taumode ? wf(dm,dk,dj,di+1) : 1.0;
      dd.d_view(9)  = taumode ? tfd(dm,dk,dj,di)  : -1.0;
      dd.d_view(10) = taumode ? tfd(dm,dk,dj,di+1): -1.0;
      dd.d_view(11) = s1v;
      dd.d_view(12) = (w1*keff(s1v) > 0.0) ? SQR(d1)/(w1*keff(s1v))*rcv*fac : -1.0;
      dd.d_view(13) = (multi_d && wa*keff(s2v) > 0.0)
                      ? SQR(d2)*s2/(wa*keff(s2v))*rcv*fac : -1.0;
      dd.d_view(14) = (three_d && wa*keff(s3v) > 0.0)
                      ? SQR(d3)*s2/(wa*keff(s3v))*rcv*fac : -1.0;
      dd.d_view(15) = ffree;
    });
    dt_diag.modify_device();
    dt_diag.sync_host();
    dt_diag_valid = true;
  }
  dtnew_prev = dtnew;

  return;
}
