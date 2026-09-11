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
#include "utils/eint_from_cons.hpp"

// <problem>/nan_report: the first x1 radiative-conduction face whose flux comes out
// non-finite, and the inputs that made it.  Pointers, not Views: a file-scope View
// outlives Kokkos::finalize.  Nothing is allocated or run unless the switch is on.
namespace {
DvceArray1D<int> *cnd_nanrep_cnt = nullptr;
DvceArray1D<Real> *cnd_nanrep_rec = nullptr;
int cnd_nanrep_lines = 0;
const int cnd_nanrep_maxlines = 400;
}  // namespace
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
      nan_report = pin->GetOrAddBoolean("problem","nan_report",false);
      rad_flux_limit = pin->GetOrAddBoolean(block,"rad_flux_limit",true);
      rad_tau_lo = pin->GetOrAddReal(block,"rad_tau_lo",0.0);
      rad_tau_hi = pin->GetOrAddReal(block,"rad_tau_hi",0.0);
      rad_tau_mode = (rad_tau_hi > 0.0);
      rad_cs_exact = pin->GetOrAddBoolean(block,"rad_cs_exact",true);
      rad_kappa_rmax = pin->GetOrAddReal(block,"rad_kappa_rmax",0.0);
      rad_kappa_above = pin->GetOrAddReal(block,"rad_kappa_above",0.0);
      rad_gate_rho = pin->GetOrAddReal(block,"rad_gate_rho",0.0);
      rad_gate_dex = pin->GetOrAddReal(block,"rad_gate_dex",0.5);
      if (rad_gate_rho > 0.0 && rad_kappa_rmax > 0.0) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_gate_rho and rad_kappa_rmax are two forms of the "
                  << "same switch: set one or the other, not both" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (rad_gate_rho > 0.0 && !(rad_gate_dex > 0.0)) {
        std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                  << std::endl << "rad_gate_dex must be positive" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      rad_blend_radial = pin->GetOrAddBoolean(block,"rad_blend_radial",true);
      rad_implicit_x1 = pin->GetOrAddBoolean(block,"rad_implicit_x1",false);
      rad_cap_ang = pin->GetOrAddReal(block,"rad_cap_ang",0.0);
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
      if (rad_implicit_x1) {
        // hydro only: the MHD inversion would have to subtract the magnetic energy too
        if (block.compare("hydro") != 0) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_implicit_x1 is implemented for <hydro> only"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        // the solve is column-local, so the whole radial extent must be in one MeshBlock
        if (pp->pmesh->mb_indcs.nx1 != pp->pmesh->mesh_indcs.nx1) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_implicit_x1 needs the whole radial extent in "
                    << "one MeshBlock: <meshblock>/nx1 = " << pp->pmesh->mb_indcs.nx1
                    << " but <mesh>/nx1 = " << pp->pmesh->mesh_indcs.nx1 << std::endl;
          std::exit(EXIT_FAILURE);
        }
        auto &indcs = pp->pmesh->mb_indcs;
        const int nmb = pp->nmb_thispack;
        const int ncells1 = indcs.nx1 + 2*indcs.ng;
        const int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
        const int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
        Kokkos::realloc(imp_wrk, nmb, nimpw, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(imp_flag, 1);
        Kokkos::realloc(imp_rec, 8);
      }
      if (rad_cap_ang > 0.0) {
        // the angular cap only makes sense once the radial direction is unconditionally
        // stable: with x1 still explicit the run dies on dt1 long before dt2/dt3 matter,
        // and dropping dt2/dt3 while dt1 still applies would be a silent half-measure
        if (!rad_implicit_x1) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_cap_ang > 0 requires rad_implicit_x1 = true"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }
        auto &indcs = pp->pmesh->mb_indcs;
        if (indcs.ng < 2) {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "rad_cap_ang needs nghost >= 2: the stiffness of a "
                    << "ghost cell reads its own outer face" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        const int nmb = pp->nmb_thispack;
        const int ncells1 = indcs.nx1 + 2*indcs.ng;
        const int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*indcs.ng) : 1;
        const int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*indcs.ng) : 1;
        Kokkos::realloc(cap_x, nmb, ncells3, ncells2, ncells1);
        Kokkos::realloc(cap_c2, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(cap_c3, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::realloc(cap_cnt, 2);
        Kokkos::realloc(cap_rec, 6);
      }
    } else if (pin->GetOrAddBoolean(block,"rad_implicit_x1",false)) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                << std::endl << "rad_implicit_x1 needs isotropic_conduction = radiative"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  Kokkos::realloc(dt_diag, ndtdiag);
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
  // the radiatively inert region above rad_kappa_rmax: tau accumulates at rad_kappa_above
  // there, so it does not grow through the corona and the blend weight stays 0
  const Real krmax = rad_kappa_rmax, kabove = rad_kappa_above;
  // the DENSITY gate (rad_gate_rho): the same inert medium, selected by what the gas is
  // rather than by where it is.  kappa_eff = G kappa_table + (1 - G) rad_kappa_above.
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  auto &x1v_t = pmy_pack->pcoord->x1v;
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
      Real kr = (krmax > 0.0 && x1v_t(m,i) > krmax)
          ? kabove
          : (ktab
             ? RosselandTable(krt, krlT, krlP, krnT, krnP, t*temp_unit,
                              krho ? rho : p*pres_unit)
             : RosselandFreedman2014(t*temp_unit, p*pres_unit, met));
      if (gaterho > 0.0) {
        const Real g = RadGate(rho, gaterho, gatedex);
        kr = g*kr + (1.0 - g)*kabove;
      }
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
//! \fn Real RadFaceKappa
//! \brief the radiative conductivity of a face in cgs from its temperature, pressure and
//! density, taking kappa_R either from the tabulated lookup (on a log p or a log rho
//! axis) or from the Freedman fit.  This is the ONE place the face conductivity is
//! defined: the explicit face flux (AddIsotropicHeatFluxRadiative) and the implicit
//! radial solve (ImplicitRadialUpdate) both call it, so the two operators cannot drift
//! apart.

KOKKOS_INLINE_FUNCTION
Real RadFaceKappa(const Real tk, const Real pcgs, const Real rhocgs, const bool ktab,
                  const DvceArray2D<Real> &krt, const DvceArray1D<Real> &krlT,
                  const DvceArray1D<Real> &krlP, const int krnT, const int krnP,
                  const bool krho, const Real met, const Real kfac) {
  return ktab
      ? RadiativeKappaKR(tk, rhocgs, kfac,
                         RosselandTable(krt, krlT, krlP, krnT, krnP, tk,
                                        krho ? rhocgs : pcgs))
      : RadiativeKappa(tk, pcgs, rhocgs, met, kfac);
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
  // no radiative diffusion above rad_kappa_rmax: the corona is transparent by fiat, and
  // an opacity of zero would divide by zero in RadiativeKappaKR, so the face is skipped
  // outright rather than handed kappa = 0
  const Real krmax = rad_kappa_rmax;
  // the DENSITY gate: kappa follows the gas, not the radius.  For a FACE the gate reads
  // the same face-averaged density the flux does, so the two operators switch off over
  // exactly the same faces.  Folded into the tau-blend weight, which is where every
  // conduction face flux is already scaled.
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  // gradn is the FACE-NORMAL temperature derivative in code units (T per length); the
  // caller forms it, which is where the grid enters
  auto face_flux = [=] (const Real tl, const Real tr, const Real pl, const Real pr,
                        const Real dl_, const Real dr_, const Real gradn) {
    const Real pf = 0.5*(pl + pr);
    if (pf < pcut) return 0.0;
    const Real tk = 0.5*(tl + tr)*temp_unit;
    const Real rhof = 0.5*(dl_ + dr_)*dens_unit;
    const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                  krnT, krnP, krho, met, kfac);
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

  // THE SAME COEFFICIENT AS A CONDUCTIVITY.  face_flux returns -K gradn; the angular
  // cap (rad_cap_ang) needs K itself, in code units, to form the explicit stiffness of a
  // cell.  Identical expressions in the same order, from the same RadFaceKappa, so the
  // capped operator and the flux it caps cannot disagree about what K is.
  auto face_kcode = [=] (const Real tl, const Real tr, const Real pl, const Real pr,
                         const Real dl_, const Real dr_, const Real gradn) {
    const Real pf = 0.5*(pl + pr);
    if (pf < pcut) return 0.0;
    const Real tk = 0.5*(tl + tr)*temp_unit;
    const Real rhof = 0.5*(dl_ + dr_)*dens_unit;
    const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                  krnT, krnP, krho, met, kfac);
    Real lf = 1.0;
    if (limit) {
      const Real f = -kap*gradn*temp_unit/len_unit;
      const Real ffree = sigma_sb*tk*tk*tk*tk;
      lf = (ffree > 0.0) ? 1.0/sqrt(1.0 + SQR(f/ffree)) : 0.0;
    }
    return kap*lf*temp_unit/len_unit/eflx_unit;
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

  // <problem>/nan_report device counter/record; untouched when the switch is off
  const bool nanrep_c = nan_report;
  if (nanrep_c && cnd_nanrep_cnt == nullptr) {
    cnd_nanrep_cnt = new DvceArray1D<int>("cnd_nanrep_cnt", 1);
    cnd_nanrep_rec = new DvceArray1D<Real>("cnd_nanrep_rec", 16);
  }
  auto cndcnt = nanrep_c ? *cnd_nanrep_cnt : DvceArray1D<int>("d", 1);
  auto cndrec = nanrep_c ? *cnd_nanrep_rec : DvceArray1D<Real>("d", 1);
  if (nanrep_c) {
    Kokkos::deep_copy(cndcnt, 0);
    Kokkos::deep_copy(cndrec, 0.0);
  }

  auto &flx1 = flx.x1f;
  // rad_implicit_x1: the INTERIOR x1 faces are handled by ImplicitRadialUpdate after the
  // RK update, so nothing is added to them here.  The two boundary faces stay explicit:
  // the inner one is the imposed wall flux (or the ghost-based gradient), the outer one
  // the ghost-based gradient, and neither is part of the tridiagonal system.
  const bool impx1 = rad_implicit_x1;
  par_for("radcond1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    // the imposed internal flux through the inner wall replaces the gradient there
    if (i == is && fin != 0.0 &&
        (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user ||
         mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::reflect)) {
      flx1(m,IEN,k,j,i) += fin;
      return;
    }
    if (impx1 && i > is && i < ie+1) return;
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = (gen ? wtemp_(m,k,j,i-1) : w0(m,IEN,k,j,i-1)/w0(m,IDN,k,j,i-1)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j,i-1) : w0(m,IEN,k,j,i-1)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? (x1v_(m,i) - x1v_(m,i-1)) : size.d_view(m).dx1;
    Real wt = (taumode && blend_r) ? wf(m,k,j,i) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k,j,i-1) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
    const Real fcnd = wt*face_flux(tl, tr, pl, pr, w0(m,IDN,k,j,i-1),
                                   w0(m,IDN,k,j,i), (tr - tl)/dl);
    flx1(m,IEN,k,j,i) += fcnd;
    // --- <problem>/nan_report: record the first face whose conduction flux, or the
    // total flux it lands in, is not finite, with the inputs that produced it.
    if (nanrep_c) {
      if (!isfinite(fcnd) || !isfinite(flx1(m,IEN,k,j,i))) {
        if (Kokkos::atomic_fetch_add(&cndcnt(0), 1) == 0) {
          const Real tk = 0.5*(tl + tr)*temp_unit;
          const Real rhof = 0.5*(w0(m,IDN,k,j,i-1) + w0(m,IDN,k,j,i))*dens_unit;
          const Real pf = 0.5*(pl + pr);
          const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                        krnT, krnP, krho, met, kfac);
          cndrec(0) = static_cast<Real>(m);
          cndrec(1) = static_cast<Real>(k);
          cndrec(2) = static_cast<Real>(j);
          cndrec(3) = static_cast<Real>(i);
          cndrec(4) = tl*temp_unit;
          cndrec(5) = tr*temp_unit;
          cndrec(6) = kap;
          cndrec(7) = wt;
          cndrec(8) = (tr - tl)/dl*temp_unit/len_unit;
          cndrec(9) = dl*len_unit;
          cndrec(10) = fcnd;
          cndrec(11) = flx1(m,IEN,k,j,i);
          cndrec(12) = w0(m,IDN,k,j,i-1);
          cndrec(13) = w0(m,IDN,k,j,i);
          cndrec(14) = pl;
          cndrec(15) = pr;
        }
      }
    }
  });
  if (nanrep_c && cnd_nanrep_lines < cnd_nanrep_maxlines) {
    auto hcc = Kokkos::create_mirror_view(cndcnt);
    Kokkos::deep_copy(hcc, cndcnt);
    if (hcc(0) > 0) {
      auto hcr = Kokkos::create_mirror_view(cndrec);
      Kokkos::deep_copy(hcr, cndrec);
      ++cnd_nanrep_lines;
      const int mb = static_cast<int>(hcr(0));
      std::cout << "### nan_report [conduction_x1_face] rank "
                << global_variable::my_rank << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time << ": " << hcc(0)
                << " bad face(s); first (m,k,j,i) = (" << mb << ","
                << static_cast<int>(hcr(1)) << "," << static_cast<int>(hcr(2)) << ","
                << static_cast<int>(hcr(3)) << ") gid = " << (pmy_pack->gids + mb)
                << " T_l = " << hcr(4) << " T_r = " << hcr(5) << " kappa = " << hcr(6)
                << " w_tau = " << hcr(7) << " dT/dx = " << hcr(8)
                << " dl = " << hcr(9) << " F_cond = " << hcr(10)
                << " flx1(IEN) = " << hcr(11) << " d_l = " << hcr(12)
                << " d_r = " << hcr(13) << " p_l = " << hcr(14)
                << " p_r = " << hcr(15) << std::endl;
    }
  }
  if (!multi_d) return;

  // ------------------------------------------------------------------------------------
  // rad_cap_ang: the CONSERVATIVE per-face cap on the explicit angular operator.  Three
  // extra kernels evaluate the frozen face coefficient A_f K_f/(dl_f sin alpha) on the
  // x2 and x3 faces and then the per-cell stiffness
  //     x_i = beta_dt alpha_i (sum over its 4 angular faces of A_f K_f/dl_f)/V_i,
  // one cell into the angular ghosts, because the first interior face needs the
  // stiffness of the cell on its other side.  The flux kernels below then multiply each
  // face by min(1, cap/max(x_i,x_j)) -- one number per face, so nothing is created or
  // destroyed, only moved more slowly than an explicit step could resolve.
  const Real capang = rad_cap_ang;
  const Real capbdt = stage_beta_dt;
  auto capx = cap_x;
  auto capc2 = cap_c2;
  auto capc3 = cap_c3;
  auto capcnt = cap_cnt;
  auto caprec = cap_rec;
  if (capang > 0.0) {
    auto &vol_ = pmy_pack->pcoord->volume;
    auto &area2_ = pmy_pack->pcoord->area.x2f;
    auto &area3_ = pmy_pack->pcoord->area.x3f;
    auto eos_ = eos;
    Kokkos::deep_copy(capcnt, 0);

    par_for("radcapc2", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+2, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      capc2(m,k,j,i) = 0.0;
      if (krmax > 0.0 && x1v_(m,i) > krmax) return;
      const Real tl = tcell(m,k,j-1,i), tr = tcell(m,k,j,i);
      const Real pl = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
      const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      const Real dl = curv ? 0.5*(dx2_(m,k,j-1,i) + dx2_(m,k,j,i)) : size.d_view(m).dx2;
      Real wt = taumode ? 0.25*(wf(m,k,j-1,i) + wf(m,k,j-1,i+1)
                                + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
      if (gaterho > 0.0) {
        wt *= RadGate(0.5*(w0(m,IDN,k,j-1,i) + w0(m,IDN,k,j,i))*dens_unit,
                      gaterho, gatedex);
      }
      Real gradn = (tr - tl)/dl;
      Real sn = 1.0;
      if (cs && three_d) {
        const Real c = 0.5*(cosc_(m,k,j-1) + cosc_(m,k,j));
        sn = 0.5*(sinc_(m,k,j-1) + sinc_(m,k,j));
        const Real ge = 0.5*((tcell(m,k+1,j-1,i) - tcell(m,k-1,j-1,i))
                              /(0.5*dx3_(m,k-1,j-1,i) + dx3_(m,k,j-1,i)
                                + 0.5*dx3_(m,k+1,j-1,i))
                           + (tcell(m,k+1,j,i) - tcell(m,k-1,j,i))
                              /(0.5*dx3_(m,k-1,j,i) + dx3_(m,k,j,i)
                                + 0.5*dx3_(m,k+1,j,i)));
        gradn = (gradn - c*ge)/sn;
      }
      // the coefficient of (T_j - T_i)/dl in the face flux is K/sin(alpha): the metric
      // normalisation of grad(xi) is part of the stiffness, exactly as it is part of the
      // dt2 the cap replaces (SQR(d2)*s2/k2 in NewTimeStep)
      const Real kc = wt*face_kcode(tl, tr, pl, pr, w0(m,IDN,k,j-1,i),
                                    w0(m,IDN,k,j,i), gradn);
      const Real af = curv ? area2_(m,k,j,i) : 1.0/size.d_view(m).dx2;
      capc2(m,k,j,i) = kc*af/(dl*sn);
    });

    if (three_d) {
      par_for("radcapc3", DevExeSpace(), 0, nmb1, ks-1, ke+2, js-1, je+1, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        capc3(m,k,j,i) = 0.0;
        if (krmax > 0.0 && x1v_(m,i) > krmax) return;
        const Real tl = tcell(m,k-1,j,i), tr = tcell(m,k,j,i);
        const Real pl = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
        const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
        const Real dl = curv ? 0.5*(dx3_(m,k-1,j,i) + dx3_(m,k,j,i))
                             : size.d_view(m).dx3;
        Real wt = taumode ? 0.25*(wf(m,k-1,j,i) + wf(m,k-1,j,i+1)
                                  + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
        if (gaterho > 0.0) {
          wt *= RadGate(0.5*(w0(m,IDN,k-1,j,i) + w0(m,IDN,k,j,i))*dens_unit,
                        gaterho, gatedex);
        }
        Real gradn = (tr - tl)/dl;
        Real sn = 1.0;
        if (cs) {
          const Real c = 0.5*(cosc_(m,k-1,j) + cosc_(m,k,j));
          sn = 0.5*(sinc_(m,k-1,j) + sinc_(m,k,j));
          const Real gx = 0.5*((tcell(m,k-1,j+1,i) - tcell(m,k-1,j-1,i))
                                /(0.5*dx2_(m,k-1,j-1,i) + dx2_(m,k-1,j,i)
                                  + 0.5*dx2_(m,k-1,j+1,i))
                             + (tcell(m,k,j+1,i) - tcell(m,k,j-1,i))
                                /(0.5*dx2_(m,k,j-1,i) + dx2_(m,k,j,i)
                                  + 0.5*dx2_(m,k,j+1,i)));
          gradn = (gradn - c*gx)/sn;
        }
        const Real kc = wt*face_kcode(tl, tr, pl, pr, w0(m,IDN,k-1,j,i),
                                      w0(m,IDN,k,j,i), gradn);
        const Real af = curv ? area3_(m,k,j,i) : 1.0/size.d_view(m).dx3;
        capc3(m,k,j,i) = kc*af/(dl*sn);
      });
    }

    par_for("radcapx", DevExeSpace(), 0, nmb1, ks-1, ke+1, js-1, je+1, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = w0(m,IDN,k,j,i);
      Real rcv = d/gm1;
      if (gen) rcv = d*eos_.SpecificHeatCv(d, w0(m,IEN,k,j,i), wtemp_(m,k,j,i));
      Real sum = capc2(m,k,j,i) + capc2(m,k,j+1,i);
      if (three_d) sum += capc3(m,k,j,i) + capc3(m,k+1,j,i);
      const Real vi = curv ? vol_(m,k,j,i) : 1.0;
      capx(m,k,j,i) = (rcv > 0.0 && vi > 0.0 && sum > 0.0) ? capbdt*sum/(rcv*vi) : 0.0;
    });

    // the diagnostic: how stiff the angular operator actually is, and where.  Active
    // cells only -- a ghost is somebody else's cell and would be double counted.
    {
      const int nx1_ = indcs.nx1, nx2_ = indcs.nx2, nx3_ = indcs.nx3;
      const int nkji_ = nx3_*nx2_*nx1_, nji_ = nx2_*nx1_;
      Kokkos::ValLocScalar<Real, int> xloc;
      int nover = 0;
      Kokkos::parallel_reduce("radcapdiag",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkji_),
      KOKKOS_LAMBDA(const int &idx, Kokkos::ValLocScalar<Real, int> &xres, int &nov) {
        const int m = idx/nkji_;
        const int k = (idx - m*nkji_)/nji_ + ks;
        const int j = (idx - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
        const int i = (idx - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
        const Real xv = capx(m,k,j,i);
        if (xv > capang) ++nov;
        if (xv > xres.val) { xres.val = xv; xres.loc = idx; }
      }, Kokkos::MaxLoc<Real, int>(xloc), nover);
      const int xl = xloc.loc;
      const Real xv = xloc.val;
      par_for("radcaploc", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
        caprec(0) = xv;
        if (xl >= 0) {
          const int m = xl/nkji_;
          const int k = (xl - m*nkji_)/nji_ + ks;
          const int j = (xl - m*nkji_ - (k - ks)*nji_)/nx1_ + js;
          const int i = (xl - m*nkji_ - (k - ks)*nji_ - (j - js)*nx1_) + is;
          caprec(1) = static_cast<Real>(m);
          caprec(2) = static_cast<Real>(k);
          caprec(3) = static_cast<Real>(j);
          caprec(4) = static_cast<Real>(i);
          caprec(5) = x1v_(m,i);
        } else {
          caprec(1) = -1.0; caprec(2) = -1.0; caprec(3) = -1.0;
          caprec(4) = -1.0; caprec(5) = -1.0;
        }
      });
      cap_diag_x = xv;
      cap_diag_over = nover;
    }
  }

  auto &flx2 = flx.x2f;
  par_for("radcond2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = (gen ? wtemp_(m,k,j-1,i) : w0(m,IEN,k,j-1,i)/w0(m,IDN,k,j-1,i)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k,j-1,i) : w0(m,IEN,k,j-1,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx2_(m,k,j-1,i) + dx2_(m,k,j,i)) : size.d_view(m).dx2;
    Real wt = taumode ? 0.25*(wf(m,k,j-1,i) + wf(m,k,j-1,i+1)
                              + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k,j-1,i) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
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
    // rad_cap_ang: one scale factor per face, applied to both of its cells
    Real wtc = wt;
    if (capang > 0.0) {
      const Real xm = fmax(capx(m,k,j-1,i), capx(m,k,j,i));
      if (xm > capang) {
        wtc = wt*(capang/xm);
        Kokkos::atomic_fetch_add(&capcnt(1), 1);
      }
    }
    flx2(m,IEN,k,j,i) += wtc*face_flux(tl, tr, pl, pr, w0(m,IDN,k,j-1,i),
                                       w0(m,IDN,k,j,i), gradn);
  });
  if (!three_d) return;

  auto &flx3 = flx.x3f;
  par_for("radcond3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (krmax > 0.0 && x1v_(m,i) > krmax) return;
    const Real tl = (gen ? wtemp_(m,k-1,j,i) : w0(m,IEN,k-1,j,i)/w0(m,IDN,k-1,j,i)*gm1);
    const Real tr = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
    const Real pl = (gen ? wder_(m,IDPR,k-1,j,i) : w0(m,IEN,k-1,j,i)*gm1);
    const Real pr = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
    const Real dl = curv ? 0.5*(dx3_(m,k-1,j,i) + dx3_(m,k,j,i)) : size.d_view(m).dx3;
    Real wt = taumode ? 0.25*(wf(m,k-1,j,i) + wf(m,k-1,j,i+1)
                              + wf(m,k,j,i) + wf(m,k,j,i+1)) : 1.0;
    if (gaterho > 0.0) {
      wt *= RadGate(0.5*(w0(m,IDN,k-1,j,i) + w0(m,IDN,k,j,i))*dens_unit,
                    gaterho, gatedex);
    }
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
    Real wtc = wt;
    if (capang > 0.0) {
      const Real xm = fmax(capx(m,k-1,j,i), capx(m,k,j,i));
      if (xm > capang) {
        wtc = wt*(capang/xm);
        Kokkos::atomic_fetch_add(&capcnt(1), 1);
      }
    }
    flx3(m,IEN,k,j,i) += wtc*face_flux(tl, tr, pl, pr, w0(m,IDN,k-1,j,i),
                                       w0(m,IDN,k,j,i), gradn);
  });

  // rad_cap_ang report: how stiff the angular operator was and how much of it had to be
  // slowed down.  Same style and line budget as the rad_implicit_x1 report.
  if (capang > 0.0 && cap_lines < 400) {
    auto hcc = Kokkos::create_mirror_view(capcnt);
    Kokkos::deep_copy(hcc, capcnt);
    const bool active = (cap_diag_over > 0 || hcc(1) > 0);
    // the first 20 calls print from EVERY rank, so the startup lines carry the global
    // maximum of x_i and not just rank 0's share; after that only a rank that actually
    // had to cap something says anything
    if (active || cap_lines < 20) {
      ++cap_lines;
      auto hcr = Kokkos::create_mirror_view(caprec);
      Kokkos::deep_copy(hcr, caprec);
      const int mb = static_cast<int>(hcr(1));
      std::cout << "### rad_cap_ang rank " << global_variable::my_rank
                << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time
                << ": max x_i = " << hcr(0) << " (cap " << capang << ") at (m,k,j,i) = ("
                << mb << "," << static_cast<int>(hcr(2)) << ","
                << static_cast<int>(hcr(3)) << "," << static_cast<int>(hcr(4))
                << ") gid = " << (pmy_pack->gids + ((mb >= 0) ? mb : 0))
                << " r = " << hcr(5) << " | cells over cap = " << cap_diag_over
                << ", faces capped = " << hcc(1) << std::endl;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Conduction::ImplicitRadialUpdate
//! \brief BACKWARD-EULER solve of the radial flux-limited radiative diffusion, one
//! tridiagonal system per (m,k,j) column.  Runs after the RK update and the source
//! terms, and before the ghost exchange, so the values it writes are the ones that are
//! communicated.
//!
//! WHY.  The explicit radial operator is limited by dt ~ dr^2 rho c_v/kappa_rad with
//! kappa_rad = 16 sigma T^3/(3 kappa_R rho).  In an evacuated cell just above the
//! photosphere (rho 60x below its shell median, tau_cell ~ 1) that limit is ~1e-4 s
//! while the run takes 30 s steps, and since dt is set from the PREVIOUS state the
//! overshoot is unbounded: the cell heats past its neighbours, kappa_rad grows as T^3,
//! and T reaches 1e10 K in a single step.  Implicit in the ENERGY removes the
//! constraint entirely: the operator below is a symmetric M-matrix, so it is
//! unconditionally stable and cannot overshoot the neighbour equilibrium.
//!
//! WHAT IS SOLVED.  With e the internal energy density, T linearised about the frozen
//! post-update state,  T_i = T*_i + (e_i - e*_i)/(rho_i c_v,i),  and the frozen face
//! conductivity K_f (the SAME RadFaceKappa, tau-blend weight and flux limiter the
//! explicit face flux uses, the limiter evaluated on the frozen gradient),
//!     F_f = -K_f (T_i - T_{i-1})/dl_f
//!     e_i = e*_i + beta_dt (A_f F_f - A_{f+1} F_{f+1})/V_i
//! with exactly the areas and volumes Hydro::RKUpdate divides by, so the update is the
//! same flux-form divergence.  The two boundary faces (is and ie+1) are NOT in the
//! system -- they are added explicitly in AddIsotropicHeatFluxRadiative -- so the
//! interior fluxes telescope and sum_i V_i (e_i - e*_i) = 0 per column to round-off.
//! Only u0(IEN) is written; w0 is rebuilt by ConToPrim later in the same stage.

void Conduction::ImplicitRadialUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                                      const Real beta_dt) {
  if (!rad_implicit_x1) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  auto &size = pmy_pack->pmb->mb_size;
  const bool curvg = pmy_pack->pmesh->use_cubed_sphere
                     || pmy_pack->pmesh->use_spherical_polar;
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  auto &cosc_ = pmy_pack->pcoord->cos_cell;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &vol_ = pmy_pack->pcoord->volume;
  auto &area1_ = pmy_pack->pcoord->area.x1f;
  auto eos_ = eos;
  const bool gen = eos.IsGeneral();
  auto &wtemp_ = pmy_pack->phydro->wtemp;
  auto &phicc_ = pmy_pack->phydro->phicc0;
  const bool etg = pmy_pack->phydro->use_etotgrav;
  const Real temp_unit = pmy_pack->punit->temperature_cgs();
  const Real pres_unit = pmy_pack->punit->pressure_cgs();
  const Real dens_unit = pmy_pack->punit->density_cgs();
  const Real len_unit  = pmy_pack->punit->length_cgs();
  const Real eflx_unit = pres_unit*pmy_pack->punit->velocity_cgs();
  const Real met = rad_met, kfac = rad_kappa_fac;
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  const bool limit = rad_flux_limit;
  const Real krmax = rad_kappa_rmax;
  // the density gate, on the same face-averaged density the explicit x1 face uses
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  const Real sigma_sb = 5.670374419e-5;
  auto &wf = rad_w;
  const bool ktab = (rad_kappa_tab && rad_kr_nT > 0);
  const bool krho = rad_kappa_rho;
  auto &krt = rad_kr_tab;
  auto &krlT = rad_kr_lT;
  auto &krlP = rad_kr_lP;
  const int krnT = rad_kr_nT, krnP = rad_kr_nP;
  auto wrk = imp_wrk;
  auto iflag = imp_flag;
  auto irec = imp_rec;
  // slot indices as plain locals: a static constexpr member would capture `this`
  const int e_ = IMPE, t_ = IMPT, al_ = IMPA, pr_ = IMPP;
  const int c_ = IMPC, cp_ = IMPCP, dp_ = IMPDP;
  Kokkos::deep_copy(iflag, 0);

  const int nkj = (ke - ks + 1)*(je - js + 1);
  const int nj = (je - js + 1);
  Real maxviol = 0.0;
  int nfail = 0;
  Kokkos::parallel_reduce("radimpx1",
  Kokkos::RangePolicy<>(DevExeSpace(), 0, (nmb1 + 1)*nkj),
  KOKKOS_LAMBDA(const int &idx, Real &mviol, int &nbad) {
    const int m = idx/nkj;
    const int k = (idx - m*nkj)/nj + ks;
    const int j = (idx - m*nkj - (k - ks)*nj) + js;
    const Real dx1c = size.d_view(m).dx1;

    // ---- the FROZEN state of every cell: internal energy, temperature, pressure and
    // 1/(rho c_v).  The internal energy is extracted exactly as ConToPrim extracts it:
    // the gravitational term (etotgrav) and, on the cubed sphere, the kinetic energy
    // formed with the non-orthogonal metric (GnomonicEquiangleRaiseVel).  A cell whose
    // internal energy is not positive is marked with a negative 1/(rho c_v) and is
    // dropped from the system rather than handed to the EOS inversion.
    for (int i=is; i<=ie; ++i) {
      const Real d = u0(m,IDN,k,j,i);
      // the one shared extraction (utils/eint_from_cons.hpp): identical arithmetic in
      // identical order to what stood here, so this refactor is bit-for-bit a no-op
      const Real ei = EintFromCons(u0, m, k, j, i, cs_ ? cosc_(m,k,j) : 0.0, cs_, etg,
                                   etg ? phicc_(m,k,j,i) : 0.0);
      wrk(m,e_,k,j,i) = ei;
      wrk(m,t_,k,j,i) = 0.0;
      wrk(m,pr_,k,j,i) = 0.0;
      wrk(m,al_,k,j,i) = -1.0;
      if (!(ei > 0.0) || !(d > 0.0) || !isfinite(ei)) continue;
      const Real tt = eos_.Temperature(d, ei, gen ? wtemp_(m,k,j,i) : -1.0);
      if (!(tt > 0.0) || !isfinite(tt)) continue;
      const Real cv = eos_.SpecificHeatCv(d, ei, tt);
      if (!(cv > 0.0) || !isfinite(cv)) continue;
      wrk(m,t_,k,j,i) = tt;
      wrk(m,pr_,k,j,i) = eos_.Pressure(d, ei, tt);
      wrk(m,al_,k,j,i) = 1.0/(d*cv);
    }

    // ---- the frozen face coefficient A_f K_f/dl_f.  The two boundary faces are
    // outside the system: they were added explicitly with the ghost states.
    wrk(m,c_,k,j,is) = 0.0;
    wrk(m,c_,k,j,ie+1) = 0.0;
    for (int i=is+1; i<=ie; ++i) {
      Real ca = 0.0;
      const Real all = wrk(m,al_,k,j,i-1), alr = wrk(m,al_,k,j,i);
      if (all > 0.0 && alr > 0.0 && !(krmax > 0.0 && x1v_(m,i) > krmax)) {
        const Real tl = wrk(m,t_,k,j,i-1), tr = wrk(m,t_,k,j,i);
        const Real pf = 0.5*(wrk(m,pr_,k,j,i-1) + wrk(m,pr_,k,j,i));
        if (!(pf < pcut)) {
          const Real dl = curvg ? (x1v_(m,i) - x1v_(m,i-1)) : dx1c;
          const Real tk = 0.5*(tl + tr)*temp_unit;
          const Real rhof = 0.5*(u0(m,IDN,k,j,i-1) + u0(m,IDN,k,j,i))*dens_unit;
          const Real kap = RadFaceKappa(tk, pf*pres_unit, rhof, ktab, krt, krlT, krlP,
                                        krnT, krnP, krho, met, kfac);
          // the flux limiter, evaluated on the FROZEN gradient and then held fixed:
          // F = -kap g/sqrt(1 + (kap g/F_free)^2) linearises to a diffusion coefficient
          // kap/sqrt(1 + s^2) at fixed s, which is what keeps the system linear
          Real lf = 1.0;
          if (limit) {
            const Real ffree = sigma_sb*tk*tk*tk*tk;
            if (ffree > 0.0) {
              const Real fu = -kap*((tr - tl)/dl)*temp_unit/len_unit;
              lf = 1.0/sqrt(1.0 + SQR(fu/ffree));
            } else {
              lf = 0.0;
            }
          }
          Real wt = (taumode && blend_r) ? wf(m,k,j,i) : 1.0;
          if (gaterho > 0.0) wt *= RadGate(rhof, gaterho, gatedex);
          const Real af = curvg ? area1_(m,k,j,i) : 1.0;
          ca = wt*kap*lf*temp_unit/len_unit/eflx_unit*af/dl;
          if (!isfinite(ca) || ca < 0.0) ca = 0.0;
        }
      }
      wrk(m,c_,k,j,i) = ca;
    }

    // ---- forward sweep of the Thomas algorithm on the energy INCREMENT
    // x_i = e_i - e*_i.  Row i:
    //   -g_i A_i C_i alpha_{i-1} x_{i-1}
    // + (1 + g_i alpha_i (A_i C_i + A_{i+1} C_{i+1})) x_i
    // - g_i A_{i+1} C_{i+1} alpha_{i+1} x_{i+1}
    // = g_i (A_i F*_i - A_{i+1} F*_{i+1}),   g_i = beta_dt/V_i
    for (int i=is; i<=ie; ++i) {
      const Real vi = curvg ? vol_(m,k,j,i) : dx1c;
      const Real g = beta_dt/vi;
      const Real cl = wrk(m,c_,k,j,i), cr = wrk(m,c_,k,j,i+1);
      const Real ac = fmax(wrk(m,al_,k,j,i), 0.0);
      const Real amm = (i > is) ? fmax(wrk(m,al_,k,j,i-1), 0.0) : 0.0;
      const Real apl = (i < ie) ? fmax(wrk(m,al_,k,j,i+1), 0.0) : 0.0;
      const Real aa = -g*cl*amm;
      const Real bb = 1.0 + g*ac*(cl + cr);
      const Real cc = -g*cr*apl;
      const Real tc = wrk(m,t_,k,j,i);
      const Real tm = (i > is) ? wrk(m,t_,k,j,i-1) : 0.0;
      const Real tp = (i < ie) ? wrk(m,t_,k,j,i+1) : 0.0;
      const Real rhs = g*(cr*(tp - tc) - cl*(tc - tm));
      if (i == is) {
        wrk(m,cp_,k,j,i) = cc/bb;
        wrk(m,dp_,k,j,i) = rhs/bb;
      } else {
        const Real den = bb - aa*wrk(m,cp_,k,j,i-1);
        wrk(m,cp_,k,j,i) = cc/den;
        wrk(m,dp_,k,j,i) = (rhs - aa*wrk(m,dp_,k,j,i-1))/den;
      }
    }

    // ---- back substitution, write-back, and the per-column conservation residual
    Real xnext = 0.0, csum = 0.0, cabs = 0.0;
    for (int i=ie; i>=is; --i) {
      Real x = wrk(m,dp_,k,j,i) - wrk(m,cp_,k,j,i)*xnext;
      const Real es = wrk(m,e_,k,j,i);
      if (!isfinite(x) || !((es + x) > 0.0)) {
        // should never fire: the system is an M-matrix and cannot undershoot below the
        // minimum of the frozen column.  Fall back to no radial conduction in the cell.
        if (x != 0.0) {
          ++nbad;
          if (Kokkos::atomic_fetch_add(&iflag(0), 1) == 0) {
            irec(0) = static_cast<Real>(m);
            irec(1) = static_cast<Real>(k);
            irec(2) = static_cast<Real>(j);
            irec(3) = static_cast<Real>(i);
            irec(4) = x1v_(m,i);
            irec(5) = wrk(m,t_,k,j,i)*temp_unit;
            irec(6) = u0(m,IDN,k,j,i)*dens_unit;
            irec(7) = x;
          }
        }
        x = 0.0;
      }
      xnext = x;
      u0(m,IEN,k,j,i) += x;
      const Real vi = curvg ? vol_(m,k,j,i) : dx1c;
      csum += vi*x;
      cabs += fabs(vi*x);
    }
    if (cabs > 0.0) mviol = fmax(mviol, fabs(csum)/cabs);
  }, Kokkos::Max<Real>(maxviol), nfail);

  // <problem>/nan_report: the conservation residual of the tridiagonal solve, and any
  // cell that had to fall back.  Nothing is printed unless the switch is on, except a
  // fallback, which is always worth a line.
  if (imp_lines < 400) {
    const bool bad = (nfail > 0);
    const bool tell = nan_report && (imp_lines < 20 || maxviol > 1.0e-10);
    if (bad || (tell && global_variable::my_rank == 0)) {
      ++imp_lines;
      std::cout << "### rad_implicit_x1 rank " << global_variable::my_rank
                << " cycle " << pmy_pack->pmesh->ncycle
                << " t = " << pmy_pack->pmesh->time
                << ": max |sum V de|/sum V|de| = " << maxviol
                << ", fallback cells = " << nfail << std::endl;
      if (bad) {
        auto hr = Kokkos::create_mirror_view(irec);
        Kokkos::deep_copy(hr, irec);
        const int mb = static_cast<int>(hr(0));
        std::cout << "    first fallback (m,k,j,i) = (" << mb << ","
                  << static_cast<int>(hr(1)) << "," << static_cast<int>(hr(2)) << ","
                  << static_cast<int>(hr(3)) << ") gid = " << (pmy_pack->gids + mb)
                  << " r = " << hr(4) << " T* = " << hr(5) << " rho = " << hr(6)
                  << " de = " << hr(7) << std::endl;
      }
    }
  }
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
  const Real pcut = rad_tau_mode ? -1.0 : rad_pcut;
  const bool taumode = rad_tau_mode;
  const bool blend_r = rad_blend_radial;
  const bool limit = rad_flux_limit;
  // the radial operator is unconditionally stable when it is solved implicitly, so it
  // carries no timestep constraint; x2/x3 are still explicit and still do
  const bool impx1 = rad_implicit_x1;
  // ...and the angular operator self-limits when rad_cap_ang caps every face, so dt2 and
  // dt3 go with it
  const bool capa = (rad_cap_ang > 0.0);
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
  // the inert region carries no flux, so it must not carry a constraint either
  const Real krmax = rad_kappa_rmax;
  // ...and neither does gas the density gate has made inert: the face fluxes it feeds
  // are all scaled by G, so the cell's constraint scales with G too
  const Real gaterho = rad_gate_rho, gatedex = rad_gate_dex;
  auto &x1v_n = pmy_pack->pcoord->x1v;

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

    Real kappa_ = kappa0;
    Real wmax = 1.0;
    if (spitzer) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      kappa_ = TempDepKappa(temp*temp_unit, limit_)/kappa_unit;
    } else if (radiative) {
      Real temp = (gen ? wtemp_(m,k,j,i) : w0(m,IEN,k,j,i)/w0(m,IDN,k,j,i)*gm1);
      Real pres = (gen ? wder_(m,IDPR,k,j,i) : w0(m,IEN,k,j,i)*gm1);
      if (pres < pcut) return;   // no flux above the cut: no constraint
      if (krmax > 0.0 && x1v_n(m,i) > krmax) return;     // radiatively inert corona
      kappa_ = (ktab
          ? RadiativeKappaKR(temp*temp_unit, w0_(m,IDN,k,j,i)*dens_unit, kfac,
                             RosselandTable(krt, krlT, krlP, krnT, krnP, temp*temp_unit,
                                            krho ? w0_(m,IDN,k,j,i)*dens_unit
                                                 : pres*pres_unit))
          : RadiativeKappa(temp*temp_unit, pres*pres_unit, w0_(m,IDN,k,j,i)*dens_unit,
                           met, kfac))/kappa_unit;
      if (gaterho > 0.0) {
        kappa_ *= RadGate(w0_(m,IDN,k,j,i)*dens_unit, gaterho, gatedex);
      }
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

    const Real d1 = (curv && radiative) ? dx1_(m,k,j,i) : size.d_view(m).dx1;
    if (!impx1) {
      const Real w1 = (taumode && blend_r) ? wmax : 1.0;
      const Real k1 = w1*keff(d1, tc(k,j,i-1), tc(k,j,i+1));
      if (k1 > 0.0) min_dt = fmin(min_dt, SQR(d1)/k1*rcv);
    }
    // on a curvilinear grid size.dx2/dx3 are ANGLES; the physical widths are pcoord's
    // cubed sphere: the exact operator's angular diffusivity is kappa/sin^2(alpha)
    const Real s2 = (cs && radiative) ? SQR(sinc_(m,k,j)) : 1.0;
    const Real wa = taumode ? wmax : 1.0;
    if (multi_d && !capa) {
      const Real d2 = (curv && radiative) ? dx2_(m,k,j,i) : size.d_view(m).dx2;
      const Real k2 = wa*keff(d2, tc(k,j-1,i), tc(k,j+1,i));
      if (k2 > 0.0) min_dt = fmin(min_dt, SQR(d2)*s2/k2*rcv);
    }
    if (three_d && !capa) {
      const Real d3 = (curv && radiative) ? dx3_(m,k,j,i) : size.d_view(m).dx3;
      const Real k3 = wa*keff(d3, tc(k-1,j,i), tc(k+1,j,i));
      if (k3 > 0.0) min_dt = fmin(min_dt, SQR(d3)*s2/k3*rcv);
    }
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
      const Real kr = (ktab ? RosselandTable(krt, krlT, krlP, krnT, krnP, temp*temp_unit,
                                             krho ? dens*dens_unit : pres*pres_unit)
                            : -1.0);
      Real kappa_ = (ktab
          ? RadiativeKappaKR(temp*temp_unit, dens*dens_unit, kfac, kr)
          : RadiativeKappa(temp*temp_unit, pres*pres_unit, dens*dens_unit, met, kfac))
          /kappa_unit;
      if (gaterho > 0.0) kappa_ *= RadGate(dens*dens_unit, gaterho, gatedex);
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
      dd.d_view(12) = (impx1 || !(w1*keff(s1v) > 0.0))
                      ? -1.0 : SQR(d1)/(w1*keff(s1v))*rcv*fac;
      dd.d_view(13) = (capa || !(multi_d && wa*keff(s2v) > 0.0))
                      ? -1.0 : SQR(d2)*s2/(wa*keff(s2v))*rcv*fac;
      dd.d_view(14) = (capa || !(three_d && wa*keff(s3v) > 0.0))
                      ? -1.0 : SQR(d3)*s2/(wa*keff(s3v))*rcv*fac;
      dd.d_view(15) = ffree;
    });
    dt_diag.template modify<DevExeSpace>();
    dt_diag.template sync<HostMemSpace>();
    dt_diag_valid = true;
  }
  dtnew_prev = dtnew;

  return;
}
