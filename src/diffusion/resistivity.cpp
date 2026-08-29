//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file resistivity.cpp
//  \brief Implements functions for Resistivity class.

#include <float.h>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string> // string

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "resistivity.hpp"
#include "current_density.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls Resistivity base class constructor

Resistivity::Resistivity(MeshBlockPack *pp, ParameterInput *pin) :
    pmy_pack(pp),
    efld_resist("efld_resist",1,1,1,1),
    jnorm("jnorm",1,1,1,1),
    eta_b("eta_b",1,1,1,1),
    u_ideal("u_ideal",1,1,1,1,1),
    b_ideal("b_ideal",1,1,1,1),
    uflx_ideal("uflx_ideal",1,1,1,1,1),
    u2("u2",1,1,1,1,1),
    b2("b2",1,1,1,1),
    efld_ideal("efld_ideal",1,1,1,1) {
  // Read parameters for Ohmic resistivity (if any)
  if (pin->DoesParameterExist("mhd","ohmic_resistivity")) {
    iso_resist_type = pin->GetString("mhd","ohmic_resistivity");
    // Check for valid type
    if (iso_resist_type.compare("constant") != 0 && iso_resist_type.compare("perna") != 0
        && iso_resist_type.compare("eos") != 0) {
      std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
                << "Invalid choice for Ohmic resistivity type '" << iso_resist_type
                << "'; expected constant, perna or eos" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    use_rkg_sts = pin->GetOrAddBoolean("mhd", "use_rkg_sts", false);
      
    // Total number of MeshBlocks on this rank to be used in array dimensioning
    int nmb = std::max((pp->nmb_thispack), (pp->pmesh->nmb_maxperrank));
    auto &indcs = pp->pmesh->mb_indcs;
    int ncells1 = indcs.nx1 + 2*(indcs.ng);
    int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
    int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
    Kokkos::realloc(efld_resist.x1e, nmb, ncells3+1, ncells2+1, ncells1);
    Kokkos::realloc(efld_resist.x2e, nmb, ncells3+1, ncells2, ncells1+1);
    Kokkos::realloc(efld_resist.x3e, nmb, ncells3, ncells2+1, ncells1+1);
    // the cubed sphere takes the curl in two passes and needs somewhere to put the
    // face-normal-frame current between them; see resistivity_gnomonic.cpp
    if (pp->pmesh->use_cubed_sphere) {
      Kokkos::realloc(jnorm.x1e, nmb, ncells3+1, ncells2+1, ncells1);
      Kokkos::realloc(jnorm.x2e, nmb, ncells3+1, ncells2, ncells1+1);
      Kokkos::realloc(jnorm.x3e, nmb, ncells3, ncells2+1, ncells1+1);
    }
    // eta_b holds the diffusivity of every cell for EVERY resistivity type, because
    // AddEMFGeneralResist() is the only EMF routine left and it reads eta_b for every
    // cell (AddEMFConstantResist is commented out). A constant resistivity therefore has
    // to fill the array too: leaving it at its 1x1x1x1 placeholder made each read run off
    // the end of the allocation, which silently returned zero and left eta_ohm_const with
    // no effect whatsoever.
    Kokkos::realloc(eta_b, nmb, ncells3, ncells2, ncells1);
    if (iso_resist_type.compare("constant") == 0) {
      // constant resistivity
      eta_ohm_const = pin->GetReal("mhd","eta_ohm_const");
      Kokkos::deep_copy(eta_b, eta_ohm_const);
    } else {
//      min_xe = pin->GetReal("mhd","min_xe");
      max_eta = pin->GetReal("mhd","max_eta");
    }
    if (use_rkg_sts) {
        int nmhd;
        // (1) construct EOS object (no default)
        std::string eqn_of_state = pin->GetString("mhd","eos");
        // ideal gas EOS
        if (eqn_of_state.compare("ideal") == 0) {
          nmhd = 5;
        // general EOS: same five conserved variables as an ideal gas, it is only their
        // relation to the primitives that differs
        } else if (eqn_of_state.compare("general") == 0) {
          nmhd = 5;
        // isothermal EOS
        } else if (eqn_of_state.compare("isothermal") == 0) {
          nmhd = 4;
        } else {
          std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__
                    << std::endl << "Super-timestepped resistivity does not support "
                    << "<mhd>/eos = " << eqn_of_state << std::endl;
          std::exit(EXIT_FAILURE);
        }
        int nscalars = pin->GetOrAddInteger("mhd","nscalars",0);
        Kokkos::realloc(u_ideal,     nmb, (nmhd+nscalars), ncells3, ncells2, ncells1);
        Kokkos::realloc(b_ideal.x1f, nmb, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(b_ideal.x2f, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(b_ideal.x3f, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::realloc(uflx_ideal.x1f, nmb, (nmhd+nscalars), ncells3, ncells2, ncells1+1);
        Kokkos::realloc(uflx_ideal.x2f, nmb, (nmhd+nscalars), ncells3, ncells2+1, ncells1);
        Kokkos::realloc(uflx_ideal.x3f, nmb, (nmhd+nscalars), ncells3+1, ncells2, ncells1);
        Kokkos::realloc(u2,     nmb, (nmhd+nscalars), ncells3, ncells2, ncells1);
        Kokkos::realloc(b2.x1f, nmb, ncells3, ncells2, ncells1+1);
        Kokkos::realloc(b2.x2f, nmb, ncells3, ncells2+1, ncells1);
        Kokkos::realloc(b2.x3f, nmb, ncells3+1, ncells2, ncells1);
        Kokkos::realloc(efld_ideal.x1e, nmb, ncells3+1, ncells2+1, ncells1);
        Kokkos::realloc(efld_ideal.x2e, nmb, ncells3+1, ncells2, ncells1+1);
        Kokkos::realloc(efld_ideal.x3e, nmb, ncells3, ncells2+1, ncells1+1);
    }
  }
}

//----------------------------------------------------------------------------------------
// Resistivity destructor

Resistivity::~Resistivity() {
}

void Resistivity::SetResistivity(const DvceArray5D<Real> &w, const EOS_Data &eos, const Real &Rgas, DvceArray4D<Real> &eta_b, const int il, const int iu, const int jl, const int ju, const int kl, const int ku) {
  int &nmb = pmy_pack->nmb_thispack;
  Real gm1 = eos.gamma - 1.0;
  // ResistivityPerna is a thermal-ionization fit and needs a temperature in KELVIN and a
  // particle number density. For an ideal gas both come from the problem's Rgas, which
  // carries a fixed mean molecular weight. A general EOS supplies temperature and
  // composition itself, so mu is asked of the EOS. The number density is then written as
  // rho pres_cgs/(mu temp_cgs k_B) rather than rho/(mu m_u): eliminating m_u in favour of
  // eos.temp_cgs = (pres_cgs/dens_cgs) m_u/k_B ties the number density to the same
  // normalisation the temperature uses, instead of introducing a second atomic mass unit
  // that would disagree with this file's k_B in the last digits.
  const bool gen = eos.IsGeneral();
  const bool use_eos_xe = (iso_resist_type.compare("eos") == 0);
  // The electron fraction has to come from somewhere: either the general EOS carries it,
  // or an ideal-gas run had the composition table built for this alone (see
  // EquationOfState::BuildElectronFractionTable). Fail here rather than silently return
  // zero and float every cell to max_eta.
  if (use_eos_xe && !eos.HasElectronFraction()) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "ohmic_resistivity = eos found no electron fraction. It needs either "
              << "eos = general with general_eos = table, or an ideal gas with a <units> "
              << "block and the eos_* composition parameters so the table can be built "
              << "for x_e alone. Either way set eos_metal_ionization = true, or the "
              << "atmosphere will have essentially no electrons below ~4000 K"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // The ideal-gas branch below forms T = p/(Rgas rho) from <problem>/Rgas, which only a
  // hot_jupiter problem sets; it is zero for anything else. Catch that here rather than
  // divide by it per cell.
  if (!gen && !(Rgas > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "the Ohmic resistivity needs a temperature, and on the ideal-gas branch "
              << "that is p/(<problem>/Rgas rho). Rgas is " << Rgas << ", so the "
              << "problem generator never set it -- only a <problem>/hot_jupiter = true "
              << "problem does." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // the temperature is read from the cache ConsToPrim filled for this same primitive
  // state, not re-derived: T(d,e) is a root find for a general EOS and doing it here
  // would repeat, per cell per stage, work that has already been done
  auto &wtemp_ = pmy_pack->pmhd->wtemp;
  // local copy: reading the member inside the kernel would capture `this`
  const Real max_eta_ = max_eta;
  par_for("mhd_resistval", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real rho = w(m,IDN,k,j,i);
    Real e = w(m,IEN,k,j,i);
    Real kb = 1.380649e-16;
    Real mh = 1.67262192369e-24;
    Real T, nn;
    if (gen) {
      Real tcode = wtemp_(m,k,j,i);
      T = tcode*eos.temp_cgs;
      nn = rho*eos.pres_cgs/(eos.MeanMolecularWeight(rho,e,tcode)*eos.temp_cgs*kb);
    } else {
      Real p = gm1*e;
      T = p/Rgas/rho;
      Real mu = kb/mh/Rgas;
      nn = rho/(mu*mh);
    }
    if (use_eos_xe) {
      // T is already in kelvin on both branches above, so the kelvin entry point is the
      // one to use -- the (d,e,T) form would scale by EOS_Data::temp_cgs, which an ideal
      // gas never sets. For a general EOS the two are the same number.
      ResistivityEOS(eos.ElectronFractionKelvin(rho, T), T, max_eta_, eta_b(m,k,j,i));
    } else {
      ResistivityPerna(nn, T, max_eta_, eta_b(m,k,j,i));
    }
//      Real lgrho = log10(rho);
//      Real lgT = log10(T);
//      lgrho = (lgrho < -7.0)? -7.0 : lgrho;
//      lgrho = (lgrho > -2.0)? -2.0 : lgrho;
//      lgT = (lgT < 3.0) ? 3.0 : lgT;
//      lgT = (lgT > 4.542780748663102) ?  4.542780748663102 : lgT;
//      lgrho = (lgrho+7.0)/(-2.0+7.0);
//      lgT = (lgT-3.0)/(4.542780748663102-3.0);
//    ResistivityKumar(lgrho, lgT, eta_b(m,k,j,i));
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddResistiveEMFs()
//! \brief Wrapper function that adds electric fields for different types of resistivity
//! Currently only Ohmic resistivity with constant coefficient is implemented.

void Resistivity::AddResistiveEMFs(const DvceFaceFld4D<Real> &b0,
    DvceEdgeFld4D<Real> &efld) {
//  AddEMFConstantResist(b0, efld);
  AddEMFDirect(efld_resist, efld);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void AddResistiveFluxes()
//! \brief Wrapper function that adds energy (Poynting) fluxes for different types of
//! resistivity.
//! Currently only Ohmic resistivity with constant coefficient is implemented.

void Resistivity::AddResistiveFluxes(const DvceFaceFld4D<Real> &b0, const DvceArray5D<Real> &bc, DvceFaceFld5D<Real> &flx) {
//  AddFluxConstantGridResist(b0, flx);
  AddFluxGeneralResist(b0, bc, flx);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn AddEMFConstantResist()
//  \brief Adds electric field from Ohmic resistivity to corner-centered electric field
//  Using Ohm's Law to compute the electric field:  E + (v x B) = \eta J, then
//    E_{inductive} = - (v x B)  [computed in the MHD Riemann solver]
//    E_{resistive} = \eta J     [computed in this function]

void Resistivity::AddEMFConstantResist(const DvceFaceFld4D<Real> &b0,
    DvceEdgeFld4D<Real> &efld) {
  // DEAD CODE TRIPWIRE. Nothing calls this (the call in AddResistiveFluxes is commented
  // out), but it reaches CurrentDensity(), whose only curvilinear branch is the
  // spherical-polar one -- so on the cubed sphere it would take the CARTESIAN curl and
  // be wrong by O(1), which is exactly the bug AddEMFGnomonicResist exists to fix. If
  // this is ever re-enabled it must be given the same two-pass treatment.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
              << "AddEMFConstantResist has no cubed-sphere form; see "
              << "diffusion/resistivity_gnomonic.cpp" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // gathered on the host; the kernels must not chase pmy_pack on the device
  const auto geom = MakeCurrentDensityGeom(pmy_pack);

  //---- 1-D problem:
  //  copy face-centered E-fields to edges and return.
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto &mbsize = pmy_pack->pmb->mb_size;
    auto eta_o = eta_ohm_const;

    int scr_level = 0;
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

    par_for_outer("ohm1", DevExeSpace(), scr_size, scr_level, 0, nmb1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

      CurrentDensity(geom, member, m, ks, js, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

      // Add E_{resistive} = \eta J to corner-centered electric fields
      par_for_inner(member, is, ie+1, [&](const int i) {
        e2(m,ks,  js  ,i) += eta_o*j2(i);
        e2(m,ke+1,js  ,i) += eta_o*j2(i);
        e3(m,ks  ,js  ,i) += eta_o*j3(i);
        e3(m,ks  ,je+1,i) += eta_o*j3(i);
      });
    });
    return;
  }

  //---- 2-D problem:
  if (pmy_pack->pmesh->two_d) {
    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto &mbsize = pmy_pack->pmb->mb_size;
    auto eta_o = eta_ohm_const;

    int scr_level = 0;
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

    par_for_outer("ohm2", DevExeSpace(), scr_size, scr_level, 0, nmb1, js, je+1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

      CurrentDensity(geom, member, m, ks, j, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

      // Add E_{resistive} = \eta J to corner-centered electric fields
      par_for_inner(member, is, ie+1, [&](const int i) {
        e1(m,ks,  j,i) += eta_o*j1(i);
        e1(m,ke+1,j,i) += eta_o*j1(i);
        e2(m,ks,  j,i) += eta_o*j2(i);
        e2(m,ke+1,j,i) += eta_o*j2(i);
        e3(m,ks  ,j,i) += eta_o*j3(i);
      });
    });
    return;
  }

  //---- 3-D problem:

  // capture class variables for the kernels
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto eta_o = eta_ohm_const;

  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

  par_for_outer("ohm3", DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke+1, js, je+1,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

    CurrentDensity(geom, member, m, k, j, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

    // Add E_{resistive} = \eta J to corner-centered electric fields
    par_for_inner(member, is, ie+1, [&](const int i) {
      e1(m,k,j,i) += eta_o*j1(i);
      e2(m,k,j,i) += eta_o*j2(i);
      e3(m,k,j,i) += eta_o*j3(i);
    });
  });

  return;
}

void Resistivity::AddEMFGeneralResist(const DvceFaceFld4D<Real> &b0,
    DvceEdgeFld4D<Real> &efld) {
  // the cubed sphere needs a two-pass curl on a non-orthogonal basis; see
  // resistivity_gnomonic.cpp for why the formulae below cannot simply be given a metric
  if (pmy_pack->pmesh->use_cubed_sphere) {
    AddEMFGnomonicResist(b0, efld);
    return;
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // gathered on the host; the kernels must not chase pmy_pack on the device
  const auto geom = MakeCurrentDensityGeom(pmy_pack);
  // local copy of the member View, so the lambdas do not capture `this`
  auto eta_b_ = eta_b;

  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;
  auto &x2v_ = pmy_pack->pcoord->x2v;
  auto &x2f_ = pmy_pack->pcoord->xx2f;
  auto &x3v_ = pmy_pack->pcoord->x3v;
  auto &x3f_ = pmy_pack->pcoord->xx3f;

  //---- 1-D problem:
  //  copy face-centered E-fields to edges and return.
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto &mbsize = pmy_pack->pmb->mb_size;

    int scr_level = 0;
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

    par_for_outer("ohm1", DevExeSpace(), scr_size, scr_level, 0, nmb1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

      CurrentDensity(geom, member, m, ks, js, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

      // Add E_{resistive} = \eta J to corner-centered electric fields.
      // The edge fields carry a k = ke+1 and a j = je+1 layer even in 1D, but eta_b is
      // cell centered and has only the ONE cell in each of those directions: indexing it
      // at ke+1 or je+1 reads past the end of the array. The two edges bound the same
      // cell, so both take that cell's eta.
      par_for_inner(member, is, ie+1, [&](const int i) {
        Real etaf = 0.5*(eta_b_(m,ks,js,i) + eta_b_(m,ks,js,i-1));
        e2(m,ks,  js  ,i) += etaf*j2(i);
        e2(m,ke+1,js  ,i) += etaf*j2(i);
        e3(m,ks  ,js  ,i) += etaf*j3(i);
        e3(m,ks  ,je+1,i) += etaf*j3(i);
      });
    });
    return;
  }

  //---- 2-D problem:
  if (pmy_pack->pmesh->two_d) {
    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto &mbsize = pmy_pack->pmb->mb_size;
    auto eta_o = eta_ohm_const;

    int scr_level = 0;
    size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

    par_for_outer("ohm2", DevExeSpace(), scr_size, scr_level, 0, nmb1, js, je+1,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j) {
      ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
      ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

      CurrentDensity(geom, member, m, ks, j, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

      // Add E_{resistive} = \eta J to corner-centered electric fields.
      // As in the 1D branch: the k = ke+1 edge layer exists but the ke+1 CELL does not,
      // so the eta of the single x3 cell is used on both x3 faces.
      par_for_inner(member, is, ie+1, [&](const int i) {
        e1(m,ks,  j,i) += 0.5*(eta_b_(m,ks,j,i)+eta_b_(m,ks,j-1,i))*j1(i);
        e1(m,ke+1,j,i) += 0.5*(eta_b_(m,ks,j,i)+eta_b_(m,ks,j-1,i))*j1(i);
        e2(m,ks,  j,i) += 0.5*(eta_b_(m,ks,j,i)+eta_b_(m,ks,j,i-1))*j2(i);
        e2(m,ke+1,j,i) += 0.5*(eta_b_(m,ks,j,i)+eta_b_(m,ks,j,i-1))*j2(i);
        e3(m,ks  ,j,i) += 0.25*(eta_b_(m,ks,j,i)+eta_b_(m,ks,j,i-1)+eta_b_(m,ks,j-1,i)+eta_b_(m,ks,j-1,i-1))*j3(i);
      });
    });
    return;
  }

  //---- 3-D problem:

  // capture class variables for the kernels
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto eta_o = eta_ohm_const;

  int scr_level = 0;
  size_t scr_size = ScrArray1D<Real>::shmem_size(ncells1) * 3;

  par_for_outer("ohm3", DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke+1, js, je+1,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray1D<Real> j1(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> j2(member.team_scratch(scr_level), ncells1);
    ScrArray1D<Real> j3(member.team_scratch(scr_level), ncells1);

    CurrentDensity(geom, member, m, k, j, is, ie+1, b0, mbsize.d_view(m), j1, j2, j3);

    // Add E_{resistive} = \eta J to corner-centered electric fields
    par_for_inner(member, is, ie+1, [&](const int i) {
      e1(m,k,j,i) += 0.25*(eta_b_(m,k,j,i)+eta_b_(m,k,j-1,i)+eta_b_(m,k-1,j,i)+eta_b_(m,k-1,j-1,i))*j1(i);
      e2(m,k,j,i) += 0.25*(eta_b_(m,k,j,i)+eta_b_(m,k,j,i-1)+eta_b_(m,k-1,j,i)+eta_b_(m,k-1,j,i-1))*j2(i);
      e3(m,k,j,i) += 0.25*(eta_b_(m,k,j,i)+eta_b_(m,k,j-1,i)+eta_b_(m,k,j,i-1)+eta_b_(m,k,j-1,i-1))*j3(i);
    });
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn AddResistiveFluxConstantResist()
//  \brief Adds Poynting flux from Ohmic resistivity to energy flux
//  Total energy equation is dE/dt = - Div(F) where F = (E X B) = \eta (J X B)

void Resistivity::AddFluxConstantGridResist(const DvceFaceFld4D<Real> &b,
                                        DvceFaceFld5D<Real> &flx) {
  // DEAD CODE TRIPWIRE, as in AddEMFConstantResist: this crosses the edge EMFs with the
  // cell-centred field with no regard for the frame each is held in.
  if (pmy_pack->pmesh->use_cubed_sphere) {
    std::cout << "### FATAL ERROR in "<< __FILE__ <<" at line " << __LINE__ << std::endl
              << "AddFluxConstantGridResist has no cubed-sphere form; see "
              << "AddFluxGeneralResist" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;
  Real qa = 0.25*eta_ohm_const;

  //------------------------------
  // energy fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("ohm_heat1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real j2k   = -(b.x3f(m,k  ,j,i) - b.x3f(m,k  ,j,i-1))/size.d_view(m).dx1;
    Real j2kp1 = -(b.x3f(m,k+1,j,i) - b.x3f(m,k+1,j,i-1))/size.d_view(m).dx1;

    Real j3j   = (b.x2f(m,k,j  ,i) - b.x2f(m,k,j  ,i-1))/size.d_view(m).dx1;
    Real j3jp1 = (b.x2f(m,k,j+1,i) - b.x2f(m,k,j+1,i-1))/size.d_view(m).dx1;

    if (multi_d) {
      j3j   -= (b.x1f(m,k,j  ,i) - b.x1f(m,k,j-1,i))/size.d_view(m).dx2;
      j3jp1 -= (b.x1f(m,k,j+1,i) - b.x1f(m,k,j  ,i))/size.d_view(m).dx2;
    }
    if (three_d) {
      j2k   += (b.x1f(m,k  ,j,i) - b.x1f(m,k-1,j,i))/size.d_view(m).dx3;
      j2kp1 += (b.x1f(m,k+1,j,i) - b.x1f(m,k  ,j,i))/size.d_view(m).dx3;
    }

    // flx1 = (E X B)_{1} =  ((\eta J) X B)_{1} = \eta (J2*B3 - J3*B2)
    flx1(m,IEN,k,j,i) += qa*(j2k  *(b.x3f(m,k  ,j  ,i) + b.x3f(m,k  ,j  ,i-1)) +
                             j2kp1*(b.x3f(m,k+1,j  ,i) + b.x3f(m,k+1,j  ,i-1)) -
                             j3j  *(b.x2f(m,k  ,j  ,i) + b.x2f(m,k  ,j  ,i-1)) -
                             j3jp1*(b.x2f(m,k  ,j+1,i) + b.x2f(m,k  ,j+1,i-1)));
  });
  if (pmy_pack->pmesh->one_d) {return;}

  //------------------------------
  // energy fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("ohm_heat2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real j1k   = (b.x3f(m,k  ,j,i) - b.x3f(m,k  ,j-1,i))/size.d_view(m).dx2;
    Real j1kp1 = (b.x3f(m,k+1,j,i) - b.x3f(m,k+1,j-1,i))/size.d_view(m).dx2;

    Real j3i   = (b.x2f(m,k,j,i  ) - b.x2f(m,k,j  ,i-1))/size.d_view(m).dx1
               - (b.x1f(m,k,j,i  ) - b.x1f(m,k,j-1,i  ))/size.d_view(m).dx2;
    Real j3ip1 = (b.x2f(m,k,j,i+1) - b.x2f(m,k,j  ,i  ))/size.d_view(m).dx1
               - (b.x1f(m,k,j,i+1) - b.x1f(m,k,j-1,i+1))/size.d_view(m).dx2;

    if (three_d) {
      j1k   -= (b.x2f(m,k  ,j,i) - b.x2f(m,k-1,j,i))/size.d_view(m).dx3;
      j1kp1 -= (b.x2f(m,k+1,j,i) - b.x2f(m,k  ,j,i))/size.d_view(m).dx3;
    }

    // E2 = \eta (J X B)_{2} = \eta (J3*B1 - J1*B3)
    flx2(m,IEN,k,j,i) += qa*(j3i  *(b.x1f(m,k  ,j,i  ) + b.x1f(m,k  ,j-1,i  )) +
                             j3ip1*(b.x1f(m,k  ,j,i+1) + b.x1f(m,k  ,j-1,i+1)) -
                             j1k  *(b.x3f(m,k  ,j,i  ) + b.x3f(m,k  ,j-1,i  )) -
                             j1kp1*(b.x3f(m,k+1,j,i  ) + b.x3f(m,k+1,j-1,i  )));
  });
  if (pmy_pack->pmesh->two_d) {return;}

  //------------------------------
  // energy fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("ohm_heat3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real j1j   = (b.x3f(m,k,j  ,i) - b.x3f(m,k  ,j-1,i))/size.d_view(m).dx2
               - (b.x2f(m,k,j  ,i) - b.x2f(m,k-1,j  ,i))/size.d_view(m).dx3;
    Real j1jp1 = (b.x3f(m,k,j+1,i) - b.x3f(m,k  ,j  ,i))/size.d_view(m).dx2
               - (b.x2f(m,k,j+1,i) - b.x2f(m,k-1,j+1,i))/size.d_view(m).dx3;

    Real j2i   = -(b.x3f(m,k,j,i  ) - b.x3f(m,k  ,j,i-1))/size.d_view(m).dx1
                + (b.x1f(m,k,j,i  ) - b.x1f(m,k-1,j,i  ))/size.d_view(m).dx3;
    Real j2ip1 = -(b.x3f(m,k,j,i+1) - b.x3f(m,k  ,j,i  ))/size.d_view(m).dx1
                + (b.x1f(m,k,j,i+1) - b.x1f(m,k-1,j,i+1))/size.d_view(m).dx3;

    // E2 = \eta (J X B)_{2} = \eta (J1*B2 - J2*B1)
    flx3(m,IEN,k,j,i) += qa*(j1j  *(b.x2f(m,k,j  ,i  ) + b.x2f(m,k-1,j  ,i  )) +
                             j1jp1*(b.x2f(m,k,j+1,i  ) + b.x2f(m,k-1,j+1,i  )) -
                             j2i  *(b.x1f(m,k,j  ,i  ) + b.x1f(m,k-1,j  ,i  )) -
                             j2ip1*(b.x1f(m,k,j  ,i+1) + b.x1f(m,k-1,j  ,i+1)));
  });

  return;
}

void Resistivity::ClearResistiveEMFs(DvceEdgeFld4D<Real> &efld) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  //---- 1-D problem:
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;

    par_for("ohm10", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int i) {
      e2(m,ks,  js  ,i) = 0.0;
      e2(m,ke+1,js  ,i) = 0.0;
      e3(m,ks  ,js  ,i) = 0.0;
      e3(m,ks  ,je+1,i) = 0.0;
    });
    return;
  }

  //---- 2-D problem:
  if (pmy_pack->pmesh->two_d) {
    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;

    par_for("ohm20", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      e1(m,ks,  j,i) = 0.0;
      e1(m,ke+1,j,i) = 0.0;
      e2(m,ks,  j,i) = 0.0;
      e2(m,ke+1,j,i) = 0.0;
      e3(m,ks  ,j,i) = 0.0;
    });
    return;
  }

  //---- 3-D problem:

  // capture class variables for the kernels
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;

  par_for("ohm30", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    e1(m,k,j,i) = 0.0;
    e2(m,k,j,i) = 0.0;
    e3(m,k,j,i) = 0.0;
  });

  return;
}

void Resistivity::AddEMFDirect(const DvceEdgeFld4D<Real> &efld_resist, DvceEdgeFld4D<Real> &efld) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  // local copy: reading the member inside the kernels would capture `this`
  const bool use_sts = use_rkg_sts;

  //---- 1-D problem:
  //  Note e2[is:ie+1,js:je,  ks:ke+1]
  //       e3[is:ie+1,js:je+1,ks:ke  ]

  if (pmy_pack->pmesh->one_d) {
    // capture class variables for the kernels
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e2r = efld_resist.x2e;
    auto e3r = efld_resist.x3e;

    par_for("ohm1add", DevExeSpace(), 0, nmb1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int i) {
      if (use_sts) {
        e2(m,ks,  js  ,i) = e2r(m,ks,  js  ,i);
        e2(m,ke+1,js  ,i) = e2r(m,ke+1,js  ,i);
        e3(m,ks  ,js  ,i) = e3r(m,ks  ,js  ,i);
        e3(m,ks  ,je+1,i) = e3r(m,ks  ,je+1,i);
      } else {
        e2(m,ks,  js  ,i) += e2r(m,ks,  js  ,i);
        e2(m,ke+1,js  ,i) += e2r(m,ke+1,js  ,i);
        e3(m,ks  ,js  ,i) += e3r(m,ks  ,js  ,i);
        e3(m,ks  ,je+1,i) += e3r(m,ks  ,je+1,i);
      }
    });
    return;
  }

  //---- 2-D problem:
  if (pmy_pack->pmesh->two_d) {
    // capture class variables for the kernels
    auto e1 = efld.x1e;
    auto e2 = efld.x2e;
    auto e3 = efld.x3e;
    auto e1r = efld_resist.x1e;
    auto e2r = efld_resist.x2e;
    auto e3r = efld_resist.x3e;

    par_for("ohm2add", DevExeSpace(), 0, nmb1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int j, const int i) {
      if (use_sts) {
        e1(m,ks,  j,i) = e1r(m,ks,  j,i);
        e1(m,ke+1,j,i) = e1r(m,ke+1,j,i);
        e2(m,ks,  j,i) = e2r(m,ks,  j,i);
        e2(m,ke+1,j,i) = e2r(m,ke+1,j,i);
        e3(m,ks  ,j,i) = e3r(m,ks  ,j,i);
      } else {
        e1(m,ks,  j,i) += e1r(m,ks,  j,i);
        e1(m,ke+1,j,i) += e1r(m,ke+1,j,i);
        e2(m,ks,  j,i) += e2r(m,ks,  j,i);
        e2(m,ke+1,j,i) += e2r(m,ke+1,j,i);
        e3(m,ks  ,j,i) += e3r(m,ks  ,j,i);
      }
    });
    return;
  }

  //---- 3-D problem:

  // capture class variables for the kernels
  auto e1 = efld.x1e;
  auto e2 = efld.x2e;
  auto e3 = efld.x3e;
  auto e1r = efld_resist.x1e;
  auto e2r = efld_resist.x2e;
  auto e3r = efld_resist.x3e;

  par_for("ohm3add", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    if (use_sts) {
      e1(m,k,j,i) = e1r(m,k,j,i);
      e2(m,k,j,i) = e2r(m,k,j,i);
      e3(m,k,j,i) = e3r(m,k,j,i);
    } else {
      e1(m,k,j,i) += e1r(m,k,j,i);
      e2(m,k,j,i) += e2r(m,k,j,i);
      e3(m,k,j,i) += e3r(m,k,j,i);
    }
  });

  return;
}

void Resistivity::AddFluxGeneralResist(const DvceFaceFld4D<Real> &b, const DvceArray5D<Real> &bc, DvceFaceFld5D<Real> &flx) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto size = pmy_pack->pmb->mb_size;
  const bool one_d = pmy_pack->pmesh->one_d;
  const bool two_d = pmy_pack->pmesh->two_d;
  const bool three_d = pmy_pack->pmesh->three_d;
    
  ClearResistiveEMFs(efld_resist);
//  if (iso_resist_type.compare("constant") == 0) {
//    AddEMFConstantResist(b, efld_resist);
//  } else {
    AddEMFGeneralResist(b, efld_resist);
//  }
  auto &e1 = efld_resist.x1e;
  auto &e2 = efld_resist.x2e;
  auto &e3 = efld_resist.x3e;
  // local copy: reading the member inside the kernels would capture `this`
  const bool use_sts = use_rkg_sts;

  auto &use_spherical_polar = pmy_pack->pmesh->use_spherical_polar;
  auto &x1v_ = pmy_pack->pcoord->x1v;
  auto &x1f_ = pmy_pack->pcoord->xx1f;
  auto &x2v_ = pmy_pack->pcoord->x2v;
  auto &x2f_ = pmy_pack->pcoord->xx2f;
  auto &x3v_ = pmy_pack->pcoord->x3v;
  auto &x3f_ = pmy_pack->pcoord->xx3f;

  // CUBED SPHERE. The Poynting flux is a cross product, so both operands have to be in
  // the SAME frame. `bc` is in the orthonormal frame {rhat, e_xi, (e_eta - c e_xi)/s}
  // (GnomonicEquiangleRaiseVelMHD) while the edge EMFs are the COVARIANT components
  // E.that along each edge, so the eta slot must be rotated, E_orth3 = (E3 - c E2)/s,
  // before the cross product and the result projected onto the face NORMAL afterwards.
  // Only the x1 and x2 faces need work: nhat_eta IS the third orthonormal axis, so the
  // x3 flux below is already (ExB).nhat_eta as written.
  const bool cs_ = pmy_pack->pmesh->use_cubed_sphere;
  auto &ccell = pmy_pack->pcoord->cos_cell;
  auto &scell = pmy_pack->pcoord->sin_cell;
  auto &cfxi = pmy_pack->pcoord->cos_face_xi;
  auto &sfxi = pmy_pack->pcoord->sin_face_xi;

  //------------------------------
  // energy fluxes in x1-direction
  auto &flx1 = flx.x1f;
  par_for("ohm_heat1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

    Real e3a, e2a;
    if (one_d) {
      e3a = e3(m,k,j,i);
    } else {
      Real fl2 = 0.5;
      Real fr2 = 0.5;
      if (use_spherical_polar) {
        Real dx2tot = x2f_(m,j+1)-x2f_(m,j);
        fl2 = (x2f_(m,j+1)-x2v_(m,j))/dx2tot;
        fr2 = 1.0 - fl2;
      }
      e3a = fl2*e3(m,k,j,i) + fr2*e3(m,k,j+1,i);
    }
    if (three_d) {
      Real fl3 = 0.5;
      Real fr3 = 0.5;
      if (use_spherical_polar) {
        Real dx3tot = x3f_(m,k+1)-x3f_(m,k);
        fl3 = (x3f_(m,k+1)-x3v_(m,k))/dx3tot;
        fr3 = 1.0 - fl3;
      }
      e2a = fl3*e2(m,k,j,i) + fr3*e2(m,k+1,j,i);
    } else {
      e2a = e2(m,k,j,i);
    }
    Real fr1 = 0.5;
    Real fl1 = 0.5;
    if (use_spherical_polar) {
      Real dx1tot = x1v_(m,i)-x1v_(m,i-1);
      fl1 = (x1v_(m,i)-x1f_(m,i))/dx1tot;
      fr1 = 1.0 - fl1;
    }
      
    // flx1 = (E X B)_{1} =  ((\eta J) X B)_{1} = \eta (J2*B3 - J3*B2)
    if (cs_) { e3a = (e3a - ccell(m,k,j)*e2a)/scell(m,k,j); }
    if (use_sts) {
      flx1(m,IEN,k,j,i) = -(fr1*bc(m,IBY,k,j,i) + fl1*bc(m,IBY,k,j,i-1))*e3a
        + (fr1*bc(m,IBZ,k,j,i) + fl1*bc(m,IBZ,k,j,i-1))*e2a;
    } else {
      flx1(m,IEN,k,j,i) += -(fr1*bc(m,IBY,k,j,i) + fl1*bc(m,IBY,k,j,i-1))*e3a
        + (fr1*bc(m,IBZ,k,j,i) + fl1*bc(m,IBZ,k,j,i-1))*e2a;
    }
  });
  if (one_d) {return;}

  //------------------------------
  // energy fluxes in x2-direction
  auto &flx2 = flx.x2f;
  par_for("ohm_heat2", DevExeSpace(), 0, nmb1, ks, ke, js, je+1, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e1a;
    Real fl1 = 0.5;
    Real fr1 = 0.5;
    if (use_spherical_polar) {
      Real dx1tot = x1f_(m,i+1)-x1f_(m,i);
      fl1 = (x1f_(m,i+1)-x1v_(m,i))/dx1tot;
      fr1 = 1.0 - fl1;
    }
    Real e3a = fl1*e3(m,k,j,i) + fr1*e3(m,k,j,i+1);
    if (three_d) {
      Real fl3 = 0.5;
      Real fr3 = 0.5;
      if (use_spherical_polar) {
        Real dx3tot = x3f_(m,k+1)-x3f_(m,k);
        fl3 = (x3f_(m,k+1)-x3v_(m,k))/dx3tot;
        fr3 = 1.0 - fl3;
      }
      e1a = fl3*e1(m,k,j,i) + fr3*e1(m,k+1,j,i);
    } else {
      e1a = e1(m,k,j,i);
    }
    Real fr2 = 0.5;
    Real fl2 = 0.5;
    if (use_spherical_polar) {
      Real dx2tot = x2v_(m,j)-x2v_(m,j-1);
      fl2 = (x2v_(m,j)-x2f_(m,j))/dx2tot;
      fr2 = 1.0 - fl2;
    }
      
    // E2 = \eta (J X B)_{2} = \eta (J3*B1 - J1*B3)
    Real fx2;
    if (cs_) {
      // the eight x2 edges around this face average to its centre
      const Real e2a = 0.125*(e2(m,k,j-1,i) + e2(m,k,j-1,i+1) + e2(m,k,j,i)
                            + e2(m,k,j,i+1) + e2(m,k+1,j-1,i) + e2(m,k+1,j-1,i+1)
                            + e2(m,k+1,j,i) + e2(m,k+1,j,i+1));
      const Real c = cfxi(m,k,j), sn = sfxi(m,k,j);
      const Real eo3 = (e3a - c*e2a)/sn;
      const Real b1 = fr2*bc(m,IBX,k,j,i) + fl2*bc(m,IBX,k,j-1,i);
      const Real b2 = fr2*bc(m,IBY,k,j,i) + fl2*bc(m,IBY,k,j-1,i);
      const Real b3 = fr2*bc(m,IBZ,k,j,i) + fl2*bc(m,IBZ,k,j-1,i);
      // (ExB).nhat_xi = s (ExB)_2 - c (ExB)_3 in the orthonormal frame
      fx2 = sn*(eo3*b1 - e1a*b3) - c*(e1a*b2 - e2a*b1);
    } else {
      fx2 = -(fr2*bc(m,IBZ,k,j,i) + fl2*bc(m,IBZ,k,j-1,i))*e1a
        + (fr2*bc(m,IBX,k,j,i) + fl2*bc(m,IBX,k,j-1,i))*e3a;
    }
    if (use_sts) {
      flx2(m,IEN,k,j,i) = fx2;
    } else {
      flx2(m,IEN,k,j,i) += fx2;
    }
  });
  if (two_d) {return;}

  //------------------------------
  // energy fluxes in x3-direction
  auto &flx3 = flx.x3f;
  par_for("ohm_heat3", DevExeSpace(), 0, nmb1, ks, ke+1, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real fl1 = 0.5;
    Real fr1 = 0.5;
    if (use_spherical_polar) {
      Real dx1tot = x1f_(m,i+1)-x1f_(m,i);
      fl1 = (x1f_(m,i+1)-x1v_(m,i))/dx1tot;
      fr1 = 1.0 - fl1;
    }
    Real e2a = fl1*e2(m,k,j,i) + fr1*e2(m,k,j,i+1);
    Real fl2 = 0.5;
    Real fr2 = 0.5;
    if (use_spherical_polar) {
      Real dx2tot = x2f_(m,j+1)-x2f_(m,j);
      fl2 = (x2f_(m,j+1)-x2v_(m,j))/dx2tot;
      fr2 = 1.0 - fl2;
    }
    Real e1a = fl2*e1(m,k,j,i) + fr2*e1(m,k,j+1,i);
    Real fr3 = 0.5;
    Real fl3 = 0.5;
    if (use_spherical_polar) {
      Real dx3tot = x3v_(m,k)-x3v_(m,k-1);
      fl3 = (x3v_(m,k)-x3f_(m,k))/dx3tot;
      fr3 = 1.0 - fl3;
    }

    // E2 = \eta (J X B)_{2} = \eta (J1*B2 - J2*B1)
    if (use_sts) {
      flx3(m,IEN,k,j,i) = -(fr3*bc(m,IBX,k,j,i) + fl3*bc(m,IBX,k-1,j,i))*e2a
        + (fr3*bc(m,IBY,k,j,i) + fl3*bc(m,IBY,k-1,j,i))*e1a;
    } else {
      flx3(m,IEN,k,j,i) += -(fr3*bc(m,IBX,k,j,i) + fl3*bc(m,IBX,k-1,j,i))*e2a
        + (fr3*bc(m,IBY,k,j,i) + fl3*bc(m,IBY,k-1,j,i))*e1a;
    }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Resistivity::NewTimeStep()
//! \brief Compute new time step for resistive MHD

void Resistivity::NewTimeStep(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  // resistive timestep on MeshBlock(s) in this pack
//  if (iso_resist_type.compare("constant") == 0 && !pmy_pack->pmesh->use_spherical_polar) {
//    NewTimeStepConstantGridResist(w0,eos_data);
//  } else {
    NewTimeStepGeneralResist(w0,eos_data);
//  }
  return;
}

void Resistivity::NewTimeStepGeneralResist(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  // resistive timestep on MeshBlock(s) in this pack
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;

  Real dt1 = std::numeric_limits<float>::max();
  Real dt2 = std::numeric_limits<float>::max();
  Real dt3 = std::numeric_limits<float>::max();
    
  auto &mbsize = pmy_pack->pmb->mb_size;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;
    
  // Coordinates::dx{1,2,3} are PHYSICAL lengths for both curvilinear systems, while
  // mb_size.dx* are index spacings -- on the cubed sphere x2 and x3 are angles on [-1,1],
  // so the Cartesian branch would overestimate the diffusive dt by orders of magnitude.
  const bool use_curvilinear = pmy_pack->pmesh->use_spherical_polar ||
                               pmy_pack->pmesh->use_cubed_sphere;
  auto &dx1_ = pmy_pack->pcoord->dx1;
  auto &dx2_ = pmy_pack->pcoord->dx2;
  auto &dx3_ = pmy_pack->pcoord->dx3;
  // local copy: reading the member inside the reduction would capture `this`
  auto eta_b_ = eta_b;

  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
  Kokkos::parallel_reduce("MHDdiffdt",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
    // compute m,k,j,i indices of thread and call function
    int m = (idx)/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    Real eta;
//    if (iso_resist_type.compare("constant") == 0) {
//      eta = eta_ohm_const;
//    } else {
      eta = eta_b_(m,k,j,i);
//    }
    if (use_curvilinear) {
      min_dt1 = fmin(fac*SQR(dx1_(m,k,j,i))/eta, min_dt1);
      min_dt2 = fmin(fac*SQR(dx2_(m,k,j,i))/eta, min_dt2);
      min_dt3 = fmin(fac*SQR(dx3_(m,k,j,i))/eta, min_dt3);
    } else {
      min_dt1 = fmin(fac*SQR(mbsize.d_view(m).dx1)/eta, min_dt1);
      min_dt2 = fmin(fac*SQR(mbsize.d_view(m).dx2)/eta, min_dt2);
      min_dt3 = fmin(fac*SQR(mbsize.d_view(m).dx3)/eta, min_dt3);
    }
  }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2),Kokkos::Min<Real>(dt3));
    
  dtnew = dt1;
  if (pmy_pack->pmesh->multi_d) { dtnew = std::min(dtnew, dt2); }
  if (pmy_pack->pmesh->three_d) { dtnew = std::min(dtnew, dt3); }
  return;
}

void Resistivity::NewTimeStepConstantGridResist(const DvceArray5D<Real> &w0, const EOS_Data &eos_data) {
  // resistive timestep on MeshBlock(s) in this pack
  dtnew = std::numeric_limits<float>::max();
  auto size = pmy_pack->pmb->mb_size;
  Real fac;
  if (pmy_pack->pmesh->three_d) {
    fac = 1.0/6.0;
  } else if (pmy_pack->pmesh->two_d) {
    fac = 0.25;
  } else {
    fac = 0.5;
  }
  for (int m=0; m<(pmy_pack->nmb_thispack); ++m) {
    dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx1)/eta_ohm_const);
    if (pmy_pack->pmesh->multi_d) {
      dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx2)/eta_ohm_const);
    }
    if (pmy_pack->pmesh->three_d) {
      dtnew = std::min(dtnew, fac*SQR(size.h_view(m).dx3)/eta_ohm_const);
    }
  }
  return;
}
