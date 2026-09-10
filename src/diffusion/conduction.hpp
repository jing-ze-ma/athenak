#ifndef DIFFUSION_CONDUCTION_HPP_
#define DIFFUSION_CONDUCTION_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file conduction.hpp
//! \brief Contains data and functions that implement various formulations for conduction.
//  Currently only isotropic conduction implemented

#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"

//----------------------------------------------------------------------------------------
//! \class Conduction
//! \brief data and functions that implement thermal conduction in Hydro and MHD

class Conduction {
 public:
  Conduction(std::string block, MeshBlockPack *pp, ParameterInput *pin);
  ~Conduction();

  // data
  Real dtnew;
  // the cell that set dtnew, for the collapse report in Mesh::NewTimeStep
  int dtnew_m = -1, dtnew_k = -1, dtnew_j = -1, dtnew_i = -1;
  // ...and its state, so the report says WHY that cell is slow rather than only where
  // it is.  Filled by one extra single-cell kernel, and only when dtnew has just
  // collapsed (or on the first call), so it costs nothing in a healthy run.
  static constexpr int ndtdiag = 16;
  DualArray1D<Real> dt_diag;
  bool dt_diag_valid = false;
  Real dtnew_prev = -1.0;
  std::string iso_cond_type; // "constant", "spitzer", "spitzer_limited", "radiative"
  // radiative conductivity kappa_rad = 16 sigma T^3/(3 kappa_R rho), kappa_R the
  // Freedman+2014 Rosseland mean (utils/rosseland.hpp), applied on faces whose pressure
  // is at or above rad_pcut (the layers the two-stream RT does not reach), flux-limited
  // to
  // sigma T^4, with an imposed heat flux through the inner x1 wall (the planet's internal
  // flux). rad_kappa_fac scales kappa_R (tests).
  Real rad_met = 0.0;          // [M/H] for the Rosseland fit
  Real rad_pcut = 0.0;         // code units; faces with p < rad_pcut get no flux
  Real rad_flux_inner = 0.0;   // code units; heat flux through the inner x1 wall
  Real rad_kappa_fac = 1.0;
  bool rad_flux_limit = true;
  // OPTICAL-DEPTH BLEND (rad_tau_hi > 0): instead of the pressure cut, each x1 face
  // carries a weight w(tau_R) rising smoothly from 0 at rad_tau_lo to 1 at rad_tau_hi,
  // tau_R the column's Rosseland depth from the top; the diffusion flux is multiplied
  // by w and a two-stream RT that reads rad_w multiplies its own by 1 - w, so the two
  // operators overlap and hand over conservatively.  rad_w and rad_tauf sit on x1 faces.
  Real rad_tau_lo = 0.0, rad_tau_hi = 0.0;
  bool rad_tau_mode = false;
  // rad_blend_radial (default true): the blend weight applies to the x1 faces too.
  // False keeps the radial diffusion at full weight everywhere and applies w only to
  // the angular faces -- for a SELF-LUMINOUS object with no two-stream above the
  // blend, where flux-limited radial diffusion is the surface treatment and the
  // horizontal exchange in the optically thin layers is what is switched off.
  bool rad_blend_radial = true;
  // rad_kappa_src = freedman (default) | table: with table, kappa_R(T,p) is a bilinear
  // lookup of log10 kappa_R over (log10 T, log10 p[cgs]) in a table the problem
  // generator hands over ONCE at start-up (deep_hot_jupiter_rt tabulates the Rosseland
  // mean of its own correlated-k table + continuum on that table's grid, so the
  // diffusion and the two-stream share one opacity). Same cost as the Freedman fit.
  // Until the table is set (rad_kr_nT == 0) the Freedman fit is used.
  bool rad_kappa_tab = false;
  // rad_kappa_src = table_rho: the same lookup with log10 rho [g/cm^3] as the second
  // axis instead of log10 p.  Stellar opacity tables (OPAL/OPLIB, AESOPUS, Ferguson) are
  // tabulated in (T, rho) -- via logR = log rho - 3 log T6 -- and a (T, p) axis would
  // need the equation of state to convert, on the table's grid, once per node.  With
  // this mode the pgen hands over (T, rho) directly and every lookup site already has
  // the face or cell density.  rad_kr_lP then holds log10 rho.
  bool rad_kappa_rho = false;
  // rad_cs_exact (default true): the exact face-normal derivative on the cubed sphere;
  // false drops the metric cross term and the 1/sin(alpha) -- DIAGNOSTIC only
  bool rad_cs_exact = true;
  int rad_kr_nT = 0, rad_kr_nP = 0;
  DvceArray2D<Real> rad_kr_tab;            // (iT, iP) log10 kappa_R [cm^2/g]
  DvceArray1D<Real> rad_kr_lT, rad_kr_lP;  // log10 T [K], log10 p [dyn/cm^2], ascending
  DvceArray4D<Real> rad_w, rad_tauf;
  // rad_w is filled by BuildRadWeights, which runs as a task inside the stage.  Until
  // it has, every weight is zero and the timestep below reads that as "no face carries
  // any diffusive flux" -- so the FIRST cycle would run at the hydro timestep with the
  // conduction operator fully on.  NewTimeStep builds the weights itself if this is
  // still false, which is the case at initialisation.
  bool rad_w_built = false;
  void BuildRadWeights(const DvceArray5D<Real> &w, const EOS_Data &eos);
  // PER-CYCLE DIAGNOSTIC (deep_hot_jupiter_rt's problem/diag_gid).  Off by default, and
  // cond_diag is not allocated until EnableDiag is called, so it costs nothing when off.
  // (m,slot,k,j,i) with slot 0/1/2 = the energy flux this module ADDS to
  // flx1/flx2/flx3(m,IEN,k,j,i) on the x1/x2/x3 face of index i/j/k, slot 3 = the
  // radiative conductivity kappa_rad at the cell exactly as NewTimeStep forms it (code
  // units), slot 4 = the flux-limited effective diffusivity keff on the x1 direction,
  // slot 5 = the cell's own conduction dt candidate, the minimum over the three
  // directions BEFORE the dimensional factor fac and before the CFL number.
  bool diag = false;
  DvceArray5D<Real> cond_diag;
  void EnableDiag(int nmb, int n3, int n2, int n1);
  Real kappa_iso;            // isotropic thermal conductivity
  Real kappa_iso_limit;      // limit to isotropic thermal conductivity

  // functions
  void AddHeatFluxes(const DvceArray5D<Real> &w, const EOS_Data &eos,
                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxConstCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                     DvceFaceFld5D<Real> &f);
  void AddIsotropicHeatFluxRadiative(const DvceArray5D<Real> &w, const EOS_Data &eos,
    DvceFaceFld5D<Real> &flx);
  void AddIsotropicHeatFluxSpitzerCond(const DvceArray5D<Real> &w, const EOS_Data &eos,
                                       DvceFaceFld5D<Real> &f);
  void NewTimeStep(const DvceArray5D<Real> &w, const EOS_Data &eos_data);

 private:
  MeshBlockPack* pmy_pack;
  // "hydro" or "mhd": identifies which physics module owns this Conduction object, so
  // that the general EOS derived-variable arrays (cached temperature and pressure) of the
  // right module can be found.
  std::string my_block;
};
//! \fn Real RosselandTable
//! \brief log-bilinear lookup of a tabulated kappa_R(T,p), indices clamped to the grid
KOKKOS_INLINE_FUNCTION
Real RosselandTable(const DvceArray2D<Real> &tab, const DvceArray1D<Real> &lT,
                    const DvceArray1D<Real> &lP, const int nT, const int nP,
                    const Real tk, const Real pcgs) {
  const Real x = log10(tk), y = log10(pcgs);
  int i = 0, j = 0;
  Real fx = 0.0, fy = 0.0;
  if (!(x > lT(0))) {
    i = 0; fx = 0.0;
  } else if (x >= lT(nT-1)) {
    i = nT-2; fx = 1.0;
  } else {
    while (i < nT-2 && lT(i+1) <= x) ++i;
    fx = (x - lT(i))/(lT(i+1) - lT(i));
  }
  if (!(y > lP(0))) {
    j = 0; fy = 0.0;
  } else if (y >= lP(nP-1)) {
    j = nP-2; fy = 1.0;
  } else {
    while (j < nP-2 && lP(j+1) <= y) ++j;
    fy = (y - lP(j))/(lP(j+1) - lP(j));
  }
  const Real lk = (1.0-fx)*((1.0-fy)*tab(i,j) + fy*tab(i,j+1))
                +      fx *((1.0-fy)*tab(i+1,j) + fy*tab(i+1,j+1));
  return pow(10.0, lk);
}

//! \fn Real RadBlendWeight
//! \brief the diffusion weight of a face: 0 for tau <= lo, 1 for tau >= hi, a raised
//! cosine in log tau between
KOKKOS_INLINE_FUNCTION
Real RadBlendWeight(const Real tau, const Real lo, const Real hi) {
  if (tau <= lo) return 0.0;
  if (tau >= hi) return 1.0;
  const Real s = log(tau/lo)/log(hi/lo);
  return 0.5*(1.0 - cos(M_PI*s));
}

#endif // DIFFUSION_CONDUCTION_HPP_
