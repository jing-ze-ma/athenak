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
  // <problem>/nan_report: record the first x1 face whose radiative heat flux comes out
  // non-finite, with the inputs that made it.  Default false; nothing runs when off.
  bool nan_report = false;
  bool rad_flux_limit = true;
  // OPTICAL-DEPTH BLEND (rad_tau_hi > 0): instead of the pressure cut, each x1 face
  // carries a weight w(tau_R) rising smoothly from 0 at rad_tau_lo to 1 at rad_tau_hi,
  // tau_R the column's Rosseland depth from the top; the diffusion flux is multiplied
  // by w and a two-stream RT that reads rad_w multiplies its own by 1 - w, so the two
  // operators overlap and hand over conservatively.  rad_w and rad_tauf sit on x1 faces.
  Real rad_tau_lo = 0.0, rad_tau_hi = 0.0;
  bool rad_tau_mode = false;
  // rad_tmax_kappa (0 = off): a CEILING on the temperature that enters kappa_rad, both
  // its T^3 and the kappa_R lookup.  kappa_rad ~ T^3/(kappa_R rho) explodes when a
  // single drained cell runs away to 1e5-1e7 K, and the explicit diffusive dt collapses
  // with it -- conduction reports the runaway rather than causing it, and capping T here
  // keeps the timestep survivable while the column refills.  The temperature GRADIENT
  // and the flux limiter's sigma T^4 are left at the true temperature.
  Real rad_tmax = 0.0;         // K; <= 0 disables the cap
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
  // rad_kappa_rmax > 0 [code length]: a RADIATIVELY INERT region above this radius.  The
  // red-giant runs put a hot hydrostatic corona above the star, and that corona must
  // neither cool nor set the timestep: it is a numerical lid, not a stellar layer, and a
  // 6e5 K gas at 1e-24 g/cm^3 has a conduction time short enough to stop the run dead.
  // Above rmax the radiative conduction flux is zero, the cell is dropped from the
  // conduction timestep, and the optical depth integrated down from the top accumulates
  // at rad_kappa_above instead of the table value -- so tau does not grow through the
  // corona and the blend weight w stays 0 there.  The two-stream reads the same two
  // numbers (see two_stream_rt.hpp) so both radiative operators go quiet together.
  Real rad_kappa_rmax = 0.0;
  // rad_kappa_above [cm^2/g]: the opacity used above rad_kappa_rmax.  0 (the default)
  // makes the corona perfectly transparent; the grey chain handles kappa = 0 exactly
  // (dtau = 0 gives e0 = 0, so both the source and the emission vanish identically).
  Real rad_kappa_above = 0.0;
  // rad_gate_rho [g/cm^3] > 0: the DENSITY form of the same idea, and the one to use
  // whenever the artificial medium can exchange gas with the star.  The radius test
  // above is a fixed shell: it makes the corona inert, but it ALSO strips the opacity
  // from stellar gas that inflates or is ejected past that radius -- which is exactly
  // what the open-top red-giant runs do from t ~ 5.7e5 on (mass above 3.6e12 grows by
  // 4.5-14x, the photosphere crosses the cutoff, and the two-stream then radiates 9000 K
  // gas as a bare blackbody at 2-4 L).  Gating on the LOCAL DENSITY instead makes the
  // opacity follow the gas: the 1e-15 corona stays inert wherever it is, and stellar
  // material stays opaque wherever it goes.  The gate is a logistic in log10 rho (see
  // RadGate) rising from 0.1 to 0.9 over rad_gate_dex decades, and the effective opacity
  // is  kappa_eff = G kappa_table + (1 - G) rad_kappa_above.  Mutually exclusive with
  // rad_kappa_rmax (fatal if both are set): one criterion or the other.
  Real rad_gate_rho = 0.0;
  Real rad_gate_dex = 0.5;
  // rad_implicit_x1 (<hydro>/ or <mhd>/rad_implicit_x1, default false): solve the RADIAL
  // radiative diffusion IMPLICITLY (backward Euler) instead of adding it to the
  // face fluxes.  The explicit radial operator has dt ~ dx^2 rho c_v/kappa_rad with
  // kappa_rad ~ T^3/(kappa_R rho): in an evacuated cell just above the photosphere
  // that limit falls below 1e-4 s while the run takes 30 s steps, and because dt is
  // evaluated from the PREVIOUS state the overshoot is unbounded -- the cell runs away
  // to 1e10 K in a single step.  With this on, the interior x1 faces carry no explicit
  // flux; a per-column tridiagonal solve (Conduction::ImplicitRadialUpdate) applies the
  // same flux-limited operator unconditionally stably, and the x1 conduction timestep
  // constraint is dropped.  The whole radial extent must live in ONE MeshBlock (the
  // solve is column-local, no MPI), which the constructor checks.  x2/x3 stay explicit.
  // In MHD the frozen internal energy also has the magnetic energy taken out of it
  // (MagEnergyCC); the operator itself is identical.
  bool rad_implicit_x1 = false;
  // scratch for the tridiagonal solve, allocated once: slots (e*, T*, 1/(rho c_v), p,
  // A_f K_f/dl_f, c', d') per cell/face of every column.  See imp_* below.
  static constexpr int nimpw = 7;
  static constexpr int IMPE = 0, IMPT = 1, IMPA = 2, IMPP = 3;
  static constexpr int IMPC = 4, IMPCP = 5, IMPDP = 6;
  DvceArray5D<Real> imp_wrk;
  DvceArray1D<int> imp_flag;   // 1 element: has the first bad cell been recorded?
  DvceArray1D<Real> imp_rec;   // 8 elements: that cell's identity and state
  int imp_lines = 0;           // lines printed so far by the debug report
  void ImplicitRadialUpdate(DvceArray5D<Real> &u0, const EOS_Data &eos,
                            const Real beta_dt);
  // rad_cap_ang (<hydro>/ or <mhd>/rad_cap_ang, default 0 = off): a CONSERVATIVE cap on
  // the explicit ANGULAR (x2/x3) radiative diffusion.  With the radial direction made
  // implicit the same evacuated-cell runaway simply migrates to the angular faces (the
  // I1 test died at t = 5.76e5 with T = 3e10 K in a cell 60x below its shell median),
  // because the angular operator is still explicit and its limit dt ~ dx^2 rho c_v/kappa
  // is evaluated from the PREVIOUS state.  Define per cell the explicit angular
  // stiffness
  //     x_i = beta_dt alpha_i (sum over its 4 angular faces of A_f K_f/dl_f)/V_i
  //         = dt/dt_cond,angular,
  // with alpha_i = 1/(rho_i c_v,i) and K_f the SAME frozen face conductivity the x2/x3
  // kernels use (RadFaceKappa x flux limiter x tau-blend weight, divided by sin(alpha)
  // on the cubed sphere, which is the coefficient of (T_j - T_i)/dl in the face flux).
  // Every angular face flux is then multiplied by
  //     f_f = min(1, rad_cap_ang/max(x_i, x_j)),
  // one number per face applied to both of its cells: exactly conservative, purely
  // local, and it makes the effective x of every cell <= rad_cap_ang, so at 0.5 the
  // operator is unconditionally stable and monotone (the checkerboard mode's
  // amplification factor stays >= 0).  Requires rad_implicit_x1; drops dt2/dt3.
  Real rad_cap_ang = 0.0;
  // beta_dt of the CURRENT stage.  The angular fluxes are added in Hydro::Fluxes, which
  // runs before the RK update, so the cap has no other way of knowing the step it is
  // capping; Hydro::Fluxes sets this immediately before calling AddHeatFluxes.
  Real stage_beta_dt = 0.0;
  DvceArray4D<Real> cap_x;    // the stiffness x_i, active cells + one angular ghost
  DvceArray4D<Real> cap_c2;   // A_f K_f/dl_f on the x2 faces
  DvceArray4D<Real> cap_c3;   // ... and on the x3 faces
  DvceArray1D<int> cap_cnt;   // 2: cells with x_i > cap, faces actually capped
  DvceArray1D<Real> cap_rec;  // 6: max x_i and the cell that carries it
  int cap_lines = 0;
  Real cap_diag_x = 0.0;      // the largest x_i of the last call (rank-local)
  int cap_diag_over = 0;      // cells with x_i > rad_cap_ang in the last call
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

//! \fn Real KappaTemp
//! \brief the temperature used in the radiative conductivity: the cell temperature,
//! capped at tmax when tmax > 0 (see Conduction::rad_tmax).  Bitwise identity when off.
KOKKOS_INLINE_FUNCTION
Real KappaTemp(const Real tk, const Real tmax) {
  return (tmax > 0.0) ? fmin(tk, tmax) : tk;
}

//! \fn Real RadGate
//! \brief the density gate G(rho): 1 where the gas is dense enough to carry its
//! tabulated opacity, 0 in the artificial low-density medium, a logistic in log10 rho
//! between.  Exactly
//!     G = 1/(1 + exp(-(log10(rho) - log10(rho_gate))/s)),   s = dex/(2 ln 9),
//! so G = 1/2 at the threshold and G runs from 0.1 to 0.9 over exactly `dex` decades
//! (G = 0.9 at u/s = ln 9, and the width in u is 2 s ln 9 = dex).  rho <= 0 gives 0; the
//! exponent is clamped so a 20-decade excursion cannot overflow.
//! ONE definition, used at every opacity site in the conduction module and the
//! two-stream, so the two operators cannot disagree about which gas is radiative.
KOKKOS_INLINE_FUNCTION
Real RadGate(const Real rho_cgs, const Real rho_gate, const Real dex) {
  if (!(rho_gate > 0.0)) return 1.0;      // gate off
  if (!(rho_cgs > 0.0)) return 0.0;
  const Real s = dex/(2.0*log(9.0));
  Real u = (log10(rho_cgs) - log10(rho_gate))/s;
  if (u > 40.0) return 1.0;
  if (u < -40.0) return 0.0;
  return 1.0/(1.0 + exp(-u));
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
