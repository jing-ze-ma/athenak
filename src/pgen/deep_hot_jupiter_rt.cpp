//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file deep_hot_jupiter.cpp
//! \brief Problem generator for the deep hot Jupiter.
//!
//! REFERENCE: Heng, Menou, Phillipps, MNRAS, 413, 2380 (2011); Deitrick, Mendonça, Schroffenegger, Grimm, Tsai, Heng, ApJS, 248, 30 (2020)

// C++ headers
#include <cmath>
#include <iostream> // cout
#include <fstream>  // ifstream, for the correlated-k table
#include <sstream>  // ostringstream, for the photosphere dump
#include <string>
#include <vector>

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/grid_stretch.hpp"  // StretchR / StretchRPoly
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/srcterms.hpp"
#include "utils/random.hpp"
#include "pgen.hpp"
#include "utils/rosseland.hpp"
#include "utils/correlated_k.hpp"
#include "utils/atm_column.hpp"
#include "utils/two_stream_rt.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "pgen_eos_utils.hpp"
#include "diffusion/resistivity.hpp"
#include "coordinates/cubed_sphere.hpp"

#include <Kokkos_Random.hpp>

// The two-stream RT solver, its configuration and scratch arrays, and the small
// geometry/gravity helpers its kernels call, live in utils/two_stream_rt.hpp and
// utils/atm_column.hpp; they were lifted from here verbatim.  These declarations
// keep every use below spelled exactly as it was.
using atm_column::CSCellAngles;
using atm_column::EffGravAt;
using atm_column::GravAccAt;
using atm_column::GravPotAt;
using atm_column::TGuess;
using atm_column::TideAccP;
using atm_column::TideAccR;
using atm_column::TideAccT;
using two_stream_rt::LimitRTSource;
using two_stream_rt::RTSourceLimiterWarn;
using two_stream_rt::ad_dump_file;
using two_stream_rt::get_Tint;
using two_stream_rt::get_albedo;
using two_stream_rt::get_kapr;
using two_stream_rt::get_picket_fence_coeff;
using two_stream_rt::par_reduce_clip3;
using two_stream_rt::par_reduce_clip4;
using two_stream_rt::picket_fence_two_stream_RT;
using two_stream_rt::rt_B_ptr;
using two_stream_rt::rt_Bb_ptr;
using two_stream_rt::rt_Fb_ptr;
using two_stream_rt::rt_Qb_ptr;
using two_stream_rt::rt_Qv_ptr;
using two_stream_rt::rt_T_ptr;
using two_stream_rt::rt_cf_ptr;
using two_stream_rt::rt_ck;
using two_stream_rt::rt_ck_pcut;
using two_stream_rt::rt_de_max;
using two_stream_rt::rt_dump_done;
using two_stream_rt::rt_dump_file;
using two_stream_rt::rt_dump_j;
using two_stream_rt::rt_dump_k;
using two_stream_rt::rt_dump_m;
using two_stream_rt::rt_icut_ptr;
using two_stream_rt::rt_int_at_cut;
using two_stream_rt::rt_kc_ptr;
using two_stream_rt::rt_nchain;
using two_stream_rt::rt_pb_ptr;
using two_stream_rt::rt_semi_implicit;
using two_stream_rt::rt_split;
using two_stream_rt::rt_srclim_warned;
using two_stream_rt::rt_star_teff;
using two_stream_rt::rt_tau_ptr;
using two_stream_rt::rt_xP_ptr;
using two_stream_rt::rt_xT_ptr;

// EOS-aware conversions shared with the other stratified problem generators
using pgen_eos::EintFromP;
using pgen_eos::PresFromEint;
using pgen_eos::TempKelvin;
using pgen_eos::PresTempFromEint;




//----------------------------------------------------------------------------------------
using pgen_eos::DensFromPT;
using pgen_eos::GradAd;

void HydrostaticEquilibrium(Mesh *pm);
void SourceFunc(Mesh *pm, Real bdt);

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt);
void double_gray_two_stream_RT(Mesh *pm, Real bdt);

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real &taulim, const Real &A, const Real &B, const Real (&C)[3], const Real (&D)[3], const Real (&E)[3], const Real (&gamv)[3], const Real &tau, Real &T);
template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr);
template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr);

void DhjPhotosphereDump(ParameterInput *pin, Mesh *pm);

KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td);
KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                    const Real &ap, const bool &grav_pmass, const int &N,
                    const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                      const Real &ap, const bool &grav_pmass, const View1D &Tarr,
                      const View1D &lgparr, const int &N, const Real &zmax,
                      View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);

//----------------------------------------------------------------------------------------
//! \fn GravAccAt / GravPotAt
//! \brief the gravity model. `<problem>/grav` is the SURFACE gravity g0 at r = ap.
//!
//! By default gravity is constant, g(r) = g0, and the potential is g0*(r - ap). That is a
//! thin-shell approximation, and this problem is not thin: the production domain reaches
//! r/ap = 1.52, where a point mass with the same g0 gives 0.43*g0. Holding g constant
//! over that range compresses the upper atmosphere -- at 1e-6 bar the level sits ~5 scale
//! heights too low and H is understated by ~2.7x, which propagates straight into
//! transmission radii and feature amplitudes.
//!
//! `<problem>/grav_point_mass = true` switches to a point mass with the SAME g0 at ap:
//!     g(r)   = g0 (ap/r)^2
//!     phi(r) = g0 ap (1 - ap/r)
//! Both reduce to the constant-g forms as r -> ap, so the potential's zero point and the
//! shallow-atmosphere limit are unchanged and the flag is a no-op at the inner boundary.
//!
//! grav_acc is NEGATIVE throughout this file (it is -g0), hence the sign juggling.
//!
//! Every call site carried a commented-out `- 0.5*SQR(omega*r*sin(theta))` after the
//! potential -- the centrifugal term, never enabled. It is recorded here once rather than
//! repeated nine times; omega*r at the top is 6 cm/s^2 against g of 400-940, so it is a
//! sub-percent correction if it is ever wanted.

//----------------------------------------------------------------------------------------
//! \fn ApplyRStretch
//! \brief map a uniform radial coordinate through whichever radial stretch is active.
//!
//! The Coordinates class already stretches x1v/x1f/dx1, so anything reading those arrays
//! needs nothing. But this file also rebuilds radii directly with CellCenterX/LeftEdgeX
//! when it lays down the initial condition and the gravitational potential, and those
//! come out UNSTRETCHED. Without this the coordinates and the initial condition would sit
//! on different grids -- silently, since both are individually self-consistent.
//!
//! StretchR/StretchRPoly are the same file-scope definitions Coordinates uses, so the two
//! cannot drift apart.

// ApplyRStretch now lives in coordinates/grid_stretch.hpp, shared with cs_test.



//! \fn TotPotAt
//! \brief the gravitational plus centrifugal potential, -Omega^2 (r sin theta)^2 / 2 on
//! top of GravPotAt (Phi increases outward here, and rotation lowers it at the equator).
//! om2 = 0 (problem/rot_potential = false) returns the gravitational potential alone.
KOKKOS_INLINE_FUNCTION
Real TotPotAt(const Real grav_acc, const Real ap, const Real r, const Real z,
              const bool pmass, const Real om2, const Real theta) {
  return GravPotAt(grav_acc, ap, r, z, pmass) - 0.5*om2*SQR(r*sin(theta));
}

//! \fn ZEffFromPot
//! \brief the height at which the gravitational potential alone equals `phi`: the
//! equivalent height of an equipotential surface.  The 1D column p(z), rho(z) integrated
//! in the gravitational potential is a function of Phi_grav; sampling it at
//! z_eff(Phi_tot)
//! makes the initial state (and the ghost column, filled by the same kernel) a barotropic
//! equilibrium of the total potential, balanced in r AND theta.
KOKKOS_INLINE_FUNCTION
Real ZEffFromPot(const Real grav_acc, const Real ap, const Real phi, const bool pmass) {
  const Real g = -grav_acc;
  if (!pmass) return phi/g;
  return ap/(1.0 - phi/(g*ap)) - ap;
}


//----------------------------------------------------------------------------------------
//! \fn TideAccR / TideAccT / TideAccP
//! \brief the host star's tidal acceleration, in the planet-centred corotating frame.
//!
//! `<problem>/stellar_tide = true` adds the term the frame has always been missing. The
//! full Hill acceleration in a frame corotating with the orbit is
//!     a = Omega^2 (3x, 0, -z),
//! with x towards the star and z along the spin axis. The block below already applies the
//! centrifugal term about the PLANET's own axis, Omega^2 (x, y, 0) -- see `oor`. What
//! is left to add is the difference,
//!     a_tide = Omega^2 (2x, -y, -z),
//! i.e. the stellar tidal tensor proper: outward at the substellar AND antistellar points
//! (+2 Omega^2 r), inward at the terminators and the poles (-Omega^2 r). A prolate bulge
//! along the star-planet axis.
//!
//! In the spherical basis, with mu = sin(theta) cos(phi) the substellar direction cosine
//! (the same `mu0` the radiative transfer uses for the stellar beam):
//!     a_r     = Omega^2 r (3 mu^2 - 1)
//!     a_theta = 3 Omega^2 r sin(theta) cos(theta) cos^2(phi)
//!     a_phi   = -3 Omega^2 r sin(theta) sin(phi) cos(phi)
//!
//! WHY IT MATTERS HERE: a_tide/g scales as r^3, so it is set by how tall the domain is.
//! At the old constant-g top (r/ap = 1.52) it is 1.5% of g at the terminator; with
//! grav_point_mass and r/ap = 2.18 it is 4.4% at the top and 2.9% at the 1e-6 bar
//! terminator isobar, and 3 Omega^2 r / g reaches 13% at the top. r_L1 = 4.28 ap, so the
//! domain now reaches half way to L1. A hydrostatic estimate against the run's own
//! measured H(r) puts the 1e-6 bar level 3.2 cells lower at the terminator and 9.3 cells
//! higher at the substellar point -- about 0.15 H on the limb observable.
//!
//! Like the centrifugal term, this is applied as an explicit momentum source and is NOT
//! folded into phicc0/phi0_x1f, so the well-balanced reconstruction does not see it. That
//! is deliberate and matches how `oor` has always been treated: the well-balanced path
//! balances the spherically symmetric radial potential, and the tidal potential is not
//! spherically symmetric. So switching this flag on leaves the initial
//! condition (built hydrostatic WITHOUT the tide) no longer exactly balanced; the
//! atmosphere adjusts over the first few dynamical times, which the initial Rayleigh drag
//! damps.
//!
//! WHERE IT IS AND IS NOT APPLIED. Applied at four places: the three momentum components
//! plus the work term in the corotating-frame block, and the radiative transfer's
//! `tau = kappa p / g` closure over the unresolved column above the domain (three call
//! sites, via EffGravAt). Deliberately NOT applied at three others:
//!   * get_wb_eos_arr / get_init_eos_arr, the 1D hydrostatic column integrators behind
//!     the initial condition and the well-balanced background. They have no angle, and
//!     the background must stay spherically symmetric for the well-balanced path.
//!   * the Maxwell outer-BC ghost clamp, which is MHD-only and uses |g| purely as a
//!     reference weight for a dimensionless ratio; the untided value is the right scale
//!     there, and theta/phi are not in scope.




//----------------------------------------------------------------------------------------
//! \fn EffGravAt
//! \brief inward gravity including the tidal correction, for the radiative transfer's
//! `tau = kappa p / g` closure over the unresolved column above the domain.
//!
//! `grav` here is POSITIVE (the RT sites pass g0, not -g0). The tidal radial acceleration
//! is outward-positive, so it is subtracted. Clamped at 10% of the untided value so a
//! domain taken close to L1 cannot divide by zero; at the production top it is
//! -8.8% substellar and +4.4% at the terminator, nowhere near the clamp.



// The correlated-k OPACITY layer now lives in utils/correlated_k.hpp so that other
// problem generators can use the same bands and tables; it was lifted from here
// verbatim.  These declarations keep every use below spelled exactly as it was.
using correlated_k::CK_DIFFUSIVITY;
using correlated_k::CK_NB;
using correlated_k::CK_NG;
using correlated_k::CK_NPF;
using correlated_k::CK_PF_TMIN;
using correlated_k::CK_PF_TMAX;
using correlated_k::CK_NCIA;
using correlated_k::CK_CIA_NTMAX;
using correlated_k::CK_NRAY;
using correlated_k::CK_NCE;
using correlated_k::ck_nq;
using correlated_k::ck_nT;
using correlated_k::ck_nP;
using correlated_k::ck_lT_ptr;
using correlated_k::ck_lP_ptr;
using correlated_k::ck_lk_ptr;
using correlated_k::ck_gw_ptr;
using correlated_k::ck_wl_ptr;
using correlated_k::ck_swf_ptr;
using correlated_k::ck_pf_ptr;
using correlated_k::ck_pf_lTmin;
using correlated_k::ck_pf_idlT;
using correlated_k::ce_nT;
using correlated_k::ce_nP;
using correlated_k::ce_lT_ptr;
using correlated_k::ce_lP_ptr;
using correlated_k::ce_ptr;
using correlated_k::cia_nT_ptr;
using correlated_k::cia_T_ptr;
using correlated_k::cia_k_ptr;
using correlated_k::ray_x_ptr;
using correlated_k::planck_fraction_below;
using correlated_k::build_planck_fractions;
using correlated_k::ck_planck_frac;
using correlated_k::ck_planck_bands;
using correlated_k::read_ck_table;
using correlated_k::read_ck_continuum;
using correlated_k::ck_tp_index;
using correlated_k::ck_kappa;
using correlated_k::ck_continuum;
using correlated_k::ck_selftest;
using correlated_k::ck_rt_selftest;
using correlated_k::ck_build_rosseland_table;


// problem/bc_outer_maxwell: whether the outer-x1 ghost extrapolation carries the
// divergence of the Maxwell stress.
//
// The extrapolation is hydrostatic, and this term was added so the ghost sees the
// magnetic force as well as gravity. It enters as e0 -= e_i*dM1mag/(rho_i*grav_acc), so
// its size relative to the hydrostatic term is (B^2/dr)/(rho g). At the outermost
// shell of the standard setup (rho = 4.0e-10, dr = 5.6e7, g = 942) that is 0.03 at 3 G,
// 0.33 at 10 G, 3.0 at 30 G and 33 at 100 G. Past O(1) it is no longer a correction: it
// swings the ghost density by its own size and the last two active radial shells go.
//
// It is now applied as an effective gravity inside the hydrostatic solve instead, and
// clamped, so it is bounded at any field strength -- see the use site. DEFAULT ON, since
// in that form it is the better physics; set false to drop the magnetic force entirely.
bool bc_outer_maxwell = true;

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  bool use_etotgrav = false;
  bool use_wellbalance_static = false;
  bool use_wellbalance_dynamic = false;
  const bool use_spherical_polar = pmy_mesh_->use_spherical_polar;
  // The cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction; see
  // CSCellAngles.  Captured by value / as a DualArray so the device lambdas below can
  // read it without touching `this`.
  const bool use_cubed_sphere_ = pmy_mesh_->use_cubed_sphere;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  // problem/photosphere_dump = <file> writes the tau = 2/3 level per band and per
  // g-point at the end of the run, using the run's own correlated-k opacity.  Needs one
  // RT evaluation first (the per-cell tables are allocated lazily), so restart and take a
  // single cycle rather than running it at nlim = 0.
  if (!pin->GetOrAddString("problem", "photosphere_dump", "").empty()) {
    pgen_final_func = DhjPhotosphereDump;
  }
  // read before anything restart-sensitive: the outer BC needs it on restarts too
  bc_outer_maxwell = pin->GetOrAddBoolean("problem","bc_outer_maxwell",true);
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: outer-x1 Maxwell-stress term in the ghost "
              << "extrapolation is " << (bc_outer_maxwell ? "ON" : "off") << std::endl;
  }
  // blocked-band RT scaling harness: sized once, before any RT call
  rt_split = pin->GetOrAddBoolean("problem","rt_split",false);
  rt_ck = pin->GetOrAddBoolean("problem","rt_ck",false);
  rt_star_teff = pin->GetOrAddReal("problem","ck_star_teff",6000.0);
  rt_dump_file = pin->GetOrAddString("problem","ck_dump_file","");
  rt_dump_m = pin->GetOrAddInteger("problem","ck_dump_m",0);
  rt_dump_j = pin->GetOrAddInteger("problem","ck_dump_j",-1);
  rt_dump_k = pin->GetOrAddInteger("problem","ck_dump_k",-1);
  rt_ck_pcut = pin->GetOrAddReal("problem","ck_pcut_bar",10.0);
  rt_de_max = pin->GetOrAddReal("problem","rt_de_max",0.5);
  rt_semi_implicit = pin->GetOrAddBoolean("problem","rt_semi_implicit",true);
  rt_int_at_cut = pin->GetOrAddBoolean("problem","ck_int_at_cut",true);
  ad_dump_file = pin->GetOrAddString("problem","ad_dump_file","");
  if (rt_ck && !rt_split) {
    // The correlated-k solver only exists inside the split path. Without this, rt_ck=true
    // with rt_split=false silently ran the GREY picket fence and looked like it worked.
    rt_split = true;
    if (global_variable::my_rank == 0) {
      std::cout << "deep_hot_jupiter_rt: problem/rt_ck implies problem/rt_split; "
                << "enabling the split RT path" << std::endl;
    }
  }
  ck_nq = pin->GetOrAddInteger("problem","ck_nquad",1);
  if (ck_nq != 1 && ck_nq != 2) {
    std::cout << "### FATAL ERROR in deep_hot_jupiter_rt: problem/ck_nquad must be 1 "
              << "(diffusivity factor) or 2 (Gauss), got " << ck_nq << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (rt_ck && ck_lk_ptr == nullptr) {
    read_ck_table(pin->GetOrAddString("problem","ck_table",
                                      "data/exo_fms_ck/ck/Premixed_1x_g8_11.txt"),
                  rt_ck_pcut);
    build_planck_fractions(rt_ck_pcut);
    read_ck_continuum(pin->GetOrAddString("problem","ck_data_dir",
                                          "data/exo_fms_ck"),
                      pin->GetOrAddString("problem","ck_swflux",
                                          "sw_band_flux_W121_11.txt"),
                      rt_star_teff);
    ck_selftest();
    ck_rt_selftest();
    // The chain set is fixed by the table once correlated-k is on: every (band, g-point)
    // for each of the two Gauss points of the two-stream angular quadrature. Note this is
    // TWICE the 88 usually quoted for an 11 x 8 scheme -- 88 counts band x g, and this
    // kernel also carries the 2-point angular quadrature the picket fence uses.
    rt_nchain = CK_NB*CK_NG*ck_nq;
    if (global_variable::my_rank == 0) {
      std::cout << "  " << CK_NB << " bands x " << CK_NG << " g-points x " << ck_nq
                << " angular point(s) = " << rt_nchain << " column solves per cell";
      if (ck_nq == 1) {
        std::cout << " (diffusivity factor " << CK_DIFFUSIVITY << ")" << std::endl;
      } else {
        std::cout << " (2-point Gauss)" << std::endl;
      }
    }
  }
  user_bcs_func = HydrostaticEquilibrium;
  // NOTE: this function must NOT return early on a restart. Only the initial condition
  // (the "probini" kernel) and the magnetic field at the end are genuinely one-off;
  // everything in between builds BACKGROUND state that lives in memory only and is
  // therefore gone after a restart -- the gravitational potential phicc0/phi0 consumed by
  // etotgrav and by the well-balanced reconstruction, and the well-balanced reference
  // atmosphere u0wb/w0wb/w0facewb. Returning here left phi identically zero, which
  // silently changed the physics at the restart point for any run with etotgrav = true.
  // The two one-off blocks are guarded with `if (!restart)` individually instead.
  if (pmy_mesh_->one_d || pmy_mesh_->two_d) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only works in 3D" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (pin->GetInteger("mesh", "nx1") != pin->GetInteger("meshblock", "nx1")) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "deep hot Jupiter problem generator only allows one meshblock in r direction for the RT to work properly" << std::endl;
    exit(EXIT_FAILURE);
  }

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // the cell's PANEL, which is what turns (x2,x3) into a direction on the cubed sphere
  auto &mbpanel_ = pmbp->pmb->mb_panel;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  int &ng = indcs.ng;
  int n1m1 = indcs.nx1 + 2*ng - 1;
  int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
  int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
    
  Real r0, r1;
  r0 = pmy_mesh_->mesh_size.x1min;
  r1 = pmy_mesh_->mesh_size.x1max;

  // Select either Hydro or MHD
  DvceArray5D<Real> u0_;
  DvceArray5D<Real> w0_;
    
  DvceArray4D<Real> phi0_x1f;
  DvceArray4D<Real> phi0_x2f;
  DvceArray4D<Real> phi0_x3f;
  DvceArray4D<Real> phicc0;
  DvceArray5D<Real> u0wb;
  DvceArray5D<Real> w0wb;
  DvceArray5D<Real> w0facewb_x1f;
  DvceArray5D<Real> w0facewb_x2f;
  DvceArray5D<Real> w0facewb_x3f;
    
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  auto &x2v_ = pmbp->pcoord->x2v;
  auto &x2f_ = pmbp->pcoord->xx2f;
  auto &x3v_ = pmbp->pcoord->x3v;
  auto &x3f_ = pmbp->pcoord->xx3f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &area2 = pmbp->pcoord->area.x2f;
  auto &area3 = pmbp->pcoord->area.x3f;
  auto &dxe1 = pmbp->pcoord->dxedge.x1e;
  auto &dxe2 = pmbp->pcoord->dxedge.x2e;
  auto &dxe3 = pmbp->pcoord->dxedge.x3e;

  Real gamma;
  if (pmbp->phydro != nullptr) {
    u0_ = pmbp->phydro->u0;
    w0_ = pmbp->phydro->w0;
    gamma = pmbp->phydro->peos->eos_data.gamma;
    use_etotgrav = pmbp->phydro->use_etotgrav;
    use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
    phi0_x1f = pmbp->phydro->phi0.x1f;
    phi0_x2f = pmbp->phydro->phi0.x2f;
    phi0_x3f = pmbp->phydro->phi0.x3f;
    phicc0 = pmbp->phydro->phicc0;
    u0wb = pmbp->phydro->u0wb;
    w0wb = pmbp->phydro->w0wb;
    w0facewb_x1f = pmbp->phydro->w0facewb.x1f;
    w0facewb_x2f = pmbp->phydro->w0facewb.x2f;
    w0facewb_x3f = pmbp->phydro->w0facewb.x3f;
  } else if (pmbp->pmhd != nullptr) {
    u0_ = pmbp->pmhd->u0;
    w0_ = pmbp->pmhd->w0;
    gamma = pmbp->pmhd->peos->eos_data.gamma;
    use_etotgrav = pmbp->pmhd->use_etotgrav;
    use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
    use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
    phi0_x1f = pmbp->pmhd->phi0.x1f;
    phi0_x2f = pmbp->pmhd->phi0.x2f;
    phi0_x3f = pmbp->pmhd->phi0.x3f;
    phicc0 = pmbp->pmhd->phicc0;
    u0wb = pmbp->pmhd->u0wb;
    w0wb = pmbp->pmhd->w0wb;
    w0facewb_x1f = pmbp->pmhd->w0facewb.x1f;
    w0facewb_x2f = pmbp->pmhd->w0facewb.x2f;
    w0facewb_x3f = pmbp->pmhd->w0facewb.x3f;
  }
  Real gm1 = gamma - 1.0;
  Real igm1 = 1.0/gm1;
  // by-value copy, capturable in the device lambdas below. Every pressure -> internal
  // energy conversion goes through it: for a general EOS e is not p/(gamma-1), and the
  // well-balanced background arrays built here feed the static scheme directly.
  auto eos = (pmbp->phydro != nullptr) ? pmbp->phydro->peos->eos_data
                                       : pmbp->pmhd->peos->eos_data;
    
//  Real Teq = 1469.0;
//  Real grav_acc = -942.0;
//  Real ap = 9.44e9;
//  Real omega = 2.06e-5;
//  Real Rgas = 4.593e7;
//  Real met = 0.0;
    
    Real Teq = pin->GetReal("problem","Teq");
    Real grav_acc = -pin->GetReal("problem","grav");
    Real ap = pin->GetReal("problem","ap");
    // Gravity model: constant g (default, historical) or a point mass with the same g0 at
    // ap. Needs a RADIAL x1, which spherical polar and the CUBED SPHERE both have -- in
    // the Cartesian branch x1v is a height in a plane-parallel box and there is no r to
    // fall off with.
    const bool grav_pmass = pin->GetOrAddBoolean("problem","grav_point_mass",false);
    if (grav_pmass && !(use_spherical_polar || use_cubed_sphere_)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "problem/grav_point_mass requires a radial x1: "
                << "mesh/use_spherical_polar or mesh/use_cubed_sphere"
                << std::endl;
      exit(EXIT_FAILURE);
    }
    // The stellar tide is NOT extended to the cubed sphere: it writes acceleration
    // components in the spherical (r, theta, phi) basis, and the gnomonic tangent pair is
    // neither that basis nor orthonormal, so the expressions would be silently wrong
    // rather than merely untested. It is deliberately off in production anyway.
    const bool tide = pin->GetOrAddBoolean("problem","stellar_tide",false);
    if (tide && !use_spherical_polar) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "problem/stellar_tide requires mesh/use_spherical_polar"
                << std::endl;
      exit(EXIT_FAILURE);
    }

    // Radial grid stretch, copied out of the Mesh so the device lambdas below capture
    // plain values. This function rebuilds radii with CellCenterX/LeftEdgeX, which are
    // UNSTRETCHED; each is passed through ApplyRStretch so the initial condition and the
    // gravitational potential land on the same grid the Coordinates arrays describe.
    const bool str_r_ = pmy_mesh_->use_grid_stretch_r;
    const bool str_rp_ = pmy_mesh_->use_grid_stretch_r_poly;
    const Real fstr_r_ = pmy_mesh_->fStretchR;
    const Real rmin_ = pmy_mesh_->mesh_size.x1min;
    const Real rmax_ = pmy_mesh_->mesh_size.x1max;
    Real cpoly_[NSTRETCH_R_POLY];
    for (int n=0; n<NSTRETCH_R_POLY; ++n) cpoly_[n] = pmy_mesh_->fStretchRPoly[n];
    Real omega = pin->GetReal("problem","omega");
    // problem/rot_potential: the centrifugal potential joins the gravitational one in
    // the column, the ghosts and the well-balanced arrays (see TotPotAt / ZEffFromPot)
    const bool rot_potential = pin->GetOrAddBoolean("problem","rot_potential",false);
    const Real om2_ = rot_potential ? SQR(omega) : 0.0;
    Real Rgas = pin->GetReal("problem","Rgas");
    Real met = pin->GetReal("problem","met");

    // <problem>/Rgas is the ideal-gas R/mu that fixes this atmosphere's mean molecular
    // weight. Under a general EOS composition lives in the EOS instead, so Rgas is unused
    // and the Kelvin conversion is eos.temp_cgs = (pres_cgs/dens_cgs) m_u/k_B, which
    // carries NO mu. The two therefore agree only at mu = 1, i.e. Rgas*temp_cgs = 1. That
    // is not required for a physical run -- but it is required for the general path to
    // reproduce the ideal one, which is how this pgen is verified, so say so loudly.
    if (eos.IsGeneral() && global_variable::my_rank == 0) {
      Real mu_implied = 1.0/(Rgas*eos.temp_cgs);
      if (fabs(mu_implied - 1.0) > 1.0e-4) {
        std::cout << std::endl << "### WARNING! in " << __FILE__ << " at line "
                  << __LINE__ << std::endl
                  << "<problem>/Rgas = " << Rgas << " implies mean molecular weight "
                  << mu_implied << ", but the general EOS supplies its own composition "
                  << "and is being asked for temperature directly." << std::endl
                  << "Rgas is now unused; this run will NOT reproduce the eos=ideal run. "
                  << "Set Rgas = " << 1.0/eos.temp_cgs << " for that comparison."
                  << std::endl << std::endl;
      }
    }

    // <problem>/met is [M/H] in dex and feeds the opacity fit in get_kapr(); the EOS
    // scales its metal electron donors by eos_metal_mh, which DEFAULTS to met and so
    // normally agrees automatically. It can only differ if it was set explicitly, and the
    // result would be an atmosphere opaque at one metallicity and conducting at another,
    // with the electron fraction -- hence the Ohmic resistivity -- wrong by
    // 10^(met - eos_metal_mh). Refuse rather than let that pass.
    if (eos.IsGeneral() && eos.MetalIonization()) {
      if (fabs(met - eos.MetalMetallicity()) > 1.0e-6) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl
                  << "<problem>/met = " << met << " but the EOS metal donors were built "
                  << "with eos_metal_mh = " << eos.MetalMetallicity()
                  << "; both are [M/H] in dex and must agree. Remove the explicit "
                  << "eos_metal_mh and it will follow met." << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }

    Real iap = 1.0/ap;
    
    Real grav = -grav_acc;
    Real Tirr = Teq*sqrt(2);
    Real Tint;
    get_Tint(Teq, Tint);
    Real mus = 1.0;
    
    const int N = 10000;
    DualArray1D<Real> zarr("zarr", N);
    DualArray1D<Real> logparr("logparr", N);
    get_wb_eos_arr(eos, Rgas, grav_acc, ap, grav_pmass, N, (r1-r0)*1.1,
                   zarr.h_view, logparr.h_view);
//    zarr.template modify<HostMemSpace>();
//    zarr.template sync<DevExeSpace>();
//    logparr.template modify<HostMemSpace>();
//    logparr.template sync<DevExeSpace>();
    zarr.modify_host();
    zarr.sync_device();
    logparr.modify_host();
    logparr.sync_device();
    
    DualArray1D<Real> Tarr_init("Tarrinit", N);
    DualArray1D<Real> lgparr_init("lgparrinit", N);
    get_picket_fence_pT_arr(eos, Rgas, gamma, Tint, Tirr, met, grav, mus, N, Tarr_init.h_view, lgparr_init.h_view);
    
    Tarr_init.modify_host();
    Tarr_init.sync_device();
    lgparr_init.modify_host();
    lgparr_init.sync_device();
    
//    DvceArray1D<Real> zinitarr("zinitarr", N);
//    DvceArray1D<Real> logpinitarr("logpinitarr", N);
    DualArray1D<Real> zarr_init("zarrinit", N);
    DualArray1D<Real> logparr_init("logparrinit", N);
    get_init_eos_arr(eos, Rgas, grav_acc, ap, grav_pmass, Tarr_init.h_view,
                     lgparr_init.h_view, N, (r1-r0)*1.1, zarr_init.h_view,
                     logparr_init.h_view);
    
    zarr_init.modify_host();
    zarr_init.sync_device();
    logparr_init.modify_host();
    logparr_init.sync_device();
  
    // one-off: the initial condition. Skipped on a restart, where u0/w0 come from file.
    if (!restart) {
    par_for("probini", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      int nx1 = indcs.nx1;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      int nx2 = indcs.nx2;
      Real &x3min = size.d_view(m).x3min;
      Real &x3max = size.d_view(m).x3max;
      int nx3 = indcs.nx3;
        
      Real x1v, x2v, x3v;
      if (use_spherical_polar) {
        x1v = x1v_(m,i);
        x2v = x2v_(m,j);
        x3v = x3v_(m,k);
      } else {
        x1v = CellCenterX(i-is, nx1, x1min, x1max);
        ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
        if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
        x2v = CellCenterX(j-js, nx2, x2min, x2max);
        x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      }
      Real r = x1v;
      if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
        
      Real lam, theta, phi;
      if (use_spherical_polar) {
        theta = x2v;
        lam = -x2v+M_PI/2.0;
        phi = x3v-M_PI;
      } else if (use_cubed_sphere_) {
        CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
      } else {
        theta = -(x2v*iap-M_PI/2.0);
        lam = x2v*iap;
        phi = x1v*iap;
      }
        
      Real pwb, denwb;
      // rot_potential: sample the column at the equivalent height of this cell's total
      // potential (barotropic in Phi_tot); otherwise at its geometric height
      const Real zs = (om2_ > 0.0) ?
          ZEffFromPot(grav_acc, ap,
                      TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta),
                      grav_pmass) : x1v;
      get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,zs,denwb,pwb);
        
      Real p, den;
      get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,zs,den,p);
//      get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//      get_init_eos(zinitarr,logpinitarr,x3v,lam,phi,den,p);
//      p = pwb;
//      den = denwb;

      u0_(m,IDN,k,j,i) = den;
      u0_(m,IM1,k,j,i) = 0.0;
      u0_(m,IM2,k,j,i) = 0.0;
      u0_(m,IM3,k,j,i) = 0.0;
      u0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
        
      w0_(m,IDN,k,j,i) = den;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVY,k,j,i) = 0.0;   // was IVX three times; harmless, ConsToPrim rebuilt w0
      w0_(m,IVZ,k,j,i) = 0.0;
      w0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
        
        Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
        if (use_etotgrav) {
            u0_(m,IEN,k,j,i) += den*phicc;
        }
    });
    }  // end of !restart guard on the initial condition

    par_for("probwb", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
        
        Real &x1min = size.d_view(m).x1min;
        Real &x1max = size.d_view(m).x1max;
        int nx1 = indcs.nx1;
        Real &x2min = size.d_view(m).x2min;
        Real &x2max = size.d_view(m).x2max;
        int nx2 = indcs.nx2;
        Real &x3min = size.d_view(m).x3min;
        Real &x3max = size.d_view(m).x3max;
        int nx3 = indcs.nx3;
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          x1v = CellCenterX(i-is, nx1, x1min, x1max);
          ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
          if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        }
        Real r = x1v;
        if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          
        Real lam, theta, phi;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -x2v+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          theta = -(x3v*iap-M_PI/2.0);
          lam = x3v*iap;
          phi = x2v*iap;
        }
          
        Real pwb, denwb;
        const Real zs = (om2_ > 0.0) ?
            ZEffFromPot(grav_acc, ap,
                        TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta),
                        grav_pmass) : x1v;
        get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,zs,denwb,pwb);
        
        Real p, den;
        get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,zs,den,p);
//        get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//        get_init_eos(zinitarr,logpinitarr,x1v,lam,phi,den,p);
//        p = pwb;
//        den = denwb;
        
        Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
        
      if (use_etotgrav || use_wellbalance_dynamic) {
          if (use_spherical_polar) {
            x1v = x1f_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1f_(m,i);   // stored centroid/face, as sp
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
          phi0_x1f(m,k,j,i) = phicc;
          if (i == ie) {
              if (use_spherical_polar) {
                x1v = x1f_(m,i+1);
              } else {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
                if (use_cubed_sphere_) x1v = x1f_(m,i+1);   // stored centroid/face, as sp
              }
              r = x1v;
              if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else if (use_cubed_sphere_) {
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
              phi0_x1f(m,k,j,i+1) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2f_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
          phi0_x2f(m,k,j,i) = phicc;
          if (j == je) {
              if (use_spherical_polar) {
                x2v = x2f_(m,j+1);
              } else {
                x2v = LeftEdgeX(j+1-js, nx2, x2min, x2max);
              }
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else if (use_cubed_sphere_) {
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
              phi0_x2f(m,k,j+1,i) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3f_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else if (use_cubed_sphere_) {
            CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
          phi0_x3f(m,k,j,i) = phicc;
          if (k == ke) {
              if (use_spherical_polar) {
                x3v = x3f_(m,k+1);
              } else {
                x3v = LeftEdgeX(k+1-ks, nx3, x3min, x3max);
              }
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else if (use_cubed_sphere_) {
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
              phi0_x3f(m,k+1,j,i) = phicc;
          }
      }
        if (use_wellbalance_static) {
            if (use_spherical_polar) {
              x1v = x1f_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1f_(m,i);   // stored centroid/face, as sp
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            Real denwb;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x1f(m,IDN,k,j,i) = denwb;
            w0facewb_x1f(m,IM1,k,j,i) = 0.0;
            w0facewb_x1f(m,IM2,k,j,i) = 0.0;
            w0facewb_x1f(m,IM3,k,j,i) = 0.0;
            w0facewb_x1f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (i == ie) {
                if (use_spherical_polar) {
                  x1v = x1f_(m,i+1);
                } else {
                  x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
                  ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
                  if (use_cubed_sphere_) x1v = x1f_(m,i+1);   // as sp
                }
                if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
                get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
                w0facewb_x1f(m,IDN,k,j,i+1) = denwb;
                w0facewb_x1f(m,IM1,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM2,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IM3,k,j,i+1) = 0.0;
                w0facewb_x1f(m,IEN,k,j,i+1) = EintFromP(eos, igm1, denwb, pwb);
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2f_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
              x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x2f(m,IDN,k,j,i) = denwb;
            w0facewb_x2f(m,IM1,k,j,i) = 0.0;
            w0facewb_x2f(m,IM2,k,j,i) = 0.0;
            w0facewb_x2f(m,IM3,k,j,i) = 0.0;
            w0facewb_x2f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (j == je) {
                w0facewb_x2f(m,IDN,k,j+1,i) = denwb;
                w0facewb_x2f(m,IM1,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM2,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IM3,k,j+1,i) = 0.0;
                w0facewb_x2f(m,IEN,k,j+1,i) = EintFromP(eos, igm1, denwb, pwb);
            }
            
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3f_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
            w0facewb_x3f(m,IDN,k,j,i) = denwb;
            w0facewb_x3f(m,IM1,k,j,i) = 0.0;
            w0facewb_x3f(m,IM2,k,j,i) = 0.0;
            w0facewb_x3f(m,IM3,k,j,i) = 0.0;
            w0facewb_x3f(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
            if (k == ke) {
                w0facewb_x3f(m,IDN,k+1,j,i) = denwb;
                w0facewb_x3f(m,IM1,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM2,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IM3,k+1,j,i) = 0.0;
                w0facewb_x3f(m,IEN,k+1,j,i) = EintFromP(eos, igm1, denwb, pwb);
            }
        }
    });
    if (use_etotgrav || use_wellbalance_dynamic) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbgrav", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            
            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            int nx1 = indcs.nx1;
            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            
            Real x1v, x2v, x3v;
            if (use_spherical_polar) {
              x1v = x1v_(m,i);
              x2v = x2v_(m,j);
              x3v = x3v_(m,k);
            } else {
              x1v = CellCenterX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
              if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            Real r = x1v;
            if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
              
            Real lam, theta, phi;
            if (use_spherical_polar) {
              theta = x2v;
              lam = -x2v+M_PI/2.0;
              phi = x3v-M_PI;
            } else if (use_cubed_sphere_) {
              CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
            } else {
              theta = -(x3v*iap-M_PI/2.0);
              lam = x3v*iap;
              phi = x2v*iap;
            }
            
            Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
            phicc0(m,k,j,i) = phicc;
            // The FACE potentials on every face, ghost faces included.  The initial
            // condition above fills only the active faces, but the well-balanced
            // background of the outermost ghost cells (which the x1 sweep reconstructs)
            // reads the ghost faces too; left at zero they made that background diverge.
            Real x1fl, x1fr;
            if (use_spherical_polar || use_cubed_sphere_) {
              x1fl = x1f_(m,i);
              x1fr = x1f_(m,i+1);
            } else {
              x1fl = LeftEdgeX(i-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1fl);
              x1fr = LeftEdgeX(i+1-is, nx1, x1min, x1max);
              ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1fr);
            }
            Real zl = x1fl, zr = x1fr;
            if (use_spherical_polar || use_cubed_sphere_) { zl -= ap; zr -= ap; }
            phi0_x1f(m,k,j,i) = TotPotAt(grav_acc, ap, x1fl, zl, grav_pmass, om2_, theta);
            if (i == n1m1) {
              phi0_x1f(m,k,j,i+1) = TotPotAt(grav_acc, ap, x1fr, zr, grav_pmass, om2_,
                                              theta);
            }
        });
//        par_for("wbgravbc", DevExeSpace(), 0, (pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1,
//        KOKKOS_LAMBDA(int m, int k, int j) {
//          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::reflect) {
//            for (int i=0; i<ng; ++i) {
//              phicc0(m,k,j,ie+i+1) = phicc0(m,k,j,ie-i);
//              u0_(m,IDN,k,j,ie+i+1) = u0_(m,IDN,k,j,ie-i);
//              u0_(m,IEN,k,j,ie+i+1) = u0_(m,IEN,k,j,ie-i);
//              phi0_x1f(m,k,j,ie+i+1) = phicc0(m,k,j,ie);
//            }
//          }
//        });
    }
    if (use_wellbalance_static) {
        int &ng = indcs.ng;
        int n1m1 = indcs.nx1 + 2*ng - 1;
        int n2m1 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng - 1) : 0;
        int n3m1 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng - 1) : 0;
        par_for("wbcc", DevExeSpace(), 0,(pmbp->nmb_thispack-1), 0, n3m1, 0, n2m1, 0, n1m1,
        KOKKOS_LAMBDA(int m, int k, int j, int i) {
            
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          int nx1 = indcs.nx1;
          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          int nx2 = indcs.nx2;
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          int nx3 = indcs.nx3;
            
          Real x1v, x2v, x3v;
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
            if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          Real r = x1v;
          if (use_spherical_polar || use_cubed_sphere_) x1v -= ap;
            
          Real pwb, denwb;
          get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
          if (use_etotgrav) {
              Real theta;
              if (use_spherical_polar) {
                theta = x2v;
              } else if (use_cubed_sphere_) {
                Real lam_, phi_;
                CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam_, phi_);
              } else {
                theta = -(x3v*iap-M_PI/2.0);
              }
              Real phicc = TotPotAt(grav_acc, ap, r, x1v, grav_pmass, om2_, theta);
              u0wb(m,IEN,k,j,i) += denwb*phicc;
          }
          w0wb(m,IDN,k,j,i) = denwb;
          w0wb(m,IM1,k,j,i) = 0.0;
          w0wb(m,IM2,k,j,i) = 0.0;
          w0wb(m,IM3,k,j,i) = 0.0;
          w0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
        });
    }

    // the cubed-sphere tangent-pair trig, for the orthonormal-frame bcc below
    auto &ccell_ic_ = pmbp->pcoord->cos_cell;
    auto &scell_ic_ = pmbp->pcoord->sin_cell;
    // initialize magnetic fields if MHD. One-off like the initial condition above: on a
    // restart b0/bcc0 are read from file.
    if (!restart && pmbp->pmhd != nullptr) {
      // Read magnetic field strength
      Real bbot = pin->GetReal("problem","bbot");
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_b0", DevExeSpace(), 0,(pmbp->nmb_thispack-1),0, n3m1, 0, n2m1, 0, n1m1, //ks,ke,js,je,is,ie,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
          if (use_spherical_polar) {
            Real x1v = x1v_(m,i);
            Real x2v = x2v_(m,j);
            Real x3v = x3v_(m,k);
            Real x1fl = x1f_(m,i);
            Real x2fl = x2f_(m,j);
            Real x3fl = x3f_(m,k);
            Real x1fr = x1f_(m,i+1);
            Real x2fr = x2f_(m,j+1);
            Real x3fr = x3f_(m,k+1);
              
            Real A1 = 0.0;
            Real A2 = 0.0;
            Real A2ip = 0.0;
            Real A2kp = 0.0;
            Real A2ipkp = 0.0;
            Real A3 = 0.5*bbot*r0*sin(x2fl)/SQR(x1fl/r0);
            Real A3ip = 0.5*bbot*r0*sin(x2fl)/SQR(x1fr/r0);
            Real A3jp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fl/r0);
            Real A3ipjp = 0.5*bbot*r0*sin(x2fr)/SQR(x1fr/r0);
//            Real A3 = 0.5*bbot*x1fl*sin(x2fl);
//            Real A3ip = 0.5*bbot*x1fr*sin(x2fl);
//            Real A3jp = 0.5*bbot*x1fl*sin(x2fr);
//            Real A3ipjp = 0.5*bbot*x1fr*sin(x2fr);
              
            b0.x1f(m,k,j,i) = (dxe3(m,k,j+1,i)*A3jp - dxe3(m,k,j,i)*A3)/area1(m,k,j,i) - (dxe2(m,k+1,j,i)*A2kp - dxe2(m,k,j,i)*A2)/area1(m,k,j,i);
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je+1)) {
              b0.x2f(m,k,j,i) = - ((x1fr*(x3fr-x3fl))*A3ip - (x1fl*(x3fr-x3fl))*A3) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0.x2f(m,k,j,i) = - (dxe3(m,k,j,i+1)*A3ip - dxe3(m,k,j,i)*A3)/area2(m,k,j,i);
            }
            b0.x3f(m,k,j,i) = (dxe2(m,k,j,i+1)*A2ip - dxe2(m,k,j,i)*A2)/area3(m,k,j,i);
            Real b0x1fip = (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j,i+1)*A3ip)/area1(m,k,j,i+1) - (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k,j,i+1)*A2ip)/area1(m,k,j,i+1);
            Real b0x2fjp;
            if ((mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::polar && j == js-1) || (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::polar && j == je)) {
              b0x2fjp = - ((x1fr*(x3fr-x3fl))*A3ipjp - (x1fl*(x3fr-x3fl))*A3jp) / (0.5*(SQR(x1fr)-SQR(x1fl))*(x3fr-x3fl));
            } else {
              b0x2fjp = - (dxe3(m,k,j+1,i+1)*A3ipjp - dxe3(m,k,j+1,i)*A3jp)/area2(m,k,j+1,i);
            }
            Real b0x3fkp = (dxe2(m,k+1,j,i+1)*A2ipkp - dxe2(m,k+1,j,i)*A2kp)/area3(m,k+1,j,i);
            if (i==n1m1) b0.x1f(m,k,j,i+1) = b0x1fip;
            if (j==n2m1) b0.x2f(m,k,j+1,i) = b0x2fjp;
            if (k==n3m1) b0.x3f(m,k+1,j,i) = b0x3fkp;
              
            Real lw, rw;
            lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
            rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
            bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0x1fip;
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0x2fjp;
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0x3fkp;
              
            u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i))+SQR(bcc0(m,IBY,k,j,i))+SQR(bcc0(m,IBZ,k,j,i)));
          } else if (use_cubed_sphere_) {
            // CUBED SPHERE.  The same dipole, built the same way: each face value is the
            // circulation of A around that face's own edges divided by the face area,
            // using the very `dxedge` and `area` arrays mhd_ct.cpp consumes.  That makes
            // div B vanish to ROUND-OFF rather than to truncation -- CT preserves the
            // divergence it is handed, so a field projected face by face would keep a
            // small permanent monopole on this grid.  Follows cs_test iprob = 12.
            //
            // A = A_phi phihat with A_phi = 0.5*bbot*r0*sin(theta)/(r/r0)^2, which is the
            // spherical-polar expression above.  In Cartesian components
            // sin(theta)*phihat is (-y,x,0)/r, so no angle is ever formed and the axis
            // is not special:
            //   A = 0.5*bbot*r0/(rf/r0)^2 * (-qhat_y, qhat_x, 0).
            // A.rhat = 0 identically, so the RADIAL edge potential vanishes and only the
            // two tangential edge components are needed.  Unlike spherical polar, where
            // phihat IS the x3 direction, both panel tangents see the field here, so the
            // x2-edge component is not zero.
            const int p = mbpanel_.d_view(m);
            // The component of A on an edge's own UNIT tangent, which is what dxedge*A
            // consumes.  `along_xi` selects the x2 edge (radial face, xi CENTRE, eta
            // face) over the x3 edge (radial face, xi face, eta CENTRE).  x2/x3 are PANEL
            // coordinates on [-1,1], not angles: the angle is pi/4 times the coordinate.
            auto Aedge = [&](const int ii, const int jj, const int kk,
                             const bool along_xi) {
              const Real rf = x1f_(m,ii);
              const Real xi = 0.25*M_PI*(along_xi ? x2v_(m,jj) : x2f_(m,jj));
              const Real et = 0.25*M_PI*(along_xi ? x3f_(m,kk) : x3v_(m,kk));
              Real qh[3], e1[3], e2[3];
              cubed_sphere::PanelToCart(p, xi, et, qh);
              cubed_sphere::PanelTangents(p, xi, et, e1, e2);
              const Real aphi = 0.5*bbot*r0/SQR(rf/r0);
              const Real ax = -aphi*qh[1];
              const Real ay =  aphi*qh[0];
              const Real *t = along_xi ? e1 : e2;
              const Real tn = sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
              return (ax*t[0] + ay*t[1])/tn;
            };
            b0.x1f(m,k,j,i) =
                (dxe3(m,k,j+1,i)*Aedge(i,j+1,k,false) - dxe3(m,k,j,i)*Aedge(i,j,k,false)
               - dxe2(m,k+1,j,i)*Aedge(i,j,k+1,true) + dxe2(m,k,j,i)*Aedge(i,j,k,true))
                /area1(m,k,j,i);
            b0.x2f(m,k,j,i) =
                -(dxe3(m,k,j,i+1)*Aedge(i+1,j,k,false)
                - dxe3(m,k,j,i)*Aedge(i,j,k,false))/area2(m,k,j,i);
            b0.x3f(m,k,j,i) =
                (dxe2(m,k,j,i+1)*Aedge(i+1,j,k,true)
               - dxe2(m,k,j,i)*Aedge(i,j,k,true))/area3(m,k,j,i);
            // the far faces of the last cell in each direction, which no cell owns
            Real b0x1fip = (dxe3(m,k,j+1,i+1)*Aedge(i+1,j+1,k,false)
                          - dxe3(m,k,j,i+1)*Aedge(i+1,j,k,false)
                          - dxe2(m,k+1,j,i+1)*Aedge(i+1,j,k+1,true)
                          + dxe2(m,k,j,i+1)*Aedge(i+1,j,k,true))/area1(m,k,j,i+1);
            Real b0x2fjp = -(dxe3(m,k,j+1,i+1)*Aedge(i+1,j+1,k,false)
                           - dxe3(m,k,j+1,i)*Aedge(i,j+1,k,false))/area2(m,k,j+1,i);
            Real b0x3fkp = (dxe2(m,k+1,j,i+1)*Aedge(i+1,j,k+1,true)
                          - dxe2(m,k+1,j,i)*Aedge(i,j,k+1,true))/area3(m,k+1,j,i);
            if (i==n1m1) b0.x1f(m,k,j,i+1) = b0x1fip;
            if (j==n2m1) b0.x2f(m,k,j+1,i) = b0x2fjp;
            if (k==n3m1) b0.x3f(m,k+1,j,i) = b0x3fkp;

            // bcc must be formed the way ConsToPrim forms it on THIS grid, or the
            // magnetic energy added here and the one subtracted on the first C2P
            // disagree.  On the cubed sphere that is GnomonicEquiangleRaiseVelMHD: the
            // radial component interpolated to the centroid x1v, and the tangential pair
            // rotated into the ORTHONORMAL frame, B.e_xi = (b_xi + c b_eta)/s with the
            // eta slot unchanged.  The plain averages used before are components on the
            // non-orthogonal tangent pair, whose sum of squares differs from |B|^2 by
            // (c^2 (b2^2 + b3^2) + 2 c b2 b3)/s^2 -- a few percent of the internal
            // energy in the top cells, where the gas energy is ~1 erg/cm^3.
            {
              const Real lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
              const Real rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
              bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0x1fip;
            }
            {
              const Real cg = ccell_ic_(m,k,j), sg = scell_ic_(m,k,j);
              const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0x2fjp);
              const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0x3fkp);
              bcc0(m,IBY,k,j,i) = (by_n + cg*bz_n)/sg;
              bcc0(m,IBZ,k,j,i) = bz_n;
            }

            u0_(m,IEN,k,j,i) += 0.5*(SQR(bcc0(m,IBX,k,j,i))
                                   + SQR(bcc0(m,IBY,k,j,i))
                                   + SQR(bcc0(m,IBZ,k,j,i)));
          }
      });

      // GATE ON THE CONSTRUCTION, printed once.  Two numbers, because div B alone is
      // satisfied by a ZERO field -- which is exactly what the cubed sphere silently had
      // before this branch existed.  B is the discrete curl of A on the very edges CT
      // uses, so sum(area*B) over a cell must vanish to ROUND-OFF at t=0, and |B|max must
      // be of order bbot.
      {
        auto b1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0.x1f);
        auto b2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0.x2f);
        auto b3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), b0.x3f);
        auto a1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmbp->pcoord->area.x1f);
        auto a2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmbp->pcoord->area.x2f);
        auto a3 = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                      pmbp->pcoord->area.x3f);
        Real dvb = 0.0, bscale = 0.0;
        for (int m=0; m<pmbp->nmb_thispack; ++m) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) {
              for (int i=is; i<=ie; ++i) {
                Real fl = b1(m,k,j,i+1)*a1(m,k,j,i+1) - b1(m,k,j,i)*a1(m,k,j,i)
                        + b2(m,k,j+1,i)*a2(m,k,j+1,i) - b2(m,k,j,i)*a2(m,k,j,i)
                        + b3(m,k+1,j,i)*a3(m,k+1,j,i) - b3(m,k,j,i)*a3(m,k,j,i);
                Real amx = fmax(a1(m,k,j,i+1), fmax(a2(m,k,j+1,i), a3(m,k+1,j,i)));
                Real bmg = fmax(fabs(b1(m,k,j,i)),
                                fmax(fabs(b2(m,k,j,i)), fabs(b3(m,k,j,i))));
                bscale = fmax(bscale, bmg);
                dvb = fmax(dvb, fabs(fl)/fmax(bmg*amx, 1.0e-300));
              }
            }
          }
        }
#if MPI_PARALLEL_ENABLED
        Real dd_[2] = {dvb, bscale};
        MPI_Allreduce(MPI_IN_PLACE, dd_, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
        dvb = dd_[0];  bscale = dd_[1];
#endif
        if (global_variable::my_rank == 0) {
          std::cout << "deep_hot_jupiter_rt: initial B -- max |sum area*B|"
                    << " / (|B| max area) = " << dvb << ", max |B| = " << bscale
                    << " G (bbot = " << bbot << ")" << std::endl;
        }
      }
    }

    // On a restart u0 and b0 are restored exactly, ghost zones included, but bcc0 is a
    // DERIVED array that the restart file does not carry, so it starts at zero. The
    // inner-x1 user boundary swaps the magnetic energy of its ghost cells by subtracting
    // 0.5*bcc0^2 and adding it back from the current face fields; with bcc0 still zero it
    // subtracts nothing and so double counts the magnetic energy on the very first step.
    // Fill bcc0 from the restored face fields here, which is exactly the state the
    // initialisation above leaves on a fresh start. The interpolation must match the one
    // ConsToPrim uses, hence the split on use_spherical_polar.
    if (restart && pmbp->pmhd != nullptr) {
      auto &b0 = pmbp->pmhd->b0;
      auto &bcc0 = pmbp->pmhd->bcc0;
      par_for("pgen_bcc_restart", DevExeSpace(), 0, (pmbp->nmb_thispack-1),
              0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(int m, int k, int j, int i) {
        if (use_spherical_polar) {
          Real lw, rw;
          lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0.x1f(m,k,j,i+1);
          lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
          rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
          bcc0(m,IBY,k,j,i) = lw*b0.x2f(m,k,j,i) + rw*b0.x2f(m,k,j+1,i);
          lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
          rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
          bcc0(m,IBZ,k,j,i) = lw*b0.x3f(m,k,j,i) + rw*b0.x3f(m,k+1,j,i);
        } else if (use_cubed_sphere_) {
          // as GnomonicEquiangleRaiseVelMHD: centroid-weighted radial, orthonormal frame
          Real lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          Real rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
          bcc0(m,IBX,k,j,i) = lw*b0.x1f(m,k,j,i) + rw*b0.x1f(m,k,j,i+1);
          const Real cg = ccell_ic_(m,k,j), sg = scell_ic_(m,k,j);
          const Real by_n = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          const Real bz_n = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
          bcc0(m,IBY,k,j,i) = (by_n + cg*bz_n)/sg;
          bcc0(m,IBZ,k,j,i) = bz_n;
        } else {
          bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
          bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        }
      });
    }

  // the correlated-k Rosseland table for the radiative diffusion (no-op unless
  // <mhd>/rad_kappa_src = table)
  ck_build_rosseland_table((pmbp->pmhd != nullptr) ? pmbp->pmhd->pcond
                                                  : pmbp->phydro->pcond);
  {
    MeshBlockPack *pmbp_ = pmy_mesh_->pmb_pack;
    Conduction *pc = (pmbp_->pmhd != nullptr) ? pmbp_->pmhd->pcond : pmbp_->phydro->pcond;
    if (pc != nullptr && pc->iso_cond_type.compare("radiative") == 0) {
      // the tau blend hands the deep column over between the diffusion and the
      // correlated-k two-stream; the grey scheme has no cut and would double count
      if (pc->rad_tau_mode && !rt_ck) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "rad_tau_hi > 0 (the optical-depth blend) needs "
                  << "problem/rt_ck = true" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      // rad_flux_inner < 0: the planet's internal flux sigma T_int^4 through the wall
      if (pc->rad_flux_inner < 0.0) {
        Real Tint;
        get_Tint(hot_jupiter_param.Teq, Tint);
        const Real fint = 5.670374419e-5*SQR(SQR(Tint));
        pc->rad_flux_inner = fint/(pmbp_->punit->pressure_cgs()
                                   *pmbp_->punit->velocity_cgs());
        if (global_variable::my_rank == 0) {
          std::cout << "  radiative diffusion: wall flux sigma T_int^4 = " << fint
                    << " erg/s/cm^2 (T_int = " << Tint << " K)" << std::endl;
        }
      }
    }
  }

  return;
}


void HydrostaticEquilibrium(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  int &ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
  int &is = indcs.is;  int &ie  = indcs.ie;
  int &js = indcs.js;  int &je  = indcs.je;
  int &ks = indcs.ks;  int &ke  = indcs.ke;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  int nmb = pm->pmb_pack->nmb_thispack;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &size = pmbp->pmb->mb_size;

    // wtemp holds the temperature ConsToPrim solved for; the outer-x1 ghost extrapolation
    // warm starts its hydrostatic solve from the last active cell's value. General EOS
    // only -- a zero-size View otherwise, captured and never read.
    DvceArray4D<Real> wtemp_;
    DvceArray5D<Real> u0_;
    DvceArray5D<Real> w0_;
    Real gamma;
    bool use_etotgrav = false;
    bool use_wellbalance_dynamic = false;
    
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    DvceArray5D<Real> bcc0;
    DvceArray4D<Real> b0_x1f;
    DvceArray4D<Real> b0_x2f;
    DvceArray4D<Real> b0_x3f;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x2f_ = pmbp->pcoord->xx2f;
    auto &x3v_ = pmbp->pcoord->x3v;
    auto &x3f_ = pmbp->pcoord->xx3f;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &area2 = pmbp->pcoord->area.x2f;
    auto &area3 = pmbp->pcoord->area.x3f;
    auto &volume = pmbp->pcoord->volume;
    auto &z_ov_rE = pmbp->pcoord->z_ov_rE;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;

    EOS_Data eos;
    if (pmbp->phydro != nullptr) {
      u0_ = pmbp->phydro->u0;
      w0_ = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
    } else if (pmbp->pmhd != nullptr) {
      u0_ = pmbp->pmhd->u0;
      w0_ = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      bcc0 = pmbp->pmhd->bcc0;
      b0_x1f = pmbp->pmhd->b0.x1f;
      b0_x2f = pmbp->pmhd->b0.x2f;
      b0_x3f = pmbp->pmhd->b0.x3f;
    }
    Real bbot = pm->pgen->hot_jupiter_param.bbot;
    
    int nvar = u0_.extent_int(1);
    
    Real igm1 = 1.0/(gamma-1.0);
    Real gigm1 = gamma*igm1;
    Real gm1ig = (gamma-1.0)/gamma;
    Real ig = 1.0/gamma;
    
    if (pmbp->pmhd != nullptr) {
      par_for("usrboundaryx1_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
      KOKKOS_LAMBDA(int m, int k, int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = is+i;
              int iex = is-i-1;
//              Real fac2 = -area2(m,k,j,iin)/area2(m,k,j,iex);
//              Real fac2p = -area2(m,k,j+1,iin)/area2(m,k,j+1,iex);
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin+1)/area1(m,k,j,iex);
              b0_x1f(m,k,j,iex) = b0_x1f(m,k,j,iin+1)*fac1;
            }
//            for (int i=0; i<ng; ++i) {
//              int iex = is-i-1;
//              Real div_rest = b0_x1f(m,k,j,iex+1)*area1(m,k,j,iex+1) + (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)-b0_x2f(m,k,j,iex)*area2(m,k,j,iex)) + (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)-b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
//              b0_x1f(m,k,j,iex) = div_rest/area1(m,k,j,iex);
//            }
          }
          if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
            for (int i=0; i<ng; ++i) {
              int iin = ie-i;
              int iex = ie+i+1;
//              Real fac2 = -area2(m,k,j,iin)/area2(m,k,j,iex);
//              Real fac2p = -area2(m,k,j+1,iin)/area2(m,k,j+1,iex);
              Real fac2 = -(SQR(x1f_(m,iin+1))-SQR(x1f_(m,iin)))/(SQR(x1f_(m,iex+1))-SQR(x1f_(m,iex)));
              Real fac2p = fac2;
              Real fac3 = -area3(m,k,j,iin)/area3(m,k,j,iex);
              Real fac3p = -area3(m,k+1,j,iin)/area3(m,k+1,j,iex);
              b0_x2f(m,k,j,iex) = b0_x2f(m,k,j,iin)*fac2;
              if (j == n2-1) {b0_x2f(m,k,j+1,iex) = b0_x2f(m,k,j+1,iin)*fac2p;}
              b0_x3f(m,k,j,iex) = b0_x3f(m,k,j,iin)*fac3;
              if (k == n3-1) {b0_x3f(m,k+1,j,iex) = b0_x3f(m,k+1,j,iin)*fac3p;}
              Real fac1 = area1(m,k,j,iin)/area1(m,k,j,iex+1);
              b0_x1f(m,k,j,iex+1) = b0_x1f(m,k,j,iin)*fac1;
            }
//            for (int i=0; i<ng; ++i) {
//              int iex = ie+i+1;
//              Real div_rest = b0_x1f(m,k,j,iex)*area1(m,k,j,iex) - (b0_x2f(m,k,j+1,iex)*area2(m,k,j+1,iex)-b0_x2f(m,k,j,iex)*area2(m,k,j,iex)) - (b0_x3f(m,k+1,j,iex)*area3(m,k+1,j,iex)-b0_x3f(m,k,j,iex)*area3(m,k,j,iex));
//              b0_x1f(m,k,j,iex+1) = div_rest/area1(m,k,j,iex+1);
//            }
//              for (int i=0; i<ng; ++i) {
//                b0_x1f(m,k,j,ie+i+2) = -b0_x1f(m,k,j,ie-i);
//                b0_x2f(m,k,j,ie+i+1) =  b0_x2f(m,k,j,ie-i);
//                if (j == n2-1) {b0_x2f(m,k,j+1,ie+i+1) = b0_x2f(m,k,j+1,ie-i);}
//                b0_x3f(m,k,j,ie+i+1) =  b0_x3f(m,k,j,ie-i);
//                if (k == n3-1) {b0_x3f(m,k+1,j,ie+i+1) = b0_x3f(m,k+1,j,ie-i);}
//              }
        }
      });
    }
    
//    if (pmbp->phydro != nullptr) {
//      if (use_etotgrav) {
//        pmbp->phydro->RemoveGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//      pmbp->phydro->peos->ConsToPrim(u0_, w0_, false, 0, is-1, 0, (n2-1), 0, (n3-1));
//      if (use_etotgrav) {
//        pmbp->phydro->AddGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//    }
//    else if (pmbp->pmhd != nullptr) {
//      auto b0 = pmbp->pmhd->b0;
//      if (use_etotgrav) {
//        pmbp->pmhd->RemoveGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//      pmbp->pmhd->peos->ConsToPrim(u0_, b0, w0_, bcc0, false, 0, is-1, 0, (n2-1), 0, (n3-1));
//      if (use_etotgrav) {
//        pmbp->pmhd->AddGravEtot(phicc0, u0_, 0, is-1, 0, (n2-1), 0, (n3-1));
//      }
//    }
    if (pmbp->phydro != nullptr) {
      if (use_etotgrav) {
        pmbp->phydro->RemoveGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
      pmbp->phydro->peos->ConsToPrim(u0_, w0_, false, ie, ie, 0, (n2-1), 0, (n3-1));
      if (use_etotgrav) {
        pmbp->phydro->AddGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
    }
    else if (pmbp->pmhd != nullptr) {
      auto b0 = pmbp->pmhd->b0;
      if (use_etotgrav) {
        pmbp->pmhd->RemoveGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
      pmbp->pmhd->peos->ConsToPrim(u0_, b0, w0_, bcc0, false, ie, ie, 0, (n2-1), 0, (n3-1));
      if (use_etotgrav) {
        pmbp->pmhd->AddGravEtot(phicc0, u0_, ie, ie, 0, (n2-1), 0, (n3-1));
      }
    }

//    par_for("usrboundaryx1", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n2-1),
//    KOKKOS_LAMBDA(int m, int n, int k, int j) {
//        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
//          for (int i=0; i<ng; ++i) {
//            u0_(m,n,k,j,ie+i+1) = u0_(m,n,k,j,ie);
//          }
//        }
//    });
    
    // Local copy of the file-scope flag. Reading the global directly from the kernel is
    // a reference to a __host__ variable in device code, which hipcc rejects outright --
    // the switch has to be captured by value like any other host state.
    const bool bc_outer_maxwell_ = bc_outer_maxwell;
    // for the metric-correct ghost kinetic energy below (cubed sphere only)
    const bool cs_bc_ = pm->use_cubed_sphere;
    auto &ccell_bc_ = pmbp->pcoord->cos_cell;
    auto &scell_bc_ = pmbp->pcoord->sin_cell;

    // The cell-centred field in the OUTER-x1 ghost zones is built here, in its own kernel,
    // and not in the extrapolation kernel below.  It used to be computed inline there, in
    // the same launch that reads bcc0 at (k,j+1) and (k+1,j) for the Maxwell stress -- cells
    // owned by OTHER threads of that launch.  Whether a thread saw its neighbour's new value
    // or the previous cycle's stale one depended on how the wavefronts happened to be
    // scheduled, which made the whole run non-reproducible on a GPU: two runs of the same
    // binary on one rank diverged within a single cycle, worst in the outer half of the
    // domain, where the two-stream RT then spread it along each radial column.  On a CPU the
    // ascending j loop always lost the race the same way, so this was deterministic (and
    // deterministically wrong) until the problem first ran on an accelerator.
    //
    // Splitting the write from the read is the whole fix: every bcc0 the stress term reads
    // is now a completed write from a previous kernel.
    //
    // The INNER-x1 ghosts are deliberately NOT hoisted: that branch subtracts the magnetic
    // energy using the OLD bcc0 and adds it back with the new one, so moving the update
    // ahead of it would change what it subtracts.  It has no cross-thread read and was
    // measured deterministic on its own.
    par_for("usrboundaryx1_bcc_outer", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user &&
            pmbp->pmhd != nullptr) {
          for (int i=0; i<ng; ++i) {
            Real lw, rw;
            lw = (x1f_(m,(ie+i+1)+1)-x1v_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
            rw = (x1v_(m,(ie+i+1))-x1f_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
            bcc0(m,IBX,k,j,(ie+i+1)) = lw*b0_x1f(m,k,j,(ie+i+1))
                                     + rw*b0_x1f(m,k,j,(ie+i+1)+1);
            lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
            bcc0(m,IBY,k,j,(ie+i+1)) = lw*b0_x2f(m,k,j,(ie+i+1))
                                     + rw*b0_x2f(m,k,j+1,(ie+i+1));
            // (the cubed-sphere frame rotation is applied after IBZ below)
            lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
            bcc0(m,IBZ,k,j,(ie+i+1)) = lw*b0_x3f(m,k,j,(ie+i+1))
                                     + rw*b0_x3f(m,k+1,j,(ie+i+1));
            if (cs_bc_) {
              // CUBED SPHERE: bcc is stored in the ORTHONORMAL frame that
              // GnomonicEquiangleRaiseVelMHD builds (B.e_xi = (b_xi + c b_eta)/s, eta slot
              // unchanged), and ConsToPrim subtracts 0.5*|bcc|^2 in that frame.  The
              // plain face averages are components on the non-orthogonal tangent pair,
              // whose sum of squares is NOT |B|^2: the two differ by the cross term
              // (c^2 (b2^2 + b3^2) + 2 c b2 b3)/s^2, of order B_tan^2 at a panel corner.
              const Real cg = ccell_bc_(m,k,j), sg = scell_bc_(m,k,j);
              bcc0(m,IBY,k,j,(ie+i+1)) = (bcc0(m,IBY,k,j,(ie+i+1))
                                          + cg*bcc0(m,IBZ,k,j,(ie+i+1)))/sg;
            }
          }
        }
    });

    par_for("usrboundaryx1_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n2-1),
    KOKKOS_LAMBDA(int m, int k, int j) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          Real rho_i = u0_(m,IDN,k,j,is);
//          Real e_i = w0_(m,IEN,k,j,is);
//          Real phi_i = phicc0(m,k,j,is);
//          Real q0_i = log(e_i);
//          Real factor_i = rho_i/e_i*igm1;
          for (int i=0; i<ng; ++i) {
            if (pmbp->pmhd != nullptr) {
              u0_(m,IEN,k,j,is-i-1) -= 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
              Real lw, rw;
              lw = (x1f_(m,(is-i-1)+1)-x1v_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              rw = (x1v_(m,(is-i-1))-x1f_(m,(is-i-1)))/(x1f_(m,(is-i-1)+1)-x1f_(m,(is-i-1)));
              bcc0(m,IBX,k,j,(is-i-1)) = lw*b0_x1f(m,k,j,(is-i-1)) + rw*b0_x1f(m,k,j,(is-i-1)+1);
              lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              bcc0(m,IBY,k,j,(is-i-1)) = lw*b0_x2f(m,k,j,(is-i-1)) + rw*b0_x2f(m,k,j+1,(is-i-1));
              lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              bcc0(m,IBZ,k,j,(is-i-1)) = lw*b0_x3f(m,k,j,(is-i-1)) + rw*b0_x3f(m,k+1,j,(is-i-1));
              if (cs_bc_) {
                // CUBED SPHERE: same orthonormal-frame rotation as the outer ghosts.  THIS
                // ONE MATTERED.  The energy update above subtracts 0.5*|bcc0|^2 with the
                // bcc0 that ConsToPrim left here (orthonormal) and adds it back with the
                // one built here (plain averages): on cs the two differ by the cross term,
                // so the ghost TOTAL energy drifted by that mismatch on every call.  The
                // ghost internal energy, uniform over (j,k) at each depth by construction,
                // had a spread of 1e-3 at rot 6, 8 % at rot 16 and ranged 3e8-2e9 around
                // 1.5e9 by rot 26 of cs_prod_mhd_rot -- bottom cells under a low-pressure
                // ghost drained at 2 km/s into the boundary, the domain lost 1 % of its
                // mass and the run went NaN at rot 32.  Spherical polar has c = 0.
                const Real cg = ccell_bc_(m,k,j), sg = scell_bc_(m,k,j);
                bcc0(m,IBY,k,j,(is-i-1)) = (bcc0(m,IBY,k,j,(is-i-1))
                                            + cg*bcc0(m,IBZ,k,j,(is-i-1)))/sg;
              }
              u0_(m,IEN,k,j,is-i-1) += 0.5*(SQR(bcc0(m,IBX,k,j,is-i-1))+SQR(bcc0(m,IBY,k,j,is-i-1))+SQR(bcc0(m,IBZ,k,j,is-i-1)));
            }
//              Real rho0_ip = u0_(m,IDN,k,j,(is-i-1));
//              u0_(m,IM2,k,j,(is-i-1)) = u0_(m,IM2,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IM3,k,j,(is-i-1)) = u0_(m,IM3,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IM1,k,j,(is-i-1)) = 0.0;//-u0_(m,IM1,k,j,is)/rho_i*rho0_ip;
//              u0_(m,IEN,k,j,(is-i-1)) = w0_(m,IEN,k,j,(is-i-1)) + 0.5*(SQR(u0_(m,IM1,k,j,(is-i-1)))+SQR(u0_(m,IM2,k,j,(is-i-1)))+SQR(u0_(m,IM3,k,j,(is-i-1))))/rho0_ip;
////            Real dphi_i = phicc0(m,k,j,(is-i-1))-phi_i;
////            Real q0_ip = q0_i - factor_i * dphi_i;
////            Real e0_ip = exp(q0_ip);
////            if (e0_ip < 0.0) e0_ip = e_i;
////            Real rho0_ip = e0_ip/e_i*rho_i;
////            u0_(m,IDN,k,j,(is-i-1)) = rho0_ip;
////            u0_(m,IM2,k,j,(is-i-1)) = u0_(m,IM2,k,j,is)/rho_i*rho0_ip;
////            u0_(m,IM3,k,j,(is-i-1)) = u0_(m,IM3,k,j,is)/rho_i*rho0_ip;
////            Real mom = u0_(m,IM1,k,j,is)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,is)/rho_i*rho0_ip); //
////            u0_(m,IM1,k,j,(is-i-1)) = mom;
////            u0_(m,IEN,k,j,(is-i-1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(is-i-1)))+SQR(u0_(m,IM2,k,j,(is-i-1)))+SQR(u0_(m,IM3,k,j,(is-i-1))))/rho0_ip;
//            if (use_etotgrav) u0_(m,IEN,k,j,(is-i-1)) += rho0_ip*phicc0(m,k,j,(is-i-1));
//            if (pmbp->pmhd != nullptr) u0_(m,IEN,k,j,(is-i-1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(is-i-1)))+SQR(bcc0(m,IBY,k,j,(is-i-1)))+SQR(bcc0(m,IBZ,k,j,(is-i-1))));
          }
        }
        if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
          Real rho_i = w0_(m,IDN,k,j,ie);
//          Real e_i = u0_(m,IEN,k,j,ie) - 0.5*(SQR(u0_(m,IM1,k,j,ie))+SQR(u0_(m,IM2,k,j,ie))+SQR(u0_(m,IM3,k,j,ie)))/rho_i;
//          if (use_etotgrav) e_i -= rho_i*phicc0(m,k,j,ie);
//          if (pmbp->pmhd != nullptr) e_i -= 0.5*(SQR(bcc0(m,IBX,k,j,ie))+SQR(bcc0(m,IBY,k,j,ie))+SQR(bcc0(m,IBZ,k,j,ie)));
          Real e_i = w0_(m,IEN,k,j,ie);
          Real phi_i = phicc0(m,k,j,ie);
          Real q0_i = log(e_i);
          Real factor_i = rho_i/e_i*igm1;
          for (int i=0; i<ng; ++i) {
            Real dM1mag = 0.0;
            if (pmbp->pmhd != nullptr) {
              // bcc0 in these ghost cells was built by usrboundaryx1_bcc_outer above.
              if (bc_outer_maxwell_) {
              // (k,j+1) and (k+1,j) are read here. They belong to other threads, which is
              // why the write had to move to its own kernel -- see the note on that kernel.
              // At the outermost ghost row j+1 and k+1 leave the cell-centred array (the
              // FACE arrays have the extra slot, the cell-centred one does not), so the
              // index is held at the edge: the correction degenerates to a one-sided
              // difference in that row rather than reading past the end, which is what the
              // previous version did.
              int jp1 = (j+1 < n2) ? (j+1) : j;
              int kp1 = (k+1 < n3) ? (k+1) : k;
              // Same one past the end in x1: on the LAST ghost cell (i = ng-1) the cell
              // centred index ie+i+2 is n1, one beyond the array. The face array b0_x1f has
              // the extra slot and is indexed as before; only bcc0 is held back.
              int ip2 = (ie+i+2 < n1) ? (ie+i+2) : (n1-1);
              Real pb = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
              Real pbp1 = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+2)))+SQR(bcc0(m,IBY,k,j,ip2))+SQR(bcc0(m,IBZ,k,j,ip2)));
              Real M11 = pb - SQR(b0_x1f(m,k,j,(ie+i+1)));
              Real M11p1 = pbp1 - SQR(b0_x1f(m,k,j,(ie+i+2)));
              Real M12 = - b0_x2f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M12p1 = - b0_x2f(m,k,j+1,(ie+i+1)) * bcc0(m,IBX,k,jp1,(ie+i+1));
              Real M13 = - b0_x3f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M13p1 = - b0_x3f(m,k+1,j,(ie+i+1)) * bcc0(m,IBX,kp1,j,(ie+i+1));
              dM1mag = -( (M11p1*area1(m,k,j,(ie+i+2))-M11*area1(m,k,j,(ie+i+1))) + (M12p1*area2(m,k,j+1,(ie+i+1))-M12*area2(m,k,j,(ie+i+1))) + (M13p1*area3(m,k+1,j,(ie+i+1))-M13*area3(m,k,j,(ie+i+1))) )/volume(m,k,j,(ie+i+1));
              dM1mag += z_ov_rE(m,k,j,(ie+i+1)) * 0.5*SQR(bcc0(m,IBX,k,j,(ie+i+1)));
              }
            }
            // Hydrostatic extrapolation into the ghost zone, at fixed temperature. The
            // closed form below is the ideal-gas isothermal background; for a general EOS
            // the same statement is WBAdvance's isothermal branch, which integrates
            // dln d/dPhi = -d/(p chi_rho) instead of assuming p = (gamma-1)e.
            Real dphi_i = phicc0(m,k,j,(ie+i+1))-phi_i;
            // The magnetic force enters as an EFFECTIVE GRAVITY, not as an additive
            // shift of the extrapolated energy. dM1mag/rho is an acceleration, so
            // dphi -> dphi*(1 - dM1mag/(rho |g|)) is the same to first order but stays
            // inside the exponential (and inside WBAdvance's EOS-consistent
            // integration), so a magnetic force comparable to gravity changes the SCALE
            // HEIGHT instead of swinging the answer linearly through zero.
            //
            // The clamp is what makes it usable at high field. gmag -> 1 is the
            // force-free limit where the atmosphere stops falling off, and gmag > 1 is
            // net outward, which an outward-decaying ghost cannot represent; gmag < -1
            // would compress the ghost without bound. Both ends are held back, so the
            // ghost degrades to "very extended" rather than to nonsense.
            if (bc_outer_maxwell_) {
              // local gravity, so the clamp means the same fraction of the real weight
              Real gmag = dM1mag/(rho_i*fabs(GravAccAt(grav_acc, ap,
                                                       x1v_(m,ie+i+1), grav_pmass)));
              gmag = fmin(fmax(gmag, -1.0), 0.9);
              dphi_i *= (1.0 - gmag);
            }
            Real e0_hyd, rho0_hyd;
            if (eos.IsGeneral()) {
              rho0_hyd = rho_i;
              e0_hyd = e_i;
              Real t_hyd = -1.0;   // WBAdvance's temperature hand-off; unused here
              // Warm start from cell ie's own temperature, which the ConsToPrim call on
              // this column above has already solved for and left in wtemp. This runs per
              // ghost cell per stage, and the isothermal branch inverts twice.
              WBAdvance(eos, 1, rho_i, e_i, dphi_i, rho0_hyd, e0_hyd, t_hyd,
                        TGuess(wtemp_, m, k, j, ie));
            } else {
              e0_hyd = exp(q0_i - factor_i * dphi_i);
              rho0_hyd = e0_hyd/e_i*rho_i;
            }
            Real e0_ip = e0_hyd;
            if (e0_ip < 0.0) e0_ip = e_i;
            // density and energy now come from the SAME hydrostatic solve, so they are
            // thermodynamically consistent by construction; the old code had to rescale
            // the density by the energy's relative shift because the magnetic term was
            // bolted on afterwards
            Real rho0_ip = rho0_hyd;
//            rho0_ip = rho_i;
//            e0_ip = e_i;
            u0_(m,IDN,k,j,(ie+i+1)) = rho0_ip;
            u0_(m,IM2,k,j,(ie+i+1)) = u0_(m,IM2,k,j,ie)/rho_i*rho0_ip;
            u0_(m,IM3,k,j,(ie+i+1)) = u0_(m,IM3,k,j,ie)/rho_i*rho0_ip;
            Real mom = u0_(m,IM1,k,j,ie)/rho_i*rho0_ip; // fmax(0.0,u0_(m,IM1,k,j,ie)/rho_i*rho0_ip); // 
            u0_(m,IM1,k,j,(ie+i+1)) = mom;
            // Kinetic energy of the ghost cell.  On the CUBED SPHERE IM2/IM3 are
            // COVARIANT components on a non-orthogonal tangent basis, so
            // KE = 0.5 m_i G^{ij} m_j / rho with G^{-1} = [[1,-c],[-c,1]]/(1-c^2) --
            // the contraction GnomonicEquiangleRaiseVelMHD uses.  The plain sum of
            // squares (kept, bit for bit, for every other grid) drops the -2c m2 m3
            // cross term, up to 50 % of the tangential KE at a panel corner where
            // |c| = 1/2, and ConsToPrim then reads that error back out of E as
            // internal energy.  Same bug class as the rotation source that used to
            // fall through into the beta-plane.
            Real ke_g;
            {
              const Real m1g = u0_(m,IM1,k,j,(ie+i+1));
              const Real m2g = u0_(m,IM2,k,j,(ie+i+1));
              const Real m3g = u0_(m,IM3,k,j,(ie+i+1));
              if (cs_bc_) {
                const Real cg = ccell_bc_(m,k,j);
                ke_g = 0.5*(m1g*m1g + (m2g*m2g - 2.0*cg*m2g*m3g + m3g*m3g)/(1.0 - cg*cg))
                       /rho0_ip;
              } else {
                ke_g = 0.5*(SQR(m1g)+SQR(m2g)+SQR(m3g))/rho0_ip;
              }
            }
            u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + ke_g;
            if (use_etotgrav) u0_(m,IEN,k,j,(ie+i+1)) += rho0_ip*phicc0(m,k,j,(ie+i+1));
            if (pmbp->pmhd != nullptr) u0_(m,IEN,k,j,(ie+i+1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
          }
        }
    });

////    par_for("usrboundaryx2", DevExeSpace(), 0,(nmb-1),0,(nvar-1),0,(n3-1),0,(n1-1),
////    KOKKOS_LAMBDA(int m, int n, int k, int i) {
////        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
////          for (int j=0; j<ng; ++j) {
////            if (n==(IVY)) {
////              u0_(m,n,k,js-j-1,i) = -u0_(m,n,k,js+j,i);
////            } else {
////              u0_(m,n,k,js-j-1,i) = u0_(m,n,k,js+j,i);
////            }
////          }
////        }
////        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
////          for (int j=0; j<ng; ++j) {
////            if (n==(IVY)) {
////              u0_(m,n,k,je+j+1,i) = -u0_(m,n,k,je-j,i);
////            } else {
////              u0_(m,n,k,je+j+1,i) = u0_(m,n,k,je-j,i);
////            }
////          }
////        }
////    });
//    par_for("usrboundaryx2", DevExeSpace(), 0,(nmb-1),0,(n3-1),0,(n1-1),
//    KOKKOS_LAMBDA(int m, int k, int i) {
//        if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//          for (int j=0; j<ng; ++j) {
//            u0_(m,IDN,k,js-j-1,i) = u0_(m,IDN,k,js,i);
//            u0_(m,IM1,k,js-j-1,i) = u0_(m,IM1,k,js,i);
//            u0_(m,IM3,k,js-j-1,i) = u0_(m,IM3,k,js,i);
//            u0_(m,IM2,k,js-j-1,i) = fmin(0.0,u0_(m,IM2,k,js,i));
//            u0_(m,IEN,k,js-j-1,i) = u0_(m,IEN,k,js,i)-0.5*SQR(u0_(m,IM2,k,js,i))+0.5*SQR(u0_(m,IM2,k,js-j-1,i));
//          }
//        }
//        if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//          for (int j=0; j<ng; ++j) {
//            u0_(m,IDN,k,je+j+1,i) = u0_(m,IDN,k,je,i);
//            u0_(m,IM1,k,je+j+1,i) = u0_(m,IM1,k,je,i);
//            u0_(m,IM3,k,je+j+1,i) = u0_(m,IM3,k,je,i);
//            u0_(m,IM2,k,je+j+1,i) = fmax(0.0,u0_(m,IM2,k,je,i));
//            u0_(m,IEN,k,je+j+1,i) = u0_(m,IEN,k,je,i)-0.5*SQR(u0_(m,IM2,k,je,i))+0.5*SQR(u0_(m,IM2,k,je+j+1,i));
//          }
//        }
//    });
//    if (pmbp->pmhd != nullptr) {
//      par_for("usrboundaryx2_bfield", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n1-1),
//      KOKKOS_LAMBDA(int m, int k, int i) {
//          if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
//                int jin = js+j;
//                int jex = js-j-1;
//                Real fac1 = -area1(m,k,jin,i)/area1(m,k,jex,i);
//                Real fac1p = -area1(m,k,jin,i+1)/area1(m,k,jex,i+1);
//                Real fac3 = -area3(m,k,jin,i)/area3(m,k,jex,i);
//                Real fac3p = -area3(m,k+1,jin,i)/area3(m,k+1,jex,i);
//                b0_x1f(m,k,jex,i) = b0_x1f(m,k,jin,i)*fac1;
//                if (i == n1-1) {b0_x1f(m,k,jex,i+1) = b0_x1f(m,k,jin,i+1)*fac1p;}
//                b0_x3f(m,k,jex,i) = b0_x3f(m,k,jin,i)*fac3;
//                if (k == n3-1) {b0_x3f(m,k+1,jex,i) = b0_x3f(m,k+1,jin,i)*fac3p;}
//                Real fac2 = area2(m,k,jin+1,i)/area2(m,k,jex,i);
//                b0_x2f(m,k,jex,i) = b0_x2f(m,k,jin+1,i)*fac2;
//              }
////              for (int j=0; j<ng; ++j) {
////                int jex = js-j-1;
////                Real div_rest = b0_x2f(m,k,jex+1,i)*area2(m,k,jex+1,i) + (b0_x1f(m,k,jex,i+1)*area1(m,k,jex,i+1)-b0_x1f(m,k,jex,i)*area1(m,k,jex,i)) + (b0_x3f(m,k+1,jex,i)*area3(m,k+1,jex,i)-b0_x3f(m,k,jex,i)*area3(m,k,jex,i));
////                b0_x2f(m,k,jex,i) = div_rest/area2(m,k,jex,i);
////              }
//          }
//          if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
//                int jin = je-j;
//                int jex = je+j+1;
//                Real fac1 = -area1(m,k,jin,i)/area1(m,k,jex,i);
//                Real fac1p = -area1(m,k,jin,i+1)/area1(m,k,jex,i+1);
//                Real fac3 = -area3(m,k,jin,i)/area3(m,k,jex,i);
//                Real fac3p = -area3(m,k+1,jin,i)/area3(m,k+1,jex,i);
//                b0_x1f(m,k,jex,i) = b0_x1f(m,k,jin,i)*fac1;
//                if (i == n1-1) {b0_x1f(m,k,jex,i+1) = b0_x1f(m,k,jin,i+1)*fac1p;}
//                b0_x3f(m,k,jex,i) = b0_x3f(m,k,jin,i)*fac3;
//                if (k == n3-1) {b0_x3f(m,k+1,jex,i) = b0_x3f(m,k+1,jin,i)*fac3p;}
//                Real fac2 = area2(m,k,jin,i)/area2(m,k,jex+1,i);
//                b0_x2f(m,k,jex+1,i) = b0_x2f(m,k,jin,i)*fac2;
//              }
////              for (int j=0; j<ng; ++j) {
////                int jex = je+j+1;
////                Real div_rest = b0_x2f(m,k,jex,i)*area2(m,k,jex,i) - (b0_x1f(m,k,jex,i+1)*area1(m,k,jex,i+1)-b0_x1f(m,k,jex,i)*area1(m,k,jex,i)) - (b0_x3f(m,k+1,jex,i)*area3(m,k+1,jex,i)-b0_x3f(m,k,jex,i)*area3(m,k,jex,i));
////                b0_x2f(m,k,jex+1,i) = div_rest/area2(m,k,jex+1,i);
////              }
//          }
//      });
//        par_for("usrboundaryx2_bfieldc", DevExeSpace(),0,(nmb-1),0,(n3-1),0,(n1-1),
//        KOKKOS_LAMBDA(int m, int k, int i) {
//            if (mb_bcs.d_view(m,BoundaryFace::inner_x2) == BoundaryFlag::user) {
//                Real lw, rw;
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,js,i) = lw*b0_x1f(m,k,js,i) + rw*b0_x1f(m,k,js,i+1);
//                lw = (x2f_(m,js+1)-x2v_(m,js))/(x2f_(m,js+1)-x2f_(m,js));
//                rw = (x2v_(m,js)-x2f_(m,js))/(x2f_(m,js+1)-x2f_(m,js));
//                bcc0(m,IBY,k,js,i) = lw*b0_x2f(m,k,js,i) + rw*b0_x2f(m,k,js+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,js,i) = lw*b0_x3f(m,k,js,i) + rw*b0_x3f(m,k+1,js,i);
//
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,je,i) = lw*b0_x1f(m,k,je,i) + rw*b0_x1f(m,k,je,i+1);
//                lw = (x2f_(m,je+1)-x2v_(m,je))/(x2f_(m,je+1)-x2f_(m,je));
//                rw = (x2v_(m,je)-x2f_(m,je))/(x2f_(m,je+1)-x2f_(m,je));
//                bcc0(m,IBY,k,je,i) = lw*b0_x2f(m,k,je,i) + rw*b0_x2f(m,k,je+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,je,i) = lw*b0_x3f(m,k,je,i) + rw*b0_x3f(m,k+1,je,i);
//
//              for (int j=0; j<ng; ++j) {
////                Real lw, rw;
////                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                bcc0(m,IBX,k,js+j,i) = lw*b0_x1f(m,k,js+j,i) + rw*b0_x1f(m,k,js+j,i+1);
////                lw = (x2f_(m,js+j+1)-x2v_(m,js+j))/(x2f_(m,js+j+1)-x2f_(m,js+j));
////                rw = (x2v_(m,js+j)-x2f_(m,js+j))/(x2f_(m,js+j+1)-x2f_(m,js+j));
////                bcc0(m,IBY,k,js+j,i) = lw*b0_x2f(m,k,js+j,i) + rw*b0_x2f(m,k,js+j+1,i);
////                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                bcc0(m,IBZ,k,js+j,i) = lw*b0_x3f(m,k,js+j,i) + rw*b0_x3f(m,k+1,js+j,i);
//                u0_(m,IEN,k,js-j-1,i) -= 0.5*(SQR(bcc0(m,IBX,k,js,i))+SQR(bcc0(m,IBY,k,js,i))+SQR(bcc0(m,IBZ,k,js,i)));
//
//                Real lw, rw;
//                lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                bcc0(m,IBX,k,js-j-1,i) = lw*b0_x1f(m,k,js-j-1,i) + rw*b0_x1f(m,k,js-j-1,i+1);
//                lw = (x2f_(m,js-j-1+1)-x2v_(m,js-j-1))/(x2f_(m,js-j-1+1)-x2f_(m,js-j-1));
//                rw = (x2v_(m,js-j-1)-x2f_(m,js-j-1))/(x2f_(m,js-j-1+1)-x2f_(m,js-j-1));
//                bcc0(m,IBY,k,js-j-1,i) = lw*b0_x2f(m,k,js-j-1,i) + rw*b0_x2f(m,k,js-j-1+1,i);
//                lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                bcc0(m,IBZ,k,js-j-1,i) = lw*b0_x3f(m,k,js-j-1,i) + rw*b0_x3f(m,k+1,js-j-1,i);
//                u0_(m,IEN,k,js-j-1,i) += 0.5*(SQR(bcc0(m,IBX,k,js-j-1,i))+SQR(bcc0(m,IBY,k,js-j-1,i))+SQR(bcc0(m,IBZ,k,js-j-1,i)));
//              }
//            }
//            if (mb_bcs.d_view(m,BoundaryFace::outer_x2) == BoundaryFlag::user) {
//              for (int j=0; j<ng; ++j) {
////                  Real lw, rw;
////                  lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                  rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
////                  bcc0(m,IBX,k,je-j,i) = lw*b0_x1f(m,k,je-j,i) + rw*b0_x1f(m,k,je-j,i+1);
////                  lw = (x2f_(m,je-j+1)-x2v_(m,je-j))/(x2f_(m,je-j+1)-x2f_(m,je-j));
////                  rw = (x2v_(m,je-j)-x2f_(m,je-j))/(x2f_(m,je-j+1)-x2f_(m,je-j));
////                  bcc0(m,IBY,k,je-j,i) = lw*b0_x2f(m,k,je-j,i) + rw*b0_x2f(m,k,je-j+1,i);
////                  lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                  rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
////                  bcc0(m,IBZ,k,je-j,i) = lw*b0_x3f(m,k,je-j,i) + rw*b0_x3f(m,k+1,je-j,i);
//                  u0_(m,IEN,k,je+j+1,i) -= 0.5*(SQR(bcc0(m,IBX,k,je,i))+SQR(bcc0(m,IBY,k,je,i))+SQR(bcc0(m,IBZ,k,je,i)));
//
//                  Real lw, rw;
//                  lw = (x1f_(m,i+1)-x1v_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                  rw = (x1v_(m,i)-x1f_(m,i))/(x1f_(m,i+1)-x1f_(m,i));
//                  bcc0(m,IBX,k,je+j+1,i) = lw*b0_x1f(m,k,je+j+1,i) + rw*b0_x1f(m,k,je+j+1,i+1);
//                  lw = (x2f_(m,je+j+1+1)-x2v_(m,je+j+1))/(x2f_(m,je+j+1+1)-x2f_(m,je+j+1));
//                  rw = (x2v_(m,je+j+1)-x2f_(m,je+j+1))/(x2f_(m,je+j+1+1)-x2f_(m,je+j+1));
//                  bcc0(m,IBY,k,je+j+1,i) = lw*b0_x2f(m,k,je+j+1,i) + rw*b0_x2f(m,k,je+j+1+1,i);
//                  lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                  rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
//                  bcc0(m,IBZ,k,je+j+1,i) = lw*b0_x3f(m,k,je+j+1,i) + rw*b0_x3f(m,k+1,je+j+1,i);
//                  u0_(m,IEN,k,je+j+1,i) += 0.5*(SQR(bcc0(m,IBX,k,je+j+1,i))+SQR(bcc0(m,IBY,k,je+j+1,i))+SQR(bcc0(m,IBZ,k,je+j+1,i)));
//              }
//            }
//        });
//    }
  return;
}



//----------------------------------------------------------------------------------------
void SourceFunc(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0, w0wb;
    DvceArray4D<Real> phi0_x1f;
    DvceArray4D<Real> phicc0;
    bool use_etotgrav = false;
    bool use_wellbalance_static = false;
    bool use_wellbalance_dynamic = false;
    DvceArray5D<Real> wbq0_;    // the well-balanced background cache (BuildWBCache)
    Real gamma;
    EOS_Data eos;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
      use_etotgrav = pmbp->phydro->use_etotgrav;
      use_wellbalance_static = pmbp->phydro->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->phydro->use_wellbalance_dynamic;
      wbq0_ = pmbp->phydro->wbq0;
      phi0_x1f = pmbp->phydro->phi0.x1f;
      phicc0 = pmbp->phydro->phicc0;
      w0wb = pmbp->phydro->w0wb;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
      use_etotgrav = pmbp->pmhd->use_etotgrav;
      use_wellbalance_static = pmbp->pmhd->use_wellbalance_static;
      use_wellbalance_dynamic = pmbp->pmhd->use_wellbalance_dynamic;
      wbq0_ = pmbp->pmhd->wbq0;
      phi0_x1f = pmbp->pmhd->phi0.x1f;
      phicc0 = pmbp->pmhd->phicc0;
      w0wb = pmbp->pmhd->w0wb;
    }
    
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
    
//    Real ap = 9.44e9;
//    Real omega = 2.06e-5;
//    Real grav_acc = -942.0;
//    Real Rgas = 4.593e7;
    
//    ParameterInput* pin;
    Real grav_acc = -pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    const bool tide = pm->pgen->hot_jupiter_param.stellar_tide;
    Real omega = pm->pgen->hot_jupiter_param.omega;
    // rot_potential: the radial centrifugal force is carried by the potential ONLY
    // through
    // the well-balanced source (the plain gravity source knows gravity alone), so the
    // explicit radial term is dropped only when that scheme is on; with etotgrav the
    // energy flux carries the potential either way
    const bool rotpot = pm->pgen->hot_jupiter_param.rot_potential;
    const bool rotpot_src = rotpot && use_wellbalance_dynamic;

    // Radial grid stretch, copied out of the Mesh so the device lambdas below capture
    // plain values. Radii rebuilt here with CellCenterX/LeftEdgeX are UNSTRETCHED, so
    // each is passed through ApplyRStretch to match the Coordinates arrays.
    const bool str_r_ = pm->use_grid_stretch_r;
    const bool str_rp_ = pm->use_grid_stretch_r_poly;
    const Real fstr_r_ = pm->fStretchR;
    const Real rmin_ = pm->mesh_size.x1min;
    const Real rmax_ = pm->mesh_size.x1max;
    Real cpoly_[NSTRETCH_R_POLY];
    for (int n=0; n<NSTRETCH_R_POLY; ++n) cpoly_[n] = pm->fStretchRPoly[n];
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    
    Real time = pm->time;
    
    picket_fence_two_stream_RT(pm, bdt);

    par_for("usrsource", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        
        Real x1v, x2v, x3v;
        if (use_spherical_polar) {
          x1v = x1v_(m,i);
          x2v = x2v_(m,j);
          x3v = x3v_(m,k);
        } else {
          // ASSIGN the x1v/x2v/x3v declared above -- do NOT redeclare them.  These three
          // lines used to read `Real x1v = ...`, which SHADOWED the outer variables and
          // threw the values away at the closing brace, leaving the outer ones
          // uninitialised for every non-spherical-polar grid.  Nothing exercised that
          // path until the cubed sphere did, and the result was the whole domain going
          // NaN on the first source-term application.
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
          ApplyRStretch(str_r_, fstr_r_, str_rp_, cpoly_, rmin_, rmax_, x1v);
          if (use_cubed_sphere_) x1v = x1v_(m,i);   // stored centroid/face, as sp

          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);

          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        }
        
        Real lam, phi, z, theta, r;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
          z = x1v-ap;
          r = x1v;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
          z = x1v-ap;
          r = x1v;
        } else {
          lam = x3v*iap;
          phi = x2v*iap;
          z = x1v;
        }
        Real rho = w0(m,IDN,k,j,i);
        Real p = eos.Pressure(w0(m,IDN,k,j,i), w0(m,IEN,k,j,i));
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,i),p);
        
        Real area_r = area1(m,k,j,i+1);
        Real area_l = area1(m,k,j,i);
        Real vol = volume(m,k,j,i);
        
        // gravity
        const Real grav_r = GravAccAt(grav_acc, ap, r, grav_pmass);
        Real src = bdt*grav_r*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
        }
        if (use_wellbalance_static) {
            src = bdt*grav_r*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        if (use_wellbalance_dynamic) {
          // The source is the background's own pressure difference across the cell, so it
          // has to come from the same entry point the reconstruction uses: getWBq0, which
          // dispatches to the ideal-gas closed forms or to the general-EOS background and
          // returns pressure directly, not an energy to be multiplied by (gamma-1).
          Real pl,pr,dum1,dum2,dum3;
          if (pmbp->phydro != nullptr) {
              WBReadCache(wbq0_, WBVar::wb_pres, m, k, j, i, dum1, pl, dum2, pr, dum3);
          } else if (pmbp->pmhd != nullptr) {
              WBReadCache(wbq0_, WBVar::wb_pres, m, k, j, i, dum1, pl, dum2, pr, dum3);
          }
          src = bdt*(area_r*(pr-p)+area_l*(p-pl))/vol;
        }
        u0(m,IM1,k,j,i) += src;
        
        // Forces in the corotating frame
        if (use_spherical_polar) {
          Real vtheta = w0(m,IVY,k,j,i);
          Real vphi = w0(m,IVZ,k,j,i);
          Real vr = w0(m,IVX,k,j,i);
          Real sine = sin(theta);
          Real cosine = cos(theta);
          Real oor = SQR(omega)*r*sine;
          Real cor = 2.0*omega*vphi;
          u0(m,IM2,k,j,i) += rho*(cor+oor)*cosine*bdt;
          u0(m,IM3,k,j,i) += -rho*2.0*omega*(vr*sine+vtheta*cosine)*bdt;
          // rot_potential: the RADIAL centrifugal force lives in the potential (the
          // well-balanced source and, with etotgrav, the energy flux carry it); only the
          // theta part and Coriolis stay explicit here
          u0(m,IM1,k,j,i) += rho*(cor + (rotpot_src ? 0.0 : oor))*sine*bdt;
//          if (!use_etotgrav)
          u0(m,IEN,k,j,i) += rho*oor*(((rotpot && use_etotgrav) ? 0.0 : vr*sine)
                                      + vtheta*cosine)*bdt;
          // The host star's tidal term, on top of the planet's own centrifugal `oor`.
          // Together they make up the Hill acceleration Omega^2 (3x, 0, -z); see
          // TideAccR/TideAccT/TideAccP. mu is the substellar direction cosine, the same
          // sin(theta)*cos(phi) the radiative transfer uses for the stellar beam.
          if (tide) {
            const Real cosph = cos(phi);
            const Real sinph = sin(phi);
            const Real mu = sine*cosph;
            const Real atr = TideAccR(omega, r, mu, tide);
            const Real att = TideAccT(omega, r, sine, cosine, cosph, tide);
            const Real atp = TideAccP(omega, r, sine, sinph, cosph, tide);
            u0(m,IM1,k,j,i) += rho*atr*bdt;
            u0(m,IM2,k,j,i) += rho*att*bdt;
            u0(m,IM3,k,j,i) += rho*atp*bdt;
            // Tidal work. Unconditional, like the centrifugal work above: the tidal
            // potential is deliberately NOT carried in phicc0, so use_etotgrav does not
            // already account for it.
            u0(m,IEN,k,j,i) += rho*(atr*vr + att*vtheta + atp*vphi)*bdt;
          }
        } else if (use_cubed_sphere_) {
          // CUBED SPHERE: the same corotating-frame forces the spherical-polar branch
          // above applies -- full Coriolis AND centrifugal -- instead of the equatorial
          // beta-plane below, which is written for the CARTESIAN BOX and which the cubed
          // sphere used to fall through into.  The coordinate block a few lines up
          // branches three ways; this one branched only two, so the omission was silent.
          //
          // The beta-plane itself is not wrong -- f = 2 omega lam IS the standard
          // equatorial form (f = beta*y with beta = 2 omega/R and y = R lam) -- but it is
          // an EQUATORIAL approximation on a GLOBAL grid: against the correct 2 omega
          // sin(lam) it is 11 % too strong at 45 deg, 21 % at 60 deg and 57 % at the
          // pole, it drops the centrifugal term entirely, and its two terms assume an
          // orthonormal (x,y) pair while IM2/IM3 here are COVARIANT components on the
          // non-orthogonal gnomonic tangent basis.  Every cubed-sphere vs spherical-polar
          // comparison of this problem was confounded by that difference.
          //
          // Done in CARTESIAN and then projected, which avoids writing the curvilinear
          // Coriolis by hand in a non-orthogonal basis:
          //     a = -2 Omega x v + Omega^2 (x, y, 0),   Omega = omega zhat,
          // and since m_i = V.e_i is COVARIANT, the momentum source is rho (a.e_i) with
          // e_i = (rhat, e_xi, e_eta).  w0's two angular velocities are CONTRAVARIANT on
          // that same pair -- see Coordinates::GnomonicEquiangleRaiseVelMHD, which builds
          // them as v2 = (m2 - c m3)/(d det) -- so V is their plain combination with the
          // basis vectors, with no metric factor.  zhat is the rotation axis on both
          // grids: CSCellAngles takes the colatitude from the Cartesian z of the same
          // chart, in the convention the spherical-polar branch uses.
          const Real xi_c = 0.25*M_PI*x2v;
          const Real eta_c = 0.25*M_PI*x3v;
          const int pnl = mbpanel_.d_view(m);
          Real qc[3], e1[3], e2[3];
          cubed_sphere::PanelToCart(pnl, xi_c, eta_c, qc);
          cubed_sphere::PanelTangents(pnl, xi_c, eta_c, e1, e2);
          const Real qn = 1.0/sqrt(qc[0]*qc[0] + qc[1]*qc[1] + qc[2]*qc[2]);
          const Real rh0 = qc[0]*qn, rh1 = qc[1]*qn, rh2 = qc[2]*qn;
          const Real v1 = w0(m,IVX,k,j,i);
          const Real v2 = w0(m,IVY,k,j,i);
          const Real v3 = w0(m,IVZ,k,j,i);
          const Real vcx = v1*rh0 + v2*e1[0] + v3*e2[0];
          const Real vcy = v1*rh1 + v2*e1[1] + v3*e2[1];
          const Real xx = r*rh0;
          const Real yy = r*rh1;
          const Real acx =  2.0*omega*vcy + SQR(omega)*xx;
          const Real acy = -2.0*omega*vcx + SQR(omega)*yy;
          // rot_potential: drop the radial projection of the centrifugal acceleration,
          // Omega^2 r sin^2(theta) = Omega^2 (x rh0 + y rh1), which the potential carries
          const Real acr_rot = rotpot_src ? SQR(omega)*(xx*rh0 + yy*rh1) : 0.0;
          u0(m,IM1,k,j,i) += rho*(acx*rh0 + acy*rh1 - acr_rot)*bdt;
          u0(m,IM2,k,j,i) += rho*(acx*e1[0] + acy*e1[1])*bdt;
          u0(m,IM3,k,j,i) += rho*(acx*e2[0] + acy*e2[1])*bdt;
          // Only the centrifugal part does work -- Coriolis is perpendicular to v -- and
          // this is the same term the spherical-polar branch adds unconditionally.
          u0(m,IEN,k,j,i) += rho*(SQR(omega)*(xx*vcx + yy*vcy)
                                  - ((rotpot && use_etotgrav) ?
                                     SQR(omega)*(xx*rh0 + yy*rh1)*v1 : 0.0))*bdt;
        } else {
          // corotating beta-plane approximation e.g. Fromang+2016
          Real omega1 = omega*lam;
          u0(m,IM2,k,j,i) += -2.0*rho*omega1*(-w0(m,IVZ,k,j,i))*bdt;
          u0(m,IM3,k,j,i) += -2.0*rho*omega1*w0(m,IVY,k,j,i)*bdt;
        }

//        // Newtonian cooling
//        Real Teq, itrad;
//        get_eq_Tp(lam, phi, p, Teq);
//        get_itrad(p,itrad);
////        Real t0 = 2.16e6;
////        Real ff = (t0-time)/t0;
////        ff = (ff < 0.0) ? 0.0 : ff;
////        Real gg = 2.0*ff;
////        itrad /= pow(10.0,gg);
//        Real Tnew = (T + Teq*itrad*bdt)/(1.0 + itrad*bdt);
//        u0(m,IEN,k,j,i) -= w0(m,IEN,k,j,i)*(Tnew-Teq)/T*itrad*bdt;

        // Rayleigh drag (initial relaxation)
        Real dyntime = 2.0*M_PI/omega;
        Real tau1 = dyntime / 10.0;
        Real tau2 = dyntime;
        Real t1   = 2.0 * dyntime;
        Real t2   = 5.0 * dyntime;
        Real itdrag, fredux;
        if(time < t1) {
          itdrag = 1.0 / tau1;
        } else if(time < t2)
          {
            Real alp = (time - t1) / (t2 - t1);
            itdrag   = 1.0 / tau1 * pow(tau1 / tau2, alp);
          }
        else {
          itdrag = 0.0;
        }
        if (time < t2) {
          fredux = itdrag*bdt;///(1.0+itdrag*bdt);
          u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
          u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
          u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        }

        // Top sponge layer
        Real bar = 1.0e6;
        Real logpl = log(1.0e-6*bar);
        Real logpt = log(1.0e-7*bar);
        Real logp = log(p);
        Real fdrag = 1.0 - (logp-logpt)/(logpl-logpt); // high p = 0, low p = 1
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/1.0e3;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
        // Bottom sponge layer
        logpl = log(1.0e2*bar);
        logpt = log(5.0e1*bar);
        logp = log(p);
        fdrag = (logp-logpt)/(logpl-logpt); // high p = 1, low p = 0
        fdrag = (fdrag < 0.0) ? 0.0 : fdrag;
        fdrag = (fdrag > 1.0) ? 1.0 : fdrag;
        itdrag = fdrag/1.0e3;
        fredux = itdrag*bdt; ///(1.0+itdrag*bdt);
//        u0(m,IM1,k,j,i) -= u0(m,IM1,k,j,i)*fredux;
        u0(m,IM2,k,j,i) -= u0(m,IM2,k,j,i)*fredux;
        u0(m,IM3,k,j,i) -= u0(m,IM3,k,j,i)*fredux;
    });

    return;
}

void double_gray_two_stream_RT(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;
        
    DvceArray5D<Real> u0, w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x1f_ = pmbp->pcoord->xx1f;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    const bool correct_spherical = false;
    const bool test_oned = false;
        
    Real gamma;
    EOS_Data eos;
    // wtemp is the temperature ConsToPrim already solved for the current w0. It is
    // allocated ONLY for a general EOS, so it is read only on the general branch of
    // TempKelvin; for an ideal gas it is a zero-size View, captured and never touched.
    DvceArray4D<Real> wtemp_;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
    }
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;

    auto area1 = pmbp->pcoord->area.x1f;
    auto volume = pmbp->pcoord->volume;
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real grav = 942.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 100.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 5;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
    int nclip = 0;
    const Real demax = rt_de_max;
    par_reduce_clip3("2stream_rt", 0, nmb1, ks, ke, js, je, nclip,
    KOKKOS_LAMBDA(const int m, const int k, const int j, int &nc) {
        constexpr int NN = 270;
        Real tau_ir_down_f[NN];
        Real F_v_down_f[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real B[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real rtop = x1f_(m,ie+1);
        Real rbot = x1f_(m,is);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = -1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real pm1 = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie),w0(m,IEN,k,j,ie));
        Real pf = exp((log(p)+log(pm1))/2.0);
        Real tau_v_f = 0.0;//pf/(grav/kap_v);
        Real tau_ir_f = 0.0;//pf/(grav/kap_ir);
        
        tau_ir_down_f[ie+1] = tau_ir_f;
        F_v_down_f[ie+1] = Fstar*mu0;//*exp(-tau_v_f/muf);
        F_v_down_f[ie+1] = (mu0 > 0.0) ? F_v_down_f[ie+1] : 0.0;
        
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
          Real kap_ir = 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = x1v_(m,i);
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
            
          Real rf = x1f_(m,i);
          Real rf1 = x1f_(m,i+1);

//          Real muf = sqrt(1.0-SQR(ap/rf)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r1));
          Real muf = (mu0 < mucr) ? mucr : mu0;
          if (test_oned) muf = mu0;
          Real dtau_v = kap_v*rho*dr;
          Real dtau_ir = kap_ir*rho*dr;
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf); // Zhang+2023
          Real trans = exp(-dtau_v/muf);
          F_v_down_f[i] = F_v_down_f[i+1]*trans*fac;
          tau_ir_down_f[i] = tau_ir_down_f[i+1] + dtau_ir;
        }
        p = PresFromEint(eos,gm1,w0(m,IDN,k,j,is-1),w0(m,IEN,k,j,is-1));
        rho = w0(m,IDN,k,j,is-1);
        T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,is-1),p);
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        I_ir_down_f[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rf = x1f_(m,i);
          Real rf1 = x1f_(m,i+1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rf1/rf);
          Real dtau = tau_ir_down_f[i]-tau_ir_down_f[i+1];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i]+B[i+1])/2.0;
          I_ir_down_f[i] = (I_ir_down_f[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        I_ir_up_f[is] = Iint;
        // up-sweep
        for (int i=is+1; i<ie+2; ++i) {
          Real rf = x1f_(m,i);
          Real rfm1 = x1f_(m,i-1);
            
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rfm1/rf);
          Real dtau = tau_ir_down_f[i-1]-tau_ir_down_f[i];
          Real trans = exp(-dtau/mu1);
          Real Bavg = (B[i-1]+B[i])/2.0;
          I_ir_up_f[i] = (I_ir_up_f[i-1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        // flux divergence
        for (int i=is; i<ie+1; ++i) {
          Real Ft = 2.0*M_PI*mu1*(I_ir_up_f[i+1]-I_ir_down_f[i+1])-F_v_down_f[i+1];
          Real Fb = 2.0*M_PI*mu1*(I_ir_up_f[i]-I_ir_down_f[i])-F_v_down_f[i];
          Real area_t = area1(m,k,j,i+1);
          Real area_b = area1(m,k,j,i);
          Real vol = volume(m,k,j,i);
          Real src = -(Ft-Fb)/dx1(m,k,j,i);
          if (correct_spherical) src = -(Ft*area_t-Fb*area_b)/vol;
          Real de = src*bdt;
          if (demax > 0.0) {
            const Real dl = LimitRTSource(de, w0(m,IEN,k,j,i), demax);
            if (dl != de) { ++nc; de = dl; }
          }
          u0(m,IEN,k,j,i) += de;
        }
        
    });
    RTSourceLimiterWarn(nclip);
    
    return;
}

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt) {
  // the cubed sphere needs the cell's PANEL to turn (x2,x3) into a direction
  const bool use_cubed_sphere_ = pm->use_cubed_sphere;
  auto &mbpanel_ = pm->pmb_pack->pmb->mb_panel;
    auto &indcs = pm->mb_indcs;
    int &ng = indcs.ng;
    int n1 = indcs.nx1 + 2*ng;
    int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
    int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
    int &is = indcs.is;  int &ie  = indcs.ie;
    int &js = indcs.js;  int &je  = indcs.je;
    int &ks = indcs.ks;  int &ke  = indcs.ke;
    auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
    int nmb1 = pm->pmb_pack->nmb_thispack - 1;
    MeshBlockPack *pmbp = pm->pmb_pack;
    auto &size = pmbp->pmb->mb_size;

    DvceArray5D<Real> u0, w0;
    const bool use_spherical_polar = pm->use_spherical_polar;
    auto &x1v_ = pmbp->pcoord->x1v;
    auto &x2v_ = pmbp->pcoord->x2v;
    auto &x3v_ = pmbp->pcoord->x3v;
    const bool correct_spherical = false;
    const bool test_oned = false;
    
    Real gamma;
    EOS_Data eos;
    // wtemp is the temperature ConsToPrim already solved for the current w0. It is
    // allocated ONLY for a general EOS, so it is read only on the general branch of
    // TempKelvin; for an ideal gas it is a zero-size View, captured and never touched.
    DvceArray4D<Real> wtemp_;
    if (pmbp->phydro != nullptr) {
      u0 = pmbp->phydro->u0;
      w0 = pmbp->phydro->w0;
      wtemp_ = pmbp->phydro->wtemp;
      gamma = pmbp->phydro->peos->eos_data.gamma;
      eos = pmbp->phydro->peos->eos_data;
    } else if (pmbp->pmhd != nullptr) {
      u0 = pmbp->pmhd->u0;
      w0 = pmbp->pmhd->w0;
      wtemp_ = pmbp->pmhd->wtemp;
      gamma = pmbp->pmhd->peos->eos_data.gamma;
      eos = pmbp->pmhd->peos->eos_data;
    }
    
    Real r0, r1;
    r0 = pm->mesh_size.x1min;
    r1 = pm->mesh_size.x1max;
    auto dx1 = pmbp->pcoord->dx1;
    
//    Real Teq = 1469.0;
//    Real ap = 9.44e9;
//    Real Rgas = 4.593e7;
//    Real grav = 942.0;
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    const bool grav_pmass = pm->pgen->hot_jupiter_param.grav_point_mass;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Fstar = boltz_sigma*SQR(SQR(Tirr));
    Real Tint = 500.0;
    Real Iint = boltz_sigma/M_PI*SQR(SQR(Tint));
    Real mu1 = 1.0/1.66; //sqrt(3);

//    size_t scr_size = ScrArray1D<Real>::shmem_size(n1) * 4;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, 0,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_ir_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> I_down(member.team_scratch(0), n1);
//        ScrArray1D<Real> B(member.team_scratch(0), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(0), n1);
        constexpr int NN = 270;
        Real tau_ir_down[NN];
        Real I_down[NN];
        Real B[NN];
        Real Q_v[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else if (use_cubed_sphere_) {
          CSCellAngles(mbpanel_.d_view(m), x2v, x3v, theta, lam, phi);
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(50.0/90.0*M_PI/2.0);
        
        // down-sweep
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kap_v = 4.0e-3;
        Real kap_ir = 1.0e-2;
        Real tau_v = 0.0;//p/(grav/kap_v);
        Real tau_ir = 0.0;//p/(grav/kap_ir);
//        if (test_oned) {
//          tau_v = p/(grav/kap_v)/mu0;
//          tau_ir = p/(grav/kap_ir)/mu1;
//        }
        tau_ir_down[ie+1] = tau_ir;
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kap_v = 4.0e-3; // Rauscher & Menou 2012; Guillot 2010
            Real kap_ir = 1.0e-2; // 2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
          Real dr = dx1(m,k,j,i);
            
          Real x1v = x1v_(m,i);
            
          Real z, r;
          if (use_spherical_polar) {
            z = x1v-ap;
            r = x1v;
          } else {
            z = x1v;
          }
////          Real mu = sqrt(1.0-SQR(r1/r)*(1.0-SQR(mu0))); // Li & Shibata 2006
          Real mucr = sqrt(1.0-SQR(r0/r));
//          Real mu = (mu0 < mucr) ? mucr : mu0;
//          if (test_oned) mu = mu0;
//          Real dtau_v = kap_v*rho*dr;
          Real delta = dr/r;
          Real drcor = r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          if (test_oned) drcor = dr/mu0;
          Real dtau_v = kap_v*rho*drcor;
          Real dtau_ir = kap_ir*rho*dr;
            
          tau_v += dtau_v;
          Real fac = 1.0;
//          if (correct_spherical) fac = SQR(r1/r); // Zhang+2023
//          Real Q_v = kap_v*rho*Fstar*fac*exp(-tau_v/mu); // Zhang+2023
//          Q_v = (mu0 > 0.0) ? Q_v : 0.0;
          Real Qv = kap_v*rho*Fstar*fac*exp(-tau_v);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
//          u0(m,IEN,k,j,i) += Q_v*bdt;
        
          tau_ir_down[i] = tau_ir;
          tau_ir += dtau_ir/mu1;
          if (i==is) tau_ir_down[i-1] = tau_ir;
        }
        p = PresFromEint(eos,gm1,w0(m,IDN,k,j,is-1),w0(m,IEN,k,j,is-1));
        rho = w0(m,IDN,k,j,is-1);
        T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,is-1),p);
        B[is-1] = boltz_sigma/M_PI*SQR(SQR(T));
        
        Real rtop = x1v_(m,ie+1);
        I_down[ie+1] = 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real r = x1v_(m,i);
          Real rp1 = x1v_(m,i+1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rp1/r);
          Real dtau = tau_ir_down[i]-tau_ir_down[i+1];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i)+B(i+1))/2.0;
          I_down[i] = (I_down[i+1]*trans + Bavg*(1.0-trans))*fac;
        }
        
        Real rbot = x1v_(m,is-1);
        Real I_up = Iint;
        // up-sweep
        for (int i=is; i<ie+1; ++i) {
          Real r = x1v_(m,i);
          Real rm1 = x1v_(m,i-1);
          Real fac = 1.0;
          if (correct_spherical) fac = SQR(rm1/r);
          Real dtau = tau_ir_down[i-1]-tau_ir_down[i];
          Real trans = exp(-dtau);
          Real Bavg = B[i];//(B(i-1)+B(i))/2.0;
          I_up = (I_up*trans + Bavg*(1.0-trans))*fac;
          Real J = (I_up+I_down[i])/2.0;
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),
                           TGuess(wtemp_, m, k, j, i),p,T);
            Real kap_ir = 1.0e-2; //2.28e-5*pow(p*cgs2Pa,0.53); // Komocek+2017
          if (test_oned) kap_ir = 1.0e-2;
//          Real Q_ir = 4.0*M_PI*kap_ir*rho*(J-B(i));
//          u0(m,IEN,k,j,i) += Q_ir*bdt;
            
          Real kk = -4.0*M_PI*kap_ir*rho*boltz_sigma/M_PI*bdt;
          Real e0, e;
          if (eos.IsGeneral()) {
            // Same implicit balance e - kk T^4 = bb, but a general EOS has no e = c_v T
            // with constant c_v, so Newton-Raphson runs on the internal energy directly:
            // F(e) = e - kk T(e)^4 - bb, with dT/de = 1/(d c_v). c_v is per unit mass in
            // CODE units while T here is in Kelvin, hence the temp_cgs factor.
            e0 = w0(m,IEN,k,j,i);
            Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
            e = e0;
            // `tc` is T in CODE temperature, carried alongside the kelvin `T` so that
            // neither EOS call has to solve for it. The two-argument SpecificHeatCv(d,e)
            // and Temperature(d,e) are the COLD-START forms -- the first is documented
            // "setup-time use only" because it inverts e(d,T) itself, and the second
            // brackets from scratch -- so using them here cost two full root finds per
            // Newton step, per cell, per stage. T and e are consistent at the top of
            // every iteration, so c_v can be evaluated at the temperature already in
            // hand, and the refresh below only needs the previous T as a warm start.
            Real tc = T/eos.temp_cgs;
            for (int n=0; n<100; ++n) {
              Real dTde = eos.temp_cgs/(rho*eos.SpecificHeatCv(rho,e,tc));
              Real de = e - kk*SQR(SQR(T)) - bb;
              e -= de / (1.0 - 4.0*kk*T*T*T*dTde);
              if (fabs(de) <= 1.0e-10*e)
                break;
              tc = eos.Temperature(rho,e,tc);
              T = tc*eos.temp_cgs;
            }
          } else {
            Real cv = Rgas*rho*igm1;
            e0 = cv*T;
            Real bb = 4.0*M_PI*kap_ir*rho*J*bdt + Q_v[i]*bdt + e0;
            // Newton-Raphson
            for (int n=0; n<100; ++n) {
              e = cv*T;
              Real de = e - kk*SQR(SQR(T)) - bb;
              T -= de / (cv - 4.0*kk*T*T*T);
              if (fabs(de) <= 1.0e-10*e)
                break;
            }
          }
          u0(m,IEN,k,j,i) += e-e0;
        }
        
    });
    
    return;
}


//----------------------------------------------------------------------------------------
//! \fn void DhjPhotosphereDump
//! \brief The tau = 2/3 PHOTOSPHERE, per band and g-point, from the run's OWN opacity.
//!
//! WHY THIS EXISTS.  An isobar is not the observed surface.  What a telescope sees is the
//! level where the vertical optical depth reaches ~2/3, and on an ultra-hot Jupiter that
//! level is not one surface: it sits DEEP in a window between water bands and HIGH in a
//! band centre, and it rises on the dayside where H- opacity climbs with temperature.
//! Presenting a single pressure level as "the observed surface" hides all three effects.
//!
//! The opacity is not re-derived here.  It is the identical lookup the longwave sweep
//! uses -- `ck_kappa(...) + rt_kc` -- so the photosphere reported is the one the run's
//! own radiative transfer saw, not a model of it.  Requires the correlated-k path and
//! at least one RT evaluation, because the per-cell tables are allocated lazily: restart
//! and take a single cycle.
//!
//! DEFINITION, stated because it is a choice.  Optical depth is accumulated VERTICALLY
//! (kappa*rho*dr, with no 1/mu slant factor) downwards from the top, starting from the
//! hydrostatic column above the domain that the sweep's own top boundary term adds.  For
//! each (band, g-point) the crossing tau = 2/3 is located by interpolating in log tau
//! against log p.  A band's reported level is the g-WEIGHTED MEAN of log p over its
//! g-points, i.e. the emission-weighted level; the spread between its lowest and
//! g-point is the window-to-line-core range, which is the point of the diagnostic.  The
//! broadband level weights bands by the local Planck fraction f_b(T).
//!
//! A g-point whose tau = 2/3 lies ABOVE the domain top is reported as p = 0: the model
//! does not contain that emitting level, which is a fact about the grid worth seeing
//! rather than clamping silently away.

void DhjPhotosphereDump(ParameterInput *pin, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const std::string fname = pin->GetOrAddString("problem", "photosphere_dump", "");
  if (fname.empty()) return;
  if (!rt_ck || ck_lk_ptr == nullptr || rt_xT_ptr == nullptr || rt_kc_ptr == nullptr) {
    if (global_variable::my_rank == 0) {
      std::cout << "### photosphere_dump: needs the correlated-k path AND at least one RT"
                << " evaluation (the per-cell tables are allocated lazily) -- restart and"
                << " run one cycle. Nothing written." << std::endl;
    }
    return;
  }
  const Real TAU_PH = 2.0/3.0;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  const int js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;

  auto xT = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_xT_ptr);
  auto xP = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_xP_ptr);
  auto pbr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_pb_ptr);
  auto Tc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_T_ptr);
  auto kc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_kc_ptr);
  auto cf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *rt_cf_ptr);
  auto lk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_lk_ptr);
  auto gw = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_gw_ptr);
  auto wl = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_wl_ptr);
  auto pf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), *ck_pf_ptr);
  auto w0h = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                 (pmbp->phydro != nullptr) ? pmbp->phydro->w0 : pmbp->pmhd->w0);
  auto dx1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->dx1);
  auto x1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto x2v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x2v);
  auto x3v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x3v);

  const HotJupiterParam &hj = pm->pgen->hot_jupiter_param;
  std::ostringstream os;
  os.precision(7);
  os << std::scientific;
  if (global_variable::my_rank == 0) {
    os << "# AthenaK tau=2/3 photosphere, from the run's own correlated-k opacity\n"
       << "# t = " << pm->time << " s   cycle " << pm->ncycle << "\n"
       << "# bands " << CK_NB << "  g-points " << CK_NG << "\n"
       << "# band edges [um], descending:";
    for (int b=0; b<=CK_NB; ++b) os << " " << wl(b);
    os << "\n# lat_deg lon_deg band lam_lo_um lam_hi_um p_eff_bar p_deep_bar p_high_bar"
       << " r_eff_cm\n"
       << "# p = 0 means tau=2/3 lies above the domain top; band index " << CK_NB
       << " is the broadband, Planck-weighted level\n";
  }

  Real tau[CK_NG], plev[CK_NG], rlev[CK_NG];
  for (int m=0; m<nmb; ++m) {
    for (int k=ks; k<=ke; ++k) {
      const Real lon = (x3v(m,k) - M_PI)*180.0/M_PI;
      for (int j=js; j<=je; ++j) {
        const Real lat = 90.0 - x2v(m,j)*180.0/M_PI;
        const Real mu0 = cf(m,k,j,3);
        Real pbb[CK_NB+1], rbb[CK_NB+1], pdeep[CK_NB+1], phigh[CK_NB+1];
        for (int b=0; b<CK_NB; ++b) {
          // start each g-point from the hydrostatic column ABOVE the domain, exactly as
          // the sweep's own top boundary term does
          const Real ptop = pbr(m,k,j,ie+1);
          const int iT0 = static_cast<int>(xT(m,k,j,ie+1));
          const int iP0 = static_cast<int>(xP(m,k,j,ie+1));
          const Real fT0 = xT(m,k,j,ie+1) - static_cast<Real>(iT0);
          const Real fP0 = xP(m,k,j,ie+1) - static_cast<Real>(iP0);
          const Real gtop = EffGravAt(hj.grav, hj.ap, x1v(m,ie+1), hj.grav_point_mass,
                                      hj.omega, mu0, hj.stellar_tide);
          for (int g=0; g<CK_NG; ++g) {
            const Real kap = ck_kappa(lk, iT0, fT0, iP0, fP0, b, g) + kc(m,b,ie+1,k,j);
            tau[g] = kap*ptop*1.0e6/gtop;
            plev[g] = 0.0;
            rlev[g] = 0.0;
          }
          for (int i=ie; i>=is; --i) {
            const int iT = static_cast<int>(xT(m,k,j,i));
            const int iP = static_cast<int>(xP(m,k,j,i));
            const Real fT = xT(m,k,j,i) - static_cast<Real>(iT);
            const Real fP = xP(m,k,j,i) - static_cast<Real>(iP);
            const Real drho = w0h(m,IDN,k,j,i)*dx1(m,k,j,i);
            for (int g=0; g<CK_NG; ++g) {
              if (plev[g] != 0.0) continue;
              const Real t0 = tau[g];
              tau[g] += (ck_kappa(lk, iT, fT, iP, fP, b, g) + kc(m,b,i,k,j))*drho;
              if (tau[g] >= TAU_PH && t0 < TAU_PH) {
                const Real f = (t0 > 0.0)
                    ? (std::log(TAU_PH) - std::log(t0))/(std::log(tau[g]) - std::log(t0))
                    : 1.0;
                const Real lp0 = std::log(pbr(m,k,j,i+1));
                const Real lp1 = std::log(pbr(m,k,j,i));
                plev[g] = std::exp(lp0 + f*(lp1 - lp0));
                rlev[g] = x1v(m,i+1) + f*(x1v(m,i) - x1v(m,i+1));
              }
            }
          }
          Real sw = 0.0, slp = 0.0, sr = 0.0, lo = 0.0, hi = 0.0;
          for (int g=0; g<CK_NG; ++g) {
            if (plev[g] <= 0.0) continue;
            sw += gw(g);
            slp += gw(g)*std::log(plev[g]);
            sr += gw(g)*rlev[g];
            lo = (lo == 0.0) ? plev[g] : fmax(lo, plev[g]);
            hi = (hi == 0.0) ? plev[g] : fmin(hi, plev[g]);
          }
          pbb[b] = (sw > 0.0) ? std::exp(slp/sw) : 0.0;
          rbb[b] = (sw > 0.0) ? sr/sw : 0.0;
          pdeep[b] = lo;
          phigh[b] = hi;
        }
        Real Bb[CK_NB], swb = 0.0, slb = 0.0, srb = 0.0;
        ck_planck_bands(pf, ck_pf_lTmin, ck_pf_idlT, 1.0, Tc(m,k,j,ie), Bb);
        for (int b=0; b<CK_NB; ++b) {
          if (pbb[b] <= 0.0) continue;
          swb += Bb[b];
          slb += Bb[b]*std::log(pbb[b]);
          srb += Bb[b]*rbb[b];
        }
        pbb[CK_NB] = (swb > 0.0) ? std::exp(slb/swb) : 0.0;
        rbb[CK_NB] = (swb > 0.0) ? srb/swb : 0.0;
        pdeep[CK_NB] = pbb[CK_NB];
        phigh[CK_NB] = pbb[CK_NB];
        for (int b=0; b<=CK_NB; ++b) {
          os << lat << " " << lon << " " << b << " "
             << ((b < CK_NB) ? wl(b+1) : 0.0) << " " << ((b < CK_NB) ? wl(b) : 0.0) << " "
             << pbb[b] << " " << pdeep[b] << " " << phigh[b] << " " << rbb[b] << "\n";
        }
      }
    }
  }

#if MPI_PARALLEL_ENABLED
  {
    // gather the per-rank text on rank 0, in rank order
    std::string mine = os.str();
    int len = static_cast<int>(mine.size());
    int nrank = global_variable::nranks;
    std::vector<int> lens(nrank), offs(nrank, 0);
    MPI_Gather(&len, 1, MPI_INT, lens.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    int tot = 0;
    if (global_variable::my_rank == 0) {
      for (int r=0; r<nrank; ++r) { offs[r] = tot; tot += lens[r]; }
    }
    std::vector<char> all((global_variable::my_rank == 0) ? tot : 1);
    MPI_Gatherv(mine.data(), len, MPI_CHAR, all.data(), lens.data(), offs.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD);
    if (global_variable::my_rank == 0) {
      std::ofstream f(fname);
      f.write(all.data(), tot);
      f.close();
      std::cout << "photosphere: written to '" << fname << "'" << std::endl;
    }
  }
#else
  {
    std::ofstream f(fname);
    f << os.str();
    f.close();
    std::cout << "photosphere: written to '" << fname << "'" << std::endl;
  }
#endif
  return;
}

KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td) {
    
    Real bar = 1.0e6;
    Real pl = 1.0e-3*bar;
    Real pt = log10(p/bar);
    Real ptl = log10(pl/bar);
    
    Real fn[13];
    Real fd[14];
    
    fn[0] = 1388.77348;
    fn[1] = 279.575848;
    fn[2] = -213.835822;
    fn[3] = 21.0010475;
    fn[4] = 100.938036;
    fn[5] = 12.7972336;
    fn[6] = -13.9266925;
    fn[7] = -3.70783272;
    fn[8] = 0.522370269;
    fn[9] = 0.320837882;
    fn[10]= 0.0451831612;
    fn[11]= 2.18195583e-3;
    fn[12]= 3.98938097e-6;

    fd[0] = 2152.06036;
    fd[1] = 29.3485512;
    fd[2] = -183.318696;
    fd[3] = 46.3893130;
    fd[4] = 19.8116485;
    fd[5] = -28.5473177;
    fd[6] = -2.52726545;
    fd[7] = 8.43627538;
    fd[8] = 2.62945375;
    fd[9] = -0.297098168;
    fd[10]= -0.286871487;
    fd[11]= -0.0590629443;
    fd[12]= -5.38679474e-3;
    fd[13]= -1.89972415e-4;
    
    Real Tnstar = 0.0;
    Real Tdstar = 0.0;
    Real Tnstar_pl = 0.0;
    Real Tdstar_pl = 0.0;
    Real ptn = 1.0;
    Real ptln = 1.0;
    for(int ilogp=0; ilogp<13; ilogp++) {
        Tnstar += fn[ilogp]*ptn;
        Tnstar_pl += fn[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }
    ptn = 1.0;
    ptln = 1.0;
    for(int ilogp=0; ilogp<14; ilogp++) {
        Tdstar += fd[ilogp]*ptn;
        Tdstar_pl += fd[ilogp]*ptln;
        ptn *= pt;
        ptln *= ptl;
    }

    Tn = Tnstar;
    Td = Tdstar;
//    if (p < pl) {
//        Tn = Tnstar_pl * exp(0.1*log10(p/pl));
//        if (Tn < 250) {
//            Tn = 250;
//        }
//        Td = Tdstar_pl * exp(0.015*log10(p/pl));
//        if (Td < 1000) {
//            Td = 1000;
//        }
//    }
//    if (Td < Tn) {
//        Td = Tn;
//    }
    Real x = log10(p / pl);

    Real Tn_new = Tnstar_pl;// * exp(0.1   * x);
    Real Td_new = Tdstar_pl;// * exp(0.015 * x);

    Tn_new = fmax(Tn_new, 250.0);
    Td_new = fmax(Td_new, 1000.0);

    // blend
    Tn = (p < pl) ? Tn_new : Tn;
    Td = (p < pl) ? Td_new : Td;

    // enforce Td >= Tn
    Td = fmax(Td, Tn);

    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T) {
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : Nt;
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
//    Real Tn, Td;
//    get_daynight_Tp(p, Tn, Td);
//
////    T = Tn;
//
//    Real Tn2 = Tn*Tn;
//    Real Tn4 = Tn2*Tn2;
//    Real Td2 = Td*Td;
//    Real Td4 = Td2*Td2;
//    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;
//
//    Real Tmid = sqrt(sqrt(Tmid4));
//    T = Tmid;
//
////    Real bar = 1.0e6;
////    Real Teq = 1469.0;
//////    Real Tirr = Teq*sqrt(2);
//////    Real Tint = 100.0;
//////    Real g = 942.0;
//////    Real mus = cos(50.0/90.0*M_PI/2.0);
//////    Real fH = 0.5;
//////    Real fK = 1.0/3.0;
//////    Real kap_v = 4.0e-3;
//////    Real kap_ir = 1.0e-2;
//////    Real gam = kap_v/kap_ir;
//////    Real tau = p/(g/kap_ir);
//////    tau = (tau < 0.0) ? 0.0 : tau;
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
////////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
//////    T = sqrt(sqrt(T4)); //2581.5574;
////    T = Teq;
    
    
    return;
}

template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T) {
    
    Real lgp = log10(p);
    Real dlgp = (lgparr(N-1)-lgparr(0))/(N-1);
    int Nt = std::floor((lgp - lgparr(0))/dlgp);
    int NN = (Nt < 0) ? 0 : Nt;
//    for (int it=Nt-2; it<Nt+3; ++it)
//    {
//        if (lgp < lgparr(it) && lgp >= lgparr(it-1)) {
//            NN = it-1;
//            break;
//        }
//    }
    T = Tarr(NN) + (Tarr(NN+1)-Tarr(NN))/(lgparr(NN+1)-lgparr(NN))*(lgp-lgparr(NN));
    T = (Nt < 0) ? Tarr(0) : T;
    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T) {
    
    Real Tn, Td;
    get_daynight_Tp(p, Tn, Td);
    
//    T = Tn;
    
    Real Tn2 = Tn*Tn;
    Real Tn4 = Tn2*Tn2;
    Real Td2 = Td*Td;
    Real Td4 = Td2*Td2;
    Real Tmid4 = 0.75*Tn4 + 0.25*Td4;

    Real Tmid = sqrt(sqrt(Tmid4));
    T = Tmid;
    
//    Real bar = 1.0e6;
//    Real Teq = 1469.0;
////    Real Tirr = Teq*sqrt(2);
////    Real Tint = 100.0;
////    Real g = 942.0;
////    Real mus = cos(50.0/90.0*M_PI/2.0);
////    Real fH = 0.5;
////    Real fK = 1.0/3.0;
////    Real kap_v = 4.0e-3;
////    Real kap_ir = 1.0e-2;
////    Real gam = kap_v/kap_ir;
////    Real tau = p/(g/kap_ir);
////    tau = (tau < 0.0) ? 0.0 : tau;
////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK) + 0.25*SQR(SQR(Tirr))*(mus/fH+SQR(mus)/gam/fK+(gam-SQR(mus)/gam/fK)*exp(-gam*tau/mus));
//////    Real T4 = 0.25*SQR(SQR(Tint))*(1.0/fH+tau/fK);
////    T = sqrt(sqrt(T4)); //2581.5574;
//    T = Teq;
    
    return;
}


template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                    const Real &ap, const bool &grav_pmass, const int &N,
                    const Real &zmax, View1D zarr, View1D logparr) {
    
//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;
        // gravity at THIS height, so the integration is consistent with the source term
        const Real gz = GravAccAt(grav_acc, ap, ap + zarr(n), grav_pmass);
        const Real fac = gz/Rgas*dz;

        get_wb_Tp(p,T);

        // n+1 is guarded: logparr holds N entries, so the last pass of this loop used to
        // integrate one step PAST THE END of it. A one-element write off the end of a
        // Kokkos allocation is silent in a Release build -- it lands in whatever follows
        // -- but it is undefined behaviour, a Kokkos Debug build aborts on it, and the
        // value written was never read. zarr(n) above still has to be set for every n.
        if (n+1 < N) {
          if (eos.IsGeneral()) {
            // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
            // branch is the same thing with rho = p/(Rgas T), kept in its original form.
            Real rho = DensFromPT(eos, Rgas, p, T);
            logparr(n+1) = logparr(n) + gz*dz*rho/p;
          } else {
            logparr(n+1) = logparr(n) + fac/T;
          }
        }
    }

    
    return;
}

KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {
    
//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_wb_Tp(p,T);
        rho = DensFromPT(eos, Rgas, p, T);
    } else {
        Real T0;
        get_wb_Tp(p0,T0);
        Real rho0 = DensFromPT(eos, Rgas, p0, T0);
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}

template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc,
                      const Real &ap, const bool &grav_pmass, const View1D &Tarr,
                      const View1D &lgparr, const int &N, const Real &zmax,
                      View1D zarr, View1D logparr) {

//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;
        // gravity at THIS height, so the integration is consistent with the source term
        const Real gz = GravAccAt(grav_acc, ap, ap + zarr(n), grav_pmass);
        const Real fac = gz/Rgas*dz;

        get_init_Tp_host(N, Tarr, lgparr, p, T);

        // n+1 is guarded: logparr holds N entries, so the last pass of this loop used to
        // integrate one step PAST THE END of it. A one-element write off the end of a
        // Kokkos allocation is silent in a Release build -- it lands in whatever follows
        // -- but it is undefined behaviour, a Kokkos Debug build aborts on it, and the
        // value written was never read. zarr(n) above still has to be set for every n.
        if (n+1 < N) {
          if (eos.IsGeneral()) {
            // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
            // branch is the same thing with rho = p/(Rgas T), kept in its original form.
            Real rho = DensFromPT(eos, Rgas, p, T);
            logparr(n+1) = logparr(n) + gz*dz*rho/p;
          } else {
            logparr(n+1) = logparr(n) + fac/T;
          }
        }
    }


    return;
}

KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p) {

//    Real Rgas = 4.593e7;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;
    Real dz = zarr(1)-zarr(0);
    Real T;

    if (z >= 0.0) {
        int Nt = std::floor(z/dz);
        Real logp = logparr(Nt) + (logparr(Nt+1)-logparr(Nt))/(zarr(Nt+1)-zarr(Nt))*(z-zarr(Nt));
        p = std::exp(logp);
        get_init_Tp(N, Tarr, lgparr, p,T);
        rho = DensFromPT(eos, Rgas, p, T);
    } else {
        Real T0;
        get_init_Tp(N, Tarr, lgparr, p0,T0);
        Real rho0 = DensFromPT(eos, Rgas, p0, T0);
//        Real grav_acc = -942.0;
        Real H0 = -p0/rho0/grav_acc;
        Real iH0 = 1.0/H0;
        p = p0 * std::exp(-z*iH0);
        rho = rho0 * std::exp(-z*iH0);
    }

    return;
}





KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real Teq = Tirr/sqrt(2);
    
    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);
    
    Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
    Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
    get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);
    
    Real R = gamir1/gamir2;
    Real gamp = gamir1 + gamir2 - SQR(gamir2)*R;
    taulim = sqrt(R/3.0)*sqrt(beta*SQR(R-1.0)-SQR(beta*(R-1.0))+R)/SQR(gamir1);
    Real At1 = SQR(gamir1)*log(1.0+1.0/(taulim*gamir1));
    Real At2 = SQR(gamir2)*log(1.0+1.0/(taulim*gamir2));
    
    Real a0 = 1.0/gamir1 + 1.0/gamir2;
    Real a1 = -1.0/(3.0*SQR(taulim))*(gamp/(1.0-gamp)*(gamir1+gamir2-2.0)/(gamir1+gamir2) + (gamir1+gamir2)*taulim - (At1+At2)*SQR(taulim));
    Real b0 = 1.0/(gamir1*gamir2/(gamir1-gamir2)*(At1-At2)/3.0 - SQR(gamir1*gamir2)/sqrt(3.0*gamp) - SQR(gamir1*gamir2)*gamir1*gamir2/(1.0-gamir1)/(1.0-gamir2)/(gamir1+gamir2));
    A = 1.0/3.0*(a0+a1*b0);
    B = -1.0/3.0*SQR(gamir1*gamir2)/gamp*b0;
    
//    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    
    for (int iv=0; iv<3; ++iv) {
        Real gamv;
        if (iv==0) gamv = gamv1;
        if (iv==1) gamv = gamv2;
        if (iv==2) gamv = gamv3;
        gamvv[iv] = gamv;
        
        Real longf = (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(gamir1+gamir2) - 3.0*gamv*(6.0*SQR(gamir1*gamir2)-SQR(gamv)*(SQR(gamir1)+SQR(gamir2)));
        Real a2 = SQR(taulim)/(gamp*SQR(gamv)) * longf/(1.0-SQR(gamv*taulim));
        Real Av1 = SQR(gamir1)*log(1.0+gamv/gamir1);
        Real Av2 = SQR(gamir2)*log(1.0+gamv/gamir2);
        Real a3 = -SQR(taulim) * (3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*(Av1+Av2) / (gamp*SQR(gamv)*gamv*(1.0-SQR(gamv*taulim)));
        Real b1 = gamir1*gamir2*(3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv))*SQR(taulim) / (gamp*SQR(gamv)*(SQR(gamv*taulim)-1.0));
        Real b2 = 3.0*(gamir1+gamir2)*SQR(gamv)*gamv / ((3.0*SQR(gamir1)-SQR(gamv))*(3.0*SQR(gamir2)-SQR(gamv)));
        Real b3 = (Av2-Av1) / (gamv*(gamir1-gamir2));
        
        C[iv] = -1.0/3.0*(b0*b1*(1.0+b2+b3)*a1 + a2 + a3);
        D[iv] = -B*b1*(1.0+b2+b3);
        E[iv] = (3.0-SQR(gamv/gamir1))*(3.0-SQR(gamv/gamir2)) / (9.0*gamv*(SQR(gamv*taulim)-1.0));
//        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C + D*exp(-tau/taulim) + E*exp(-gamv*tau));
    }
    
//    T = sqrt(sqrt(T4));
    return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real taulim, const Real &A, const Real &B, const Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamv)[3], const Real &tau, Real &T) {
    
    Real Tirr4 = SQR(SQR(Tirr));
    Real Tint4 = SQR(SQR(Tint));
    Real T4 = 3.0/4.0*Tint4*(tau + A + B*exp(-tau/taulim));
    for (int iv=0; iv<3; ++iv) {
        T4 += 3.0/4.0*1.0/3.0*Tirr4*mus*(C[iv] + D[iv]*exp(-tau/taulim) + E[iv]*exp(-gamv[iv]*tau));
    }
    T = sqrt(sqrt(T4));
    return;
}

template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr) {
    Real bar = 1.0e6;
    Real tautop = 1.0e-6;
    Real Ttop;
    Real taulim, A, B, C[3], D[3], E[3], gamv[3];
    get_picket_fence_Ttau_coeff(Tint, Tirr, met, grav, mus, taulim, A, B, C, D, E, gamv);
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    Real kapr, ptop;
    ptop = 2.0; // 2e-6 bar
    get_kapr(Ttop, ptop, met, kapr);
    tautop = kapr/grav*ptop;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tautop, Ttop);
    
    Real lgptop = log10(ptop);
    Real pbot = 300.0*bar;
    Real lgpbot = log10(pbot);
    Real dlgp = (lgpbot-lgptop)/(N-1);
    Real lgp = lgptop;
    Real p = pow(10.0,lgp);
    Real tau = tautop;
    Real lgtau = log10(tau);
    Real T;
    get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
    get_kapr(Ttop, p, met, kapr);
    for (int ip=0; ip<N; ++ ip) {
        Tarr(ip) = T;
        lgparr(ip) = lgp;
        Real K = p*kapr/grav*log(10.0);
        tau += K*dlgp;
        lgp += dlgp;
        p = pow(10.0,lgp);
//        lgp += dlgp;
//        Real pp = pow(10.0,lgp);
//        Real dp = pp-p;
//        Real K = kapr/grav;
//        tau += K*dp;
//        p = pp;
        get_picket_fence_Ttau(Tint, Tirr, mus, taulim, A, B, C, D, E, gamv, tau, T);
        get_kapr(T, p, met, kapr);
    }
    
    // Snapshot before the convective adjustment so the dump can show what it changed.
    std::vector<double> Tpre;
    if (!ad_dump_file.empty()) {
      Tpre.resize(N);
      for (int ip=0; ip<N; ++ip) Tpre[ip] = Tarr(ip);
    }

    adjust_ad_pT_arr(eos, Rgas, gamma, N, Tarr, lgparr);

    if (!ad_dump_file.empty() && global_variable::my_rank == 0) {
      // index 0 is the TOP (lowest pressure); gradients are forward differences in
      // log10 p, which is uniform here.  nabla > nabla_ad means convectively unstable.
      std::ofstream f(ad_dump_file);
      f.precision(10);
      f << std::scientific;
      f << "# deep_hot_jupiter_rt initial (p,T) profile and convective stability\n"
        << "# N = " << N << ", index 0 = top.  adjust_ad_pT_arr uses 0.9*grad_ad.\n"
        << "# nabla > nabla_ad_used  =>  CONVECTIVELY UNSTABLE at that level\n"
        << "# ip  p[bar]  T_before[K]  T_after[K]  nabla_before  nabla_after  "
        << "grad_ad  nabla_ad_used  unstable_before  unstable_after\n"
        << "# 'unstable' means nabla exceeds nabla_ad_used by more than 1 % of grad_ad\n";
      const double bar = 1.0e6;
      int nub = 0, nua = 0;
      for (int ip=0; ip<N-1; ++ip) {
        const double dlgp = lgparr(ip+1) - lgparr(ip);
        const double nb = (log10(Tpre[ip+1]) - log10(Tpre[ip]))/dlgp;
        const double na = (log10(Tarr(ip+1)) - log10(Tarr(ip)))/dlgp;
        const double pp = pow(10.0, lgparr(ip));
        const double gad = GradAd(eos, gamma, Rgas, pp, Tarr(ip));
        const double gus = 0.9*gad;
        // TOLERANCE, not a bare `>`: the adjustment sets every level it touched to
        // exactly 0.9*grad_ad, so an exact comparison reports round-off on all of them
        // as instability -- 1169 spurious levels of 9999 on the shipped setup. 1 % of
        // grad_ad is far below any real super-adiabaticity (the genuine band before the
        // adjustment exceeds it by up to 0.25) and far above the noise.
        const double tol = 0.01*gad;
        const int ub = (nb - gus > tol) ? 1 : 0;
        const int ua = (na - gus > tol) ? 1 : 0;
        nub += ub; nua += ua;
        f << ip << " " << pp/bar << " " << Tpre[ip] << " " << Tarr(ip) << " "
          << nb << " " << na << " " << gad << " " << gus << " "
          << ub << " " << ua << "\n";
      }
      f.close();
      std::cout << "deep_hot_jupiter_rt: wrote initial-profile stability dump to '"
                << ad_dump_file << "'; convectively unstable levels: "
                << nub << " before the adjustment, " << nua << " after (of " << N-1
                << ")" << std::endl;
    }

    return;
}

template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr) {

  // --- find convective boundary (search from bottom) ---
  int ic = -1;

  for (int ip = 0; ip < N-1; ++ip) {
    int i_inv = N - 1 - ip;
    int i_inv_p1 = N - 1 - (ip + 1);

    Real lgT  = log10(Tarr(i_inv));
    Real lgT1 = log10(Tarr(i_inv_p1));

    Real lgp  = lgparr(i_inv);
    Real lgp1 = lgparr(i_inv_p1);

    Real nabla = (lgT - lgT1) / (lgp - lgp1);

    Real T = Tarr(i_inv);
    Real nabla_ad = 0.9*GradAd(eos, gamma, Rgas, pow(10.0,lgp), T);

    if (nabla < nabla_ad) {
      ic = i_inv_p1;  // map back to original indexing
      break;
    }
  }

  // fallback if no crossing found
  if (ic > 0) {
    // --- enforce adiabat downward ---
    for (int ip = ic; ip < N-1; ++ip) {
      Real T = Tarr(ip);
      Real nabla_ad = 0.9*GradAd(eos, gamma, Rgas, pow(10.0,lgparr(ip)), T);

      Tarr(ip+1) = pow(10.0,
        nabla_ad * (lgparr(ip+1) - lgparr(ip)) + log10(T)
      );
    }
  }
    
  return;
}


