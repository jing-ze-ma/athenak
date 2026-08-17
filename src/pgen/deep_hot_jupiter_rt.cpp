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

// Athena++ headers
#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "srcterms/srcterms.hpp"
#include "utils/random.hpp"
#include "pgen.hpp"
#include "pgen_eos_utils.hpp"
#include "diffusion/resistivity.hpp"

#include <Kokkos_Random.hpp>

// EOS-aware conversions shared with the other stratified problem generators
using pgen_eos::EintFromP;
using pgen_eos::PresFromEint;
using pgen_eos::TempKelvin;
using pgen_eos::PresTempFromEint;
using pgen_eos::DensFromPT;
using pgen_eos::GradAd;

void HydrostaticEquilibrium(Mesh *pm);
void SourceFunc(Mesh *pm, Real bdt);

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt);
void double_gray_two_stream_RT(Mesh *pm, Real bdt);

void picket_fence_two_stream_RT(Mesh *pm, Real bdt);
KOKKOS_INLINE_FUNCTION
void get_albedo(const Real &Teff, const Real &gg, Real &A);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2, Real &gamv3, Real &beta, Real &gamir1, Real &gamir2);
KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr);
KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau_coeff(const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, Real &taulim, Real &A, Real &B, Real (&C)[3], Real (&D)[3], Real (&E)[3], Real (&gamvv)[3]);
KOKKOS_INLINE_FUNCTION
void get_picket_fence_Ttau(const Real &Tint, const Real &Tirr, const Real &mus, const Real &taulim, const Real &A, const Real &B, const Real (&C)[3], const Real (&D)[3], const Real (&E)[3], const Real (&gamv)[3], const Real &tau, Real &T);
template <typename View1D>
void get_picket_fence_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const Real &Tint, const Real &Tirr, const Real &met, const Real &grav, const Real &mus, const int &N, View1D Tarr, View1D lgparr);
template <typename View1D>
void adjust_ad_pT_arr(const EOS_Data &eos, const Real &Rgas, const Real &gamma, const int &N, View1D Tarr, View1D lgparr);

KOKKOS_INLINE_FUNCTION
void get_daynight_Tp(const Real &p, Real &Tn, Real &Td);
KOKKOS_INLINE_FUNCTION
void get_wb_Tp(const Real &p, Real &T);
KOKKOS_INLINE_FUNCTION
void get_init_Tp(const int &N, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_init_Tp_host(const int &N, const View1D &Tarr, const View1D &lgparr, const Real &p, Real &T);
template <typename View1D>
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_wb_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);
template <typename View1D>
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const View1D &Tarr, const View1D &lgparr, const int &N, const Real &zmax, View1D zarr, View1D logparr);
KOKKOS_INLINE_FUNCTION
void get_init_eos(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const DvceArray1D<Real> &Tarr, const DvceArray1D<Real> &lgparr, const int &N, const DvceArray1D<Real> &zarr, const DvceArray1D<Real> &logparr, const Real &z, Real &rho, Real &p);


namespace {
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
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//  \brief Problem Generator for the shallow hot Jupiter test

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  bool use_etotgrav = false;
  bool use_wellbalance_static = false;
  bool use_wellbalance_dynamic = false;
  const bool use_spherical_polar = pmy_mesh_->use_spherical_polar;
  bool user_srcs = pin->GetOrAddBoolean("problem","user_srcs",false);
  if (user_srcs) user_srcs_func = SourceFunc;
  // read before anything restart-sensitive: the outer BC needs it on restarts too
  bc_outer_maxwell = pin->GetOrAddBoolean("problem","bc_outer_maxwell",true);
  if (global_variable::my_rank == 0) {
    std::cout << "deep_hot_jupiter_rt: outer-x1 Maxwell-stress term in the ghost "
              << "extrapolation is " << (bc_outer_maxwell ? "ON" : "off") << std::endl;
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
    Real omega = pin->GetReal("problem","omega");
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
    get_wb_eos_arr(eos, Rgas, grav_acc, N, (r1-r0)*1.1, zarr.h_view, logparr.h_view);
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
    get_init_eos_arr(eos, Rgas, grav_acc, Tarr_init.h_view, lgparr_init.h_view, N, (r1-r0)*1.1, zarr_init.h_view, logparr_init.h_view);
    
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
        x2v = CellCenterX(j-js, nx2, x2min, x2max);
        x3v = CellCenterX(k-ks, nx3, x3min, x3max);
      }
      Real r = x1v;
      if (use_spherical_polar) x1v -= ap;
        
      Real lam, theta, phi;
      if (use_spherical_polar) {
        theta = x2v;
        lam = -x2v+M_PI/2.0;
        phi = x3v-M_PI;
      } else {
        theta = -(x2v*iap-M_PI/2.0);
        lam = x2v*iap;
        phi = x1v*iap;
      }
        
      Real pwb, denwb;
      get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
      Real p, den;
      get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,x1v,den,p);
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
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IVX,k,j,i) = 0.0;
      w0_(m,IEN,k,j,i) = EintFromP(eos, igm1, den, p);
        
        Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
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
          x2v = CellCenterX(j-js, nx2, x2min, x2max);
          x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        }
        Real r = x1v;
        if (use_spherical_polar) x1v -= ap;
          
        Real lam, theta, phi;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -x2v+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          theta = -(x3v*iap-M_PI/2.0);
          lam = x3v*iap;
          phi = x2v*iap;
        }
          
        Real pwb, denwb;
        get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
        
        Real p, den;
        get_init_eos(eos, Rgas, grav_acc, Tarr_init.d_view,lgparr_init.d_view,N,zarr_init.d_view,logparr_init.d_view,x1v,den,p);
//        get_init_eos_arr(lam, phi, N, zinitarr, logpinitarr);
//        get_init_eos(zinitarr,logpinitarr,x1v,lam,phi,den,p);
//        p = pwb;
//        den = denwb;
        
        Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
        
      if (use_etotgrav || use_wellbalance_dynamic) {
          if (use_spherical_polar) {
            x1v = x1f_(m,i);
            x2v = x2v_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = LeftEdgeX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
          phi0_x1f(m,k,j,i) = phicc;
          if (i == ie) {
              if (use_spherical_polar) {
                x1v = x1f_(m,i+1);
              } else {
                x1v = LeftEdgeX(i+1-is, nx1, x1min, x1max);
              }
              r = x1v;
              if (use_spherical_polar) x1v -= ap;
              if (use_spherical_polar) {
                theta = x2v;
                lam = -x2v+M_PI/2.0;
                phi = x3v-M_PI;
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x1f(m,k,j,i+1) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2f_(m,j);
            x3v = x3v_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
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
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
              phi0_x2f(m,k,j+1,i) = phicc;
          }
          
          if (use_spherical_polar) {
            x1v = x1v_(m,i);
            x2v = x2v_(m,j);
            x3v = x3f_(m,k);
          } else {
            x1v = CellCenterX(i-is, nx1, x1min, x1max);
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
          }
          r = x1v;
          if (use_spherical_polar) x1v -= ap;
          if (use_spherical_polar) {
            theta = x2v;
            lam = -x2v+M_PI/2.0;
            phi = x3v-M_PI;
          } else {
            theta = -(x3v*iap-M_PI/2.0);
            lam = x3v*iap;
            phi = x2v*iap;
          }
          phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
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
              } else {
                theta = -(x3v*iap-M_PI/2.0);
                lam = x3v*iap;
                phi = x2v*iap;
              }
              phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
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
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
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
                }
                if (use_spherical_polar) x1v -= ap;
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
              x2v = LeftEdgeX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
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
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = LeftEdgeX(k-ks, nx3, x3min, x3max);
            }
            if (use_spherical_polar) x1v -= ap;
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
              x2v = CellCenterX(j-js, nx2, x2min, x2max);
              x3v = CellCenterX(k-ks, nx3, x3min, x3max);
            }
            Real r = x1v;
            if (use_spherical_polar) x1v -= ap;
              
            Real lam, theta, phi;
            if (use_spherical_polar) {
              theta = x2v;
              lam = -x2v+M_PI/2.0;
              phi = x3v-M_PI;
            } else {
              theta = -(x3v*iap-M_PI/2.0);
              lam = x3v*iap;
              phi = x2v*iap;
            }
            
            Real phicc = - grav_acc * x1v;// - 0.5*SQR(omega*r*sin(theta));
            phicc0(m,k,j,i) = phicc;
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
            x2v = CellCenterX(j-js, nx2, x2min, x2max);
            x3v = CellCenterX(k-ks, nx3, x3min, x3max);
          }
          Real r = x1v;
          if (use_spherical_polar) x1v -= ap;
            
          Real pwb, denwb;
          get_wb_eos(eos, Rgas, grav_acc, zarr.d_view,logparr.d_view,x1v,denwb,pwb);
          u0wb(m,IDN,k,j,i) = denwb;
          u0wb(m,IM1,k,j,i) = 0.0;
          u0wb(m,IM2,k,j,i) = 0.0;
          u0wb(m,IM3,k,j,i) = 0.0;
          u0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
          if (use_etotgrav) {
              Real phicc = - grav_acc * x1v;
              u0wb(m,IEN,k,j,i) += denwb*phicc;
          }
          w0wb(m,IDN,k,j,i) = denwb;
          w0wb(m,IM1,k,j,i) = 0.0;
          w0wb(m,IM2,k,j,i) = 0.0;
          w0wb(m,IM3,k,j,i) = 0.0;
          w0wb(m,IEN,k,j,i) = EintFromP(eos, igm1, denwb, pwb);
        });
    }

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
          }
      });
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
        } else {
          bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i) + b0.x1f(m,k,j,i+1));
          bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i) + b0.x2f(m,k,j+1,i));
          bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i) + b0.x3f(m,k+1,j,i));
        }
      });
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
//                u0_(m,IEN,k,j,(ie+i+1)) -=  0.5*(SQR(bcc0(m,IBX,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
              Real lw, rw;
              lw = (x1f_(m,(ie+i+1)+1)-x1v_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
              rw = (x1v_(m,(ie+i+1))-x1f_(m,(ie+i+1)))/(x1f_(m,(ie+i+1)+1)-x1f_(m,(ie+i+1)));
              bcc0(m,IBX,k,j,(ie+i+1)) = lw*b0_x1f(m,k,j,(ie+i+1)) + rw*b0_x1f(m,k,j,(ie+i+1)+1);
              lw = (x2f_(m,j+1)-x2v_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              rw = (x2v_(m,j)-x2f_(m,j))/(x2f_(m,j+1)-x2f_(m,j));
              bcc0(m,IBY,k,j,(ie+i+1)) = lw*b0_x2f(m,k,j,(ie+i+1)) + rw*b0_x2f(m,k,j+1,(ie+i+1));
              lw = (x3f_(m,k+1)-x3v_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              rw = (x3v_(m,k)-x3f_(m,k))/(x3f_(m,k+1)-x3f_(m,k));
              bcc0(m,IBZ,k,j,(ie+i+1)) = lw*b0_x3f(m,k,j,(ie+i+1)) + rw*b0_x3f(m,k+1,j,(ie+i+1));
//                u0_(m,IEN,k,j,(ie+i+1)) +=  0.5*(SQR(bcc0(m,IBX,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));

              // The cell-centred ghost field above is always needed. What follows is the
              // Maxwell-stress correction to the hydrostatic extrapolation, which is
              // optional -- see bc_outer_maxwell.
              if (bc_outer_maxwell_) {
              Real pb = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+1)))+SQR(bcc0(m,IBY,k,j,(ie+i+1)))+SQR(bcc0(m,IBZ,k,j,(ie+i+1))));
              Real pbp1 = 0.5*(SQR(b0_x1f(m,k,j,(ie+i+2)))+SQR(bcc0(m,IBY,k,j,(ie+i+2)))+SQR(bcc0(m,IBZ,k,j,(ie+i+2))));
              Real M11 = pb - SQR(b0_x1f(m,k,j,(ie+i+1)));
              Real M11p1 = pbp1 - SQR(b0_x1f(m,k,j,(ie+i+2)));
              Real M12 = - b0_x2f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M12p1 = - b0_x2f(m,k,j+1,(ie+i+1)) * bcc0(m,IBX,k,j+1,(ie+i+1));
              Real M13 = - b0_x3f(m,k,j,(ie+i+1)) * bcc0(m,IBX,k,j,(ie+i+1));
              Real M13p1 = - b0_x3f(m,k+1,j,(ie+i+1)) * bcc0(m,IBX,k+1,j,(ie+i+1));
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
              Real gmag = dM1mag/(rho_i*fabs(grav_acc));
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
                        wtemp_(m,k,j,ie));
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
            u0_(m,IEN,k,j,(ie+i+1)) = e0_ip + 0.5*(SQR(u0_(m,IM1,k,j,(ie+i+1)))+SQR(u0_(m,IM2,k,j,(ie+i+1)))+SQR(u0_(m,IM3,k,j,(ie+i+1))))/rho0_ip;
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


void SourceFunc(Mesh *pm, Real bdt) {
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
    Real omega = pm->pgen->hot_jupiter_param.omega;
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
          Real &x1min = size.d_view(m).x1min;
          Real &x1max = size.d_view(m).x1max;
          Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
        
          Real &x2min = size.d_view(m).x2min;
          Real &x2max = size.d_view(m).x2max;
          Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
        
          Real &x3min = size.d_view(m).x3min;
          Real &x3max = size.d_view(m).x3max;
          Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
        }
        
        Real lam, phi, z, theta, r;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
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
        Real src = bdt*grav_acc*w0(m,IDN,k,j,i);
        if (!use_etotgrav) {
            u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
        }
        if (use_wellbalance_static) {
            src = bdt*grav_acc*(w0(m,IDN,k,j,i)-w0wb(m,IDN,k,j,i));
        }
        if (use_wellbalance_dynamic) {
          // The source is the background's own pressure difference across the cell, so it
          // has to come from the same entry point the reconstruction uses: getWBq0, which
          // dispatches to the ideal-gas closed forms or to the general-EOS background and
          // returns pressure directly, not an energy to be multiplied by (gamma-1).
          Real pl,pr,dum1,dum2,dum3;
          if (pmbp->phydro != nullptr) {
              pmbp->phydro->getWBq0(eos, WBVar::wb_pres,
                w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
                w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
                phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
                dum1,pl,dum2,pr,dum3);
          } else if (pmbp->pmhd != nullptr) {
              pmbp->pmhd->getWBq0(eos, WBVar::wb_pres,
                w0(m,IDN,k,j,i-1),w0(m,IDN,k,j,i),w0(m,IDN,k,j,i+1),
                w0(m,IEN,k,j,i-1),w0(m,IEN,k,j,i),w0(m,IEN,k,j,i+1),
                phicc0(m,k,j,i-1),phi0_x1f(m,k,j,i),phicc0(m,k,j,i),phi0_x1f(m,k,j,i+1),phicc0(m,k,j,i+1),
                dum1,pl,dum2,pr,dum3);
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
          u0(m,IM1,k,j,i) += rho*(cor+oor)*sine*bdt;
//          if (!use_etotgrav)
          u0(m,IEN,k,j,i) += rho*oor*(vr*sine+vtheta*cosine)*bdt;
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
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
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
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
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
          u0(m,IEN,k,j,i) += src*bdt;
        }
        
    });
    
    return;
}

void double_gray_two_stream_RT_source(Mesh *pm, Real bdt) {
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
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
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
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
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
void get_wb_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const int &N, const Real &zmax, View1D zarr, View1D logparr) {
    
//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_wb_Tp(p,T);

        if (eos.IsGeneral()) {
          // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
          // branch is the same thing with rho = p/(Rgas T), kept in its original form.
          Real rho = DensFromPT(eos, Rgas, p, T);
          logparr(n+1) = logparr(n) + grav_acc*dz*rho/p;
        } else {
          logparr(n+1) = logparr(n) + fac/T;
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
void get_init_eos_arr(const EOS_Data &eos, const Real &Rgas, const Real &grav_acc, const View1D &Tarr, const View1D &lgparr, const int &N, const Real &zmax, View1D zarr, View1D logparr) {

//    Real Rgas = 4.593e7;
//    Real grav_acc = -942.0;
    Real bar = 1.0e6;
    Real p0 = 250.0*bar;

    Real zmin = 0.0;
//    Real zmax = 1.2e9;
    Real dz = (zmax-zmin)/N;
    Real fac = grav_acc/Rgas*dz;
    logparr(0) = std::log(p0);

    for(int n=0; n<N; n++) {
        Real T;
        Real p = exp(logparr(n));
        zarr(n) = zmin + n*dz;

        get_init_Tp_host(N, Tarr, lgparr, p, T);

        if (eos.IsGeneral()) {
          // dln p/dz = rho g/p, closed with the EOS's (p,T) -> rho inversion. The ideal
          // branch is the same thing with rho = p/(Rgas T), kept in its original form.
          Real rho = DensFromPT(eos, Rgas, p, T);
          logparr(n+1) = logparr(n) + grav_acc*dz*rho/p;
        } else {
          logparr(n+1) = logparr(n) + fac/T;
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
void get_albedo(const Real &Teff, const Real &gg, Real &A) {
  // Parmentier+2015
  Real X = log10(Teff);
  Real g = gg*0.01;
  Real a,b;
  if (Teff < 250.0) {
    a = -0.335*pow(g,0.07);
    b = 0.0;
  } else if (Teff < 750.0) {
    a = -0.335*pow(g,0.07) + 2.149*pow(g,0.135);
    b = -0.896*pow(g,0.135);
  } else if (Teff < 1250.0) {
    a = -0.335*pow(g,0.07) - 0.428*pow(g,0.135);
    b = 0.0;
  } else {
    a = 16.947 - 3.174*pow(g,0.07) - 4.051*pow(g,0.135);
    b = -5.472 + 0.917*pow(g,0.07) + 1.170*pow(g,0.135);
  }
  Real log10A = a + b*X;
  A = pow(10.0,log10A);
  return;
}

KOKKOS_INLINE_FUNCTION
void get_picket_fence_coeff(const Real &Teq, const Real &Teff, Real &gamv1, Real &gamv2, Real &gamv3, Real &beta, Real &gamir1, Real &gamir2) {
  // Parmentier & Giollot 2014; Parmentier+2015; Roth+2024
  Real X = log10(Teff);
  Real a3, a2, a1, b3, b2, b1, ab, bb;
  
  Real ap = -2.36;
  Real bp = 13.92;
  Real cp = -19.38;
  if (Teff >= 1400.0 && Teq < 1800.0) {
    ap = -12.45;
    bp = 82.25;
    cp = -134.42;
  }
  
  if (Teff < 2000.0) {
    ab = 0.84;
    bb = 0.0;
  } else {
    ab = 6.21;
    bb = -1.63;
  }
  if (Teff >= 1400.0 && Teq < 1800.0) {
    ab = 3.0;
    bb = -0.69;
  }
  
  if (Teff < 200.0) {
    a3 = -3.03;
    b3 = -0.2;
    a2 = -7.37;
    b2 = 2.53;
    a1 = -5.51;
    b1 = 2.48;
  } else if (Teff < 300.0) {
    a3 = -13.87;
    b3 = 4.51;
    a2 = 13.99;
    b2 = -6.75;
    a1 = 1.23;
    b1 = -0.45;
  } else if (Teff < 600.0) {
    a3 = -11.95;
    b3 = 3.74;
    a2 = -15.18;
    b2 = 5.02;
    a1 = 8.65;
    b1 = -3.45;
  } else if (Teff < 1400.0) {
    a3 = -6.97;
    b3 = 1.94;
    a2 = -10.41;
    b2 = 3.31;
    a1 = -12.96;
    b1 = 4.33;
  } else if (Teff < 2000.0) {
    a3 = -3.65;
    b3 = 0.89;
    a2 = -19.95;
    b2 = 6.34;
    a1 = -23.75;
    b1 = 7.76;
  } else {
    a3 = -6.02;
    b3 = 1.61;
    a2 = 13.56;
    b2 = -3.81;
    a1 = 12.65;
    b1 = -3.27;
  }
  if (Teff >= 1400.0 && Teq < 1800.0) {
    if (Teff < 2000.0) {
      a3 = 0.02;
      b3 = -0.28;
      a2 = 6.96;
      b2 = -2.21;
      a1 = -1.68;
      b1 = 0.75;
    } else {
      a3 = -16.54;
      b3 = 4.74;
      a2 = -2.4;
      b2 = 0.62;
      a1 = 10.37;
      b1 = -2.91;
    }
  }
  Real log10gamv1 = a1 + b1*X;
  Real log10gamv2 = a2 + b2*X;
  Real log10gamv3 = a3 + b3*X;
  Real log10gamp = ap*SQR(X) + bp*X + cp;
  beta = ab + bb*X;
    
  gamv1 = pow(10.0,log10gamv1);
  gamv2 = pow(10.0,log10gamv2);
  gamv3 = pow(10.0,log10gamv3);
  Real gamp = pow(10.0,log10gamp);
  Real dum = (gamp-1.0)/(2.0*beta*(1.0-beta));
  Real R = 1.0 + dum + sqrt(SQR(dum)+dum);
  gamir1 = beta + R - beta*R;
  gamir2 = gamir1/R;
    
  return;
}

KOKKOS_INLINE_FUNCTION
void get_kapr(const Real &T, const Real &p, const Real &met, Real &kapr) {
  // Freedman+2014
  Real T1 = T;
  Real p1 = p;
  if (T > 4000.0) T1 = 4000.0;
  if (T < 75.0) T1 = 75.0;
  if (p > 3.0e8) p1 = 3.0e8;
  if (p < 1.0) p1 = 1.0;
  Real lgT = log10(T1);
  Real lgp = log10(p1);
  Real c1 = 10.602;
  Real c2 = 2.882;
  Real c3 = 6.09e-15;
  Real c4 = 2.954;
  Real c5 = -2.526;
  Real c6 = 0.843;
  Real c7 = -5.490;
  Real c8, c9, c10, c11, c12, c13;
  if (T1 < 800.0) {
    c8 = -14.051;
    c9 = 3.055;
    c10 = 0.024;
    c11 = 1.877;
    c12 = -0.445;
    c13 = 0.8321;
  } else {
    c8 = 82.241;
    c9 = -55.456;
    c10 = 8.754;
    c11 = 0.7048;
    c12 = -0.0414;
    c13 = 0.8321;
  }
    
  Real lgkl = c1*atan(lgT-c2) - c3/(lgp+c4)*exp(SQR(lgT-c5)) + c6*met + c7;
  Real lgkh = c8 + c9*lgT + c10*SQR(lgT) + lgp*(c11+c12*lgT) + c13*met*(0.5+1.0/M_PI*atan((lgT-2.5)/0.2));
  Real kl = pow(10.0,lgkl);
  Real kh = pow(10.0,lgkh);
  kapr = kl + kh;
  return;
}

KOKKOS_INLINE_FUNCTION
void get_Tint(const Real &Teq, Real &Tint) {
  // Thorngren+2019 +Erratum
  Real boltz_sigma = 5.6704e-5;
  Real F = 4.0*boltz_sigma*SQR(SQR(Teq));
  Tint = 0.39*Teq*exp(-SQR(log10(F)-9.0-0.14)/1.095);
  Tint = (Tint < 100.0) ? 100.0 : Tint;
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
    
    adjust_ad_pT_arr(eos, Rgas, gamma, N, Tarr, lgparr);
    
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

void picket_fence_two_stream_RT(Mesh *pm, Real bdt) {
    // Noti+2023; Lee+2021
    
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
//    Real met = 0.0;
    
//    ParameterInput* pin;
    Real grav = pm->pgen->hot_jupiter_param.grav;
    Real ap = pm->pgen->hot_jupiter_param.ap;
    Real Rgas = pm->pgen->hot_jupiter_param.Rgas;
    Real Teq = pm->pgen->hot_jupiter_param.Teq;
    Real met = pm->pgen->hot_jupiter_param.met;
    
    Real iap = 1.0/ap;
    Real gm1 = gamma-1.0;
    Real igm1 = 1.0/gm1;
    Real cgs2Pa = 0.1;
    Real boltz_sigma = 5.6704e-5;
    
    Real Tirr = Teq*sqrt(2);
    Real Tirr4 = SQR(SQR(Tirr));
    Real Fstar = boltz_sigma*Tirr4;
    Real Tint;
    get_Tint(Teq, Tint);
    Real Tint4 = SQR(SQR(Tint));
    Real Iint = boltz_sigma/M_PI*Tint4;
    
    Real mug[2];
    Real wg[2];
    mug[0] = 0.21132487;
    mug[1] = 0.78867513;
    wg[0] = 0.5;
    wg[1] = 0.5;
    
    Real Teff0 = sqrt(sqrt(Tint4+Tirr4/sqrt(3.0)));
    Real albedo;
    get_albedo(Teff0,grav,albedo);

//    size_t scr_size = 8 * ScrArray1D<Real>::shmem_size(n1);
//    int scr_level = 0;
//    par_for_outer("2stream_rt", DevExeSpace(), scr_size, scr_level,
//                  0, nmb1, ks, ke, js, je,
//    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    par_for("2stream_rt", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
//        ScrArray1D<Real> tau_down_r_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_v_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> B(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_down_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> I_ir_up_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> F_ir_f(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> Q_v(member.team_scratch(scr_level), n1);
//        ScrArray1D<Real> kapJ_ir(member.team_scratch(scr_level), n1);
        constexpr int NN = 270;
        Real tau_down_r_f[NN];
//        Real F_v_down_f[NN];
        Real B[NN];
        Real I_ir_down_f[NN];
        Real I_ir_up_f[NN];
        Real F_ir_f[NN];
        Real Q_v[NN];
//        Real kapJ_ir[NN];
        
        Real x2v = x2v_(m,j);
        Real x3v = x3v_(m,k);
        
        Real rtop = x1v_(m,ie+1);
        Real rbot = x1v_(m,is);
        
        Real lam, phi, theta;
        if (use_spherical_polar) {
          theta = x2v;
          lam = -theta+M_PI/2.0;
          phi = x3v-M_PI;
        } else {
          lam = x3v*iap;
          theta = -lam+M_PI/2.0;
          phi = x2v*iap;
        }
        Real ex = sin(theta)*cos(phi);
        Real ex0 = 1.0;
        Real mu0 = ex*ex0;
        if (test_oned) mu0 = cos(85.0/90.0*M_PI/2.0);
        
        Real mus = (mu0 > 0.0) ? mu0 : 0.0;
        Real Teff = sqrt(sqrt(Tint4+(1.0-albedo)*mus*Tirr4));
        Real gamv1, gamv2, gamv3, beta, gamir1, gamir2;
        get_picket_fence_coeff(Teq, Teff, gamv1, gamv2, gamv3, beta, gamir1, gamir2);
        
        // 3 V Bands
        // top
        Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,ie+1),w0(m,IEN,k,j,ie+1));
        Real rho = w0(m,IDN,k,j,ie+1);
        Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,ie+1),p);
        B[ie+1] = boltz_sigma/M_PI*SQR(SQR(T));
        Real kapr;
        get_kapr(T, p, met, kapr);
        Real tau_r_f = kapr*p/grav;
        tau_down_r_f[ie+1] = tau_r_f;
        Real drtop = tau_r_f/(kapr*rho);
        Real delta = drtop/rtop;
        Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
        fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
        Real tausl = tau_r_f*fac;
        Real trans1 = exp(-gamv1*tau_down_r_f[ie+1]*fac);
        Real trans2 = exp(-gamv2*tau_down_r_f[ie+1]*fac);
        Real trans3 = exp(-gamv3*tau_down_r_f[ie+1]*fac);
//        F_v_down_f[ie+1] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
//        F_v_down_f(ie+1) = (mu0 > 0.0)? F_v_down_f(ie+1) : 0.0;
        // down-sweep
        for (int i=ie; i>is-1; --i) {
          Real rho = w0(m,IDN,k,j,i);
          Real p, T;
          PresTempFromEint(eos,gm1,Rgas,rho,w0(m,IEN,k,j,i),wtemp_(m,k,j,i),p,T);
          B[i] = boltz_sigma/M_PI*SQR(SQR(T));
          Real kapr;
          get_kapr(T, p, met, kapr);
          Real dr = dx1(m,k,j,i);
          tau_down_r_f[i] = tau_down_r_f[i+1] + kapr*rho*dr;
          Real r = x1f_(m,i);
////          Real delta = (drtop+(rtop-r))/r;
//          Real delta = dr/r;
//          Real fac = (sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0)/delta;
//          tausl += kapr*rho*r*(sqrt(SQR(mu0)+2.0*delta+SQR(delta)) - mu0);
          Real fac = (mu0 > 0.1) ? (1.0/mu0) : (1.0/0.1);
          Real trans1 = exp(-gamv1*tau_down_r_f[i]*fac);
          Real trans2 = exp(-gamv2*tau_down_r_f[i]*fac);
          Real trans3 = exp(-gamv3*tau_down_r_f[i]*fac);
//          Real trans1 = exp(-gamv1*tausl);
//          Real trans2 = exp(-gamv2*tausl);
//          Real trans3 = exp(-gamv3*tausl);
//          F_v_down_f[i] = (1.0-albedo)*Fstar*mus*1.0/3.0*(trans1+trans2+trans3);
          Real mucr = 0.0; //sqrt(1.0-SQR(r0/r));
          Real Qv = kapr*rho*(1.0-albedo)*Fstar*1.0/3.0*(gamv1*trans1+gamv2*trans2+gamv3*trans3);
          Q_v[i] = (mu0 > -mucr) ? Qv : 0.0;
        }
        
        // 2 IR Bands
        for (int i=is; i<ie+2; ++i) {
          F_ir_f[i] = 0.0;
//          kapJ_ir[i] = 0.0;
        }
        // two quadrature
        for (int n=0; n<2; ++n) {
          Real mugg = mug[n];
          Real wgg = wg[n];
          for (int vir=0; vir<2; ++vir) {
            Real gamir, fb;
            if (vir == 0) {
              gamir = gamir1;
              fb = beta;
            } else {
              gamir = gamir2;
              fb = 1.0-beta;
            }

            // top
            Real dtauir = gamir*tau_down_r_f[ie+1];
            Real trans = exp(-dtauir/mugg);
            I_ir_down_f[ie+1] = (1.0-trans)*(fb*B[ie+1]);
            // down-sweep
            for (int i=ie; i>is-1; --i) {
              Real dtauir = gamir*(tau_down_r_f[i]-tau_down_r_f[i+1]);
//              Real trans = exp(-dtauir/mugg);
//              Real e0 = -expm1(-dtauir/mugg);
//              Real e1 = dtauir/mugg - e0;
//              Real alpa = 0.5*e0*(fb*B[i+1]+fb*B[i])/(fb*B[i+1]);
//              Real alp = (dtauir > 1.0e-3) ? (e0 - e1/(dtauir/mugg)) : alpa;
//              Real bet = (dtauir > 1.0e-3) ? (e1/(dtauir/mugg)) : 0.0;
              Real x = dtauir/mugg;
              Real e0 = -expm1(-x);
              Real alp = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              I_ir_down_f[i] = (1.0-e0)*I_ir_down_f[i+1] + alp*fb*B[i+1] + bet*fb*B[i];
            }
              
            // bottom
            I_ir_up_f[is] = Iint + I_ir_down_f[is];
            // up-sweep
            for (int i=is+1; i<ie+2; ++i) {
              Real dtauir = gamir*(tau_down_r_f[i-1]-tau_down_r_f[i]);
//              Real trans = exp(-dtauir/mugg);
//              Real e0 = -expm1(-dtauir/mugg);
//              Real e1 = dtauir/mugg - e0;
//              Real beto = 0.5*e0*(fb*B[i]+fb*B[i-1])/(fb*B[i]);
//              Real bet = (dtauir > 1.0e-3) ? (e1/(dtauir/mugg)) : beto;
//              Real gam = (dtauir > 1.0e-3) ? (e0 - e1/(dtauir/mugg)) : 0.0;
              Real x = dtauir/mugg;
              Real e0 = -expm1(-x);
              Real bet = (x > 1.0e-3) ? (1.0 - e0/x) : (x/2.0-SQR(x)/6.0);
              Real gam = (x > 1.0e-3) ? (e0 - 1.0 + e0/x) : (x/2.0-SQR(x)/3.0);
              I_ir_up_f[i] = (1.0-e0)*I_ir_up_f[i-1] + bet*fb*B[i] + gam*fb*B[i-1];
            }
              
            for (int i=is; i<ie+2; ++i) {
              Real F_ir_down_f = 2.0*M_PI*wgg*mugg*I_ir_down_f[i];
              Real F_ir_up_f = 2.0*M_PI*wgg*mugg*I_ir_up_f[i];
              F_ir_f[i] += F_ir_up_f - F_ir_down_f;
            }
//            for (int i=is; i<ie+1; ++i) {
//              Real J = wgg*I_ir_down_f[i] + wgg*I_ir_up_f[i];
//              Real p = w0(m,IEN,k,j,i)*gm1;
//              Real rho = w0(m,IDN,k,j,i);
//              Real T = p/Rgas/rho;
//              Real kapr;
//              get_kapr(T, p, met, kapr);
//              kapJ_ir[i] += gamir*kapr*J;
//            }
          }
        }
        
//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();
        
//        par_for_inner(member, is, ie, [&](const int i) {
        for (int i=is; i<ie+1; ++i) {
          // source term as flux divergence
          Real area_t = area1(m,k,j,i+1);
          Real area_b = area1(m,k,j,i);
          Real vol = volume(m,k,j,i);
            Real Ft = F_ir_f[i+1];//-F_v_down_f(i+1);
            Real Fb = F_ir_f[i];//-F_v_down_f(i);
          Real src = -(Ft-Fb)/dx1(m,k,j,i);
          if (correct_spherical) {
            src = -(Ft*area_t-Fb*area_b)/vol;
          }
            src += Q_v[i];
          Real du_flux = src*bdt;
            
//          // source term semi-implicit
//          Real p = PresFromEint(eos,gm1,w0(m,IDN,k,j,i),w0(m,IEN,k,j,i));
//          Real rho = w0(m,IDN,k,j,i);
//          Real T = TempKelvin(eos,Rgas,rho,w0(m,IEN,k,j,i),p);
//          Real kapr;
//          get_kapr(T, p, met, kapr);
//          Real cv = Rgas*rho*igm1;
//          Real e0 = eos.IsGeneral() ? w0(m,IEN,k,j,i) : cv*T;
//          Real kk = 0.0;
//          Real bb = du_flux + e0;
////          Real bb = Q_v(i)*bdt + e0;
//          for (int vir=0; vir<2; ++vir) {
//            Real gamir, fb;
//            if (vir == 0) {
//              gamir = gamir1;
//              fb = beta;
//            } else {
//              gamir = gamir2;
//              fb = 1.0-beta;
//            }
//            kk += -4.0*M_PI*gamir*kapr*rho*fb*boltz_sigma/M_PI*bdt;
//            bb += 4.0*M_PI*gamir*kapr*rho*fb*B[i]*bdt;
//          }
////          bb += 4.0*M_PI*rho*kapJ_ir(i)*bdt;
//          int ierr=0;
//          Real e = e0;
//          // Newton-Raphson. A general EOS has no e = c_v T with constant c_v, so the
//          // iteration runs on the internal energy directly rather than on T:
//          // F(e) = e - kk T(e)^4 - bb, with dT/de = temp_cgs/(d c_v) since T is in K.
//          for (int n=0; n<100; ++n) {
//            Real de;
//            if (eos.IsGeneral()) {
//              Real dTde = eos.temp_cgs/(rho*eos.SpecificHeatCv(rho,e));
//              de = e - kk*SQR(SQR(T)) - bb;
//              e -= de / (1.0 - 4.0*kk*T*T*T*dTde);
//              T = eos.Temperature(rho,e)*eos.temp_cgs;
//            } else {
//              e = cv*T;
//              de = e - kk*SQR(SQR(T)) - bb;
//              T -= de / (cv - 4.0*kk*T*T*T);
//            }
//            if (T < 0.0) {
//              e = e0;
//              ierr = 1;
//              break;
//            }
//            if (fabs(de) <= 1.0e-10*e)
//              break;
//          }
//          Real du_src = e-e0;
//
//          Real du = (fabs(du_flux) < e0 && ierr == 1) ? du_flux : du_src;
          Real du = du_flux;
          u0(m,IEN,k,j,i) += du;
        }
//        });
        
//        // Sync all threads in the team so that scratch memory is consistent
//        member.team_barrier();
        
    });
    
    return;
}

