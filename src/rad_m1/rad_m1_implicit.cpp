//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_implicit.cpp
//! \brief IMPLICIT transport along x1 columns, <rad_m1>/transport = implicit_x1
//! (milestone 3a of docs/dev/rad_m1_implicit_design.md).
//!
//! ONE backward-Euler solve per hydro step, at the TRUE speed of light, no sub-cycling
//! and no PD-ARS.  The normal flux lives on the x1 FACES and is eliminated there, which
//! turns the coupled (E, F) system into one tridiagonal M-matrix system for E per column:
//!
//!   F0'_f = theta_f [ F0^n_f - chat c dt (w_i E'_i - w_{i-1} E'_{i-1})/dx
//!                     - chat dt v_f g0_f ],    theta_f = 1/(1 + chat dt (rho k_t)_f)
//!   E'_i + (dt/dx)(chat/c)[ (F0' + A')_{i+1/2} - (F0' + A')_{i-1/2} ]
//!        = E^n_i + dt chat (rho kappa_P a T'^4 - rho kappa_E E0'_i)
//!
//! with w = P_11/E from the LAGGED closure (= chi in 1-D), A = a E the enthalpy flux
//! (a = v1 (1 + w)) upwinded with the face velocity, and the emission term linearised in
//! T and eliminated into the diagonal (design sect. 2).  Both the diffusion and the
//! upwind-advection off-diagonals are non-positive and the source adds
//! dt chat rho kappa_E rho c_v/B >= 0 to the diagonal, so the matrix is an M-matrix and
//! E' > 0 at any dt.
//!
//! OUTER LOOP: Picard on (w, theta, a, de0, g0, the upwind directions, optionally the
//! opacities) plus the safeguarded scalar root find for T' at fixed E'.  At convergence
//! the linearised emission correction vanishes identically, so the converged state
//! satisfies the NONLINEAR backward-Euler equations (see the note above
//! M1ImplTemperature).
//!
//! RESTRICTIONS of 3a, all fatal:
//!   * exactly ONE MeshBlock along x1 (each column is solved by a plain Thomas sweep
//!     inside one block; the partitioned line solve of two_stream_column_partition.hpp
//!     is NOT used);
//!   * nx2 = nx3 = 1 unless <rad_m1>/implicit_allow_multid = true, in which case the
//!     columns are INDEPENDENT: there is no transport along x2/x3 in this mode and
//!     F_2 = F_3 = 0 always;
//!   * no SMR/AMR.

#include <float.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {

namespace {
//----------------------------------------------------------------------------------------
//! \fn ImplBCFromString
//! \brief one x1 boundary type of the implicit solve, from the input string, or from the
//! mesh boundary flag when the input does not name one.

int ImplBCFromString(const std::string &s, const BoundaryFlag mbc) {
  if (s.compare("auto") == 0) {
    if (mbc == BoundaryFlag::periodic) return M1_IBC_PERIODIC;
    if (mbc == BoundaryFlag::reflect) return M1_IBC_REFLECT;
    return M1_IBC_MARSHAK;
  }
  if (s.compare("marshak") == 0) return M1_IBC_MARSHAK;
  if (s.compare("flux") == 0) return M1_IBC_FLUX;
  if (s.compare("reflect") == 0) return M1_IBC_REFLECT;
  if (s.compare("periodic") == 0) return M1_IBC_PERIODIC;
  if (s.compare("efix") == 0) return M1_IBC_EFIX;
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
    << std::endl << "<rad_m1>/implicit_bc_x1min|max = '" << s << "' is not a valid "
    << "choice (auto | marshak | flux | reflect | periodic | efix)" << std::endl;
  std::exit(EXIT_FAILURE);
  return M1_IBC_MARSHAK;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitInit
//! \brief read the <rad_m1> parameters of the implicit solver, check the restrictions of
//! 3a and allocate the face array and the work array.  Called at the END of the
//! constructor, from the delimited hook there, and a no-op in explicit mode.

void RadiationM1::ImplicitInit(ParameterInput *pin) {
  if (transport != M1_TRANSPORT_IMPLICIT_X1) return;

  impl_cfl = pin->GetOrAddReal("rad_m1","implicit_cfl",-1.0);
  impl_tol = pin->GetOrAddReal("rad_m1","implicit_tol",1.0e-8);
  impl_maxit = pin->GetOrAddInteger("rad_m1","implicit_maxit",30);
  impl_opac_update = pin->GetOrAddBoolean("rad_m1","implicit_opac_update",false);
  impl_allow_multid = pin->GetOrAddBoolean("rad_m1","implicit_allow_multid",false);
  marshak_q = pin->GetOrAddReal("rad_m1","marshak_q",0.5);
  if (!(impl_tol > 0.0) || impl_maxit < 1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/implicit_tol must be positive and implicit_maxit >= 1"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &mbcs = pmy_pack->pmesh->mesh_bcs;
  ibc_x1min = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1min","auto"),
      mbcs[static_cast<int>(BoundaryFace::inner_x1)]);
  ibc_x1max = ImplBCFromString(
      pin->GetOrAddString("rad_m1","implicit_bc_x1max","auto"),
      mbcs[static_cast<int>(BoundaryFace::outer_x1)]);
  iflux_x1min = pin->GetOrAddReal("rad_m1","implicit_flux_x1min",0.0);
  iflux_x1max = pin->GetOrAddReal("rad_m1","implicit_flux_x1max",0.0);
  iebath_x1min = pin->GetOrAddReal("rad_m1","implicit_ebath_x1min",0.0);
  iebath_x1max = pin->GetOrAddReal("rad_m1","implicit_ebath_x1max",0.0);
  if ((ibc_x1min == M1_IBC_PERIODIC) != (ibc_x1max == M1_IBC_PERIODIC)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> implicit x1 boundaries: periodic must be set on BOTH "
      << "ends or neither" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ---- the restrictions of 3a
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  auto &mindcs = pmy_pack->pmesh->mesh_indcs;
  if (mindcs.nx1 != indcs.nx1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 (milestone 3a) needs exactly ONE "
      << "MeshBlock along x1: <meshblock>/nx1 must equal <mesh>/nx1 (" << indcs.nx1
      << " vs " << mindcs.nx1 << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if ((indcs.nx2 > 1 || indcs.nx3 > 1) && !impl_allow_multid) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does NOT transport along x2/x3. "
      << "Set <rad_m1>/implicit_allow_multid = true to run a set of INDEPENDENT x1 "
      << "columns (F_2 = F_3 = 0 everywhere)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pmy_pack->pmesh->multilevel) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/transport = implicit_x1 does not support SMR/AMR"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // one solve per hydro step: no sub-cycling, one stage
  subcycle = false;
  nstage = 1;
  bool other_sets_dt = (pin->DoesBlockExist("hydro") || pin->DoesBlockExist("mhd") ||
                        pin->DoesBlockExist("z4c") || pin->DoesBlockExist("particles"));
  sets_mesh_dt = (impl_cfl > 0.0) || (!other_sets_dt);

  // ---- arrays
  int nmb = std::max((pmy_pack->nmb_thispack), (pmy_pack->pmesh->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(f0x1, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1, 0.0);
  Kokkos::realloc(f0x1n, nmb, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(f0x1n, 0.0);
  Kokkos::realloc(iw, nmb, M1_NIW, ncells3, ncells2, ncells1);
  Kokkos::deep_copy(iw, 0.0);

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: transport=implicit_x1 (milestone 3a) implicit_cfl="
              << impl_cfl << " tol=" << impl_tol << " maxit=" << impl_maxit
              << " opac_update=" << (impl_opac_update ? "true" : "false")
              << " marshak_q=" << marshak_q << std::endl;
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::SetImplicitX1BC
//! \brief let a problem generator name the x1 boundary types (and the imposed fluxes) of
//! the implicit solve, for a mesh whose x1 flags are `user` and whose own BC routine the
//! module cannot interpret.  Alternative to <rad_m1>/implicit_bc_x1min|max.

void RadiationM1::SetImplicitX1BC(int lo_type, Real lo_flux, int hi_type, Real hi_flux) {
  ibc_x1min = lo_type;
  iflux_x1min = lo_flux;
  ibc_x1max = hi_type;
  iflux_x1max = hi_flux;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitReport
//! \brief one line at the end of the run with the Picard statistics

void RadiationM1::ImplicitReport() {
  if (transport != M1_TRANSPORT_IMPLICIT_X1) return;
  // the Picard iteration count is MPI_MAX-reduced every step (see ImplicitSolve), so
  // every rank holds the same three numbers and no reduction is needed here
  if (global_variable::my_rank != 0) return;
  Real mean = (impl_nstep > 0.0) ? (impl_itsum/impl_nstep) : 0.0;
  std::cout << "<rad_m1> implicit transport: solves=" << impl_nstep
            << " Picard iterations mean=" << mean << " max=" << impl_itmax
            << " NON-CONVERGED=" << impl_nfail << std::endl;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::ImplicitSolve
//! \brief the whole backward-Euler step: the Picard loop, the tridiagonal column solves,
//! the write-back into u0 and into the gas.

TaskStatus RadiationM1::ImplicitSolve(Driver *pdrive, int stage) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto u0_ = u0;
  auto iw_ = iw;
  auto f0_ = f0x1;
  // F0^n, the face flux at the START of the step: the Picard loop overwrites f0x1 with
  // each new iterate, so the backward-Euler right-hand side needs its own copy
  auto f0n_ = f0x1n;
  Kokkos::deep_copy(DevExeSpace(), f0x1n, f0x1);
  auto opac_ = opac;
  auto &mbsize = pmy_pack->pmb->mb_size;
  Real cl = c_light;
  Real ch = chat;
  Real efl = e_floor;
  Real ar = arad;
  Real dt = dt_sub;
  bool edd = eddington;
  bool ovc = source_ovc;
  bool feedback = gas_feedback;
  // the DEBUG switches of the explicit coupling (dbg_gas_force / dbg_gas_heat /
  // opac_freeze) are not on this branch point; the implicit path runs the production
  // values unconditionally.
  const bool dbgf = true;
  const bool dbgh = true;
  bool fref = (force_ref == M1_FREF_WB_ARAD);
  auto aref_ = arad_ref;
  Real mq = marshak_q;
  int bclo = ibc_x1min, bchi = ibc_x1max;
  Real fxlo = iflux_x1min, fxhi = iflux_x1max;
  Real eblo = iebath_x1min, ebhi = iebath_x1max;
  bool cyclic = (bclo == M1_IBC_PERIODIC);

  const bool have_hydro = (pmy_pack->phydro != nullptr);
  const bool src_on = have_hydro && coupling && dbgh && !opac_zero;
  auto uh = have_hydro ? pmy_pack->phydro->u0 : u0;
  const bool etg = have_hydro ? pmy_pack->phydro->use_etotgrav : false;
  auto phicc = have_hydro ? pmy_pack->phydro->phicc0 : arad_ref;

  //-------------------------------------------------------------------------- start state
  par_for("m1_impl_i0", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = fmax(u0_(m,M1_E,k,j,i), efl);
    iw_(m,M1_IW_EN,k,j,i) = e;
    iw_(m,M1_IW_EP,k,j,i) = e;
    iw_(m,M1_IW_F1,k,j,i) = u0_(m,M1_F1,k,j,i);
    iw_(m,M1_IW_V1,k,j,i) = 0.0;
    iw_(m,M1_IW_SRCB,k,j,i) = 0.0;
    iw_(m,M1_IW_SRCR,k,j,i) = 0.0;
    iw_(m,M1_IW_DE0,k,j,i) = 0.0;
    iw_(m,M1_IW_G0,k,j,i) = 0.0;
    iw_(m,M1_IW_TP,k,j,i) = 0.0;
    iw_(m,M1_IW_EGN,k,j,i) = 0.0;
  });

  if (have_hydro) {
    auto eos = pmy_pack->phydro->peos->eos_data;
    par_for("m1_impl_i1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dd = uh(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      Real ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                       SQR(uh(m,IM3,k,j,i)))*idd;
      Real egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      Real eg = uh(m,IEN,k,j,i) - ekin - egrv;
      iw_(m,M1_IW_EGN,k,j,i) = eg;
      iw_(m,M1_IW_TP,k,j,i) = eos.Temperature(dd, fmax(eg, 1.0e-300));
      iw_(m,M1_IW_V1,k,j,i) = uh(m,IM1,k,j,i)*idd;
    });
  }

  //--------------------------------------------------------------------- the Picard loop
  int it = 0;
  Real resid = 0.0;
  bool converged = false;
  for (it = 0; it < impl_maxit && !converged; ++it) {
    // (a) optional opacity re-evaluation at the current temperature iterate
    if (impl_opac_update && it > 0 && have_hydro && !opac_zero) {
      int otype = opacity_type;
      Real kp = kappa_p, kev = kappa_e, kf = kappa_f, kscat = kappa_s;
      Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
      M1OpacTab ot = otab;
      par_for("m1_impl_opac", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real d = uh(m,IDN,k,j,i);
        Real t = iw_(m,M1_IW_TP,k,j,i);
        Real op, oe, of, os;
        if (otype == M1_OPAC_TABLE) {
          M1TableOpacities(ot, d, t, op, oe, of, os);
        } else {
          M1Opacities(otype, d, t, kp, kev, kf, kscat, rref, tref, aa, bb, op, oe,
                      of, os);
        }
        opac_(m,M1_OP_P,k,j,i) = d*op;
        opac_(m,M1_OP_E,k,j,i) = d*oe;
        opac_(m,M1_OP_T,k,j,i) = d*(of + os);
      });
    }

    // (b) the lagged closure, the enthalpy-flux coefficient, de0 and g0
    par_for("m1_impl_lag", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real e = fmax(iw_(m,M1_IW_EP,k,j,i), efl);
      Real f1 = iw_(m,M1_IW_F1,k,j,i);
      Real rf = f1/(cl*e);
      if (rf > 1.0) {rf = 1.0;}
      if (rf < -1.0) {rf = -1.0;}
      Real chi = edd ? (1.0/3.0) : M1Chi(fabs(rf));
      Real v1 = iw_(m,M1_IW_V1,k,j,i);
      iw_(m,M1_IW_WCHI,k,j,i) = chi;
      iw_(m,M1_IW_ADV,k,j,i) = v1*(1.0 + chi);
      // E0 - E with F_2 = F_3 = 0: P_11 = chi E, P_22 = P_33 = (1-chi) E/2
      Real b1 = v1/cl;
      Real de0 = ovc ? (-2.0*b1*f1/cl) : (b1*b1*e - 2.0*b1*f1/cl + b1*b1*chi*e);
      iw_(m,M1_IW_DE0,k,j,i) = de0;
      Real rkev = opac_(m,M1_OP_E,k,j,i);
      Real rkpv = opac_(m,M1_OP_P,k,j,i);
      Real tp = iw_(m,M1_IW_TP,k,j,i);
      Real t2 = tp*tp;
      iw_(m,M1_IW_G0,k,j,i) = rkev*(e + de0) - rkpv*ar*t2*t2;
    });

    // (c) the emission/absorption source, linearised in T about the iterate
    if (src_on) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      par_for("m1_impl_src", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real rkpv = opac_(m,M1_OP_P,k,j,i);
        Real rkev = opac_(m,M1_OP_E,k,j,i);
        if (rkpv == 0.0 && rkev == 0.0) {
          iw_(m,M1_IW_SRCB,k,j,i) = 0.0;
          iw_(m,M1_IW_SRCR,k,j,i) = 0.0;
          return;
        }
        Real dd = uh(m,IDN,k,j,i);
        Real tk = iw_(m,M1_IW_TP,k,j,i);
        Real ee, pp, cr, ct, cv;
        eos.ThermoAt(dd, tk, ee, pp, cr, ct, cv);
        Real t3 = tk*tk*tk;
        Real t4 = t3*tk;
        Real de0 = iw_(m,M1_IW_DE0,k,j,i);
        Real bk = dd*cv + 4.0*cl*dt*rkpv*ar*t3;
        Real rk = iw_(m,M1_IW_EGN,k,j,i) - ee - cl*dt*rkpv*ar*t4 + cl*dt*rkev*de0;
        Real emis = dt*ch*rkpv*ar;
        Real kk = (bk > 0.0) ? (emis*4.0*t3*cl*dt*rkev/bk) : 0.0;
        iw_(m,M1_IW_SRCB,k,j,i) = dt*ch*rkev - kk;
        iw_(m,M1_IW_SRCR,k,j,i) = emis*t4 - dt*ch*rkev*de0
                                  + ((bk > 0.0) ? (emis*4.0*t3*rk/bk) : 0.0);
      });
    }

    // (d) assemble the tridiagonal system of every column
    par_for("m1_impl_asm", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx = mbsize.d_view(m).dx1;
      Real nu = dt/dx;
      Real cr = ch/cl;
      Real wi = iw_(m,M1_IW_WCHI,k,j,i);
      Real ai = iw_(m,M1_IW_ADV,k,j,i);
      Real vi = iw_(m,M1_IW_V1,k,j,i);
      Real aa = 0.0, bb = 1.0, cc = 0.0;
      Real rr = iw_(m,M1_IW_EN,k,j,i) + iw_(m,M1_IW_SRCR,k,j,i);
      bb += iw_(m,M1_IW_SRCB,k,j,i);

      // ---- face i+1/2
      if (i < ie || cyclic) {
        int ip = (i < ie) ? (i+1) : is;
        Real ktf = 0.5*(opac_(m,M1_OP_T,k,j,i) + opac_(m,M1_OP_T,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = th*ch*ch*dt/dx;
        bb += nu*df*wi;
        cc -= nu*df*iw_(m,M1_IW_WCHI,k,j,ip);
        Real vf = 0.5*(vi + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,i) + iw_(m,M1_IW_G0,k,j,ip));
        rr -= nu*cr*th*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f);
        if (vf > 0.0) {
          bb += nu*cr*ai;
        } else {
          cc += nu*cr*iw_(m,M1_IW_ADV,k,j,ip);
        }
      } else if (bchi == M1_IBC_MARSHAK) {
        bb += nu*ch*mq;
        rr += nu*ch*mq*ebhi;
      } else if (bchi == M1_IBC_FLUX) {
        rr -= nu*cr*fxhi;
      }

      // ---- face i-1/2
      if (i > is || cyclic) {
        int im = (i > is) ? (i-1) : ie;
        Real ktf = 0.5*(opac_(m,M1_OP_T,k,j,im) + opac_(m,M1_OP_T,k,j,i));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = th*ch*ch*dt/dx;
        bb += nu*df*wi;
        aa -= nu*df*iw_(m,M1_IW_WCHI,k,j,im);
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + vi);
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,i));
        rr += nu*cr*th*(f0n_(m,k,j,i) - ch*dt*vf*g0f);
        if (vf > 0.0) {
          aa -= nu*cr*iw_(m,M1_IW_ADV,k,j,im);
        } else {
          bb -= nu*cr*ai;
        }
      } else if (bclo == M1_IBC_MARSHAK) {
        bb += nu*ch*mq;
        rr += nu*ch*mq*eblo;
      } else if (bclo == M1_IBC_FLUX) {
        rr += nu*cr*fxlo;
      }

      // a Dirichlet end cell: the whole row is replaced, which keeps the matrix an
      // M-matrix and anchors the level of E (see M1_IBC_EFIX)
      if (!cyclic && ((i == is && bclo == M1_IBC_EFIX) ||
                      (i == ie && bchi == M1_IBC_EFIX))) {
        aa = 0.0;
        bb = 1.0;
        cc = 0.0;
        rr = iw_(m,M1_IW_EN,k,j,i);
      }
      iw_(m,M1_IW_TA,k,j,i) = aa;
      iw_(m,M1_IW_TB,k,j,i) = bb;
      iw_(m,M1_IW_TC,k,j,i) = cc;
      iw_(m,M1_IW_TR,k,j,i) = rr;
    });

    // (e) one Thomas sweep per column (cyclic: Sherman-Morrison)
    par_for("m1_impl_thomas", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (!cyclic) {
        Real bet = iw_(m,M1_IW_TB,k,j,is);
        iw_(m,M1_IW_S2,k,j,is) = iw_(m,M1_IW_TR,k,j,is)/bet;
        for (int i=is+1; i<=ie; ++i) {
          iw_(m,M1_IW_S1,k,j,i) = iw_(m,M1_IW_TC,k,j,i-1)/bet;
          bet = iw_(m,M1_IW_TB,k,j,i) - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S1,k,j,i);
          iw_(m,M1_IW_S2,k,j,i) = (iw_(m,M1_IW_TR,k,j,i)
                                   - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S2,k,j,i-1))/bet;
        }
        for (int i=ie-1; i>=is; --i) {
          iw_(m,M1_IW_S2,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S2,k,j,i+1);
        }
      } else {
        // cyclic tridiagonal, Sherman-Morrison (Press et al. `cyclic`).  alpha is the
        // BOTTOM-LEFT corner, c(ie) (row ie coupling to cell is), and beta the TOP-RIGHT
        // one, a(is) (row is coupling to cell ie) -- that is the convention u and v below
        // are built for, and swapping them is invisible on a symmetric matrix but wrong
        // as soon as upwind advection makes the two off-diagonals differ: it produced a
        // spurious dipole across the seam (-14 % in the first cell, +11 % in the last)
        // on the very first step of T4b, where the exact answer is "nothing moves".
        Real alpha = iw_(m,M1_IW_TC,k,j,ie);
        Real beta = iw_(m,M1_IW_TA,k,j,is);
        Real gam = -iw_(m,M1_IW_TB,k,j,is);
        Real bb0 = iw_(m,M1_IW_TB,k,j,is) - gam;
        Real bbn = iw_(m,M1_IW_TB,k,j,ie) - alpha*beta/gam;
        // solve A' y = r and A' z = u with u = (gam,0,...,0,alpha)
        Real bet = bb0;
        iw_(m,M1_IW_S2,k,j,is) = iw_(m,M1_IW_TR,k,j,is)/bet;
        iw_(m,M1_IW_S3,k,j,is) = gam/bet;
        for (int i=is+1; i<=ie; ++i) {
          Real bd = (i == ie) ? bbn : iw_(m,M1_IW_TB,k,j,i);
          iw_(m,M1_IW_S1,k,j,i) = iw_(m,M1_IW_TC,k,j,i-1)/bet;
          bet = bd - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S1,k,j,i);
          iw_(m,M1_IW_S2,k,j,i) = (iw_(m,M1_IW_TR,k,j,i)
                                   - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S2,k,j,i-1))/bet;
          Real uu = (i == ie) ? alpha : 0.0;
          iw_(m,M1_IW_S3,k,j,i) = (uu
                                   - iw_(m,M1_IW_TA,k,j,i)*iw_(m,M1_IW_S3,k,j,i-1))/bet;
        }
        for (int i=ie-1; i>=is; --i) {
          iw_(m,M1_IW_S2,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S2,k,j,i+1);
          iw_(m,M1_IW_S3,k,j,i) -= iw_(m,M1_IW_S1,k,j,i+1)*iw_(m,M1_IW_S3,k,j,i+1);
        }
        // x = y - z (v.y)/(1 + v.z),  v = (1,0,...,0,beta/gam)
        Real vy = iw_(m,M1_IW_S2,k,j,is) + (beta/gam)*iw_(m,M1_IW_S2,k,j,ie);
        Real vz = iw_(m,M1_IW_S3,k,j,is) + (beta/gam)*iw_(m,M1_IW_S3,k,j,ie);
        Real fac = vy/(1.0 + vz);
        for (int i=is; i<=ie; ++i) {
          iw_(m,M1_IW_S2,k,j,i) -= fac*iw_(m,M1_IW_S3,k,j,i);
        }
      }
    });

    // (f) accept E', solve for T' and measure the Picard residual
    if (src_on) {
      auto eos = pmy_pack->phydro->peos->eos_data;
      par_for("m1_impl_tsolve", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efl);
        Real eold = iw_(m,M1_IW_EP,k,j,i);
        Real rkpv = opac_(m,M1_OP_P,k,j,i);
        Real rkev = opac_(m,M1_OP_E,k,j,i);
        Real told = iw_(m,M1_IW_TP,k,j,i);
        Real tnew = told;
        if (rkpv > 0.0 || rkev > 0.0) {
          Real dd = uh(m,IDN,k,j,i);
          Real de0 = iw_(m,M1_IW_DE0,k,j,i);
          bool ok = true;
          (void) M1ImplTemperature(eos, dd, told, iw_(m,M1_IW_EGN,k,j,i),
                                   cl*dt*rkpv*ar, cl*dt*rkev*(enew + de0), tnew, ok);
          if (!ok) {tnew = told;}
        }
        iw_(m,M1_IW_EP,k,j,i) = enew;
        iw_(m,M1_IW_TP,k,j,i) = tnew;
        Real re = fabs(enew - eold)/fmax(fabs(enew), 1.0e-300);
        Real rt = fabs(tnew - told)/fmax(fabs(tnew), 1.0e-300);
        iw_(m,M1_IW_RES,k,j,i) = fmax(re, rt);
      });
    } else {
      par_for("m1_impl_accept", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efl);
        Real eold = iw_(m,M1_IW_EP,k,j,i);
        iw_(m,M1_IW_EP,k,j,i) = enew;
        iw_(m,M1_IW_RES,k,j,i) = fabs(enew - eold)/fmax(fabs(enew), 1.0e-300);
      });
    }

    // (g) the face fluxes of the new iterate, and the derived cell flux
    par_for("m1_impl_face", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx = mbsize.d_view(m).dx1;
      bool lo = (i == is), hi = (i == ie+1);
      if ((lo || hi) && !cyclic) {
        int bc = lo ? bclo : bchi;
        int ic = lo ? is : ie;
        Real sgn = lo ? -1.0 : 1.0;
        Real fb = 0.0;
        if (bc == M1_IBC_MARSHAK) {
          fb = sgn*cl*mq*(iw_(m,M1_IW_EP,k,j,ic) - (lo ? eblo : ebhi));
        } else if (bc == M1_IBC_FLUX) {
          fb = lo ? fxlo : fxhi;
        } else if (bc == M1_IBC_EFIX) {
          // the Dirichlet row does not define a face flux; the adjacent interior face is
          // copied so that the cell flux and the momentum deposit stay finite.  The
          // boundary cell is held fixed anyway, so nothing downstream depends on it.
          fb = f0_(m,k,j,lo ? (is+1) : ie);
        }
        f0_(m,k,j,i) = fb;
      } else {
        int im = (i == is) ? ie : (i-1);
        int ip = (i == ie+1) ? is : i;
        Real ktf = 0.5*(opac_(m,M1_OP_T,k,j,im) + opac_(m,M1_OP_T,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,ip));
        Real gr = (iw_(m,M1_IW_WCHI,k,j,ip)*iw_(m,M1_IW_EP,k,j,ip)
                   - iw_(m,M1_IW_WCHI,k,j,im)*iw_(m,M1_IW_EP,k,j,im))/dx;
        Real fn = th*(f0n_(m,k,j,(i == ie+1 && cyclic) ? is : i)
                      - ch*cl*dt*gr - ch*dt*vf*g0f);
        f0_(m,k,j,i) = fn;
      }
    });
    if (cyclic) {
      // the wrapped face is stored twice; keep the two copies identical
      par_for("m1_impl_facewrap", DevExeSpace(), 0, nmb1, ks, ke, js, je,
      KOKKOS_LAMBDA(const int m, const int k, const int j) {
        f0_(m,k,j,ie+1) = f0_(m,k,j,is);
      });
    }
    par_for("m1_impl_f1", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      iw_(m,M1_IW_F1,k,j,i) = 0.5*(f0_(m,k,j,i) + f0_(m,k,j,i+1))
                              + iw_(m,M1_IW_ADV,k,j,i)*iw_(m,M1_IW_EP,k,j,i);
    });

    // (h) convergence
    resid = 0.0;
    Kokkos::parallel_reduce("m1_impl_res",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = iw_(m,M1_IW_RES,k,j,i);
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(resid));
    converged = (resid < impl_tol);
  }
#if MPI_PARALLEL_ENABLED
  {
    // a column never crosses a rank here, but the ITERATION COUNT must be global or the
    // ranks would take different numbers of Picard passes and diverge
    int ilocal = it, iglob;
    MPI_Allreduce(&ilocal, &iglob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    it = iglob;
  }
#endif
  impl_nstep += 1.0;
  impl_itsum += static_cast<Real>(it);
  impl_itmax = std::max(impl_itmax, static_cast<Real>(it));
  if (!converged) {impl_nfail += 1.0;}

  //------------------------------------------------------------------------- write back
  par_for("m1_impl_wb", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real ep = iw_(m,M1_IW_EP,k,j,i);
    Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
    Real fp1 = 0.5*(fl + fr) + iw_(m,M1_IW_ADV,k,j,i)*ep;

    Real work = 0.0, dm1 = 0.0, dmref = 0.0, eg = 0.0, ekin = 0.0, egrv = 0.0;
    Real dd = 0.0, v1 = 0.0;
    if (have_hydro) {
      dd = uh(m,IDN,k,j,i);
      Real idd = 1.0/fmax(dd, 1.0e-300);
      v1 = uh(m,IM1,k,j,i)*idd;
      ekin = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                  SQR(uh(m,IM3,k,j,i)))*idd;
      egrv = etg ? (dd*phicc(m,k,j,i)) : 0.0;
      eg = iw_(m,M1_IW_EGN,k,j,i);
      // (a) the energy the RADIATION gained from the gas over the step.  It is taken
      // from the ASSEMBLED row, q = SRCR - SRCB E', and not from rho kappa_P a T'^4:
      // the two differ by the Picard remainder of the linearisation, and only the first
      // is the amount the solved E actually received.  Setting the gas energy from it
      // makes e_gas + (c/chat) E change by the face fluxes and the work term ALONE, to
      // round-off -- measured: the T'^4 form drifted 2.8e-11 over 2000 steps of T5, this
      // one 0.  At convergence the two agree, so T' stays the consistent temperature.
      if (coupling && dbgh) {
        Real qq = iw_(m,M1_IW_SRCR,k,j,i) - iw_(m,M1_IW_SRCB,k,j,i)*ep;
        eg -= (cl/ch)*qq;
      }
      // (b) MOMENTUM.  Each x1 face hands dt (rho k_t)_f F0'_f/c to the gas, half to
      // each of its two cells (a physical boundary face gives all of it to its one
      // interior cell), which is exactly what the implicit face source removed from the
      // radiation: sum_cells dm = sum_faces dt (rho k_t)_f F0'_f/c.
      if (coupling && dbgf) {
        Real ktl, ktr;
        Real wl = 0.5, wr = 0.5;
        if (i == is && !cyclic) {
          ktl = opac_(m,M1_OP_T,k,j,i);
          wl = 1.0;
        } else {
          int im = (i == is) ? ie : (i-1);
          ktl = 0.5*(opac_(m,M1_OP_T,k,j,im) + opac_(m,M1_OP_T,k,j,i));
        }
        if (i == ie && !cyclic) {
          ktr = opac_(m,M1_OP_T,k,j,i);
          wr = 1.0;
        } else {
          int ip = (i == ie) ? is : (i+1);
          ktr = 0.5*(opac_(m,M1_OP_T,k,j,i) + opac_(m,M1_OP_T,k,j,ip));
        }
        dm1 = (dt/cl)*(wl*ktl*fl + wr*ktr*fr);
        dmref = fref ? (dt*dd*aref_(m,k,j,i)) : 0.0;
        if (feedback) {
          Real w1 = (uh(m,IM1,k,j,i) + dm1)/fmax(dd, 1.0e-300);
          work = 0.5*(v1 + w1)*dm1;
          ep -= (ch/cl)*work;
        }
      }
    }

    Real f2 = 0.0, f3 = 0.0;
    M1ApplyLimits(cl, efl, ep, fp1, f2, f3);
    u0_(m,M1_E,k,j,i) = ep;
    u0_(m,M1_F1,k,j,i) = fp1;
    u0_(m,M1_F2,k,j,i) = 0.0;
    u0_(m,M1_F3,k,j,i) = 0.0;
    if (have_hydro && feedback) {
      uh(m,IM1,k,j,i) = uh(m,IM1,k,j,i) + dm1 - dmref;
      uh(m,IEN,k,j,i) = eg + ekin + egrv + work;
    }
  });

  return TaskStatus::complete;
}

} // namespace radm1
