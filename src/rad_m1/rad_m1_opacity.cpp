//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_opacity.cpp
//! \brief fills RadiationM1::opac with rho*kappa_{P,E,transport} once per stage, before
//! the fluxes and including the ghost cells.  The state is read from the conserved u0
//! (the coupling updates u0 inside the stage, so w0 lags).

#include <math.h>

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_parfor.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "utils/eint_from_cons.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {
//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationM1::Opacity

TaskStatus RadiationM1::Opacity(Driver *pdrive, int stage) {
  // first gas reader of both M1 chains: refresh |B|^2/2 of the stage's field
  if (fl_mhd && !emag_hold) EmagBuild(false);
  if (opac_zero) return TaskStatus::complete;
  if (!fl_on) return TaskStatus::complete;
  // <rad_m1>/opac_freeze: fill the array once and keep it
  if (opac_freeze && opac_frozen) return TaskStatus::complete;
  if (opacity_type == M1_OPAC_TABLE && otab.nT <= 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/opacity = table, but no tables were handed to the "
      << "module: the problem generator must call RadiationM1::SetOpacityTables"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int n1 = indcs.nx1 + 2*(indcs.ng);
  int n2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int n3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb1 = pmy_pack->nmb_thispack - 1;

  auto opac_ = opac;
  FluidRef fl = FluidRef::Get(pmy_pack);
  auto uh = fl.u0;
  auto eos = fl.eos;
  // with <hydro>/etotgrav the conserved energy includes rho*Phi: subtract it as well
  const bool etg = fl.use_etotgrav;
  auto phicc = fl.phicc0;
  const bool mhd = fl_mhd;      // m1-mhd: and |B|^2/2
  auto emag_ = emag0;
  // cubed sphere: the momentum is covariant, KE = 0.5 m_a v^a (gnomonic_raisevel.hpp)
  const bool csk = pmy_pack->pmesh->use_cubed_sphere;
  auto ccl = pmy_pack->pcoord->cos_cell;
  int otype = opacity_type;
  Real kp = kappa_p, ke = kappa_e, kf = kappa_f, ks = kappa_s;
  Real rref = opac_rho_ref, tref = opac_t_ref, aa = opac_a, bb = opac_b;
  bool need_t = (otype != M1_OPAC_CONST);
  // copy the tables (held by value) into the closure instead of capturing this
  M1OpacTab ot = otab;
  const bool osc = vscat;   // vet_scatter: rho kappa_e -> opac(M1_OP_S)
  if (ot.scat == 2 && !eos.tbl.active) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1>/vet_scatter_kappa_e = eos needs a general (tabulated) "
      << "EOS: the electron fraction comes from its Saha composition" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  par_for_lb("m1_opacity", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real d = uh(m,IDN,k,j,i);
    Real t = 0.0;
    if (need_t) {
      Real ke_dens = 0.5*(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                          SQR(uh(m,IM3,k,j,i)))/fmax(d, 1.0e-300);
      if (csk) {
        const Real c = ccl(m,k,j);
        const Real m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
        ke_dens = CsKeMom(uh(m,IM1,k,j,i), m2, m3, c)
                  /fmax(d, 1.0e-300);
      }
      Real eint = uh(m,IEN,k,j,i) - ke_dens;
      if (etg) eint -= d*phicc(m,k,j,i);
      if (mhd) eint -= emag_(m,k,j,i);
      t = eos.Temperature(d, fmax(eint, 0.0));
    }
    Real op, oe, of, os;
    if (otype == M1_OPAC_TABLE) {
      M1TableOpacities(ot, d, t, op, oe, of, os);
      M1ScatterEos(ot, eos, d, t, of, os);
    } else {
      M1Opacities(otype, d, t, kp, ke, kf, ks, rref, tref, aa, bb, op, oe, of, os);
    }
    opac_(m,M1_OP_P,k,j,i) = d*op;
    opac_(m,M1_OP_E,k,j,i) = d*oe;
    opac_(m,M1_OP_T,k,j,i) = d*(of + os);
    if (osc) {opac_(m,M1_OP_S,k,j,i) = d*os;}
  });

  // <rad_m1>/dbg_opac_patch (diagnostic): opacities times f inside an x1-x2 box
  if (dbg_opac_patch != 1.0) {
    const Real pf = dbg_opac_patch;
    const Real x1a = dbg_opac_x1lo, x1b = dbg_opac_x1hi;
    const Real x2a = dbg_opac_x2lo, x2b = dbg_opac_x2hi;
    const int is = indcs.is, js = indcs.js;
    const int nx1 = indcs.nx1, nx2 = indcs.nx2;
    auto &size = pmy_pack->pmb->mb_size;
    par_for("m1_opac_patch", DevExeSpace(), 0, nmb1, 0, n3-1, 0, n2-1, 0, n1-1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real x1v = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
      Real x2v = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
      if (x1v >= x1a && x1v <= x1b && x2v >= x2a && x2v <= x2b) {
        opac_(m,M1_OP_P,k,j,i) *= pf;
        opac_(m,M1_OP_E,k,j,i) *= pf;
        opac_(m,M1_OP_T,k,j,i) *= pf;
        if (osc) {opac_(m,M1_OP_S,k,j,i) *= pf;}
      }
    });
  }

  opac_frozen = true;
  if (!rsla_done) RSLACheck();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::RSLACheck
//! \brief reduced-speed-of-light admissibility check, run once on the first filled
//! opacity array: prints K = chat/(v_max tau_max), with tau_max the largest x1 column
//! optical depth over the active cells and v_max = <rad_m1>/rsla_vmax if positive, else
//! the measured max |v|

void RadiationM1::RSLACheck() {
  rsla_done = true;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  auto opac_ = opac;
  auto uh = FluidRef::Get(pmy_pack).u0;
  auto &size = pmy_pack->pmb->mb_size;

  Real taumax = 0.0;
  Kokkos::parallel_reduce("m1_rsla_tau",
  Kokkos::MDRangePolicy<Kokkos::Rank<3>>(DevExeSpace(), {0,ks,js}, {nmb1+1,ke+1,je+1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, Real &lmax) {
    Real tau = 0.0;
    for (int i=is; i<=ie; ++i) {tau += opac_(m,M1_OP_T,k,j,i)*size.d_view(m).dx1;}
    lmax = (tau > lmax) ? tau : lmax;
  }, Kokkos::Max<Real>(taumax));

  Real vmax = 0.0;
  Kokkos::parallel_reduce("m1_rsla_v",
  Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                         {nmb1+1,ke+1,je+1,ie+1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
    Real d = fmax(uh(m,IDN,k,j,i), 1.0e-300);
    Real v = sqrt(SQR(uh(m,IM1,k,j,i)) + SQR(uh(m,IM2,k,j,i)) +
                  SQR(uh(m,IM3,k,j,i)))/d;
    lmax = (v > lmax) ? v : lmax;
  }, Kokkos::Max<Real>(vmax));

#if MPI_PARALLEL_ENABLED
  {
    Real loc[2] = {taumax, vmax}, glob[2];
    MPI_Allreduce(loc, glob, 2, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    taumax = glob[0];
    vmax = glob[1];
  }
#endif
  Real vmeas = vmax;
  if (rsla_vmax > 0.0) vmax = rsla_vmax;
  // chat = c: nothing to check
  if (chat_over_c >= 1.0) {
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_m1> RSLA check: chat = c, no reduced speed of light (tau_max="
                << taumax << ", measured v_max=" << vmeas << ")" << std::endl;
    }
    return;
  }
  Real denom = vmax*((taumax > 1.0) ? taumax : 1.0);
  Real ratio = (chat > 0.0) ? (denom/chat) : 0.0;
  Real kk = (denom > 0.0) ? (chat/denom) : 0.0;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1> RSLA check: tau_max=" << taumax
              << " v_max=" << vmax << " (measured " << vmeas << ")"
              << " chat=" << chat << " -> K = chat/(v_max tau_max) = " << kk
              << ", v_max tau_max/chat = " << ratio << std::endl;
  }
  if (ratio > 1.0 && !rsla_force) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_m1> reduced speed of light is inadmissible: "
      << "v_max*tau_max/chat = " << ratio << " > 1 (design sect. 2).  Raise "
      << "chat_over_c, or set <rad_m1>/rsla_force = true to run anyway" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ratio > rsla_warn && global_variable::my_rank == 0) {
    std::cout << "### WARNING: <rad_m1> v_max*tau_max/chat = " << ratio
              << " exceeds <rad_m1>/rsla_warn = " << rsla_warn
              << "; the reduced speed of light is marginal (design sect. 2)"
              << std::endl;
  }
  return;
}

} // namespace radm1
