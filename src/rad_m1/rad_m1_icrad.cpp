//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_icrad.cpp
//! \brief icrad-1011: an opt-in start of the radiation field that is consistent with the
//! VET transport (implicit transport, closure = vet_col with vet_gd).  FRESH STARTS
//! ONLY: Driver::Initialize calls IcRadRun when res_flag is false, after the first dt
//! and before the cycle-0 outputs; a restart never runs it (the keys are then inert).
//!
//! <rad_m1>/ic_rad (default fs_relax for implicit vet_col + vet_gd with a fluid, none for
//! every other configuration; none = bitwise the old start):
//!   fs        (a) formal-solution start: ic_rad_fs_iter passes of E := J, J the vet_gd
//!             formal solution of the current state (twin-corrected under the fused
//!             twin: J *= <J_twin>/J_twin, the same correction as vgd_hr's J), the IC
//!             flux F and the gas kept; the tensor is rebuilt after each pass, so the
//!             last pass leaves E = J(E_{n-1}) and the tensor of that E.  Gas and hydro
//!             frozen; no implicit solve (the pass returns after the build).
//!   relax     (b) radiation-only relaxation: up to ic_rad_relax_n implicit radiation
//!             steps of dt_relax = ic_rad_relax_dtfac x dt with the hydro frozen (no
//!             hydro stages, the time and the cycle do not advance) and no radiative
//!             force on the gas (dbg_gas_force(_trans) = false for the IC steps only:
//!             the velocity stays the IC's), the gas-radiation ENERGY exchange on (the
//!             envelope gas relaxes to the radiation temperature with it).  Stops when
//!             |L_top/L_top(prev) - 1| < ic_rad_relax_tol (default 1e-4; at most
//!             ic_rad_relax_n = 40 steps; deterministic sums).
//!   fs_relax  (a) then (b).
//! ic_rad_diag = N (default 0): the shell means of E, J (twin-corrected) and F1 at every
//!   IC pass and after the vet_col build of the first N solves of the run, appended to
//!   icrad_ej.txt (rank 0, run directory).  Read-only: also usable with ic_rad = none.
//! After the IC passes the predictor history is dropped (pred_ok = pred2_ok = false), so
//! the first step starts cold as on any fresh start.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/m1_fluid.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radm1 {
namespace {
void IcrFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in rad_m1_icrad.cpp" << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

// rank-ordered sum of per-rank host vectors (reproducible run to run, as vet_gd_twin_det)
void IcrDetSum(std::vector<Real> &v) {
#if MPI_PARALLEL_ENABLED
  const int nr = global_variable::nranks;
  if (nr > 1) {
    const int n = static_cast<int>(v.size());
    std::vector<Real> all(static_cast<size_t>(n)*nr);
    MPI_Allgather(v.data(), n, MPI_ATHENA_REAL, all.data(), n, MPI_ATHENA_REAL,
                  MPI_COMM_WORLD);
    for (int i = 0; i < n; ++i) {
      Real a = 0.0;
      for (int r = 0; r < nr; ++r) {a += all[static_cast<size_t>(r)*n + i];}
      v[i] = a;
    }
  }
#endif
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::VsrcNoes
//! \brief vet_source_noesrc in force for the current solve: every backward-Euler solve
//! (time_scheme = be; under hesdirk2 the BESTORE / fallback / ic_rad solves unless
//! vet_source_noesrc_h2be = false); never a hesdirk2 STAGE solve (its build reads iw EN,
//! where time2_vet_col = predict puts the predicted state)

bool RadiationM1::VsrcNoes() const {
  if (!vsrc_noes) return false;
  if (time_scheme == M1_TIME_BE) return true;
  return vsrc_noes_h2be && (t2_solve != M1_T2S_STAGE1) && (t2_solve != M1_T2S_STAGE2);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::IcRadInit
//! \brief the ic_rad keys (read only when named) and their restrictions

void RadiationM1::IcRadInit(ParameterInput *pin) {
  // icrad-1011 default (user 10-11): fs_relax on a FRESH start of an implicit vet_col +
  // vet_gd run with a fluid (the configuration it was tested on); none for every other
  // problem (no M1 / explicit / no vet_gd: the hook is not even armed).  A restart never
  // runs it (Driver::Initialize, res_flag), whatever its input says.  Opt out: none.
  // prodcand-1011: and only in a he_star_m1 build (the stellar-envelope problems: AG Car,
  // He giant, BSG; CMake defines M1_ICRAD_DEFAULT_FS_RELAX for this file).  Every other
  // pgen (the built-in tests: a time-dependent pulse must keep its IC radiation) defaults
  // to none.  An explicit ic_rad key wins everywhere.
#if defined(M1_ICRAD_DEFAULT_FS_RELAX)
  const bool elig = (transport >= M1_TRANSPORT_IMPLICIT_X1) && vet_col && vgd_on && fl_on;
#else
  const bool elig = false;
#endif
  icr_mode = elig ? 3 : 0;
  if (pin->DoesParameterExist("rad_m1","ic_rad")) {
    const std::string s = pin->GetString("rad_m1","ic_rad");
    if (s == "none") {
      icr_mode = 0;
    } else if (s == "fs") {
      icr_mode = 1;
    } else if (s == "relax") {
      icr_mode = 2;
    } else if (s == "fs_relax") {
      icr_mode = 3;
    } else {
      IcrFatal("<rad_m1>/ic_rad must be none, fs, relax or fs_relax");
    }
  }
  if (pin->DoesParameterExist("rad_m1","ic_rad_diag")) {
    icr_diag = pin->GetInteger("rad_m1","ic_rad_diag");
  }
  if (icr_mode == 0 && icr_diag <= 0) return;
  icr_fs_iter = pin->GetOrAddInteger("rad_m1","ic_rad_fs_iter",3);
  icr_relax_n = pin->GetOrAddInteger("rad_m1","ic_rad_relax_n",40);
  icr_relax_tol = pin->GetOrAddReal("rad_m1","ic_rad_relax_tol",1.0e-4);
  icr_relax_dtfac = pin->GetOrAddReal("rad_m1","ic_rad_relax_dtfac",1.0);
  if (transport < M1_TRANSPORT_IMPLICIT_X1 || !vet_col || !vgd_on) {
    IcrFatal("<rad_m1>/ic_rad and ic_rad_diag need transport = implicit, closure = "
             "vet_col and vet_gd = true");
  }
  if ((icr_mode & 1) && icr_fs_iter < 1) {IcrFatal("ic_rad_fs_iter must be >= 1");}
  if ((icr_mode & 2) && (icr_relax_n < 1 || !(icr_relax_dtfac > 0.0))) {
    IcrFatal("ic_rad_relax_n must be >= 1 and ic_rad_relax_dtfac > 0");
  }
  if (global_variable::my_rank == 0) {
    const char *nm[4] = {"none", "fs", "relax", "fs_relax"};
    std::cout << "<rad_m1> ic_rad = " << nm[icr_mode] << " (fresh starts only)";
    if (icr_mode & 1) {std::cout << " fs_iter = " << icr_fs_iter;}
    if (icr_mode & 2) {
      std::cout << " relax_n = " << icr_relax_n << " relax_tol = " << icr_relax_tol
                << " relax_dtfac = " << icr_relax_dtfac;
    }
    std::cout << " diag = " << icr_diag
              << (global_variable::restart_run ? " [RESTART: the IC passes are skipped]"
                                               : "")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::IcRadShellSum
//! \brief global sums over the active lateral cells (m, k, j) of components 0..nc-1 of
//! a(m, c, k, j, i), per shell i: one thread per (c, i) in a fixed (m, k, j) order, then
//! the rank-ordered sum (reproducible).  out[c*c1 + i] for component c0 + c

void RadiationM1::IcRadShellSum(DvceArray5D<Real> &a, const int c0, const int nc,
                                std::vector<Real> &out) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int nm = pmy_pack->nmb_thispack;
  DvceArray2D<Real> s("icr_ss", nc, c1);
  auto a_ = a;
  par_for("icr_shsum", DevExeSpace(), 0, nc - 1, is, ie,
  KOKKOS_LAMBDA(const int c, const int i) {
    Real acc = 0.0;
    for (int m = 0; m < nm; ++m) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {acc += a_(m,c0+c,k,j,i);}
      }
    }
    s(c,i) = acc;
  });
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), s);
  out.assign(static_cast<size_t>(nc)*c1, 0.0);
  for (int c = 0; c < nc; ++c) {
    for (int i = 0; i < c1; ++i) {out[static_cast<size_t>(c)*c1 + i] = h(c,i);}
  }
  IcrDetSum(out);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::IcRadJ
//! \brief J of the last vet_gd build per cell (E units, as m1_fs), shells above the cut;
//! under the fused twin divided by the twin's J/<J_twin> (vgd_hr's J correction).
//! jc(m, 0, k, j, i) = J (0 below the cut); component 1 is scratch.

void RadiationM1::IcRadJ(DvceArray5D<Real> &jc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int ilo = is + vlat_icut;
  const int og = vgd_w - indcs.ng;
  const int nd = vgd_n;
  const bool tw = vgd_twin && vgd_twfuse && !vgd_twseq &&
                  (vgd_itw.d.extent_int(0) > 1);
  auto vi_ = vgd_i;
  auto vt_ = vgd_itw;
  auto dir_ = vgd_dir;
  Kokkos::deep_copy(jc, 0.0);
  par_for("icr_j", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real jr = 0.0, jt = 0.0;
    for (int d = 0; d < nd; ++d) {
      const Real w = dir_(d,3);
      jr += w*vi_(m,d,k+og,j+og,i);
      if (tw) {jt += w*vt_(m,d,k+og,j+og,i);}
    }
    jc(m,0,k,j,i) = jr;
    jc(m,1,k,j,i) = jt;
  });
  if (!tw) {
    if (!icr_warned && global_variable::my_rank == 0) {
      std::cout << "<rad_m1> ic_rad: no fused twin intensities: J is the raw formal "
                << "solution (ray noise not removed)" << std::endl;
    }
    icr_warned = true;
    return;
  }
  // <J_twin> per shell (all ranks, reproducible), then J *= <J_twin>/J_twin
  std::vector<Real> sm;
  IcRadShellSum(jc, 1, 1, sm);
  const Real ncol = static_cast<Real>(pmy_pack->pmesh->mesh_indcs.nx2)*
                    pmy_pack->pmesh->mesh_indcs.nx3;
  DvceArray1D<Real> mt("icr_jtm", c1);
  auto mh = Kokkos::create_mirror_view(mt);
  for (int i = 0; i < c1; ++i) {mh(i) = sm[i]/ncol;}
  Kokkos::deep_copy(mt, mh);
  par_for("icr_jtw", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real jt = jc(m,1,k,j,i);
    if (jt > 0.0 && mt(i) > 0.0) {jc(m,0,k,j,i) *= mt(i)/jt;}
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::IcRadProfile
//! \brief the shell means of E, J (jc component 0) and F1 -> icrad_ej.txt (when
//! ic_rad_diag > 0), and one summary line on stdout: the max over the gd shells above
//! the first of |<E>/<J> - 1| and <E>/<J> at the top shell

void RadiationM1::IcRadProfile(const std::string &tag, DvceArray5D<Real> &jc) {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  const int n3 = jc.extent_int(2), n2 = jc.extent_int(3);
  DvceArray5D<Real> a("icr_prof", nmb1 + 1, 3, n3, n2, c1);
  auto u0_ = u0;
  par_for("icr_pf", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    a(m,0,k,j,i) = u0_(m,M1_E,k,j,i);
    a(m,1,k,j,i) = jc(m,0,k,j,i);
    a(m,2,k,j,i) = u0_(m,M1_F1,k,j,i);
  });
  std::vector<Real> s;
  IcRadShellSum(a, 0, 3, s);
  if (global_variable::my_rank != 0) return;
  const int ilo = is + vlat_icut;
  Real emax = 0.0, rtop = 0.0;
  for (int i = ilo + 1; i <= ie; ++i) {
    const Real ej = (s[c1 + i] > 0.0) ? s[i]/s[c1 + i] : 0.0;
    emax = std::max(emax, std::fabs(ej - 1.0));
    if (i == ie) {rtop = ej;}
  }
  std::printf("<rad_m1> ic_rad %-12s cycle=%d max|<E>/<J>-1| (gd shells) = %.4e "
              "<E>/<J>(top) = %.5f\n", tag.c_str(), pmy_pack->pmesh->ncycle, emax, rtop);
  if (icr_diag <= 0) return;
  auto x1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmy_pack->pcoord->x1v);
  const Real ncol = static_cast<Real>(pmy_pack->pmesh->mesh_indcs.nx2)*
                    pmy_pack->pmesh->mesh_indcs.nx3;
  FILE *fp = std::fopen("icrad_ej.txt", "a");
  if (fp == nullptr) return;
  for (int i = is; i <= ie; ++i) {
    std::fprintf(fp, "%s %d %d %.10e %.10e %.10e %.10e\n", tag.c_str(),
                 pmy_pack->pmesh->ncycle, i - is, x1v(0,i), s[i]/ncol, s[c1 + i]/ncol,
                 s[2*c1 + i]/ncol);
  }
  std::fclose(fp);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationM1::IcRadLtop
//! \brief the luminosity through the outer x1 face, sum f0x1(ie+1) area (as the
//! he_star_m1 history's L_top), reproducible (fixed order, rank-ordered)

Real RadiationM1::IcRadLtop() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nm = pmy_pack->nmb_thispack;
  auto f0_ = f0x1;
  auto ar_ = pmy_pack->pcoord->area.x1f;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  DvceArray1D<Real> s("icr_lt", 1);
  par_for("icr_ltop", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int) {
    Real acc = 0.0;
    for (int m = 0; m < nm; ++m) {
      if (mbbcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::block) continue;
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {acc += f0_(m,k,j,ie+1)*ar_(m,k,j,ie+1);}
      }
    }
    s(0) = acc;
  });
  auto h = Kokkos::create_mirror_view_and_copy(HostMemSpace(), s);
  std::vector<Real> v(1, h(0));
  IcrDetSum(v);
  return v[0];
}

//----------------------------------------------------------------------------------------
//! \fn bool RadiationM1::IcRadHook
//! \brief called by ImplicitSolve right after the start-of-step vet_col (vet_gd) build:
//! (a) the E := J passes when pending; the E/J diagnostic otherwise.  Returns true when
//! the solve is to be skipped (the formal-solution pass of IcRadRun).

bool RadiationM1::IcRadHook() {
  const bool fs = icr_fs_pending;
  const bool diag = (icr_diag > 0) && (icr_active || icr_ndiag_done < icr_diag);
  if (!fs && !diag) return false;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int c1 = indcs.nx1 + 2*indcs.ng;
  const int n3 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*indcs.ng) : 1;
  const int n2 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*indcs.ng) : 1;
  const int nmb1 = pmy_pack->nmb_thispack - 1;
  DvceArray5D<Real> jc("icr_jc", nmb1 + 1, 2, n3, n2, c1);
  if (!fs) {
    IcRadJ(jc);
    IcRadProfile(icr_active ? icr_tag : (std::string("run_c") +
                 std::to_string(pmy_pack->pmesh->ncycle)), jc);
    if (!icr_active) {icr_ndiag_done += 1;}
    return false;
  }
  // ---- (a) E := J, ic_rad_fs_iter passes, the tensor rebuilt after each
  Kokkos::fence();
  Kokkos::Timer tm;
  icr_fs_pending = false;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const Real efl = e_floor;
  const bool eso = esrc_on;
  auto es_ = esrc;
  const Real dtes = dt_sub*chat/c_light;
  auto u0_ = u0;
  auto iw_ = iw;
  for (int it = 0; it < icr_fs_iter; ++it) {
    IcRadJ(jc);
    IcRadProfile("fs" + std::to_string(it), jc);
    // the first gd shell keeps its E (its inflow is the approximate start condition)
    const int ilo = is + vlat_icut + 1;
    par_for("icr_e2j", DevExeSpace(), 0, nmb1, ks, ke, js, je, ilo, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real jj = jc(m,0,k,j,i);
      if (!(jj > 0.0)) return;
      const Real e = fmax(jj, efl);
      u0_(m,M1_E,k,j,i) = e;
      iw_(m,M1_IW_EP,k,j,i) = e;
      iw_(m,M1_IW_EN,k,j,i) = e + (eso ? dtes*es_(m,k,j,i) : 0.0);
    });
    if (trans_on) {
      ImplicitTransverseHalo(M1_NHALO_T);
    } else {
      ImplicitX1Halo(true);
    }
    VetColBuild();   // the tensor (and J) of the new E
  }
  IcRadJ(jc);
  IcRadProfile("fs_end", jc);
  Kokkos::fence();
  icr_tfs += tm.seconds();
  return icr_fs_only;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::IcRadRun
//! \brief the IC passes (fresh start; called by Driver::Initialize after the first dt).
//! Each pass is one ordinary m1 task-list execution (closure limits, opacity, the
//! implicit solve with its vet_col build, the fluid refresh, the ghost exchange) with
//! t2_solve = none; the driver re-initializes the boundaries and the dt afterwards.

void RadiationM1::IcRadRun(Driver *pdrive) {
  if (icr_mode == 0 || global_variable::restart_run) return;
  Mesh *pm = pmy_pack->pmesh;
  Kokkos::fence();
  Kokkos::Timer tall;
  const bool sgf = dbg_gas_force, sgft = dbg_gas_force_trans, svimp = impl_vimp;
  const int st2 = t2_solve;
  const Real dt0 = pm->dt;
  dbg_gas_force = false;        // the gas momentum is frozen during the IC passes
  dbg_gas_force_trans = false;
  impl_vimp = false;            // (vimp IS the radiative force on the gas)
  t2_solve = M1_T2S_NONE;       // plain backward Euler, no hesdirk2 slope stored
  icr_active = true;
  auto pass = [&]() {
    const int nsub = SetSubsteps(pm->dt);
    for (int s = 0; s < nsub; ++s) {
      for (int stage = 1; stage <= nstage; ++stage) {
        pdrive->ExecuteTaskList(pm, "m1_before_stagen", stage);
        pdrive->ExecuteTaskList(pm, "m1_stagen", stage);
        pdrive->ExecuteTaskList(pm, "m1_after_stagen", stage);
      }
      (void) ApplyClosureLimits(pdrive, 1);
    }
  };
  const Real l0 = IcRadLtop();
  if (global_variable::my_rank == 0) {
    std::printf("<rad_m1> ic_rad: start L_top = %.10e dt = %.6e\n", l0, dt0);
  }
  if (icr_mode & 1) {
    icr_tag = "fs";
    icr_fs_pending = true;
    icr_fs_only = true;
    pass();
    icr_fs_only = false;
    icr_fs_pending = false;
    if (global_variable::my_rank == 0) {
      std::printf("<rad_m1> ic_rad fs: %d passes, %.3f s (L_top keeps the IC flux)\n",
                  icr_fs_iter, icr_tfs);
    }
  }
  if (icr_mode & 2) {
    pm->dt = icr_relax_dtfac*dt0;
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int nmb1 = pmy_pack->nmb_thispack - 1;
    DvceArray4D<Real> eold("icr_eold", nmb1 + 1, u0.extent_int(2), u0.extent_int(3),
                           u0.extent_int(4));
    auto u0_ = u0;
    Real lp = IcRadLtop();
    for (int n = 0; n < icr_relax_n; ++n) {
      Kokkos::fence();
      Kokkos::Timer t1;
      par_for("icr_eo", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        eold(m,k,j,i) = u0_(m,M1_E,k,j,i);
      });
      icr_tag = "relax" + std::to_string(n);
      pass();
      Real dem = 0.0;
      Kokkos::parallel_reduce("icr_dem",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                             {nmb1+1,ke+1,je+1,ie+1}),
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lm) {
        const Real r = fabs(u0_(m,M1_E,k,j,i)/fmax(eold(m,k,j,i), 1.0e-300) - 1.0);
        lm = (r > lm) ? r : lm;
      }, Kokkos::Max<Real>(dem));
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, &dem, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
      const Real l = IcRadLtop();
      const Real dl = std::fabs(l/lp - 1.0);
      Kokkos::fence();
      icr_trx += t1.seconds();
      icr_nrx += 1;
      if (global_variable::my_rank == 0) {
        std::printf("<rad_m1> ic_rad relax %2d: L_top = %.10e (L/L_IC = %.6f) "
                    "|dL/L| = %.3e max|dE/E| = %.3e %.3f s\n", n, l, l/l0, dl, dem,
                    t1.seconds());
      }
      lp = l;
      if (dl < icr_relax_tol && n >= 1) break;
    }
    pm->dt = dt0;
  }
  // dic-1013 (icr_fmeas, set by the pgen; default false = nothing): one pass at the run's
  // dt with the radiative force on the gas (as the first step will apply it); the net
  // x1 momentum it deposits, dm - dt rho a_ref, is kept per cell in icr_dm, then the gas
  // momentum is restored and its kinetic-energy change taken back out of IEN (the gas
  // keeps the internal energy of the pass)
  FluidRef flr = FluidRef::Get(pmy_pack);
  if (icr_fmeas && flr.on) {
    auto &indcs = pmy_pack->pmesh->mb_indcs;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int nmb1 = pmy_pack->nmb_thispack - 1;
    auto uh = flr.u0;
    DvceArray5D<Real> sv("icr_fm_sv", nmb1 + 1, 3, uh.extent_int(2), uh.extent_int(3),
                         uh.extent_int(4));
    par_for("icr_fm_save", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      sv(m,0,k,j,i) = uh(m,IM1,k,j,i);
      sv(m,1,k,j,i) = uh(m,IM2,k,j,i);
      sv(m,2,k,j,i) = uh(m,IM3,k,j,i);
    });
    Kokkos::realloc(icr_dm, nmb1 + 1, uh.extent_int(2), uh.extent_int(3),
                    uh.extent_int(4));
    Kokkos::deep_copy(icr_dm, 0.0);
    pm->dt = dt0;
    dbg_gas_force = sgf;
    dbg_gas_force_trans = sgft;
    icr_tag = "fmeas";
    pass();
    dbg_gas_force = false;
    dbg_gas_force_trans = false;
    auto dm_ = icr_dm;
    par_for("icr_fm_rest", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real d = uh(m,IDN,k,j,i);
      const Real m1 = uh(m,IM1,k,j,i), m2 = uh(m,IM2,k,j,i), m3 = uh(m,IM3,k,j,i);
      const Real o1 = sv(m,0,k,j,i), o2 = sv(m,1,k,j,i), o3 = sv(m,2,k,j,i);
      dm_(m,k,j,i) = m1 - o1;
      uh(m,IEN,k,j,i) -= 0.5*((m1*m1 + m2*m2 + m3*m3) - (o1*o1 + o2*o2 + o3*o3))/d;
      uh(m,IM1,k,j,i) = o1;
      uh(m,IM2,k,j,i) = o2;
      uh(m,IM3,k,j,i) = o3;
    });
    icr_dm_dt = dt0;
    if (global_variable::my_rank == 0) {
      std::printf("<rad_m1> ic_rad: force pass (dic-1013) at dt = %.6e, L_top = "
                  "%.10e\n", dt0, IcRadLtop());
    } else {
      (void) IcRadLtop();
    }
  }
  icr_active = false;
  dbg_gas_force = sgf;
  dbg_gas_force_trans = sgft;
  impl_vimp = svimp;
  t2_solve = st2;
  pred_ok = false;    // the first step starts cold, as on any fresh start
  pred2_ok = false;
  (void) SetSubsteps(dt0);
  const Real lend = IcRadLtop();   // (collective)
  Kokkos::fence();
  if (global_variable::my_rank == 0) {
    std::printf("<rad_m1> ic_rad: done, %.3f s (fs %.3f s, relax %d steps %.3f s), "
                "L_top = %.10e (L/L_IC = %.6f)\n", tall.seconds(), icr_tfs, icr_nrx,
                icr_trx, lend, lend/l0);
  }
}

} // namespace radm1
