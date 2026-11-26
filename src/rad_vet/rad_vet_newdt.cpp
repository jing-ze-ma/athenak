//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_vet_newdt.cpp
//! \brief time step, the reduced-light-speed check, the radiation inputs of the hydro's
//! signal speed.
//!
//!  called from: rad_vet_opacity.cpp, driver/driver.cpp, hydro/hydro_newdt.cpp,
//!    mhd/mhd_newdt.cpp, outputs/restart.cpp, pgen/pgen.cpp
//!  calls into: rad_vet_opacity.cpp

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "tasklist/task_list.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "coordinates/coordinates.hpp"
#include "driver/driver.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "rad_vet/rad_vet_args.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

int RadiationVet::RssRstMode() const {return RssMode() + 16;}

//----------------------------------------------------------------------------------------
//! \fn int RadiationVet::RssMode
//! \brief the closure tensor <hydro>/rad_signal_speed reads (hydro_newdt.cpp)

int RadiationVet::RssMode() const {
  // clean: its own Eddington tensor cln_s, once CleanInit ran
  if (cln_s.extent(0) > 0 && cln_init) return 5;
  return (vet_sc && vet_full && vet_cell.extent(0) > 0) ? 4 : 0;
}

int RadiationVet::RssRstNch(int mode) const {
  // extended restart state (mode + 16): also the other opacity components.  Modes 2, 3
  // (files of removed closure tensors): their channels are counted and read past
  if (mode >= 16) return RssRstNch(mode - 16) + (opac.extent_int(1) - 1);
  if (mode == 2 || mode == 3) return 5;
  if (mode == 4) return 4;
  if (mode == 5) return 7;   // opac_T and the six components of cln_s's D
  return 1;
}

void RadiationVet::RssSlot(int mode, int n, DvceArray5D<Real> *&a, int &c) {
  if (mode >= 16) {   // base slots, then opac (all but T)
    const int b = mode - 16, nb = RssRstNch(b);
    if (n < nb) {RssSlot(b, n, a, c); return;}
    const int k = n - nb;
    a = &opac;
    c = (k < VET_OP_T) ? k : (k + 1);
    return;
  }
  if (n == 0) {a = &opac; c = VET_OP_T; return;}
  if (mode == 5) {a = &cln_s; c = CLN_D11 + n - 1; return;}
  a = &vet_cell;
  c = VETRAY_D11 + n - 1;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RssRstPack
//! \brief signal-speed inputs of MeshBlocks m0 .. m0+nmb-1, channel by channel, into
//! dst (nmb, nch, k, j, i)

void RadiationVet::RssRstPack(DvceArray5D<Real> &dst, int nmb, int mode, int m0) {
  const int nch = RssRstNch(mode);
  for (int n=0; n<nch; ++n) {
    DvceArray5D<Real> *a = nullptr;
    int c = 0;
    RssSlot(mode, n, a, c);
    Kokkos::deep_copy(Kokkos::subview(dst, Kokkos::ALL, n, Kokkos::ALL, Kokkos::ALL,
                                      Kokkos::ALL),
                      Kokkos::subview(*a, std::make_pair(m0, m0 + nmb), c, Kokkos::ALL,
                                      Kokkos::ALL, Kokkos::ALL));
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RssPrime
//! \brief before the first dt with <hydro>/rad_signal_speed: put back the signal-speed
//! inputs a restart file carries, else fill opac from the current state (the closure
//! tensor starts from what the module has)

void RadiationVet::RssPrime(bool restart) {
  if (rss_stage_mode >= 16) {
    // extended record; the tensor slots only into the same closure
    const int mode = rss_stage_mode, b = mode - 16, nb = RssRstNch(b);
    const bool tens = (b == 5 && RssMode() == 5) ||
                      (b == 4 && vet_sc && vet_cell.extent(0) > 0 && vet_full);
    const int nmb = static_cast<int>(rss_stage.extent(0));
    const int nall = static_cast<int>(rss_stage.extent(1));
    const int nop = opac.extent_int(1) - 1;
    for (int n=0; n<nall; ++n) {
      const bool isop = (n == 0) || (n >= nb && n < nb + nop);
      if (!isop && !tens) continue;
      DvceArray5D<Real> *a = nullptr;
      int c = 0;
      RssSlot(mode, n, a, c);
      if (a->extent_int(1) <= c) continue;
      Kokkos::deep_copy(Kokkos::subview(*a, std::make_pair(0,nmb), c, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL),
                        Kokkos::subview(rss_stage, Kokkos::ALL, n, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
    }
    rss_stage = DvceArray5D<Real>();
    rss_stage_mode = -1;
    return;
  }
  if (rss_stage_mode >= 0) {
    const int mode = rss_stage_mode;
    // the tensor is put back only into the same kind of closure; opac always
    int want = 1;
    if ((mode == 4 && vet_sc && vet_cell.extent(0) > 0 && vet_full) ||
        (mode == 5 && RssMode() == 5)) {
      want = RssRstNch(mode);
    } else if (mode >= 2 && global_variable::my_rank == 0) {
      std::cout << "### WARNING: <rad_vet> the restart file's signal-speed closure "
                   "tensor "
                << "does not match this run's closure; only the opacity is restored"
                << std::endl;
    }
    const int nmb = static_cast<int>(rss_stage.extent(0));
    for (int n=0; n<want; ++n) {
      DvceArray5D<Real> *a = nullptr;
      int c = 0;
      RssSlot(mode, n, a, c);
      Kokkos::deep_copy(Kokkos::subview(*a, std::make_pair(0,nmb), c, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL),
                        Kokkos::subview(rss_stage, Kokkos::ALL, n, Kokkos::ALL,
                                        Kokkos::ALL, Kokkos::ALL));
    }
    rss_stage = DvceArray5D<Real>();
    rss_stage_mode = -1;
    return;
  }
  if (restart && global_variable::my_rank == 0) {
    std::cout << "### WARNING: restart file has no <hydro|mhd>/rad_signal_speed inputs; "
              << "the opacity is rebuilt from the restored state and this restart is not "
              << "bitwise." << std::endl;
  }
  (void) Opacity(nullptr, 1);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::RSLACheck
//! \brief reduced-speed-of-light check, once on the first opacity fill: prints
//! K = chat/(v_max tau_max); tau_max: largest x1 column depth, v_max: rsla_vmax or max|v|

void RadiationVet::RSLACheck() {
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
  Kokkos::parallel_reduce("vet_rsla_tau",
  Kokkos::MDRangePolicy<Kokkos::Rank<3>>(DevExeSpace(), {0,ks,js}, {nmb1+1,ke+1,je+1}),
  KOKKOS_LAMBDA(const int m, const int k, const int j, Real &lmax) {
    Real tau = 0.0;
    for (int i=is; i<=ie; ++i) {
      tau += opac_(m,VET_OP_T,k,j,i)*size.d_view(m).dx1;
    }
    lmax = (tau > lmax) ? tau : lmax;
  }, Kokkos::Max<Real>(taumax));

  Real vmax = 0.0;
  Kokkos::parallel_reduce("vet_rsla_v",
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
      std::cout << "<rad_vet> RSLA check: chat = c, no reduced speed of light (tau_max="
                << taumax << ", measured v_max=" << vmeas << ")" << std::endl;
    }
    return;
  }
  Real denom = vmax*((taumax > 1.0) ? taumax : 1.0);
  Real ratio = (chat > 0.0) ? (denom/chat) : 0.0;
  Real kk = (denom > 0.0) ? (chat/denom) : 0.0;
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> RSLA check: tau_max=" << taumax
              << " v_max=" << vmax << " (measured " << vmeas << ")"
              << " chat=" << chat << " -> K = chat/(v_max tau_max) = " << kk
              << ", v_max tau_max/chat = " << ratio << std::endl;
  }
  if (ratio > 1.0 && !rsla_force) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<rad_vet> reduced speed of light is inadmissible: "
      << "v_max*tau_max/chat = " << ratio << " > 1 (design sect. 2).  Raise "
      << "chat_over_c, or set <rad_vet>/rsla_force = true to run anyway" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (ratio > rsla_warn && global_variable::my_rank == 0) {
    std::cout << "### WARNING: <rad_vet> v_max*tau_max/chat = " << ratio
              << " exceeds <rad_vet>/rsla_warn = " << rsla_warn
              << "; the reduced speed of light is marginal (design sect. 2)"
              << std::endl;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus RadiationVet::NewTimeStep

TaskStatus RadiationVet::NewTimeStep(Driver *pdriver, int stage) {
  // no light-speed CFL unless implicit_cfl > 0 (dt = implicit_cfl dx/c)
  dtnew = static_cast<Real>(std::numeric_limits<float>::max());
  if (!(impl_cfl > 0.0)) return TaskStatus::complete;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int nx1 = indcs.nx1;
  int nx2 = indcs.nx2;
  int nx3 = indcs.nx3;

  Real dt1 = std::numeric_limits<float>::max();
  Real dt2 = std::numeric_limits<float>::max();
  Real dt3 = std::numeric_limits<float>::max();

  auto &size = pmy_pack->pmb->mb_size;
  const int nmkji = (pmy_pack->nmb_thispack)*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;

  Kokkos::parallel_reduce("RadVetNudt", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
    int m = (idx)/nkji;
    min_dt1 = fmin((size.d_view(m).dx1), min_dt1);
    min_dt2 = fmin((size.d_view(m).dx2), min_dt2);
    min_dt3 = fmin((size.d_view(m).dx3), min_dt3);
  }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2), Kokkos::Min<Real>(dt3));

  Real dxmin = dt1;
  if (pmy_pack->pmesh->multi_d) {dxmin = std::min(dxmin, dt2);}
  if (pmy_pack->pmesh->three_d) {dxmin = std::min(dxmin, dt3);}
  if (sph_geom) {
    // spherical: mb_size dx2/dx3 are angles, so take the smallest physical width
    // (dr, r dtheta, r sin(theta) dphi) over the active cells
    auto cx1 = pmy_pack->pcoord->dx1;
    auto cx2 = pmy_pack->pcoord->dx2;
    auto cx3 = pmy_pack->pcoord->dx3;
    const int is = indcs.is, js = indcs.js, ks = indcs.ks;
    const int nji = nx2*nx1;
    const bool md = pmy_pack->pmesh->multi_d, td = pmy_pack->pmesh->three_d;
    Real smin = std::numeric_limits<float>::max();
    Kokkos::parallel_reduce("RadVetNudtSph",
                            Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &lmin) {
      int m = idx/nkji;
      int r = idx - m*nkji;
      int k = r/nji;
      r -= k*nji;
      int j = r/nx1;
      int i = r - j*nx1;
      k += ks; j += js; i += is;
      Real d = cx1(m,k,j,i);
      if (md) {d = fmin(d, cx2(m,k,j,i));}
      if (td) {d = fmin(d, cx3(m,k,j,i));}
      lmin = fmin(d, lmin);
    }, Kokkos::Min<Real>(smin));
    dxmin = smin;
  }

  dtnew = impl_cfl*dxmin/c_light;
  return TaskStatus::complete;
}

}  // namespace radvet
