#ifndef UTILS_FOFC_REPORT_HPP_
#define UTILS_FOFC_REPORT_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file fofc_report.hpp
//! \brief <hydro|mhd>/fofc_report = N (DIAGNOSTIC, default off, read only when named):
//! count the active cells where the first-order flux correction FIRES, per RK stage, in N
//! bins of the GLOBAL x1 index (r on the spherical-polar / cubed-sphere grids; the cell
//! index -> r map of a run is in its x1grid file) and split by whether the cell's density
//! entering the stage is at the floor (w0 IDN <= 1.000001 dfloor) or not.  One block of
//! lines per <hydro|mhd>/fofc_report_dt (read only when named; default the <output1> dt,
//! 0 = every cycle) on rank 0, summed over all ranks.  The counts are integers stored
//! in doubles (exact, order-independent).  Reads the flags only: the run is unchanged.

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"

struct FofcReport {
  int nb = 0;                  // number of x1 bins; 0 = off
  Real dt = 0.0, next = 0.0;   // report interval and the next report time
  Real nstage = 0.0;           // RK stages accumulated since the last report
  std::string who;
  DvceArray2D<Real> cnt;       // (nb, 2): [0] density at the floor, [1] not
  DvceArray1D<int> ioff;       // per MeshBlock: global x1 index of its first active cell

  void Init(ParameterInput *pin, const std::string &blk) {
    if (!pin->DoesParameterExist(blk, "fofc_report")) {return;}
    nb = pin->GetInteger(blk, "fofc_report");
    if (nb <= 0) {nb = 0; return;}
    who = blk;
    if (pin->DoesParameterExist(blk, "fofc_report_dt")) {
      dt = pin->GetReal(blk, "fofc_report_dt");
    } else if (pin->DoesParameterExist("output1", "dt")) {
      dt = pin->GetReal("output1", "dt");
    }
    next = -1.0;   // set from the mesh time on the first call
    cnt = DvceArray2D<Real>("fofc_rep", nb, 2);
  }

  // count the flagged active cells of this stage (call after the last flag writer)
  void Count(MeshBlockPack *pp, const DvceArray4D<bool> &flag,
             const DvceArray5D<Real> &w0, const Real dfloor) {
    if (nb == 0) {return;}
    Mesh *pm = pp->pmesh;
    auto &indcs = pm->mb_indcs;
    const int nmb = pp->nmb_thispack;
    if (ioff.extent_int(0) != nmb) {ioff = DvceArray1D<int>("fofc_rep_ioff", nmb);}
    auto hio = Kokkos::create_mirror_view(ioff);
    for (int m = 0; m < nmb; ++m) {
      hio(m) = static_cast<int>(pm->lloc_eachmb[pp->gids + m].lx1)*indcs.nx1;
    }
    Kokkos::deep_copy(ioff, hio);
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    const int nbin = nb, ntot = pm->mesh_indcs.nx1;
    const Real dfl = 1.000001*dfloor;
    auto c_ = cnt;
    auto io_ = ioff;
    par_for("fofc_report", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (!flag(m,k,j,i)) {return;}
      int b = static_cast<int>((static_cast<int64_t>(io_(m) + i - is)*nbin)/ntot);
      b = (b < 0) ? 0 : ((b >= nbin) ? nbin - 1 : b);
      const int c = (w0(m,IDN,k,j,i) <= dfl) ? 0 : 1;
      Kokkos::atomic_add(&c_(b,c), 1.0);   // integer-valued: exact in any order
    });
    nstage += 1.0;
    Report(pm, false);
  }

  void Report(Mesh *pm, const bool force) {
    if (nb == 0) {return;}
    const Real t = pm->time;
    if (next < 0.0) {next = t + dt;}
    if (!force && dt > 0.0 && t < next) {return;}
    while (dt > 0.0 && next <= t) {next += dt;}
    auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cnt);
    std::vector<Real> v(2*nb);
    for (int b = 0; b < nb; ++b) {v[2*b] = h(b,0); v[2*b + 1] = h(b,1);}
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, v.data(), 2*nb, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
#endif
    if (global_variable::my_rank == 0) {
      const int ntot = pm->mesh_indcs.nx1;
      Real s0 = 0.0, s1 = 0.0;
      for (int b = 0; b < nb; ++b) {s0 += v[2*b]; s1 += v[2*b + 1];}
      std::printf("<%s> fofc_report t=%.6e cycle=%d stages=%d fired cell-stages: "
                  "at-floor=%.0f other=%.0f | per x1 bin [i0-i1] floor/other:",
                  who.c_str(), t, pm->ncycle, static_cast<int>(nstage), s0, s1);
      for (int b = 0; b < nb; ++b) {
        if (v[2*b] == 0.0 && v[2*b + 1] == 0.0) {continue;}
        std::printf(" [%d-%d] %.0f/%.0f", (b*ntot)/nb, ((b + 1)*ntot)/nb - 1,
                    v[2*b], v[2*b + 1]);
      }
      std::printf("\n");
      std::fflush(stdout);
    }
    Kokkos::deep_copy(cnt, 0.0);
    nstage = 0.0;
  }
};

#endif  // UTILS_FOFC_REPORT_HPP_
