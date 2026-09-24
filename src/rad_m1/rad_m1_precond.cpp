//========================================================================================
// AthenaK astrophysical fluid dynamics & numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_precond.cpp
//! \brief implicit M1: preconditioner study (tests_m1/runs_5m_precond).
//!
//! <rad_m1>/implicit_dump_op = N (debug, read only when named): at the first implicit
//! solve of cycle N the stored stencil, b and x0 of the frozen system are written, per
//! rank, to m1op.c<N>.r<rank>.bin for the offline study (runs_5m_precond/study.py).

#include <cstdio>
#include <string>
#include <vector>
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitDumpOp
//! \brief file layout (all int32 then float64): nmb, nst, n1, n2, n3; per block gid,
//! lx1, lx2, lx3, level; then per block st(0..nst-1), b, x0 over the active cells
//! (k,j,i order, i fastest)

void RadiationM1::ImplicitDumpOp() {
  auto *pm = pmy_pack->pmesh;
  if (impl_dump_cyc < 0 || pm->ncycle != impl_dump_cyc || impl_dump_done) {return;}
  impl_dump_done = true;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1 = indcs.nx1, n2 = indcs.nx2, n3 = indcs.nx3;
  const int nmb = pmy_pack->nmb_thispack;
  const int nst = static_cast<int>(ost.extent(1));
  auto hst = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ost);
  auto hiw = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), iw);
  char fn[128];
  std::snprintf(fn, sizeof(fn), "m1op.c%d.r%d.bin", impl_dump_cyc,
                global_variable::my_rank);
  std::FILE *f = std::fopen(fn, "wb");
  if (f == nullptr) {return;}
  int hd[5] = {nmb, nst, n1, n2, n3};
  std::fwrite(hd, sizeof(int), 5, f);
  pmy_pack->pmb->mb_gid.sync_host();
  for (int m = 0; m < nmb; ++m) {
    int g = pmy_pack->pmb->mb_gid.h_view(m);
    auto &ll = pm->lloc_eachmb[g];
    int b[5] = {g, ll.lx1, ll.lx2, ll.lx3, ll.level};
    std::fwrite(b, sizeof(int), 5, f);
  }
  std::vector<double> buf(static_cast<size_t>(n1)*n2*n3);
  for (int m = 0; m < nmb; ++m) {
    for (int q = 0; q < nst + 2; ++q) {
      size_t p = 0;
      for (int k = ks; k < ks + n3; ++k) {
        for (int j = js; j < js + n2; ++j) {
          for (int i = is; i < is + n1; ++i) {
            if (q < nst) {
              buf[p++] = hst(m,q,k,j,i);
            } else {
              buf[p++] = hiw(m, (q == nst) ? M1_IW_KB : M1_IW_EP, k, j, i);
            }
          }
        }
      }
      std::fwrite(buf.data(), sizeof(double), buf.size(), f);
    }
  }
  std::fclose(f);
}
