//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file hydro_ringavg.cpp
//! \brief conservative polar ring averaging (see hydro_ringavg.hpp).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/coordinates.hpp"
#include "hydro/hydro_ringavg.hpp"

namespace hydro {
namespace {
[[noreturn]] void RavgFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
            << "<mesh>/polar_ring_avg: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}
}  // namespace

RingAvg::RingAvg(MeshBlockPack *ppack, ParameterInput *pin) : pmy_pack(ppack) {
  Mesh *pm = ppack->pmesh;
  if (!pm->use_spherical_polar || !pm->use_polar_boundary) {
    RavgFatal("needs use_spherical_polar and use_polar_boundary");
  }
  if (!pm->three_d) RavgFatal("needs a 3-D mesh");
  if (pm->multilevel) RavgFatal("SMR/AMR not supported");
  auto &ic = pm->mb_indcs;
  const int nrow = pm->mesh_indcs.nx2;   // theta rows of the whole mesh
  const int nphi = pm->mesh_indcs.nx3;

  // the group size of a row doubles while 2 n sin(theta) <= tol
  const Real tol = 1.05;

  const int nmb = ppack->nmb_thispack;
  const int n2 = ic.nx2 + 2*ic.ng;
  Kokkos::realloc(nfac, nmb, n2);
  auto hfac = Kokkos::create_mirror_view(nfac);
  auto hx2v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ppack->pcoord->x2v);
  std::vector<int> vm, vj, vn, vjg;
  std::vector<int> nrep(nrow, 1);   // per global row (report)
  for (int m = 0; m < nmb; ++m) {
    const int lx2 = static_cast<int>(pm->lloc_eachmb[ppack->gids + m].lx2);
    for (int j = 0; j < n2; ++j) {
      hfac(m,j) = 1.0;
      if (j < ic.js || j > ic.je) continue;
      const int jg = lx2*ic.nx2 + (j - ic.js);
      const int jp = std::min(jg, nrow - 1 - jg);   // row from the nearer pole
      // the doubling is capped by the half-ring (a mesh property), not by the
      // MeshBlock nx3, so the group size does not depend on the phi decomposition; a
      // block nx3 that the mesh-wide group size does not divide is refused below
      int n = 1;
      const Real st = std::sin(hx2v(m,j));
      while (2*n*st <= tol && ((nphi/2) % (2*n)) == 0) n *= 2;
      if ((ic.nx3 % n) != 0 || ((nphi/2) % n) != 0) {
        RavgFatal("group size " + std::to_string(n) + " (row " + std::to_string(jp) +
                  " from the pole) must divide the MeshBlock nx3 = " +
                  std::to_string(ic.nx3) + " and the half-ring nx3/2 = " +
                  std::to_string(nphi/2));
      }
      if (n > 1) {
        hfac(m,j) = static_cast<Real>(n);
        vm.push_back(m); vj.push_back(j); vn.push_back(n); vjg.push_back(jg);
        nrep[jg] = n;
      }
    }
  }
  Kokkos::deep_copy(nfac, hfac);
  nent = static_cast<int>(vm.size());
  if (nent > 0) {
    Kokkos::realloc(em, nent);
    Kokkos::realloc(ej, nent);
    Kokkos::realloc(en, nent);
    auto hm = Kokkos::create_mirror_view(em);
    auto hj = Kokkos::create_mirror_view(ej);
    auto hn = Kokkos::create_mirror_view(en);
    for (int e = 0; e < nent; ++e) {
      hm(e) = vm[e];
      hj(e) = vj[e];
      hn(e) = vn[e];
    }
    Kokkos::deep_copy(em, hm);
    Kokkos::deep_copy(ej, hj);
    Kokkos::deep_copy(en, hn);
  }
  if (global_variable::my_rank == 0) {
    std::cout << "RINGAVG polar ring averaging: group size per row from the pole:";
    for (int jp = 0; jp < nrow/2 && nrep[jp] > 1; ++jp) std::cout << " " << nrep[jp];
    std::cout << std::endl;
  }
  // <mesh>/polar_ring_avg_filter, read only when named
  if (pin->DoesParameterExist("mesh", "polar_ring_avg_filter")) {
    const std::string f = pin->GetString("mesh", "polar_ring_avg_filter");
    if (f == "fourier") {
      fourier = true;
      FourierInit(pin, vjg);
    } else if (f != "average") {
      RavgFatal("polar_ring_avg_filter must be average or fourier");
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn RingAvg::Average
//! \brief replace each phi group of the averaged rows by its group average (see the file
//! header). One thread per (group, i); sums in a fixed order: deterministic.

void RingAvg::Average(DvceArray5D<Real> &u, const int nvar, const DvceArray4D<Real> &phi,
                      const bool grav) {
  if (fourier) {
    Filter(u, nvar, phi, grav, nullptr, 0, 0.0);
    return;
  }
  if (nent == 0) return;
  auto &ic = pmy_pack->pmesh->mb_indcs;
  const int is = ic.is, ie = ic.ie, ks = ic.ks, ke = ic.ke;
  auto u_ = u;
  auto phi_ = phi;
  auto em_ = em, ej_ = ej, en_ = en;
  par_for("ring_avg", DevExeSpace(), 0, nent-1, ks, ke, is, ie,
  KOKKOS_LAMBDA(const int e, const int k, const int i) {
    const int n = en_(e);
    if (((k - ks) % n) != 0) return;
    const int m = em_(e), j = ej_(e);
    const Real rn = 1.0/static_cast<Real>(n);
    for (int v = 0; v < nvar; ++v) {
      Real s = 0.0;
      for (int kk = k; kk < k + n; ++kk) {
        s += u_(m,v,kk,j,i);
      }
      s *= rn;
      for (int kk = k; kk < k + n; ++kk) {
        u_(m,v,kk,j,i) = s;
      }
    }
    // etotgrav: IEN carries rho Phi with the cell's own Phi
    if (grav) {
      Real pb = 0.0;
      for (int kk = k; kk < k + n; ++kk) {
        pb += phi_(m,kk,j,i);
      }
      pb *= rn;
      const Real db = u_(m,IDN,k,j,i);
      for (int kk = k; kk < k + n; ++kk) {
        u_(m,IEN,kk,j,i) += db*(phi_(m,kk,j,i) - pb);
      }
    }
  });
}

}  // namespace hydro
