//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file sc_cart_setup.cpp
//! \brief Cartesian short characteristics: keys, direction sets, buffers of several
//! MeshBlocks.
//!
//!  called from: rad_vet.cpp, rad_vet_closure_faces.cpp, rad_vet_closure_ring.cpp,
//!    rad_vet_params_implicit.cpp, rays/sc_cart_sweep.cpp, rays/sc_cart_sweep_mb.cpp,
//!    rays/sc_cs_sweep.cpp, rays/sc_sph_sweep.cpp
//!  calls into: rad_vet_geometry.cpp, rad_vet_params_closure.cpp, rays/sc_directions.cpp

#include <math.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rays/sc_cart_exchange.hpp"
#include "utils/recipe.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace radvet {

namespace {

// Gauss-Legendre nodes and weights on [0,1] (weights sum to 1)
void GaussLegendre01(const int n, std::vector<double> &x, std::vector<double> &w) {
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(M_PI*(i + 0.75)/(n + 0.5));
    double pp = 1.0;
    for (int it = 0; it < 200; ++it) {
      double p1 = 1.0, p2 = 0.0;
      for (int jj = 1; jj <= n; ++jj) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0*jj - 1.0)*z*p2 - (jj - 1.0)*p3)/jj;
      }
      pp = n*(z*p1 - p2)/(z*z - 1.0);
      double z1 = z;
      z = z1 - p1/pp;
      if (std::fabs(z - z1) < 1.0e-15) break;
    }
    x[i] = 0.5*(1.0 + z);
    w[i] = 1.0/((1.0 - z*z)*pp*pp);
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VetBathLo
//! \brief isotropic intensity entering the sweep at the inner x1 face when vet_bc_inner
//! resolves to emit: Marshak bath, 4 F/c of a flux face, 0 at outflow; -1 = diffusion

Real RadiationVet::VetBathLo() const {
  if (!VetBcEmit()) {return -1.0;}
  if (ibc_x1min == VET_IBC_OUTFLOW) return 0.0;
  return (ibc_x1min == VET_IBC_FLUX) ? 4.0*iflux_x1min/c_light : iebath_x1min;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetInit
//! \brief read the vet_sc keys, build the angle set and the sweep buffers

void RadiationVet::VetInit(ParameterInput *pin) {
  Mesh *pm = pmy_pack->pmesh;
  if (!trans_on) {
    VetFatal("<rad_vet>/closure = sc needs transport = implicit on a MULTI-D mesh");
  }
  // several MeshBlocks or ranks: the banded sweep VetSweepMB
  const bool mbpath = (pm->nmb_total != 1 || global_variable::nranks != 1);
  // a non-periodic x2 (x3) boundary is a vacuum face of the formal solution (I = 0
  // entering); single MeshBlock only
  vet_x2open = (pm->mesh_bcs[BoundaryFace::inner_x2] != BoundaryFlag::periodic ||
                pm->mesh_bcs[BoundaryFace::outer_x2] != BoundaryFlag::periodic);
  vet_x3open = trans_x3 &&
               (pm->mesh_bcs[BoundaryFace::inner_x3] != BoundaryFlag::periodic ||
                pm->mesh_bcs[BoundaryFace::outer_x3] != BoundaryFlag::periodic);
  // the moment solve's transverse faces must be the same vacuum (marshak | outflow)
  auto vac = [](const int b) {return b == VET_IBC_MARSHAK || b == VET_IBC_OUTFLOW;};
  if ((vet_x2open && !(vac(ibc_x2min) && vac(ibc_x2max))) ||
      (vet_x3open && !(vac(ibc_x3min) && vac(ibc_x3max)))) {
    VetFatal("<rad_vet>/closure = sc on a non-periodic x2 (x3) mesh needs "
             "implicit_bc_x2min = implicit_bc_x2max (x3min, x3max) = marshak: the "
             "vacuum face of both the formal solution and the moment solve");
  }
  if ((vet_x2open || vet_x3open) && mbpath) {
    VetFatal("<rad_vet>/closure = sc with non-periodic x2 (x3) boundaries (vacuum "
             "faces of the formal solution) needs a single MeshBlock on one rank");
  }
  vet_nmu = pin->GetOrAddInteger("rad_vet", "vet_nmu", 4);
  vet_nphi = pin->GetOrAddInteger("rad_vet", "vet_nphi", 8);
  // vet_bc_inner / vet_bc_reflect (VetBcInit); vet_bc_bath: an x1 end with
  // implicit_bc = marshak gives its entering rays the bath implicit_ebath
  VetBcInit(pin, false);
  vet_bc_bath = pin->GetOrAddBoolean("rad_vet", "vet_bc_bath",
      recipe::Fix(pin, "rad_vet", "vet_bc_bath", "true", "false"));
  // vet_bc_diff: E of the bottom diffusion intensity from the cell, or moved to the face
  // by the source's linear extrapolation, E + (S_f - S_0), S_f = 1.5 S_0 - 0.5 S_1
  {
    std::string bd = pin->GetOrAddString("rad_vet", "vet_bc_diff",
        recipe::Fix(pin, "rad_vet", "vet_bc_diff", "face", "cell") ? "face" : "cell");
    if (bd.compare("face") == 0) {
      vet_bc_face = true;
    } else if (bd.compare("cell") == 0) {
      vet_bc_face = false;
    } else {
      VetFatal("<rad_vet>/vet_bc_diff must be cell or face");
    }
  }
  {std::string ax = pin->GetOrAddString("rad_vet", "vet_axis", "flux");
  if (ax.compare("flux") == 0) {
    vet_axis_flux = true;
  } else if (ax.compare("eigen") == 0) {
    vet_axis_flux = false;
  } else {
    VetFatal("<rad_vet>/vet_axis = '" + ax + "' not implemented (flux | eigen)");
  }
  }
  vet_full = true;
  if (vet_nmu < 1 || vet_nphi < 1) {VetFatal("<rad_vet>/vet_nmu, vet_nphi must be >= 1");}
  // vet_x1_periodic: periodic x1 sweep, up to vet_x1_npass passes, each entering from
  // the exit plane of the pass before; only the last pass accumulates
  vet_x1per = false;
  vet_x1npass = 1;
  // vet_x1_npass is also the pass count of the mirrored cavity (vet_mirror_npass, >= 2)
  if (pin->DoesParameterExist("rad_vet", "vet_x1_npass")) {
    vet_mirror_npass = std::max(2, pin->GetInteger("rad_vet", "vet_x1_npass"));
  }
  // default true on a fresh run with a periodic x1 mesh: the open sweep is not periodic
  const bool x1meshper =
      (pm->mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic);
  if (pin->DoesParameterExist("rad_vet", "vet_x1_periodic")) {
    vet_x1per = pin->GetBoolean("rad_vet", "vet_x1_periodic");
  } else if (x1meshper && !global_variable::restart_run) {
    vet_x1per = true;
    pin->SetBoolean("rad_vet", "vet_x1_periodic", vet_x1per);
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_sc: periodic x1 mesh, vet_x1_periodic defaulted to true"
                << std::endl;
    }
  }
  if (vet_x1per) {
    if (pm->mesh_bcs[BoundaryFace::inner_x1] != BoundaryFlag::periodic) {
      VetFatal("<rad_vet>/vet_x1_periodic needs periodic x1 boundaries");
    }
    const bool np_named = pin->DoesParameterExist("rad_vet", "vet_x1_npass");
    vet_x1npass = pin->GetOrAddInteger("rad_vet", "vet_x1_npass", 64);
    if (vet_x1npass < 1) {VetFatal("<rad_vet>/vet_x1_npass must be >= 1");}
    // vet_x1_ptol: adaptive passes.  d_p = max|I_p - I_{p-1}|/max I_p of the exit plane,
    // q = d_p/d_{p-1}; stop when d_p/(1 - q) < vet_x1_ptol (at most vet_x1_npass).
    // 0 = exactly vet_x1_npass passes.  Later builds start from the last exit plane.
    const bool fixed = np_named || global_variable::restart_run || mbpath;
    vet_x1ptol = pin->GetOrAddReal("rad_vet", "vet_x1_ptol", fixed ? 0.0 : 1.0e-6);
    if (vet_x1ptol < 0.0) {VetFatal("<rad_vet>/vet_x1_ptol must be >= 0");}
    if (global_variable::my_rank == 0) {
      std::cout << "<rad_vet> vet_sc periodic x1 sweep: "
                << ((vet_x1ptol > 0.0) ? "adaptive, vet_x1_ptol=" : "fixed passes")
                << ((vet_x1ptol > 0.0) ? std::to_string(vet_x1ptol) : std::string())
                << ", vet_x1_npass=" << vet_x1npass
                << ((vet_x1ptol > 0.0) ? " (cap)" : "") << std::endl;
    }
    // vet_x1_periodic_inflow: entry intensity of the first pass.
    //   source:   both directions enter with the face source S (uniform state: I = S,
    //             H = 0 exactly);
    //   boundary: the ordinary boundary intensities (one-sided: a spurious H ~ J at
    //             small npass).
    std::string vin = pin->GetOrAddString("rad_vet", "vet_x1_periodic_inflow",
        recipe::Fix(pin, "rad_vet", "vet_x1_periodic_inflow", "source", "boundary")
            ? "source" : "boundary");
    if (vin.compare("source") == 0) {
      vet_x1src = true;
    } else if (vin.compare("boundary") == 0) {
      vet_x1src = false;
    } else {
      VetFatal("<rad_vet>/vet_x1_periodic_inflow = '" + vin
               + "' is not a choice (source | boundary)");
    }
  }
  vet_nray = 2*vet_nmu*vet_nphi;
  // vet_ring = N: the 2-D ring of N in-plane directions theta = (l + 1/2) 2 pi/N from +x1
  // toward +x2 (N a multiple of 4), equal weights.  vet_ring_groups = 2 adds the tilted
  // rings sqrt(1/3) (cos, sin) theta, +-sqrt(2/3) in x3 (weights 1/2 : 1/4 : 1/4), giving
  // f = 1/3 for an isotropic field.
  int vring = 0, vrgrp = 1;
  if (pin->DoesParameterExist("rad_vet", "vet_ring")) {
    vring = pin->GetInteger("rad_vet", "vet_ring");
    vrgrp = pin->GetOrAddInteger("rad_vet", "vet_ring_groups", 1);
    if (vring < 4 || vring % 4 != 0)
      VetFatal("<rad_vet>/vet_ring must be a multiple of 4");
    if (vrgrp != 1 && vrgrp != 2) VetFatal("<rad_vet>/vet_ring_groups must be 1 or 2");
    vet_nray = vring*((vrgrp == 2) ? 3 : 1);
  }

  vet_ring_n = vring;
  vet_ring_g = vrgrp;
  VetBuildAngles(0.0);
  auto &indcs = pm->mb_indcs;
  int nmb = std::max((pmy_pack->nmb_thispack), (pm->nmb_maxperrank));
  int ncells1 = indcs.nx1 + 2*(indcs.ng);
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(vet_cell, nmb, VETRAY_NC + VETRAY_NKP +
                  (cln_consist ? VETRAY_NFS : 0), ncells3, ncells2, ncells1);
  Kokkos::deep_copy(vet_cell, 0.0);
  if (mbpath) {
    VetMBInit(pin);   // the banded plane buffers live in vet_mbs
  } else {
    Kokkos::realloc(vet_ipl, nmb, 2, vet_nray, ncells3, ncells2);
    Kokkos::deep_copy(vet_ipl, 0.0);
    // device: one thread per (block, ray, cell), then per-cell sums; host: one thread per
    // cell over the rays; same bits
    vsc_rayk = !std::is_same<DevExeSpace, Kokkos::DefaultHostExecutionSpace>::value;
    if (vsc_rayk) {Kokkos::realloc(vet_rayg, nmb, 5*vet_nray, ncells3, ncells2);}
    if (vet_x1per && vet_x1ptol > 0.0) {
      Kokkos::realloc(vet_x1ent, nmb, vet_nray, ncells3, ncells2);
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << "<rad_vet> closure = sc: short characteristics, " << vet_nray
              << " rays (double-Gauss nmu=" << vet_nmu << " x nphi=" << vet_nphi
              << "), lagged per hydro step"
              << (vet_full ? ", FULL tensor D = K/J (eigenvalue floor " : "")
              << (vet_full ? std::to_string(vet_eig_min) + ")" : "")
              << ", vet_bc_diff = " << (vet_bc_face ? "face" : "cell")
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetBuildAngles
//! \brief the vet_sc angle set (mu1, mu2, mu3, weight), up hemisphere first, turned by
//! beta about x1 (beta = 0: the fixed set)

void RadiationVet::VetBuildAngles(const Real beta) {
  const int vring = vet_ring_n, vrgrp = vet_ring_g;
  // the angle set: (mu1, mu2, mu3, weight), upward hemisphere first
  std::vector<double> xg, wg;
  GaussLegendre01(vet_nmu, xg, wg);
  if (vet_ang.extent_int(0) != vet_nray) {Kokkos::realloc(vet_ang, vet_nray, 4);}
  auto ah = Kokkos::create_mirror_view(vet_ang);
  int r = 0;
  for (int h = 0; h < ((vring > 0) ? 0 : 2); ++h) {
    for (int a = 0; a < vet_nmu; ++a) {
      double m1 = (h == 0) ? xg[a] : -xg[a];
      double st = std::sqrt(1.0 - m1*m1);
      for (int l = 0; l < vet_nphi; ++l) {
        double ph = 2.0*M_PI*(l + 0.5)/vet_nphi;
        if (beta != 0.0) {ph += beta;}
        ah(r, 0) = m1;
        ah(r, 1) = st*std::cos(ph);
        ah(r, 2) = st*std::sin(ph);
        ah(r, 3) = 0.5*wg[a]/vet_nphi;
        ++r;
      }
    }
  }
  if (vring > 0) {   // the 2-D ring (up rays first, ray r + nray/2 its
                     // x1 mirror)
    std::vector<std::array<double, 4>> up, dn;
    for (int g = 0; g < ((vrgrp == 2) ? 3 : 1); ++g) {
      const double sz = (g == 0) ? 0.0 : std::sqrt(2.0/3.0)*((g == 1) ? 1.0 : -1.0);
      const double sp = (g == 0) ? 1.0 : std::sqrt(1.0/3.0);
      const double wg0 = (vrgrp == 1) ? 1.0 : ((g == 0) ? 0.5 : 0.25);
      for (int l = 0; l < vring; ++l) {
        double th = 2.0*M_PI*(l + 0.5)/vring;
        if (beta != 0.0) {th += beta;}
        if (std::cos(th) > 0.0) {
          up.push_back({sp*std::cos(th), sp*std::sin(th), sz, wg0/vring});
        } else {
          dn.push_back({sp*std::cos(th), sp*std::sin(th), sz, wg0/vring});
        }
      }
    }
    const int nhu = static_cast<int>(up.size());
    if (beta != 0.0 && static_cast<int>(dn.size()) != nhu) {
      VetFatal("<rad_vet>/vet_ring: the rotated ring lost its up/down balance");
    }
    for (int q = 0; q < nhu; ++q) {
      for (int c = 0; c < 4; ++c) {
        ah(q, c) = up[q][c];
        // fixed ring: down ray = x1 mirror of the up ray; rotated ring: its own down rays
        ah(q + nhu, c) = (beta != 0.0) ? dn[q][c] : ((c == 0) ? -up[q][c] : up[q][c]);
      }
    }
  }
  Kokkos::deep_copy(vet_ang, ah);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRotUpdate
//! \brief vet_sc set of the current rotation block (VetGdBlock, the step's midpoint);
//! double-Gauss: azimuth q 2 pi/nphi; ring: 0.8 (q - 1/2) 2 pi/N (rays >= 0.1 ring step
//! from +-x2).  Rebuilt only when the angle changes.

void RadiationVet::VetRotUpdate() {
  if (!vsc_rot) {return;}
  const Real beta = VetRotBeta(VetGdBlock(pmy_pack->pmesh->ncycle));
  if (beta == vsc_beta) {return;}
  VetRotSet(beta);
}

//----------------------------------------------------------------------------------------
//! \fn Real RadiationVet::VetRotBeta
//! \brief vet_sc set rotation of block blk

Real RadiationVet::VetRotBeta(const int64_t blk) const {
  // a mirrored x1 end needs x1 mirror pairs: a turned ring has none, so it stays fixed
  if (vet_ring_n > 0 && !vet_x1per && ((vet_bc_mirror && ibc_x1min == VET_IBC_REFLECT) ||
      (vet_bc_mirror_hi && ibc_x1max == VET_IBC_REFLECT))) {
    return 0.0;
  }
  // stratified phase with vet_rotate_avg (RotPhase)
  const Real q = vra_on ? static_cast<Real>(RotPhase(blk))
                        : VetGdAngleBlk(blk)/(0.5*M_PI);
  return (vet_ring_n > 0) ? 0.8*(q - 0.5)*2.0*M_PI/vet_ring_n : q*2.0*M_PI/vet_nphi;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetRotSet
//! \brief vet_sc set turned by beta

void RadiationVet::VetRotSet(const Real beta) {
  vsc_beta = beta;
  VetBuildAngles(beta);
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetScRotAvgSetup
//! \brief register the vet_sc fields averaged by vet_rotate_avg: J (ratio) and, as X/J,
//! K_ab, H_a, the half-range sums and the C1 sum; the projection follows from the
//! averaged K/J

void RadiationVet::VetScRotAvgSetup() {
  RotAvgTab &t = vra_tab;
  t.nf = 0;
  vra_a0 = vet_cell;
  vra_a1 = vet_cell;
  vra_c = DvceArray3D<Real>();
  auto add = [&](const int comp, const int kind, const Real lo, const Real hi) {
    t.arr[t.nf] = 0; t.comp[t.nf] = comp; t.kind[t.nf] = kind;
    t.lo[t.nf] = lo; t.hi[t.nf] = hi; ++t.nf;
  };
  t.jf = 0;
  t.jw = 0;   // the sets weighted by their J (union of the sets)
  add(VETRAY_J, 2, 0.0, 0.0);
  for (int c = 0; c < 3; ++c) {add(VETRAY_K11 + c, 3, 0.0, 1.0);}
  for (int c = 3; c < 6; ++c) {add(VETRAY_K11 + c, 3, -1.0, 1.0);}
  for (int c = 0; c < 3; ++c) {add(VETRAY_H1 + c, 3, -1.0, 1.0);}
  if (impl_beam_hr) {
    for (int c = 0; c < 3; ++c) {add(VETRAY_HP1 + c, 3, 0.0, 1.0);}
  }
  if (thsw_c1) {add(VETRAY_NL, 3, 0.0, 1.0e300);}
  t.farr = -1;
  t.fcomp = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetMBInit
//! \brief several MeshBlocks / ranks: band widths, neighbour tables, buffers (uniform
//! mesh, periodic x2 (x3), every block at least 2 cells in x1)

void RadiationVet::VetMBInit(ParameterInput *pin) {
  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  auto &mindcs = pm->mesh_indcs;
  vet_mbs = new VetMBState;
  VetMBState &st = *vet_mbs;
  const int nmb = pmy_pack->nmb_thispack;
  const int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
  const bool thrd = trans_x3;
  if ((mindcs.nx1 % nx1) != 0 || (mindcs.nx2 % nx2) != 0 ||
      (thrd && (mindcs.nx3 % nx3) != 0)) {
    VetFatal("<rad_vet>/closure = sc on several MeshBlocks needs a uniform mesh");
  }
  st.nbx1 = mindcs.nx1/nx1;
  const int nbx2 = mindcs.nx2/nx2;
  const int nbx3 = thrd ? (mindcs.nx3/nx3) : 1;
  st.nx1g = mindcs.nx1;
  if (nx1 < 2) {
    VetFatal("<rad_vet>/closure = sc needs <meshblock>/nx1 >= 2");
  }
  const int me = global_variable::my_rank;
  // sweep blocks: the rank's own; rays: all, up rays (mu1 > 0) first
  st.gsb0 = pm->gids_eachrank[me];
  st.nsb = nmb;
  st.nrl = vet_nray;
  st.nhl = vet_nray/2;

  // band: the largest horizontal reach of a ray per layer, in cells
  auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
  const Real dx1 = MeshBlock0Dx(0), dx2 = MeshBlock0Dx(1), dx3 = MeshBlock0Dx(2);
  Real s2 = 0.0, s3 = 0.0;
  for (int r = 0; r < vet_nray; ++r) {
    s2 = std::max(s2, dx1*std::fabs(ah(r,1))/(std::fabs(ah(r,0))*dx2));
    s3 = std::max(s3, dx1*std::fabs(ah(r,2))/(std::fabs(ah(r,0))*dx3));
  }
  if (vsc_rot) {
    // the band covers every rotated set: |mu_perp|/|mu1| bounds both axes
    // (double-Gauss), the ring's steepest ray at either end of its rotation
    for (int r = 0; r < vet_nray; ++r) {
      const Real tp = std::sqrt(ah(r,1)*ah(r,1) + ah(r,2)*ah(r,2))/std::fabs(ah(r,0));
      if (vet_ring_n == 0) {
        s2 = std::max(s2, dx1*tp/dx2);
        s3 = std::max(s3, dx1*tp/dx3);
      }
    }
    if (vet_ring_n > 0) {
      const Real tm = std::tan(0.5*M_PI - 0.1*2.0*M_PI/vet_ring_n);
      s2 = std::max(s2, dx1*tm/dx2);
      s3 = std::max(s3, dx1*tm/dx3);
    }
  }
  // one cell more for the PLM foot slopes (VetPlmFoot), and one more for a reach within
  // round-off of an integer
  const int wp = 1;
  auto wband = [wp](const Real s) {
    return static_cast<int>(std::floor(s)) + 1 + wp +
           ((std::fabs(s - std::round(s)) < 1.0e-9) ? 1 : 0);
  };
  st.w2 = wband(s2);
  st.w3 = thrd ? wband(s3) : 0;
  st.b2 = st.w2;
  st.b3 = st.w3;
  st.npb = std::max((st.b2 + nx2 - 1)/nx2, thrd ? (st.b3 + nx3 - 1)/nx3 : 1);
  st.n2w = nx2 + 2*st.b2;
  st.n3w = thrd ? (nx3 + 2*st.b3) : 1;
  const int nx3e = thrd ? nx3 : 1;
  if (thrd) {
    st.nof = 8;
    for (int t = 0, o = 0; t < 9; ++t) {
      if (t == 4) continue;
      st.oj[o] = t%3 - 1;
      st.ok[o] = t/3 - 1;
      ++o;
    }
  } else {
    st.nof = 2;
    st.oj[0] = -1;
    st.oj[1] = 1;
    st.ok[0] = st.ok[1] = 0;
  }
  for (int o = st.nof; o < 8; ++o) {
    st.oj[o] = st.ok[o] = 0;
  }
  st.maxreg = thrd ? std::max(std::max(st.b2*nx3, nx2*st.b3), st.b2*st.b3) : st.b2;
  // per-ray band: cells a ray reads beyond the active area per layer, below and above in
  // x2 (x3), one more within round-off of an integer.  With rotation the band of each
  // slot is the union over the whole rotation range (sampled, padded by the largest
  // offset change between samples, capped at w).
  {
    Kokkos::realloc(st.rdep, st.nrl, 4);
    auto rdh = Kokkos::create_mirror_view(st.rdep);
    for (int rl = 0; rl < st.nrl; ++rl) {
      const int r = rl;
      for (int d = 0; d < 2; ++d) {
        rdh(rl,2*d) = 0;
        rdh(rl,2*d + 1) = 0;
        if (d == 1 && !thrd) continue;
        const Real dxp = (d == 0) ? dx2 : dx3;
        const Real sh = -dx1*ah(r,1 + d)/(std::fabs(ah(r,0))*dxp);
        const int o = static_cast<int>(std::floor(sh));
        const int eps = (std::fabs(sh - std::round(sh)) < 1.0e-9) ? 1 : 0;
        rdh(rl,2*d) = std::max(0, -o) + eps + wp;
        rdh(rl,2*d + 1) = std::max(0, o + 1) + eps + wp;
        if (rdh(rl,2*d) > ((d == 0) ? st.w2 : st.w3) ||
            rdh(rl,2*d + 1) > ((d == 0) ? st.w2 : st.w3)) {
          VetFatal("<rad_vet>/closure = sc: per-ray band wider than the band");
        }
      }
    }
    if (vsc_rot) {
      // every set VetRotBeta can give (q in [0, 1])
      const int nsam = 512;
      // a real copy (on a host build create_mirror_view would alias vet_ang)
      auto keep = Kokkos::create_mirror(Kokkos::HostSpace(), vet_ang);
      Kokkos::deep_copy(keep, vet_ang);
      std::vector<Real> smn(2*st.nrl, 1.0e300), smx(2*st.nrl, -1.0e300);
      std::vector<Real> sprev(2*st.nrl, 0.0), dmax(2*st.nrl, 0.0);
      for (int qs = 0; qs <= nsam; ++qs) {
        const Real q = static_cast<Real>(qs)/nsam;
        const Real beta = (vet_ring_n > 0) ? 0.8*(q - 0.5)*2.0*M_PI/vet_ring_n
                                           : q*2.0*M_PI/vet_nphi;
        VetBuildAngles(beta);
        auto bh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
        for (int rl = 0; rl < st.nrl; ++rl) {
          const int r = rl;
          for (int d = 0; d < 2; ++d) {
            if (d == 1 && !thrd) continue;
            const Real dxp = (d == 0) ? dx2 : dx3;
            const Real sh = -dx1*bh(r,1 + d)/(std::fabs(bh(r,0))*dxp);
            const int x = 2*rl + d;
            smn[x] = std::min(smn[x], sh);
            smx[x] = std::max(smx[x], sh);
            if (qs > 0) {dmax[x] = std::max(dmax[x], std::fabs(sh - sprev[x]));}
            sprev[x] = sh;
          }
        }
      }
      Kokkos::deep_copy(vet_ang, keep);   // the current set, bit for bit
      for (int rl = 0; rl < st.nrl; ++rl) {
        for (int d = 0; d < 2; ++d) {
          if (d == 1 && !thrd) continue;
          const int x = 2*rl + d;
          const int w = (d == 0) ? st.w2 : st.w3;
          const int olo = static_cast<int>(std::floor(smn[x] - dmax[x] - 1.0e-9));
          const int ohi = static_cast<int>(std::floor(smx[x] + dmax[x] + 1.0e-9));
          rdh(rl,2*d) = std::min(w, std::max(rdh(rl,2*d), std::max(0, -olo) + wp));
          rdh(rl,2*d + 1) = std::min(w, std::max(rdh(rl,2*d + 1),
                                                 std::max(0, ohi + 1) + wp));
        }
      }
    }
    Kokkos::deep_copy(st.rdep, rdh);
    Kokkos::realloc(st.offs, st.nof, st.nrl + 1);
    Kokkos::realloc(st.offr, st.nof, st.nrl + 1);
    auto osh = Kokkos::create_mirror_view(st.offs);
    auto orh = Kokkos::create_mirror_view(st.offr);
    st.cnts.assign(st.nof, 0);
    st.cntr.assign(st.nof, 0);
    for (int o = 0; o < st.nof; ++o) {
      for (int sd = 0; sd < 2; ++sd) {   // sd = 0: sent (direction -o), 1: received
        const int dj = (sd == 0) ? -st.oj[o] : st.oj[o];
        const int dk = (sd == 0) ? -st.ok[o] : st.ok[o];
        int off = 0;
        for (int r = 0; r < st.nrl; ++r) {
          ((sd == 0) ? osh : orh)(o,r) = off;
          int j0, nj, k0, nk;
          VetRayRegion(dj, nx2, st.b2, rdh(r,0), rdh(r,1), j0, nj);
          VetRayRegion(dk, nx3e, st.b3, rdh(r,2), rdh(r,3), k0, nk);
          off += nj*nk;
        }
        ((sd == 0) ? osh : orh)(o,st.nrl) = off;
        ((sd == 0) ? st.cnts : st.cntr)[o] = off;
      }
    }
    Kokkos::deep_copy(st.offs, osh);
    Kokkos::deep_copy(st.offr, orh);
  }

  // neighbour tables from the global LogicalLocation list over the sweep blocks
  // (gid gsb0 + s); a neighbour on this rank is local (its sweep index)
  const int nbt = st.nbx1*nbx2*nbx3;
  std::vector<int> gof(nbt, -1);
  for (int g = 0; g < pm->nmb_total; ++g) {
    LogicalLocation &l = pm->lloc_eachmb[g];
    gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + l.lx1] = g;
  }
  for (int t = 0; t < nbt; ++t) {
    if (gof[t] < 0) {
      VetFatal("<rad_vet>/closure = sc: incomplete MeshBlock grid");
    }
  }
  const int nsb = st.nsb, g0 = st.gsb0;
  // owner of a block: rank and its sweep index there
  auto cpart = [&](const int g, int &rk, int &lid) {
    rk = pm->rank_eachmb[g];
    lid = g - pm->gids_eachrank[rk];
  };
  Kokkos::realloc(st.lx1, nsb);
  Kokkos::realloc(st.hloc, 8*nsb);
  Kokkos::realloc(st.xloc, 2*nsb);
  st.hrank.assign(8*nsb, me);
  st.hlid.assign(8*nsb, -1);
  st.xrank.assign(2*nsb, -1);
  st.xlid.assign(2*nsb, -1);
  for (int m = 0; m < nsb; ++m) {
    LogicalLocation &l = pm->lloc_eachmb[g0 + m];
    st.lx1.h_view(m) = l.lx1;
    for (int o = 0; o < 8; ++o) {
      st.hloc.h_view(8*m + o) = -1;
      if (o >= st.nof) continue;
      const int l2 = (l.lx2 + st.oj[o] + nbx2) % nbx2;
      const int l3 = (l.lx3 + st.ok[o] + nbx3) % nbx3;
      const int g = gof[(l3*nbx2 + l2)*st.nbx1 + l.lx1];
      if (pm->rank_eachmb[g] == me) {
        st.hloc.h_view(8*m + o) = g - g0;
      } else {
        cpart(g, st.hrank[8*m + o], st.hlid[8*m + o]);
        st.hmpi = true;
      }
    }
    for (int s = 0; s < 2; ++s) {
      int lx = l.lx1 + ((s == 0) ? -1 : 1);
      st.xloc.h_view(2*m + s) = -2;
      if (vet_x1per) lx = (lx + st.nbx1) % st.nbx1;   // the periodic image
      if (lx < 0 || lx >= st.nbx1) continue;
      const int g = gof[(l.lx3*nbx2 + l.lx2)*st.nbx1 + lx];
      if (pm->rank_eachmb[g] == me) {
        st.xrank[2*m + s] = me;
        st.xloc.h_view(2*m + s) = g - g0;
      } else {
        cpart(g, st.xrank[2*m + s], st.xlid[2*m + s]);
        st.xloc.h_view(2*m + s) = -1;
        st.xmpi = true;
      }
    }
  }
  st.lx1.modify_host();
  st.lx1.sync_device();
  st.hloc.modify_host();
  st.hloc.sync_device();
  st.xloc.modify_host();
  st.xloc.sync_device();

  // message plans: a piece is one (sweep block, slot) pair with a remote neighbour,
  // ordered by the receiver's (block, slot) = (hlid, nof-1-o) on the sender, (m, o) on
  // the receiver
  auto plan = [&](VetMBState::Agg &A, const std::vector<int> &csz,
                  const std::vector<int> &crz) {
    std::vector<std::array<int, 4>> ps, pr;   // (rank, key1, key2, index 8m+o)
    for (int m = 0; m < nsb; ++m) {
      for (int o = 0; o < st.nof; ++o) {
        if (st.hloc.h_view(8*m + o) >= 0) continue;
        ps.push_back({st.hrank[8*m + o], st.hlid[8*m + o], st.nof - 1 - o, 8*m + o});
        pr.push_back({st.hrank[8*m + o], m, o, 8*m + o});
      }
    }
    std::sort(ps.begin(), ps.end());
    std::sort(pr.begin(), pr.end());
    Kokkos::realloc(A.soff, 8*nsb);
    Kokkos::realloc(A.roff, 8*nsb);
    auto soh = Kokkos::create_mirror_view(A.soff);
    auto roh = Kokkos::create_mirror_view(A.roff);
    for (int t = 0; t < 8*nsb; ++t) {
      soh(t) = -1;
      roh(t) = -1;
    }
    A.prk.clear();
    A.sdsp.clear();
    A.scnt.clear();
    A.rdsp.clear();
    A.rcnt.clear();
    for (const auto &x : ps) {
      if (A.prk.empty() || A.prk.back() != x[0]) {
        A.prk.push_back(x[0]);
        A.sdsp.push_back(A.stot);
        A.scnt.push_back(0);
      }
      soh(x[3]) = A.stot;
      const int c = csz[x[3]%8];
      A.stot += c;
      A.scnt.back() += c;
    }
    // receive side: the same peers (the neighbour relation is symmetric per slot)
    A.rdsp.assign(A.prk.size(), 0);
    A.rcnt.assign(A.prk.size(), 0);
    size_t ip = 0;
    int last = -1;
    for (const auto &x : pr) {
      while (ip < A.prk.size() && A.prk[ip] != x[0]) ++ip;
      if (ip == A.prk.size()) {
        VetFatal("<rad_vet>/closure = sc: asymmetric neighbour ranks");
      }
      if (static_cast<int>(ip) != last) {
        A.rdsp[ip] = A.rtot;
        last = static_cast<int>(ip);
      }
      roh(x[3]) = A.rtot;
      const int c = crz[x[3]%8];
      A.rtot += c;
      A.rcnt[ip] += c;
    }
    Kokkos::deep_copy(A.soff, soh);
    Kokkos::deep_copy(A.roff, roh);
  };
  {
    std::vector<int> csP(8, 0), crP(8, 0), csB(8, 0), crB(8, 0);
    const int nab = 2*(nx1 + 4);
    for (int o = 0; o < st.nof; ++o) {
      csP[o] = st.cnts[o];
      crP[o] = st.cntr[o];
      const int nj = (st.oj[o] == 0) ? nx2 : st.b2, nk = (st.ok[o] == 0) ? nx3e : st.b3;
      csB[o] = nab*nj*nk;
      crB[o] = nab*nj*nk;
    }
    plan(st.aP, csP, crP);
    plan(st.aB, csB, crB);
    const int ns = std::max(std::max(st.aP.stot, st.aB.stot), 1);
    const int nr = std::max(std::max(st.aP.rtot, st.aB.rtot), 1);
    Kokkos::realloc(st.sfl, ns);
    Kokkos::realloc(st.rfl, nr);
  }

  Kokkos::realloc(st.csw, nsb, nx1 + 4, 2, st.n3w, st.n2w);
  Kokkos::deep_copy(st.csw, 0.0);
  Kokkos::realloc(st.ipl, nsb, st.nring, st.nrl, st.n3w, st.n2w);
  Kokkos::deep_copy(st.ipl, 0.0);
  if (cln_consist) Kokkos::realloc(st.spl, nsb, 4, VETRAY_NFS, st.n3w, st.n2w);
  Kokkos::realloc(vet_ipl, 1, 1, 1, 1, 1);
  if (vet_full) {
    const int n1 = nx1 + 2*indcs.ng;
    const int n2 = nx2 + 2*indcs.ng;
    const int n3 = thrd ? (nx3 + 2*indcs.ng) : 1;
    Kokkos::realloc(st.dsc, nmb, 6, n3, n2, n1);
    Kokkos::deep_copy(st.dsc, 0.0);
    Kokkos::realloc(st.dsc_c, nmb, 6, 1, 1, 1);
    st.pbd = new MeshBoundaryValuesCC(pmy_pack, pin, false);
    st.pbd->InitializeBuffers(6);
    // the D halo follows <mesh>/cs_seam_blockend (read here: pbval_u does not exist yet)
    st.pbd->cs_blockend_ext_ = CsSeamBlockendExtend(pmy_pack, pin);
    // the six D_ab are scalars to the exchange (no vector pairs)
    st.pbd->SetVectorPairs(6, {});
  }
  if (vet_x1per && vet_x1ptol > 0.0) {
    VetFatal("<rad_vet>/vet_x1_periodic on several MeshBlocks needs a fixed "
             "vet_x1_npass");
  }
#if MPI_PARALLEL_ENABLED
  MPI_Comm_dup(MPI_COMM_WORLD, &st.comm);
#endif
  if (global_variable::my_rank == 0) {
    const double mb = 8.0e-6*(static_cast<double>(st.csw.size()) + st.ipl.size()
                              + st.dsc.size() + st.sfl.size() + st.rfl.size());
    std::cout << "<rad_vet> vet_sc on " << pm->nmb_total << " MeshBlocks (" << st.nbx1
              << " x " << nbx2 << " x " << nbx3 << "), " << global_variable::nranks
              << " rank(s): exact banded sweep, reach " << st.w2 << " x " << st.w3
              << " cells, cell-parallel kernel, one message per rank, " << st.nof
              << " horizontal slots, " << mb << " MB of sweep buffers on this rank"
              << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::VetFree

void RadiationVet::VetFree() {
  if (vet_mbs == nullptr) return;
  if (vet_mbs->pbd != nullptr) {
    delete vet_mbs->pbd;
  }
#if MPI_PARALLEL_ENABLED
  MPI_Comm_free(&(vet_mbs->comm));
#endif
  delete vet_mbs;
  vet_mbs = nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationVet::ImplBeamInit
//! \brief boundary-injected collimated beams (Jiang et al. 2012 / Jiang 2021 beam and
//! shadow tests).  implicit_beam_n = 1 or 2; beam b: implicit_beam<b>_e (E_b),
//! _ang (degrees from +x1 toward +x2, |ang| < 90), _y0 / _y1 (x2 patch of the x1min
//! face).  Carried by the up ray nearest the beam direction (2-D: plus its mu3 mirror),
//! entering intensity frac E_b/(n w_r) added to the face's own; the moment solve takes
//! F_in = c frac E_b mu1 as an imposed flux there.  Needs vet_sc, implicit_bc_x1min =
//! marshak, no periodic x1 sweep.

void RadiationVet::ImplBeamInit(ParameterInput *pin) {
  ibm_n = 0;
  if (!pin->DoesParameterExist("rad_vet", "implicit_beam_n")) return;
  ibm_n = pin->GetInteger("rad_vet", "implicit_beam_n");
  if (ibm_n <= 0) {ibm_n = 0; return;}
  if (ibm_n > 2) VetFatal("<rad_vet>/implicit_beam_n must be 0, 1 or 2");
  if (!vet_sc) VetFatal("<rad_vet>/implicit_beam_n needs closure = sc");
  if (ibc_x1min != VET_IBC_MARSHAK)
    VetFatal("<rad_vet>/implicit_beam_n needs implicit_bc_x1min = marshak");
  if (vet_x1per) VetFatal("<rad_vet>/implicit_beam_n is not for a periodic x1 sweep");
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int nx2 = indcs.nx2, nx3 = trans_x3 ? indcs.nx3 : 1;
  const int nmb = pmy_pack->nmb_thispack;
  auto ah = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), vet_ang);
  const int nray = vet_nray;
  ibm_nr = 0;
  Kokkos::realloc(ibm_i, nmb, 4, nx3, nx2);
  Kokkos::realloc(ibm_f, nmb, nx3, nx2);
  Kokkos::realloc(ibm_s, nmb, 3, nx3, nx2);
  auto hi = Kokkos::create_mirror_view(ibm_i);
  auto hf = Kokkos::create_mirror_view(ibm_f);
  auto hs = Kokkos::create_mirror_view(ibm_s);
  Kokkos::deep_copy(hi, 0.0);
  Kokkos::deep_copy(hf, 0.0);
  Kokkos::deep_copy(hs, 0.0);
  auto &msz = pmy_pack->pmb->mb_size;
  for (int b = 0; b < ibm_n; ++b) {
    const std::string p = "implicit_beam" + std::to_string(b) + "_";
    const Real eb = pin->GetReal("rad_vet", p + "e");
    const Real ang = pin->GetReal("rad_vet", p + "ang")*M_PI/180.0;
    const Real y0 = pin->GetReal("rad_vet", p + "y0");
    const Real y1 = pin->GetReal("rad_vet", p + "y1");
    int rb = -1;
    Real best = 1.0e30;
    for (int r = 0; r < nray; ++r) {
      if (!(ah(r,0) > 0.0)) continue;
      // ray nearest the in-plane beam direction (cos ang, sin ang, 0)
      const Real d = -(ah(r,0)*std::cos(ang) + ah(r,1)*std::sin(ang))
                     + 1.0e-9*std::fabs(ah(r,2));
      if (d < best - 1.0e-12) {best = d; rb = r;}
    }
    int rs[2] = {rb, -1};
    int nr = 1;
    if (!trans_x3) {
      for (int r = 0; r < nray; ++r) {
        if (r != rb && std::fabs(ah(r,0) - ah(rb,0)) < 1.0e-12 &&
            std::fabs(ah(r,1) - ah(rb,1)) < 1.0e-12 &&
            std::fabs(ah(r,2) + ah(rb,2)) < 1.0e-12) {
          rs[1] = r;
          nr = 2;
          break;
        }
      }
    }
    for (int q = 0; q < nr; ++q) {
      const int qq = ibm_nr++;
      ibm_ray[qq] = rs[q];
      const Real wr = ah(rs[q],3);
      for (int m = 0; m < nmb; ++m) {
        // beams enter through the x1min face of the mesh
        auto bf = pmy_pack->pmb->mb_bcs.h_view(m,BoundaryFace::inner_x1);
        if (bf == BoundaryFlag::block) continue;
        const Real x2a = msz.h_view(m).x2min, dx2 = msz.h_view(m).dx2;
        for (int k = 0; k < nx3; ++k) {
          for (int j = 0; j < nx2; ++j) {
            const Real ya = x2a + j*dx2, yb = ya + dx2;
            const Real fr = std::fmax(std::fmin(yb, y1) - std::fmax(ya, y0), 0.0)/dx2;
            if (!(fr > 0.0)) continue;
            hi(m,qq,k,j) += fr*eb/(nr*wr);
            hf(m,k,j) += c_light*wr*ah(rs[q],0)*fr*eb/(nr*wr);
            for (int a = 0; a < 3; ++a) {
              hs(m,a,k,j) += c_light*wr*ah(rs[q],a)*ah(rs[q],0)*fr*eb/(nr*wr);
            }
          }
        }
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "  implicit_beam" << b << ": E_b=" << eb << " requested ang="
                << ang*180.0/M_PI << " deg -> ray " << rb << " (mu1 mu2 mu3) = ("
                << ah(rb,0) << " " << ah(rb,1) << " " << ah(rb,2) << "), in-plane ang="
                << std::atan2(ah(rb,1), ah(rb,0))*180.0/M_PI << " deg, " << nr
                << " ray(s), x2 patch [" << y0 << ", " << y1 << "]" << std::endl;
    }
  }
  Kokkos::deep_copy(ibm_i, hi);
  Kokkos::deep_copy(ibm_f, hf);
  Kokkos::deep_copy(ibm_s, hs);
}

}  // namespace radvet
