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
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "reconstruct/plm.hpp"
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

//----------------------------------------------------------------------------------------
//! \fn ImplFatal
//! \brief one fatal-error exit with a message, used by the 3a2 option parsers

void ImplFatal(const std::string &msg) {
  std::cout << "### FATAL ERROR in " << __FILE__ << std::endl << msg << std::endl;
  std::exit(EXIT_FAILURE);
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
  // ---- milestone 3a2 options.  All three default to the 3a behaviour, so an input file
  // that does not name them reproduces RESULTS.txt of runs_3a exactly.
  std::string sfx = pin->GetOrAddString("rad_m1","implicit_flux","central");
  if (sfx.compare("central") == 0) {
    impl_flux = M1_IFLUX_CENTRAL;
  } else if (sfx.compare("ap_hll") == 0) {
    impl_flux = M1_IFLUX_APHLL;
  } else if (sfx.compare("berthon") == 0) {
    impl_flux = M1_IFLUX_BERTHON;
  } else if (sfx.compare("blend") == 0) {
    impl_flux = M1_IFLUX_BLEND;
  } else {
    ImplFatal("<rad_m1>/implicit_flux = '" + sfx
              + "' is not a choice (central | ap_hll | berthon | blend)");
  }
  // ---- milestone 3c.  The weight of implicit_flux = blend.  Inert for every other
  // flux, and the two ends of the blend are BITWISE central and berthon.
  std::string sbl = pin->GetOrAddString("rad_m1","implicit_blend","tau_f");
  if (sbl.compare("tau") == 0) {
    impl_blend = M1_IBLEND_TAU;
  } else if (sbl.compare("f") == 0) {
    impl_blend = M1_IBLEND_F;
  } else if (sbl.compare("tau_f") == 0) {
    impl_blend = M1_IBLEND_TAUF;
  } else {
    ImplFatal("<rad_m1>/implicit_blend = '" + sbl
              + "' is not a choice (tau | f | tau_f)");
  }
  std::string sbm = pin->GetOrAddString("rad_m1","implicit_blend_fmode","max");
  if (sbm.compare("max") == 0) {
    impl_blend_fmode = M1_IBFM_MAX;
  } else if (sbm.compare("mean") == 0) {
    impl_blend_fmode = M1_IBFM_MEAN;
  } else {
    ImplFatal("<rad_m1>/implicit_blend_fmode = '" + sbm
              + "' is not a choice (max | mean)");
  }
  std::string sbw = pin->GetOrAddString("rad_m1","implicit_blend_mode","flux");
  if (sbw.compare("flux") == 0) {
    impl_blend_mode = M1_IBMODE_FLUX;
  } else if (sbw.compare("dissipation") == 0) {
    impl_blend_mode = M1_IBMODE_DISSIP;
  } else {
    ImplFatal("<rad_m1>/implicit_blend_mode = '" + sbw
              + "' is not a choice (flux | dissipation)");
  }
  impl_blend_tau0 = pin->GetOrAddReal("rad_m1","implicit_blend_tau0",1.0);
  impl_blend_flo = pin->GetOrAddReal("rad_m1","implicit_blend_flo",0.6);
  impl_blend_fhi = pin->GetOrAddReal("rad_m1","implicit_blend_fhi",0.9);
  if (!(impl_blend_tau0 > 0.0) || !(impl_blend_fhi > impl_blend_flo)) {
    ImplFatal("<rad_m1>: implicit_blend_tau0 must be positive and implicit_blend_fhi "
              "must exceed implicit_blend_flo");
  }
  std::string srn = pin->GetOrAddString("rad_m1","implicit_recon","dc");
  if (srn.compare("dc") == 0) {
    impl_recon = M1_IRECON_DC;
  } else if (srn.compare("plm_dc") == 0) {
    impl_recon = M1_IRECON_PLMDC;
  } else {
    ImplFatal("<rad_m1>/implicit_recon = '" + srn + "' is not a choice (dc | plm_dc)");
  }
  // LIMIT 4 of the 3a findings is NOT implemented in 3a2: a column still has to live
  // inside one MeshBlock along x1 (the fatal below).  The option is parsed so that the
  // input files and the gate scripts can already name it, and `gather` fatals rather
  // than silently doing something else.
  std::string spt = pin->GetOrAddString("rad_m1","implicit_partition","none");
  if (spt.compare("none") == 0) {
    impl_part = M1_IPART_NONE;
  } else if (spt.compare("gather") == 0) {
    impl_part = M1_IPART_GATHER;
  } else {
    ImplFatal("<rad_m1>/implicit_partition = '" + spt
              + "' is not a choice (none | gather)");
  }
  impl_recon_w = pin->GetOrAddReal("rad_m1","implicit_recon_w",-1.0);
  impl_res_floor = pin->GetOrAddReal("rad_m1","implicit_res_floor",0.0);
  std::string slg = pin->GetOrAddString("rad_m1","implicit_recon_lag","picard");
  impl_recon_freeze = (slg.compare("step") == 0);
  // MILESTONE 3c: freeze the deferred correction, and with it the plm limiter's choice,
  // after this many Picard passes.  The 3a2 finding is that the limiter keeps switching
  // on a handful of cells and the iteration is a small limit cycle that never meets
  // implicit_tol, so plm_dc costs implicit_maxit passes per step; freezing the
  // correction after a few passes leaves an ordinary linear system to converge.
  // <= 0 (the default) never freezes, i.e. reproduces 3a2.
  impl_recon_npass = pin->GetOrAddInteger("rad_m1","implicit_recon_npass",-1);
  if (!impl_recon_freeze && slg.compare("picard") != 0) {
    ImplFatal("<rad_m1>/implicit_recon_lag = '" + slg
              + "' is not a choice (step | picard)");
  }
  // LIMIT 3 of the 3a findings: the DEFAULT is now `true`.  It is a bug fix, not a
  // tuning knob (the boundary face used to hand its whole momentum to one interior cell,
  // which gave that cell 1.5 face-shares of radiative force: He column bottom-cell |v1|
  // 10.3 -> 0.47 v_MLT).  The key is kept so that `false` reproduces runs_3a/RESULTS.txt.
  impl_bmom_half = pin->GetOrAddBoolean("rad_m1","implicit_bmom_half",true);
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
  part_nblk = 1;
  if (mindcs.nx1 != indcs.nx1) {
    if (impl_part != M1_IPART_GATHER) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "<rad_m1>/transport = implicit_x1 with implicit_partition = none "
        << "needs exactly ONE MeshBlock along x1: <meshblock>/nx1 must equal <mesh>/nx1 ("
        << indcs.nx1 << " vs " << mindcs.nx1 << ").  Set <rad_m1>/implicit_partition = "
        << "gather for the line solve partitioned over MeshBlocks and ranks" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if ((mindcs.nx1 % indcs.nx1) != 0) {
      ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nx1 to be an exact "
                "multiple of <meshblock>/nx1 (uniform mesh only)");
    }
    part_nblk = mindcs.nx1/indcs.nx1;
  }
  if (part_nblk > 1 && (ibc_x1min == M1_IBC_PERIODIC)) {
    ImplFatal("<rad_m1>/implicit_partition = gather does not support PERIODIC x1 across "
              "more than one MeshBlock (the cyclic Thomas sweep of 3a wraps inside one "
              "block).  Use one MeshBlock along x1, or a non-periodic x1 boundary pair");
  }
  if (part_nblk > 1 && indcs.ng < 2) {
    ImplFatal("<rad_m1>/implicit_partition = gather needs <mesh>/nghost >= 2");
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
  Kokkos::realloc(ifw, nmb, M1_NIFW, ncells3, ncells2, ncells1+1);
  Kokkos::deep_copy(ifw, 0.0);
  ImplicitPartitionInit();

  if (global_variable::my_rank == 0) {
    std::cout << "<rad_m1>: transport=implicit_x1 (milestone 3a) implicit_cfl="
              << impl_cfl << " tol=" << impl_tol << " maxit=" << impl_maxit
              << " opac_update=" << (impl_opac_update ? "true" : "false")
              << " marshak_q=" << marshak_q << std::endl;
    std::cout << "         implicit_flux="
              << ((impl_flux == M1_IFLUX_APHLL) ? "ap_hll" :
                  ((impl_flux == M1_IFLUX_BERTHON) ? "berthon" :
                   ((impl_flux == M1_IFLUX_BLEND) ? "blend" : "central")))
              << ((impl_flux == M1_IFLUX_BLEND) ? (":" + sbl + "/" + sbm + "/" + sbw)
                                                : std::string(""))
              << " implicit_recon="
              << ((impl_recon == M1_IRECON_PLMDC) ? "plm_dc" : "dc")
              << " implicit_partition="
              << ((impl_part == M1_IPART_GATHER) ? "gather" : "none") << std::endl;
    std::cout << "         x1 boundaries: min=" << ibc_x1min << " max=" << ibc_x1max
              << " (0 marshak, 1 flux, 2 reflect, 3 periodic) flux_min=" << iflux_x1min
              << " flux_max=" << iflux_x1max << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitPartitionInit
//! \brief milestone 3b, LIMIT 4: build the topology of the x1 stacks and allocate the
//! gather/scatter and halo buffers.  A "stack" is the set of MeshBlocks that share the
//! same (x2,x3) footprint, ordered by their x1 logical location; its ROOT is the block
//! with the lowest one.  Uniform mesh only (the caller has already fatalled on SMR/AMR
//! and on a non-integer block count along x1), so the stack of a block is found by a
//! scan of the global LogicalLocation list.
//!
//! Nothing is allocated and nothing is communicated when part_nblk == 1: the 3a path is
//! then bitwise what it was.

void RadiationM1::ImplicitPartitionInit() {
  part_nroot = 0;
  part_nx1g = 0;
  part_nlay = 0;
  part_nqa = M1_NHALO_A;
  part_any_mpi = false;

  Mesh *pm = pmy_pack->pmesh;
  auto &indcs = pm->mb_indcs;
  int nmb = pmy_pack->nmb_thispack;
  int g0 = pmy_pack->gids;
  // the three per-block tables are ALWAYS allocated: ImplicitSolve reads part_pos in
  // every kernel, and with one MeshBlock per column the defaults (position 0 of a stack
  // of one) select exactly the 3a branches.
  Kokkos::realloc(part_pos, nmb);
  Kokkos::realloc(part_slot, nmb);
  Kokkos::realloc(part_nbr, 2*nmb);
  for (int m=0; m<nmb; ++m) {
    part_pos.h_view(m) = 0;
    part_slot.h_view(m) = -1;
    part_nbr.h_view(2*m) = -1;
    part_nbr.h_view(2*m+1) = -1;
  }
  if (part_nblk <= 1) {
    part_pos.modify_host();
    part_pos.sync_device();
    part_slot.modify_host();
    part_slot.sync_device();
    part_nbr.modify_host();
    part_nbr.sync_device();
    return;
  }
  part_nx1g = part_nblk*indcs.nx1;
  part_nlay = std::min(indcs.ng, 2);

  // (1) for every LOCAL block: its position in the stack, the global ids of the stack
  // members and of its two x1 neighbours.
  part_rootgid.assign(nmb, -1);
  part_rootrank.assign(nmb, -1);
  part_nbrrank.assign(2*nmb, -1);
  part_nbrgid.assign(2*nmb, -1);
  std::vector<int> stack(part_nblk);
  for (int m=0; m<nmb; ++m) {
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {stack[p] = -1;}
    for (int g=0; g<pm->nmb_total; ++g) {
      LogicalLocation &lg = pm->lloc_eachmb[g];
      if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level) {
        if (lg.lx1 >= 0 && lg.lx1 < part_nblk) {stack[lg.lx1] = g;}
      }
    }
    for (int p=0; p<part_nblk; ++p) {
      if (stack[p] < 0) {
        ImplFatal("<rad_m1>/implicit_partition = gather: the x1 stack of a MeshBlock is "
                  "incomplete (a non-uniform mesh?)");
      }
    }
    int pos = static_cast<int>(lm.lx1);
    part_pos.h_view(m) = pos;
    part_rootgid[m] = stack[0];
    part_rootrank[m] = pm->rank_eachmb[stack[0]];
    if (pos > 0) {
      part_nbrgid[2*m] = stack[pos-1];
      part_nbrrank[2*m] = pm->rank_eachmb[stack[pos-1]];
    }
    if (pos < part_nblk-1) {
      part_nbrgid[2*m+1] = stack[pos+1];
      part_nbrrank[2*m+1] = pm->rank_eachmb[stack[pos+1]];
    }
    for (int s=0; s<2; ++s) {
      int gn = part_nbrgid[2*m+s];
      bool loc = (gn >= 0) && (part_nbrrank[2*m+s] == global_variable::my_rank);
      part_nbr.h_view(2*m+s) = loc ? (gn - g0) : -1;
      if (gn >= 0 && part_nbrrank[2*m+s] != global_variable::my_rank) {
        part_any_mpi = true;
      }
    }
    if (part_rootrank[m] != global_variable::my_rank) {part_any_mpi = true;}
  }

  // (2) the local ROOTS and their member lists
  part_mrank.clear();
  part_mgid.clear();
  std::vector<int> rootmb;
  for (int m=0; m<nmb; ++m) {
    if (part_pos.h_view(m) != 0) continue;
    rootmb.push_back(m);
    LogicalLocation &lm = pm->lloc_eachmb[g0+m];
    for (int p=0; p<part_nblk; ++p) {
      int gfound = -1;
      for (int g=0; g<pm->nmb_total; ++g) {
        LogicalLocation &lg = pm->lloc_eachmb[g];
        if (lg.lx2 == lm.lx2 && lg.lx3 == lm.lx3 && lg.level == lm.level && lg.lx1 == p) {
          gfound = g;
        }
      }
      part_mgid.push_back(gfound);
      part_mrank.push_back(pm->rank_eachmb[gfound]);
      if (pm->rank_eachmb[gfound] != global_variable::my_rank) {part_any_mpi = true;}
    }
  }
  part_nroot = static_cast<int>(rootmb.size());
  for (int m=0; m<nmb; ++m) {
    int sl = -1;
    if (part_rootrank[m] == global_variable::my_rank) {
      for (int s=0; s<part_nroot; ++s) {
        if (rootmb[s] + g0 == part_rootgid[m]) {sl = s;}
      }
    }
    part_slot.h_view(m) = sl;
  }
  part_pos.modify_host();
  part_pos.sync_device();
  part_slot.modify_host();
  part_slot.sync_device();
  part_nbr.modify_host();
  part_nbr.sync_device();

  // (3) buffers
  int ncells2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int ncells3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  Kokkos::realloc(part_sys, std::max(part_nroot,1), 6, ncells3, ncells2, part_nx1g);
  Kokkos::deep_copy(part_sys, 0.0);
  int nrow = 4*ncells3*ncells2*indcs.nx1;
  Kokkos::realloc(part_sbuf, nmb, nrow);
  Kokkos::realloc(part_rbuf, std::max(part_nroot*part_nblk,1), nrow);
  Kokkos::realloc(part_sbuf_h, nmb, nrow);
  Kokkos::realloc(part_rbuf_h, std::max(part_nroot*part_nblk,1), nrow);
  int nhal = M1_NHALO_A*part_nlay*ncells3*ncells2;
  Kokkos::realloc(part_hbuf, nmb, 4, nhal);
  Kokkos::realloc(part_hbuf_h, nmb, 4, nhal);

  if (global_variable::my_rank == 0) {
    std::cout << "         implicit_partition=gather: " << part_nblk
              << " MeshBlocks per x1 column, " << part_nx1g << " rows per gathered line, "
              << part_nlay << " halo layers" << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitX1Halo
//! \brief exchange the x1 ghost layers of the work array `iw` between the MeshBlocks of
//! one x1 stack.  `eponly` picks the one-quantity set (M1_IW_EP, the new iterate, needed
//! by the face update and by the next pass) instead of the six LAGGED quantities
//! (M1HaloCompA).
//!
//! What makes the partitioned solve BITWISE identical to a single-block solve is that
//! the ghost value a block reads is the VERY NUMBER its neighbour computed, copied, and
//! never a quantity recomputed from a hydro ghost: this routine moves nothing else.

void RadiationM1::ImplicitX1Halo(bool eponly) {
  if (part_nblk <= 1) return;
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  int nl = part_nlay;
  int nq = eponly ? 1 : M1_NHALO_A;
  const bool ep1 = eponly;
  auto iw_ = iw;
  auto nbr_ = part_nbr;
  int nj = je - js + 1, nk = ke - ks + 1;

  // (1) same-rank neighbours: a plain device copy, neighbour ACTIVE -> my GHOST
  par_for("m1_impl_halo_loc", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int mlo = nbr_.d_view(2*m);
    int mhi = nbr_.d_view(2*m+1);
    if (mlo >= 0) {iw_(m,nc,k,j,is-1-l) = iw_(mlo,nc,k,j,ie-l);}
    if (mhi >= 0) {iw_(m,nc,k,j,ie+1+l) = iw_(mhi,nc,k,j,is+l);}
  });
  if (!part_any_mpi) return;

#if MPI_PARALLEL_ENABLED
  // (2) remote neighbours.  Buffer slots: 0 = send to lo, 1 = send to hi,
  // 2 = recv from lo, 3 = recv from hi.
  auto hb_ = part_hbuf;
  par_for("m1_impl_halo_pack", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int idx = (((n*nl + l)*nk + (k-ks))*nj + (j-js));
    hb_(m,0,idx) = iw_(m,nc,k,j,is+l);
    hb_(m,1,idx) = iw_(m,nc,k,j,ie-l);
  });
  int nbuf = nq*nl*nk*nj;
  Kokkos::deep_copy(part_hbuf_h, part_hbuf);
  std::vector<MPI_Request> req;
  // the tag identifies (receiving local block, receiving side, which halo set); the
  // source rank is named in the receive, so it need only be unique per rank pair.
  int te = eponly ? 1 : 0;
  int *gr = pmy_pack->pmesh->gids_eachrank;
  for (int m=0; m<nmb; ++m) {
    for (int s=0; s<2; ++s) {
      int rk = part_nbrrank[2*m+s];
      if (rk < 0 || rk == global_variable::my_rank) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(&part_hbuf_h(m,2+s,0), nbuf, MPI_ATHENA_REAL, rk,
                4*m + 2*s + te, MPI_COMM_WORLD, &req.back());
    }
  }
  for (int m=0; m<nmb; ++m) {
    for (int s=0; s<2; ++s) {
      int rk = part_nbrrank[2*m+s];
      if (rk < 0 || rk == global_variable::my_rank) continue;
      // the neighbour receives this message into ITS slot 2+(1-s)
      int lidn = part_nbrgid[2*m+s] - gr[rk];
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(&part_hbuf_h(m,s,0), nbuf, MPI_ATHENA_REAL, rk,
                4*lidn + 2*(1-s) + te, MPI_COMM_WORLD, &req.back());
    }
  }
  MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
  Kokkos::deep_copy(part_hbuf, part_hbuf_h);
  auto nbrk = part_nbr;
  par_for("m1_impl_halo_unpack", DevExeSpace(), 0, nmb-1, 0, nq-1, 0, nl-1,
          ks, ke, js, je,
  KOKKOS_LAMBDA(const int m, const int n, const int l, const int k, const int j) {
    int nc = ep1 ? M1_IW_EP : M1HaloCompA(n);
    int idx = (((n*nl + l)*nk + (k-ks))*nj + (j-js));
    if (nbrk.d_view(2*m) < 0) {iw_(m,nc,k,j,is-1-l) = hb_(m,2,idx);}
    if (nbrk.d_view(2*m+1) < 0) {iw_(m,nc,k,j,ie+1+l) = hb_(m,3,idx);}
  });
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::ImplicitGatherSolve
//! \brief gather the assembled rows (a,b,c,r) of every x1 stack onto its root block, run
//! there the IDENTICAL serial Thomas sweep the single-block solve runs, and scatter the
//! solution back into iw(M1_IW_S2).  Two messages per Picard iteration per non-root
//! block (a gather of 4*nx1 reals per column and a scatter of nx1).
//!
//! Cost: the root sweeps part_nblk*nx1 rows serially per column.  That is acceptable
//! here (the columns are independent and the kernel is parallel over (root,k,j), and
//! these columns are 84-512 cells long), but it is the serial bottleneck of the scheme;
//! the scalable successor is Schur condensation to the block-interface unknowns (one
//! reduced tridiagonal system of part_nblk rows per column, solved after two local
//! sweeps) or cyclic reduction, neither of which can be bitwise.

void RadiationM1::ImplicitGatherSolve() {
  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmy_pack->nmb_thispack;
  int nx1 = indcs.nx1, nx1g = part_nx1g;
  int nj = je - js + 1, nk = ke - ks + 1;
  auto iw_ = iw;
  auto sys_ = part_sys;
  auto pos_ = part_pos;
  auto slot_ = part_slot;

  // (1) same-rank members: copy their rows straight into the root's gathered system
  par_for("m1_impl_gth_loc", DevExeSpace(), 0, nmb-1, 0, 3, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
    int sl = slot_.d_view(m);
    if (sl < 0) return;
    sys_(sl,n,k,j,pos_.d_view(m)*(ie-is+1) + (i-is)) = iw_(m,M1_IW_TA+n,k,j,i);
  });

#if MPI_PARALLEL_ENABLED
  int nrow = 4*nk*nj*nx1;
  std::vector<MPI_Request> req;
  if (part_any_mpi) {
    auto sb_ = part_sbuf;
    par_for("m1_impl_gth_pack", DevExeSpace(), 0, nmb-1, 0, 3, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int n, const int k, const int j, const int i) {
      if (slot_.d_view(m) >= 0) return;
      sb_(m,((n*nk + (k-ks))*nj + (j-js))*(ie-is+1) + (i-is)) = iw_(m,M1_IW_TA+n,k,j,i);
    });
    Kokkos::deep_copy(part_sbuf_h, part_sbuf);
    int *gr = pmy_pack->pmesh->gids_eachrank;
    for (int s=0; s<part_nroot; ++s) {
      for (int p=0; p<part_nblk; ++p) {
        int rk = part_mrank[s*part_nblk+p];
        if (rk == global_variable::my_rank) continue;
        req.push_back(MPI_REQUEST_NULL);
        MPI_Irecv(&part_rbuf_h(s*part_nblk+p,0), nrow, MPI_ATHENA_REAL, rk,
                  2*(part_mgid[s*part_nblk+p] - gr[rk]), MPI_COMM_WORLD, &req.back());
      }
    }
    for (int m=0; m<nmb; ++m) {
      if (part_slot.h_view(m) >= 0) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Isend(&part_sbuf_h(m,0), nrow, MPI_ATHENA_REAL, part_rootrank[m],
                2*m, MPI_COMM_WORLD, &req.back());
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    req.clear();
    if (part_nroot > 0) {
      Kokkos::deep_copy(part_rbuf, part_rbuf_h);
      auto rb_ = part_rbuf;
      // which (slot,position) pairs are remote: encoded as a host loop over kernels
      for (int s=0; s<part_nroot; ++s) {
        for (int p=0; p<part_nblk; ++p) {
          if (part_mrank[s*part_nblk+p] == global_variable::my_rank) continue;
          const int ss = s, pp = p, row = s*part_nblk + p;
          par_for("m1_impl_gth_unp", DevExeSpace(), 0, 3, ks, ke, js, je, 0, nx1-1,
          KOKKOS_LAMBDA(const int n, const int k, const int j, const int i) {
            sys_(ss,n,k,j,pp*nx1 + i) =
                rb_(row,((n*nk + (k-ks))*nj + (j-js))*nx1 + i);
          });
        }
      }
    }
  }
#endif

  // (2) the gathered Thomas sweep: the same arithmetic, over nx1g rows
  if (part_nroot > 0) {
    par_for("m1_impl_gth_thomas", DevExeSpace(), 0, part_nroot-1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int s, const int k, const int j) {
      Real bet = sys_(s,1,k,j,0);
      sys_(s,5,k,j,0) = sys_(s,3,k,j,0)/bet;
      for (int i=1; i<nx1g; ++i) {
        sys_(s,4,k,j,i) = sys_(s,2,k,j,i-1)/bet;
        bet = sys_(s,1,k,j,i) - sys_(s,0,k,j,i)*sys_(s,4,k,j,i);
        sys_(s,5,k,j,i) = (sys_(s,3,k,j,i) - sys_(s,0,k,j,i)*sys_(s,5,k,j,i-1))/bet;
      }
      for (int i=nx1g-2; i>=0; --i) {
        sys_(s,5,k,j,i) -= sys_(s,4,k,j,i+1)*sys_(s,5,k,j,i+1);
      }
    });
  }

  // (3) scatter
#if MPI_PARALLEL_ENABLED
  if (part_any_mpi) {
    int nsol = nk*nj*nx1;
    if (part_nroot > 0) {
      auto rb_ = part_rbuf;
      for (int s=0; s<part_nroot; ++s) {
        for (int p=0; p<part_nblk; ++p) {
          if (part_mrank[s*part_nblk+p] == global_variable::my_rank) continue;
          const int ss = s, pp = p, row = s*part_nblk + p;
          par_for("m1_impl_sct_pack", DevExeSpace(), ks, ke, js, je, 0, nx1-1,
          KOKKOS_LAMBDA(const int k, const int j, const int i) {
            rb_(row,((k-ks)*nj + (j-js))*nx1 + i) = sys_(ss,5,k,j,pp*nx1 + i);
          });
        }
      }
      Kokkos::deep_copy(part_rbuf_h, part_rbuf);
    }
    int *gr = pmy_pack->pmesh->gids_eachrank;
    for (int m=0; m<nmb; ++m) {
      if (part_slot.h_view(m) >= 0) continue;
      req.push_back(MPI_REQUEST_NULL);
      MPI_Irecv(&part_sbuf_h(m,0), nsol, MPI_ATHENA_REAL, part_rootrank[m],
                2*m + 1, MPI_COMM_WORLD, &req.back());
    }
    for (int s=0; s<part_nroot; ++s) {
      for (int p=0; p<part_nblk; ++p) {
        int rk = part_mrank[s*part_nblk+p];
        if (rk == global_variable::my_rank) continue;
        req.push_back(MPI_REQUEST_NULL);
        MPI_Isend(&part_rbuf_h(s*part_nblk+p,0), nsol, MPI_ATHENA_REAL, rk,
                  2*(part_mgid[s*part_nblk+p] - gr[rk]) + 1, MPI_COMM_WORLD,
                  &req.back());
      }
    }
    MPI_Waitall(static_cast<int>(req.size()), req.data(), MPI_STATUSES_IGNORE);
    Kokkos::deep_copy(part_sbuf, part_sbuf_h);
    auto sb_ = part_sbuf;
    par_for("m1_impl_sct_unp", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (slot_.d_view(m) >= 0) return;
      iw_(m,M1_IW_S2,k,j,i) = sb_(m,((k-ks)*nj + (j-js))*(ie-is+1) + (i-is));
    });
  }
#endif
  par_for("m1_impl_sct_loc", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    int sl = slot_.d_view(m);
    if (sl < 0) return;
    iw_(m,M1_IW_S2,k,j,i) = sys_(sl,5,k,j,pos_.d_view(m)*(ie-is+1) + (i-is));
  });
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
  auto ifw_ = ifw;
  const bool aphll = (impl_flux != M1_IFLUX_CENTRAL);
  const bool berth = (impl_flux == M1_IFLUX_BERTHON) || (impl_flux == M1_IFLUX_BLEND);
  // MILESTONE 3c: the smooth per-face convex blend.  `bkind` etc. are read only when
  // `blend` is true, and `blend` is false for every 3a/3a2 flux, so those paths keep the
  // arithmetic they had (the weight enters as an exact multiplication by 1.0).
  const bool blend = (impl_flux == M1_IFLUX_BLEND);
  const int bkind = impl_blend, bfm = impl_blend_fmode;
  const bool bdis = blend && (impl_blend_mode == M1_IBMODE_DISSIP);
  const Real btau0 = impl_blend_tau0;
  const Real bflo = impl_blend_flo, bfhi = impl_blend_fhi;
  const bool plmdc = (impl_recon == M1_IRECON_PLMDC);
  const Real rwin = impl_recon_w;
  const Real rfl_ = impl_res_floor;
  const bool rfreeze = impl_recon_freeze;
  const int rnpass = impl_recon_npass;
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
  // LIMIT 3 of the 3a findings.  A physical boundary face hands its whole
  // dt (rho k_t)_f F0_f/c to its ONE interior cell in 3a, so that cell receives 1.5
  // face-shares of radiative force where every interior cell receives 1.0; the residual
  // is a steady force the well-balanced reference a_rad_ref does not carry, and it drives
  // the 10.3 v_MLT bottom-cell flow of I8.  With implicit_bmom_half the boundary face
  // gives HALF, like any other face: the cell-averaged radiative force is then
  // (rho k_t F/c) with F the mean of the cell's two faces, everywhere.  The other half
  // leaves the domain with the radiation, which is where it physically goes.
  const bool bmhalf = impl_bmom_half;
  bool fref = (force_ref == M1_FREF_WB_ARAD);
  auto aref_ = arad_ref;
  Real mq = marshak_q;
  int bclo = ibc_x1min, bchi = ibc_x1max;
  Real fxlo = iflux_x1min, fxhi = iflux_x1max;
  Real eblo = iebath_x1min, ebhi = iebath_x1max;
  bool cyclic = (bclo == M1_IBC_PERIODIC);
  // MILESTONE 3b, LIMIT 4.  With more than one MeshBlock along x1 a block is at a
  // PHYSICAL x1 boundary only when it sits at the corresponding end of its stack; the
  // faces it shares with a stack neighbour are ordinary interior faces whose other cell
  // is a GHOST cell, filled by ImplicitX1Halo with the very numbers the neighbour
  // computed.  With part_nblk == 1 every block is both ends and nothing below changes.
  const int nblkx1 = part_nblk;
  auto pos_ = part_pos;
  const int nlay_ = (part_nblk > 1) ? part_nlay : 0;

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
    iw_(m,M1_IW_KT,k,j,i) = opac_(m,M1_OP_T,k,j,i);
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

  // The SCALE of the Picard convergence test.  With implicit_res_floor = 0 (the 3a
  // default) the residual is the pure relative change |dE|/E, which in a run with a large
  // dynamic range is dominated by cells many orders below the peak: on gate I6 with
  // implicit_recon = plm_dc the tail cells sit 9 orders under the maximum and keep the
  // reported residual above the tolerance for ever, although the SOLUTION is converged
  // (maxit 30 and maxit 100 give the same amplitude to five digits).  A positive
  // implicit_res_floor scales those cells by the column peak instead,
  // res = |dE|/max(E, implicit_res_floor*max(E)).
  Real emax0 = 0.0;
  if (rfl_ > 0.0) {
    Kokkos::parallel_reduce("m1_impl_emax",
    Kokkos::MDRangePolicy<Kokkos::Rank<4>>(DevExeSpace(), {0,ks,js,is},
                                           {nmb1+1,ke+1,je+1,ie+1}),
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i, Real &lmax) {
      Real r = iw_(m,M1_IW_EN,k,j,i);
      lmax = (r > lmax) ? r : lmax;
    }, Kokkos::Max<Real>(emax0));
#if MPI_PARALLEL_ENABLED
    {Real g;
    MPI_Allreduce(&emax0, &g, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    emax0 = g;}
#endif
  }
  const Real escale = rfl_*emax0;

  // the x1 ghost layers the partitioned solve reads.  E of the start-of-step state is
  // sent once here; the six lagged quantities and the new E are sent inside the loop.
  ImplicitX1Halo(true);

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
        iw_(m,M1_IW_KT,k,j,i) = d*(of + os);
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
      // The COMOVING reduced flux of the iterate, which is what the HLL part of the
      // ap_hll flux lags (the HLL acts on F0 only; the enthalpy flux A is added back
      // upwinded, exactly as in the explicit advective split).  The LAB f above still
      // drives the closure chi, as it does in the explicit scheme.
      //
      // It is NOT F0_cell/(c E) with F0_cell the arithmetic mean of the two faces.  That
      // is design risk R4 and it is fatal here: in free streaming the upwind face flux is
      // c E_{i-1}, so the cell mean is c (E_{i-1}+E_i)/2 and the derived f is
      // (1 + E_{i-1}/E_i)/2, i.e. 0.5 rather than 1 on the steep side of a pulse.  The
      // wave speeds then reopen to +-c/sqrt(3), the HLL flux turns CENTRED, and the
      // I6 pulse is flattened to its box mean in one crossing (amplitude ratio 0.0014).
      // Each FACE flux is therefore normalised by the E of the cell it comes FROM, which
      // is exactly 1 for an upwind free-streaming face, and the cell value is the mean of
      // the two face ratios.  On a cold start (f0x1 is zero-initialised and the problem
      // generator's state lives in u0) the cell flux is used instead.
      Real r0;
      Real fl = f0_(m,k,j,i), fr = f0_(m,k,j,i+1);
      if (fabs(fl) + fabs(fr) > 0.0) {
        int iml = (i > is) ? (i-1)
                  : (cyclic ? ie : ((pos_.d_view(m) > 0) ? (is-1) : is));
        int ipr = (i < ie) ? (i+1)
                  : (cyclic ? is : ((pos_.d_view(m) < nblkx1-1) ? (ie+1) : ie));
        Real eul = fmax((fl > 0.0) ? iw_(m,M1_IW_EP,k,j,iml) : e, efl);
        Real eur = fmax((fr > 0.0) ? e : iw_(m,M1_IW_EP,k,j,ipr), efl);
        r0 = 0.5*(fl/(cl*eul) + fr/(cl*eur));
      } else {
        r0 = (f1 - iw_(m,M1_IW_ADV,k,j,i)*e)/(cl*e);
      }
      if (r0 > 1.0) {r0 = 1.0;}
      if (r0 < -1.0) {r0 = -1.0;}
      iw_(m,M1_IW_RF0,k,j,i) = r0;
    });

    // the x1 halo of the LAGGED quantities (w, a, g0, v1, the comoving reduced flux and
    // the transport opacity).  Every face of the stack is then assembled by both of its
    // blocks from bit-identical numbers.
    ImplicitX1Halo(false);

    // (b2) the FACE coefficients of the asymptotic-preserving HLL blend.  Nothing here
    // runs under implicit_flux = central, where ifw stays identically zero and the row
    // assembled below is bitwise the 3a one.
    // The deferred correction is by default evaluated ONCE per step, at the start-of-step
    // state (implicit_recon_lag = step).  Recomputing it every Picard pass
    // (= picard) makes the loop a limit cycle: the plm limiter keeps switching on a few
    // cells and the strict tolerance is never reached, at 30 iterations per step against
    // 2, although the answer is the same to five digits.  The correction is a lagged,
    // explicit term in any case, so evaluating it at E^n costs nothing in order.
    const int iter = it;
    const bool doface = aphll && (it == 0 || !rfreeze);
    if (doface) {
      const bool dodg = plmdc && (it == 0
                                  || (!rfreeze && (rnpass <= 0 || it < rnpass)));
      par_for("m1_impl_aphll", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        int ipos = pos_.d_view(m);
        bool phys = (((i == is) && (ipos == 0)) ||
                     ((i == ie+1) && (ipos == nblkx1-1))) && !cyclic;
        if (phys) {
          // a physical boundary face: the flux is IMPOSED there (flux / Marshak /
          // reflect / efix), so there is no Riemann problem and no blend.
          ifw_(m,M1_IFW_AL,k,j,i) = 0.0;
          ifw_(m,M1_IFW_HCL,k,j,i) = 0.0;
          ifw_(m,M1_IFW_HCR,k,j,i) = 0.0;
          ifw_(m,M1_IFW_DG,k,j,i) = 0.0;
          return;
        }
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        Real dx = mbsize.d_view(m).dx1;
        Real rfl = iw_(m,M1_IW_RF0,k,j,im);
        Real rfr = iw_(m,M1_IW_RF0,k,j,ip);
        // closed-form M1 wave speeds of the two LAGGED states (1-D: mu = sign f)
        Real bl, br;
        if (edd) {
          br = ch/sqrt(3.0);
          bl = -br;
        } else {
          Real lml, lpl, lmr, lpr;
          M1WaveSpeeds(fabs(rfl), (rfl >= 0.0) ? 1.0 : -1.0, lml, lpl);
          M1WaveSpeeds(fabs(rfr), (rfr >= 0.0) ? 1.0 : -1.0, lmr, lpr);
          bl = ch*fmin(fmin(lml, lmr), 0.0);
          br = ch*fmax(fmax(lpl, lpr), 0.0);
        }
        // alpha: Bloch et al. (2021) eq. 25 with the (1-f^2) guard and the arithmetic
        // face mean of the CELL optical depth.  lp*lm <= 0, so den >= 1 and alpha <= 1.
        Real tauf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,ip))*dx;
        Real al = 1.0;
        if (tauf > 0.0) {
          Real fbar = 0.5*(fabs(rfl) + fabs(rfr));
          Real guard = fmax(1.0 - fbar*fbar, 0.0);
          Real lp = br/ch, lm = bl/ch;
          Real den = 1.0 - 3.0*tauf*guard*lp*lm/(lp - lm + 1.0e-300);
          al = 1.0/fmax(den, 1.0);
        }
        // F_HLL = [b_R c_h f_L E'_L - b_L c_h f_R E'_R + b_R b_L (E'_R - E'_L)]/(b_R-b_L)
        // is linear in E'.  Split into the ADVECTIVE part (the two physical fluxes) and
        // the DISSIPATION (the jump term), because the two carry different weights:
        //
        //   F = alpha F_adv + alpha^2 F_dis + (1 - alpha) F_diff.
        //
        // The dissipation must carry alpha^2 and not alpha.  At piecewise-constant states
        // -- which is what the MATRIX is built from, whatever implicit_recon says -- the
        // HLL dissipation IS the physical diffusion (b_R b_L dE/(b_R-b_L) -> -c dE/3),
        // so weighting it alpha and adding (1-alpha) F_diff on top counts the diffusion
        // TWICE: measured on gate I1 at tau_cell = 1e3, d(sigma^2)/dt came out 2.0022 x
        // the analytic 2D at every CFL.  alpha^2 ~ 1/tau^2 kills it against F_diff
        // ~ 1/tau and leaves the thin limit (alpha -> 1) exactly the plain HLL flux.
        // This is the `alpha2` form of the explicit scheme (rad_m1_closure.hpp), reached
        // here for the same reason.
        //
        // The two fmax()/fmin() are the M-MATRIX GUARDS.  At alpha = 1 they are provably
        // inactive (the HLL consistency condition b_L <= c_h f <= b_R holds on the M1
        // admissible set), but alpha < 1 rescales the two parts differently and the
        // E'_R coefficient can turn positive; the clamp then drops it to zero, which
        // only makes the face flux more upwind and leaves conservation exact (it is one
        // number per face, used with opposite signs by the two cells).
        Real invb = 1.0/(br - bl + 1.0e-300);
        Real adl = br*ch*rfl*invb;        // E'_L coefficient of F_adv
        Real adr = -bl*ch*rfr*invb;       // E'_R coefficient of F_adv
        Real dk = -br*bl*invb;            // >= 0, the dissipation coefficient
        //
        // implicit_flux = berthon drops F_diff altogether and weights BOTH parts of the
        // HLL flux by alpha, which is what alpha was constructed for: alpha (F_adv +
        // F_dis) is the physical diffusion to first order in 1/tau, the advective part
        // supplying the 1/(0.866 tau) that the dissipation alone is short of.  The
        // F_diff weight is then zero, which is what storing AL = 1 below means.
        Real wdis = berth ? al : (al*al);
        // MILESTONE 3c.  w_f in [0,1] from the LAGGED face quantities: the whole face
        // flux is (1 - w_f) F_central + w_f F_berthon (implicit_blend_mode = flux), or
        // the FULL central flux plus w_f times the HLL DISSIPATION alone
        // (= dissipation).  w_f = 1 with mode = flux reproduces `berthon` bitwise and
        // w_f = 0 reproduces `central` bitwise: the multiplications below are by exactly
        // 1.0 or exactly 0.0.  Every input of the weight lives in iw, which the x1 halo
        // of the partitioned solve already carries.
        Real wf = 1.0;
        if (blend) {
          wf = M1BlendWeight(bkind, bfm, tauf, btau0, rfl, rfr, bflo, bfhi);
        }
        // In `dissipation` mode the advective part of the HLL flux is NOT added (the
        // central flux already carries the transport); only the jump term is.
        Real wadv = bdis ? 0.0 : al;
        Real wdsq = bdis ? al : wdis;
        Real ccl = wf*fmax(wadv*adl + wdsq*dk, 0.0);
        Real ccr = wf*fmin(wadv*adr - wdsq*dk, 0.0);
        // how much of the CENTRAL (face-eliminated) flux the row keeps is 1 - AL.
        ifw_(m,M1_IFW_AL,k,j,i) = bdis ? 0.0 : (berth ? wf : al);
        ifw_(m,M1_IFW_HCL,k,j,i) = ccl;
        ifw_(m,M1_IFW_HCR,k,j,i) = ccr;
        // the plm DEFERRED CORRECTION: the difference between the plm and the dc HLL
        // flux at the PREVIOUS iterate.  It goes to the right-hand side, so the matrix
        // stays the low-order M-matrix.  Both E and the comoving reduced flux are
        // reconstructed, with the same limiter the explicit scheme uses; a face whose
        // 4-cell stencil leaves the block falls back to dc (zero correction).
        Real dg = 0.0;
        if (dodg && al > 0.0) {
          int ilo = is - ((ipos > 0) ? nlay_ : 0);
          int ihi = ie + ((ipos < nblkx1-1) ? nlay_ : 0);
          int imm = (im > ilo) ? (im-1) : (cyclic ? ie : -1);
          int ipp = (ip < ihi) ? (ip+1) : (cyclic ? is : -1);
          if (imm >= 0 && ipp >= 0) {
            Real dum;
            Real elp, erp, flp, frp;
            PLM(iw_(m,M1_IW_EP,k,j,imm), iw_(m,M1_IW_EP,k,j,im),
                iw_(m,M1_IW_EP,k,j,ip), elp, dum);
            PLM(iw_(m,M1_IW_EP,k,j,im), iw_(m,M1_IW_EP,k,j,ip),
                iw_(m,M1_IW_EP,k,j,ipp), dum, erp);
            PLM(iw_(m,M1_IW_RF0,k,j,imm), iw_(m,M1_IW_RF0,k,j,im),
                iw_(m,M1_IW_RF0,k,j,ip), flp, dum);
            PLM(iw_(m,M1_IW_RF0,k,j,im), iw_(m,M1_IW_RF0,k,j,ip),
                iw_(m,M1_IW_RF0,k,j,ipp), dum, frp);
            Real ecl = iw_(m,M1_IW_EP,k,j,im), ecr = iw_(m,M1_IW_EP,k,j,ip);
            // the deferred correction applies to the UPWIND part only, so it carries the
            // same weight w_f the upwind part carries (3c).
            Real gp = wf*(wadv*(br*ch*flp*elp - bl*ch*frp*erp)*invb
                          + wdsq*(-dk)*(erp - elp));
            Real gc = ccl*ecl + ccr*ecr;
            // ADMISSIBILITY of the corrected face flux against the DONOR cell.  The
            // reconstructed face energy may exceed the donor cell's own (plm puts
            // E_i (r-1)/(r+1) on top of E_i for a geometric ratio r), and c times that is
            // then faster than the donor can physically emit: the cell drains below what
            // it receives and, on an exponentially falling background, the drain
            // cascades.  Measured on I6 before this clamp: the 1e-4 background of the
            // free-streaming pulse collapsed onto the floor and the peak grew 7x
            // (amplitude ratio 6.98, Picard never converging).  The low-order flux gc
            // already satisfies this bound, so the clamp never removes the whole
            // correction, only the inadmissible part of it.
            Real gmax = wf*ch*ecl, gmin = -wf*ch*ecr;
            gp = fmin(fmax(gp, gmin), gmax);
            // The DEFERRED-CORRECTION WEIGHT.  A deferred correction is a fixed-point
            // iteration x <- A_low^-1 (b + (A_low - A_high) x), and for advection its
            // contraction factor is ~ 2 nu/(1 + nu) with nu = chat dt/dx: it converges
            // only below CFL ~ 1 and DIVERGES above it (measured: Picard never converges
            // at implicit_cfl = 10 and the I6 pulse amplitude comes out 3.85).  The
            // correction is therefore weighted by w = 1/(1 + nu) unless
            // <rad_m1>/implicit_recon_w names a fixed value.  That makes the contraction
            // factor 2 nu/(1 + nu)^2 <= 1/2 at EVERY CFL, and the fixed point a convex
            // blend of the dc and plm fluxes -- still a monotone flux, second-order where
            // w -> 1 (nu << 1, which is where a propagating front is resolved in time at
            // all) and dc where the step is so long that the front is not resolved.
            Real wdc = (rwin > 0.0) ? rwin : (1.0/(1.0 + ch*dt/dx));
            dg = wdc*(gp - gc);
            // UNDER-RELAXATION across the Picard passes.  The plm limiter keeps switching
            // on a handful of cells and the un-relaxed iteration is a small-amplitude
            // limit cycle that never meets the tolerance (30 passes per step against 2,
            // with the answer already right to five digits).  Averaging with the previous
            // pass leaves the fixed point untouched and breaks the cycle.
            if (iter > 0) {dg = 0.5*(dg + ifw_(m,M1_IFW_DG,k,j,i));}
          }
        }
        ifw_(m,M1_IFW_DG,k,j,i) = dg;
      });
    }

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
      int ipos = pos_.d_view(m);
      bool botb = (ipos == 0), topb = (ipos == nblkx1-1);
      Real wi = iw_(m,M1_IW_WCHI,k,j,i);
      Real ai = iw_(m,M1_IW_ADV,k,j,i);
      Real vi = iw_(m,M1_IW_V1,k,j,i);
      Real aa = 0.0, bb = 1.0, cc = 0.0;
      Real rr = iw_(m,M1_IW_EN,k,j,i) + iw_(m,M1_IW_SRCR,k,j,i);
      bb += iw_(m,M1_IW_SRCB,k,j,i);

      // ---- face i+1/2
      if (i < ie || cyclic || !topb) {
        int ip = (i < ie) ? (i+1) : (cyclic ? is : (ie+1));
        Real om = 1.0 - ifw_(m,M1_IFW_AL,k,j,i+1);
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = om*th*ch*ch*dt/dx;
        bb += nu*df*wi;
        cc -= nu*df*iw_(m,M1_IW_WCHI,k,j,ip);
        Real vf = 0.5*(vi + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,i) + iw_(m,M1_IW_G0,k,j,ip));
        rr -= nu*cr*om*th*(f0n_(m,k,j,i+1) - ch*dt*vf*g0f);
        // the HLL part: its E'_L coefficient is >= 0 (diagonal) and its E'_R coefficient
        // <= 0 (upper off-diagonal), so the blend keeps the M-matrix.
        bb += nu*ifw_(m,M1_IFW_HCL,k,j,i+1);
        cc += nu*ifw_(m,M1_IFW_HCR,k,j,i+1);
        rr -= nu*ifw_(m,M1_IFW_DG,k,j,i+1);
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
      if (i > is || cyclic || !botb) {
        int im = (i > is) ? (i-1) : (cyclic ? ie : (is-1));
        Real om = 1.0 - ifw_(m,M1_IFW_AL,k,j,i);
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real df = om*th*ch*ch*dt/dx;
        bb += nu*df*wi;
        aa -= nu*df*iw_(m,M1_IW_WCHI,k,j,im);
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + vi);
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,i));
        rr += nu*cr*om*th*(f0n_(m,k,j,i) - ch*dt*vf*g0f);
        aa -= nu*ifw_(m,M1_IFW_HCL,k,j,i);
        bb -= nu*ifw_(m,M1_IFW_HCR,k,j,i);
        rr += nu*ifw_(m,M1_IFW_DG,k,j,i);
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
      if (!cyclic && ((i == is && botb && bclo == M1_IBC_EFIX) ||
                      (i == ie && topb && bchi == M1_IBC_EFIX))) {
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

    // (e) one Thomas sweep per column (cyclic: Sherman-Morrison).  When the column
    // spans several MeshBlocks the rows are gathered onto its root block and swept there
    // by the identical serial recurrence (ImplicitGatherSolve).
    if (nblkx1 > 1) {
      ImplicitGatherSolve();
    } else {
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
    }

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
        Real re = fabs(enew - eold)/fmax(fmax(fabs(enew), escale), 1.0e-300);
        Real rt = fabs(tnew - told)/fmax(fabs(tnew), 1.0e-300);
        iw_(m,M1_IW_RES,k,j,i) = fmax(re, rt);
      });
    } else {
      par_for("m1_impl_accept", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        Real enew = fmax(iw_(m,M1_IW_S2,k,j,i), efl);
        Real eold = iw_(m,M1_IW_EP,k,j,i);
        iw_(m,M1_IW_EP,k,j,i) = enew;
        iw_(m,M1_IW_RES,k,j,i) = fabs(enew - eold)
                                 /fmax(fmax(fabs(enew), escale), 1.0e-300);
      });
    }

    // the new iterate's E has to reach the ghost cells before the FACE update, or the
    // two blocks that share a face would build it from different states.
    ImplicitX1Halo(true);

    // (g) the face fluxes of the new iterate, and the derived cell flux
    par_for("m1_impl_face", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      Real dx = mbsize.d_view(m).dx1;
      int ipos = pos_.d_view(m);
      bool lo = (i == is) && (ipos == 0);
      bool hi = (i == ie+1) && (ipos == nblkx1-1);
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
        int im = (cyclic && i == is) ? ie : (i-1);
        int ip = (cyclic && i == ie+1) ? is : i;
        Real ktf = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,ip));
        Real th = 1.0/(1.0 + ch*dt*ktf);
        Real vf = 0.5*(iw_(m,M1_IW_V1,k,j,im) + iw_(m,M1_IW_V1,k,j,ip));
        Real g0f = 0.5*(iw_(m,M1_IW_G0,k,j,im) + iw_(m,M1_IW_G0,k,j,ip));
        Real gr = (iw_(m,M1_IW_WCHI,k,j,ip)*iw_(m,M1_IW_EP,k,j,ip)
                   - iw_(m,M1_IW_WCHI,k,j,im)*iw_(m,M1_IW_EP,k,j,im))/dx;
        Real fn = th*(f0n_(m,k,j,(i == ie+1 && cyclic) ? is : i)
                      - ch*cl*dt*gr - ch*dt*vf*g0f);
        // the SAME blend the row was assembled with, so that the stored comoving face
        // flux (which the restart file carries, which the momentum deposit uses and which
        // the next iterate's reduced flux is built from) is the flux the solve applied.
        // The F0^n MEMORY term sits entirely in the F_diff branch: in the thin limit
        // alpha -> 1 and it must NOT survive, or the lagged flux would fight the upwind
        // HLL flux and the front would be damped exactly as it is under `central`.
        // (3c) the test is on the FLUX FORM, not on al: with implicit_blend_mode =
        // dissipation the central weight is 1 (al = 0) and the HLL part is still there.
        Real al = ifw_(m,M1_IFW_AL,k,j,i);
        if (aphll) {
          Real gh = ifw_(m,M1_IFW_HCL,k,j,i)*iw_(m,M1_IW_EP,k,j,im)
                    + ifw_(m,M1_IFW_HCR,k,j,i)*iw_(m,M1_IW_EP,k,j,ip)
                    + ifw_(m,M1_IFW_DG,k,j,i);
          fn = (1.0 - al)*fn + (cl/ch)*gh;
        }
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
#if MPI_PARALLEL_ENABLED
    // the convergence test must be GLOBAL: with a partitioned column the ranks would
    // otherwise take different numbers of Picard passes and the gather would deadlock,
    // and even with rank-local columns a per-rank test makes the answer depend on the
    // decomposition.  One MPI_MAX of one double per pass.
    {Real rg;
    MPI_Allreduce(&resid, &rg, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    resid = rg;}
#endif
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
        int ipos = pos_.d_view(m);
        if (i == is && ipos == 0 && !cyclic) {
          ktl = iw_(m,M1_IW_KT,k,j,i);
          wl = bmhalf ? 0.5 : 1.0;
        } else {
          int im = (cyclic && i == is) ? ie : (i-1);
          ktl = 0.5*(iw_(m,M1_IW_KT,k,j,im) + iw_(m,M1_IW_KT,k,j,i));
        }
        if (i == ie && ipos == nblkx1-1 && !cyclic) {
          ktr = iw_(m,M1_IW_KT,k,j,i);
          wr = bmhalf ? 0.5 : 1.0;
        } else {
          int ip = (cyclic && i == ie) ? is : (i+1);
          ktr = 0.5*(iw_(m,M1_IW_KT,k,j,i) + iw_(m,M1_IW_KT,k,j,ip));
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
