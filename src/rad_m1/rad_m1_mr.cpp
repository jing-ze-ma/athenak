//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_mr.cpp
//! \brief <rad_m1>/implicit_mr_every = k: MULTI-RATE implicit radiation (m1-fast4,
//! tests_m1/runs_5r_fast4).
//!
//! The implicit radiation (transport + gas coupling) is unconditionally stable, and where
//! its state follows the hydro quasi-statically nothing requires it to be solved on
//! every hydro step.  With k > 1 the hydro (Heun) takes every step ALONE, and the
//! radiation R advances the gas + radiation over a whole window Delta of k hydro steps,
//! operator split:
//!   H^{k/2} R(Delta) H^k R(Delta) H^k ...   =   (H(Delta/2) R(Delta) H(Delta/2))^N,
//! i.e. R sits at the window CENTRES (the first one after k/2 steps, with Delta twice
//! the elapsed time), which is Strang splitting: second order in time when R is.  At a
//! time that is a window END (t = N Delta from the start) the state is the Strang one;
//! in between, it is off by a non-accumulating half-window splitting term.
//!
//! R(Delta) is the two-stage, L-stable, stiffly accurate SDIRK2 (g = 1 - 1/sqrt 2):
//!   Y_A = Y_0 + g Delta K_A,  Y_B = Y_0 + (1-g) Delta K_A + g Delta K_B,  Y_1 = Y_B,
//! built from the hesdirk2 stage machinery (rad_m1_time2.cpp): stage A is a STAGE1
//! solve (dt = g Delta, t2inc = 0, K_A stored in t2k2), stage B a STAGE2 solve from
//! Y_A with t2inc = (1 - 2g) Delta K_A.  The closure tensor is built once per window at
//! Y_0, as the hesdirk2 step builds it at U^n.  A stage that fails puts Y_0 back and
//! takes R(Delta) as backward Euler (counted).
//!
//! Needs time_scheme = hesdirk2 (its stage arrays).  Default k = 1: nothing here runs,
//! and every other path is bitwise unchanged.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "hydro/hydro.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/m1_fluid.hpp"

namespace radm1 {

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MRInit
//! \brief read <rad_m1>/implicit_mr_every (read only when named; 1 = off)

void RadiationM1::MRInit(ParameterInput *pin) {
  // vet_sc_every (read only when named): see VscSkip below
  // unify-sc-gd-1012 (#13): the canonical key vet_every, vet_sc_every its alias
  vsc_every = 1;
  const std::string kev = M1KeyName(pin, "vet_every", {"vet_sc_every"});
  if (vet_sc && pin->DoesParameterExist("rad_m1", kev)) {
    vsc_every = std::max(1, pin->GetInteger("rad_m1", kev));
    if (vsc_every > 1) {
      if (time_scheme != M1_TIME_HESDIRK2) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "<rad_m1>/vet_sc_every > 1 needs time_scheme = hesdirk2"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      const int nmb = pmy_pack->nmb_thispack;
      Kokkos::realloc(vsc_d0, nmb, M1_T2_NVET, static_cast<int>(vet_prev.extent(2)),
                      static_cast<int>(vet_prev.extent(3)),
                      static_cast<int>(vet_prev.extent(4)));
      Kokkos::realloc(vsc_d1, nmb, M1_T2_NVET, static_cast<int>(vet_prev.extent(2)),
                      static_cast<int>(vet_prev.extent(3)),
                      static_cast<int>(vet_prev.extent(4)));
      if (global_variable::my_rank == 0) {
        std::cout << "<rad_m1> vet_sc_every = " << vsc_every << ": formal solution on "
                  << "every " << vsc_every << "-th tensor build, linear extrapolation "
                  << "in time in between" << std::endl;
      }
    }
  }
  mr_every = 1;
  if (pin->DoesParameterExist("rad_m1", "dbg_cell_lo")) {
    dbg_cell_lo = pin->GetInteger("rad_m1", "dbg_cell_lo");
    dbg_cell_hi = pin->GetInteger("rad_m1", "dbg_cell_hi");
  }
  if (pin->DoesParameterExist("rad_m1", "dbg_hydro_skip_t2")) {   // raycheap2-1109
    dbg_hydro_skip_t2 = pin->GetBoolean("rad_m1", "dbg_hydro_skip_t2");
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadiationM1::MRReport
//! \brief one line at the end of the run (rank 0; called from ImplicitReport)

void RadiationM1::MRReport() {
  if (vsc_every > 1) {
    std::cout << "<rad_m1> vet_sc_every = " << vsc_every << ": builds=" << vsc_cnt
              << " extrapolated=" << vsc_nskip << " clipped cells=" << vsc_nclip
              << std::endl;
  }
  if (mr_every <= 1) return;
  std::cout << "<rad_m1> implicit_mr_every = " << mr_every << ": radiation steps="
            << mr_nr << " mean window=" << ((mr_nr > 0.0) ? mr_ksum/mr_nr : 0.0)
            << " backward-Euler fallbacks=" << mr_nfall;
  if (mr_theta > 0.0) {
    std::cout << " theta max=" << mr_thmax << " (guard " << mr_theta << ")";
  }
  std::cout << std::endl;
}

} // namespace radm1
