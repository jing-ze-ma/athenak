//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_wball.cpp
//! \brief wball-1104: <mhd>/wb_dirs = all -- MHD's members of the shared all-direction
//! well-balanced routines (hydro/hydro_wball.hpp/.cpp).  The background is the GAS
//! pressure only (the field is part of the deviation, as in the x1 scheme); the AM
//! correction of the background fluxes is taken field-free.

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro_wball.hpp"
#include "mhd.hpp"

namespace mhd {

hydro::WbAllCtx MHD::WbCtx() {
  hydro::WbAllCtx c;
  c.pmbp = pmy_pack; c.eos = peos->eos_data; c.wbo = wb_option;
  c.wgt = wb_weight; c.wbgi = wb_guard_ideal; c.wbs = use_wellbalance_static;
  c.r0 = wbw_r0; c.r1 = wbw_r1; c.m0 = wbw_m0; c.m1 = wbw_m1; c.eta = wbw_eta;
  c.w0 = w0; c.u0 = u0; c.wder = wder; c.wbq0 = wbq0; c.wbq2 = wbq2; c.wbq3 = wbq3;
  c.wbw = wbw; c.phicc = phicc0; c.wtemp = wtemp; c.pwb = pwb; c.phi1f = phi0.x1f; c.phi2f = phi0.x2f; c.phi3f = phi0.x3f; c.pam = pam;
  c.wbam_u = &wbam_u; c.wbam_w = &wbam_w; c.wbam_f = &wbam_f;
  c.phimax = wb_phimax; c.phiwb = phicc_wb;
  return c;
}
void MHD::BuildWBCacheAll() {hydro::WbAllCtx c = WbCtx(); hydro::WbAllBuild(c);}
void MHD::AddWbGravity(const Real bdt, const bool with_x1) {
  hydro::WbAllCtx c = WbCtx(); hydro::WbAllGravity(c, bdt, with_x1);
}

}  // namespace mhd
