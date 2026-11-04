//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gresho_mhd.cpp
//! \brief Problem generator for a magnetized Gresho vortex at low Mach number.
//!
//! The Gresho vortex is a time-independent solution of the Euler equations: the
//! centrifugal force of a rotating ring of gas is balanced by a radial pressure
//! gradient.  Adding a UNIFORM field leaves the equilibrium untouched (a uniform field
//! exerts no force), so the exact solution is stationary in MHD as well, and any decay
//! of the kinetic energy is numerical dissipation.  Lowering the Mach number (i.e.
//! raising the background pressure at fixed peak rotation speed) makes this a standard
//! test of how a Riemann solver behaves in the low Mach number limit: the dissipation of
//! a standard Godunov solver grows like 1/M, and the vortex is destroyed in a turnover.
//!
//! The azimuthal velocity and pressure are (Gresho & Chan 1990; Liska & Wendroff 2003)
//!   u_phi = 5 r,                 p = p0 + 12.5 r^2                      (r < 0.2)
//!   u_phi = 2 - 5 r,             p = p0 + 12.5 r^2 + 4 - 20 r + 4 ln(5r)(0.2 < r < 0.4)
//!   u_phi = 0,                   p = p0 - 2 + 4 ln 2                    (r > 0.4)
//! with rho = d0, so that max(u_phi) = 1 at r = 0.2 and the core rotates rigidly with
//! angular frequency 5 (one turnover in t = 2 pi/5 = 1.2566).  The background pressure
//! follows from the requested Mach number, p0 = d0/(gamma M^2), and the uniform field
//! along x1 from the requested plasma beta, B0 = sqrt(2 p0/beta).
//!
//! REFERENCES:
//! - P. M. Gresho & S. T. Chan, "On the theory of semi-implicit projection methods...",
//!   Int. J. Numer. Methods Fluids, 11, 621 (1990)
//! - T. Minoshima & T. Miyoshi, "A low-dissipation HLLD approximate Riemann solver for a
//!   very wide range of Mach numbers", JCP, 446, 110639 (2021)

// C++ headers
#include <math.h>
#include <cstdint>
#include <iostream>   // endl

// Athena++ headers
#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "pgen/pgen.hpp"
#include "outputs/outputs.hpp"

//----------------------------------------------------------------------------------------
//! \fn Real GreshoVphi(const Real r)
//! \brief azimuthal velocity of the Gresho vortex

KOKKOS_INLINE_FUNCTION
Real GreshoVphi(const Real r) {
  if (r < 0.2) {
    return 5.0*r;
  } else if (r < 0.4) {
    return 2.0 - 5.0*r;
  }
  return 0.0;
}

//----------------------------------------------------------------------------------------
//! \fn Real GreshoPres(const Real r, const Real p0)
//! \brief pressure of the Gresho vortex, in hydrostatic balance with GreshoVphi

KOKKOS_INLINE_FUNCTION
Real GreshoPres(const Real r, const Real p0) {
  if (r < 0.2) {
    return p0 + 12.5*r*r;
  } else if (r < 0.4) {
    return p0 + 12.5*r*r + 4.0 - 20.0*r + 4.0*std::log(5.0*r);
  }
  return p0 - 2.0 + 4.0*std::log(2.0);
}

//----------------------------------------------------------------------------------------
//! \fn Real GreshoVint(const Real r)
//! \brief alfmach-1102: integral of GreshoVphi from 0 to r (flux function of b_phi)

KOKKOS_INLINE_FUNCTION
Real GreshoVint(const Real r) {
  if (r < 0.2) {
    return 2.5*r*r;
  } else if (r < 0.4) {
    return 0.1 + 2.0*(r - 0.2) - 2.5*(r*r - 0.04);
  }
  return 0.2;
}

namespace {
// lowmachpu-1103: state for the <problem>/user_hist = true history (GreshoPUHistory)
Real gpu_p0_ = 1.0, gpu_d0_ = 1.0, gpu_c0_ = 1.0, gpu_vamp_ = 1.0;
Real gpu_ubg1_ = 0.0, gpu_ubg2_ = 0.0;
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void GreshoPUHistory(HistoryData *pdata, Mesh *pm)
//! \brief lowmachpu-1103 (TEST, <problem>/user_hist = true): low-Mach pressure-velocity
//! coupling diagnostics, all volume integrals over the mesh (domain area 1 = means):
//!  pcb_ij  = int (p/p0 - 1)(-1)^(gi+gj)      odd-even pressure mode
//!  pdev    = int (p/p0 - 1)                   mean pressure drift
//!  pdev2   = int (p/p0 - 1)^2
//!  uxcb_i  = int (vx - ubg1)(-1)^gi / c0      odd-even x-velocity (divergence) mode
//!  uycb_j  = int (vy - ubg2)(-1)^gj / c0      odd-even y-velocity (divergence) mode
//!  ukern   = int (vx - ubg1)(-1)^(gi+gj) / c0
//!  dpv1/dpv2/dpvm = int |p - p_G|, (p - p_G)^2, (p - p_G) with p_G the (vortex_amp-scaled)
//!                   exact Gresho pressure (meaningful with the vortex at rest)
//!  fRk/fIk = int (p/p0 - 1) cos/sin(2 pi k x1/L1), k = 1, 2, 4, 8 (acoustic pulse)
void GreshoPUHistory(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &w0_ = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  // w0(IEN) is the internal energy density e = p/(gamma - 1)
  const Real gm1 = (is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                           : pmbp->phydro->peos->eos_data.gamma) - 1.0;
  pdata->nhist = 18;
  const char *lab[18] = {"pcb_ij", "pdev", "pdev2", "uxcb_i", "uycb_j", "ukern",
                         "dpv1", "dpv2", "dpvm", "fR1", "fI1", "fR2", "fI2", "fR4", "fI4",
                         "fR8", "fI8", "vol"};
  for (int n=0; n<18; ++n) pdata->label[n] = lab[n];
  const Real p0 = gpu_p0_, c0 = gpu_c0_, vamp = gpu_vamp_;
  const Real ubg1 = gpu_ubg1_, ubg2 = gpu_ubg2_;
  const Real x1mn = pm->mesh_size.x1min, l1 = pm->mesh_size.x1max - pm->mesh_size.x1min;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, nx1 = indcs.nx1, js = indcs.js, nx2 = indcs.nx2;
  const int ks = indcs.ks, nx3 = indcs.nx3;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1, nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto &size = pmbp->pmb->mb_size;
  array_sum::GlobalSum sum_this_mb;
  Kokkos::parallel_reduce("GreshoPUHist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
    int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks; j += js;
    const Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
    const Real x1v = CellCenterX(i-is, nx1, size.d_view(m).x1min, size.d_view(m).x1max);
    const Real x2v = CellCenterX(j-js, nx2, size.d_view(m).x2min, size.d_view(m).x2max);
    // same global parity rule as the initial condition
    const int gi = static_cast<int>(std::floor((x1v + 1000.0)*nx1/
                                    (size.d_view(m).x1max - size.d_view(m).x1min)));
    const int gj = static_cast<int>(std::floor((x2v + 1000.0)*nx2/
                                    (size.d_view(m).x2max - size.d_view(m).x2min)));
    const Real si = (gi % 2 == 0) ? 1.0 : -1.0;
    const Real sj = (gj % 2 == 0) ? 1.0 : -1.0;
    const Real sij = ((gi + gj) % 2 == 0) ? 1.0 : -1.0;
    const Real p = gm1*w0_(m,IEN,k,j,i);
    const Real dp = p/p0 - 1.0;
    const Real ux = w0_(m,IVX,k,j,i) - ubg1, uy = w0_(m,IVY,k,j,i) - ubg2;
    const Real rad = std::sqrt(x1v*x1v + x2v*x2v);
    const Real pg = p0 + vamp*vamp*(GreshoPres(rad, p0) - p0);
    const Real ph = 2.0*M_PI*(x1v - x1mn)/l1;
    array_sum::GlobalSum h;
    for (int n=0; n<NHISTORY_VARIABLES; ++n) h.the_array[n] = 0.0;
    h.the_array[0] = vol*dp*sij;
    h.the_array[1] = vol*dp;
    h.the_array[2] = vol*dp*dp;
    h.the_array[3] = vol*ux*si/c0;
    h.the_array[4] = vol*uy*sj/c0;
    h.the_array[5] = vol*ux*sij/c0;
    h.the_array[6] = vol*std::fabs(p - pg);
    h.the_array[7] = vol*(p - pg)*(p - pg);
    h.the_array[8] = vol*(p - pg);
    h.the_array[9]  = vol*dp*std::cos(ph);     h.the_array[10] = vol*dp*std::sin(ph);
    h.the_array[11] = vol*dp*std::cos(2.0*ph); h.the_array[12] = vol*dp*std::sin(2.0*ph);
    h.the_array[13] = vol*dp*std::cos(4.0*ph); h.the_array[14] = vol*dp*std::sin(4.0*ph);
    h.the_array[15] = vol*dp*std::cos(8.0*ph); h.the_array[16] = vol*dp*std::sin(8.0*ph);
    h.the_array[17] = vol;
    mb_sum += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this_mb));
  for (int n=0; n<pdata->nhist; ++n) pdata->hdata[n] = sum_this_mb.the_array[n];
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::GreshoMHD(ParameterInput *pin, const bool restart)
//! \brief Problem generator for the magnetized Gresho vortex.  Assumes the domain is
//! [-0.5,0.5] x [-0.5,0.5], so that the vortex is centered on the origin.

void ProblemGenerator::GreshoMHD(ParameterInput *pin, const bool restart) {
  if (restart) return;

  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  // lowmach-noise-1030: also runs in hydro (<hydro> block, no field)
  const bool is_mhd = (pmbp->pmhd != nullptr);
  if (!is_mhd && pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Gresho vortex test needs an <mhd> or a <hydro> block" << std::endl;
    exit(EXIT_FAILURE);
  }

  EOS_Data &eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  Real gm1 = eos.gamma - 1.0;

  Real d0 = pin->GetOrAddReal("problem", "d0", 1.0);
  Real mach = pin->GetOrAddReal("problem", "mach", 0.01);
  Real beta = pin->GetOrAddReal("problem", "beta", 100.0);
  // peak rotation speed is 1, so M = 1/c_s fixes the background pressure
  Real p0 = d0/(eos.gamma*mach*mach);
  Real b0mag = is_mhd ? std::sqrt(2.0*p0/beta) : 0.0;
  // lowmach-noise-1030 (TEST keys, read only when set): uniform field at angle b_angle
  // (degrees) from x1 in the x1-x2 plane (default 0 = along x1, as before), and a
  // grid-scale velocity seed vpert_grid*(-1)^(i+j) (in units of the peak speed 1) on
  // both velocity components, or with vpert_rand > 0 a white-noise seed of that rms
  // (cell-hashed, reproducible)
  Real bang = 0.0, vgrid = 0.0, vrand = 0.0;
  if (pin->DoesParameterExist("problem", "b_angle")) {
    bang = pin->GetReal("problem", "b_angle")*M_PI/180.0;
  }
  if (pin->DoesParameterExist("problem", "vpert_grid")) {
    vgrid = pin->GetReal("problem", "vpert_grid");
  }
  if (pin->DoesParameterExist("problem", "vpert_rand")) {
    vrand = pin->GetReal("problem", "vpert_rand");
  }
  // lowmach2-1102 (TEST keys, read only when set): vortex_amp scales the vortex
  // (velocity x a, pressure perturbation x a^2; 0 = no vortex), ubg1/ubg2 add a uniform
  // flow, and three grid-scale seeds with the global cell parity gi, gj:
  //   cb_p:    p  *= 1 + cb_p (-1)^gi        (pressure checkerboard along x1)
  //   cb_un:   vx += cb_un (-1)^gi           (normal-velocity checkerboard along x1)
  //   cb_kern: vx += cb_kern (-1)^(gi+gj)    (the corner-divergence kernel mode)
  Real vamp = 1.0, ubg1 = 0.0, ubg2 = 0.0, cbp = 0.0, cbun = 0.0, cbkern = 0.0;
  if (pin->DoesParameterExist("problem", "vortex_amp")) {
    vamp = pin->GetReal("problem", "vortex_amp");
  }
  if (pin->DoesParameterExist("problem", "ubg1")) ubg1 = pin->GetReal("problem", "ubg1");
  if (pin->DoesParameterExist("problem", "ubg2")) ubg2 = pin->GetReal("problem", "ubg2");
  if (pin->DoesParameterExist("problem", "cb_p")) cbp = pin->GetReal("problem", "cb_p");
  if (pin->DoesParameterExist("problem", "cb_un")) cbun = pin->GetReal("problem", "cb_un");
  if (pin->DoesParameterExist("problem", "cb_kern")) {
    cbkern = pin->GetReal("problem", "cb_kern");
  }
  // alfmach-1102 (TEST key, read only when set): cb_ut: vy += cb_ut (-1)^gi, a 2dx
  // TRANSVERSE-velocity checkerboard along x1 (an Alfvenic 2dx mode when B has an x1 part)
  Real cbut = 0.0;
  if (pin->DoesParameterExist("problem", "cb_ut")) cbut = pin->GetReal("problem", "cb_ut");
  // lowmachpu-1103 (TEST keys, read only when set):
  //   cb_pij: p  *= 1 + cb_pij (-1)^(gi+gj)  (2-D odd-even pressure checkerboard)
  //   cb_uy:  vy += cb_uy (-1)^gj            (with cb_un: the odd-even divergence mode)
  //   pulse_amp, pulse_w, pulse_x0, pulse_dir: a plane Gaussian sound pulse along x1,
  //     dp = pulse_amp p0 exp(-(dx/pulse_w)^2) (nearest periodic image), isentropic
  //     drho = dp/c0^2, and dvx = pulse_dir dp/(d0 c0) (+1/-1 = simple right/left wave,
  //     0 = pure pressure pulse that splits in two)
  Real cbpij = 0.0, cbuy = 0.0, pamp = 0.0, pw = 0.05, px0 = 0.0, pdir = 1.0;
  if (pin->DoesParameterExist("problem", "cb_pij")) cbpij = pin->GetReal("problem", "cb_pij");
  if (pin->DoesParameterExist("problem", "cb_uy")) cbuy = pin->GetReal("problem", "cb_uy");
  if (pin->DoesParameterExist("problem", "pulse_amp")) {
    pamp = pin->GetReal("problem", "pulse_amp");
    pw = pin->GetOrAddReal("problem", "pulse_w", 0.05);
    px0 = pin->GetOrAddReal("problem", "pulse_x0", 0.0);
    pdir = pin->GetOrAddReal("problem", "pulse_dir", 1.0);
  }
  // lowmachpu-1103 <problem>/pulse_fast (TEST key, read only when named; MHD with the
  // field along x2, b_angle = 90): the pulse is the exact linear FAST magnetosonic simple
  // wave along x1 (k perpendicular to B), dp = pulse_amp p0 g, drho = dp/c0^2,
  // dB2 = B0 drho/d0, dvx = pulse_dir c_f drho/d0, c_f^2 = c0^2 + B0^2/d0
  bool pfast = false;
  if (pin->DoesParameterExist("problem", "pulse_fast")) {
    pfast = pin->GetBoolean("problem", "pulse_fast");
  }
  const Real c0s = std::sqrt(eos.gamma*p0/d0);
  const Real l1m = pmy_mesh_->mesh_size.x1max - pmy_mesh_->mesh_size.x1min;
  gpu_p0_ = p0; gpu_d0_ = d0; gpu_c0_ = c0s; gpu_vamp_ = 1.0;
  const Real cfs = std::sqrt(SQR(c0s) + (is_mhd ? b0mag*b0mag/d0 : 0.0));
  const Real pfb = (pfast && is_mhd) ? b0mag : 0.0;   // dB2 = pfb drho/d0
  if (pfast && is_mhd && std::fabs(bang - 0.5*M_PI) > 1.0e-12) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<problem>/pulse_fast needs b_angle = 90 (field along x2)" << std::endl;
    exit(EXIT_FAILURE);
  }
  const bool lm2 = (cbpij != 0.0) || (cbuy != 0.0) || (pamp != 0.0) || (cbut != 0.0) || (vamp != 1.0) || (ubg1 != 0.0) || (ubg2 != 0.0) || (cbp != 0.0) ||
                   (cbun != 0.0) || (cbkern != 0.0);
  Real b1 = b0mag*std::cos(bang), b2 = b0mag*std::sin(bang), b3 = 0.0;
  // ... and with b_incl (degrees from x3) set, a field out of the x1-x2 plane, e.g.
  // b_incl = 0: B along x3, perpendicular to every face normal of the 2-D mesh
  if (pin->DoesParameterExist("problem", "b_incl")) {
    const Real binc = pin->GetReal("problem", "b_incl")*M_PI/180.0;
    b1 = b0mag*std::sin(binc)*std::cos(bang);
    b2 = b0mag*std::sin(binc)*std::sin(bang);
    b3 = b0mag*std::cos(binc);
  }

  // alfmach-1102 <problem>/b_phi (TEST key, read only when named): in place of the
  // uniform field, an AZIMUTHAL field ALIGNED with the flow, B_phi = B0 u_phi(r) (peak
  // B0 = sqrt(2 p0/beta) at r = 0.2, zero for r > 0.4), so the Alfven Mach number
  // M_A = u_phi/v_A = sqrt(d0)/B0 is the same at every radius of the vortex.  u x B = 0,
  // so the induction equation is stationary, and the pressure takes the hoop stress:
  //   dp_T/dr = (d0 - B0^2) u_phi^2/r,  p = p_Gresho(r) + B0^2 (G(0.4) - G(r))
  //                                          - B0^2 u_phi^2/2,   G = p_Gresho - p0
  // (d0 = 1, as the hydro profile above), with the outside pressure unchanged.  The face
  // field is the curl of A_z = -B0 int_0^r u_phi dr' (divergence-free to round-off).
  // An exact stationary state for any M and M_A (Leidi et al. 2022 use the same idea).
  bool bphi = false;
  if (pin->DoesParameterExist("problem", "b_phi")) {
    bphi = pin->GetBoolean("problem", "b_phi");
  }
  if (bphi) { b1 = 0.0; b2 = 0.0; b3 = 0.0; }
  const Real bphi0 = bphi ? b0mag : 0.0;
  const Real gout = GreshoPres(1.0, 0.0);

  gpu_vamp_ = vamp; gpu_ubg1_ = ubg1; gpu_ubg2_ = ubg2;
  if (user_hist) user_hist_func = GreshoPUHistory;

  // capture variables for kernel
  auto &indcs = pmy_mesh_->mb_indcs;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;

  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &size = pmbp->pmb->mb_size;

  par_for("pgen_gresho", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real rad = std::sqrt(x1v*x1v + x2v*x2v);
    Real vphi = GreshoVphi(rad);
    Real vx = 0.0, vy = 0.0;
    if (rad > 0.0) {
      vx = -vphi*x2v/rad;
      vy =  vphi*x1v/rad;
    }
    if (vgrid != 0.0) {
      // global parity from the cell centre (same on every MeshBlock layout)
      const int gi = static_cast<int>(std::floor((x1v + 1000.0)*nx1/(x1max - x1min)));
      const int gj = static_cast<int>(std::floor((x2v + 1000.0)*nx2/(x2max - x2min)));
      const Real sgn = ((gi + gj) % 2 == 0) ? 1.0 : -1.0;
      vx += vgrid*sgn;
      vy += vgrid*sgn;
    }
    if (vrand != 0.0) {
      // cell-hashed uniform deviates in [-1,1), rms 1/sqrt(3)
      const uint64_t gi = static_cast<uint64_t>((x1v + 1000.0)*nx1/(x1max - x1min));
      const uint64_t gj = static_cast<uint64_t>((x2v + 1000.0)*nx2/(x2max - x2min));
      uint64_t h = (gi*0x9E3779B97F4A7C15ULL)
                   ^ ((gj + 0x632BE59BD9B4E019ULL)*0xBF58476D1CE4E5B9ULL);
      h ^= h >> 31; h *= 0x94D049BB133111EBULL; h ^= h >> 29;
      const Real r1 = static_cast<Real>(h >> 11)*(1.0/9007199254740992.0);
      h *= 0xD6E8FEB86659FD93ULL; h ^= h >> 32;
      const Real r2 = static_cast<Real>(h >> 11)*(1.0/9007199254740992.0);
      vx += vrand*std::sqrt(3.0)*(2.0*r1 - 1.0);
      vy += vrand*std::sqrt(3.0)*(2.0*r2 - 1.0);
    }

    Real pres = GreshoPres(rad,p0);
    Real emag = 0.5*b0mag*b0mag;
    if (bphi0 != 0.0) {
      pres += bphi0*bphi0*(gout - GreshoPres(rad, 0.0) - 0.5*SQR(GreshoVphi(rad)));
      // cell-centred field = mean of the two face fields, from A_z at the 4 corners
      const Real x1l = LeftEdgeX(i-is, nx1, x1min, x1max);
      const Real x1r = LeftEdgeX(i-is+1, nx1, x1min, x1max);
      const Real x2l = LeftEdgeX(j-js, nx2, x2min, x2max);
      const Real x2r = LeftEdgeX(j-js+1, nx2, x2min, x2max);
      auto az = [=](Real x, Real y) { return -bphi0*GreshoVint(std::sqrt(x*x + y*y)); };
      const Real dx1 = x1r - x1l, dx2 = x2r - x2l;
      const Real bxl = (az(x1l, x2r) - az(x1l, x2l))/dx2;
      const Real bxr = (az(x1r, x2r) - az(x1r, x2l))/dx2;
      const Real byl = -(az(x1r, x2l) - az(x1l, x2l))/dx1;
      const Real byr = -(az(x1r, x2r) - az(x1l, x2r))/dx1;
      emag = 0.5*(SQR(0.5*(bxl + bxr)) + SQR(0.5*(byl + byr)));
    }
    Real dens = d0;
    if (lm2) {
      const int gi = static_cast<int>(std::floor((x1v + 1000.0)*nx1/(x1max - x1min)));
      const int gj = static_cast<int>(std::floor((x2v + 1000.0)*nx2/(x2max - x2min)));
      const Real si = (gi % 2 == 0) ? 1.0 : -1.0;
      const Real sj = (gj % 2 == 0) ? 1.0 : -1.0;
      const Real sij = ((gi + gj) % 2 == 0) ? 1.0 : -1.0;
      vx = vamp*vx + ubg1 + cbun*si + cbkern*sij;
      vy = vamp*vy + ubg2 + cbut*si + cbuy*sj;
      pres = (p0 + vamp*vamp*(pres - p0))*(1.0 + cbp*si)*(1.0 + cbpij*sij);
      if (pamp != 0.0) {
        Real dxp = x1v - px0;
        dxp -= l1m*std::round(dxp/l1m);
        const Real dpp = pamp*p0*std::exp(-SQR(dxp/pw));
        pres += dpp;
        dens += dpp/SQR(c0s);
        if (pfb != 0.0) {
          vx += pdir*cfs*dpp/(SQR(c0s)*d0);
          const Real by = b2 + pfb*dpp/(SQR(c0s)*d0);
          emag = 0.5*(b1*b1 + by*by + b3*b3);
        } else {
          vx += pdir*dpp/(d0*c0s);
        }
      }
    }
    u0(m,IDN,k,j,i) = dens;
    u0(m,IM1,k,j,i) = dens*vx;
    u0(m,IM2,k,j,i) = dens*vy;
    u0(m,IM3,k,j,i) = 0.0;
    // total energy, with the uniform field contributing 0.5*B0^2 everywhere
    u0(m,IEN,k,j,i) = pres/gm1 + 0.5*dens*(vx*vx + vy*vy)
                      + emag;
  });
  if (!is_mhd) return;

  // uniform field at angle b_angle from x1 in the x1-x2 plane (divergence-free)
  auto &b0 = pmbp->pmhd->b0;
  par_for("pgen_gresho_b", DevExeSpace(), 0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m,k,j,i) = b1;
    b0.x2f(m,k,j,i) = b2;
    b0.x3f(m,k,j,i) = b3;
    if (i==ie) {
      b0.x1f(m,k,j,i+1) = b1;
    }
    if (j==je) {
      b0.x2f(m,k,j+1,i) = b2;
    }
    if (pfb != 0.0) {
      // x2-face field depends on x1 only: divergence-free; value at the cell centre x1
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      const Real x1v = CellCenterX(i-is, indcs.nx1, x1min, x1max);
      Real dxp = x1v - px0;
      dxp -= l1m*std::round(dxp/l1m);
      const Real by = b2 + pfb*pamp*p0*std::exp(-SQR(dxp/pw))/(SQR(c0s)*d0);
      b0.x2f(m,k,j,i) = by;
      if (j==je) b0.x2f(m,k,j+1,i) = by;
    }
    if (bphi0 != 0.0) {
      Real &x1min = size.d_view(m).x1min;
      Real &x1max = size.d_view(m).x1max;
      Real &x2min = size.d_view(m).x2min;
      Real &x2max = size.d_view(m).x2max;
      const int nx1 = indcs.nx1, nx2 = indcs.nx2;
      const Real x1l = LeftEdgeX(i-is, nx1, x1min, x1max);
      const Real x1r = LeftEdgeX(i-is+1, nx1, x1min, x1max);
      const Real x2l = LeftEdgeX(j-js, nx2, x2min, x2max);
      const Real x2r = LeftEdgeX(j-js+1, nx2, x2min, x2max);
      auto az = [=](Real x, Real y) { return -bphi0*GreshoVint(std::sqrt(x*x + y*y)); };
      b0.x1f(m,k,j,i) = (az(x1l, x2r) - az(x1l, x2l))/(x2r - x2l);
      b0.x2f(m,k,j,i) = -(az(x1r, x2l) - az(x1l, x2l))/(x1r - x1l);
      if (i==ie) b0.x1f(m,k,j,i+1) = (az(x1r, x2r) - az(x1r, x2l))/(x2r - x2l);
      if (j==je) b0.x2f(m,k,j+1,i) = -(az(x1r, x2r) - az(x1l, x2r))/(x1r - x1l);
    }
    if (k==ke) {
      b0.x3f(m,k+1,j,i) = b3;
    }
  });

  return;
}
