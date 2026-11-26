//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file star_vet.cpp
//! \brief Problem generator for stellar envelopes with the VET radiation module
//! (<rad_vet>), on spherical-polar wedges and the cubed sphere, used with
//! `cmake -D PROBLEM=star_vet`.  Inputs: inputs/stars/star_*.athinput,
//! star3d_*.athinput.
//!
//! It sets up any star from a 1-D column (problem/star_ic_file) and opacity tables:
//! gas and radiation in the column's balance, gravity (star_gravity, star_gm*), walls
//! (star_bc_inner, star_bc_outer, star_wall_*), seeds (star_seed*, star_ic_profile*),
//! a top sponge (star_sponge_*), a frozen or adaptive MLT++ scaffold (mlt_*) and an
//! optional time-dependent base (star_drive.hpp).  Hydro only (<mhd> is refused):
//!   * the gas-only tabulated EOS (<hydro>/eos_radiation = false, enforced by the module)
//!     or an ideal gas (eos = ideal; <units>/mu sets its kelvin scale, which
//!     <rad_vet>/temp_unit_kelvin must equal; no EOS range check then)
//!     and the Rosseland + Planck tables (problem/star_opac_table, star_planck_table);
//!   * the initial column star_ic_file, columns  r rho eint F_r  (cgs, ascending r):
//!     F_r is the DIFFUSIVE part of the luminosity flux, F_r = L/(4 pi r^2) - F_MLT, so
//!     the radiation force kappa_R F_r/c that the gas is in balance with is known
//!     everywhere; E = a T(rho,eint)^4 is set here from the run's own EOS;
//!   * force_reference = wb_arad with the reference acceleration a_ref = kappa_t F_r/c
//!     of that column, delivered through the EFFECTIVE potential of the x1
//!     well-balanced pair Phi_eff = Phi - int a_ref dr (<rad_vet> then applies only the
//!     residual).  The pgen also gives the reference WORK (fref_wsplit_ok), which is what
//!     lets <rad_vet>/force_reference_work = split be accepted;
//!   * x1 boundaries (ix1_bc = ox1_bc = user): inner = closed hydro wall carrying the
//!     scaled initial profile, the luminosity entering through the face flux
//!     (<rad_vet>/implicit_bc_x1min = flux, implicit_flux_x1min = L/(4 pi r_in^2));
//!     outer = hydro outflow with NO inflow (mass may leave through the top, none
//!     enters), Marshak top (implicit_bc_x1max = marshak);
//!   * a startup CHECK that the whole initial column, ghost margin included, lies inside
//!     the EOS table (rho, T) and inside the Rosseland/Planck table grid (fatal);
//!   * history columns (user_hist): face luminosities, L_in, E_rad, the total energy of
//!     gas + radiation, the mass through both x1 faces, the interior kinetic sums, the
//!     mean v_r^2 of the wall cells and the mean Picard passes per implicit solve since
//!     the previous history output;
//!   * an optional multi-mode entropy/temperature seed (star_seed_nk, star_seed_rad,
//!     star_seed_signs, star_seed_kdist).
//!   * problem/star_seed_band = n2 (default) | old: n2 places the fresh-start seed where
//!     the initial column is most unstable (Brunt-Vaisala N^2 < 0, gas + radiation,
//!     general EOS), real Y_lm with l from pi r/H_p clipped to the grid, amplitude
//!     star_seed default 1e-3; old = the seed of the keys as named (star_seed default 0),
//!     bitwise as before.  See StarN2Band().
//!   * problem/star_gm_column (default off): the gravity of the enclosed
//!     mass of the initial column instead of the point mass: G m(r) = star_gm + G int 4
//!     pi r^2 rho_IC dr from r_in (star_gm = G m(r_in)), a static monopole frozen at t =
//!     0.
//!   * problem/star_base_heat (default off): a wall inside a convective
//!     zone: the scaffold's wall-face convective luminosity is kept as a base heating
//!     layer while the frozen scaffold is ramped off.
//!
//! RESTARTS: the pgen is not skipped on a restart: the column, the tables, the
//! potentials, the reference acceleration and the hooks are rebuilt from the input file;
//! what cannot be (Picard window, star_ic_equil a_ref, automatic mlt_ramp, adaptive MLT
//! scaffold, drive ledger) travels in the file's pgen block.
//!
//! CUDA-safe: no lambdas inside kernels, no host reads of device Views (the fine column
//! is built on the host and copied), no class members inside kernels.

#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "hydro/hydro_am.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "utils/wb_background.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/cell_locations.hpp"
#include "rad_vet/rad_vet.hpp"
#include "rad_vet/rad_vet_closure.hpp"
#include "rad_vet/rad_vet_implicit.hpp"
#include "rad_vet/rad_vet_opacity.hpp"
#include "rad_vet/rays/sc_segment.hpp"
#include "rad_vet/rad_vet_fluid.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"
#include "pgen/star_drive.hpp"
#include "driver/dt_retry.hpp"
#include "utils/eint_from_cons.hpp"
#include "utils/imex_sdirk2.hpp"
#include "utils/repro_sum.hpp"
#include "utils/p_ledger.hpp"
#include "utils/recipe.hpp"

// <hydro>/fofc defaults to TRUE for this problem generator (every
// stellar/box input ran it explicitly); see hydro::fofc_pgen_default (hydro.hpp).
// <mhd>/fofc likewise (mhd::fofc_pgen_default, mhd.hpp).
namespace {
[[maybe_unused]] const bool kFofcPgenDefault = (hydro::fofc_pgen_default = true);
[[maybe_unused]] const bool kStarPgen = (recipe::star_pgen = true);
[[maybe_unused]] const bool kMhdFofcPgenDefault = (mhd::fofc_pgen_default = true);
// this pgen fills the potential and carries rho*Phi (both gravity modes need
// <hydro>/etotgrav = true): the new etotgrav default applies
// (hydro::etotgrav_pgen_default)
[[maybe_unused]] const bool kEtgPgenDefault = (hydro::etotgrav_pgen_default = true);
[[maybe_unused]] const bool kMhdEtgPgenDefault = (mhd::etotgrav_pgen_default = true);
}  // namespace

namespace {
// the initial column on a uniform fine radial grid (device), read by the x1 BCs
DvceArray1D<Real> star_rho_, star_eint_;
Real star_rlo_ = 0.0, star_dr_ = 1.0;
int star_nf_ = 0;
Real star_gm_ = 0.0, star_rin_ = 1.0, star_rint_ = 0.0, star_fin_ = 0.0;
// problem/star_gm_column (read only when named, default false = the point
// mass GM = star_gm, bit for bit): the gravity of the ENCLOSED MASS of the initial
// column, G m(r) = star_gm + G int_{r_in}^{r} 4 pi r'^2 rho_IC dr' (star_gm = G m(r_in),
// the mass inside the wall), FROZEN at t = 0 (a static monopole, rebuilt from
// star_ic_file on every start).  star_phig_ = Phi(r) = int_{r_in}^{r} G m/r'^2 dr' and
// star_gmr_ = G m(r) on the uniform fine grid of the column (star_rlo_, star_dr_,
// star_nf_).
bool star_gmc_ = false;
DvceArray1D<Real> star_phig_, star_gmr_;
// problem/star_base_heat [cm] (read only when named, default none): the wall
// lies inside a CONVECTIVE zone, so the frozen MLT scaffold's wall-face flux F_MLT(r_in)
// is the convective luminosity entering the domain.  While the scaffold is ramped off
// (mlt_ramp_start/_time, weight w) that luminosity is kept as a BASE HEATING layer: the
// deposit is w (-div A F_MLT) + (1 - w) (-div A F_b), F_b(r) = F_MLT(r_in) (r_in/r)^2
// max(0, 1 - (r - r_in)/star_base_heat) on the x1 faces, so A F_b = L_MLT(r_in) at the
// wall face and 0 above r_in + star_base_heat.  Frozen closure only.  (Injecting it as
// VET face flux instead puts a radiation force of ~Gamma_F ~ 100 g on the wall cells.)
Real star_bheat_ = 0.0;
DvceArray1D<Real> star_fb_;
Real star_sp_rate_ = 0.0, star_sp_r0_ = 0.0, star_rtop_ = 1.0;
bool star_zflux_ = true;
DvceArray3D<Real> star_min_;   // (the top inflow mass flux removed, StarWallFlux)
// problem/star_wb_grav_window (read only when named, default false): the x1
// WB gravity source honours the fluid's wb_rmax / wb_phimax window -- the cells
// the reconstruction treats as plain (Hydro::wb_rmax) take the isodensity walk
// p -/+ rho dPhi_wb instead of the cached hydrostatic walk (see StarGravity)
bool star_wbwin_ = false;
// he-wind-bc: problem/star_bc_inner = inflow (star_binf_) and star_bc_outer = outflow
// (star_bout_); defaults wall / noinflow = the he-presn-m1 behaviour
bool star_binf_ = false, star_bout_ = false;
// the initial column's rho, eint at r_in (scales of the table ghosts)
Real star_drv_d0_ = 1.0, star_drv_e0_ = 1.0;
// problem/star_bc_outer = hse (star_bhse_): the top ghosts are the HYDROSTATIC
// continuation of the last active cell (StarBC); star_bc_hse_gmax clamps the Eddington
// factor used there, star_bc_hse_tgrad continues T along the interior's dT/dPhi_eff
// (polytropic WB walk) instead of holding it (isothermal).  problem/star_sponge_mode =
// radial (star_sp_rad_) makes the top sponge damp v_r only (default all = v).
// problem/star_sponge_dmax (default 0 = off): the sponge also acts, at its full rate, on
// every cell with rho < star_sponge_dmax (the floor gas above the photosphere, which can
// not be hydrostatic and falls freely)
// problem/star_sponge_form (default lawson on a fresh start, recorded; a
// restart file without the key keeps exp, NOTE): exp = u0 *= exp(-bdt rate w)
// after each RK stage (1st order in time); lawson = the integrating-factor form of the
// SSP-RK stages (exact for the decay, 2nd order, weights >= 0): the stage's u^n part
// decays over c_s dt, the rest over (c_s - c_{s-1}) dt; w and the rho < dmax flag from
// u^n (u1), the same on every stage
bool star_bhse_ = false, star_bhse_tg_ = false;
bool star_sp_rad_ = false, star_sp_law_ = false;
// problem/star_bc_hse_flux = face (default on a fresh start) | cell: the radial flux in
// the hse top's Gamma.  face = the transported comoving face flux (f0x1, mean of the top
// cell's two x1 faces); cell = the old clipped lab-frame cell F1.  A restart whose
// embedded input lacks the key keeps cell (no silent change mid-run).
bool star_bhse_face_ = true;
Real star_bhse_gmax_ = 0.9, star_sp_dmax_ = 0.0;
// problem/star_ic_balance: the discretely balanced initial cells (all x1 cells incl.
// ghosts of the ONE MeshBlock along x1; every block has the same radial grid)
bool star_bal_ = false;
DvceArray1D<Real> star_bd_, star_be_;
// problem/star_ic_rad_atm_post (StarIcPost): the thin atmosphere's first
// mesh cell above R_ph (x1 index incl. ghosts) and the column's floor density
bool star_atp_ = false;
ParameterInput *star_pin_ = nullptr;   // the run's input (StarIcPost records a_ref)
int star_atp_ia_ = -1;
Real star_atp_dfl_ = 0.0;
// (Ma, Bildsten & Jiang 2026 setup, all three default off):
// problem/star_gravity = pointmass (star_pm_): NO well-balanced pair; the point-mass
//   momentum source in the Athena++ form rho GM <1/r>_vol/x1v (etotgrav keeps the
//   energy);
// problem/star_bc_inner = zerocontact (star_bzc_): every inner ghost holds the UNSCALED
//   initial column (rho, eint), v = 0, except the first ghost layer whose v_r is solved
//   so that the HLLC contact speed of the code's own PLM states at the wall face is 0;
// problem/star_bc_outer = copy (star_bcopy_): the edge cell's rho, eint, v_t copied,
//   v_r copied where it points out and 0 where it points in.
bool star_pm_ = false, star_bzc_ = false, star_bcopy_ = false;
// problem/star_bc_outer_r2 (read only when named; default false = bitwise;
// needs star_bc_outer = copy): bsg_star.cpp Outflow_X1 exactly: edge v_r < 0 -> rho
// copied, v_r = 0; else (v_r >= 0) ghost rho = previous ghost rho (r_{g-1}/r_g)^2 and v_r
// copied; P (here eint, ideal gas) and v_t copied in both branches
bool star_bcr2_ = false;
// problem/mlt_flux_frozen: the frozen MLT flux of the IC on the x1 faces (index = face)
bool star_mlt_ = false;
DvceArray1D<Real> star_fm_;
// problem/mlt_ramp_start, mlt_ramp_time: the time-only weight w(t) on the frozen MLT
// deposit (1 before the start, cosine to 0 over the ramp time, 0 after; start < 0 = no
// ramp, w = 1 always); star_mw_ = the w now in <rad_vet>/esrc; star_lmx_ = max over the
// x1 faces of the UNSCALED L_MLT = F_MLT r^2 Omega (the wedge's solid angle)
Real star_mrs_ = -1.0, star_mrt_ = 0.0, star_mw_ = 1.0, star_lmx_ = 0.0;
// mlt_ramp_start/_time = auto (default) resolved from the IC column's
// convective growth rates at the fresh start; star_mr_store_ = the two numbers travel in
// the restart file (trailer kStarRmpTag, the LAST one: start, time, tag)
bool star_mr_store_ = false;
constexpr char kStarRmpTag[8] = {'H', 'S', 'R', 'A', 'M', 'P', '0', '1'};
// the IC column's per-cell N^2 profile (StarN2Compute), computed once
struct StarN2Col {
  bool done = false;
  int is = 0, iex = -1, ipk = 0, iu0 = 0, iu1 = 0, nu = 0;
  Real smax = 0.0, med = 0.0;
  std::vector<Real> r, sg, dm, pt, dd, gg, tt, g1;
};
StarN2Col star_n2_;
// problem/star_esrc_const, star_esrc_rmax (the esrc test, kept when w is re-applied)
Real star_esc_ = 0.0, star_esrmx_ = 0.0;
// problem/mlt_closure = adaptive (default frozen): the ADAPTIVE SHELL-MEAN DEFICIT
// closure.  The applied subgrid flux F_sub (x1 faces) relaxes toward
// s(r) max(0, L - L_rad - L_conv,res)/A over problem/mlt_relax_time, updated every
// cycle at its START (before any stage), from
// shell sums over all ranks: L_rad = sum A F_0 (the VET comoving face flux), L_conv,res =
// the FLUCTUATION flux A (<v_r'(e+p)_gas'> + (4/3) <v_r'E'>) of the cell shells averaged
// to the face.  s = clamp(fmlt_IC/problem/mlt_closure_taper, 0, 1): the FeCZ of the IC
// (F_MLT > 0) with its edges tapered.  F_sub(t=0) = F_MLT of the IC.  star_fsub_ and the
// update clock star_ad_tl_ are RESTART STATE (the pgen block of pgen.hpp).
bool star_ad_ = false;
// problem/mlt_closure = pressure (opt-in): the same F_sub machinery (restart
// state, relaxation over mlt_relax_time, every cycle), but the target
// is the IC's F_MLT/F as a function of the SHELL-MEAN total pressure, F_sub =
// fm_IC(P_face) x F_tot,IC(face) with F_tot,IC = L/(4 pi r^2) of the IC: the scaffold is
// frozen in P instead of r and moves with the envelope.  star_plnp_/star_pfm_: ln P_tot,
// fm of the IC (ascending in ln P), star_ftot_: F_r + F_MLT of the IC on the faces.
bool star_adp_ = false;
// problem/star_ad_dead_skip (read only when named; default true;
// bitwise):
// under mesh/layout_invariant the per-shell serial sums of StarAdaptiveUpdate (a thread
// per radial index over every column, with the EOS pressure) are overwritten by the
// order-independent repro::ShellSums: skip them
bool star_ad_dskip_ = false;
std::vector<Real> star_plnp_, star_pfm_, star_ftot_;
DvceArray2D<Real> star_shp_;
Real star_ad_tau_ = 1500.0, star_ad_tl_ = 0.0;
Real star_ad_dout_ = 0.0, star_ad_tout_ = 0.0, star_lsub_ = 0.0;
int star_ad_cyc_ = -1;
DvceArray1D<Real> star_fsub_;
// problem/mlt_src_time = stage (default; StarMltC) | start (the old form, 1st order in
// time): the ramp w and the adaptive esrc at the time of each implicit stage, t + c dt,
// c = a21 + g and 1 (hesdirk2 1, 1 with the FSAL slope K1 at 0; sdirk2 g, 1; the
// backward-Euler solves 1), 2nd order: sum b c = 1/2, and the stiffly accurate last stage
// has the source of t + dt.  adaptive: F_sub^n = F^{n-1} + a(t^n - t^{n-1})
// ((tgt^{n-1} + tgt^n)/2 - F^{n-1}) (exponential trapezoid; star_tgp_ = tgt^{n-1}, restart
// state behind F_sub, hdr[1] = 1) and the stage's F = F^n + a(c dt) (tgt^n - F^n)
bool star_mtc_ = false;
int star_tgp_ok_ = 0;
DvceArray1D<Real> star_tgp_;
DvceArray2D<Real> star_shs_;
std::vector<Real> star_tap_, star_r2o_;
char star_ad_file_[1024] = {0};
// Picard counters at the previous history output (StarHist)
Real star_pic_n0_ = 0.0, star_pic_s0_ = 0.0;
// the open Picard window travels in the pgen restart block as a tagged
// trailer (tag, solves, passes since the last history output), so the first user.hst row
// after a restart is the continuous run's
constexpr char kStarPicTag[8] = {'H', 'S', 'P', 'I', 'C', 'W', '0', '1'};
Mesh *star_mesh_ = nullptr;   // the mesh (StarAdaptiveRstWrite: the rad_vet counters)
// problem/star_ic_equil (StarIcEquil; read only when named, default false =
// nothing): the discrete equilibrium of the initial state on the code's own operators.
// star_eq_ on; star_eq_it_ the hook calls so far; star_eq_dfl_ the column's floor
// density; star_eq_a_ a_ref of the active cells is..ie (the restart trailer kStarDicTag);
// star_sfac_ the star_seed factor of every cell, applied after the equilibrium (deferred
// seed)
bool star_eq_ = false, star_eq_srad_ = false;
bool star_eq_seed_ = false, star_eq_last_ = false;
int star_eq_it_ = 0, star_eq_maxit_ = 30;
Real star_eq_tol_ = 1.0e-9, star_eq_dfl_ = 0.0;
std::vector<Real> star_eq_a_;
DvceArray4D<Real> star_sfac_;
constexpr char kStarDicTag[8] = {'H', 'S', 'D', 'I', 'C', 'E', '0', '1'};
// (the removed fixed-profile inner face's wall-luminosity trailer: stripped from old
// restart files)
constexpr char kStarLbcTag[8] = {'H', 'S', 'L', 'B', 'A', 'T', 'H', '1'};

//! the frozen-MLT weight w(t): a function of the time only (restart-safe)
Real StarMltW(const Real t) {
  if (star_mrs_ < 0.0 || t < star_mrs_) return 1.0;
  if (star_mrt_ <= 0.0 || t >= star_mrs_ + star_mrt_) return 0.0;
  return 0.5*(1.0 + cos(M_PI*(t - star_mrs_)/star_mrt_));
}
//! the stage time of the MLT source, c in t + c dt (mlt_src_time), of the hydro stage bdt
Real StarMltC(const Mesh *pm, const Real bdt) {
  if (!star_mtc_) return 0.0;
  const auto *pvet = pm->pmb_pack->pradvet;
  return (pvet->time_scheme == radvet::VET_TIME_HESDIRK2 && pvet->t2_tab != 0 &&
          bdt == pvet->fref_hbeta[0]*pm->dt) ? pvet->t2_a21 + imex_sdirk2::Gamma() : 1.0;
}

//! log-linear interpolation on the fine grid, clamped to its end nodes
KOKKOS_INLINE_FUNCTION
Real StarLogInterp(const DvceArray1D<Real> &a, const Real rlo, const Real dr,
                   const int nf, const Real r) {
  Real x = (r - rlo)/dr;
  int i = static_cast<int>(floor(x));
  i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
  Real w = x - i;
  w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
  return exp((1.0 - w)*log(a(i)) + w*log(a(i+1)));
}

//! linear interpolation on the fine grid, clamped
KOKKOS_INLINE_FUNCTION
Real StarLinInterp(const DvceArray1D<Real> &a, const Real rlo, const Real dr,
                   const int nf, const Real r) {
  Real x = (r - rlo)/dr;
  int i = static_cast<int>(floor(x));
  i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
  Real w = x - i;
  w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
  return (1.0 - w)*a(i) + w*a(i+1);
}

//! (problem/star_ic_profile): a uniform deviate U in [0, 1) that depends only
//! on the GLOBAL cell index (and the key's seed): splitmix64 of the index.  Identical on
//! every layout, device and site (the per-block srand(gid)/rand() stream of the Athena++
//! pgen cannot be reproduced; only its statistics are).
KOKKOS_INLINE_FUNCTION
Real StarCellUniform(const std::uint64_t idx, const std::uint64_t seed) {
  std::uint64_t z = idx + seed*0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL;
  z = (z ^ (z >> 30))*0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27))*0x94D049BB133111EBULL;
  z = z ^ (z >> 31);
  return static_cast<Real>(z >> 11)*(1.0/9007199254740992.0);
}

[[noreturn]] void StarFatal(const std::string &msg, const int line) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << line << std::endl
            << "star_vet: " << msg << std::endl;
  // (_Exit: no static destructors behind live Kokkos views, as ScFatalAt)
  std::_Exit(EXIT_FAILURE);
}

//! the problem key `k` (star_ic_profile*), or its deprecated name star_ic_linhao* when
//! only that one is in the input (older decks and restart files)
std::string StarIcKey(ParameterInput *pin, const std::string &k) {
  const std::string old = "star_ic_linhao" + k.substr(std::strlen("star_ic_profile"));
  return (!pin->DoesParameterExist("problem", k) &&
          pin->DoesParameterExist("problem", old)) ? old : k;
}

//! a restart file without problem/<key>: the old default v applies (NOTE)
void StarOldDefault(ParameterInput *pin, const char *key, const char *v) {
  if (global_variable::restart_run && !pin->DoesParameterExist("problem", key) &&
      global_variable::my_rank == 0) {
    std::cout << "star_vet: NOTE problem/" << key << " = " << v << " (a restart file "
              << "without the key: the old default)" << std::endl;
  }
}

//! the merged opacity table format of box_convection.cpp's ReadOpacityTable (comment
//! lines, one of them "# nT nD lTmin dlT lDmin dlD", then nT*nD values of log10 kappa,
//! T slowest), filled on the host and deep-copied to the device; the grid ranges are
//! returned for the startup range check
void StarReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                          DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD,
                          Real &lt_lo, Real &lt_hi, Real &ld_lo, Real &ld_hi,
                          const Real ld_ext = 1.0e300) {
  std::ifstream f(fname);
  if (!f.good()) StarFatal("cannot open opacity table '" + fname + "'", __LINE__);
  std::string line;
  Real lt0 = 0.0, dlt = 0.0, ld0 = 0.0, dld = 0.0;
  bool have_grid = false;
  std::vector<Real> vals;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      if (!have_grid) {
        std::istringstream ss(line.substr(1));
        int a, b;
        Real c, d, e, g;
        if (ss >> a >> b >> c >> d >> e >> g) {
          nT = a; nD = b; lt0 = c; dlt = d; ld0 = e; dld = g;
          have_grid = true;
        }
      }
      continue;
    }
    vals.push_back(std::stod(line));
  }
  if (!have_grid || static_cast<int>(vals.size()) != nT*nD) {
    StarFatal("opacity table '" + fname + "' has no grid line or the wrong value count",
              __LINE__);
  }
  // problem/star_opac_logd_min (he-wind-bc, default none): the grid is EXTENDED to lower
  // density by nadd columns, log10 kappa continued linearly in log10 rho with the slope
  // of the table's first density interval clamped to [0, 1] (an absorption/scattering
  // opacity per unit mass does not grow as rho falls, and falls at most like rho, the
  // free-free/bound-free limit).  Where the slope lies inside [0, 1] the extension is
  // continuous in value AND slope; the stellar tables are edge-filled (constant in rho)
  // for log R < -8, i.e. flat (slope 0) at the edge for log T >= 4.2.
  if (ld_ext < ld0 - 1.0e-9*dld) {
    const int nadd = static_cast<int>(ceil((ld0 - ld_ext)/dld - 1.0e-9));
    const int nDn = nD + nadd;
    std::vector<Real> vn(static_cast<std::size_t>(nT)*nDn);
    int nclip = 0;
    for (int i=0; i<nT; ++i) {
      const Real v0 = vals[i*nD], v1 = vals[i*nD + 1];
      Real sl = (v1 - v0)/dld;
      if (sl < 0.0 || sl > 1.0) ++nclip;
      sl = (sl < 0.0) ? 0.0 : ((sl > 1.0) ? 1.0 : sl);
      for (int j=0; j<nadd; ++j) vn[i*nDn + j] = v0 - sl*(nadd - j)*dld;
      for (int j=0; j<nD; ++j) vn[i*nDn + nadd + j] = vals[i*nD + j];
    }
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: opacity table '" << fname << "' extended from log10 rho "
                << ld0 << " to " << ld0 - nadd*dld << " (" << nadd << " columns, slope "
                << "clamped to [0,1] in " << nclip << " of " << nT << " T rows)"
                << std::endl;
    }
    vals.swap(vn);
    nD = nDn;
    ld0 -= nadd*dld;
  }
  Kokkos::realloc(tab, nT, nD);
  Kokkos::realloc(lT, nT);
  Kokkos::realloc(lD, nD);
  auto htab = Kokkos::create_mirror_view(tab);
  auto hlT = Kokkos::create_mirror_view(lT);
  auto hlD = Kokkos::create_mirror_view(lD);
  for (int i=0; i<nT; ++i) {
    hlT(i) = lt0 + i*dlt;
    for (int j=0; j<nD; ++j) htab(i,j) = vals[i*nD + j];
  }
  for (int j=0; j<nD; ++j) hlD(j) = ld0 + j*dld;
  Kokkos::deep_copy(tab, htab);
  Kokkos::deep_copy(lT, hlT);
  Kokkos::deep_copy(lD, hlD);
  lt_lo = lt0; lt_hi = lt0 + (nT - 1)*dlt;
  ld_lo = ld0; ld_hi = ld0 + (nD - 1)*dld;
}

//----------------------------------------------------------------------------------------
// (problem/star_ic_rad = sph1d): helpers of the 1-D spherical grey
// radiative-equilibrium solve of the initial column (see UserProblem)

// kappa_t = of + os (the total, = kappa_R: the rays' extinction and the pgen's
// own ckt) of the module's tables at (d, t) [code units], one device thread per point
void StarKtEval(const radvet::VetOpacTab &ot, const std::vector<Real> &d,
                const std::vector<Real> &t, std::vector<Real> &k) {
  const int n = static_cast<int>(d.size());
  DvceArray1D<Real> dd("star_icr_d", n), dt("star_icr_t", n), dk("star_icr_k", n);
  auto hd = Kokkos::create_mirror_view(dd);
  auto ht = Kokkos::create_mirror_view(dt);
  for (int i=0; i<n; ++i) {
    hd(i) = d[i];
    ht(i) = t[i];
  }
  Kokkos::deep_copy(dd, hd);
  Kokkos::deep_copy(dt, ht);
  par_for("star_icr_kt", DevExeSpace(), 0, n-1, KOKKOS_LAMBDA(const int i) {
    Real op, oe, of, os;
    radvet::VetTableOpacities(ot, dd(i), dt(i), op, oe, of, os);
    dk(i) = of + os;
  });
  auto hk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), dk);
  k.resize(n);
  for (int i=0; i<n; ++i) k[i] = hk(i);
}

// mean of chi over an interval with chi log-linear in the path (the arithmetic mean
// where the two ends are equal or one vanishes)
inline Real StarLogMean(const Real c0, const Real c1) {
  if (c0 > 0.0 && c1 > 0.0 && std::fabs(c1 - c0) > 1.0e-8*c0) {
    return (c1 - c0)/std::log(c1/c0);
  }
  return 0.5*(c0 + c1);
}

// one first-order short-characteristic step (the positive weights VetSegW,
// source linear in tau): I = I_up e^-dtau + w0 S_0 + wu S_up
inline Real StarScStep(const Real iup, const Real dtau, const Real su, const Real s0) {
  Real ex, w0, wu;
  VetSegW(dtau, ex, w0, wu);
  return iup*ex + w0*s0 + wu*su;
}

// ONE VET pass on the solve grid rg (ascending; rg[0] = r_b, the inner boundary):
// (a) the p-z formal solution of n.grad I = chi (S - I), S = E (radiative equilibrium,
//     E units: I 4 pi/c), one ray per shell tangent plus ncore core rays p < r_b uniform
//     in p^2, no incoming intensity at the top, the diffusion intensity E_b + 3 mu F_b/c
//     at r_b; J, H, K by the trapezoid in mu (the rays reach every shell in descending
//     mu, so the sums are accumulated ray by ray);
// (b) the moment equation d(f E)/dr + (3f - 1) E/r = -chi F/c with f = K/J, integrated
//     inward (trapezoid; chi log-linear) from E_top = F_top/(c g_top), g = H/J.
// E is replaced by the new E; jout = J of the formal solution (for diagnostics).
void StarSph1dPass(const std::vector<Real> &rg, const std::vector<Real> &chi,
                   const std::vector<Real> &F, const int ncore, const Real cl,
                   std::vector<Real> &E, std::vector<Real> &jout, Real &gtop,
                   Real &ftop) {
  const int N = static_cast<int>(rg.size());
  std::vector<Real> J(N, 0.0), H(N, 0.0), K(N, 0.0), mup(N, -1.0), ap(N), bp(N);
  std::vector<Real> im(N), ip(N), z(N), dtau(N);
  const Real rb = rg[0];
  const int nr = ncore + N;
  for (int q=0; q<nr; ++q) {
    const bool core = (q < ncore);
    const Real p = core ? rb*std::sqrt(static_cast<Real>(q)/ncore) : rg[q - ncore];
    const int k0 = core ? 0 : (q - ncore);
    for (int k=k0; k<N; ++k) {
      z[k] = (!core && k == k0) ? 0.0
                                : std::sqrt(std::fmax((rg[k] - p)*(rg[k] + p), 0.0));
    }
    for (int k=k0; k<N-1; ++k) dtau[k] = (z[k+1] - z[k])*StarLogMean(chi[k], chi[k+1]);
    im[N-1] = 0.0;
    for (int k=N-2; k>=k0; --k) im[k] = StarScStep(im[k+1], dtau[k], E[k+1], E[k]);
    if (core) {
      const Real mub = std::sqrt(std::fmax(1.0 - (p/rb)*(p/rb), 0.0));
      ip[k0] = E[0] + 3.0*mub*F[0]/cl;
    } else {
      ip[k0] = im[k0];
    }
    for (int k=k0+1; k<N; ++k) ip[k] = StarScStep(ip[k-1], dtau[k-1], E[k-1], E[k]);
    for (int k=k0; k<N; ++k) {
      const Real mu = z[k]/rg[k], a = ip[k], b = im[k];
      if (mup[k] >= 0.0) {
        const Real dm = mup[k] - mu;
        J[k] += 0.25*dm*((ap[k] + bp[k]) + (a + b));
        H[k] += 0.25*dm*((ap[k] - bp[k])*mup[k] + (a - b)*mu);
        K[k] += 0.25*dm*((ap[k] + bp[k])*mup[k]*mup[k] + (a + b)*mu*mu);
      }
      mup[k] = mu;
      ap[k] = a;
      bp[k] = b;
    }
  }
  std::vector<Real> f(N), y(N);
  for (int k=0; k<N; ++k) f[k] = K[k]/J[k];
  gtop = H[N-1]/J[N-1];
  ftop = f[N-1];
  y[N-1] = f[N-1]*F[N-1]/(cl*gtop);
  for (int k=N-2; k>=0; --k) {
    const Real h = rg[k+1] - rg[k];
    const Real a0 = (3.0*f[k] - 1.0)/(f[k]*rg[k]);
    const Real a1 = (3.0*f[k+1] - 1.0)/(f[k+1]*rg[k+1]);
    const Real b = StarLogMean(chi[k], chi[k+1])*h*0.5*(F[k] + F[k+1])/cl;
    y[k] = (y[k+1]*(1.0 + 0.5*h*a1) + b)/(1.0 - 0.5*h*a0);
  }
  for (int k=0; k<N; ++k) E[k] = y[k]/f[k];
  jout = J;
}
}  // namespace

namespace {
//! the kinetic energy of a cell from the CONSERVED momenta.  On the cubed
//! sphere the transverse momenta are COVARIANT on the panel basis (angle cos c between
//! e_xi and e_eta), so KE = (m1^2 + (m2^2 + m3^2 - 2 c m2 m3)/(1 - c^2))/(2 d); the plain
//! sum of squares (the sp/Cartesian form, kept textually below) is wrong by the cross
//! term at every off-midline cell with transverse motion
KOKKOS_INLINE_FUNCTION
Real StarKeCs(const Real m1, const Real m2, const Real m3, const Real d, const Real c) {
  return CsKeMom(m1, m2, m3, c)/d;
}
void StarGravity(Mesh *pm, const Real bdt);
void StarWallFlux(Mesh *pm);
void StarBC(Mesh *pm);
void StarHist(HistoryData *pdata, Mesh *pm);
void StarApplyMltW(Mesh *pm, const Real w, const Real h);
void StarAdaptiveUpdate(Mesh *pm);
std::vector<char> StarAdaptiveRstWrite();
void StarFinal(ParameterInput *pin, Mesh *pm);
void StarIcPost(Mesh *pm);
void StarIcEquil(Mesh *pm);
void StarRemapLoad(Mesh *pm, const std::string &fn, const std::int64_t off,
                   const int *n);
void StarAtmPhiEff(const int ia, const int ie, const std::vector<Real> &a,
                   const std::vector<Real> &pc0, const std::vector<Real> &pf0,
                   const std::vector<Real> &p0c, const std::vector<Real> &p0f,
                   const std::vector<Real> &x1v, const std::vector<Real> &x1f,
                   std::vector<Real> &pc, std::vector<Real> &pf);
void StarAtmRadial(MeshBlockPack *pmbp, std::vector<Real> &pc0, std::vector<Real> &pf0,
                   std::vector<Real> &p0c, std::vector<Real> &p0f, std::vector<Real> &x1v,
                   std::vector<Real> &x1f, std::vector<Real> &ar0);
void StarAtmWrite(MeshBlockPack *pmbp, const int ia, const std::vector<Real> &a,
                  const std::vector<Real> &pc, const std::vector<Real> &pf);
void StarN2Compute(MeshBlockPack *pmbp, const EOS_Data &eos, const Real ar);
bool StarN2Band(MeshBlockPack *pmbp, const EOS_Data &eos, const Real ar, const Real tk,
                const bool hcs, const Real frac, const int nmin, const int ngap,
                const Real ppw, const int lcap, Real &rlo, Real &rhi, int &lmin,
                int &lmax);
Real StarN2Median(MeshBlockPack *pmbp, const EOS_Data &eos, const Real ar);
}  // namespace
// the L budget history of pgen/tests/rad_vet_tests2.cpp (am_hist)
void RadVetAmHistory(HistoryData *pdata, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief see the file header.

void StarStripBackground(ParameterInput *pin, Mesh *pm);   // (star_strip_bg.cpp)
void StarRestartPerturb(ParameterInput *pin, Mesh *pm);    // (star_strip_bg.cpp)
void StarEosReinitFromT(ParameterInput *pin, Mesh *pm);    // (star_strip_bg.cpp)

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto *pvet = pmbp->pradvet;
  star_pin_ = pin;
  // STAGE CS3 (m1-cs-implicit): the cubed sphere as well (x1 = r on both meshes; the
  // column, the potentials, the WB pair and the x1 BCs are radial).  On cs the solid
  // angle is 4 pi and the multi-mode seed is laid down in the global (theta, phi) of the
  // cell centre (see the seed below).
  const bool hcs = pmy_mesh_->use_cubed_sphere;
  if (pvet == nullptr || pmbp->phydro == nullptr ||
      !(pmy_mesh_->use_spherical_polar || hcs)) {
    StarFatal("needs <hydro>, <rad_vet> and mesh/use_spherical_polar or use_cubed_sphere",
              __LINE__);
  }
  if (pmbp->pmhd != nullptr) {
    StarFatal("<mhd> is not supported by this pgen (hydro only)", __LINE__);
  }
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto *ph = &fl;
  const Real cl = pvet->c_light, efl = pvet->e_floor, ar = pvet->arad;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto eos = ph->peos->eos_data;
  // kelvin per unit code temperature: the tabulated general EOS keeps its own scale
  // (EOS_Data::temp_cgs, mu_ref = 1); an ideal gas has T_code = p/d, whose scale is the
  // <units> temperature with <units>/mu folded in (EOS_Data::temp_cgs stays 1 there)
  const Real star_tk = eos.tbl.active ? eos.temp_cgs :
                       ((pmbp->punit != nullptr) ? pmbp->punit->temperature_cgs() : 1.0);
  const bool etg = ph->use_etotgrav;
  const bool wbdyn = ph->use_wellbalance_dynamic;
  {
    const std::string hg = pin->GetOrAddString("problem","star_gravity","wb");
    if (hg != "wb" && hg != "pointmass") {
      StarFatal("problem/star_gravity must be wb|pointmass", __LINE__);
    }
    star_pm_ = (hg == "pointmass");
  }
  if (star_pm_) {
    // plain (not well-balanced) hydro with the full VET force
    if (wbdyn || !etg) {
      StarFatal("problem/star_gravity = pointmass needs <hydro>/etotgrav = true and "
                "wellbalance_dynamic = false", __LINE__);
    }
    if (pvet->force_ref != radvet::VET_FREF_NONE) {
      StarFatal("problem/star_gravity = pointmass needs <rad_vet>/force_reference = none "
                "(no well-balanced pair to carry a reference force)", __LINE__);
    }
    if (pin->DoesParameterExist("problem","star_gm_column") &&
        pin->GetBoolean("problem","star_gm_column")) {
      StarFatal("problem/star_gravity = pointmass does not support star_gm_column",
                __LINE__);
    }
  } else {
    if (!(wbdyn && ph->use_wb_x1 && etg)) {
      StarFatal("needs <hydro>/etotgrav = true, wellbalance_dynamic = true and wb_x1 = "
                "true (the gravity source reads the x1 well-balanced cache)", __LINE__);
    }
    if (pvet->force_ref != radvet::VET_FREF_WB_ARAD) {
      StarFatal("needs <rad_vet>/force_reference = wb_arad (the well-balanced pair "
                "carries the reference force of the initial column)", __LINE__);
    }
  }
  if (pvet->opacity_type != radvet::VET_OPAC_TABLE) {
    StarFatal("needs <rad_vet>/opacity = table", __LINE__);
  }

  // ---- parameters
  star_gm_ = pin->GetReal("problem","star_gm");
  star_gmc_ = pin->DoesParameterExist("problem","star_gm_column") &&
              pin->GetBoolean("problem","star_gm_column");
  if (pin->DoesParameterExist("problem","star_base_heat")) {
    star_bheat_ = pin->GetReal("problem","star_base_heat");
    if (!(star_bheat_ > 0.0)) StarFatal("problem/star_base_heat must be > 0 [cm]",
                                        __LINE__);
  }
  star_rin_ = pmy_mesh_->mesh_size.x1min;
  const Real rtop = pmy_mesh_->mesh_size.x1max;
  star_fin_ = pin->GetReal("rad_vet","implicit_flux_x1min");
  const std::string fn = pin->GetString("problem","star_ic_file");
  const int nf = pin->GetOrAddInteger("problem","star_nfine",16384);
  const Real tau_int = pin->GetOrAddReal("problem","star_tau_int",1.0);
  // star_ic_cols = 5: the file also carries the radiation energy density E (a state that
  // was relaxed by a run); 4 (default): E = a T(rho,eint)^4
  const int ncols = pin->GetOrAddInteger("problem","star_ic_cols",4);
  if (ncols != 4 && ncols != 5) StarFatal("problem/star_ic_cols must be 4 or 5",
                                          __LINE__);
  // problem/mlt_flux_frozen (default false; in 3-D only as the start-up scaffold that is
  // ramped off by problem/mlt_ramp_start/_time while convection grows, or relaxed by
  // problem/mlt_closure = adaptive): the MLT flux of the IC, F_MLT = F_r fmlt/(1 - fmlt)
  // with fmlt = F_MLT/F the file's column 7 (make_ic mlt: r rho eint F_r E T fmlt), so
  // that F_r + F_MLT = L/(4 pi r^2) exactly, is deposited as the conservative energy
  // source -(A F_MLT|_{i+1/2} - A F_MLT|_{i-1/2})/V into the RADIATION energy through
  // <rad_vet>/esrc (the implicit solve's old vector, every stage with its weight).
  // Rebuilt from the file on every start (restarts).
  star_mlt_ = pin->GetOrAddBoolean("problem","mlt_flux_frozen",false);
  if (star_mlt_ && ncols != 5) {
    StarFatal("problem/mlt_flux_frozen = true needs star_ic_cols = 5 and the 7-column "
              "IC file (columns 6, 7: T, F_MLT/F)", __LINE__);
  }
  // ONE MeshBlock along x1 (the margins below, the history's face picks and the balanced
  // IC read the radial grid of block 0)
  if (pmy_mesh_->mesh_indcs.nx1 != indcs.nx1) {
    StarFatal("needs ONE MeshBlock along x1 (meshblock/nx1 = mesh/nx1)", __LINE__);
  }
  // the fine grid spans the OUTERMOST GHOST FACES (the code's own x1f, stretched or
  // not) plus two end-cell widths on each side: on a uniform grid this is the old
  // (ng + 2) dx margin exactly; on a coarse stretched grid the old "at least the mean
  // width" margin reached rho < 1e-14 above the top (v2/runs/P03 fatal)
  Real rlo_g = star_rin_, rhi_g = rtop, dxlo = 0.0, dxhi = 0.0;
  {
    auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
    rlo_g = hx1f(0,0);
    rhi_g = hx1f(0,n1m1+1);
    dxlo = hx1f(0,indcs.is+1) - hx1f(0,indcs.is);
    dxhi = hx1f(0,indcs.ie+1) - hx1f(0,indcs.ie);
    // problem/star_grid_dump (default false): rank 0 writes the code's own radial grid of
    // block 0, "i x1f(i) x1v(i) x1f(i+1)-x1f(i)" for every cell incl. ghosts (i = 0 is
    // the first ghost; active cells is..ie), to <basename>.x1grid.txt
    if (pin->DoesParameterExist("problem","star_grid_dump") &&
        pin->GetBoolean("problem","star_grid_dump") &&
        global_variable::my_rank == 0) {
      auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
      const std::string gf = pin->GetString("job","basename") + ".x1grid.txt";
      FILE *fp = fopen(gf.c_str(), "w");
      if (fp != nullptr) {
        fprintf(fp, "# is = %d ie = %d; i x1f x1v dx1\n", indcs.is, indcs.ie);
        for (int i=0; i<=n1m1; ++i) {
          fprintf(fp, "%d %.17e %.17e %.17e\n", i, hx1f(0,i), hx1v(0,i),
                  hx1f(0,i+1) - hx1f(0,i));
        }
        fprintf(fp, "%d %.17e\n", n1m1 + 1, hx1f(0,n1m1+1));
        fclose(fp);
      }
    }
  }
  star_rlo_ = rlo_g - 2.0*dxlo;
  const Real rhi = rhi_g + 2.0*dxhi;
  star_nf_ = nf;
  star_dr_ = (rhi - star_rlo_)/(nf - 1);
  std::vector<Real> hr(nf), hd(nf), he(nf), hF(nf), hE(nf, -1.0), hM(nf, 0.0);
  for (int n=0; n<nf; ++n) hr[n] = star_rlo_ + n*star_dr_;

  // problem/star_wind_ic (he-wind-bc, default false): above the radius r_c where the
  // column's density falls below that of a beta-law wind (star_wind_mdot [g/s, full
  // sphere], star_wind_vinf, star_wind_beta, star_wind_v0 [cm/s], star_wind_r0 [cm]: v =
  // max(v0, vinf (1 - r0/r)^beta), rho = Mdot/(4 pi r^2 v)), the IC is that wind: F_r =
  // F_c (r_c/ r)^2, E = (F_r/c) (1 + (c E_c/F_c - 1)(r_c/r)^2) (free streaming far out),
  // T_gas = T_c (E/E_c)^(1/4) (continuous at r_c), v_r = the beta law.  The file then
  // need not cover the mesh top (it is extrapolated in log up to r_c, which must lie
  // inside it).
  const bool wic = pin->GetOrAddBoolean("problem","star_wind_ic",false);
  std::vector<Real> hV(nf, 0.0), hdc, hec;
  Real wrc = 0.0;
  // problem/star_ic_eint_from_t (default false = bitwise the old IC): set the internal
  // energy of the column from (rho, T) with the code's OWN EOS, T = the file's column 6
  // [K] (the 7-column make_ic files), instead of taking the file's eint (column 3).  The
  // file's eint was computed by the IC builder with its own interpolant of the EOS (for
  // the general EOS: a cubic spline of the code's table dump, not the code's bicubic
  // Hermite patch), so T(rho, eint) in the code missed the column T by the difference of
  // the two interpolants (agcar/geos B: 5.8e-6, vs 4.6e-7 for the ideal gas); E = aT^4 of
  // the file is then out of equilibrium with the gas by 4x that.  With this key the
  // code's T equals the column T to round-off (the log-linear resampling aside).
  // Default true on a fresh run whose file carries the T
  // column (star_ic_cols = 5 and >= 7 columns on its first data line); otherwise
  // (4-column files, a restart whose file lacks the key) the old default false.
  // Recorded.
  bool eft_def = false;
  if (!global_variable::restart_run && ncols == 5) {
    std::ifstream f0(fn);
    std::string l0;
    while (std::getline(f0, l0)) {
      if (l0.empty() || l0[0] == '#') continue;
      std::istringstream s0(l0);
      int nc = 0;
      Real v0;
      while (s0 >> v0) ++nc;
      eft_def = (nc >= 7);
      break;
    }
  }
  const bool eft = pin->GetOrAddBoolean("problem","star_ic_eint_from_t",eft_def);
  if (eft && ncols != 5) {
    StarFatal("problem/star_ic_eint_from_t = true needs star_ic_cols = 5 and the "
              "7-column file (column 6 = T [K])", __LINE__);
  }
  bool tcol_read = false;
  std::vector<Real> hTcol(nf, -1.0);   // the file's T column [K] on the fine grid
  std::vector<Real> hfr_keep, hfM_keep;  // star_ic_rad = sph1d: the file's r, fmlt (dump)
  // ---- the column file: r rho eint F_r, ascending r, resampled in log
  {
    std::ifstream f(fn);
    if (!f.good()) StarFatal("cannot open problem/star_ic_file '" + fn + "'", __LINE__);
    std::vector<Real> fr, fd, fe, fF, fE, fM, fT;
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      Real a, b, c, d, e = -1.0;
      if (!(ss >> a >> b >> c >> d)) continue;
      if (ncols == 5 && !(ss >> e)) continue;
      Real tcol = 0.0, fm = 0.0;
      if ((star_mlt_ || eft) && !(ss >> tcol >> fm)) {
        StarFatal("problem/mlt_flux_frozen / star_ic_eint_from_t: star_ic_file line "
                  "without columns 6-7 (T, fmlt)", __LINE__);
      }
      fM.push_back(fm);
      fT.push_back(tcol);
      fr.push_back(a); fd.push_back(b); fe.push_back(c); fF.push_back(d);
      fE.push_back(e);
    }
    if (fr.size() < 2 || fr.front() > hr.front() || (!wic && fr.back() < hr.back())) {
      std::ostringstream os;
      os << "problem/star_ic_file '" << fn << "' must cover r = " << hr.front() << " .. "
         << hr.back() << " (the mesh plus the ghost margin), it has "
         << (fr.empty() ? 0.0 : fr.front()) << " .. " << (fr.empty() ? 0.0 : fr.back());
      StarFatal(os.str(), __LINE__);
    }
    std::size_t kk = 0;
    for (int n=0; n<nf; ++n) {
      const Real r = hr[n];
      while (kk + 2 < fr.size() && fr[kk+1] < r) ++kk;
      const Real w = (r - fr[kk])/(fr[kk+1] - fr[kk]);
      hd[n] = exp((1.0 - w)*log(fd[kk]) + w*log(fd[kk+1]));
      he[n] = exp((1.0 - w)*log(fe[kk]) + w*log(fe[kk+1]));
      hF[n] = exp((1.0 - w)*log(fF[kk]) + w*log(fF[kk+1]));
      if (ncols == 5) hE[n] = exp((1.0 - w)*log(fE[kk]) + w*log(fE[kk+1]));
      if (star_mlt_) hM[n] = (1.0 - w)*fM[kk] + w*fM[kk+1];
      if (star_mlt_ || eft) hTcol[n] = exp((1.0 - w)*log(fT[kk]) + w*log(fT[kk+1]));
    }
    tcol_read = (star_mlt_ || eft);
    hfr_keep = fr;
    hfM_keep = fM;
    if (eft) {
      for (int n=0; n<nf; ++n) {
        Real ee, pp, cr, ct, cv;
        eos.ThermoAt(hd[n], hTcol[n]/star_tk, ee, pp, cr, ct, cv);
        he[n] = ee;
      }
    }
    if (wic) {
      const Real mdot = pin->GetReal("problem","star_wind_mdot");
      const Real vinf = pin->GetReal("problem","star_wind_vinf");
      const Real wbeta = pin->GetOrAddReal("problem","star_wind_beta",1.0);
      const Real wv0 = pin->GetOrAddReal("problem","star_wind_v0",1.0e6);
      const Real wr0 = pin->GetReal("problem","star_wind_r0");
      auto vw = [&](const Real r) {
        const Real v = (r > wr0) ? vinf*pow(1.0 - wr0/r, wbeta) : 0.0;
        return std::max(v, wv0);
      };
      int nc = -1;
      for (int n=0; n<nf; ++n) {
        if (hr[n] >= wr0 && hd[n] <= mdot/(4.0*M_PI*SQR(hr[n])*vw(hr[n]))) {
          nc = n;
          break;
        }
      }
      if (nc < 1 || hr[nc] > fr.back()) {
        StarFatal("problem/star_wind_ic: the wind density never exceeds the column's "
                  "inside the file's range above star_wind_r0", __LINE__);
      }
      wrc = hr[nc];
      hdc = hd;   // the column itself (extrapolated), for the balance anchor
      hec = he;
      const Real tc = eos.Temperature(hd[nc], he[nc]);
      const Real ec = (hE[nc] > 0.0) ? hE[nc] : ar*SQR(SQR(tc));
      const Real fc = hF[nc];
      for (int n=nc; n<nf; ++n) {
        const Real r = hr[n], x2 = SQR(wrc/r);
        const Real v = vw(r);
        const Real d = mdot/(4.0*M_PI*SQR(r)*v);
        const Real fr_ = fc*x2;
        const Real er = (fr_/cl)*(1.0 + (cl*ec/fc - 1.0)*x2);
        const Real t = tc*pow(er/ec, 0.25);
        Real e, p, cr, ct, cv;
        eos.ThermoAt(d, t, e, p, cr, ct, cv);
        hd[n] = d;
        he[n] = e;
        hF[n] = fr_;
        hE[n] = er;
        hM[n] = 0.0;
        hV[n] = v;
      }
      if (global_variable::my_rank == 0) {
        std::cout << "star_vet: star_wind_ic: Mdot = " << mdot << " g/s, vinf = " << vinf
                  << ", beta = " << wbeta << ", r0 = " << wr0 << "; the wind starts at "
                  << "r_c = " << wrc << " (rho " << hd[nc] << ", T [K] "
                  << tc*star_tk << ", v "
                  << hV[nc] << ")" << std::endl;
      }
    }
  }

  // ---- the Rosseland + Planck tables (one shared grid), handed to the module
  Real tab_lt_lo = 0.0, tab_lt_hi = 0.0, tab_ld_lo = 0.0, tab_ld_hi = 0.0;
  {
    const std::string rt = pin->GetString("problem","star_opac_table");
    const std::string pt = pin->GetString("problem","star_planck_table");
    DvceArray2D<Real> krt, kpt;
    DvceArray1D<Real> mlT, mlD, plT, plD;
    int mnT = 0, mnD = 0, pnT = 0, pnD = 0;
    Real a1, a2, a3, a4, b1, b2, b3, b4;
    // problem/star_opac_logd_min (default none): extend both tables to this log10 rho
    const Real ldx = pin->DoesParameterExist("problem","star_opac_logd_min") ?
                     pin->GetReal("problem","star_opac_logd_min") : 1.0e300;
    StarReadOpacityTable(rt, krt, mlT, mlD, mnT, mnD, tab_lt_lo, tab_lt_hi, tab_ld_lo,
                         tab_ld_hi, ldx);
    StarReadOpacityTable(pt, kpt, plT, plD, pnT, pnD, a1, a2, a3, a4, ldx);
    b1 = tab_lt_lo; b2 = tab_lt_hi; b3 = tab_ld_lo; b4 = tab_ld_hi;
    if (mnT != pnT || mnD != pnD || fabs(a1 - b1) + fabs(a2 - b2) + fabs(a3 - b3) +
        fabs(a4 - b4) > 1.0e-9) {
      StarFatal("the Rosseland and Planck tables are not on the same grid", __LINE__);
    }
    pvet->SetOpacityTables(krt, kpt, mlT, mlD, mnT, mnD);
    // the lookup's units against the EOS's own code temperature and <units> density
    const Real tcgs = star_tk;
    const Real dcgs = (pmbp->punit != nullptr) ? pmbp->punit->density_cgs() : 1.0;
    if (std::fabs(pvet->otab.tunit/tcgs - 1.0) > 1.0e-5 ||
        std::fabs(pvet->otab.dunit/dcgs - 1.0) > 1.0e-5) {
      std::ostringstream os;
      os << "<rad_vet>/temp_unit_kelvin = " << pvet->otab.tunit << ", rho_unit_cgs = "
         << pvet->otab.dunit << " but this run's units are (" << tcgs << ", " << dcgs
         << ")";
      StarFatal(os.str(), __LINE__);
    }
  }

  // ---- problem/star_ic_rad (file = the old IC; the default is set below).
  // sph1d: above the depth tau_R = star_ic_rad_taub (from the mesh top, the run's
  // own kappa_t), the radiation energy density of the column is replaced by the solution
  // of the 1-D spherical grey transfer problem in radiative equilibrium (S = J = E, the
  // the rays' source in RE, the vet_scatter split included) with the column's own
  // radiative flux F_r: a VET iteration of a p-z formal solution (f = K/J, g_top = H/J)
  // and the radial moment equation (StarSph1dPass).  The file's E above R_ph is an
  // isothermal-atmosphere guess (2.1-2.6x the relaxed E at 1.0-1.2 R_ph for AG Car A).
  // E = w E_sph1d + (1 - w) E_col, w a smoothstep LINEAR in tau, 1 at star_ic_rad_taua, 0
  // at star_ic_rad_taub.  star_ic_rad_gas = column (old default): T of the column; re:
  // T^4
  // = E at fixed rho (blended with the same w; eint from the code's EOS); re_atm: as re
  // but the floor gas (rho at the column's floor) keeps the column's T, so its a_ref = kt
  // F/c is the file IC's (with re its cold floor gas falls 3x faster and the implicit
  // solve then takes ~6 instead of ~4 Picard iterations, AG Car A).  star_ic_rad_hse
  // (needs gas = re | re_atm): rho re-integrated outward from r_b at fixed T, dP_gas/dr
  // = -rho (g - kt F/c) (the continuum balance of the WB pair, whose a_ref = kt F/c is
  // built below from this column), floored at the column's floor; E solve and HSE
  // iterated together.
  // F_r and the MLT fraction are not changed.  Deterministic (host solve, per-point
  // device opacities) and rebuilt on every start, as the column itself.
  // star_ic_rad_dump = <file>: rank 0 writes the new column in the IC file format (r rho
  // eint F_r E T fmlt; usable with star_ic_rad = file; none = no dump).  NOTE-icgen.md.
  // DEFAULT sph1d (with star_ic_rad_gas = re_atm, star_ic_rad_hse = true) on
  // a FRESH start of every star_vet run where it applies: opacity = table, no
  // star_wind_ic, star_ic_rad_thin not set.  Elsewhere the default stays file (the old
  // IC), silently.  A restart whose file lacks the key keeps file (the column is rebuilt
  // on every start, so this matters for the restart too).  The resolved value is echoed
  // (later restarts keep it); explicit input overrides.  To keep the old IC on a fresh
  // start name problem/star_ic_rad = file (and rad_vet/ic_rad = none for the old
  // radiation start).
  // (star_ic_rad_atm_force = module): the factor phi = a_module/a_continuum of
  // the radiative acceleration of each mesh cell of the thin atmosphere (hatm_c, index
  // = the x1 cell index of block 0 incl. ghosts) and on the fine grid (hatm_f); 1 (and
  // unused) elsewhere / when off
  bool hatm_on = false;
  std::vector<Real> hatm_c, hatm_f;
  const bool icrdef = !global_variable::restart_run &&
      (pvet->opacity_type == radvet::VET_OPAC_TABLE) && !wic &&
      !(pin->DoesParameterExist("problem","star_ic_rad_thin") &&
        pin->GetBoolean("problem","star_ic_rad_thin"));
  const std::string icr0 = pin->GetOrAddString("problem","star_ic_rad",
                                               icrdef ? "sph1d" : "file");
  if (icr0 != "file") {
    const std::string icr = icr0;
    if (icr != "sph1d") StarFatal("problem/star_ic_rad must be file | sph1d", __LINE__);
    if (wic) {
      StarFatal("problem/star_ic_rad = sph1d does not combine with star_wind_ic",
                __LINE__);
    }
    if (pin->DoesParameterExist("problem","star_ic_rad_thin") &&
        pin->GetBoolean("problem","star_ic_rad_thin")) {
      StarFatal("problem/star_ic_rad = sph1d replaces star_ic_rad_thin (set one)",
                __LINE__);
    }
    Kokkos::Timer tw0;
    const Real taua = pin->GetOrAddReal("problem","star_ic_rad_taua",30.0);
    const Real taub = pin->GetOrAddReal("problem","star_ic_rad_taub",3000.0);
    const Real dtau0 = pin->GetOrAddReal("problem","star_ic_rad_dtau",0.05);
    const Real dlnr = pin->GetOrAddReal("problem","star_ic_rad_dlnr",2.0e-3);
    const int ncore = pin->GetOrAddInteger("problem","star_ic_rad_ncore",64);
    // fresh-start defaults re_atm + hse (the recommended IC, NOTE-icgen);
    // a restart whose file lacks the keys keeps column / false
    const bool icgdef = !global_variable::restart_run;
    const std::string gasm = pin->GetOrAddString("problem","star_ic_rad_gas",
                                                 icgdef ? "re_atm" : "column");
    const bool hse = pin->GetOrAddBoolean("problem","star_ic_rad_hse",icgdef);
    constexpr int maxit = 60;
    const Real tol = pin->GetOrAddReal("problem","star_ic_rad_tol",1.0e-8);
    const Real htol = pin->GetOrAddReal("problem","star_ic_rad_hse_tol",1.0e-6);
    // HSE outer loop on a near-Eddington column:
    // (1) Gamma = kt F/(c g) of the HSE march is clamped at star_ic_rad_hse_gmax (as
    //     star_bc_hse_gmax).  Where Gamma > 1 the fixed-T march makes rho grow outward,
    //     kappa_t(rho) grows with it and the E-solve <-> HSE loop runs away.  Clamped
    //     cells keep a net outward force (Gamma - gmax) g.
    // (2) At the first rise of the outer residual max |dln rho| the update is damped,
    //     ln rho <- (1 - w) ln rho + w ln rho_HSE, w = star_ic_rad_hse_damp, from then
    // on. (3) The outer loop has its own limit star_ic_rad_hse_maxit and warns when it
    //     does not converge.
    // A column with Gamma <= gmax and a monotone residual is not affected.
    const Real hgmax = pin->GetOrAddReal("problem","star_ic_rad_hse_gmax",0.9);
    const Real hsdamp = pin->GetOrAddReal("problem","star_ic_rad_hse_damp",0.5);
    const int hsmaxo = pin->GetOrAddInteger("problem","star_ic_rad_hse_maxit",400);
    if (!(taua > 0.0 && taub > taua && dtau0 > 0.0 && dlnr > 0.0 && ncore >= 4 &&
          tol > 0.0 && htol > 0.0 && hgmax > 0.0 && hgmax < 1.0 &&
          hsdamp > 0.0 && hsdamp <= 1.0 && hsmaxo >= 1)) {
      StarFatal("problem/star_ic_rad: needs 0 < taua < taub, dtau, dlnr, tol > 0, ncore "
                ">= 4, 0 < hse_gmax < 1, 0 < hse_damp <= 1, hse_maxit >= 1",
                __LINE__);
    }
    if (gasm != "column" && gasm != "re" && gasm != "re_atm") {
      StarFatal("problem/star_ic_rad_gas must be column | re | re_atm", __LINE__);
    }
    // re_atm: radiative equilibrium only where rho is above the column's floor (the floor
    // gas keeps the column's T, hence its a_ref = kt F/c, as with star_ic_rad = file)
    const bool gre = (gasm == "re" || gasm == "re_atm");
    const bool gatm = (gasm == "re_atm");
    if (hse && !gre) {
      StarFatal("problem/star_ic_rad_hse = true needs star_ic_rad_gas = re | re_atm",
                __LINE__);
    }
    // the column as read: T [code] (the file's T with star_ic_eint_from_t, else the
    // code's T(rho, eint)), E (the file's, else a T^4), rho
    std::vector<Real> d0(hd), t0(nf), e0(nf);
    for (int n=0; n<nf; ++n) {
      t0[n] = eft ? hTcol[n]/star_tk : eos.Temperature(hd[n], he[n]);
      e0[n] = (hE[n] > 0.0) ? hE[n] : ar*SQR(SQR(t0[n]));
    }
    std::vector<Real> tcur(t0), ef(e0), kt, tauf(nf, 0.0);
    radvet::VetOpacTab otr = pvet->otab;
    auto ftau = [&]() {   // kappa_t and tau from the mesh top (as the pgen's tau below)
      StarKtEval(otr, hd, tcur, kt);
      tauf[nf-1] = 0.0;
      for (int n=nf-2; n>=0; --n) {
        const Real a0 = (hr[n] < rtop) ? hd[n]*kt[n] : 0.0;
        const Real a1 = (hr[n+1] < rtop) ? hd[n+1]*kt[n+1] : 0.0;
        tauf[n] = tauf[n+1] + 0.5*star_dr_*(a0 + a1);
      }
    };
    auto frph = [&]() {   // tau_R = 2/3
      for (int n=nf-2; n>=0; --n) {
        if (tauf[n] >= 2.0/3.0 && tauf[n+1] < 2.0/3.0) {
          return hr[n] + star_dr_*(tauf[n] - 2.0/3.0)/(tauf[n] - tauf[n+1]);
        }
      }
      return -1.0;
    };
    auto flog = [&](const std::vector<Real> &a, const Real r) {   // fine grid, log-linear
      Real x = (r - star_rlo_)/star_dr_;
      int i = static_cast<int>(floor(x));
      i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
      Real w = x - i;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return exp((1.0 - w)*log(a[i]) + w*log(a[i+1]));
    };
    auto glog = [](const std::vector<Real> &x, const std::vector<Real> &y, const Real r) {
      const int n = static_cast<int>(x.size());
      int i = static_cast<int>(std::upper_bound(x.begin(), x.end(), r) - x.begin()) - 1;
      i = (i < 0) ? 0 : ((i > n - 2) ? (n - 2) : i);
      Real w = (r - x[i])/(x[i+1] - x[i]);
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return exp((1.0 - w)*log(y[i]) + w*log(y[i+1]));
    };
    ftau();
    const Real rph0 = frph();
    int nb = 0;
    for (int n=nf-1; n>=0; --n) {
      if (tauf[n] >= taub) {
        nb = n;
        break;
      }
    }
    if (tauf[0] < taub && global_variable::my_rank == 0) {
      std::cout << "star_vet: star_ic_rad = sph1d: WARNING the column's bottom tau ("
                << tauf[0] << ") is below star_ic_rad_taub; r_b = the fine-grid bottom"
                << std::endl;
    }
    Real dfl = 1.0e300;   // the column's floor: its minimum density at and above r_b
    for (int n=nb; n<nf; ++n) dfl = std::min(dfl, d0[n]);
    // gravity of the HSE march: the point mass, or (star_gm_column) the static monopole
    // of the current column, as built below
    Real gcode = 0.0;
    if (star_gmc_) {
      if (pmbp->punit == nullptr) StarFatal("problem/star_gm_column needs <units>",
                                            __LINE__);
      const Real lu = pmbp->punit->length_cgs(), tu = pmbp->punit->time_cgs();
      gcode = units::Units::grav_constant_cgs*pmbp->punit->mass_cgs()*tu*tu/(lu*lu*lu);
    }
    std::vector<Real> rg, dg, tg, fg, eg, cg, jg, kg;
    Real gtop = 0.0, ftop = 0.0, dlr = 0.0, rmx = 0.0;
    int nit = 0, nout = 0, nflo = -1, ncl = 0, ngcl = -1;
    Real hsom = 1.0, dlrp = 1.0e300;
    int hsodmp = -1;
    for (int outer=0; outer<(hse ? hsmaxo : 1); ++outer) {
      ++nout;
      // the solve grid: uniform in s = int dr/(dlnr r) + dtau/(dtau0 max(1, tau)^0.7)
      std::vector<Real> sm(nf, 0.0);
      for (int n=nb+1; n<nf; ++n) {
        const Real a0 = (hr[n-1] < rtop) ? hd[n-1]*kt[n-1] : 0.0;
        const Real a1 = (hr[n] < rtop) ? hd[n]*kt[n] : 0.0;
        sm[n] = sm[n-1] + star_dr_*(1.0/(dlnr*hr[n]) +
                             0.5*(a0 + a1)/(dtau0*pow(std::max(1.0, tauf[n]), 0.7)));
      }
      const int ngr = std::min(std::max(16, static_cast<int>(ceil(sm[nf-1])) + 1), 50000);
      rg.assign(ngr, 0.0);
      {
        int n = nb + 1;
        for (int g=0; g<ngr; ++g) {
          const Real s = sm[nf-1]*g/(ngr - 1);
          while (n < nf - 1 && sm[n] < s) ++n;
          rg[g] = hr[n-1] + star_dr_*(s - sm[n-1])/(sm[n] - sm[n-1]);
        }
        rg[0] = hr[nb];
        rg[ngr-1] = hr[nf-1];
      }
      dg.resize(ngr); tg.resize(ngr); fg.resize(ngr); eg.resize(ngr); cg.resize(ngr);
      for (int g=0; g<ngr; ++g) {
        dg[g] = flog(hd, rg[g]);
        tg[g] = flog(tcur, rg[g]);
        fg[g] = flog(hF, rg[g]);
        eg[g] = flog(ef, rg[g]);
      }
      for (int it=0; it<maxit; ++it) {
        if (gre) {
          for (int g=0; g<ngr; ++g) {
            tg[g] = (gatm && dg[g] <= dfl*(1.0 + 1.0e-6)) ? flog(t0, rg[g])
                                                         : pow(eg[g]/ar, 0.25);
          }
        }
        StarKtEval(otr, dg, tg, kg);
        for (int g=0; g<ngr; ++g) cg[g] = (rg[g] < rtop) ? dg[g]*kg[g] : 0.0;
        const std::vector<Real> ep(eg);
        StarSph1dPass(rg, cg, fg, ncore, cl, eg, jg, gtop, ftop);
        ++nit;
        Real de = 0.0;
        for (int g=0; g<ngr; ++g) de = std::max(de, std::fabs(eg[g]/ep[g] - 1.0));
        if (!(eg[0] > 0.0) || !std::isfinite(de)) {
          StarFatal("problem/star_ic_rad = sph1d: the VET iteration failed", __LINE__);
        }
        if (de < tol) break;
      }
      // the blended column on the fine grid
      for (int n=nb; n<nf; ++n) {
        const Real x = (taub - tauf[n])/(taub - taua);
        const Real w = (x >= 1.0) ? 1.0 : ((x <= 0.0) ? 0.0 : x*x*(3.0 - 2.0*x));
        ef[n] = w*glog(rg, eg, hr[n]) + (1.0 - w)*e0[n];
        if (gre) tcur[n] = w*pow(ef[n]/ar, 0.25) + (1.0 - w)*t0[n];
        if (gatm && hd[n] <= dfl*(1.0 + 1.0e-9)) tcur[n] = t0[n];
      }
      if (!hse) {
        if (gre) ftau();
        break;
      }
      // hydrostatic balance at fixed T, outward from r_b (kappa_t at the previous rho):
      // P(rho_{n+1}) + dr/2 rho_{n+1} g_{n+1} = P(rho_n) - dr/2 rho_n g_n
      StarKtEval(otr, hd, tcur, kt);
      std::vector<Real> gm(nf, star_gm_);
      if (star_gmc_) {
        std::vector<Real> cm(nf, 0.0);
        for (int n=1; n<nf; ++n) {
          cm[n] = cm[n-1] + 0.5*star_dr_*4.0*M_PI*(hr[n-1]*hr[n-1]*hd[n-1] +
                                                   hr[n]*hr[n]*hd[n]);
        }
        const Real x = (star_rin_ - star_rlo_)/star_dr_;
        const int i0 = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
        const Real cm0 = (1.0 - (x - i0))*cm[i0] + (x - i0)*cm[i0+1];
        for (int n=0; n<nf; ++n) gm[n] = star_gm_ + gcode*(cm[n] - cm0);
      }
      std::vector<Real> dn(hd);
      bool flo = false;
      nflo = -1;
      ncl = 0;
      ngcl = -1;
      rmx = 0.0;
      for (int n=nb; n<nf-1; ++n) {
        if (flo) {
          dn[n+1] = dfl;
          continue;
        }
        Real g0 = gm[n]/SQR(hr[n]) - kt[n]*hF[n]/cl;
        Real g1 = gm[n+1]/SQR(hr[n+1]) - kt[n+1]*hF[n+1]/cl;
        if (kt[n]*hF[n]/cl > hgmax*gm[n]/SQR(hr[n])) {
          g0 = (1.0 - hgmax)*gm[n]/SQR(hr[n]);
        }
        if (kt[n+1]*hF[n+1]/cl > hgmax*gm[n+1]/SQR(hr[n+1])) {
          g1 = (1.0 - hgmax)*gm[n+1]/SQR(hr[n+1]);
          if (n + 1 > ngcl) ngcl = n + 1;
          ++ncl;
        }
        Real ee, p0, cr, ct, cv;
        eos.ThermoAt(dn[n], tcur[n], ee, p0, cr, ct, cv);
        const Real rhs = p0 - 0.5*star_dr_*dn[n]*g0;
        auto res = [&](const Real lx) {
          Real e1, p1, c1, c2, c3;
          eos.ThermoAt(exp(lx), tcur[n+1], e1, p1, c1, c2, c3);
          return (p1 + 0.5*star_dr_*exp(lx)*g1)/rhs - 1.0;
        };
        Real x0 = log(dn[n]) + log(hd[n+1]/hd[n]), x1 = x0 + 1.0e-6;
        Real f0 = res(x0), f1 = res(x1);
        for (int it=0; it<60 && fabs(f1) > 1.0e-14 && f1 != f0; ++it) {
          const Real x2 = x1 - f1*(x1 - x0)/(f1 - f0);
          x0 = x1; f0 = f1; x1 = x2;
          f1 = res(x1);
        }
        rmx = std::max(rmx, fabs(f1));
        if (!(rhs > 0.0) || !std::isfinite(x1) || !(exp(x1) > dfl)) {
          flo = true;
          nflo = n + 1;
          dn[n+1] = dfl;
        } else {
          dn[n+1] = exp(x1);
        }
      }
      dlr = 0.0;
      for (int n=nb; n<nf; ++n) dlr = std::max(dlr, std::fabs(log(dn[n]/hd[n])));
      // the outer loop oscillates (the first rise of the residual): damp from now on
      if (outer > 0 && dlr > dlrp && hsom == 1.0) {
        hsom = hsdamp;
        hsodmp = outer;
      }
      dlrp = dlr;
      if (hsom < 1.0) {
        for (int n=nb; n<nf; ++n) {
          dn[n] = exp((1.0 - hsom)*log(hd[n]) + hsom*log(dn[n]));
        }
      }
      hd = dn;
      ftau();
      if (dlr < htol) break;
    }
    // problem/star_ic_rad_atm_force = module | continuum (needs hse and gas =
    // re_atm; default continuum = bitwise the old IC; recorded).  In the optically thin
    // atmosphere above R_ph the flux is fixed by the luminosity (F0 = L/(4 pi r^2) to
    // 3e-4 after ic_rad, measured) and <rad_vet> gives cell i the force
    // [(rho k)_{i-1/2} F_{i-1/2} + (rho k)_{i+1/2} F_{i+1/2}]/(2c), (rho k)_face the
    // distance-weighted mean of the two cells (the matched deposit, and the steady-state
    // conservative one), not rho_i kt_i F/c: with ~2 cells per density scale height the
    // face means exceed the cell value by phi = (1 + cosh(dln(rho k)))/2 ~ 1.06 (AG Car
    // A), an uncompensated outward acceleration (phi - 1) Gamma g ~ 1-4 % g that lifted
    // the atmosphere (NOTE-vimpbc 4a).  Here, after the HSE loop, rho is re-marched
    // outward from the lower face of the first mesh cell above R_ph (tau_R 2/3 of the
    // converged column) with the acceleration phi(r) kt F/c, phi from the code's own mesh
    // cells (piecewise linear between their centres) at the current rho, T held at the
    // column's (RE where rho is above the floor, the column's in the floor, as re_atm),
    // iterated to a fixed point; a_ref of those cells (and Phi_eff above that face) take
    // the same phi, so the WB pair and the module's force agree.  E, F and everything at
    // and below that face (fine grid, a_ref, Phi_eff, cells) are bitwise unchanged; the
    // change reaches below R_ph only through <rad_vet>/ic_rad's relaxation (|de/e| <=
    // 5e-4 at 0.99 R_ph, <= 1e-6 below 0.95 R_ph at cycle 0, AG Car A), hence opt-in.
    // See also star_ic_rad_atm_post (StarIcPost).
    // star_ic_rad_atm_relax w (fresh-start default 0.5; a restart whose file lacks the
    // key keeps 1 = the old march, bitwise): log phi under-relaxed by w, and every face
    // solve of the march starts from the node below.  The old guess, the last pass's rho
    // at the node, is the floor once a node fell to it; the secant then misses the root
    // there, the node stays on the floor: a one-cell rho cliff with phi ~ 20-60 whose
    // cell alternates between passes (not converged, the result
    // roundoff-dependent).  With w < 1, not converged in
    // star_ic_rad_atm_maxit passes (fresh default 1000) is fatal.
    {
      const std::string afm = pin->GetOrAddString("problem","star_ic_rad_atm_force",
                                                   "continuum");
      if (afm != "module" && afm != "continuum") {
        StarFatal("problem/star_ic_rad_atm_force must be module | continuum", __LINE__);
      }
      if (afm == "module" && !(hse && gatm)) {
        StarFatal("problem/star_ic_rad_atm_force = module needs star_ic_rad_hse = true "
                  "and star_ic_rad_gas = re_atm", __LINE__);
      }
      if (afm == "module") {
        const bool rst = global_variable::restart_run;
        const int amaxit = pin->GetOrAddInteger("problem","star_ic_rad_atm_maxit",
                                                rst ? 200 : 1000);
        const Real arel = pin->GetOrAddReal("problem","star_ic_rad_atm_relax",
                                            rst ? 1.0 : 0.5);
        if (amaxit < 1 || !(arel > 0.0 && arel <= 1.0)) {
          StarFatal("problem/star_ic_rad_atm_maxit >= 1 and 0 < star_ic_rad_atm_relax <= "
                    "1", __LINE__);
        }
        const Real atol = pin->GetOrAddReal("problem","star_ic_rad_atm_tol",1.0e-12);
        auto ax1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                        pmbp->pcoord->x1v);
        auto ax1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(),
                                                        pmbp->pcoord->xx1f);
        const int n1 = n1m1 + 1, aie = indcs.ie;
        const Real rpa = frph();
        int ia = -1;
        for (int i=0; i<n1; ++i) {
          if (ax1v(0,i) >= rpa) {
            ia = i;
            break;
          }
        }
        if (!(rpa > 0.0) || ia < 1 || ia > aie) {
          StarFatal("problem/star_ic_rad_atm_force = module: R_ph not inside the mesh",
                    __LINE__);
        }
        const int n0 = static_cast<int>(floor((ax1f(0,ia) - star_rlo_)/star_dr_));
        auto lin = [&](const std::vector<Real> &a, const Real r) {   // as StarLinInterp
          Real x = (r - star_rlo_)/star_dr_;
          int n = static_cast<int>(floor(x));
          n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
          Real w = x - n;
          w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
          return (1.0 - w)*a[n] + w*a[n+1];
        };
        std::vector<Real> cf(n1 + 1), cd(n1), ctt(n1), ckk, cdx(n1), ph0(n1, 1.0);
        for (int i=0; i<=n1; ++i) cf[i] = lin(hF, ax1f(0,i));
        for (int i=0; i<n1; ++i) cdx[i] = ax1f(0,i+1) - ax1f(0,i);
        hatm_c.assign(n1, 1.0);
        hatm_f.assign(nf, 1.0);
        int ait = 0, natm = 0;
        Real dph = 0.0, rmxa = 0.0, gpre = 0.0, rgpre = 0.0, gpost = 0.0;
        for (int it=0; it<amaxit; ++it) {
          ++ait;
          for (int i=0; i<n1; ++i) {
            cd[i] = flog(hd, ax1v(0,i));
            ctt[i] = flog(tcur, ax1v(0,i));
          }
          StarKtEval(otr, cd, ctt, ckk);
          dph = 0.0;
          natm = 0;
          for (int i=ia; i<=aie; ++i) {
            Real ph = 1.0;
            if (cd[i] > dfl*(1.0 + 1.0e-6)) {
              const Real rkc = cd[i]*ckk[i];
              const Real rkl = (cdx[i]*cd[i-1]*ckk[i-1] + cdx[i-1]*rkc)
                               /(cdx[i-1] + cdx[i]);
              const Real rkr = (i == aie) ? rkc
                  : (cdx[i+1]*rkc + cdx[i]*cd[i+1]*ckk[i+1])/(cdx[i] + cdx[i+1]);
              ph = (rkl*cf[i] + rkr*cf[i+1])/(rkc*(cf[i] + cf[i+1]));
              ++natm;
            }
            if (it > 0 && arel < 1.0) {
              ph = exp((1.0 - arel)*log(hatm_c[i]) + arel*log(ph));
            }
            dph = std::max(dph, std::fabs(ph/hatm_c[i] - 1.0));
            hatm_c[i] = ph;
          }
          if (it == 0) {   // the old IC's uncompensated acceleration (phi - 1) Gamma
            ph0 = hatm_c;
            for (int i=ia; i<=aie; ++i) {
              const Real gg = (hatm_c[i] - 1.0)*ckk[i]*0.5*(cf[i] + cf[i+1])/cl
                              /(star_gm_/SQR(ax1v(0,i)));
              if (gg > gpre) {
                gpre = gg;
                rgpre = ax1v(0,i);
              }
            }
          }
          if (it > 0 && dph < atol) break;
          // phi on the fine grid above the face x1f(ia): constant up to x1v(ia), then
          // linear between the cell centres (1 above the last active cell)
          for (int n=n0+1; n<nf; ++n) {
            const Real r = hr[n];
            Real p = 1.0;
            if (r <= ax1v(0,ia)) {
              p = hatm_c[ia];
            } else if (r < ax1v(0,aie)) {
              int i = ia;
              while (i < aie && ax1v(0,i+1) <= r) ++i;
              const Real w = (r - ax1v(0,i))/(ax1v(0,i+1) - ax1v(0,i));
              p = (1.0 - w)*hatm_c[i] + w*hatm_c[i+1];
            } else {
              p = hatm_c[aie];
            }
            hatm_f[n] = p;
          }
          // the outward march at fixed T from fine node n0 (as the HSE loop's, the
          // radiative acceleration times phi; the same Gamma clamp)
          StarKtEval(otr, hd, tcur, kt);
          bool flo = false;
          rmxa = 0.0;
          for (int n=n0; n<nf-1; ++n) {
            if (flo) {
              hd[n+1] = dfl;
              continue;
            }
            const Real gg0 = star_gm_/SQR(hr[n]), gg1 = star_gm_/SQR(hr[n+1]);
            const Real ar0 = hatm_f[n]*kt[n]*hF[n]/cl;
            const Real ar1 = hatm_f[n+1]*kt[n+1]*hF[n+1]/cl;
            const Real g0 = (ar0 > hgmax*gg0) ? (1.0 - hgmax)*gg0 : (gg0 - ar0);
            const Real g1 = (ar1 > hgmax*gg1) ? (1.0 - hgmax)*gg1 : (gg1 - ar1);
            Real ee, p0, cr, ct, cv;
            eos.ThermoAt(hd[n], tcur[n], ee, p0, cr, ct, cv);
            const Real rhs = p0 - 0.5*star_dr_*hd[n]*g0;
            auto res = [&](const Real lx) {
              Real e1, p1, c1, c2, c3;
              eos.ThermoAt(exp(lx), tcur[n+1], e1, p1, c1, c2, c3);
              return (p1 + 0.5*star_dr_*exp(lx)*g1)/rhs - 1.0;
            };
            Real x0 = log(hd[n]);
            if (arel == 1.0) x0 += log(std::max(hd[n+1], dfl)/hd[n]);
            Real x1 = x0 + 1.0e-6;
            Real f0 = res(x0), f1 = res(x1);
            for (int k=0; k<60 && fabs(f1) > 1.0e-14 && f1 != f0; ++k) {
              const Real x2 = x1 - f1*(x1 - x0)/(f1 - f0);
              x0 = x1; f0 = f1; x1 = x2;
              f1 = res(x1);
            }
            rmxa = std::max(rmxa, fabs(f1));
            if (!(rhs > 0.0) || !std::isfinite(x1) || !(exp(x1) > dfl)) {
              flo = true;
              hd[n+1] = dfl;
            } else {
              hd[n+1] = exp(x1);
            }
          }
          // T as re_atm: the RE T (the blend weight is 1 this far out) above the floor,
          // the column's in it
          for (int n=n0+1; n<nf; ++n) {
            const Real x = (taub - tauf[n])/(taub - taua);
            const Real w = (x >= 1.0) ? 1.0 : ((x <= 0.0) ? 0.0 : x*x*(3.0 - 2.0*x));
            tcur[n] = (hd[n] <= dfl*(1.0 + 1.0e-9)) ? t0[n]
                      : (w*pow(ef[n]/ar, 0.25) + (1.0 - w)*t0[n]);
          }
        }
        // the converged cells: the remaining uncompensated acceleration (0 by
        // construction up to the fixed-point tolerance; printed)
        for (int i=ia; i<=aie; ++i) {
          const Real gg = (hatm_c[i] - 1.0)*ckk[i]*0.5*(cf[i] + cf[i+1])/cl
                          /(star_gm_/SQR(ax1v(0,i)));
          gpost = std::max(gpost, gg);
        }
        ftau();
        hatm_on = true;
        star_atp_ia_ = ia;
        star_atp_dfl_ = dfl;
        if (global_variable::my_rank == 0) {
          Real pmx = 1.0, pmx0 = 1.0;
          for (int i=ia; i<=aie; ++i) {
            pmx = std::max(pmx, hatm_c[i]);
            pmx0 = std::max(pmx0, ph0[i]);
          }
          std::cout << "star_vet: star_ic_rad_atm_force = module: atmosphere from cell "
                    << ia << " (lower face " << ax1f(0,ia) << ", R_ph " << rpa << " -> "
                    << frph() << "), " << natm << " cells above the floor, " << ait
                    << " passes (last max |dphi| " << dph << ", march residual " << rmxa
                    << "); phi max " << pmx0 << " (old IC) -> " << pmx
                    << "; the old IC's (phi - 1) Gamma max " << gpre << " at r/R_ph "
                    << rgpre/rpa << "; floor from r = ";
          Real rfl = -1.0;
          for (int n=n0; n<nf; ++n) {
            if (hd[n] <= dfl*(1.0 + 1.0e-9)) {
              rfl = hr[n];
              break;
            }
          }
          std::cout << rfl << " (" << rfl/rpa << " R_ph)" << std::endl;
          if (!(dph < atol) && arel == 1.0) {
            std::cout << "star_vet: star_ic_rad_atm_force: WARNING not converged in "
                      << amaxit << " iterations" << std::endl;
          }
          for (int i=ia; i<=aie; ++i) {
            if (!(cd[i] > dfl*(1.0 + 1.0e-6))) break;
            std::cout << "star_vet:   atm cell " << i << " r/R_ph " << ax1v(0,i)/rpa
                      << " rho " << cd[i] << " T [K] " << ctt[i]*star_tk << " phi "
                      << ph0[i] << " -> " << hatm_c[i] << ", Gamma "
                      << ckk[i]*0.5*(cf[i] + cf[i+1])/cl/(star_gm_/SQR(ax1v(0,i)))
                      << std::endl;
          }
        }
        (void) gpost;
        if (arel < 1.0 && !(dph < atol && rmxa < 1.0e-10)) {
          StarFatal("problem/star_ic_rad_atm_force = module: the atmosphere march did "
                    "not converge (star_ic_rad_atm_maxit, star_ic_rad_atm_relax)",
                    __LINE__);
        }
      }
    }
    // write back: E, and with gas = re the T (file column) and eint of the code's EOS
    for (int n=nb; n<nf; ++n) {
      hE[n] = ef[n];
      if (gre) {
        hTcol[n] = tcur[n]*star_tk;
        Real ee, pp, cr, ct, cv;
        eos.ThermoAt(hd[n], tcur[n], ee, pp, cr, ct, cv);
        he[n] = ee;
      }
    }
    const Real rph1 = frph();
    const double wsec = tw0.seconds();
    if (global_variable::my_rank == 0) {
      const Real x1a = pmy_mesh_->mesh_size.x1min;
      Real emx = 0.0, tmx = 0.0, dmx = 0.0, rem = 0.0, rtm = 0.0, rdm = 0.0;
      for (int n=nb; n<nf; ++n) {
        if (hr[n] < x1a || hr[n] > rtop) continue;
        if (fabs(ef[n]/e0[n] - 1.0) > emx) {
          emx = fabs(ef[n]/e0[n] - 1.0);
          rem = hr[n];
        }
        if (fabs(tcur[n]/t0[n] - 1.0) > tmx) {
          tmx = fabs(tcur[n]/t0[n] - 1.0);
          rtm = hr[n];
        }
        if (fabs(hd[n]/d0[n] - 1.0) > dmx) {
          dmx = fabs(hd[n]/d0[n] - 1.0);
          rdm = hr[n];
        }
      }
      std::cout << "star_vet: star_ic_rad = sph1d (gas " << gasm << ", hse "
                << (hse ? "on" : "off") << "): r_b = " << hr[nb] << " (tau "
                << tauf[nb] << "), blend tau " << taua << " .. " << taub << ", "
                << rg.size() << " solve points + " << ncore << " core rays, " << nit
                << " VET passes, " << nout << " outer, g_top = H/J " << gtop
                << ", f_top = K/J " << ftop << ", " << wsec << " s" << std::endl;
      std::cout << "star_vet: star_ic_rad = sph1d: R_ph (tau_R 2/3, run kappa) "
                << rph0 << " -> " << rph1 << "; max |E/E_col - 1| " << emx << " at r = "
                << rem << ", |T/T_col - 1| " << tmx << " at " << rtm << ", |rho/rho_col"
                << " - 1| " << dmx << " at " << rdm;
      if (hse) {
        std::cout << "; HSE last max |dln rho| " << dlr << ", residual " << rmx
                  << ", floor " << dfl << " from r = "
                  << ((nflo >= 0) ? hr[nflo] : -1.0) << ", Gamma > " << hgmax
                  << " (clamped) in " << ncl << " fine cells"
                  << ((ncl > 0) ? ", outermost r = " + std::to_string(hr[ngcl]) : "")
                  << ", damped (" << hsdamp << ") from outer "
                  << ((hsodmp >= 0) ? std::to_string(hsodmp) : std::string("none"));
        if (!(dlr < htol)) {
          std::cout << std::endl << "star_vet: star_ic_rad = sph1d: WARNING the HSE "
                    << "outer loop did NOT converge in star_ic_rad_hse_maxit = " << hsmaxo
                    << " (max |dln rho| " << dlr << " > star_ic_rad_hse_tol " << htol
                    << ")";
        }
      }
      std::cout << std::endl;
      const Real xs[] = {0.9, 0.95, 0.98, 0.99, 0.995, 1.0, 1.003, 1.006, 1.01, 1.02,
                         1.05, 1.1, 1.2, 1.5, 2.0, 2.5, 2.9};
      for (const Real xx : xs) {
        const Real r = xx*rph0;
        if (r < hr[nb] || r > rtop) continue;
        std::cout << "star_vet:   r/R_ph " << xx << ": E/E_col "
                  << flog(ef, r)/flog(e0, r)
                  << ", T [K] " << flog(tcur, r)*star_tk << " (/T_col " << flog(tcur, r)/
                  flog(t0, r) << "), rho/rho_col " << flog(hd, r)/flog(d0, r)
                  << ", F/(cE) " << flog(hF, r)/(cl*flog(ef, r)) << std::endl;
      }
    }
    // star_ic_rad_dump: the new column on the union of the file radii and the solve grid
    if (pin->DoesParameterExist("problem","star_ic_rad_dump") &&
        pin->GetString("problem","star_ic_rad_dump") != "none" &&
        global_variable::my_rank == 0) {
      const std::string dfn = pin->GetString("problem","star_ic_rad_dump");
      // (the fine grid's span: this mesh plus its ghost margin, ends included)
      std::vector<Real> ru = {hr[0], hr[nf-1]};
      for (const Real r : hfr_keep) if (r > hr[0] && r < hr[nf-1]) ru.push_back(r);
      for (const Real r : rg) ru.push_back(r);
      std::sort(ru.begin(), ru.end());
      ru.erase(std::unique(ru.begin(), ru.end()), ru.end());
      std::vector<Real> tk(nf);
      for (int n=0; n<nf; ++n) tk[n] = tcur[n]*star_tk;
      FILE *fp = fopen(dfn.c_str(), "w");
      if (fp == nullptr) StarFatal("cannot write problem/star_ic_rad_dump", __LINE__);
      {
        std::ifstream f0(fn);
        std::string l0;
        while (std::getline(f0, l0)) {
          if (!l0.empty() && l0[0] == '#') fprintf(fp, "%s\n", l0.c_str());
        }
      }
      fprintf(fp, "# star_ic_rad = sph1d (star_vet): E (and with gas = re T, "
              "eint; with hse rho) above tau_R %g replaced; gas %s hse %d taua %g dtau "
              "%g dlnr %g ncore %d; R_ph %.10e -> %.10e; resampled from the fine grid\n",
              taub, gasm.c_str(), hse ? 1 : 0, taua, dtau0, dlnr, ncore, rph0, rph1);
      for (const Real r : ru) {
        Real fm = 0.0;
        if (hfr_keep.size() > 1) {
          const int nk = static_cast<int>(hfr_keep.size());
          int i = static_cast<int>(std::upper_bound(hfr_keep.begin(), hfr_keep.end(), r) -
                                   hfr_keep.begin()) - 1;
          i = (i < 0) ? 0 : ((i > nk - 2) ? (nk - 2) : i);
          Real w = (r - hfr_keep[i])/(hfr_keep[i+1] - hfr_keep[i]);
          w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
          fm = (1.0 - w)*hfM_keep[i] + w*hfM_keep[i+1];
        }
        fprintf(fp, "%.10e %.10e %.10e %.10e %.10e %.10e %.6e\n", r, flog(hd, r),
                flog(he, r), flog(hF, r), flog(ef, r), flog(tk, r), fm);
      }
      fclose(fp);
      std::cout << "star_vet: star_ic_rad = sph1d: column written to " << dfn << " ("
                << ru.size() << " lines)" << std::endl;
    }
  }

  // problem/star_ic_equil (read only when named; default false = nothing): see
  // StarIcEquil.  Its floor density is the column's minimum (the sph1d floor)
  star_eq_ = pin->DoesParameterExist("problem","star_ic_equil") &&
             pin->GetBoolean("problem","star_ic_equil");
  if (star_eq_) {
    star_eq_dfl_ = 1.0e300;
    for (int n=0; n<nf; ++n) star_eq_dfl_ = std::min(star_eq_dfl_, hd[n]);
  }
  // ---- the column on the device: T, E = a T^4 and kappa_t
  Kokkos::realloc(star_rho_, nf);
  Kokkos::realloc(star_eint_, nf);
  DvceArray1D<Real> cE("star_E", nf), cT("star_T", nf), ckt("star_kt", nf);
  DvceArray1D<Real> cF("star_F", nf);
  {
    auto h1 = Kokkos::create_mirror_view(star_rho_);
    auto h2 = Kokkos::create_mirror_view(star_eint_);
    auto h4 = Kokkos::create_mirror_view(cF);
    auto h5 = Kokkos::create_mirror_view(cE);
    for (int n=0; n<nf; ++n) {
      h1(n) = hd[n];
      h2(n) = he[n];
      h4(n) = hF[n];
      h5(n) = hE[n];
    }
    Kokkos::deep_copy(cE, h5);
    Kokkos::deep_copy(star_rho_, h1);
    Kokkos::deep_copy(star_eint_, h2);
    Kokkos::deep_copy(cF, h4);
  }
  {
    auto crho = star_rho_, ceint = star_eint_;
    radvet::VetOpacTab ot = pvet->otab;
    par_for("star_col", DevExeSpace(), 0, nf-1, KOKKOS_LAMBDA(const int n) {
      const Real d = crho(n);
      const Real t = eos.Temperature(d, ceint(n));
      Real op, oe, of, os;
      radvet::VetTableOpacities(ot, d, t, op, oe, of, os);
      cT(n) = t;
      ckt(n) = of + os;
      if (!(cE(n) > 0.0)) cE(n) = ar*t*t*t*t;
    });
  }
  auto hT = Kokkos::create_mirror_view(cT);
  auto hkt = Kokkos::create_mirror_view(ckt);
  Kokkos::deep_copy(hT, cT);
  Kokkos::deep_copy(hkt, ckt);
  // diagnostic (output only): the code's T(rho, eint) of the column against the file's T
  // (column 6), over the mesh's radial range; skipped with star_wind_ic (T replaced
  // there)
  if (tcol_read && !wic && global_variable::my_rank == 0) {
    const Real x1a = pmy_mesh_->mesh_size.x1min, x1b = pmy_mesh_->mesh_size.x1max;
    Real tdm = 0.0, rdm = 0.0;
    for (int n=0; n<nf; ++n) {
      if (hr[n] < x1a || hr[n] > x1b) continue;
      const Real dv = fabs(hT(n)*star_tk/hTcol[n] - 1.0);
      if (dv > tdm) {
        tdm = dv;
        rdm = hr[n];
      }
    }
    std::cout << "star_vet: IC column T(rho,eint)/T_col - 1: max " << tdm << " at r = "
              << rdm << " (star_ic_eint_from_t = " << (eft ? "true" : "false") << ")"
              << std::endl;
  }

  // ---- startup check (fatal): the column, ghost margin included, is inside the EOS
  // table (tabulated EOS only; an ideal gas has no range) and the opacity table grid,
  // in (rho [g/cm^3], T [K])
  {
    const Real tk = pvet->otab.tunit, dk = pvet->otab.dunit;
    const bool etab = eos.tbl.active;
    const Real elo_d = etab ? pin->GetReal("hydro","eos_logd_min") : -1.0e300;
    const Real ehi_d = etab ? pin->GetReal("hydro","eos_logd_max") : 1.0e300;
    const Real elo_t = etab ? pin->GetReal("hydro","eos_logt_min") : -1.0e300;
    const Real ehi_t = etab ? pin->GetReal("hydro","eos_logt_max") : 1.0e300;
    Real ldmin = 1.0e300, ldmax = -1.0e300, ltmin = 1.0e300, ltmax = -1.0e300;
    for (int n=0; n<nf; ++n) {
      const Real ld = log10(hd[n]*dk), lt = log10(hT(n)*tk);
      ldmin = std::min(ldmin, ld); ldmax = std::max(ldmax, ld);
      ltmin = std::min(ltmin, lt); ltmax = std::max(ltmax, lt);
    }
    std::ostringstream os;
    bool bad = false;
    if (ldmin < elo_d || ldmax > ehi_d || ltmin < elo_t || ltmax > ehi_t) {
      os << "the initial column leaves the EOS table: log10 rho " << ldmin << " .. "
         << ldmax << " (table " << elo_d << " .. " << ehi_d << "), log10 T " << ltmin
         << " .. " << ltmax << " (table " << elo_t << " .. " << ehi_t << ")  ";
      bad = true;
    }
    if (ldmin < tab_ld_lo || ldmax > tab_ld_hi || ltmin < tab_lt_lo ||
        ltmax > tab_lt_hi) {
      os << "the initial column leaves the Rosseland/Planck table grid: log10 rho "
         << ldmin << " .. " << ldmax << " (grid " << tab_ld_lo << " .. " << tab_ld_hi
         << "), log10 T " << ltmin << " .. " << ltmax << " (grid " << tab_lt_lo
         << " .. " << tab_lt_hi << ")  ";
      bad = true;
    }
    if (bad) StarFatal(os.str(), __LINE__);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: column inside the tables: log10 rho " << ldmin << " .. "
                << ldmax << ", log10 T " << ltmin << " .. " << ltmax << std::endl;
    }
  }
  // the imposed bottom flux must be the column's own F_r at r_in
  {
    const Real x = (star_rin_ - star_rlo_)/star_dr_;
    const int i = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w = x - i;
    const Real fcol = (1.0 - w)*hF[i] + w*hF[i+1];
    if (fabs(star_fin_/fcol - 1.0) > 1.0e-3) {
      std::ostringstream os;
      os << "<rad_vet>/implicit_flux_x1min = " << star_fin_ << " but the column has F_r("
         << star_rin_ << ") = " << fcol;
      StarFatal(os.str(), __LINE__);
    }
  }
  // tau from the top face, r_int (tau = star_tau_int), Phi_eff = Phi - int a_ref dr
  std::vector<Real> tau(nf, 0.0), aref(nf), cum(nf, 0.0);
  for (int n=0; n<nf; ++n) aref[n] = hkt(n)*hF[n]/cl;
  // the module form in the thin atmosphere (hatm_f = 1 elsewhere: bitwise)
  if (hatm_on) {
    for (int n=0; n<nf; ++n) aref[n] *= hatm_f[n];
  }
  for (int n=nf-2; n>=0; --n) {
    const Real a0 = (hr[n] < rtop) ? hd[n]*hkt(n) : 0.0;
    const Real a1 = (hr[n+1] < rtop) ? hd[n+1]*hkt(n+1) : 0.0;
    tau[n] = tau[n+1] + 0.5*star_dr_*(a0 + a1);
  }
  star_rint_ = star_rin_;
  for (int n=nf-1; n>=0; --n) {
    if (tau[n] >= tau_int) {
      star_rint_ = std::max(hr[n], star_rin_);
      break;
    }
  }
  // F7 problem/star_ic_rad_thin (default false = bitwise): the initial
  // radiation above the photosphere (tau_R < 2/3 from the mesh top, this column) is the
  // uniformly bright sphere's, E = F_r/(c f), f = (1 + mu_c)/2, mu_c = sqrt(1 -
  // (R_ph/r)^2), instead of the column's a T^4; F_r is the column's own (= L/(4 pi r^2)
  // there; printed).  Below R_ph E is blended into the column's: E = w 2F/c + (1 - w)
  // E_col with w a smoothstep in tau from 1 at 2/3 to 0 at star_ic_rad_thin_taub (default
  // 2), so E is continuous at R_ph (f = 1/2 there).  Gas T, rho, F, a_ref unchanged
  // (the thin gas is decoupled; changing T would break the column's balance).
  if (pin->GetOrAddBoolean("problem","star_ic_rad_thin",false)) {
    const Real taub = pin->GetOrAddReal("problem","star_ic_rad_thin_taub",2.0);
    const Real t23 = 2.0/3.0;
    if (!(taub > t23)) StarFatal("problem/star_ic_rad_thin_taub must exceed 2/3",
                                 __LINE__);
    Real rph = -1.0;
    for (int n=nf-2; n>=0; --n) {
      if (tau[n] >= t23 && tau[n+1] < t23) {
        rph = hr[n] + star_dr_*(tau[n] - t23)/(tau[n] - tau[n+1]);
        break;
      }
    }
    if (!(rph > 0.0)) StarFatal("problem/star_ic_rad_thin: tau_R = 2/3 not on the column",
                                __LINE__);
    auto hEc = Kokkos::create_mirror_view(cE);
    Kokkos::deep_copy(hEc, cE);
    const Real lph = 4.0*M_PI*rph*rph;
    Real fdev = 0.0, lum = 0.0;
    {
      const Real x = (rph - star_rlo_)/star_dr_;
      const int i = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
      const Real w = x - i;
      lum = lph*((1.0 - w)*hF[i] + w*hF[i+1]);
    }
    int nmod = 0;
    for (int n=0; n<nf; ++n) {
      if (tau[n] >= taub) continue;
      const Real r = hr[n], f = hF[n];
      Real e;
      if (r >= rph) {
        const Real mu = sqrt(fmax(1.0 - SQR(rph/r), 0.0));
        e = f/(cl*0.5*(1.0 + mu));
        if (r <= rtop) fdev = std::max(fdev, std::fabs(4.0*M_PI*r*r*f/lum - 1.0));
      } else {
        const Real x = (taub - tau[n])/(taub - t23);   // 0 at taub, 1 at 2/3
        const Real w = (x >= 1.0) ? 1.0 : x*x*(3.0 - 2.0*x);
        e = w*2.0*f/cl + (1.0 - w)*hEc(n);
      }
      hEc(n) = e;
      ++nmod;
    }
    Kokkos::deep_copy(cE, hEc);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: star_ic_rad_thin: R_ph (tau_R = 2/3) = " << rph
                << ", L = " << lum << ", blend to tau " << taub << ", column points "
                << nmod << ", max |4 pi r^2 F/L - 1| above R_ph (mesh) = " << fdev
                << std::endl;
    }
  }
  for (int n=1; n<nf; ++n) cum[n] = cum[n-1] + 0.5*star_dr_*(aref[n-1] + aref[n]);
  Real cin = 0.0;
  {
    const Real x = (star_rin_ - star_rlo_)/star_dr_;
    int i = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w = x - i;
    cin = (1.0 - w)*cum[i] + w*cum[i+1];
  }
  // star_gm_column: G m(r) and Phi(r) of the column on the fine grid (trapezoid sums,
  // both zero-referenced at r_in like the point-mass Phi = GM (1/r_in - 1/r))
  std::vector<Real> hphig;
  if (star_gmc_) {
    if (pmbp->punit == nullptr) StarFatal("problem/star_gm_column needs <units>",
                                          __LINE__);
    const Real lu = pmbp->punit->length_cgs(), tu = pmbp->punit->time_cgs();
    const Real gcode = units::Units::grav_constant_cgs*pmbp->punit->mass_cgs()*tu*tu/
                       (lu*lu*lu);
    std::vector<Real> cm(nf, 0.0), hgm(nf), cp(nf, 0.0);
    for (int n=1; n<nf; ++n) {
      cm[n] = cm[n-1] + 0.5*star_dr_*4.0*M_PI*(hr[n-1]*hr[n-1]*hd[n-1] +
                                               hr[n]*hr[n]*hd[n]);
    }
    const Real x = (star_rin_ - star_rlo_)/star_dr_;
    const int i0 = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w0 = x - i0;
    const Real cm0 = (1.0 - w0)*cm[i0] + w0*cm[i0+1];
    for (int n=0; n<nf; ++n) hgm[n] = star_gm_ + gcode*(cm[n] - cm0);
    for (int n=1; n<nf; ++n) {
      cp[n] = cp[n-1] + 0.5*star_dr_*(hgm[n-1]/(hr[n-1]*hr[n-1]) + hgm[n]/(hr[n]*hr[n]));
    }
    const Real cp0 = (1.0 - w0)*cp[i0] + w0*cp[i0+1];
    hphig.resize(nf);
    for (int n=0; n<nf; ++n) hphig[n] = cp[n] - cp0;
    Kokkos::realloc(star_phig_, nf);
    Kokkos::realloc(star_gmr_, nf);
    auto h1 = Kokkos::create_mirror_view(star_phig_);
    auto h2 = Kokkos::create_mirror_view(star_gmr_);
    for (int n=0; n<nf; ++n) {
      h1(n) = hphig[n];
      h2(n) = hgm[n];
    }
    Kokkos::deep_copy(star_phig_, h1);
    Kokkos::deep_copy(star_gmr_, h2);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: star_gm_column: G m(r_in) = " << star_gm_
                << ", G m(r_top) = " << hgm[nf-1] << " (G = " << gcode
                << " code), static monopole of the initial column" << std::endl;
    }
  }
  DvceArray1D<Real> cphi("star_phieff", nf);
  {
    auto hp = Kokkos::create_mirror_view(cphi);
    for (int n=0; n<nf; ++n) {
      if (star_gmc_) {
        hp(n) = hphig[n] - (cum[n] - cin);
      } else {
        hp(n) = star_gm_*(1.0/star_rin_ - 1.0/hr[n]) - (cum[n] - cin);
      }
    }
    Kokkos::deep_copy(cphi, hp);
  }

  // ---- potentials: the TRUE one (etotgrav, the WB default) ...
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  const Real gm = star_gm_, rin = star_rin_;
  if (star_gmc_) {
    auto phicc = ph->phicc0;
    auto ph1 = ph->phi0->x1f, ph2 = ph->phi0->x2f, ph3 = ph->phi0->x3f;
    auto pg = star_phig_;
    const Real rlo0 = star_rlo_, dr0 = star_dr_;
    par_for("star_phi_col", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real pc = StarLinInterp(pg, rlo0, dr0, nf, x1v(m,i));
      phicc(m,k,j,i) = pc;
      ph1(m,k,j,i) = StarLinInterp(pg, rlo0, dr0, nf, x1f(m,i));
      if (i == n1m1) ph1(m,k,j,i+1) = StarLinInterp(pg, rlo0, dr0, nf, x1f(m,i+1));
      ph2(m,k,j,i) = pc;
      ph3(m,k,j,i) = pc;
      if (j == n2m1) ph2(m,k,j+1,i) = pc;
      if (k == n3m1) ph3(m,k+1,j,i) = pc;
    });
  } else {
    auto phicc = ph->phicc0;
    auto ph1 = ph->phi0->x1f, ph2 = ph->phi0->x2f, ph3 = ph->phi0->x3f;
    par_for("star_phi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real pc = gm*(1.0/rin - 1.0/x1v(m,i));
      phicc(m,k,j,i) = pc;
      ph1(m,k,j,i) = gm*(1.0/rin - 1.0/x1f(m,i));
      if (i == n1m1) ph1(m,k,j,i+1) = gm*(1.0/rin - 1.0/x1f(m,i+1));
      ph2(m,k,j,i) = pc;
      ph3(m,k,j,i) = pc;
      if (j == n2m1) ph2(m,k,j+1,i) = pc;
      if (k == n3m1) ph3(m,k+1,j,i) = pc;
    });
  }
  // ... and the EFFECTIVE one of the x1 well-balanced pair
  const Real rlo = star_rlo_, dr = star_dr_;
  // (star_gravity = pointmass: no WB pair, no reference force; phicc_wb aliases phicc0)
  if (!star_pm_) {
    ph->EnableWBEffectivePotential();
    auto pwc = ph->phicc_wb, pwf = ph->phi_wb_x1f;
    par_for("star_phieff", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pwc(m,k,j,i) = StarLinInterp(cphi, rlo, dr, nf, x1v(m,i));
      pwf(m,k,j,i) = StarLinInterp(cphi, rlo, dr, nf, x1f(m,i));
      if (i == n1m1) pwf(m,k,j,i+1) = StarLinInterp(cphi, rlo, dr, nf, x1f(m,i+1));
    });
  }
  // the reference acceleration, face form: kt(cell) (F(r_l) + F(r_r))/(2c), from the
  // INITIAL column (identical on a restart)
  if (!star_pm_) {
    DvceArray4D<Real> aref_d("star_aref", nmb1+1, n3m1+1, n2m1+1, n1m1+1);
    par_for("star_aref", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real kt = StarLinInterp(ckt, rlo, dr, nf, x1v(m,i));
      const Real rl = x1f(m,i), rr = x1f(m,i+1);
      aref_d(m,k,j,i) = 0.5*kt*(StarLinInterp(cF, rlo, dr, nf, rl)
                                + StarLinInterp(cF, rlo, dr, nf, rr))/cl;
    });
    // (star_ic_rad_atm_force = module): the thin-atmosphere cells take the
    // module's face-mean form, phi times the continuum value (phi = 1 elsewhere)
    if (hatm_on) {
      DvceArray1D<Real> aphi("star_aphi", n1m1+1);
      auto hp = Kokkos::create_mirror_view(aphi);
      for (int i=0; i<=n1m1; ++i) hp(i) = hatm_c[i];
      Kokkos::deep_copy(aphi, hp);
      par_for("star_aref_atm", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        aref_d(m,k,j,i) *= aphi(i);
      });
    }
    // StarGravity gives the reference work (split); set before SetForceReference,
    // which resolves force_reference_work = auto
    pvet->fref_wsplit_ok = true;
    pvet->SetForceReference(aref_d);
  }

  // ---- DISCRETE hydrostatic balance of the initial cells (problem/star_ic_balance,
  // default false).  The x1 well-balanced scheme is exactly static when, at every
  // interior face, the pressure of the local background of the cell on the left
  // (walked from its centre to the face along Phi_eff, WBBackgroundStencil) equals that
  // of the cell on the right: the Riemann states are then the backgrounds themselves
  // (the deviations vanish) and the WB source cancels the flux divergence.  For
  // wb_option = polytropic the walk of cell i depends only on (d_i, T_i) and on the
  // stencil's dT/dPhi, i.e. on T_{i-1}, T_{i+1}: with T FIXED (the column's, set by the
  // radiation) each face condition PR_i(d_i) = PL_{i+1}(d_{i+1}) is one equation for d_i.
  // It is solved from the first outer ghost (anchor: the column) inward to the first
  // active cell, with the code's own EOS and walk (host), so it is exact to round-off.
  // The inner ghosts are the column scaled to the edge cell, as StarBC makes them;
  // their T enters the stencil of the edge cell, which is iterated.  The inner wall face
  // (wall_closed_ix1, mirror state) is then balanced by construction.  E and F_r keep the
  // column's values.  Rebuilt on a restart (the BCs read the balanced ghost profiles).
  // With star_ic_profile the column is the profile file's (rho, T, F_r incl. the inner
  // ghosts); its seed then multiplies the balanced rho (isobaric, in its block below).
  star_bal_ = pin->GetOrAddBoolean("problem","star_ic_balance",false);
  // problem/star_ic_profile (read only when named; default false = this
  // IC bitwise, no key added to the input; deprecated name star_ic_linhao) replaces the
  // initial state by the construction of the BSG reference Athena++ setup bsg_star.cpp
  // (the block at the end of this function).  Its profile is read here:
  // problem/star_ic_balance re-marches it (also on a restart).
  const bool hlin = pin->DoesParameterExist("problem",StarIcKey(pin,"star_ic_profile")) &&
                    pin->GetBoolean("problem",StarIcKey(pin,"star_ic_profile"));
  std::string lfn;
  std::vector<Real> pr, pT, pd, pM, pF;
  bool gxi = false;
  Real lrmax = 0.0, ltmax = 0.0, ldmax = 0.0, lgmax = 0.0;
  if (hlin && (star_bal_ || !restart)) {
    lfn = pin->GetString("problem",StarIcKey(pin,"star_ic_profile_file"));
    // Athena++ Constants::grav_const_cgs (bsg_star.cpp's grav_rmax: G, the file's M)
    const Real lG = pin->GetOrAddReal("problem",StarIcKey(pin,"star_ic_profile_gconst"),
                                      6.67408e-8);
    // star_ic_profile_gext = r2 (default, bsg_star.cpp:429 g = g(rmax)
    // (r/rmax)^+2, as coded) | invsq (the run's point mass GM/r^2, hydrostatic; see
    // below)
    const std::string gxs = pin->GetOrAddString("problem",
                                              StarIcKey(pin,"star_ic_profile_gext"),"r2");
    if (gxs != "r2" && gxs != "invsq") {
      StarFatal("problem/star_ic_profile_gext must be r2 or invsq", __LINE__);
    }
    gxi = (gxs == "invsq");
    if (gxi && star_gmc_) {
      StarFatal("problem/star_ic_profile_gext = invsq needs the point-mass gravity",
                __LINE__);
    }
    {
      std::ifstream f(lfn);
      if (!f.good()) StarFatal("cannot open problem/star_ic_profile_file '" + lfn + "'",
                               __LINE__);
      std::string line;
      while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        Real a, b, c, d, e, g;
        if (!(ss >> a >> b >> c >> d >> e >> g)) continue;
        pr.push_back(a); pT.push_back(b); pd.push_back(c); pM.push_back(d);
        pF.push_back(g);
      }
    }
    if (pr.size() < 2) StarFatal("problem/star_ic_profile_file: fewer than 2 rows",
                                 __LINE__);
    lrmax = pr.back(); ltmax = pT.back(); ldmax = pd.back();
    lgmax = lG*pM.back()/(lrmax*lrmax);
  }
  // the profile at x1: r <= rmax rho, T [K], F_r LINEAR in r between the bracketing
  // rows (bsg_star.cpp ProblemGenerator's lleft/lright search, verbatim); r > rmax
  // the atmosphere (see the star_ic_profile block), F_r = 0
  auto lprof = [&](const Real x1, Real &rho, Real &tem, Real &flx) {
    const int np = static_cast<int>(pr.size());
    flx = 0.0;
    if (x1 > lrmax) {
      tem = ltmax;
      const Real tc = tem/star_tk;   // P/rho of this EOS at T(rmax)
      if (gxi) {
        // star_ic_profile_gext = invsq: the isothermal hydrostatic atmosphere of the
        // run's own point mass, g = GM/r^2 (= g(rmax) (rmax/r)^2 with g(rmax) =
        // GM/rmax^2): rho = rho(rmax) exp(-(GM/(P/rho)) (1/rmax - 1/r))
        rho = ldmax*std::exp(-(star_gm_/tc)*(1.0/lrmax - 1.0/x1));
      } else {
        const Real gl = lgmax*std::pow(x1/lrmax, 2.0);
        rho = ldmax*std::exp(-gl*(x1 - lrmax)/tc);
      }
      rho = std::max(rho, eos.dfloor);
      return;
    }
    int lleft = 0, lright = 1;
    while ((x1 > pr[lright]) && (lright < np - 1)) lright = lright + 1;
    if (lright - lleft > 1) lleft = lright - 1;
    rho = pd[lleft]
          + (x1 - pr[lleft])*(pd[lright] - pd[lleft])/(pr[lright] - pr[lleft]);
    tem = pT[lleft]
          + (x1 - pr[lleft])*(pT[lright] - pT[lleft])/(pr[lright] - pr[lleft]);
    flx = pF[lleft]
          + (x1 - pr[lleft])*(pF[lright] - pF[lleft])/(pr[lright] - pr[lleft]);
  };
  if (star_bal_ && star_pm_) {
    StarFatal("problem/star_ic_balance = true needs the well-balanced pair (star_gravity "
              "= wb)", __LINE__);
  }
  if (star_bal_) {
    if (ph->wb_option != WBOption::polytropic) {
      StarFatal("problem/star_ic_balance = true needs <hydro>/wb_option = polytropic",
                __LINE__);
    }
    const int n1 = n1m1 + 1;
    const int is = indcs.is, ie = indcs.ie, k0 = indcs.ks, j0 = indcs.js;
    auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1v);
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ph->phicc_wb);
    auto hpf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ph->phi_wb_x1f);
    std::vector<Real> dc(n1), ec(n1), tc(n1), bd(n1), be(n1);
    auto logint = [&](const std::vector<Real> &a, const Real r) {
      Real x = (r - star_rlo_)/star_dr_;
      int i = static_cast<int>(floor(x));
      i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
      Real w = x - i;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return exp((1.0 - w)*log(a[i]) + w*log(a[i+1]));
    };
    for (int i=0; i<n1; ++i) {
      dc[i] = logint(hd, hx1v(0,i));
      ec[i] = logint(he, hx1v(0,i));
      if (hlin) {   // star_ic_profile: the profile (rho, T) is the column
        Real tl, fl, p, cr, ct, cv;
        lprof(hx1v(0,i), dc[i], tl, fl);
        eos.ThermoAt(dc[i], tl/star_tk, ec[i], p, cr, ct, cv);
      }
      tc[i] = eos.Temperature(dc[i], ec[i]);
      bd[i] = dc[i];
      be[i] = ec[i];
    }
    const WBOption wbo = ph->wb_option;
    auto walk = [&](const int i, const Real di, Real &pimh, Real &piph) {
      Real ei, p, cr, ct, cv;
      eos.ThermoAt(di, tc[i], ei, p, cr, ct, cv);
      WBState s0, s1, s2, s3, s4;
      WBBackgroundStencil(eos, wbo, bd[i-1], di, bd[i+1], be[i-1], ei, be[i+1],
                          hpc(0,k0,j0,i-1), hpf(0,k0,j0,i), hpc(0,k0,j0,i),
                          hpf(0,k0,j0,i+1), hpc(0,k0,j0,i+1), s0, s1, s2, s3, s4,
                          tc[i-1], tc[i], tc[i+1]);
      pimh = s1.p;
      piph = s3.p;
    };
    Real rmax = 0.0;
    auto solve = [&](const int i, const Real target) {  // PR_i(d_i) = target
      Real pl, pr;
      Real x0 = log(bd[i]), x1 = x0 + 1.0e-6;
      walk(i, exp(x0), pl, pr);
      Real f0 = pr/target - 1.0;
      walk(i, exp(x1), pl, pr);
      Real f1 = pr/target - 1.0;
      for (int it=0; it<60 && fabs(f1) > 1.0e-15 && f1 != f0; ++it) {
        const Real x2 = x1 - f1*(x1 - x0)/(f1 - f0);
        x0 = x1; f0 = f1; x1 = x2;
        walk(i, exp(x1), pl, pr);
        f1 = pr/target - 1.0;
      }
      rmax = std::max(rmax, fabs(f1));
      bd[i] = exp(x1);
      Real p, cr, ct, cv;
      eos.ThermoAt(bd[i], tc[i], be[i], p, cr, ct, cv);
    };
    // THE WALL CELL IN THE MODULE'S OWN DISCRETE FORM.  In the IC state <rad_vet> gives
    // cell i the radiative force [(rho k)_{i-1/2} F_{i-1/2} + (rho k)_{i+1/2} F_{i+1/2}]
    // /(2c), (rho k)_face = the mean of the two cells but, at a PHYSICAL x1 face, the
    // cell's own value (implicit_bmom_half).  Inside, this equals the continuum
    // reference kt F/c of the column to O(dx^2) (~1e-4 rho g, measured); in the two
    // edge cells it differs by -(rho k)' dx/4: -6.6e-3 rho g in the wall cell of the
    // mlt column (v2/runs/mE03, predicted -6.59e-3).  So in the WALL cells a_ref is
    // set to the module's value and Phi_eff is rebuilt across them with that
    // acceleration (the ghosts' Phi_eff shifted with the wall face), so that the WB
    // pair and the module's force agree there too.  NOT done: (a) EVERY cell (in the
    // Gamma ~ 1 layers the balance then depends on a_ref to 1e-3 and the rho <-> kappa
    // fixed point diverges); (b) the TOP cell: the march is anchored at the top ghost
    // and, at fixed T, hydrostatic balance is homogeneous in rho, so a change there
    // rescales rho of the WHOLE column (-1.2 %, v2/runs/mG2) against the reference.
    auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1f);
    auto linF = [&](const Real r) {
      Real x = (r - star_rlo_)/star_dr_;
      int n = static_cast<int>(floor(x));
      n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
      Real w = x - n;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return (1.0 - w)*hF[n] + w*hF[n+1];
    };
    auto phig = [&](const Real r) {   // star_gm_column: the column's Phi, linear
      Real x = (r - star_rlo_)/star_dr_;
      int n = static_cast<int>(floor(x));
      n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
      Real w = x - n;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return (1.0 - w)*hphig[n] + w*hphig[n+1];
    };
    auto seg = [&](const Real a_, const Real ra, const Real rb) {   // Phi_eff(rb) - (ra)
      if (star_gmc_) return phig(rb) - phig(ra) - a_*(rb - ra);
      return star_gm_*(1.0/ra - 1.0/rb) - a_*(rb - ra);
    };
    // star_wind_ic: the march starts below r_c (the wind cells keep their IC)
    int ianc = ie;
    if (wic) {
      while (ianc > is && !(hx1v(0,ianc) < wrc)) --ianc;
    }
    // problem/star_ic_balance_rmax [cm] (default none = the top ghost): the march starts
    // below this radius and the cells above keep the column state.  For a column that
    // ends in a constant density floor (make_ic_mlt_star.py --atm iso): the floor is not
    // hydrostatic, and anchored there the march would rescale the whole column
    if (pin->DoesParameterExist("problem","star_ic_balance_rmax")) {
      const Real rbm = pin->GetReal("problem","star_ic_balance_rmax");
      while (ianc > is && !(hx1v(0,ianc) < rbm)) --ianc;
    }
    // the anchor cell ianc+1 takes the COLUMN's own state for the march (the wind
    // there is not hydrostatic; its density would rescale the whole column)
    if (wic && ianc < ie) {
      const int ia = ianc + 1;
      dc[ia] = logint(hdc, hx1v(0,ia));
      ec[ia] = logint(hec, hx1v(0,ia));
      tc[ia] = eos.Temperature(dc[ia], ec[ia]);
      bd[ia] = dc[ia];
      be[ia] = ec[ia];
    }
    auto march = [&]() {
      for (int i=ianc; i>=is; --i) {
        Real pl, pr;
        walk(i+1, bd[i+1], pl, pr);
        solve(i, pl);
      }
      for (int it=0; it<8; ++it) {
        for (int g=0; g<is; ++g) {
          bd[g] = dc[g]*bd[is]/dc[is];
          be[g] = ec[g]*be[is]/ec[is];
          tc[g] = eos.Temperature(bd[g], be[g]);
        }
        Real pl, pr;
        walk(is+1, bd[is+1], pl, pr);
        solve(is, pl);
      }
    };
    // the wall cell AND the next NW-1 cells take the module form: moving rho of the wall
    // cell moves (rho k) of the face above it, i.e. the force of cell is+1 (+6.3e-4
    // rho g with the wall cell alone, v2/runs/mG3); the kink dies out within a few cells.
    // Iterated with the balance (a depends on rho); the face x1f(is+NW) keeps its
    // continuum Phi_eff, the cells below are rebuilt downward from it.
    // problem/star_ic_balance_aref = file (default; read only when named) | dpr (implied
    // by star_ic_profile): a_ref = -(dPr/dr)/rho of the column (Pr = a T^4/3 at the
    // faces, the column's rho in the cell), the discrete force of the steady start at the
    // column's T, in every cell (no module form in the wall cells); Phi_eff rebuilt
    // upward from the first face with it.  file = kt F_file/c (above)
    const bool adpr = hlin ||
                      (pin->DoesParameterExist("problem","star_ic_balance_aref") &&
                       pin->GetString("problem","star_ic_balance_aref") == "dpr");
    auto tfc = [&](const Real r) {   // the column's code T at a face
      Real dl, tl, fl;
      if (hlin) {lprof(r, dl, tl, fl); return tl/star_tk;}
      return eos.Temperature(logint(hd, r), logint(he, r));
    };
    if (adpr) {
      auto har = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pvet->arad_ref);
      std::vector<Real> af(n1);
      Real prl = ar*SQR(SQR(tfc(hx1f(0,0))))/3.0;
      for (int i=0; i<n1; ++i) {
        const Real prr = ar*SQR(SQR(tfc(hx1f(0,i+1))))/3.0;
        af[i] = -(prr - prl)/(dc[i]*(hx1f(0,i+1) - hx1f(0,i)));
        prl = prr;
      }
      for (int m=0; m<=nmb1; ++m) {
        for (int k=0; k<=n3m1; ++k) {
          for (int j=0; j<=n2m1; ++j) {
            for (int i=0; i<n1; ++i) {
              hpc(m,k,j,i) = hpf(m,k,j,i) + seg(af[i], hx1f(0,i), hx1v(0,i));
              hpf(m,k,j,i+1) = hpf(m,k,j,i) + seg(af[i], hx1f(0,i), hx1f(0,i+1));
              har(m,k,j,i) = af[i];
            }
          }
        }
      }
      Kokkos::deep_copy(pvet->arad_ref, har);
    }
    const int NW = adpr ? 0 : std::min(4, ie - is);
    DvceArray1D<Real> kd_d("star_kd", NW + 1), kt_d("star_kt", NW + 1);
    DvceArray1D<Real> kk_d("star_kk", NW + 1);
    auto hkd = Kokkos::create_mirror_view(kd_d);
    auto hkt = Kokkos::create_mirror_view(kt_d);
    radvet::VetOpacTab ot1 = pvet->otab;
    std::vector<Real> aw(NW, 0.0);
    const Real pf_top = hpf(0,k0,j0,is+NW);
    const Real pf_old = hpf(0,k0,j0,is);
    std::vector<Real> gpc(is), gpf(is);
    for (int i=0; i<is; ++i) {
      gpc[i] = hpc(0,k0,j0,i) - pf_old;
      gpf[i] = hpf(0,k0,j0,i) - pf_old;
    }
    Real wchg = 0.0;
    for (int outer=0; outer<8; ++outer) {
      for (int n=0; n<=NW; ++n) {
        hkd(n) = bd[is+n];
        hkt(n) = tc[is+n];
      }
      Kokkos::deep_copy(kd_d, hkd);
      Kokkos::deep_copy(kt_d, hkt);
      par_for("star_kap1", DevExeSpace(), 0, NW, KOKKOS_LAMBDA(const int n) {
        Real op, oe, of, os;
        radvet::VetTableOpacities(ot1, kd_d(n), kt_d(n), op, oe, of, os);
        kk_d(n) = of + os;
      });
      auto hkk = Kokkos::create_mirror_view_and_copy(HostMemSpace(), kk_d);
      wchg = 0.0;
      for (int n=0; n<NW; ++n) {
        const int i = is + n;
        const Real rkc = bd[i]*hkk(n), rku = bd[i+1]*hkk(n+1);
        const Real rkl = (n == 0) ? rkc : 0.5*(rkc + bd[i-1]*hkk(n-1));
        const Real anew = 0.5*(rkl*linF(hx1f(0,i)) + 0.5*(rkc + rku)*linF(hx1f(0,i+1)))
                          /(cl*bd[i]);
        if (outer > 0) wchg = std::max(wchg, fabs(anew/aw[n] - 1.0));
        aw[n] = anew;
      }
      // Phi_eff of the wall cells, downward from the face x1f(is+NW); ghosts follow
      std::vector<Real> pf(NW + 1), pc(NW);
      pf[NW] = pf_top;
      for (int n=NW-1; n>=0; --n) {
        const int i = is + n;
        pf[n] = pf[n+1] - seg(aw[n], hx1f(0,i), hx1f(0,i+1));
        pc[n] = pf[n] + seg(aw[n], hx1f(0,i), hx1v(0,i));
      }
      for (int m=0; m<=nmb1; ++m) {
        for (int k=0; k<=n3m1; ++k) {
          for (int j=0; j<=n2m1; ++j) {
            for (int i=0; i<is; ++i) {
              hpc(m,k,j,i) = gpc[i] + pf[0];
              hpf(m,k,j,i) = gpf[i] + pf[0];
            }
            for (int n=0; n<NW; ++n) {
              hpf(m,k,j,is+n) = pf[n];
              hpc(m,k,j,is+n) = pc[n];
            }
          }
        }
      }
      rmax = 0.0;
      march();
      if (outer > 0 && wchg < 1.0e-13) break;
    }
    {
      auto har = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pvet->arad_ref);
      for (int m=0; m<=nmb1; ++m) {
        for (int k=0; k<=n3m1; ++k) {
          for (int j=0; j<=n2m1; ++j) {
            for (int n=0; n<NW; ++n) har(m,k,j,is+n) = aw[n];
          }
        }
      }
      Kokkos::deep_copy(ph->phicc_wb, hpc);
      Kokkos::deep_copy(ph->phi_wb_x1f, hpf);
      Kokkos::deep_copy(pvet->arad_ref, har);
    }
    Real dmax = 0.0;
    int imax = is;
    for (int i=is; i<=ie; ++i) {
      if (fabs(bd[i]/dc[i] - 1.0) > dmax) {
        dmax = fabs(bd[i]/dc[i] - 1.0);
        imax = i;
      }
    }
    if ((tcol_read || hlin) && !wic && global_variable::my_rank == 0) {   // diagnostic
      Real tdm = 0.0, rdm = 0.0;
      for (int i=is; i<=ie; ++i) {
        Real tr = 0.0, dl, fl;
        if (hlin) lprof(hx1v(0,i), dl, tr, fl);
        const Real dv = fabs(eos.Temperature(bd[i], be[i])*star_tk
                             /(hlin ? tr : logint(hTcol, hx1v(0,i))) - 1.0);
        if (dv > tdm) {
          tdm = dv;
          rdm = hx1v(0,i);
        }
      }
      std::cout << "star_vet: star_ic_balance cells: max |T/T_col - 1| = " << tdm
                << " at r = " << rdm << std::endl;
    }
    if (wic && ianc < ie) {   // the anchor cell is a wind cell again
      bd[ianc+1] = logint(hd, hx1v(0,ianc+1));
      be[ianc+1] = logint(he, hx1v(0,ianc+1));
    }
    Kokkos::realloc(star_bd_, n1);
    Kokkos::realloc(star_be_, n1);
    auto h1 = Kokkos::create_mirror_view(star_bd_);
    auto h2 = Kokkos::create_mirror_view(star_be_);
    for (int i=0; i<n1; ++i) {
      h1(i) = bd[i];
      h2(i) = be[i];
    }
    Kokkos::deep_copy(star_bd_, h1);
    Kokkos::deep_copy(star_be_, h2);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: star_ic_balance: max |rho/rho_col - 1| = " << dmax
                << " at r = " << hx1v(0,imax) << ", max face residual |PR/PL - 1| = "
                << rmax << ", wall cells a_ref (module form), last change " << wchg
                << "; wall ghost/edge rho " << bd[is-1]/bd[is] << " T " << tc[is-1]/tc[is]
                << std::endl;
    }
  }

  // ---- a RESTART of a run whose fresh start re-balanced the atmosphere
  // (StarIcPost) applies the a_ref it recorded in the input, and the Phi_eff rebuilt
  // from it above the face x1f(ia), exactly as the fresh start did
  if (restart && pin->DoesParameterExist("problem","star_ic_rad_atm_na")) {
    const int ia = pin->GetInteger("problem","star_ic_rad_atm_ia");
    const int na = pin->GetInteger("problem","star_ic_rad_atm_na");
    if (ia < indcs.is + 2 || na < 0 || ia + na - 1 > indcs.ie) {
      StarFatal("problem/star_ic_rad_atm_ia, _na do not match this radial grid",
                __LINE__);
    }
    std::vector<Real> pc0, pf0, p0c, p0f, x1v, x1f, ar0, pc, pf;
    StarAtmRadial(pmbp, pc0, pf0, p0c, p0f, x1v, x1f, ar0);
    std::vector<Real> a(indcs.ie - ia + 1);
    for (int i=ia; i<=indcs.ie; ++i) a[i - ia] = ar0[i];
    char key[40];
    for (int n=0; n<na; ++n) {
      snprintf(key, sizeof(key), "star_ic_rad_atm_a%03d", n);
      a[n] = std::strtod(pin->GetString("problem", key).c_str(), nullptr);
    }
    StarAtmPhiEff(ia, indcs.ie, a, pc0, pf0, p0c, p0f, x1v, x1f, pc, pf);
    StarAtmWrite(pmbp, ia, a, pc, pf);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: star_ic_rad_atm_a*: a_ref (cells " << ia << " .. "
                << ia + na - 1 << ") and Phi_eff from the fresh start's re-balance"
                << std::endl;
    }
  }

  // ---- the frozen MLT flux on the x1 faces (problem/mlt_flux_frozen)
  if (star_mlt_) {
    auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1f);
    const int nfc = n1m1 + 2;
    Kokkos::realloc(star_fm_, nfc);
    auto hfm = Kokkos::create_mirror_view(star_fm_);
    Real fmax = 0.0, rmx = 0.0;
    for (int i=0; i<nfc; ++i) {
      const Real r = hx1f(0,i);
      Real x = (r - star_rlo_)/star_dr_;
      int n = static_cast<int>(floor(x));
      n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
      Real w = x - n;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      const Real fm = (1.0 - w)*hM[n] + w*hM[n+1];
      const Real fr = exp((1.0 - w)*log(hF[n]) + w*log(hF[n+1]));
      hfm(i) = (fm > 0.0) ? fr*fm/(1.0 - fm) : 0.0;
      if (fm > fmax) {
        fmax = fm;
        rmx = r;
      }
    }
    Kokkos::deep_copy(star_fm_, hfm);
    // star_base_heat: the base-heating flux F_b on the same faces (see star_bheat_)
    if (star_bheat_ > 0.0) {
      if (pin->GetOrAddString("problem","mlt_closure","frozen") != "frozen") {
        StarFatal("problem/star_base_heat needs mlt_closure = frozen", __LINE__);
      }
      Kokkos::realloc(star_fb_, nfc);
      auto hfb = Kokkos::create_mirror_view(star_fb_);
      const Real fw = hfm(indcs.is);
      for (int i=0; i<nfc; ++i) {
        const Real r = hx1f(0,i);
        Real x = 1.0 - (r - star_rin_)/star_bheat_;
        x = (x < 0.0) ? 0.0 : ((x > 1.0) ? 1.0 : x);
        hfb(i) = fw*SQR(star_rin_/r)*x;
      }
      Kokkos::deep_copy(star_fb_, hfb);
      if (global_variable::my_rank == 0) {
        std::cout << "star_vet: star_base_heat: F_MLT(r_in) = " << fw << " kept as base "
                  << "heating over r_in .. r_in + " << star_bheat_ << " when w(t) < 1"
                  << std::endl;
      }
    }
    {
      const auto &ms = pmy_mesh_->mesh_size;
      const Real omg = hcs ? (4.0*M_PI)
                           : (cos(ms.x2min) - cos(ms.x2max))*(ms.x3max - ms.x3min);
      star_lmx_ = 0.0;
      star_r2o_.assign(nfc, 0.0);
      for (int i=0; i<nfc; ++i) {
        star_r2o_[i] = SQR(hx1f(0,i))*omg;
        star_lmx_ = std::max(star_lmx_, hfm(i)*star_r2o_[i]);
      }
    }
    // the adaptive closure's radial taper, from the IC's F_MLT/F on the faces
    star_tap_.assign(nfc, 0.0);
    {
      const Real tp = pin->GetOrAddReal("problem","mlt_closure_taper",0.01);
      for (int i=0; i<nfc; ++i) {
        const Real r = hx1f(0,i);
        Real x = (r - star_rlo_)/star_dr_;
        int n = static_cast<int>(floor(x));
        n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
        Real w = x - n;
        w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
        const Real fm = (1.0 - w)*hM[n] + w*hM[n+1];
        Real s = (tp > 0.0) ? fm/tp : ((fm > 0.0) ? 1.0 : 0.0);
        s = (s < 0.0) ? 0.0 : ((s > 1.0) ? 1.0 : s);
        star_tap_[i] = (hfm(i) > 0.0) ? s : 0.0;
      }
    }
    // the RADIATION takes the deposit (<rad_vet>/esrc, added to the implicit solve's old
    // vector with the stage weights): beta = Pg/P ~ 0.01, so the gas holds ~1 % of the
    // heat capacity, and an explicit gas-side deposit was -0.5 e_gas per step at cfl 0.3
    // (v2/runs/mF03: force 0.94 rho g).  The implicit exchange gives the gas its share.
    Kokkos::realloc(pvet->esrc, nmb1+1, n3m1+1, n2m1+1, n1m1+1);
    Kokkos::deep_copy(pvet->esrc, 0.0);
    auto es = pvet->esrc;
    auto fmd = star_fm_;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &volume = pmbp->pcoord->volume;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    par_for("star_esrc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      es(m,k,j,i) = -(area1(m,k,j,i+1)*fmd(i+1) - area1(m,k,j,i)*fmd(i))/volume(m,k,j,i);
    });
    pvet->esrc_on = true;
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: mlt_flux_frozen: max F_MLT/F on the faces = " << fmax
                << " at r = " << rmx << "; F_MLT at the bottom/top faces = "
                << hfm(indcs.is) << " / " << hfm(indcs.ie + 1) << std::endl;
    }
  }

  // TEST of <rad_vet>/esrc (default off): a constant source star_esrc_const (erg/cm^3/s)
  // in the cells with r < star_esrc_rmax, ON TOP of the frozen MLT deposit if that is on
  const Real esc = pin->GetOrAddReal("problem","star_esrc_const",0.0);
  star_esc_ = esc;
  if (esc != 0.0) {
    const Real rmx = pin->GetOrAddReal("problem","star_esrc_rmax",star_rin_);
    star_esrmx_ = rmx;
    if (!pvet->esrc_on) {
      Kokkos::realloc(pvet->esrc, nmb1+1, n3m1+1, n2m1+1, n1m1+1);
      Kokkos::deep_copy(pvet->esrc, 0.0);
      pvet->esrc_on = true;
    }
    auto es = pvet->esrc;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    par_for("star_esrc_t", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (x1v(m,i) < rmx) es(m,k,j,i) += esc;
    });
  }

  // problem/mlt_ramp_start (s, default -1 = no ramp) and mlt_ramp_time (s): the weight
  // w(t) on the frozen MLT deposit, w = 1 for t < start, 0.5 (1 + cos(pi (t - start)/
  // time)) during the ramp, 0 after.  A function of t only, re-applied on a restart (both
  // keys may be overridden on the restart's command line).  Unused: esrc bitwise.
  star_mw_ = 1.0;
  if (star_mrs_ >= 0.0 && !star_mlt_) {
    StarFatal("problem/mlt_ramp_start >= 0 needs problem/mlt_flux_frozen = true",
              __LINE__);
  }
  // the Picard window trailer of the pgen restart block (stripped here,
  // so the adaptive state below reads the block as before); the rad_vet counters restart
  // at 0, so the window base is minus the open window's sums
  star_mesh_ = pmy_mesh_;
  // problem/star_ic_remap_file (read only when named; absent = nothing,
  // bitwise the old code): a FRESH start whose state is a remapped restart of another
  // run (tools: remap_e1w.py).  The file carries the parent's time, its pgen restart
  // block (Picard window + the star_ic_equil a_ref trailer, applied below exactly as a
  // restart applies them) and u0, the VET u0 and f0x1 of this mesh including ghosts.
  // <time>/start_time must equal the parent's time (bitwise; the generated deck sets it)
  // so that w_mlt(t), the vet_gd rotation clock (physical) and the outputs continue.  No
  // ic_rad, no star_ic_equil, no seed; everything else is a fresh start (hesdirk2 starts
  // with time2_start, the implicit F0 faces from implicit_f0_init, empty rotavg history).
  const bool star_rmp = !restart &&
                        pin->DoesParameterExist("problem","star_ic_remap_file");
  std::string star_rmp_fn;
  std::int64_t star_rmp_data = 0;   // byte offset of the arrays
  int star_rmp_n[4] = {0, 0, 0, 0};  // ng, nx1, nx2, nx3
  if (star_rmp) {
    star_rmp_fn = pin->GetString("problem","star_ic_remap_file");
    std::ifstream rf(star_rmp_fn, std::ios::binary);
    if (!rf) StarFatal("cannot open problem/star_ic_remap_file " + star_rmp_fn, __LINE__);
    char mg[8];
    std::int32_t hi[8];
    double hd[7];
    std::uint64_t npg = 0;
    rf.read(mg, 8);
    rf.read(reinterpret_cast<char*>(hi), sizeof(hi));
    rf.read(reinterpret_cast<char*>(hd), sizeof(hd));
    rf.read(reinterpret_cast<char*>(&npg), sizeof(npg));
    if (!rf || std::memcmp(mg, "HSREMAP1", 8) != 0 || hi[0] != 1 || hi[5] != 5 ||
        hi[6] != 4 || npg > (1u << 30)) {
      StarFatal("problem/star_ic_remap_file: not a version-1 HSREMAP1 file", __LINE__);
    }
    auto &ms = pmy_mesh_->mesh_size;
    auto &mi = pmy_mesh_->mesh_indcs;
    auto rel = [](Real a, Real b) {return std::abs(a - b) <= 1.0e-12*std::abs(b) +
                                          1.0e-14;};
    if (hi[1] != mi.ng || hi[2] != mi.nx1 || hi[3] != mi.nx2 || hi[4] != mi.nx3 ||
        !rel(hd[1], ms.x1min) || !rel(hd[2], ms.x1max) || !rel(hd[3], ms.x2min) ||
        !rel(hd[4], ms.x2max) || !rel(hd[5], ms.x3min) || !rel(hd[6], ms.x3max)) {
      StarFatal("problem/star_ic_remap_file: its mesh (ng, nx1..3, extents) differs from "
                "<mesh>", __LINE__);
    }
    if (hd[0] != pmy_mesh_->time) {
      char buf[128];
      snprintf(buf, sizeof(buf), "%.17g", hd[0]);
      StarFatal(std::string("problem/star_ic_remap_file: set <time>/start_time = ") +
                buf + " (the parent's time; the remap deck from remap_e1w.py does)",
                __LINE__);
    }
    pgen_rststate.resize(npg);
    if (npg > 0) rf.read(pgen_rststate.data(), npg);
    if (!rf) StarFatal("problem/star_ic_remap_file: truncated pgen block", __LINE__);
    star_rmp_data = 8 + sizeof(hi) + sizeof(hd) + sizeof(npg) + npg;
    for (int n=0; n<4; ++n) star_rmp_n[n] = hi[n+1];
    const std::int64_t c1 = hi[2] + 2*hi[1], c2 = hi[3] + 2*hi[1], c3 = hi[4] + 2*hi[1];
    const std::int64_t want = star_rmp_data + 8*(9*c1*c2*c3 + (c1 + 1)*c2*c3);
    rf.seekg(0, std::ios::end);
    if (static_cast<std::int64_t>(rf.tellg()) != want) {
      StarFatal("problem/star_ic_remap_file: file size does not match its header",
                __LINE__);
    }
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: star_ic_remap_file " << star_rmp_fn << ": t = " << hd[0]
                << ", pgen block " << npg << " bytes (remapped fresh start)" << std::endl;
    }
  }
  // Both keys default to auto.  A number is used as given, mlt_ramp_start = -1 = no ramp.
  // auto on a FRESH start (mlt_flux_frozen only, else no ramp): start = 0, time = the
  // convective saturation time t_sat = mlt_ramp_nefold/(mlt_ramp_sigma_frac sigma_max)
  // (defaults 10.5 = ln(v_MLT/v_seed) and 0.61 = sigma_eff/sigma_max, measured on a red
  // supergiant; read only when a ramp key is auto), sigma_max from the IC column
  // (StarN2Compute).  The
  // mass-weighted median 1/sigma (StarN2Median) is printed for reference.  Computed once
  // and kept in the restart file (kStarRmpTag), which a
  // restart reuses (never recomputed).  A restart (or star_ic_remap_file start) with auto
  // and no stored numbers keeps NO ramp (an old restart file carries the old binary's
  // explicit default -1 anyway).
  // the drive's trailer is the LAST one: strip it first (star_drive.cpp)
  if (restart || star_rmp) {stardrive::RestartStrip(&pgen_rststate);}
  {
    const std::string ms = pin->GetOrAddString("problem","mlt_ramp_start","auto");
    const std::string mt = pin->GetOrAddString("problem","mlt_ramp_time","auto");
    const bool sa = (ms == "auto"), ta = (mt == "auto");
    star_mrs_ = sa ? -1.0 : static_cast<Real>(atof(ms.c_str()));
    star_mrt_ = ta ? 0.0 : static_cast<Real>(atof(mt.c_str()));
    Real st[2] = {-1.0, 0.0};
    bool have = false;
    if (restart || star_rmp) {   // strip the trailer (the last one) whatever the keys say
      auto &blk = pgen_rststate;
      const std::size_t nt = sizeof(kStarRmpTag) + 2*sizeof(Real);
      if (blk.size() >= nt &&
          std::memcmp(blk.data() + blk.size() - sizeof(kStarRmpTag), &(kStarRmpTag[0]),
                      sizeof(kStarRmpTag)) == 0) {
        std::memcpy(&(st[0]), blk.data() + blk.size() - nt, 2*sizeof(Real));
        blk.resize(blk.size() - nt);
        have = true;
      }
    }
    const bool on = !(!sa && star_mrs_ < 0.0);   // an explicit negative start: no ramp
    if ((sa || ta) && on) {
      std::string how;
      if (restart || star_rmp) {
        if (have) {
          if (sa) star_mrs_ = st[0];
          if (ta) star_mrt_ = st[1];
          star_mr_store_ = true;
          how = "restart file (computed at the fresh start)";
        } else {
          star_mrs_ = -1.0;
          how = "no stored numbers in the restart file: NO ramp";
        }
      } else if (!star_mlt_) {
        star_mrs_ = -1.0;
        how = "mlt_flux_frozen = false: NO ramp";
      } else {
        StarN2Compute(pmbp, eos, ar);
        if (star_n2_.smax > 0.0) {
          const Real nef = pin->GetOrAddReal("problem","mlt_ramp_nefold",10.5);
          const Real sfr = pin->GetOrAddReal("problem","mlt_ramp_sigma_frac",0.61);
          if (!(nef > 0.0) || !(sfr > 0.0)) {
            StarFatal("problem/mlt_ramp_nefold and mlt_ramp_sigma_frac must be > 0",
                      __LINE__);
          }
          StarN2Median(pmbp, eos, ar);
          if (sa) star_mrs_ = 0.0;
          if (ta) star_mrt_ = nef/(sfr*star_n2_.smax);
          if (global_variable::my_rank == 0) {
            std::cout << "star_vet: mlt_ramp auto rule v2: t_sat = n_efold " << nef
                      << " / (sigma_frac " << sfr << " sigma_max) = "
                      << nef/(sfr*star_n2_.smax) << " s" << std::endl;
          }
          star_mr_store_ = true;
          how = "IC column N^2";
        } else {
          star_mrs_ = -1.0;
          how = "no unstable cell below r_int: NO ramp";
        }
      }
      if (global_variable::my_rank == 0) {
        std::cout << "star_vet: mlt_ramp auto (start " << ms << ", time " << mt
                  << ") from " << how << ": start = " << star_mrs_ << " s ("
                  << star_mrs_/86400.0 << " d), time = " << star_mrt_ << " s ("
                  << star_mrt_/86400.0 << " d)";
        if (star_n2_.done && star_n2_.smax > 0.0) {
          std::cout << "; sigma_max " << star_n2_.smax << " /s (1/sigma "
                    << 1.0/star_n2_.smax/86400.0 << " d) at r/r_int "
                    << star_n2_.r[star_n2_.ipk]/star_rint_ << ", unstable cells r/r_int "
                    << star_n2_.r[star_n2_.iu0]/star_rint_ << " - "
                    << star_n2_.r[star_n2_.iu1]/star_rint_ << " (" << star_n2_.nu
                    << " cells), mass-weighted median 1/sigma " << star_n2_.med/86400.0
                    << " d";
        }
        std::cout << std::endl;
      }
    }
  }
  // (star_ic_equil): the equilibrium's a_ref of the active cells travels in the
  // pgen restart block as the LAST tagged trailer (a_ref..., Real n, tag); a restart
  // re-applies it with the Phi_eff rebuilt from it above the face x1f(is), exactly as the
  // fresh start's StarIcEquil left them (stripped here, before the Picard trailer)
  if (restart) {
    // the removed fixed-profile inner face's wall-luminosity trailer (the last one) is
    // dropped
    auto &blk = pgen_rststate;
    const std::size_t nt = sizeof(kStarLbcTag) + sizeof(Real);
    if (blk.size() >= nt &&
        std::memcmp(blk.data() + blk.size() - sizeof(kStarLbcTag), &(kStarLbcTag[0]),
                    sizeof(kStarLbcTag)) == 0) {
      blk.resize(blk.size() - nt);
    }
  }
  // (also from a star_ic_remap_file's copy of the parent's block)
  if (restart || star_rmp) {
    auto &blk = pgen_rststate;
    const std::size_t nt = sizeof(kStarDicTag) + sizeof(Real);
    if (blk.size() >= nt &&
        std::memcmp(blk.data() + blk.size() - sizeof(kStarDicTag), &(kStarDicTag[0]),
                    sizeof(kStarDicTag)) == 0) {
      Real rn;
      std::memcpy(&rn, blk.data() + blk.size() - nt, sizeof(Real));
      const int na = static_cast<int>(rn);
      const std::size_t nb = nt + static_cast<std::size_t>(na)*sizeof(Real);
      if (na != indcs.nx1 || blk.size() < nb) {
        StarFatal("the restart file's star_ic_equil trailer does not match this radial "
                  "grid", __LINE__);
      }
      star_eq_a_.assign(na, 0.0);
      std::memcpy(star_eq_a_.data(), blk.data() + blk.size() - nb, na*sizeof(Real));
      blk.resize(blk.size() - nb);
      std::vector<Real> pc0, pf0, p0c, p0f, x1v, x1f, ar0, pc, pf;
      StarAtmRadial(pmbp, pc0, pf0, p0c, p0f, x1v, x1f, ar0);
      StarAtmPhiEff(indcs.is, indcs.ie, star_eq_a_, pc0, pf0, p0c, p0f, x1v, x1f, pc, pf);
      StarAtmWrite(pmbp, indcs.is, star_eq_a_, pc, pf);
      if (global_variable::my_rank == 0) {
        std::cout
            << "star_vet: star_ic_equil: a_ref (cells is .. ie) and Phi_eff from the "
                  << "fresh start's equilibrium (restart trailer)" << std::endl;
      }
    } else if (star_eq_) {
      StarFatal("problem/star_ic_equil = true but the restart file carries no "
                "star_ic_equil trailer", __LINE__);
    }
  }
  if (restart || star_rmp) {
    auto &blk = pgen_rststate;
    const std::size_t nt = sizeof(kStarPicTag) + 2*sizeof(Real);
    if (blk.size() >= nt &&
        std::memcmp(blk.data() + blk.size() - nt, &(kStarPicTag[0]),
                    sizeof(kStarPicTag)) == 0) {
      Real w[2];
      std::memcpy(&(w[0]), blk.data() + blk.size() - 2*sizeof(Real), 2*sizeof(Real));
      // a remapped start opens a new window (the trailer is only stripped)
      if (restart) {
        star_pic_n0_ = -w[0];
        star_pic_s0_ = -w[1];
      }
      blk.resize(blk.size() - nt);
    }
    // an adaptive-closure remap keeps the adaptive state at the block's
    // front (read below as a restart reads it); anything else left is refused
    const bool rmp_ad = star_rmp && (pin->GetOrAddString("problem","mlt_closure",
                                                         "frozen") != "frozen");
    if (star_rmp && !blk.empty() && !rmp_ad) {
      StarFatal("problem/star_ic_remap_file: the parent's pgen block carries state this "
                "start cannot take (adaptive MLT?)", __LINE__);
    }
  }
  // problem/mlt_closure: frozen (default: the IC's F_MLT, times w(t)) or adaptive (see
  // star_ad_ above; w(t) caps it as an overall factor)
  {
    const std::string mc = pin->GetOrAddString("problem","mlt_closure","frozen");
    if (mc != "frozen" && mc != "adaptive" && mc != "pressure") {
      StarFatal("problem/mlt_closure must be frozen, adaptive or pressure", __LINE__);
    }
    star_ad_ = (mc == "adaptive" || mc == "pressure");
    star_adp_ = (mc == "pressure");
    if (star_mlt_) {
      recipe::Fix(pin, "problem", "mlt_src_time", "stage", "start");
      recipe::OldStored(pin, "problem", "mlt_src_time", "start", "stage");
    }
    const std::string tc = pin->GetOrAddString("problem","mlt_src_time","stage");
    if (tc != "start" && tc != "stage") {
      StarFatal("problem/mlt_src_time must be start|stage", __LINE__);
    }
    star_mtc_ = (tc == "stage");
  }
  // problem/star_ad_dead_skip (default true; acts only under mesh/layout_invariant)
  star_ad_dskip_ = true;
  if (pin->DoesParameterExist("problem","star_ad_dead_skip")) {
    star_ad_dskip_ = pin->GetBoolean("problem","star_ad_dead_skip");
  }
  if (star_ad_) {
    if (!star_mlt_) {
      StarFatal("problem/mlt_closure = adaptive needs problem/mlt_flux_frozen = true "
                "(the IC's F_MLT is its start profile and its FeCZ mask)", __LINE__);
    }
    if (pvet->f0x1.extent_int(0) <= 0) {
      StarFatal("problem/mlt_closure = adaptive needs the implicit VET face flux f0x1",
                __LINE__);
    }
    star_ad_tau_ = pin->GetOrAddReal("problem","mlt_relax_time",1500.0);
    star_ad_dout_ = pin->GetOrAddReal("problem","mlt_closure_dout",0.0);
    if (!(star_ad_tau_ > 0.0)) {
      StarFatal("problem/mlt_relax_time must be > 0", __LINE__);
    }
    const int nfc = n1m1 + 2;
    Kokkos::realloc(star_fsub_, nfc);
    Kokkos::deep_copy(star_fsub_, star_fm_);
    if (star_mtc_) Kokkos::realloc(star_tgp_, nfc);
    Kokkos::realloc(star_shs_, nfc, 8);
    if (star_adp_) {
      // the IC's F_MLT/F against ln P_tot (P_gas(rho, T) + E/3 of the start-up column,
      // the fine grid between x1min and x1max), and F_tot = F_r/(1 - fm) on the faces
      if (!tcol_read) {
        StarFatal("problem/mlt_closure = pressure needs the 7-column file", __LINE__);
      }
      std::vector<std::pair<Real, Real>> pf;
      for (int n=0; n<nf; ++n) {
        if (hr[n] < star_rin_ || hr[n] > rtop) continue;
        Real ee, pp, cr, ct, cv;
        eos.ThermoAt(hd[n], hTcol[n]/star_tk, ee, pp, cr, ct, cv);
        const Real er = (hE[n] > 0.0) ? hE[n] : ar*SQR(SQR(hTcol[n]/star_tk));
        pf.emplace_back(log(pp + er/3.0), hM[n]);
      }
      std::sort(pf.begin(), pf.end());
      star_plnp_.clear();
      star_pfm_.clear();
      for (const auto &q : pf) {
        if (!star_plnp_.empty() && !(q.first > star_plnp_.back())) continue;
        star_plnp_.push_back(q.first);
        star_pfm_.push_back(q.second);
      }
      auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
      star_ftot_.assign(nfc, 0.0);
      for (int i=0; i<nfc; ++i) {
        const Real r = hx1f(0,i);
        Real x = (r - star_rlo_)/star_dr_;
        int n = static_cast<int>(floor(x));
        n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
        Real w = x - n;
        w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
        const Real fm = (1.0 - w)*hM[n] + w*hM[n+1];
        const Real fr = exp((1.0 - w)*log(hF[n]) + w*log(hF[n+1]));
        star_ftot_[i] = (fm < 1.0) ? fr/(1.0 - fm) : fr;
      }
      Kokkos::realloc(star_shp_, nfc, 3);
      if (star_plnp_.size() < 2) {
        StarFatal("problem/mlt_closure = pressure: no IC pressure table", __LINE__);
      }
    }
    star_ad_tl_ = pmy_mesh_->time;
    star_ad_cyc_ = -1;
    // a star_ic_remap_file start takes F_sub from the file's adaptive state
    // (the parent's), and its update clock restarts at this start's time
    if (restart || star_rmp) {
      const auto &blk = pgen_rststate;
      if (blk.size() < 2*sizeof(std::int32_t) + sizeof(Real)) {
        if (global_variable::my_rank == 0) {
          std::cout << "### star_vet: this restart file carries no adaptive MLT "
                       "state; "
                    << "F_sub re-seeded with the IC's F_MLT" << std::endl;
        }
      } else {
        std::int32_t hdr[2];
        std::memcpy(&(hdr[0]), blk.data(), sizeof(hdr));
        const int ntg = (hdr[1] == 1) ? nfc : 0;   // tgt^{n-1} (mlt_src_time)
        if (hdr[0] != nfc || (hdr[1] != 0 && hdr[1] != 1) ||
            blk.size() != sizeof(hdr) + (nfc + 1 + ntg)*sizeof(Real)) {
          StarFatal("the restart file's adaptive MLT profile does not match this radial "
                    "grid", __LINE__);
        }
        std::vector<Real> tmp(nfc + 1 + ntg);
        std::memcpy(tmp.data(), blk.data() + sizeof(hdr), tmp.size()*sizeof(Real));
        star_ad_tl_ = star_rmp ? pmy_mesh_->time : tmp[0];
        auto hf = Kokkos::create_mirror_view(star_fsub_);
        for (int i=0; i<nfc; ++i) hf(i) = tmp[i+1];
        Kokkos::deep_copy(star_fsub_, hf);
        if (star_mtc_ && ntg > 0 && !star_rmp) {
          auto hp = Kokkos::create_mirror_view(star_tgp_);
          for (int i=0; i<nfc; ++i) hp(i) = tmp[nfc+1+i];
          Kokkos::deep_copy(star_tgp_, hp);
          star_tgp_ok_ = 1;
        }
      }
      if (star_mtc_ && star_tgp_ok_ == 0 && global_variable::my_rank == 0) {
        std::cout << "star_vet: NOTE mlt_src_time = stage: no tgt^{n-1} in the "
                     "restart "
                  << "file, the first update takes tgt^n for both ends" << std::endl;
      }
    }
    {
      auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_fsub_);
      star_lsub_ = 0.0;
      for (int i=0; i<nfc; ++i) star_lsub_ = std::max(star_lsub_, hf(i)*star_r2o_[i]);
    }
    if (star_ad_dout_ > 0.0) {
      snprintf(star_ad_file_, sizeof(star_ad_file_), "%s.fsub.txt",
               pin->GetString("job","basename").c_str());
      star_ad_tout_ = star_ad_dout_*floor(pmy_mesh_->time/star_ad_dout_ + 1.0);
      if (!restart) star_ad_tout_ = 0.0;
    }
    // esrc from F_sub (bitwise the start-up esrc when F_sub = F_MLT and w = 1)
    StarApplyMltW(pmy_mesh_, StarMltW(pmy_mesh_->time), 0.0);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: mlt_closure = " << (star_adp_ ? "pressure" : "adaptive")
                << ", relax time " << star_ad_tau_
                << " s, every cycle, taper "
                << pin->GetReal("problem","mlt_closure_taper") << ", update clock "
                << star_ad_tl_ << ", max L_sub = " << star_lsub_ << std::endl;
    }
  }
  if (star_mlt_) {
    const Real w = StarMltW(pmy_mesh_->time);
    if (w != star_mw_) StarApplyMltW(pmy_mesh_, w, 0.0);
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: mlt ramp start = " << star_mrs_ << ", time = " << star_mrt_
                << ", w(t = " << pmy_mesh_->time << ") = " << star_mw_
                << ", max L_MLT (unscaled) = " << star_lmx_ << std::endl;
    }
  }

  // top sponge (default off) and hooks
  star_zflux_ = pin->GetOrAddBoolean("problem","star_wall_zero_flux",true);
  // default true; with the window keys absent it is bitwise the plain branch
  star_wbwin_ = true;
  if (pin->DoesParameterExist("problem","star_wb_grav_window")) {
    star_wbwin_ = pin->GetBoolean("problem","star_wb_grav_window");
  }
  // the inner wall face sees the exact mirror of the interior-side state
  // (hydro_fluxes.cpp, wall_closed_ix1, as deep_hot_jupiter_rt problem/wall_closed): zero
  // mass flux to round-off.  Measured: it does not change the first-cell force
  // (-2.1e-3 rho g at 640 zones, first order in dx), which is not the wall flux.
  // problem/star_bc_inner: wall (default: the closed wall above) | inflow:
  // the ghosts hold the (balanced) initial column at fixed rho and eint, UNSCALED, with
  // v_r copied from the first active cell (zero gradient), so mass may enter or leave;
  // no wall-face correction.  problem/star_bc_outer: noinflow (default) | outflow: the
  // ghosts continue the edge cell as a constant-velocity wind, rho ~ r^-2, the edge's
  // specific internal energy (the no-inflow face correction of star_wall_zero_flux stays)
  // | hse: the ghosts are the edge cell walked outward in hydrostatic balance
  // with the effective gravity of the code's own forces (see StarBC), T held (or
  // continued along the interior's dT/dPhi_eff with star_bc_hse_tgrad), v_r copied where
  // it points out and zero where it points in, v_t copied (box_convection bc_mode_top 4).
  {
    const std::string bi = pin->GetOrAddString("problem","star_bc_inner","wall");
    const std::string bo = pin->GetOrAddString("problem","star_bc_outer","noinflow");
    if ((bi != "wall" && bi != "inflow" && bi != "zerocontact") ||
        (bo != "noinflow" && bo != "outflow" && bo != "hse" && bo != "copy")) {
      StarFatal("problem/star_bc_inner must be wall|inflow|zerocontact, star_bc_outer "
                "noinflow|outflow|hse|copy", __LINE__);
    }
    star_binf_ = (bi == "inflow");
    star_bzc_ = (bi == "zerocontact");
    star_bout_ = (bo == "outflow");
    star_bhse_ = (bo == "hse");
    star_bcopy_ = (bo == "copy");
    if (pin->DoesParameterExist("problem","star_bc_outer_r2")) {
      star_bcr2_ = pin->GetBoolean("problem","star_bc_outer_r2");
      if (star_bcr2_ && !star_bcopy_) {
        StarFatal("problem/star_bc_outer_r2 = true needs star_bc_outer = copy", __LINE__);
      }
    }
    if (star_bzc_) {
      // the root find replays the plain (non-WB) spherical-polar
      // PLM (Hydro::PLM_nonuniform) and the ideal-gas HLLC contact speed
      auto *phy = pmbp->phydro;
      if (phy->use_wellbalance_dynamic || phy->rsolver_method != Hydro_RSolver::hllc ||
          !pmy_mesh_->use_spherical_polar || eos.IsGeneral() || star_bal_ ||
          indcs.nx1 < 2 || ng < 2) {
        StarFatal("problem/star_bc_inner = zerocontact needs spherical polar, an ideal "
                  "gas, <hydro>/rsolver = hllc, no well-balanced scheme, "
                  "star_ic_balance = false, nx1 >= 2 and nghost >= 2", __LINE__);
      }
    }
    // problem/star_bc_linhao (the ghost-built fixed-profile inner face of the old solver)
    // is removed
    if (pin->DoesParameterExist("problem","star_bc_linhao") &&
        pin->GetBoolean("problem","star_bc_linhao")) {
      StarFatal("problem/star_bc_linhao was removed with the old solver", __LINE__);
    }
    if (star_bhse_) {
      star_bhse_gmax_ = pin->GetOrAddReal("problem","star_bc_hse_gmax",0.9);
      star_bhse_tg_ = pin->GetOrAddBoolean("problem","star_bc_hse_tgrad",false);
      std::string hf;
      if (restart && !pin->DoesParameterExist("problem","star_bc_hse_flux")) {
        hf = "cell";
        pin->SetString("problem","star_bc_hse_flux",hf);
        if (global_variable::my_rank == 0) {
          std::cout << "### star_vet: restart input has no problem/star_bc_hse_flux; "
                    << "keeping the old cell-flux Gamma (star_bc_hse_flux = cell)"
                    << std::endl;
        }
      } else {
        hf = pin->GetOrAddString("problem","star_bc_hse_flux","face");
      }
      if (hf != "face" && hf != "cell") {
        StarFatal("problem/star_bc_hse_flux must be face|cell", __LINE__);
      }
      star_bhse_face_ = (hf == "face");
      if (star_bhse_face_ && (pmbp->pradvet == nullptr ||
                              pmbp->pradvet->f0x1.extent_int(0) == 0)) {
        StarFatal("problem/star_bc_hse_flux = face needs the implicit rad_vet transport "
                  "(face fluxes f0x1); use cell", __LINE__);
      }
      if (!(star_bhse_gmax_ >= 0.0 && star_bhse_gmax_ < 1.0)) {
        StarFatal("problem/star_bc_hse_gmax must be in [0,1)", __LINE__);
      }
    }
  }
  pmbp->phydro->wall_closed_ix1 = pin->GetOrAddBoolean("problem","star_wall_closed",
                                                       !(star_binf_ || star_bzc_));
  if ((star_binf_ || star_bzc_) && pmbp->phydro->wall_closed_ix1) {
    StarFatal("problem/star_bc_inner = inflow|zerocontact needs star_wall_closed = false",
              __LINE__);
  }
  star_sp_rate_ = pin->GetOrAddReal("problem","star_sponge_rate",0.0);
  star_sp_r0_ = pin->GetOrAddReal("problem","star_sponge_r0",star_rint_);
  if (pin->DoesParameterExist("problem","star_sponge_mode")) {
    const std::string sm = pin->GetString("problem","star_sponge_mode");
    if (sm != "all" && sm != "radial") {
      StarFatal("problem/star_sponge_mode must be all|radial", __LINE__);
    }
    star_sp_rad_ = (sm == "radial");
  }
  if (pin->DoesParameterExist("problem","star_sponge_dmax")) {
    star_sp_dmax_ = pin->GetReal("problem","star_sponge_dmax");
  }
  {StarOldDefault(pin, "star_sponge_form", "exp");
  const std::string sf = pin->GetOrAddString("problem","star_sponge_form",
                                  global_variable::restart_run ? "exp" : "lawson");
  if (sf != "exp" && sf != "lawson") {
    StarFatal("problem/star_sponge_form must be exp|lawson", __LINE__);
  }
  star_sp_law_ = (sf == "lawson");}
  star_rtop_ = rtop;
  user_srcs_func = StarGravity;
  if (star_zflux_) {
    user_hydro_flux_func = StarWallFlux;
    Kokkos::realloc(star_min_, pmbp->nmb_thispack, indcs.nx3 + 2*indcs.ng,
                    indcs.nx2 + 2*indcs.ng);
  }
  // fin51-r-fix: the per-cycle state of StarGravity joins the <time>/dt_retry
  // checkpoint (a retried cycle redoes the ramp weight / adaptive closure update exactly
  // as the straight cycle; esrc itself is in the m1 minimal list)
  dtretry::RegisterValue("star_mw", &star_mw_);
  dtretry::RegisterValue("star_ad_cyc", &star_ad_cyc_);
  dtretry::RegisterValue("star_ad_tl", &star_ad_tl_);
  dtretry::RegisterValue("star_ad_tout", &star_ad_tout_);
  dtretry::RegisterValue("star_lsub", &star_lsub_);
  dtretry::RegisterView("star_fsub", &star_fsub_);
  if (star_mtc_ && star_ad_) {
    dtretry::RegisterView("star_tgp", &star_tgp_);
    dtretry::RegisterValue("star_tgp_ok", &star_tgp_ok_);
  }
  user_bcs_func = StarBC;
  user_hist_func = StarHist;
  // (read only when named): problem/am_hist = true takes the history to
  // the energy / angular-momentum budget RadVetAmHistory (pgen/tests/rad_vet_tests2.cpp:
  // Eg, Er, gas L, radiation L, am_vector spin, frame L) instead
  if (pin->DoesParameterExist("problem", "am_hist") &&
      pin->GetBoolean("problem", "am_hist")) {
    user_hist_func = RadVetAmHistory;
  }
  pgen_rst_write_func = StarAdaptiveRstWrite;   // adaptive state + the Picard window
  pgen_final_func = StarFinal;
  // problem/star_ic_rad_atm_post (default true with star_ic_rad_atm_force =
  // module; a fresh start whose <rad_vet>/ic_rad runs; the x1 WB pair): after the ic_rad
  // passes the thin atmosphere's gas sits at the code's own radiative-equilibrium T
  // (0.2-1.8 % in e/rho above the pgen's sph1d RE T, AG Car A; the same shift happens in
  // the first ~10 cycles with ic_rad = none), which leaves it out of hydrostatic balance
  // by ~0.5 % g.  StarIcPost re-marches rho of the atmosphere cells outward from the
  // first cell above R_ph in the DISCRETE balance of the x1 WB pair (the face condition
  // of star_ic_balance, walk of the code's own WBBackgroundStencil along Phi_eff) at the
  // relaxed T (shell mean; each cell keeps its own T), until rho reaches the floor.
  // a_ref and Phi_eff are NOT changed (restart-exact: a restart rebuilds them from the
  // input); the cells below R_ph are not touched.
  {
    const bool atp = pin->GetOrAddBoolean("problem","star_ic_rad_atm_post",hatm_on);
    if (atp && !hatm_on) {
      StarFatal("problem/star_ic_rad_atm_post = true needs star_ic_rad_atm_force = "
                "module", __LINE__);
    }
    star_atp_ = atp && !global_variable::restart_run && wbdyn && !star_pm_ && !star_rmp;
    if (star_atp_) user_icpost_func = StarIcPost;
  }
  // problem/star_ic_equil (read only when named; default false = nothing, bitwise
  // the old start): the DISCRETE EQUILIBRIUM of the initial state on the code's own
  // operators, iterated with <rad_vet>/ic_rad (StarIcEquil, called after every ic_rad
  // run; Driver::Initialize repeats ic_rad + hook while user_icpost_again).  Each round:
  // (1) ic_rad relaxes E, F and the gas T at fixed rho (the relax steps of
  //     star_ic_equil_dtfac x dt: the steady state of the transport, the frozen MLT
  //     scaffold (w = 1) and the gas-radiation exchange, i.e. L_rad + L_MLT = L_in at
  //     every face), then ONE pass with the radiative force on (rad_vet icr_fmeas) gives
  //     the module's own net momentum deposit dm - dt rho a_ref of every cell;
  // (2) the force per volume f = dm/dt + rho a_ref (radiation pressure gradient form,
  //     independent of rho at fixed E) gives a_ref = f/rho; Phi_eff is rebuilt above the
  //     face x1f(is) from these a_ref (StarAtmPhiEff) and rho is re-marched at the
  //     relaxed T in the discrete balance of the x1 WB pair (the star_ic_balance face
  //     condition) OUTWARD from the wall cell (its rho kept) to the floor, a_i =
  //     f_i/rho_i solved with rho_i in each face solve; e from (rho, T).
  // It stops (no update) when both residuals -- the measured module deposit and the WB
  // face condition, as accelerations / g -- are below star_ic_equil_tol, or after
  // star_ic_equil_maxit rounds.  The star_seed factor is applied after it (deferred).
  // a_ref travels in the restart files (kStarDicTag).  Needs star_ic_balance, ic_rad, the
  // WB pair, a fresh start; replaces star_ic_rad_atm_post.
  if (star_eq_ && !restart && !star_rmp) {
    if (!star_bal_ || star_pm_ || (pvet->icr_mode == 0 && !pvet->cln_steady) || star_atp_)
        {
      StarFatal("problem/star_ic_equil = true needs star_ic_balance = true, star_gravity "
                "= wb, <rad_vet>/ic_rad != none (scheme = clean: steady_init = steady) "
                "and star_ic_rad_atm_post = false", __LINE__);
    }
    // (scheme = clean: each round's relax is SteadyStart, E and the face F at the
    // scheme's own fixed point with the gas T, and its force pass the last steady pass;
    // the x1 shells of the force, its radial component: spherical polar or the cubed
    // sphere)
    if (!(pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere)) {
      StarFatal("problem/star_ic_equil with scheme = clean needs spherical polar or the "
                "cubed sphere", __LINE__);
    }
    if (wic) StarFatal("problem/star_ic_equil does not combine with star_wind_ic",
                       __LINE__);
    star_eq_maxit_ = pin->GetOrAddInteger("problem","star_ic_equil_maxit",30);
    // star_ic_equil_last_march (default false = the last round only checks): the round
    // star_ic_equil_maxit still updates (march at the relaxed T, a_i = f_i/rho_i) and the
    // equilibrium stops there WITHOUT a further ic_rad relax: E, F, T are the relaxed
    // ones of that round, rho/a_ref/Phi_eff the march's (module deposit and WB faces
    // balanced at fixed E; the face luminosities then differ by the rho change)
    star_eq_last_ = pin->GetOrAddBoolean("problem","star_ic_equil_last_march",false);
    star_eq_tol_ = pin->GetOrAddReal("problem","star_ic_equil_tol",1.0e-9);
    const Real eqdt = pin->GetOrAddReal("problem","star_ic_equil_dtfac",1.0e6);
    const Real eqrt = pin->GetOrAddReal("problem","star_ic_equil_relax_tol",1.0e-12);
    if (star_eq_maxit_ < 1 || !(star_eq_tol_ > 0.0) || !(eqdt > 0.0) || !(eqrt > 0.0)) {
      StarFatal("problem/star_ic_equil_maxit >= 1, _tol, _dtfac, _relax_tol > 0",
                __LINE__);
    }
    pvet->icr_fmeas = true;
    star_eq_it_ = 0;
    user_icpost_func = StarIcEquil;
    user_icpost_again = false;
  }
  if (global_variable::my_rank == 0) {
    const int i0 = static_cast<int>((star_rin_ - star_rlo_)/star_dr_);
    std::cout << "star_vet: GM = " << star_gm_ << ", F_in = " << star_fin_
              << ", tau(r_in) = " << tau[i0] << ", r_int (tau = " << tau_int
              << ") = " << star_rint_ << ", r_top = " << rtop
              << ", sponge rate = " << star_sp_rate_ << " above r = " << star_sp_r0_
              << std::endl;
  }
  // (star_drive.hpp; every key read only when named, absent = bitwise):
  // the inflow table, the live monopole gravity and the boundary ledgers
  {
    stardrive::SetupArgs da;
    if (pmbp->punit != nullptr) {
      const Real lu = pmbp->punit->length_cgs(), tu = pmbp->punit->time_cgs();
      da.gcode = units::Units::grav_constant_cgs*pmbp->punit->mass_cgs()*tu*tu/
                 (lu*lu*lu);
    }
    da.tk = star_tk;
    da.gm = star_gm_;
    da.fconst = star_fin_;
    da.binf = star_binf_;
    da.pointmass = star_pm_;
    da.bc_hse = star_bhse_;
    if (pin->DoesParameterExist("problem","star_bc_inner_table")) {
      DvceArray1D<Real> c2("star_drv_c", 2);
      auto crho = star_rho_, ceint = star_eint_;
      const Real rlo = star_rlo_, dr = star_dr_, rin = star_rin_;
      const int nf2 = star_nf_;
      par_for("star_drv_col", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int n) {
        c2(0) = StarLogInterp(crho, rlo, dr, nf2, rin);
        c2(1) = StarLogInterp(ceint, rlo, dr, nf2, rin);
      });
      auto hc2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), c2);
      da.rho_rin = hc2(0);
      da.eint_rin = hc2(1);
      star_drv_d0_ = hc2(0);
      star_drv_e0_ = hc2(1);
    }
    stardrive::Setup(pin, pmy_mesh_, restart, da);
    if (stardrive::On()) {user_cycle_func = stardrive::Cycle;}
  }
  // (star_strip_bg.cpp): one-shot restart-time background strip, off unless
  // problem/strip_bg_on_restart = true (absent: bitwise)
  if (restart) {
    stardrive::RestartApply(pmy_mesh_);   // (absent: nothing)
    StarEosReinitFromT(pin, pmy_mesh_);   // (absent: bitwise)
    StarStripBackground(pin, pmy_mesh_);
    StarRestartPerturb(pin, pmy_mesh_);   // (absent: bitwise)
    return;
  }
  // a fresh start builds e with this run's EOS table already, so
  // eos_reinit_from_t has nothing to do; mark it done so that the restarts of this run
  // never fire it (see StarEosReinitFromT).  The state is unchanged.
  if (pin->DoesParameterExist("problem", "eos_reinit_from_t") &&
      pin->GetBoolean("problem", "eos_reinit_from_t")) {
    pin->SetBoolean("problem", "eos_reinit_from_t", false);
    pin->SetBoolean("problem", "eos_reinit_done", true);
  }
  if (star_rmp) {
    // no <rad_vet>/ic_rad (the state is an evolved one), then the arrays
    pvet->icr_mode = 0;
    StarRemapLoad(pmy_mesh_, star_rmp_fn, star_rmp_data, star_rmp_n);
    return;
  }

  // ---- the initial state (active + ghost cells)
  // problem/star_seed_band = n2 (default) | old.
  // old = the seed exactly as the keys below name it (star_seed default 0 = no seed).  n2
  // = the radial band, the Y_lm l range and (star_seed default 1e-3) the amplitude come
  // from the convective growth rate of the initial column (StarN2Band);
  // star_seed_rlo/rhi, _ylm/_lmin/_lmax are then ignored, star_seed_rad/_rng are used as
  // named.  Fresh starts only (restarts and remaps returned above).
  const std::string sband = pin->GetOrAddString("problem","star_seed_band","n2");
  if (sband != "n2" && sband != "old") {
    StarFatal("star_seed_band must be n2 or old", __LINE__);
  }
  const bool sn2 = (sband == "n2");
  Real seed = pin->GetOrAddReal("problem","star_seed",sn2 ? 1.0e-3 : 0.0);
  if (hlin && seed != 0.0) {
    StarFatal("problem/star_ic_profile = true replaces the star_seed perturbation by its "
              "own (set star_seed = 0)", __LINE__);
  }
  const Real seedk = pin->GetOrAddReal("problem","star_seed_k",4.0);
  Real srlo = pin->GetOrAddReal("problem","star_seed_rlo",star_rin_);
  Real srhi = pin->GetOrAddReal("problem","star_seed_rhi",star_rint_);
  int n2lmin = 1, n2lmax = 1;
  if (sn2 && seed != 0.0) {
    const Real n2fr = pin->GetOrAddReal("problem","star_seed_n2_frac",0.3);
    const int n2nm = pin->GetOrAddInteger("problem","star_seed_n2_nmin",4);
    const int n2gp = pin->GetOrAddInteger("problem","star_seed_n2_gap",2);
    const Real n2pw = pin->GetOrAddReal("problem","star_seed_n2_ppw",8.0);
    if (!(n2fr > 0.0 && n2fr < 1.0) || n2nm < 1 || n2gp < 0 || !(n2pw >= 2.0)) {
      StarFatal("star_seed_n2_frac in (0,1), _nmin >= 1, _gap >= 0, _ppw >= 2", __LINE__);
    }
    if (!StarN2Band(pmbp, eos, ar, star_tk, hcs, n2fr, n2nm, n2gp, n2pw, 64, srlo, srhi,
                    n2lmin, n2lmax)) {
      seed = 0.0;   // no unstable band: no seed (said in the log)
    }
  }
  // multi-mode seed (star_seed_nk > 0, default 0 = the single mode above): nk random
  // lateral modes sin(2 pi (k2 (theta - theta_a)/L_theta + k3 (phi - phi_a)/L_phi) + ph),
  // integer k2, k3 drawn uniformly in [star_seed_kmin, star_seed_kmax] (periodic in the
  // wedge), amplitudes normalised to sum a_n^2 = 1 (field rms 1/sqrt 2), under the radial
  // envelope sin(pi (r - rlo)/(rhi - rlo)); box_convection's vpert_var = eint recipe.
  // star_seed_rad = true scales E by (1 + d)^4 with eint by (1 + d): a TEMPERATURE seed
  // at fixed rho, gas and radiation in equilibrium (a gas-only eint seed in a radiation-
  // dominated layer is removed by the stiff exchange in the first step).
  const int snk = pin->GetOrAddInteger("problem","star_seed_nk",0);
  const bool srad = pin->GetOrAddBoolean("problem","star_seed_rad",false);
  const int skmin = pin->GetOrAddInteger("problem","star_seed_kmin",1);
  const int skmax = pin->GetOrAddInteger("problem","star_seed_kmax",2);
  const int srng = pin->GetOrAddInteger("problem","star_seed_rng",1234);
  // star_seed_signs = positive (default, the historical draw: k2, k3 >= 0, so every
  // wavefront tilts the same diagonal way) | random (an independent random sign for k2
  // and for k3 of each mode).  star_seed_kdist = box (default: k2, k3 independently
  // uniform in [kmin, kmax]) | shell (|k| = sqrt(k2^2 + k3^2) uniform in [kmin, kmax],
  // direction uniform, both rounded to integers, k = 0 rejected; with signs = positive
  // the direction is limited to the first quadrant).  Both extra draws come AFTER the
  // historical per-mode draws (a separate pass on the same RNG stream), so the default
  // keys reproduce the old seed bitwise.
  const std::string ssign = pin->GetOrAddString("problem","star_seed_signs","positive");
  const std::string skdist = pin->GetOrAddString("problem","star_seed_kdist","box");
  if (ssign != "positive" && ssign != "random") {
    StarFatal("star_seed_signs must be positive or random", __LINE__);
  }
  if (skdist != "box" && skdist != "shell") {
    StarFatal("star_seed_kdist must be box or shell", __LINE__);
  }
  const bool srsign = (ssign == "random"), sshell = (skdist == "shell");
  if (snk > 0 && (skmin < 0 || skmax < skmin || !(srhi > srlo))) {
    StarFatal("need 0 <= star_seed_kmin <= star_seed_kmax and star_seed_rhi > "
              "star_seed_rlo", __LINE__);
  }
  if (snk > 0 && sshell && skmax < 1) {
    StarFatal("star_seed_kdist = shell needs star_seed_kmax >= 1", __LINE__);
  }
  DualArray2D<Real> smd("star_seed_modes", std::max(snk, 1), 4);   // k2, k3, amp, phase
  {
    std::mt19937 rng(srng);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const Real nkr = static_cast<Real>(skmax - skmin + 1);
    Real norm = 0.0;
    for (int n=0; n<snk; ++n) {
      smd.h_view(n,0) = skmin + std::floor(nkr*u01(rng));
      smd.h_view(n,1) = skmin + std::floor(nkr*u01(rng));
      smd.h_view(n,2) = u01(rng) + 0.25;
      smd.h_view(n,3) = 2.0*M_PI*u01(rng);
      norm += SQR(smd.h_view(n,2));
    }
    norm = (norm > 0.0) ? 1.0/std::sqrt(norm) : 1.0;
    for (int n=0; n<snk; ++n) smd.h_view(n,2) *= norm;
    // opt-in wave-vector redraws (separate pass: the default path draws nothing here)
    if (sshell) {
      const Real amax = srsign ? 2.0*M_PI : 0.5*M_PI;
      for (int n=0; n<snk; ++n) {
        Real k2 = 0.0, k3 = 0.0;
        while (k2 == 0.0 && k3 == 0.0) {
          const Real kk = skmin + (skmax - skmin)*u01(rng);
          const Real an = amax*u01(rng);
          k2 = std::round(kk*std::cos(an));
          k3 = std::round(kk*std::sin(an));
        }
        smd.h_view(n,0) = k2;
        smd.h_view(n,1) = k3;
      }
    } else if (srsign) {
      for (int n=0; n<snk; ++n) {
        if (u01(rng) < 0.5) smd.h_view(n,0) = -smd.h_view(n,0);
        if (u01(rng) < 0.5) smd.h_view(n,1) = -smd.h_view(n,1);
      }
    }
    if (snk > 0 && global_variable::my_rank == 0) {
      std::cout << "star_vet: seed " << snk << " modes, amplitude " << seed
                << (srad ? " (temperature: eint and E)" : " (eint)") << ", k in ["
                << skmin << "," << skmax << "], signs " << ssign << ", kdist "
                << skdist << ", r " << srlo << " .. " << srhi << std::endl;
      for (int n=0; n<snk; ++n) {
        std::cout << "  seed mode " << n << ": k2 " << smd.h_view(n,0) << " k3 "
                  << smd.h_view(n,1) << " a " << smd.h_view(n,2) << std::endl;
      }
    }
  }
  smd.modify_host();
  smd.sync_device();
  auto smd_d = smd.d_view;
  // problem/star_seed_ylm = true (read
  // only when named) replaces the lateral modes by a random combination of REAL spherical
  // harmonics Y_lm(theta, phi) in the GLOBAL polar angles of the cell centre (cubed
  // sphere and spherical polar alike), l in [star_seed_lmin, star_seed_lmax], every m =
  // -l..l, each coefficient Gaussian (own mt19937 stream seeded by star_seed_rng),
  // normalised so the lateral field has rms 1/sqrt 2 over the full sphere (the multi-mode
  // convention), under the same radial envelope.  Smooth at the poles (no taper) and with
  // no preferred azimuthal period: the 90-degree wedge periods of the modes above are
  // exactly 4-fold symmetric in phi, which is the cubed sphere's own panel symmetry.
  // (star_seed_band = n2 forces Y_lm with its own l range; cap 32 -> 64)
  const bool sylm = sn2 || (pin->DoesParameterExist("problem","star_seed_ylm") &&
                            pin->GetBoolean("problem","star_seed_ylm"));
  // (n2: an explicit star_seed_lmin / star_seed_lmax in the deck overrides the l window)
  const int slmin = sylm ? pin->GetOrAddInteger("problem","star_seed_lmin",
                                                sn2 ? n2lmin : 1) : 0;
  const int slmax = sylm ? pin->GetOrAddInteger("problem","star_seed_lmax",
                                                sn2 ? n2lmax : 8) : 0;
  if (sylm && (slmin < 1 || slmax < slmin || slmax > 64 || !(srhi > srlo))) {
    StarFatal("star_seed_ylm needs 1 <= star_seed_lmin <= star_seed_lmax <= 64 and "
              "star_seed_rhi > star_seed_rlo", __LINE__);
  }
  const int nylm = sylm ? (slmax + 1)*(slmax + 1) : 1;
  DualArray1D<Real> sylc("star_seed_ylm", nylm);   // index l*l + l + m: c_lm N_lm
  for (int q = 0; q < nylm; ++q) {sylc.h_view(q) = 0.0;}
  if (sylm) {
    std::mt19937 rngy(static_cast<unsigned>(srng) + 0x9e37u);
    std::normal_distribution<double> g01(0.0, 1.0);
    double ss = 0.0;
    for (int l = slmin; l <= slmax; ++l) {
      for (int mm = -l; mm <= l; ++mm) {
        const double a = g01(rngy);
        sylc.h_view(l*l + l + mm) = a;
        ss += a*a;
      }
    }
    // orthonormal real Y (int Y^2 dOmega = 1): <field^2> = sum c^2/(4 pi) = 1/2
    const double cn = std::sqrt(2.0*M_PI/ss);
    for (int l = slmin; l <= slmax; ++l) {
      for (int mm = -l; mm <= l; ++mm) {
        const int am = std::abs(mm);
        const double lnf = std::lgamma(l - am + 1.0) - std::lgamma(l + am + 1.0);
        double nlm = std::sqrt((2.0*l + 1.0)/(4.0*M_PI)*std::exp(lnf));
        if (mm != 0) {nlm *= std::sqrt(2.0);}
        sylc.h_view(l*l + l + mm) *= cn*nlm;
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "star_vet: seed = real spherical harmonics l " << slmin << ".."
          << slmax
                << " (" << (slmax + 1)*(slmax + 1) - slmin*slmin << " terms, all m), "
                << "amplitude " << seed
                << (srad ? " (temperature: eint and E)" : " (eint)")
                << ", r " << srlo << " .. " << srhi << ", rng " << srng << std::endl;
    }
  }
  sylc.modify_host();
  sylc.sync_device();
  auto sylc_d = sylc.d_view;
  // STAGE CS3: on the cubed sphere the modes are laid down in the GLOBAL polar angles of
  // the cell centre, with the 90 x 90 degree wedge's periods (theta_a = pi/4, L_theta =
  // L_phi = pi/2: the same k2, k3 give the same wavelengths as on the production wedge;
  // phi-periodic over 2 pi since 4 k3 is an integer), tapered by sin(theta) so the seed
  // is continuous at the poles (two panel centres)
  const Real x2a = hcs ? (0.25*M_PI) : pmy_mesh_->mesh_size.x2min;
  const Real x3a = hcs ? 0.0 : pmy_mesh_->mesh_size.x3min;
  const Real lth = hcs ? (0.5*M_PI) : (pmy_mesh_->mesh_size.x2max - x2a);
  const Real lph = hcs ? (0.5*M_PI) : (pmy_mesh_->mesh_size.x3max - x3a);
  auto mbpan = pmbp->pmb->mb_panel;
  auto &msz = pmbp->pmb->mb_size;
  const int hjs = indcs.js, hks = indcs.ks, hnx2 = indcs.nx2, hnx3 = indcs.nx3;
  auto uh = ph->u0;
  auto ur = pvet->u0;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto crho = star_rho_, ceint = star_eint_;
  const bool gmc = star_gmc_;
  auto pgc = star_phig_;
  const bool bal = star_bal_;
  auto cbd = star_bd_, cbe = star_be_;
  DvceArray1D<Real> cwv("star_wv", wic ? nf : 1);
  if (wic) {
    auto hw = Kokkos::create_mirror_view(cwv);
    for (int n=0; n<nf; ++n) hw(n) = hV[n];
    Kokkos::deep_copy(cwv, hw);
  }
  // (star_ic_equil): the seed factor is stored and applied after the equilibrium
  const bool sdef = star_eq_ && (seed != 0.0);
  star_eq_seed_ = sdef;
  star_eq_srad_ = srad;
  if (sdef) {
    Kokkos::realloc(star_sfac_, nmb1 + 1, n3m1 + 1, n2m1 + 1, n1m1 + 1);
    Kokkos::deep_copy(star_sfac_, 1.0);
  }
  auto sfc = star_sfac_;
  par_for("star_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const Real d = bal ? cbd(i) : StarLogInterp(crho, rlo, dr, nf, r);
    Real e = bal ? cbe(i) : StarLogInterp(ceint, rlo, dr, nf, r);
    // seed (default off): one smooth lateral mode in eint, wave number star_seed_k
    Real er = StarLogInterp(cE, rlo, dr, nf, r);
    if (seed != 0.0 && r >= srlo && r <= srhi) {
      Real fac;
      Real th2 = x2v(m,j), ph3 = x3v(m,k), tap = 1.0;
      if (hcs) {
        const Real xi = 0.25*M_PI*CellCenterX(j - hjs, hnx2, msz.d_view(m).x2min,
                                               msz.d_view(m).x2max);
        const Real eta = 0.25*M_PI*CellCenterX(k - hks, hnx3, msz.d_view(m).x3min,
                                                msz.d_view(m).x3max);
        Real q[3];
        cubed_sphere::PanelToCart(mbpan.d_view(m), xi, eta, q);
        th2 = acos(fmin(fmax(q[2], -1.0), 1.0));
        ph3 = atan2(q[1], q[0]);
        tap = sin(th2);
      }
      if (sylm) {
        // sum over m of P_l^m(cos theta) [c_lm cos(m phi) + c_l,-m sin(m phi)], the
        // associated Legendre functions by the stable upward recurrence in l (no
        // Condon-Shortley phase: the coefficients are random)
        const Real x = cos(th2), sth = sin(th2);
        Real amp = 0.0, pmm = 1.0;
        for (int mm = 0; mm <= slmax; ++mm) {
          if (mm > 0) {pmm *= (2.0*mm - 1.0)*sth;}
          const Real cm = cos(mm*ph3), sm = sin(mm*ph3);
          Real pl2 = 0.0, pl1 = pmm;   // P_{l-2}, P_{l-1}; start at l = mm
          for (int l = mm; l <= slmax; ++l) {
            Real pl;
            if (l == mm) {
              pl = pmm;
            } else if (l == mm + 1) {
              pl = x*(2.0*mm + 1.0)*pmm;
            } else {
              pl = (x*(2.0*l - 1.0)*pl1 - (l + mm - 1.0)*pl2)/(l - mm);
            }
            if (l > mm) {pl2 = pl1; pl1 = pl;}
            if (l >= slmin) {
              const Real cp = sylc_d(l*l + l + mm);
              amp += (mm == 0) ? cp*pl : pl*(cp*cm + sylc_d(l*l + l - mm)*sm);
            }
          }
        }
        fac = 1.0 + seed*sin(M_PI*(r - srlo)/(srhi - srlo))*amp;
      } else if (snk > 0) {
        Real amp = 0.0;
        for (int n=0; n<snk; ++n) {
          amp += smd_d(n,2)*sin(2.0*M_PI*(smd_d(n,0)*(th2 - x2a)/lth
                                          + smd_d(n,1)*(ph3 - x3a)/lph)
                                + smd_d(n,3));
        }
        if (hcs) {amp *= tap;}   // (the sp expression below is the pre-CS3 one)
        fac = 1.0 + seed*sin(M_PI*(r - srlo)/(srhi - srlo))*amp;
      } else {
        fac = 1.0 + seed*sin(2.0*M_PI*seedk*(th2 - x2a)/lth + 0.3)
                        *sin(2.0*M_PI*seedk*(ph3 - x3a)/lph + 1.1);
        if (hcs) {fac = 1.0 + (fac - 1.0)*tap;}
      }
      if (sdef) {
        sfc(m,k,j,i) = fac;
      } else {
        e *= fac;
        if (srad) er *= SQR(SQR(fac));
      }
    }
    const Real vr = wic ? StarLinInterp(cwv, rlo, dr, nf, r) : 0.0;
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = d*vr;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    if (gmc) {
      uh(m,IEN,k,j,i) = e + d*StarLinInterp(pgc, rlo, dr, nf, r);
    } else {
      uh(m,IEN,k,j,i) = e + d*gm*(1.0/rin - 1.0/r);
    }
    if (wic) uh(m,IEN,k,j,i) += 0.5*d*vr*vr;
    ur(m,radvet::VET_E,k,j,i) = fmax(er, efl);
    ur(m,radvet::VET_F1,k,j,i) = StarLinInterp(cF, rlo, dr, nf, r);
    ur(m,radvet::VET_F2,k,j,i) = 0.0;
    ur(m,radvet::VET_F3,k,j,i) = 0.0;
  });
  // the comoving x1 face fluxes of the implicit transport (allocated only there)
  auto ff0 = pvet->f0x1;
  if (ff0.extent_int(0) >= nmb1 + 1) {
    const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
    const int e1 = ff0.extent_int(3) - 1;
    par_for("star_icf", DevExeSpace(), 0, nmb1, 0, e3, 0, e2, 0, e1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      ff0(m,k,j,i) = StarLinInterp(cF, rlo, dr, nf, x1f(m,i));
    });
  }
  if (!hlin) return;

  // ---- problem/star_ic_profile = true.  The initial state of the BSG reference
  // Athena++ pgen bsg_star.cpp (ProblemGenerator), in cgs (star_ic_profile_file = its
  // MESA_profile_combined profile: r T rho M kappa F_r per line, cgs, ascending r; that
  // code reads the same table in code units, a pure rescaling):
  //  * r <= rmax (the file's last r): rho, T, F_r LINEAR in r between the bracketing rows
  //    (its lleft/lright search and formula, verbatim);
  //  * r > rmax: T = T(rmax), rho = max(rho_rmax exp(-g (r - rmax)/(P/rho)), dfloor),
  //    F_r = 0, with g = G M_rmax/rmax^2 (r/rmax)^2 EXACTLY as that code (the factor
  //    grows outward; the paper's "local gravity" would be (rmax/r)^2: reproduced, not
  //    fixed), P/rho = this EOS's code temperature T/T_unit (gas pressure only);
  //  * gas energy from the UNPERTURBED rho and T (its IEN = rho T/(gamma - 1)), E = a
  // T^4, F1 = F_r (its two-zone intensities have exactly these moments), F2 = F3 = 0, v =
  // 0; * the perturbation rho *= 1 + amp (U - 0.5) where the unperturbed rho > rho_thr,
  //    U uniform per cell from StarCellUniform(global cell index) (there: srand(gid),
  //    rand()): isobaric (P, E unchanged).  The potential energy of etotgrav uses the
  // perturbed rho. * the comoving x1 face fluxes f0x1 = the profile's F_r at the faces.
  // Fresh starts only (a restart returned above).  The ghost cells are filled the same
  // way and then overwritten by the (unchanged) boundary conditions.
  {
    if (hcs || wic) {
      StarFatal("problem/star_ic_profile needs the spherical-polar wedge and no "
                "star_wind_ic", __LINE__);
    }
    if (pmbp->punit != nullptr && (std::fabs(pmbp->punit->length_cgs() - 1.0) > 1.0e-12
        || std::fabs(pmbp->punit->density_cgs() - 1.0) > 1.0e-12)) {
      StarFatal("problem/star_ic_profile needs <units> length_cgs = mass_cgs = 1 (cgs r, "
                "rho)", __LINE__);
    }
    const Real lamp = pin->GetOrAddReal("problem",StarIcKey(pin,"star_ic_profile_amp"),
                                        0.05);
    // the reference threshold rho > 0.01 in ITS code units: rho_thr [cgs] = 0.01
    // rho_unit (of its input file, not supplied).  No default: the run states its
    // assumption.
    const Real lthr = pin->GetReal("problem",StarIcKey(pin,"star_ic_profile_rho_thr"));
    const int lrng = pin->GetOrAddInteger("problem",StarIcKey(pin,"star_ic_profile_rng"),
                                         0);
    const int np = static_cast<int>(pr.size());
    const int nm = nmb1 + 1, nc1 = n1m1 + 1;
    DualArray2D<Real> lD("star_lin_d", nm, nc1), lE("star_lin_eg", nm, nc1);
    DualArray2D<Real> lR("star_lin_er", nm, nc1), lF("star_lin_f", nm, nc1);
    DualArray2D<Real> lFf("star_lin_ff", nm, nc1 + 1);
    DualArray2D<std::int64_t> lof("star_lin_off", nm, 3);
    auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
    auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
    const Real tk = star_tk;
    for (int m=0; m<nm; ++m) {
      const LogicalLocation &ll = pmy_mesh_->lloc_eachmb[pmbp->gids + m];
      lof.h_view(m,0) = static_cast<std::int64_t>(ll.lx1)*indcs.nx1;
      lof.h_view(m,1) = static_cast<std::int64_t>(ll.lx2)*indcs.nx2;
      lof.h_view(m,2) = static_cast<std::int64_t>(ll.lx3)*indcs.nx3;
      for (int i=0; i<nc1; ++i) {
        const Real x1 = hx1v(m,i);
        Real rho, tem, flx;
        lprof(x1, rho, tem, flx);
        Real ee, pp, cr, ct, cv;
        eos.ThermoAt(rho, tem/tk, ee, pp, cr, ct, cv);
        const Real tc = tem/tk;
        lD.h_view(m,i) = rho;
        lE.h_view(m,i) = ee;
        lR.h_view(m,i) = ar*tc*tc*tc*tc;
        lF.h_view(m,i) = flx;
      }
      for (int i=0; i<=nc1; ++i) {
        const Real x1 = hx1f(m,i);
        Real rho, tem, flx;
        lprof(x1, rho, tem, flx);
        lFf.h_view(m,i) = flx;
      }
    }
    for (auto *a : {&lD, &lE, &lR, &lF, &lFf}) {a->modify_host(); a->sync_device();}
    lof.modify_host(); lof.sync_device();
    auto ld = lD.d_view, le = lE.d_view, lr = lR.d_view, lf = lF.d_view;
    auto lff = lFf.d_view;
    auto lo = lof.d_view;
    const int lis = indcs.is, ljs = indcs.js, lks = indcs.ks;
    const std::int64_t gn1 = pmy_mesh_->mesh_indcs.nx1 + 2*ng;
    const std::int64_t gn2 = pmy_mesh_->mesh_indcs.nx2 + 2*ng;
    const std::uint64_t useed = static_cast<std::uint64_t>(lrng);
    const bool m3 = (indcs.nx3 > 1), m2 = (indcs.nx2 > 1);
    // problem/star_ic_profile_seed_group = true (read only when named): one
    // deviate per hydro merge_lat group, so a merged 128^2 shell has the seed of the real
    // coarse cells (same amp and rms) instead of the average of 4 fine deviates.
    DualArray2D<int> lfg("star_lin_fg", nm, nc1);
    {
      std::vector<Real> mrc;
      std::vector<int> mfs;
      const std::string ksg = StarIcKey(pin,"star_ic_profile_seed_group");
      if (pin->DoesParameterExist("problem",ksg) && pin->GetBoolean("problem",ksg)) {
        if (!pin->DoesParameterExist("hydro","merge_lat_rc")) {
          StarFatal("problem/star_ic_profile_seed_group needs hydro/merge_lat_rc",
                    __LINE__);
        }
        std::istringstream sr(pin->GetString("hydro","merge_lat_rc"));
        std::istringstream sf(pin->GetString("hydro","merge_lat_f"));
        std::string tk1;
        while (std::getline(sr, tk1, ',')) mrc.push_back(std::stod(tk1));
        while (std::getline(sf, tk1, ',')) mfs.push_back(std::stoi(tk1));
        if (mrc.size() != mfs.size()) {
          StarFatal("merge_lat_rc and merge_lat_f sizes differ", __LINE__);
        }
      }
      for (int m=0; m<nm; ++m) {
        for (int i=0; i<nc1; ++i) {
          int fg = 1;
          for (std::size_t n=0; n<mrc.size(); ++n) {
            if (hx1v(m,i) < mrc[n]) {fg = mfs[n]; break;}
          }
          lfg.h_view(m,i) = fg;
        }
      }
      lfg.modify_host(); lfg.sync_device();
    }
    auto lfgd = lfg.d_view;
    par_for("star_ic_profile", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real r = x1v(m,i);
      const Real d0 = bal ? cbd(i) : ld(m,i), eg = bal ? cbe(i) : le(m,i);
      Real d = d0;
      if (d0 > lthr) {
        const std::int64_t gi = lo(m,0) + (i - lis) + ng;
        const int fg = lfgd(m,i);
        const std::int64_t gj = m2 ? (lo(m,1) + ((j - ljs)/fg)*fg + ng) : 0;
        const std::int64_t gk = m3 ? (lo(m,2) + ((k - lks)/fg)*fg + ng) : 0;
        const std::uint64_t idx = static_cast<std::uint64_t>((gk*gn2 + gj)*gn1 + gi);
        d = d0*(1.0 + lamp*(StarCellUniform(idx, useed) - 0.5));
      }
      uh(m,IDN,k,j,i) = d;
      uh(m,IM1,k,j,i) = 0.0;
      uh(m,IM2,k,j,i) = 0.0;
      uh(m,IM3,k,j,i) = 0.0;
      if (gmc) {
        uh(m,IEN,k,j,i) = eg + d*StarLinInterp(pgc, rlo, dr, nf, r);
      } else {
        uh(m,IEN,k,j,i) = eg + d*gm*(1.0/rin - 1.0/r);
      }
      ur(m,radvet::VET_E,k,j,i) = fmax(lr(m,i), efl);
      ur(m,radvet::VET_F1,k,j,i) = lf(m,i);
      ur(m,radvet::VET_F2,k,j,i) = 0.0;
      ur(m,radvet::VET_F3,k,j,i) = 0.0;
    });
    if (ff0.extent_int(0) >= nmb1 + 1) {
      const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
      const int e1 = std::min(ff0.extent_int(3) - 1, nc1);
      par_for("star_icf_profile", DevExeSpace(), 0, nmb1, 0, e3, 0, e2, 0, e1,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        ff0(m,k,j,i) = lff(m,i);
      });
    }
    if (global_variable::my_rank == 0) {
      int nthr = 0;
      for (int n=0; n<np; ++n) {
        if (pd[n] > lthr) nthr = n;
      }
      std::cout << "star_vet: star_ic_profile: " << lfn << " (" << np << " rows, r "
                << pr[0] << " .. rmax " << lrmax << "), T(rmax) " << ltmax
                << " K, rho(rmax) "
                << ldmax << ", g(rmax) " << lgmax
                << (gxi ? " (gext invsq: hydrostatic in GM/r^2 above, GM = star_gm)"
                        : " (x (r/rmax)^2 above, as bsg_star.cpp)")
                << "; perturbation amp " << lamp << " where rho > " << lthr
                << " (last profile row above: r " << pr[nthr] << "), rng " << lrng
                << "; star_seed off" << std::endl;
    }
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn void StarGravity()
//! \brief the point-mass source with the x1 well-balanced pressure form (etotgrav), the
//! reference work of the WB kick (force_reference_work = split), the wall/top mass-flux
//! correction and the optional top sponge.

void StarGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto *ph = &fl;
  auto u0 = ph->u0;
  auto w0 = ph->w0;
  auto eos = ph->peos->eos_data;
  auto wbq0 = ph->wbq0;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  // the adaptive closure: F_sub once per cycle, at its first stage (before any VET solve
  // of the cycle), esrc every stage (StarMltC)
  const Real tc = StarMltC(pm, bdt)*pm->dt;
  const Real wc = StarMltW(pm->time + tc);
  const bool ad = star_ad_ && (pm->ncycle != star_ad_cyc_);
  // the closure is skipped once the ramp has ended (w(t) = 0 AND w = 0 already applied
  // to esrc): F_sub then enters nothing (esrc = w F_sub deposit, history slot 21 =
  // w L_sub), so F_sub and its clock stay FROZEN at their last values (also in the
  // restart files), except when a problem/mlt_closure_dout profile is due (diagnostic
  // cadence kept: the update runs and writes its lines)
  if (ad) {
    star_ad_cyc_ = pm->ncycle;
    const bool off = (star_mrs_ >= 0.0 && wc == 0.0 && star_mw_ == 0.0);
    const bool dout = (star_ad_dout_ > 0.0 && pm->time >= star_ad_tout_);
    if (!off || dout) StarAdaptiveUpdate(pm);
  }
  if (star_mlt_ && (wc != star_mw_ || (star_ad_ && wc > 0.0 && (ad || star_mtc_)))) {
    StarApplyMltW(pm, wc, tc);
  }
  auto htsum = [&]() {
    Real s = 0.0;
    const int nx1 = ie - is + 1, nx2 = je - js + 1, nx3 = ke - ks + 1;
    Kokkos::parallel_reduce("star_htally", Kokkos::RangePolicy<>(DevExeSpace(), 0,
                            (nmb1 + 1)*nx3*nx2*nx1),
    KOKKOS_LAMBDA(const int n, Real &acc) {
      const int m = n/(nx3*nx2*nx1);
      const int k = (n/(nx2*nx1))%nx3 + ks, j = (n/nx1)%nx2 + js, i = n%nx1 + is;
      acc += u0(m,IEN,k,j,i)*volume(m,k,j,i);
    }, Kokkos::Sum<Real>(s));
#if MPI_PARALLEL_ENABLED
    Real g = 0.0;
    MPI_Allreduce(&s, &g, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    s = g;
#endif
    return s;
  };
  if (star_pm_) {
    // star_gravity = pointmass: the Athena++ point-mass source
    // (srcterms PointMass): d(rho v_r) = -dt rho GM coord_src1/x1v, coord_src1 =
    // (r_+^2 - r_-^2)/2 / ((r_+^3 - r_-^3)/3).  The energy needs no source: etotgrav
    // carries rho Phi in the total energy and its flux.
    auto &x1vg = pmbp->pcoord->x1v;
    auto &x1fg = pmbp->pcoord->xx1f;
    const Real gmp = star_gm_;
    par_for("star_grav_pm", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real rm = x1fg(m,i), rp = x1fg(m,i+1);
      const Real src1 = 1.5*(rp*rp - rm*rm)/(rp*rp*rp - rm*rm*rm);
      u0(m,IM1,k,j,i) -= bdt*w0(m,IDN,k,j,i)*gmp*src1/x1vg(m,i);
    });
  } else if (star_wbwin_) {
  // problem/star_wb_grav_window: outside the fluid's WB window (x1v > wb_rmax > 0 or
  // phicc_wb > wb_phimax: exactly the cells whose x1
  // reconstruction is plain, hydro_fluxes.cpp / mhd_fluxes.cpp) the cached hydrostatic
  // walk describes no supported atmosphere (H << dr above a star): its face pressures
  // fall by orders of magnitude across the cell and the pressure form over-pulls the gas.
  // There the face pressures are walked at the CELL density, pl = p + d (Phi_c - Phi_l),
  // pr = p - d (Phi_r - Phi_c): the source becomes -d [A_r (Phi_r - Phi_c) + A_l (Phi_c -
  // Phi_l)]/V, the plain -rho dPhi_wb/dr in the same face form and on the same effective
  // potential (so the a_ref split and the energy flux are untouched).  The switch is per
  // cell and coincides with the reconstruction's; the two forms differ by O(dr/H) of
  // the force where H is resolved, so the transition at wb_rmax is continuous to that
  // order.  Window absent (both at their defaults): bitwise the branch below.
  const Real wrmax = fl.mhd ? fl.pm->wb_rmax : fl.ph->wb_rmax;
  const Real wphmax = fl.mhd ? fl.pm->wb_phimax : fl.ph->wb_phimax;
  auto &x1vw = pmbp->pcoord->x1v;
  auto pwcw = ph->phicc_wb, pwfw = ph->phi_wb_x1f;
  par_for("star_grav_win", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    Real pl, pr, d1, d2, d3;
    WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
    const Real e_ = w0(m,IEN,k,j,i);
    const Real p = (e_ > 0.0) ? eos.Pressure(d, e_) : 0.5*(pl + pr);
    const Real r = x1vw(m,i);
    if ((wrmax > 0.0 && r > wrmax) || pwcw(m,k,j,i) > wphmax) {
      pl = p + d*(pwcw(m,k,j,i) - pwfw(m,k,j,i));
      pr = p - d*(pwfw(m,k,j,i+1) - pwcw(m,k,j,i));
    }
    u0(m,IM1,k,j,i) += bdt*(area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))
                       /volume(m,k,j,i);
  });
  } else {
  par_for("star_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    Real pl, pr, d1, d2, d3;
    WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
    const Real e_ = w0(m,IEN,k,j,i);
    const Real p = (e_ > 0.0) ? eos.Pressure(d, e_) : 0.5*(pl + pr);
    u0(m,IM1,k,j,i) += bdt*(area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))
                       /volume(m,k,j,i);
  });
  }
  // <rad_vet>/force_reference_work = split: the work of the rho a_ref part of the WB kick
  // at the stage-start velocity; the radiation pays it in its solve
  auto *pvet = pmbp->pradvet;
  if (pvet != nullptr && pvet->fref_wsplit) {
    auto aref = pvet->arad_ref;
    if (pvet->FrefWaccOn()) {
      // (time_scheme = be): the same increment also goes to the
      // accumulator the be solve makes E pay (acc = gam0 acc + increment)
      auto wacc = pvet->FrefWacc();
      const Real g0 = pvet->FrefWaccGam0(bdt, pm->dt);
      par_for("star_grav_fws_x", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real w = bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
        u0(m,IEN,k,j,i) += w;
        wacc(m,k,j,i) = g0*wacc(m,k,j,i) + w;
      });
    } else {
      par_for("star_grav_fws", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        u0(m,IEN,k,j,i) += bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
      });
    }
  }
  // the RK-exact boundary/sponge ledgers of this stage (star_drive.cpp)
  const bool hdrv = stardrive::On();
  if (star_sp_rate_ <= 0.0) {
    if (hdrv) {stardrive::StageLedger(pm, bdt, star_binf_, star_zflux_, 0.0);}
    return;
  }
  const Real hdsp0 = hdrv ? htsum() : 0.0;
  // top sponge: the velocity relaxes to zero at the rate star_sponge_rate*w,
  // w = ((r - r0)/(r_top - r0))^2 above r0; the kinetic energy removed leaves the total
  // energy (the heat is not kept, so that the sponge adds no buoyancy).
  // star_sponge_mode = radial: only v_r (the transverse momenta are untouched)
  auto &x1v = pmbp->pcoord->x1v;
  const Real rate = star_sp_rate_, r0 = star_sp_r0_, rt = star_rtop_, spd = star_sp_dmax_;
  // star_sponge_form = lawson: this stage's abscissa c (c_s = gam0_s c_{s-1} + beta_s)
  // from the integrator's weights (found by beta_s dt = bdt): u0 = fa (u0 - g1 u^n) +
  // fb g1 u^n, fa = exp(-ca rate w), fb = exp(-cb rate w)
  const bool law = star_sp_law_;
  auto u1 = ph->u1;
  Real ca = bdt, cb = bdt, g1 = 0.0;
  if (law) {
    Real c = 0.0, cp = 0.0;
    int s = -1;
    for (int q = 0; q < pvet->fref_hnst && s < 0; ++q) {
      cp = c;
      c = pvet->fref_hgam0[q]*c + pvet->fref_hbeta[q];
      if (std::abs(bdt - pvet->fref_hbeta[q]*pm->dt) <= 1.0e-12*bdt) s = q;
    }
    if (s < 0) StarFatal("problem/star_sponge_form = lawson: no rk1/rk2/rk3 stage has "
                         "this beta dt", __LINE__);
    ca = (c - cp)*pm->dt;
    cb = c*pm->dt;
    g1 = 1.0 - pvet->fref_hgam0[s];
  }
  if (star_sp_rad_) {
    par_for("star_sponge_r", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real r = x1v(m,i);
      const bool low = ((law ? u1(m,IDN,k,j,i) : u0(m,IDN,k,j,i)) < spd);
      if (r <= r0 && !low) return;
      const Real w = low ? 1.0 : SQR((r - r0)/(rt - r0));
      const Real f = exp(-ca*rate*w);
      const Real ke1 = 0.5*SQR(u0(m,IM1,k,j,i))/u0(m,IDN,k,j,i);
      if (law) {
        u0(m,IM1,k,j,i) = f*u0(m,IM1,k,j,i) + (exp(-cb*rate*w) - f)*g1*u1(m,IM1,k,j,i);
        u0(m,IEN,k,j,i) -= ke1 - 0.5*SQR(u0(m,IM1,k,j,i))/u0(m,IDN,k,j,i);
        return;
      }
      u0(m,IM1,k,j,i) *= f;
      u0(m,IEN,k,j,i) -= (1.0 - f*f)*ke1;
    });
    if (hdrv) {stardrive::StageLedger(pm, bdt, star_binf_, star_zflux_, htsum() - hdsp0);}
    return;
  }
  const bool hcsg = pm->use_cubed_sphere;   // metric KE (StarKeCs)
  auto hccl = pmbp->pcoord->cos_cell;
  par_for("star_sponge", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const bool low = ((law ? u1(m,IDN,k,j,i) : u0(m,IDN,k,j,i)) < spd);
    if (r <= r0 && !low) return;
    const Real w = low ? 1.0 : SQR((r - r0)/(rt - r0));
    const Real f = exp(-ca*rate*w);
    const Real d = u0(m,IDN,k,j,i);
    const Real ke0 = hcsg ? StarKeCs(u0(m,IM1,k,j,i), u0(m,IM2,k,j,i), u0(m,IM3,k,j,i), d,
                                     hccl(m,k,j))
                          : 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                          + SQR(u0(m,IM3,k,j,i)))/d;
    if (law) {
      const Real fu = (exp(-cb*rate*w) - f)*g1;
      const Real m1 = f*u0(m,IM1,k,j,i) + fu*u1(m,IM1,k,j,i);
      const Real m2 = f*u0(m,IM2,k,j,i) + fu*u1(m,IM2,k,j,i);
      const Real m3 = f*u0(m,IM3,k,j,i) + fu*u1(m,IM3,k,j,i);
      u0(m,IM1,k,j,i) = m1;
      u0(m,IM2,k,j,i) = m2;
      u0(m,IM3,k,j,i) = m3;
      u0(m,IEN,k,j,i) -= ke0 - (hcsg ? StarKeCs(m1, m2, m3, d, hccl(m,k,j))
                                     : 0.5*(SQR(m1) + SQR(m2) + SQR(m3))/d);
      return;
    }
    u0(m,IM1,k,j,i) *= f;
    u0(m,IM2,k,j,i) *= f;
    u0(m,IM3,k,j,i) *= f;
    u0(m,IEN,k,j,i) -= (1.0 - f*f)*ke0;
  });
  if (hdrv) {stardrive::StageLedger(pm, bdt, star_binf_, star_zflux_, htsum() - hdsp0);}
}

//----------------------------------------------------------------------------------------
//! \fn void StarWallFlux()
//! \brief star_wall_zero_flux (user_hydro_flux_func): the scaled-profile ghosts are not a
//! mirror image of the edge cell, so the Riemann flux through a "closed" x1 wall carries
//! mass.  The inner wall's fluxes are set to 0, the pressure (IM1) flux stays; the outer
//! face keeps OUTFLOW and is closed only where its mass flux is inflow or exactly 0 here.
//! Before the update: the source terms and am_transport see that flux (L_z exact), as do
//! the scalars.  dfloor_flux_limit runs later (Hydro::RKUpdate): a top face that it makes
//! massless (th = 0) it closes by the same rule (DflFace: all but IM1 set to 0).

void StarWallFlux(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie;
  auto flx1 = radvet::FluidRef::Get(pmbp).uflx->x1f;
  const int nv = flx1.extent_int(1);
  auto &mbbcs = pmbp->pmb->mb_bcs;
  const bool binf = star_binf_;
  auto mn = star_min_;
  par_for("star_wallflux", DevExeSpace(), 0, pmbp->nmb_thispack - 1, indcs.ks, indcs.ke,
          indcs.js, indcs.je,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    const bool z[2] = {
      !binf && mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user,
      mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user &&
      !(flx1(m,IDN,k,j,ie+1) > 0.0)};
    mn(m,k,j) = z[1] ? flx1(m,IDN,k,j,ie+1) : 0.0;   // (Min_top)
    for (int n = 0; n < nv; ++n) {
      if (n == IM1) continue;
      if (z[0]) flx1(m,n,k,j,is) = 0.0;
      if (z[1]) flx1(m,n,k,j,ie+1) = 0.0;
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn Real StarZcAm()
//! \brief star_bc_inner = zerocontact: the HLLC contact speed S_M at the inner wall face
//! for a trial v_r = vg of the first ghost (index is-1), exactly as the x1 sweep forms
//! it: the plain spherical-polar PLM (Hydro::PLM_nonuniform, Mignone 2014) of the left
//! state from cells (is-2, is-1, is) and of the right state from (is-1, is, is+1), then
//! hllc_hyd.hpp steps 2-5 (ideal gas).  dl/el, dr/er: the reconstructed density and
//! internal energy (independent of vg); v2g: v_r of the second ghost (0); va, vb: v_r of
//! the cells is, is+1; cl* / cr*: the PLM grid factors (dxL, dxR, dxLh, dxRh) of the
//! cells is-1 and is.

KOKKOS_INLINE_FUNCTION
Real StarZcAm(const EOS_Data &eos, const Real vg, const Real v2g, const Real va,
              const Real vb, const Real dl, const Real el, const Real dr, const Real er,
              const Real *cl4, const Real *cr4) {
  Real vl, vr, dum;
  hydro::Hydro::PLM_nonuniform(v2g, vg, va, cl4[0], cl4[1], cl4[2], cl4[3], vl, dum);
  hydro::Hydro::PLM_nonuniform(vg, va, vb, cr4[0], cr4[1], cr4[2], cr4[3], dum, vr);
  const Real g = eos.gamma, alpha = (g + 1.0)/(2.0*g);
  const Real pl = eos.IdealGasPressure(el), pr = eos.IdealGasPressure(er);
  const Real qa = eos.SoundSpeedFromP(dl, pl, g);
  const Real qb = eos.SoundSpeedFromP(dr, pr, g);
  const Real qc = 0.25*(dl + dr)*(qa + qb);
  const Real pm = 0.5*(pl + pr + (vl - vr)*qc);
  const Real ql = (pm <= pl) ? 1.0 : sqrt(1.0 + alpha*((pm/pl) - 1.0));
  const Real qr = (pm <= pr) ? 1.0 : sqrt(1.0 + alpha*((pm/pr) - 1.0));
  const Real al = vl - qa*ql;
  const Real ar = vr + qb*qr;
  const Real vxl = vl - al, vxr = vr - ar;
  const Real tl = pl + vxl*dl*vl, tr = pr + vxr*dr*vr;
  const Real ml = dl*vxl, mr = -(dr*vxr);
  return (tl - tr)/(ml + mr);
}

//----------------------------------------------------------------------------------------
//! \fn void StarBC()
//! \brief x1 ghosts carrying the scaled initial profile (rho and eint scaled by the ratio
//! of the adjacent active cell to its own initial value), the radial velocity mirrored
//! (wall) or, at the outer edge where it points outward, copied (outflow).  VET ghosts by
//! VetFillGhost (the implicit solve imposes its own face BCs).

void StarBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto &x1v = pmbp->pcoord->x1v;
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto *ph = &fl;
  auto uh = ph->u0;
  auto phicc = ph->phicc0;
  auto crho = star_rho_, ceint = star_eint_;
  const Real rlo = star_rlo_, dr = star_dr_;
  const int nf = star_nf_;
  const bool bal = star_bal_;
  auto cbd = star_bd_, cbe = star_be_;
  const bool binf = star_binf_, bout = star_bout_;
  const bool bhse = star_bhse_, bhtg = star_bhse_tg_;
  auto ur = pmbp->pradvet->u0;
  const Real cl = pmbp->pradvet->c_light;
  const Real efl = pmbp->pradvet->e_floor;
  // star_bc_outer = hse
  auto eos = ph->peos->eos_data;
  auto pwc = ph->phicc_wb;
  auto aref = pmbp->pradvet->arad_ref;
  auto opac = pmbp->pradvet->opac;
  const bool haveop = (opac.extent_int(0) > 0);
  const Real gm = star_gm_, gmax = star_bhse_gmax_;
  const bool gmc = star_gmc_;
  auto cgm = star_gmr_;
  const Real dfl = eos.dfloor;
  const bool bhfc = bhse && star_bhse_face_;
  auto ff1 = pmbp->pradvet->f0x1;
  // star_bc_inner = zerocontact, star_bc_outer = copy
  const bool bzc = star_bzc_, bcopy = star_bcopy_;
  const bool bcr2 = star_bcr2_;
  auto &x1fb = pmbp->pcoord->xx1f;
  const bool hcsb = pm->use_cubed_sphere;   // metric KE (StarKeCs)
  auto hcclb = pmbp->pcoord->cos_cell;
  // problem/star_bc_inner_table (star_bc_inner = inflow): the table's rho,
  // T (code units) and v_r at r_in at the BC time
  const bool itab = stardrive::InflowOn();
  Real tbd = 0.0, tbt = 0.0, tbv = 0.0;
  if (itab) {stardrive::InflowAt(pm, &tbd, &tbt, &tbv);}
  const Real tbd0 = star_drv_d0_, tbe0 = star_drv_e0_, tbr2 = SQR(star_rin_);
  par_for("star_bc", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    for (int side=0; side<2; ++side) {
      const bool lo = (side == 0);
      if (mbbcs.d_view(m, lo ? BoundaryFace::inner_x1 : BoundaryFace::outer_x1)
          != BoundaryFlag::user) continue;
      const int ia = lo ? is : ie;
      const Real da = uh(m,IDN,k,j,ia);
      const Real kea = hcsb ? StarKeCs(uh(m,IM1,k,j,ia), uh(m,IM2,k,j,ia),
                                       uh(m,IM3,k,j,ia), da, hcclb(m,k,j))
                            : 0.5*(SQR(uh(m,IM1,k,j,ia)) + SQR(uh(m,IM2,k,j,ia))
                            + SQR(uh(m,IM3,k,j,ia)))/da;
      const Real ea = uh(m,IEN,k,j,ia) - kea - da*phicc(m,k,j,ia);
      const bool inflo = lo && binf, outfl = !lo && bout;
      if (lo && bzc) {
        // ZERO-CONTACT INNER WALL (Ma, Bildsten & Jiang 2026, sect. 2.2.3).  All inner
        // ghosts hold the UNSCALED initial column (rho, eint at their radii) and v = 0,
        // except the first ghost (is-1): its v_r is the root vg of the HLLC contact speed
        // S_M(vg) = 0 at the wall face, with the left/right states the x1 sweep will
        // reconstruct (plain PLM, StarZcAm), v_t = 0 there.  S_M = 0 makes the HLLC mass,
        // transverse-momentum and energy fluxes vanish; the momentum flux is p*.
        // S_M increases with vg; bracketed Illinois iteration; no bracket -> mirror v.
        const int g0 = is - 1, g1 = is - 2, ib = is + 1;
        const Real d0 = StarLogInterp(crho, rlo, dr, nf, x1v(m,g0));
        const Real e0 = StarLogInterp(ceint, rlo, dr, nf, x1v(m,g0));
        const Real d1 = StarLogInterp(crho, rlo, dr, nf, x1v(m,g1));
        const Real e1 = StarLogInterp(ceint, rlo, dr, nf, x1v(m,g1));
        const Real db = uh(m,IDN,k,j,ib);
        const Real vb = uh(m,IM1,k,j,ib)/db;
        const Real eb = uh(m,IEN,k,j,ib) - 0.5*(SQR(uh(m,IM1,k,j,ib))
                        + SQR(uh(m,IM2,k,j,ib)) + SQR(uh(m,IM3,k,j,ib)))/db
                        - db*phicc(m,k,j,ib);
        const Real va = uh(m,IM1,k,j,is)/da;
        Real cl4[4], cr4[4];
        cl4[0] = x1v(m,g0) - x1v(m,g1);
        cl4[1] = x1v(m,is) - x1v(m,g0);
        cl4[2] = x1v(m,g0) - x1fb(m,g0);
        cl4[3] = x1fb(m,is) - x1v(m,g0);
        cr4[0] = x1v(m,is) - x1v(m,g0);
        cr4[1] = x1v(m,ib) - x1v(m,is);
        cr4[2] = x1v(m,is) - x1fb(m,is);
        cr4[3] = x1fb(m,ib) - x1v(m,is);
        Real dlf, elf, drf, erf, dum;
        hydro::Hydro::PLM_nonuniform(d1, d0, da, cl4[0], cl4[1], cl4[2], cl4[3],
                                     dlf, dum);
        hydro::Hydro::PLM_nonuniform(e1, e0, ea, cl4[0], cl4[1], cl4[2], cl4[3],
                                     elf, dum);
        hydro::Hydro::PLM_nonuniform(d0, da, db, cr4[0], cr4[1], cr4[2], cr4[3],
                                     dum, drf);
        hydro::Hydro::PLM_nonuniform(e0, ea, eb, cr4[0], cr4[1], cr4[2], cr4[3],
                                     dum, erf);
        const Real csa = sqrt(eos.gamma*eos.IdealGasPressure(ea)/da);
        Real w = 2.0*(csa + fabs(va) + fabs(vb));
        Real xa = -w, xb = w;
        Real fa = StarZcAm(eos, xa, 0.0, va, vb, dlf, elf, drf, erf, cl4, cr4);
        Real fb = StarZcAm(eos, xb, 0.0, va, vb, dlf, elf, drf, erf, cl4, cr4);
        for (int it=0; it<30 && !(fa < 0.0 && fb > 0.0); ++it) {
          w *= 2.0;
          if (!(fa < 0.0)) {
            xa = -w;
            fa = StarZcAm(eos, xa, 0.0, va, vb, dlf, elf, drf, erf, cl4, cr4);
          }
          if (!(fb > 0.0)) {
            xb = w;
            fb = StarZcAm(eos, xb, 0.0, va, vb, dlf, elf, drf, erf, cl4, cr4);
          }
        }
        Real vg = -va;
        if (fa < 0.0 && fb > 0.0) {
          // Illinois (modified regula falsi), stopped at |S_M| <= 1e-14 c_s
          const Real ftol = 1.0e-14*csa;
          int side = 0;
          Real xc = xa, fc = fa;
          for (int it=0; it<100; ++it) {
            xc = (xa*fb - xb*fa)/(fb - fa);
            fc = StarZcAm(eos, xc, 0.0, va, vb, dlf, elf, drf, erf, cl4, cr4);
            if (fabs(fc) <= ftol || !(xb - xa > 1.0e-15*w)) break;
            if (fc < 0.0) {
              xa = xc;
              fa = fc;
              if (side == -1) fb *= 0.5;
              side = -1;
            } else {
              xb = xc;
              fb = fc;
              if (side == 1) fa *= 0.5;
              side = 1;
            }
          }
          vg = xc;
        }
        for (int g=0; g<ng; ++g) {
          const int ig = is - 1 - g;
          const Real rg = x1v(m,ig);
          const Real dg = (g == 0) ? d0 : ((g == 1) ? d1 :
                                           StarLogInterp(crho, rlo, dr, nf, rg));
          const Real eg = (g == 0) ? e0 : ((g == 1) ? e1 :
                                           StarLogInterp(ceint, rlo, dr, nf, rg));
          const Real v1 = (g == 0) ? vg : 0.0;
          uh(m,IDN,k,j,ig) = dg;
          uh(m,IM1,k,j,ig) = dg*v1;
          uh(m,IM2,k,j,ig) = 0.0;
          uh(m,IM3,k,j,ig) = 0.0;
          uh(m,IEN,k,j,ig) = eg + 0.5*dg*v1*v1 + dg*phicc(m,k,j,ig);
          radvet::VetFillGhost(ur, m, k, j, ig, k, j, ia, 1, 0, -1.0, cl, efl);
        }
        continue;
      }
      if (!lo && bcopy) {
        // Ma+2026 OUTFLOW TOP: rho, eint and v_t of the edge cell copied into every
        // ghost, v_r copied where it points out and 0 where it points in
        const Real v1e = uh(m,IM1,k,j,ia)/da;
        Real v1 = (v1e > 0.0) ? v1e : 0.0;
        const Real v2 = uh(m,IM2,k,j,ia)/da;
        const Real v3 = uh(m,IM3,k,j,ia)/da;
        // star_bc_outer_r2: its branch test is v_r < 0 (v_r = 0 -> outflow)
        const bool r2o = bcr2 && !(v1e < 0.0);
        if (bcr2) v1 = r2o ? v1e : 0.0;
        Real dg = da;
        for (int g=0; g<ng; ++g) {
          const int ig = ie + 1 + g;
          if (r2o) dg *= SQR(x1v(m,ig-1)/x1v(m,ig));
          uh(m,IDN,k,j,ig) = dg;
          uh(m,IM1,k,j,ig) = dg*v1;
          uh(m,IM2,k,j,ig) = dg*v2;
          uh(m,IM3,k,j,ig) = dg*v3;
          uh(m,IEN,k,j,ig) = ea + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,j,ig);
          radvet::VetFillGhost(ur, m, k, j, ig, k, j, ia, 1, 2, 1.0, cl, efl);
        }
        continue;
      }
      if (!lo && bhse) {
        // THE HYDROSTATIC TOP (star_bc_outer = hse).  Every ghost layer continues the
        // LAST ACTIVE cell a = ie (box_convection bc_mode_top 4, ported): with the WB
        // pair's effective potential Phi_eff (gravity minus the reference radiation force
        // a_ref of the column, what the x1 well-balanced stencil walks along) plus the
        // RESIDUAL radiation force the VET module applies on top of it,
        //   a_rad = opac_T F_1/(c rho) of cell a  (the code's own force; the reference
        //           a_ref when the opacity is not yet filled, i.e. on the first step),
        //           F_1 per star_bc_hse_flux: face (default) = 0.5 (f0x1(ie) +
        //           f0x1(ie+1)), the comoving flux the transport carries and the gas
        //           force uses; cell = the clipped lab-frame cell F1 (old; noisy, ~0.76
        //           of the face flux in the 3-D BSG top where |F| = cE clips)
        //   Gamma = a_rad/g(r_a) clamped to [0, star_bc_hse_gmax], scaled ~ r^-2,
        // the ghost is the cell walked by dPhi = Phi_eff(r_g) - Phi_eff(r_a)
        //   - (Gamma g_a - a_ref(a)) r_a^2 (1/r_a - 1/r_g)
        // with the WB closure (utils/wb_background.hpp WBAdvance): ISOTHERMAL (T_g = T_a,
        // exact for the ideal gas: rho_g = rho_a exp(-dPhi rho_a/p_a)) or, with
        // star_bc_hse_tgrad, polytropic along the interior's one-sided dT/dPhi_eff.  So
        // the total effective gravity over the ghosts is g (1 - Gamma), and with
        // a_rad = a_ref the ghost IS the WB background of the edge cell: the top-face
        // Riemann problem then carries no hydrostatic residual (no one-signed mass flux).
        // Not finite / not positive -> copy of the edge cell; floors after.
        const int ia1 = ie - 1;
        const Real ta = eos.Temperature(da, ea);
        const Real ra = x1v(m,ia);
        const Real ga = gmc ? StarLinInterp(cgm, rlo, dr, nf, ra)/SQR(ra) : gm/SQR(ra);
        const Real ar0 = aref(m,k,j,ia);
        Real arad = ar0;
        if (haveop && opac(m,radvet::VET_OP_T,k,j,ia) > 0.0) {
          // face (star_bc_hse_flux): the comoving radial flux the transport carries, the
          // mean of the cell's two x1 faces = the unclipped cell F1 minus its lab-frame
          // (v + v.D) E; cell: the clipped lab-frame cell F1 (old)
          const Real f1 = bhfc ? 0.5*(ff1(m,k,j,ia) + ff1(m,k,j,ia+1))
                               : ur(m,radvet::VET_F1,k,j,ia);
          arad = opac(m,radvet::VET_OP_T,k,j,ia)*f1/(cl*da);
        }
        Real gam = arad/ga;
        gam = (gam < 0.0) ? 0.0 : ((gam > gmax) ? gmax : gam);
        const Real dres = gam*ga - ar0;
        Real dltdphi = 0.0;
        if (bhtg) {
          const Real d1 = uh(m,IDN,k,j,ia1);
          const Real e1 = uh(m,IEN,k,j,ia1) - (hcsb ? StarKeCs(uh(m,IM1,k,j,ia1),
                          uh(m,IM2,k,j,ia1), uh(m,IM3,k,j,ia1), d1, hcclb(m,k,j))
                          : 0.5*(SQR(uh(m,IM1,k,j,ia1))
                          + SQR(uh(m,IM2,k,j,ia1)) + SQR(uh(m,IM3,k,j,ia1)))/d1)
                          - d1*phicc(m,k,j,ia1);
          const Real dph1 = pwc(m,k,j,ia) - pwc(m,k,j,ia1);
          if (dph1 != 0.0) dltdphi = (ta - eos.Temperature(d1, e1, ta))/dph1;
        }
        const Real v1e = uh(m,IM1,k,j,ia)/da;
        const Real v1 = (v1e > 0.0) ? v1e : 0.0;
        const Real v2 = uh(m,IM2,k,j,ia)/da;
        const Real v3 = uh(m,IM3,k,j,ia)/da;
        for (int g=0; g<ng; ++g) {
          const int ig = ie + 1 + g;
          const Real rg = x1v(m,ig);
          const Real dphi = (pwc(m,k,j,ig) - pwc(m,k,j,ia))
                            - dres*SQR(ra)*(1.0/ra - 1.0/rg);
          Real dg = da, eg = ea, tg = ta;
          WBAdvance(eos, bhtg ? 3 : 1, da, ea, dphi, dg, eg, tg, ta, dltdphi, ta, ta);
          if (!(Kokkos::isfinite(dg) && (dg > 0.0) && Kokkos::isfinite(eg)
                && (eg > 0.0))) {
            dg = da;
            eg = ea;
          }
          if (dg < dfl) {
            eg *= dfl/dg;
            dg = dfl;
          }
          uh(m,IDN,k,j,ig) = dg;
          uh(m,IM1,k,j,ig) = dg*v1;
          uh(m,IM2,k,j,ig) = dg*v2;
          uh(m,IM3,k,j,ig) = dg*v3;
          uh(m,IEN,k,j,ig) = eg + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,j,ig);
          radvet::VetFillGhost(ur, m, k, j, ig, k, j, ia, 1, 2, 1.0, cl, efl);
        }
        continue;
      }
      if (inflo && itab) {
        // the initial column's shape scaled to the table's (rho,
        // e(rho, T)) at r_in; v_r carries the table's mass flux rho v r^2 through every
        // ghost; v_t = 0 (non-rotating interior)
        Real tbe, tbp, tcr, tct, tcv;
        eos.ThermoAt(tbd, tbt, tbe, tbp, tcr, tct, tcv);
        const Real sdt = tbd/tbd0, set = tbe/tbe0;
        for (int g=0; g<ng; ++g) {
          const int ig = is - 1 - g;
          const Real rg = x1v(m,ig);
          const Real dg = sdt*StarLogInterp(crho, rlo, dr, nf, rg);
          const Real eg = set*StarLogInterp(ceint, rlo, dr, nf, rg);
          const Real v1 = tbv*tbr2*tbd/(SQR(rg)*dg);
          uh(m,IDN,k,j,ig) = dg;
          uh(m,IM1,k,j,ig) = dg*v1;
          uh(m,IM2,k,j,ig) = 0.0;
          uh(m,IM3,k,j,ig) = 0.0;
          uh(m,IEN,k,j,ig) = eg + 0.5*dg*v1*v1 + dg*phicc(m,k,j,ig);
          radvet::VetFillGhost(ur, m, k, j, ig, k, j, ia, 1, 0, -1.0, cl, efl);
        }
        continue;
      }
      const Real sd = inflo ? 1.0 :
                      da/(bal ? cbd(ia) : StarLogInterp(crho, rlo, dr, nf, x1v(m,ia)));
      const Real se = inflo ? 1.0 :
                      ea/(bal ? cbe(ia) : StarLogInterp(ceint, rlo, dr, nf, x1v(m,ia)));
      for (int g=0; g<ng; ++g) {
        const int ig = lo ? (is - 1 - g) : (ie + 1 + g);
        const int im = lo ? (is + g) : (ie - g);        // the mirror cell
        const Real rg = x1v(m,ig);
        const Real dg = outfl ? da*SQR(x1v(m,ia)/rg) :
                        sd*(bal ? cbd(ig) : StarLogInterp(crho, rlo, dr, nf, rg));
        const Real eg = outfl ? ea*(dg/da) :
                        se*(bal ? cbe(ig) : StarLogInterp(ceint, rlo, dr, nf, rg));
        const Real dm = uh(m,IDN,k,j,im);
        // wall: v1 mirrored; outer edge: v1 of the edge cell where it flows out
        const Real v1e = uh(m,IM1,k,j,ia)/da;
        const Real v1 = inflo ? v1e : ((!lo && v1e > 0.0) ? v1e : -uh(m,IM1,k,j,im)/dm);
        const Real v2 = uh(m,IM2,k,j,im)/dm;
        const Real v3 = uh(m,IM3,k,j,im)/dm;
        uh(m,IDN,k,j,ig) = dg;
        uh(m,IM1,k,j,ig) = dg*v1;
        uh(m,IM2,k,j,ig) = dg*v2;
        uh(m,IM3,k,j,ig) = dg*v3;
        uh(m,IEN,k,j,ig) = eg + 0.5*dg*(v1*v1 + v2*v2 + v3*v3) + dg*phicc(m,k,j,ig);
        radvet::VetFillGhost(ur, m, k, j, ig, k, j, ia, 1, lo ? 0 : 2,
                            lo ? -1.0 : 1.0, cl, efl);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void StarHist()
//! \brief history columns, all plain sums (see the file header).  Etot = gas total energy
//! (conserved IEN incl. kinetic and the static potential) + radiation energy: with closed
//! walls its rate is L_in - L_top(lab) (the comoving L_top differs from it by O(v/c)).
//! Mdot_top is the mass really leaving through the top (the inflow part is removed by
//! StarWallFlux), Min_top the inflow that was removed.
//! M_hflx, E_hflx, Lz_hflx: the hydro's mass, energy (IEN) and L_z leaving through the
//! physical x1 faces over the step, as rates (am_bflx, hydro_am.hpp: the RK-weighted
//! fluxes); Lz_tot: the L_z of gas + radiation (LzTotal, incl. the L_z ledger), so
//! that d M_tot/dt = -M_hflx and d Lz_tot/dt = -Lz_hflx close in the history.

void StarHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  // w_mlt, L_MLT (= w max_faces F_MLT r^2 Omega) only with mlt_flux_frozen
  // F4: + P_esrc (power the explicit source puts into E, (chat/c) sum
  // esrc V) and P_sponge (kinetic-energy loss rate of the top sponge at this state), the
  // two non-flux terms of d Etot/dt = L_in - L_top + P_esrc - P_sponge
  const int nb = star_mlt_ ? 22 : 20;
  // L_midE, L_intE: the clean energy flux Phi of the step's last solve (t2pf1, what
  // moved E; c/chat its stored chat/c) through the L_mid / L_int faces; L_mid / L_int
  // are the faces' F (clean staggered: the face momentum, not the energy flux).
  // L_midC, L_intC: the comoving part Phi - (4/3) v E of the same flux (v, E the face
  // means; thick D = 1/3), what old's L_mid / L_int sum (old f0x1 = F - a E). L_bot is
  // that comoving flux too (old's L_bot, face is+1; the face F adds the wall's 4/3 v E)
  pdata->nhist = nb + 10;
  const char *lab[32] = {"L_bot", "L_mid", "L_int", "L_top", "L_in", "E_rad", "e_gas",
                         "M_int", "KE_int", "KEr_int", "Mr_int", "PV_int", "V_int",
                         "M_tot", "Mdot_top", "Mdot_bot", "Etot", "Min_top",
                         "v1sq_wall", "Picard", "w_mlt", "L_MLT", "P_esrc", "P_sponge",
                         "M_hflx", "E_hflx", "Lz_hflx", "Lz_tot", "L_midE", "L_intE",
                         "L_midC", "L_intC"};
  for (int n=0; n<nb; ++n) pdata->label[n] = lab[n];
  for (int n=0; n<10; ++n) pdata->label[nb+n] = lab[22+n];
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2, ks = indcs.ks, nx3 = indcs.nx3;
  const int imid = is + nx1/2;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto *ph = &fl;
  auto w0 = ph->w0;
  auto u0 = ph->u0;
  auto fl1 = ph->uflx->x1f;
  const bool hmn = star_zflux_;
  auto mn = star_min_;
  auto eos = ph->peos->eos_data;
  auto *pvet = pmbp->pradvet;
  auto ur = pvet->u0;
  auto ff0 = pvet->f0x1;
  const bool haveff = (ff0.extent_int(0) > 0);
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  const Real rint = star_rint_, fin = star_fin_;   // hg262-drive
  // clean staggered: L_in / L_top = the energy the step's solves put through the x1
  // faces over dt (BoundaryEnergyAcc), the flux that moved E
  const bool lac = pvet->cln_ldt > 0.0 && pvet->cln_lacc.extent_int(0) > 0;
  const bool hpf = pvet->t2pf1_ok && pvet->t2pf1.extent_int(0) > 0;
  auto pf1 = pvet->t2pf1;
  const Real icc = pvet->c_light/pvet->chat;
  auto lacc = pvet->cln_lacc;
  const Real ildt = lac ? 1.0/pvet->cln_ldt : 0.0;
  // v1sq_wall: mean v_r^2 over the wall (first active) cells of the whole mesh
  const Real inwall = 1.0/(static_cast<Real>(pm->mesh_indcs.nx2)*pm->mesh_indcs.nx3);
  const bool eso = pvet->esrc_on;
  auto es_ = pvet->esrc;
  const Real chc = pvet->chat/pvet->c_light;
  const Real sprate = star_sp_rate_, spr0 = star_sp_r0_;
  const Real sprt = star_rtop_, spd = star_sp_dmax_;
  const bool sprad = star_sp_rad_;
  const bool hcsh = pm->use_cubed_sphere;   // metric KE
  auto hccl_ = pmbp->pcoord->cos_cell;
  array_sum::GlobalSum sum_this;
  repro::HistSums("star_hist", pmbp->nmb_thispack, nkji, nx1,
  KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &msum) {
    const int m = idx/nkji;
    int k = (idx - m*nkji)/nji;
    int j = (idx - m*nkji - k*nji)/nx1;
    const int i = (idx - m*nkji - k*nji - j*nx1) + is;
    k += ks;
    j += js;
    array_sum::GlobalSum h;
    for (int n=0; n<NREDUCTION_VARIABLES; ++n) h.the_array[n] = 0.0;
    const Real vol = volume(m,k,j,i);
    if (haveff) {
      if (i == is && !hpf) h.the_array[0] = ff0(m,k,j,i+1)*area1(m,k,j,i+1);
      if (i == imid) h.the_array[1] = ff0(m,k,j,i)*area1(m,k,j,i);
      if (x1f(m,i) <= rint && x1f(m,i+1) > rint) {
        h.the_array[2] = ff0(m,k,j,i)*area1(m,k,j,i);
      }
      if (i == ie) h.the_array[3] = ff0(m,k,j,i+1)*area1(m,k,j,i+1);
    }
    const Real cof = (w0(m,IVX,k,j,i-1) + w0(m,IVX,k,j,i))*(ur(m,radvet::VET_E,k,j,i-1) +
                     ur(m,radvet::VET_E,k,j,i))*area1(m,k,j,i)/3.0;
    if (hpf && i == imid) h.the_array[24] = icc*pf1(m,0,k,j,i)*area1(m,k,j,i);
    if (hpf && i == imid) h.the_array[26] = h.the_array[24] - cof;
    if (hpf && i == is + 1) h.the_array[0] = icc*pf1(m,0,k,j,i)*area1(m,k,j,i) - cof;
    if (hpf && x1f(m,i) <= rint && x1f(m,i+1) > rint) {
      h.the_array[25] = icc*pf1(m,0,k,j,i)*area1(m,k,j,i);
      h.the_array[27] = h.the_array[25] - cof;
    }
    if (i == ie && lac) h.the_array[3] = lacc(m,1,k,j)*ildt;
    if (i == is) h.the_array[4] = lac ? lacc(m,0,k,j)*ildt : fin*area1(m,k,j,i);
    h.the_array[5] = ur(m,radvet::VET_E,k,j,i)*vol;
    h.the_array[6] = w0(m,IEN,k,j,i)*vol;
    h.the_array[13] = w0(m,IDN,k,j,i)*vol;
    if (i == ie) {
      const Real fm = fl1(m,IDN,k,j,i+1)*area1(m,k,j,i+1);
      h.the_array[14] = (fm > 0.0) ? fm : 0.0;
      h.the_array[17] = hmn ? mn(m,k,j)*area1(m,k,j,i+1) : ((fm > 0.0) ? 0.0 : fm);
    }
    if (i == is) {
      h.the_array[15] = fl1(m,IDN,k,j,i)*area1(m,k,j,i);
      h.the_array[18] = SQR(w0(m,IVX,k,j,i))*inwall;
    }
    h.the_array[16] = (u0(m,IEN,k,j,i) + ur(m,radvet::VET_E,k,j,i))*vol;
    if (eso) h.the_array[22] = chc*es_(m,k,j,i)*vol;
    if (sprate > 0.0) {
      const Real r = x1v(m,i);
      const Real d = u0(m,IDN,k,j,i);
      const bool low = (d < spd);
      if (r > spr0 || low) {
        const Real wsp = low ? 1.0 : SQR((r - spr0)/(sprt - spr0));
        const Real kes = sprad ? 0.5*SQR(u0(m,IM1,k,j,i))/d
                               : (hcsh ? StarKeCs(u0(m,IM1,k,j,i), u0(m,IM2,k,j,i),
                                                  u0(m,IM3,k,j,i), d, hccl_(m,k,j))
                               : 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                                      + SQR(u0(m,IM3,k,j,i)))/d);
        h.the_array[23] = 2.0*sprate*wsp*kes*vol;
      }
    }
    if (x1v(m,i) <= rint) {
      const Real d = w0(m,IDN,k,j,i);
      const Real v1 = w0(m,IVX,k,j,i), v2 = w0(m,IVY,k,j,i), v3 = w0(m,IVZ,k,j,i);
      h.the_array[7] = d*vol;
      // w0 holds the CONTRAVARIANT v^a on the cubed sphere
      h.the_array[8] = hcsh ? 0.5*d*(v1*v1 + v2*v2 + v3*v3
                                     + 2.0*hccl_(m,k,j)*v2*v3)*vol
                            : 0.5*d*(v1*v1 + v2*v2 + v3*v3)*vol;
      h.the_array[9] = 0.5*d*v1*v1*vol;
      h.the_array[10] = d*v1*vol;
      h.the_array[11] = eos.Pressure(d, w0(m,IEN,k,j,i))*vol;
      h.the_array[12] = vol;
    }
    msum += h;
  }, sum_this);
  for (int n=0; n<20; ++n) pdata->hdata[n] = sum_this.the_array[n];
  pdata->hdata[nb] = sum_this.the_array[22];
  pdata->hdata[nb+1] = sum_this.the_array[23];
  // M_hflx, E_hflx, Lz_hflx: the step's am_bflx over its dt (0 without the module)
  pdata->hdata[nb+6] = sum_this.the_array[24];
  pdata->hdata[nb+7] = sum_this.the_array[25];
  pdata->hdata[nb+8] = sum_this.the_array[26];
  pdata->hdata[nb+9] = sum_this.the_array[27];
  for (int n=0; n<3; ++n) pdata->hdata[nb+2+n] = 0.0;
  hydro::AmModule *pam = hydro::AmModuleOf(pmbp);
  if (pam != nullptr && pam->am_bflx.extent_int(0) > 0 && pam->am_bdt > 0.0) {
    auto bf = pam->am_bflx;
    const int nkj = nx3*nx2;
    array_sum::GlobalSum fsum;
    repro::HistSums("star_hflx", pmbp->nmb_thispack, nkj, 1,
    KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &msum) {
      const int m = idx/nkj, k = (idx - m*nkj)/nx2 + ks, j = idx%nx2 + js;
      array_sum::GlobalSum h;
      for (int n=0; n<NREDUCTION_VARIABLES; ++n) h.the_array[n] = 0.0;
      for (int n=0; n<3; ++n) h.the_array[n] = bf(m,n,0,k,j) + bf(m,n,1,k,j);
      msum += h;
    }, fsum);
    for (int n=0; n<3; ++n) pdata->hdata[nb+2+n] = fsum.the_array[n]/pam->am_bdt;
  }
  // (a global value: rank 0 contributes)
  const Real lzt = pvet->LzTotal();
  pdata->hdata[nb+5] = (global_variable::my_rank == 0) ? lzt : 0.0;
  if (star_mlt_) {
    // host values, summed over ranks by the history output: rank 0 contributes
    const bool r0 = (global_variable::my_rank == 0);
    pdata->hdata[20] = r0 ? star_mw_ : 0.0;
    pdata->hdata[21] = r0 ? star_mw_*(star_ad_ ? star_lsub_ : star_lmx_) : 0.0;
  }
  // Picard: mean passes per implicit solve since the previous history output (the
  // counters are MPI_MAX-reduced, identical on every rank: rank 0 contributes)
  pdata->hdata[19] = 0.0;
  if (global_variable::my_rank == 0) {
    const Real dn = pvet->impl_nstep - star_pic_n0_, ds = pvet->impl_itsum - star_pic_s0_;
    pdata->hdata[19] = (dn > 0.0) ? ds/dn : 0.0;
    star_pic_n0_ = pvet->impl_nstep;
    star_pic_s0_ = pvet->impl_itsum;
  }
  p_ledger::Hist(pdata, pm);   // <hydro>/p_ledger: the momentum ledger's columns
}

//----------------------------------------------------------------------------------------
//! \fn void StarApplyMltW()
//! \brief rebuild <rad_vet>/esrc = w (frozen MLT deposit) + the star_esrc_const test
//! source, with exactly the operations of the start-up kernels (w = 1: bitwise the same
//! esrc).  h: the time past F_sub's (mlt_src_time = stage: c dt).

void StarApplyMltW(Mesh *pm, const Real w, const Real h) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto es = pmbp->pradvet->esrc;
  auto fmd = star_ad_ ? star_fsub_ : star_fm_;
  auto ftg = (star_ad_ && star_mtc_) ? star_tgp_ : fmd;
  const Real a = 1.0 - exp(-h/star_ad_tau_);
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  auto &x1v = pmbp->pcoord->x1v;
  const Real esc = star_esc_, rmx = star_esrmx_;
  const bool bh = (star_bheat_ > 0.0);
  auto fbd = star_fb_;
  par_for("star_esrc_w", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real fl = fmd(i) + a*(ftg(i) - fmd(i));
    const Real fr = fmd(i+1) + a*(ftg(i+1) - fmd(i+1));
    Real e = -(area1(m,k,j,i+1)*fr - area1(m,k,j,i)*fl)/volume(m,k,j,i);
    if (bh) {
      const Real eb = -(area1(m,k,j,i+1)*fbd(i+1) - area1(m,k,j,i)*fbd(i))
                      /volume(m,k,j,i);
      e = w*e + (1.0 - w)*eb;
    } else {
      e = w*e;
    }
    if (esc != 0.0 && x1v(m,i) < rmx) e += esc;
    es(m,k,j,i) = e;
  });
  star_mw_ = w;
}

//----------------------------------------------------------------------------------------
//! \fn void StarAdaptiveUpdate()
//! \brief the adaptive shell-mean deficit closure (star_ad_): shell sums (deterministic:
//! one thread per radial index, serial over the columns; MPI sum), the target, the
//! relaxation F_sub += (1 - exp(-(t - t_last)/tau)) (target - F_sub), clip >= 0.

// the rank-ordered sum of a host array over all ranks (MPI_Allgather, then
// the ranks summed 0, 1, ... in that order on every rank).  An MPI_Allreduce(MPI_SUM)
// does not fix the order in which the ranks' doubles are added: with HCOLL the order
// changes with the communication history, and a restart is not bitwise.
// Same pattern as IcrDetSum / vet_gd_twin_det.
void StarDetSum(Real *v, const int n) {
#if MPI_PARALLEL_ENABLED
  const int nr = global_variable::nranks;
  if (nr <= 1 || n <= 0) return;
  std::vector<Real> mine(v, v + n);
  std::vector<Real> all(static_cast<size_t>(n)*nr);
  MPI_Allgather(mine.data(), n, MPI_ATHENA_REAL, all.data(), n, MPI_ATHENA_REAL,
                MPI_COMM_WORLD);
  for (int i = 0; i < n; ++i) {
    Real a = 0.0;
    for (int r = 0; r < nr; ++r) {a += all[static_cast<size_t>(r)*n + i];}
    v[i] = a;
  }
#else
  (void) v; (void) n;
#endif
}

void StarAdaptiveUpdate(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto w0 = fl.w0;
  auto eos = fl.peos->eos_data;
  auto *pvet = pmbp->pradvet;
  auto ur = pvet->u0;
  auto f0 = pvet->f0x1;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  auto sh = star_shs_;
  // star_ad_dead_skip: the serial sums only when they are used
  const bool dead = star_ad_dskip_ && pm->layout_invariant;
  if (!dead) {
  par_for("star_ad_shell", DevExeSpace(), is, ie+1, KOKKOS_LAMBDA(const int i) {
    Real s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0, s4 = 0.0, s5 = 0.0, s6 = 0.0, s7 = 0.0;
    for (int m=0; m<nmb; ++m) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          s0 += f0(m,k,j,i)*area1(m,k,j,i);
          s1 += area1(m,k,j,i);
          if (i <= ie) {
            const Real vol = volume(m,k,j,i);
            const Real d = w0(m,IDN,k,j,i), e = w0(m,IEN,k,j,i);
            const Real h = e + eos.Pressure(d, e);
            const Real v = w0(m,IVX,k,j,i), er = ur(m,radvet::VET_E,k,j,i);
            s2 += vol;
            s3 += vol*v;
            s4 += vol*h;
            s5 += vol*er;
            s6 += vol*v*h;
            s7 += vol*v*er;
          }
        }
      }
    }
    sh(i,0) = s0; sh(i,1) = s1; sh(i,2) = s2; sh(i,3) = s3;
    sh(i,4) = s4; sh(i,5) = s5; sh(i,6) = s6; sh(i,7) = s7;
  });
  }
  auto hsh = Kokkos::create_mirror_view(HostMemSpace(), sh);
  if (!dead) {Kokkos::deep_copy(hsh, sh);}
  const int nfc = hsh.extent_int(0);
  if (pm->layout_invariant) {
    // the same per-cell terms, summed order-independently over all ranks
    std::vector<Real> rs;
    repro::ShellSums("star_ad_shell_r", 8, nfc, nmb, ks, ke, js, je, is, ie + 1,
                     KOKKOS_LAMBDA(const int q, const int m, const int k, const int j,
                                   const int i) {
      if (q == 0) {return f0(m,k,j,i)*area1(m,k,j,i);}
      if (q == 1) {return area1(m,k,j,i);}
      if (i > ie) {return 0.0;}
      const Real vol = volume(m,k,j,i);
      if (q == 2) {return vol;}
      const Real d = w0(m,IDN,k,j,i), e = w0(m,IEN,k,j,i);
      const Real h = e + eos.Pressure(d, e);
      const Real v = w0(m,IVX,k,j,i), er = ur(m,radvet::VET_E,k,j,i);
      if (q == 3) {return vol*v;}
      if (q == 4) {return vol*h;}
      if (q == 5) {return vol*er;}
      if (q == 6) {return vol*v*h;}
      return vol*v*er;
    }, rs);
    for (int i=is; i<=ie+1; ++i) {
      for (int q=0; q<8; ++q) {hsh(i,q) = rs[static_cast<size_t>(q)*nfc + i];}
    }
  } else {
    StarDetSum(&hsh(is,0), 8*(ie + 2 - is));   // rank-ordered (restart-exact)
  }
  // luminosities of the wedge (the units of the history's L_in = F_in A(r_in))
  std::vector<Real> lc(nfc, 0.0), lrad(nfc, 0.0), lcf(nfc, 0.0), tgt(nfc, 0.0);
  for (int i=is; i<=ie; ++i) {
    const Real sv = hsh(i,2);
    const Real vm = hsh(i,3)/sv;
    const Real fc = (hsh(i,6)/sv - vm*hsh(i,4)/sv)
                    + (4.0/3.0)*(hsh(i,7)/sv - vm*hsh(i,5)/sv);
    lc[i] = fc*0.5*(hsh(i,1) + hsh(i+1,1));
  }
  const Real lin = star_fin_*hsh(is,1);
  for (int i=is; i<=ie+1; ++i) {
    lrad[i] = hsh(i,0);
    lcf[i] = (i > is && i <= ie) ? 0.5*(lc[i-1] + lc[i]) : 0.0;
    const Real def = lin - lrad[i] - lcf[i];
    tgt[i] = star_tap_[i]*((def > 0.0) ? def : 0.0)/hsh(i,1);
  }
  if (star_adp_) {
    // mlt_closure = pressure: target fm_IC(ln P_face) F_tot,IC(face); P_face =
    // the log-mean of the adjacent shell-mean P_gas + E/3 (the edge cell's at is, ie+1)
    auto sp = star_shp_;
    auto ww = w0;
    auto eo = eos;
    auto uu = ur;
    auto vv = volume;
    if (!dead) {
    par_for("star_adp_shell", DevExeSpace(), is, ie, KOKKOS_LAMBDA(const int i) {
      Real s0 = 0.0, s1 = 0.0, s2 = 0.0;
      for (int m=0; m<nmb; ++m) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real vol = vv(m,k,j,i);
            s0 += vol;
            s1 += vol*eo.Pressure(ww(m,IDN,k,j,i), ww(m,IEN,k,j,i));
            s2 += vol*uu(m,radvet::VET_E,k,j,i);
          }
        }
      }
      sp(i,0) = s0; sp(i,1) = s1; sp(i,2) = s2;
    });
    }
    auto hsp = Kokkos::create_mirror_view(HostMemSpace(), sp);
    if (!dead) {Kokkos::deep_copy(hsp, sp);}
    if (pm->layout_invariant) {   // order-independent
      const int nsp = hsp.extent_int(0);
      std::vector<Real> rs;
      repro::ShellSums("star_adp_shell_r", 3, nsp, nmb, ks, ke, js, je, is, ie,
                       KOKKOS_LAMBDA(const int q, const int m, const int k, const int j,
                                     const int i) {
        const Real vol = vv(m,k,j,i);
        if (q == 0) {return vol;}
        if (q == 1) {return vol*eo.Pressure(ww(m,IDN,k,j,i), ww(m,IEN,k,j,i));}
        return vol*uu(m,radvet::VET_E,k,j,i);
      }, rs);
      for (int i=is; i<=ie; ++i) {
        for (int q=0; q<3; ++q) {hsp(i,q) = rs[static_cast<size_t>(q)*nsp + i];}
      }
    } else {
      StarDetSum(&hsp(is,0), 3*(ie + 1 - is));   // rank-ordered
    }
    std::vector<Real> lp(nfc, 0.0);
    for (int i=is; i<=ie; ++i) lp[i] = log((hsp(i,1) + hsp(i,2)/3.0)/hsp(i,0));
    const int np = static_cast<int>(star_plnp_.size());
    for (int i=is; i<=ie+1; ++i) {
      const Real x = (i == is) ? lp[is] : ((i == ie + 1) ? lp[ie]
                                           : 0.5*(lp[i-1] + lp[i]));
      Real fm;
      if (x <= star_plnp_[0]) {
        fm = star_pfm_[0];
      } else if (x >= star_plnp_[np-1]) {
        fm = star_pfm_[np-1];
      } else {
        const int k = static_cast<int>(std::upper_bound(star_plnp_.begin(),
                                                        star_plnp_.end(), x)
                                       - star_plnp_.begin()) - 1;
        const Real a = (x - star_plnp_[k])/(star_plnp_[k+1] - star_plnp_[k]);
        fm = (1.0 - a)*star_pfm_[k] + a*star_pfm_[k+1];
      }
      tgt[i] = fm*star_ftot_[i];
    }
  }
  const Real t = pm->time;
  const Real a = 1.0 - exp(-(t - star_ad_tl_)/star_ad_tau_);
  star_ad_tl_ = t;
  auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_fsub_);
  // mlt_src_time = stage: the trapezoid of tgt^{n-1}, tgt^n; hp keeps tgt^n
  auto hp = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_tgp_);
  if (star_mtc_) {
    for (int i=0; i<nfc; ++i) {
      const Real tn = tgt[i];
      if (star_tgp_ok_) tgt[i] = 0.5*(hp(i) + tn);
      hp(i) = tn;
    }
    Kokkos::deep_copy(star_tgp_, hp);
    star_tgp_ok_ = 1;
  }
  star_lsub_ = 0.0;
  for (int i=0; i<nfc; ++i) {
    Real f = hf(i) + a*(tgt[i] - hf(i));
    hf(i) = (f > 0.0) ? f : 0.0;
    star_lsub_ = std::max(star_lsub_, hf(i)*star_r2o_[i]);
  }
  Kokkos::deep_copy(star_fsub_, hf);
  // profiles (problem/mlt_closure_dout > 0): rank 0 appends, per output, three lines
  // "t <kind> values over the faces is..ie+1" in units of L_in: fsub (w A F_sub),
  // lrad, lconv; and one "t fmlt" line (A F_MLT) at the first output
  if (star_ad_dout_ > 0.0 && t >= star_ad_tout_ && global_variable::my_rank == 0) {
    auto hm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_fm_);
    FILE *fp = fopen(star_ad_file_, "a");
    if (fp != nullptr) {
      const char *kind[4] = {"fsub", "lrad", "lconv", "fmlt"};
      for (int q=0; q<4; ++q) {
        if (q == 3 && star_ad_tout_ > 0.0) break;
        fprintf(fp, "%.10e %s", t, kind[q]);
        for (int i=is; i<=ie+1; ++i) {
          Real v = 0.0;
          if (q == 0) v = StarMltW(t)*hf(i)*hsh(i,1);
          if (q == 1) v = lrad[i];
          if (q == 2) v = lcf[i];
          if (q == 3) v = hm(i)*hsh(i,1);
          fprintf(fp, " %.6e", v/lin);
        }
        fprintf(fp, "\n");
      }
      // the shell-mean radiation energy density per CELL is..ie
      // (code units, not normalised; a trailing 0 keeps the face count)
      fprintf(fp, "%.10e emean", t);
      for (int i=is; i<=ie; ++i) {fprintf(fp, " %.6e", hsh(i,5)/hsh(i,2));}
      fprintf(fp, " 0\n");
      fclose(fp);
    }
    star_ad_tout_ = star_ad_dout_*floor(t/star_ad_dout_ + 1.0);
  } else if (star_ad_dout_ > 0.0 && t >= star_ad_tout_) {
    star_ad_tout_ = star_ad_dout_*floor(t/star_ad_dout_ + 1.0);
  }
}

//----------------------------------------------------------------------------------------
//! \fn std::vector<char> StarAdaptiveRstWrite()
//! \brief the adaptive closure's restart state: int32 nface, int32 0, Real t_last (the
//! update clock), nface Reals F_sub.

std::vector<char> StarAdaptiveRstWrite() {
  std::vector<char> out;
  // the open Picard window (solves, passes since the last history
  // output) as a tagged trailer behind the adaptive state (if any)
  auto pic = [&out]() {
    Real w[2] = {0.0, 0.0};
    const radvet::RadiationVet *pvet =
        (star_mesh_ != nullptr) ? star_mesh_->pmb_pack->pradvet : nullptr;
    if (pvet != nullptr) {
      w[0] = pvet->impl_nstep - star_pic_n0_;
      w[1] = pvet->impl_itsum - star_pic_s0_;
    }
    const std::size_t o = out.size();
    out.resize(o + sizeof(kStarPicTag) + sizeof(w));
    std::memcpy(out.data() + o, &(kStarPicTag[0]), sizeof(kStarPicTag));
    std::memcpy(out.data() + o + sizeof(kStarPicTag), &(w[0]), sizeof(w));
  };
  // the star_ic_equil a_ref as the last trailer (only when it ran / was restored)
  auto dic = [&out]() {
    if (star_eq_a_.empty()) return;
    const std::size_t na = star_eq_a_.size();
    const Real rn = static_cast<Real>(na);
    const std::size_t o = out.size();
    out.resize(o + (na + 1)*sizeof(Real) + sizeof(kStarDicTag));
    std::memcpy(out.data() + o, star_eq_a_.data(), na*sizeof(Real));
    std::memcpy(out.data() + o + na*sizeof(Real), &rn, sizeof(Real));
    std::memcpy(out.data() + o + (na + 1)*sizeof(Real), &(kStarDicTag[0]),
                sizeof(kStarDicTag));
  };
  // the auto ramp's two numbers as the very last trailer
  auto rmp = [&out]() {
    if (!star_mr_store_) return;
    const Real w[2] = {star_mrs_, star_mrt_};
    const std::size_t o = out.size();
    out.resize(o + sizeof(w) + sizeof(kStarRmpTag));
    std::memcpy(out.data() + o, &(w[0]), sizeof(w));
    std::memcpy(out.data() + o + sizeof(w), &(kStarRmpTag[0]), sizeof(kStarRmpTag));
  };
  if (!star_ad_ || star_fsub_.extent_int(0) <= 0) {
    pic(); dic(); rmp();
    stardrive::RstAppend(&out);   // the very last trailer (if on)
    return out;
  }
  const int nfc = star_fsub_.extent_int(0);
  auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_fsub_);
  // hdr[1] = 1: tgt^{n-1} follows F_sub (mlt_src_time = stage, once it is known)
  const int ntg = (star_mtc_ && star_tgp_ok_) ? nfc : 0;
  const std::int32_t hdr[2] = {static_cast<std::int32_t>(nfc), (ntg > 0) ? 1 : 0};
  std::vector<Real> tmp(nfc + 1 + ntg);
  tmp[0] = star_ad_tl_;
  for (int i=0; i<nfc; ++i) tmp[i+1] = hf(i);
  if (ntg > 0) {
    auto hp = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_tgp_);
    for (int i=0; i<nfc; ++i) tmp[nfc+1+i] = hp(i);
  }
  out.resize(sizeof(hdr) + tmp.size()*sizeof(Real));
  std::memcpy(out.data(), &(hdr[0]), sizeof(hdr));
  std::memcpy(out.data() + sizeof(hdr), tmp.data(), tmp.size()*sizeof(Real));
  pic();
  dic();
  rmp();
  stardrive::RstAppend(&out);   // the very last trailer (if on)
  return out;
}

//----------------------------------------------------------------------------------------
//! \fn void StarFinal()
//! \brief release the namespace-scope Views before Kokkos::finalize.

//----------------------------------------------------------------------------------------
//! \fn void StarAtmPhiEff
//! \brief the x1 WB pair's Phi_eff rebuilt above the face x1f(ia) from the
//! per-cell reference accelerations a[i - ia] (cells ia..ie), each segment taking the
//! true potential's own difference minus a dr; the ghosts above ie are shifted with the
//! face x1f(ie+1).  Radial vectors (block 0): pc0/pf0 the pgen's Phi_eff at the cells /
//! x1 faces, p0c/p0f the true potential, x1v/x1f the grid; out pc (n1), pf (n1 + 1),
//! equal to pc0/pf0 at and below the face x1f(ia).

void StarAtmPhiEff(const int ia, const int ie, const std::vector<Real> &a,
                   const std::vector<Real> &pc0, const std::vector<Real> &pf0,
                   const std::vector<Real> &p0c, const std::vector<Real> &p0f,
                   const std::vector<Real> &x1v, const std::vector<Real> &x1f,
                   std::vector<Real> &pc, std::vector<Real> &pf) {
  const int n1 = static_cast<int>(pc0.size());
  pc = pc0;
  pf = pf0;
  for (int i=ia; i<=ie; ++i) {
    const Real ai = a[i - ia];
    pf[i+1] = pf[i] + (p0f[i+1] - p0f[i]) - ai*(x1f[i+1] - x1f[i]);
    pc[i] = pf[i] + (p0c[i] - p0f[i]) - ai*(x1v[i] - x1f[i]);
  }
  const Real dsh = pf[ie+1] - pf0[ie+1];
  for (int i=ie+1; i<n1; ++i) {
    pc[i] = pc0[i] + dsh;
    pf[i+1] = pf0[i+1] + dsh;
  }
}

//! \fn void StarAtmRadial
//! \brief the radial (block 0, k = ks, j = js) profiles StarAtmPhiEff reads (device)
void StarAtmRadial(MeshBlockPack *pmbp, std::vector<Real> &pc0, std::vector<Real> &pf0,
                   std::vector<Real> &p0c, std::vector<Real> &p0f, std::vector<Real> &x1v,
                   std::vector<Real> &x1f, std::vector<Real> &ar0) {
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ks = indcs.ks, js = indcs.js;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phicc_wb);
  auto hpf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phi_wb_x1f);
  auto hp0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phicc0);
  auto hq0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phi0->x1f);
  auto hxv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto hxf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
  auto har = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pradvet->arad_ref);
  pc0.resize(n1); pf0.resize(n1 + 1); p0c.resize(n1); p0f.resize(n1 + 1);
  x1v.resize(n1); x1f.resize(n1 + 1); ar0.resize(n1);
  for (int i=0; i<n1; ++i) {
    pc0[i] = hpc(0,ks,js,i);
    p0c[i] = hp0(0,ks,js,i);
    x1v[i] = hxv(0,i);
    ar0[i] = har(0,ks,js,i);
  }
  for (int i=0; i<=n1; ++i) {
    pf0[i] = hpf(0,ks,js,i);
    p0f[i] = hq0(0,ks,js,i);
    x1f[i] = hxf(0,i);
  }
}

//! \fn void StarAtmWrite
//! \brief Phi_eff (cells >= ia, x1 faces >= ia + 1, ghosts included) and a_ref (cells
//! ia..ie) of every block from the radial vectors
void StarAtmWrite(MeshBlockPack *pmbp, const int ia, const std::vector<Real> &a,
                  const std::vector<Real> &pc, const std::vector<Real> &pf) {
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ie = indcs.ie;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phicc_wb);
  auto hpf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phi_wb_x1f);
  auto har = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pradvet->arad_ref);
  for (int m=0; m<hpc.extent_int(0); ++m) {
    for (int k=0; k<hpc.extent_int(1); ++k) {
      for (int j=0; j<hpc.extent_int(2); ++j) {
        for (int i=ia; i<n1; ++i) {
          hpc(m,k,j,i) = pc[i];
          hpf(m,k,j,i+1) = pf[i+1];
        }
        for (int i=ia; i<=ie; ++i) har(m,k,j,i) = a[i - ia];
      }
    }
  }
  Kokkos::deep_copy(fl.phicc_wb, hpc);
  Kokkos::deep_copy(fl.phi_wb_x1f, hpf);
  Kokkos::deep_copy(pmbp->pradvet->arad_ref, har);
}

//----------------------------------------------------------------------------------------
//! \fn void StarIcPost
//! \brief problem/star_ic_rad_atm_post: the thin atmosphere's discrete
//! hydrostatic re-balance at the gas T left by <rad_vet>/ic_rad, with a_ref = the
//! module's own face-mean force at that state (see the enrolment in UserProblem).
//! Called once, on a fresh start, by Driver::Initialize after IcRadRun.

void StarIcPost(Mesh *pm) {
  if (!star_atp_) return;
  MeshBlockPack *pmbp = pm->pmb_pack;
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto *pvet = pmbp->pradvet;
  auto eos = fl.eos;
  const WBOption wbo = fl.wb_option;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int ia = star_atp_ia_;
  const Real dfl = star_atp_dfl_, cl = pvet->c_light;
  if (ia < is + 2 || ia > ie) {
    StarFatal("star_ic_rad_atm_post: the first atmosphere cell is not inside the mesh",
              __LINE__);
  }
  std::vector<Real> pc0, pf0, p0c, p0f, x1v, x1f, ar0;
  StarAtmRadial(pmbp, pc0, pf0, p0c, p0f, x1v, x1f, ar0);
  auto hu = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.u0);
  auto hp0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phicc0);
  auto hf0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pvet->f0x1);
  // the cell T (host, the code's EOS) of every cell i >= ia - 2, its shell sums
  const int i0 = ia - 2;
  const int nr = ie - i0 + 1;
  const int nk = ke - ks + 1, nj = je - js + 1;
  std::vector<Real> tcell(static_cast<std::size_t>(nmb)*nk*nj*nr);
  // per radial index: sum T, sum rho, count, sum F0 (x1 face i); min T, max T
  std::vector<Real> sh(4*nr, 0.0), tmn(nr, 1.0e300), tmx(nr, -1.0e300);
  auto cid = [&](const int m, const int k, const int j, const int i) {
    return ((static_cast<std::size_t>(m)*nk + (k - ks))*nj + (j - js))*nr + (i - i0);
  };
  Real ff = 0.0;   // F0 of the top face x1f(ie+1), summed
  // <mesh>/layout_invariant: the shell sums order-independent (ReproSum)
  const bool lay = pm->layout_invariant;
  repro::ReproVec rv(lay ? 4*nr + 1 : 0);
  for (int m=0; m<nmb; ++m) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=i0; i<=ie; ++i) {
          const Real d = hu(m,IDN,k,j,i);
          const Real ek = 0.5*(SQR(hu(m,IM1,k,j,i)) + SQR(hu(m,IM2,k,j,i)) +
                               SQR(hu(m,IM3,k,j,i)))/d;
          const Real e = hu(m,IEN,k,j,i) - ek - d*hp0(m,k,j,i);
          const Real t = eos.Temperature(d, e);
          tcell[cid(m,k,j,i)] = t;
          if (lay) {
            rv.Add(4*(i - i0), t);
            rv.Add(4*(i - i0) + 1, d);
            rv.Add(4*(i - i0) + 2, 1.0);
            rv.Add(4*(i - i0) + 3, hf0(m,k,j,i));
          } else {
            sh[4*(i - i0)] += t;
            sh[4*(i - i0) + 1] += d;
            sh[4*(i - i0) + 2] += 1.0;
            sh[4*(i - i0) + 3] += hf0(m,k,j,i);
          }
          tmn[i - i0] = std::min(tmn[i - i0], t);
          tmx[i - i0] = std::max(tmx[i - i0], t);
        }
        if (lay) {
          rv.Add(4*nr, hf0(m,k,j,ie+1));
        } else {
          ff += hf0(m,k,j,ie+1);
        }
      }
    }
  }
  if (lay) {
    rv.Merge();
    for (int q=0; q<4*nr; ++q) {sh[q] = rv.Value(q);}
    ff = rv.Value(4*nr);
  }
#if MPI_PARALLEL_ENABLED
  if (!lay) {
    MPI_Allreduce(MPI_IN_PLACE, sh.data(), 4*nr, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &ff, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  }
  MPI_Allreduce(MPI_IN_PLACE, tmn.data(), nr, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, tmx.data(), nr, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  // the 1-D column on x1 indices (shell means): T, rho, the face flux F0
  std::vector<Real> bd(n1, 0.0), be(n1, 0.0), tc(n1, 0.0), d0(n1, 0.0), f0(n1 + 1, 0.0);
  for (int i=i0; i<=ie; ++i) {
    const Real cnt = sh[4*(i - i0) + 2];
    tc[i] = sh[4*(i - i0)]/cnt;
    bd[i] = sh[4*(i - i0) + 1]/cnt;
    f0[i] = sh[4*(i - i0) + 3]/cnt;
    d0[i] = bd[i];
  }
  f0[ie+1] = ff/sh[4*(ie - i0) + 2];
  tc[ie+1] = tc[ie];   // the outflow ghost: the edge cell's T and rho
  bd[ie+1] = bd[ie];
  d0[ie+1] = bd[ie+1];
  for (int i=i0; i<=ie+1; ++i) {
    Real p, cr, ct, cv;
    eos.ThermoAt(bd[i], tc[i], be[i], p, cr, ct, cv);
  }
  std::vector<Real> pc, pf;
  auto walk = [&](const int i, const Real di, Real &pimh, Real &piph) {
    Real ei, p, cr, ct, cv;
    eos.ThermoAt(di, tc[i], ei, p, cr, ct, cv);
    WBState s0, s1, s2, s3, s4;
    WBBackgroundStencil(eos, wbo, bd[i-1], di, bd[i+1], be[i-1], ei, be[i+1],
                        pc[i-1], pf[i], pc[i], pf[i+1], pc[i+1], s0, s1, s2, s3, s4,
                        tc[i-1], tc[i], tc[i+1]);
    pimh = s1.p;
    piph = s3.p;
  };
  // a_ref of the cells ia..ie: start from the pgen's
  std::vector<Real> a(ie - ia + 1), an(ie - ia + 1);
  for (int i=ia; i<=ie; ++i) a[i - ia] = ar0[i];
  int ilast = ia - 1, npass = 0, iprev = -1;
  Real rmx = 0.0, chg = 0.0;
  for (int pass=0; pass<60; ++pass) {
    ++npass;
    iprev = ilast;
    StarAtmPhiEff(ia, ie, a, pc0, pf0, p0c, p0f, x1v, x1f, pc, pf);
    // outward: PL_i(d_i) = PR_{i-1}(d_{i-1}), from ia; stop where rho reaches the floor
    for (int i=ia; i<=ie+1; ++i) {
      bd[i] = d0[i];
      Real p, cr, ct, cv;
      eos.ThermoAt(bd[i], tc[i], be[i], p, cr, ct, cv);
    }
    ilast = ia - 1;
    rmx = 0.0;
    for (int i=ia; i<=ie; ++i) {
      Real pl, pr;
      walk(i-1, bd[i-1], pl, pr);
      const Real target = pr;
      auto res = [&](const Real lx) {
        Real q1, q2;
        walk(i, exp(lx), q1, q2);
        return q1/target - 1.0;
      };
      Real x0 = log(bd[i]), x1 = x0 + 1.0e-6;
      Real fa = res(x0), fb = res(x1);
      for (int it=0; it<60 && fabs(fb) > 1.0e-15 && fb != fa; ++it) {
        const Real x2 = x1 - fb*(x1 - x0)/(fb - fa);
        x0 = x1; fa = fb; x1 = x2;
        fb = res(x1);
      }
      if (!std::isfinite(x1) || !(exp(x1) > dfl*(1.0 + 1.0e-6)) || fabs(fb) > 1.0e-10) {
        break;
      }
      rmx = std::max(rmx, fabs(fb));
      bd[i] = exp(x1);
      Real p, cr, ct, cv;
      eos.ThermoAt(bd[i], tc[i], be[i], p, cr, ct, cv);
      ilast = i;
    }
    if (ilast < ia) break;
    // the module's face-mean force at this state (cells ia..ilast; the floor above keeps
    // the pgen's a_ref)
    std::vector<Real> kd, kt, kk;
    for (int i=ia-1; i<=ie+1; ++i) {
      kd.push_back(bd[i]);
      kt.push_back(tc[i]);
    }
    StarKtEval(pvet->otab, kd, kt, kk);
    chg = 0.0;
    for (int i=ia; i<=ie; ++i) {
      if (i > ilast) {
        an[i - ia] = ar0[i];
        continue;
      }
      const int q = i - (ia - 1);
      const Real dl = x1f[i] - x1f[i-1], dc = x1f[i+1] - x1f[i];
      const Real du = x1f[i+2] - x1f[i+1];
      const Real rkc = bd[i]*kk[q];
      const Real rkl = (dc*bd[i-1]*kk[q-1] + dl*rkc)/(dl + dc);
      const Real rkr = (i == ie) ? rkc : (du*rkc + dc*bd[i+1]*kk[q+1])/(dc + du);
      an[i - ia] = 0.5*(rkl*f0[i] + rkr*f0[i+1])/(cl*bd[i]);
      chg = std::max(chg, fabs(an[i - ia]/a[i - ia] - 1.0));
    }
    if (chg < 1.0e-12 && ilast == iprev) break;
    a = an;
  }
  // the cells above ilast keep the pgen's a_ref (what a restart rebuilds; a no-op after
  // a converged pass) and Phi_eff follows a
  for (int i=ilast+1; i<=ie; ++i) a[i - ia] = ar0[i];
  StarAtmPhiEff(ia, ie, a, pc0, pf0, p0c, p0f, x1v, x1f, pc, pf);
  // write back: rho of the column, e from (rho, the cell's own T), v kept; a_ref and
  // Phi_eff of the pass that produced this rho
  Real dmx = 0.0, rdm = 0.0, tsp = 0.0, amx = 0.0, ram = 0.0;
  for (int i=ia; i<=ilast; ++i) {
    if (fabs(bd[i]/d0[i] - 1.0) > dmx) {
      dmx = fabs(bd[i]/d0[i] - 1.0);
      rdm = x1v[i];
    }
    if (fabs(a[i - ia] - ar0[i])/(star_gm_/SQR(x1v[i])) > amx) {
      amx = fabs(a[i - ia] - ar0[i])/(star_gm_/SQR(x1v[i]));
      ram = x1v[i];
    }
    tsp = std::max(tsp, tmx[i - i0]/tmn[i - i0] - 1.0);
  }
  for (int m=0; m<nmb; ++m) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=ia; i<=ilast; ++i) {
          const Real dold = hu(m,IDN,k,j,i), d = bd[i];
          const Real v1 = hu(m,IM1,k,j,i)/dold, v2 = hu(m,IM2,k,j,i)/dold;
          const Real v3 = hu(m,IM3,k,j,i)/dold;
          Real e, p, cr, ct, cv;
          eos.ThermoAt(d, tcell[cid(m,k,j,i)], e, p, cr, ct, cv);
          hu(m,IDN,k,j,i) = d;
          hu(m,IM1,k,j,i) = d*v1;
          hu(m,IM2,k,j,i) = d*v2;
          hu(m,IM3,k,j,i) = d*v3;
          hu(m,IEN,k,j,i) = e + 0.5*d*(v1*v1 + v2*v2 + v3*v3) + d*hp0(m,k,j,i);
        }
      }
    }
  }
  Kokkos::deep_copy(fl.u0, hu);
  StarAtmWrite(pmbp, ia, a, pc, pf);
  // recorded in the input (and so in every restart file): a restart rebuilds the pgen's
  // a_ref and Phi_eff and then applies these (UserProblem), bit for bit.  One short key
  // per re-balanced cell (a long value would pad every line of the restart file's
  // parameter dump past its 40 KB header limit); the cells above keep the pgen's a_ref
  {
    const int na = std::max(ilast - ia + 1, 0);
    star_pin_->SetString("problem", "star_ic_rad_atm_ia", std::to_string(ia));
    star_pin_->SetString("problem", "star_ic_rad_atm_na", std::to_string(na));
    char key[40], buf[40];
    for (int n=0; n<na; ++n) {
      snprintf(key, sizeof(key), "star_ic_rad_atm_a%03d", n);
      snprintf(buf, sizeof(buf), "%.17g", a[n]);
      star_pin_->SetString("problem", key, buf);
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << "star_vet: star_ic_rad_atm_post: cells " << ia << " .. " << ilast
              << " (r " << x1v[ia] << " .. " << x1v[std::max(ilast, ia)] << ") "
              << "re-balanced at the ic_rad gas T, " << npass << " passes (last max "
              << "|da/a| " << chg << ", face residual " << rmx << "): max |rho/rho_ic - "
              << "1| " << dmx << " at r = " << rdm << ", max |a_ref - a_ref,pgen|/g "
              << amx << " at r = " << ram << ", lateral T spread " << tsp << std::endl;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void StarIcEquil
//! \brief problem/star_ic_equil: one round of the discrete equilibrium of the
//! initial state (see the enrolment in UserProblem).  Called by Driver::Initialize after
//! every <rad_vet>/ic_rad run (fresh start); sets user_icpost_again while it wants
//! another round.  The x1 structure is taken from the shell means (the initial state is l
//! = 0); each cell keeps its own T and its rho/<rho> ratio.

void StarIcEquil(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  radvet::FluidRef fl = radvet::FluidRef::Get(pmbp);
  auto *pvet = pmbp->pradvet;
  auto eos = fl.eos;
  const WBOption wbo = fl.wb_option;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const Real dfl = star_eq_dfl_;
  // the march's anchor: the WALL cell (its rho kept), marched OUTWARD.  At fixed force
  // per volume f (fixed E) the hydrostatic march is stable only outward: a rho error
  // grows inward relative to rho like (rho_in/rho_out)^(Gamma/(1 - Gamma)) (from the
  // star_ic_balance anchor at 342 Rsun inward: 1 %/cell in the CZ, then NaN, run d4)
  const int ianc = is;
  ++star_eq_it_;
  std::vector<Real> pc0, pf0, p0c, p0f, x1v, x1f, ar0;
  StarAtmRadial(pmbp, pc0, pf0, p0c, p0f, x1v, x1f, ar0);
  auto hu = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.u0);
  auto hp0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), fl.phicc0);
  auto hf0 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pvet->f0x1);
  auto har = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->area.x1f);
  const bool hdm = (pvet->icr_dm.extent_int(0) == nmb);
  if (!hdm || !(pvet->icr_dm_dt > 0.0)) {
    StarFatal("star_ic_equil: no force pass from <rad_vet>/ic_rad", __LINE__);
  }
  auto hdmv = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pvet->icr_dm);
  const Real dtm = pvet->icr_dm_dt;
  // shell sums per x1 index: rho, T, the net deposit dm, count; per x1 face: F0 area
  std::vector<Real> sh(4*n1, 0.0), lf(n1 + 1, 0.0), af(n1 + 1, 0.0);
  // <mesh>/layout_invariant: the shell sums order-independent (ReproSum);
  // slots [0, 4 n1) sh, [4 n1, 5 n1 + 1) lf, [5 n1 + 1, 6 n1 + 2) af
  const bool lay = pm->layout_invariant;
  repro::ReproVec rv(lay ? 6*n1 + 2 : 0);
  const int olf = 4*n1, oaf = 5*n1 + 1;
  for (int m=0; m<nmb; ++m) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          const Real d = hu(m,IDN,k,j,i);
          const Real ek = 0.5*(SQR(hu(m,IM1,k,j,i)) + SQR(hu(m,IM2,k,j,i)) +
                               SQR(hu(m,IM3,k,j,i)))/d;
          const Real e = hu(m,IEN,k,j,i) - ek - d*hp0(m,k,j,i);
          if (lay) {
            rv.Add(4*i, d);
            rv.Add(4*i + 1, eos.Temperature(d, e));
            rv.Add(4*i + 2, hdmv(m,k,j,i));
            rv.Add(4*i + 3, 1.0);
          } else {
            sh[4*i] += d;
            sh[4*i + 1] += eos.Temperature(d, e);
            sh[4*i + 2] += hdmv(m,k,j,i);
            sh[4*i + 3] += 1.0;
          }
        }
        for (int i=is; i<=ie+1; ++i) {
          if (lay) {
            rv.Add(olf + i, hf0(m,k,j,i)*har(m,k,j,i));
            rv.Add(oaf + i, har(m,k,j,i));
          } else {
            lf[i] += hf0(m,k,j,i)*har(m,k,j,i);
            af[i] += har(m,k,j,i);
          }
        }
      }
    }
  }
  if (lay) {
    rv.Merge();
    for (int q=0; q<4*n1; ++q) {sh[q] = rv.Value(q);}
    for (int i=0; i<=n1; ++i) {
      lf[i] = rv.Value(olf + i);
      af[i] = rv.Value(oaf + i);
    }
  }
#if MPI_PARALLEL_ENABLED
  if (!lay) {
    MPI_Allreduce(MPI_IN_PLACE, sh.data(), 4*n1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, lf.data(), n1 + 1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, af.data(), n1 + 1, MPI_ATHENA_REAL, MPI_SUM,
                  MPI_COMM_WORLD);
  }
#endif
  // the 1-D column on x1 indices: rho, T, e, the force per volume f = dm/dt + rho a_ref
  std::vector<Real> bd(n1, 0.0), be(n1, 0.0), tc(n1, 0.0), d0(n1, 0.0), fv(n1, 0.0);
  std::vector<Real> resm(n1, 0.0), resh(n1, 0.0);
  for (int i=is; i<=ie; ++i) {
    const Real cnt = sh[4*i + 3];
    bd[i] = sh[4*i]/cnt;
    tc[i] = sh[4*i + 1]/cnt;
    d0[i] = bd[i];
    const Real dmi = sh[4*i + 2]/cnt;
    fv[i] = dmi/dtm + bd[i]*ar0[i];
    resm[i] = std::fabs(dmi/dtm)/(bd[i]*star_gm_/SQR(x1v[i]));
  }
  // ghosts as StarBC makes them: inner = the balanced profile scaled to the edge cell,
  // outer = the edge cell (the floor; the march stops below it)
  std::vector<Real> gbd(n1), gbe(n1);
  {
    auto h1 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_bd_);
    auto h2 = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_be_);
    for (int i=0; i<n1; ++i) {
      gbd[i] = h1(i);
      gbe[i] = h2(i);
    }
  }
  auto setcell = [&](const int i) {
    Real p, cr, ct, cv;
    eos.ThermoAt(bd[i], tc[i], be[i], p, cr, ct, cv);
  };
  for (int i=is; i<=ie; ++i) setcell(i);
  auto ghosts = [&]() {
    for (int g=0; g<is; ++g) {
      bd[g] = gbd[g]*bd[is]/gbd[is];
      be[g] = gbe[g]*be[is]/gbe[is];
      tc[g] = eos.Temperature(bd[g], be[g]);
    }
    for (int g=ie+1; g<n1; ++g) {
      bd[g] = bd[ie];
      be[g] = be[ie];
      tc[g] = tc[ie];
    }
  };
  ghosts();
  std::vector<Real> pc, pf;
  auto walk = [&](const int i, const Real di, Real &pimh, Real &piph) {
    Real ei, p, cr, ct, cv;
    eos.ThermoAt(di, tc[i], ei, p, cr, ct, cv);
    WBState s0, s1, s2, s3, s4;
    WBBackgroundStencil(eos, wbo, bd[i-1], di, bd[i+1], be[i-1], ei, be[i+1],
                        pc[i-1], pf[i], pc[i], pf[i+1], pc[i+1], s0, s1, s2, s3, s4,
                        tc[i-1], tc[i], tc[i+1]);
    pimh = s1.p;
    piph = s3.p;
  };
  // the floor: the cells at the column's floor density are not hydrostatic (sponge)
  auto above = [&](const int i) { return bd[i] > dfl*(1.0 + 1.0e-6); };
  // residual 2: the WB face condition at the current state with the current Phi_eff,
  // as an acceleration (Delta p / (rho_face dr)) / g, interior faces below the floor
  pc = pc0;
  pf = pf0;
  Real rmm = 0.0, rmh = 0.0, rrm = 0.0, rrh = 0.0;
  for (int i=is; i<ie; ++i) {
    if (!above(i) || !above(i+1)) continue;
    Real pl0, pr0, pl1, pr1;
    walk(i, bd[i], pl0, pr0);
    walk(i+1, bd[i+1], pl1, pr1);
    const Real g = star_gm_/SQR(x1f[i+1]);
    resh[i] = std::fabs(pr0 - pl1)/(0.5*(bd[i] + bd[i+1])*(x1v[i+1] - x1v[i])*g);
    if (resh[i] > rmh) {
      rmh = resh[i];
      rrh = x1f[i+1];
    }
  }
  for (int i=is; i<=ie; ++i) {
    if (!above(i)) continue;
    if (resm[i] > rmm) {
      rmm = resm[i];
      rrm = x1v[i];
    }
  }
  // L_rad + L_MLT (w = 1) through every x1 face against L_in (the wall face)
  std::vector<Real> lm(n1 + 1, 0.0);
  if (star_mlt_ && star_fm_.extent_int(0) >= n1 + 1) {
    auto hfm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_fm_);
    for (int i=is; i<=ie+1; ++i) lm[i] = hfm(i)*af[i];
  }
  const Real lin = lf[is] + lm[is];
  Real ldv = 0.0, rld = 0.0;
  for (int i=is; i<=ie+1; ++i) {
    const Real dv = std::fabs((lf[i] + lm[i])/lin - 1.0);
    if (dv > ldv) {
      ldv = dv;
      rld = x1f[i];
    }
  }
  const bool conv = (std::max(rmm, rmh) < star_eq_tol_);
  const bool lastm = !conv && star_eq_last_ && (star_eq_it_ >= star_eq_maxit_);
  const bool done = conv || (star_eq_it_ >= star_eq_maxit_);
  if (global_variable::my_rank == 0) {
    std::printf("star_vet: star_ic_equil round %d: module residual max |a_net|/g = %.3e "
                "(r %.5e), WB face residual %.3e (r %.5e), max |L(face)/L_in - 1| = "
                "%.3e (r %.5e), L_in %.10e, L_top %.10e%s\n", star_eq_it_, rmm, rrm, rmh,
                rrh, ldv, rld, lin, lf[ie+1] + lm[ie+1],
                done ? ((std::max(rmm, rmh) < star_eq_tol_) ? " CONVERGED" :
                        (lastm ? " (star_ic_equil_maxit: last round, march then STOP)" :
                         " (star_ic_equil_maxit: STOP, not converged)")) : "");
    // the profile of this round (appended): i r rho T a_ref f/rho resm resh L/L_in
    FILE *fp = std::fopen("star_ic_equil.txt", "a");
    if (fp != nullptr) {
      for (int i=is; i<=ie; ++i) {
        std::fprintf(fp, "%d %d %.12e %.12e %.12e %.12e %.12e %.6e %.6e %.12e\n",
                     star_eq_it_, i - is, x1v[i], bd[i], tc[i], ar0[i], fv[i]/bd[i],
                     resm[i], resh[i], (lf[i] + lm[i])/lin);
      }
      std::fclose(fp);
    }
  }
  // the end of the equilibrium: the deferred seed, a_ref kept for the restarts (af0 = the
  // a_ref of the cells is..ie now in place), no further ic_rad
  auto finish = [&](const std::vector<Real> &af0) {
    // the deferred star_seed: eint (and with star_seed_rad E) of every cell times its
    // factor
    if (star_eq_seed_) {
      auto hsf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_sfac_);
      auto hr = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pvet->u0);
      for (int m=0; m<nmb; ++m) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              const Real fac = hsf(m,k,j,i);
              if (fac == 1.0) continue;
              const Real d = hu(m,IDN,k,j,i);
              const Real ek = 0.5*(SQR(hu(m,IM1,k,j,i)) + SQR(hu(m,IM2,k,j,i)) +
                                   SQR(hu(m,IM3,k,j,i)))/d;
              const Real e = hu(m,IEN,k,j,i) - ek - d*hp0(m,k,j,i);
              hu(m,IEN,k,j,i) = e*fac + ek + d*hp0(m,k,j,i);
              if (star_eq_srad_) hr(m,radvet::VET_E,k,j,i) *= SQR(SQR(fac));
            }
          }
        }
      }
      Kokkos::deep_copy(fl.u0, hu);
      Kokkos::deep_copy(pvet->u0, hr);
    }
    star_eq_a_.assign(af0.begin(), af0.end());
    pvet->icr_fmeas = false;
    pm->pgen->user_icpost_again = false;
    if (global_variable::my_rank == 0) {
      std::printf("star_vet: star_ic_equil: done after %d rounds%s\n", star_eq_it_,
                  star_eq_seed_ ? " (deferred star_seed applied)" : "");
    }
  };
  // not converged in star_ic_equil_maxit rounds is fatal, unless the last round marches
  // (star_ic_equil_last_march); then that march must converge (below)
  if (done && !conv && !lastm) {
    StarFatal("problem/star_ic_equil: not converged in star_ic_equil_maxit rounds",
              __LINE__);
  }
  if (done && !lastm) {
    finish(std::vector<Real>(ar0.begin() + is, ar0.begin() + ie + 1));
    return;
  }
  // ---- the update: a = f/rho, Phi_eff from a above x1f(is), rho re-marched
  std::vector<Real> a(ie - is + 1), an(ie - is + 1);
  for (int i=is; i<=ie; ++i) a[i - is] = above(i) ? fv[i]/bd[i] : ar0[i];
  std::vector<Real> ap(a);   // a of the previous pass (convergence of the passes)
  int ilast = ie, npass = 0, nstop = 0, istop = -1;   // (nstop: passes stopped early)
  Real chg = 0.0, rmx = 0.0, rstop = 0.0;
  // Phi_eff of cell i for a trial a_i = f_i/d_i, the march's known face kept: inward
  // (face x1f(i+1) known) pc[i], pf[i] and pc[i-1] shifted with pf[i]; outward (face
  // x1f(i) known) pc[i], pf[i+1] and pc[i+1] shifted with pf[i+1].  Only differences
  // of Phi_eff enter the walk, so the frame below (inward) / above (outward) is carried
  // by these shifts; the final StarAtmPhiEff rebuilds the same differences.  Solving a_i
  // WITH d_i is what makes the march converge: lagged (a = f/rho of the previous pass)
  // it amplifies a global rho error by ~ sum Gamma dr/H_p > 1 (diverged, run d1)
  auto setphi = [&](const int i, const Real ai, const bool outward) {
    if (outward) {
      const Real sh0 = pf[i+1];
      pc[i] = pf[i] + (p0c[i] - p0f[i]) - ai*(x1v[i] - x1f[i]);
      pf[i+1] = pf[i] + (p0f[i+1] - p0f[i]) - ai*(x1f[i+1] - x1f[i]);
      pc[i+1] += pf[i+1] - sh0;
    } else {
      const Real sh0 = pf[i];
      pc[i] = pf[i+1] - (p0f[i+1] - p0c[i]) + ai*(x1f[i+1] - x1v[i]);
      pf[i] = pc[i] - (p0c[i] - p0f[i]) + ai*(x1v[i] - x1f[i]);
      pc[i-1] += pf[i] - sh0;
    }
  };
  auto solve = [&](const int i, const Real target, const bool outward) {
    // outward: PL_i(d_i) = target; inward: PR_i(d_i) = target; a_i = f_i/d_i
    auto res = [&](const Real lx) {
      Real q1, q2;
      const Real d = exp(lx);
      setphi(i, fv[i]/d, outward);
      walk(i, d, q1, q2);
      return (outward ? q1 : q2)/target - 1.0;
    };
    Real x0 = log(bd[i]), x1 = x0 + 1.0e-6;
    Real fa = res(x0), fb = res(x1);
    for (int it=0; it<60 && fabs(fb) > 1.0e-15 && fb != fa; ++it) {
      const Real x2 = x1 - fb*(x1 - x0)/(fb - fa);
      x0 = x1; fa = fb; x1 = x2;
      fb = res(x1);
    }
    setphi(i, fv[i]/exp(x1), outward);   // the Phi_eff of the solution
    a[i - is] = fv[i]/exp(x1);
    return std::make_pair(x1, fb);
  };
  for (int pass=0; pass<200; ++pass) {
    ++npass;
    StarAtmPhiEff(is, ie, a, pc0, pf0, p0c, p0f, x1v, x1f, pc, pf);
    rmx = 0.0;
    ghosts();
    // outward from the wall cell to the floor
    ilast = ianc;
    for (int i=ianc+1; i<=ie; ++i) {
      if (!(d0[i] > dfl*(1.0 + 1.0e-6))) break;
      Real pl, pr;
      walk(i-1, bd[i-1], pl, pr);
      const Real asv = a[i - is];
      auto r = solve(i, pr, true);
      if (!std::isfinite(r.first) || !(exp(r.first) > dfl*(1.0 + 1.0e-6)) ||
          fabs(r.second) > 1.0e-10) {
        a[i - is] = asv;
        if (!std::isfinite(r.first) || fabs(r.second) > 1.0e-10) {   // (not the floor)
          ++nstop;
          istop = i;
          rstop = r.second;
        }
        break;
      }
      rmx = std::max(rmx, fabs(r.second));
      bd[i] = exp(r.first);
      setcell(i);
      ilast = i;
    }
    for (int i=ilast+1; i<=ie; ++i) {
      bd[i] = d0[i];
      setcell(i);
    }
    ghosts();
    chg = 0.0;
    for (int i=is; i<=ie; ++i) {
      an[i - is] = (i <= ilast) ? fv[i]/bd[i] : ar0[i];
      chg = std::max(chg, fabs(an[i - is]/ap[i - is] - 1.0));
    }
    a = an;
    ap = an;
    if (chg < 1.0e-14) break;
  }
  StarAtmPhiEff(is, ie, a, pc0, pf0, p0c, p0f, x1v, x1f, pc, pf);
  // write back: rho of every cell times its shell's ratio, e from (rho, the cell's own
  // T), v kept; a_ref and Phi_eff of the last pass
  Real dmx = 0.0, rdm = 0.0;
  for (int i=is; i<=ie; ++i) {
    if (fabs(bd[i]/d0[i] - 1.0) > dmx) {
      dmx = fabs(bd[i]/d0[i] - 1.0);
      rdm = x1v[i];
    }
  }
  for (int m=0; m<nmb; ++m) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          const Real q = bd[i]/d0[i];
          if (q == 1.0) continue;
          const Real dold = hu(m,IDN,k,j,i);
          const Real v1 = hu(m,IM1,k,j,i)/dold, v2 = hu(m,IM2,k,j,i)/dold;
          const Real v3 = hu(m,IM3,k,j,i)/dold;
          const Real eo = hu(m,IEN,k,j,i) - 0.5*dold*(v1*v1 + v2*v2 + v3*v3)
                          - dold*hp0(m,k,j,i);
          const Real t = eos.Temperature(dold, eo);
          const Real d = dold*q;
          Real e, p, cr, ct, cv;
          eos.ThermoAt(d, t, e, p, cr, ct, cv);
          hu(m,IDN,k,j,i) = d;
          hu(m,IM1,k,j,i) = d*v1;
          hu(m,IM2,k,j,i) = d*v2;
          hu(m,IM3,k,j,i) = d*v3;
          hu(m,IEN,k,j,i) = e + 0.5*d*(v1*v1 + v2*v2 + v3*v3) + d*hp0(m,k,j,i);
        }
      }
    }
  }
  Kokkos::deep_copy(fl.u0, hu);
  StarAtmWrite(pmbp, is, a, pc, pf);
  pm->pgen->user_icpost_again = true;
  if (lastm) finish(a);
  if (global_variable::my_rank == 0) {
    std::printf("star_vet: star_ic_equil round %d update: anchor cell %d (r %.5e), "
                "cells %d .. %d re-marched, %d passes (last max |da/a| %.3e, face "
                "residual %.3e), max |drho/rho| %.3e at r %.5e\n", star_eq_it_, ianc - is,
                x1v[ianc], 0, ilast - is, npass, chg, rmx, dmx, rdm);
  }
  // (a march stopped below the floor: its cells above keep the column's rho and a_ref)
  if (nstop > 0 && global_variable::my_rank == 0) {
    std::printf("star_vet: WARNING star_ic_equil round %d: the march stopped at cell %d "
                "(r %.5e) in %d of %d passes, no face solve (last residual %.3e); the "
                "cells above keep the column\n", star_eq_it_, istop - is, x1v[istop],
                nstop, npass, rstop);
  }
  if (lastm && !(chg < 1.0e-12)) {
    StarFatal("problem/star_ic_equil: the last round's march did not converge", __LINE__);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void StarFinal()
//! \brief release the namespace-scope Views before Kokkos::finalize.

//----------------------------------------------------------------------------------------
//! \fn void StarRemapLoad()
//! \brief problem/star_ic_remap_file: hydro u0 (5), VET u0 (4) and f0x1 of
//! every MeshBlock of this rank, ghosts included, from the file's global arrays
//! [nx3+2ng][nx2+2ng][nx1+2ng(+1)] (header checked by the caller; uniform mesh only)

void StarRemapLoad(Mesh *pm, const std::string &fn, const std::int64_t off,
                   const int *n) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto *pvet = pmbp->pradvet;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  if (pm->multilevel) StarFatal("star_ic_remap_file needs a uniform mesh", __LINE__);
  const int nout1 = indcs.nx1 + 2*ng, nout2 = indcs.nx2 + 2*ng, nout3 = indcs.nx3 + 2*ng;
  const std::int64_t c1 = n[1] + 2*n[0], c2 = n[2] + 2*n[0], c3 = n[3] + 2*n[0];
  const std::int64_t cc = c1*c2*c3;
  std::ifstream rf(fn, std::ios::binary);
  if (!rf) StarFatal("cannot open problem/star_ic_remap_file " + fn, __LINE__);
  auto uh = Kokkos::create_mirror_view(pmbp->phydro->u0);
  auto ur = Kokkos::create_mirror_view(pvet->u0);
  const bool hf0 = (pvet->f0x1.extent_int(0) >= pmbp->nmb_thispack);
  auto uf = Kokkos::create_mirror_view(pvet->f0x1);
  std::vector<double> row(nout1 + 1);
  auto rd = [&](const std::int64_t pos, const int len) {
    rf.seekg(pos);
    rf.read(reinterpret_cast<char*>(row.data()), len*sizeof(double));
    if (!rf) StarFatal("problem/star_ic_remap_file: read failed", __LINE__);
  };
  for (int m=0; m<pmbp->nmb_thispack; ++m) {
    const LogicalLocation &ll = pm->lloc_eachmb[pmbp->gids + m];
    const std::int64_t o1 = static_cast<std::int64_t>(ll.lx1)*indcs.nx1;
    const std::int64_t o2 = static_cast<std::int64_t>(ll.lx2)*indcs.nx2;
    const std::int64_t o3 = static_cast<std::int64_t>(ll.lx3)*indcs.nx3;
    for (int v=0; v<9; ++v) {
      for (int k=0; k<nout3; ++k) {
        for (int j=0; j<nout2; ++j) {
          rd(off + 8*(v*cc + ((o3 + k)*c2 + (o2 + j))*c1 + o1), nout1);
          for (int i=0; i<nout1; ++i) {
            if (v < 5) {
              uh(m,v,k,j,i) = row[i];
            } else {
              ur(m,v-5,k,j,i) = row[i];
            }
          }
        }
      }
    }
    if (hf0) {
      for (int k=0; k<nout3; ++k) {
        for (int j=0; j<nout2; ++j) {
          rd(off + 8*(9*cc + ((o3 + k)*c2 + (o2 + j))*(c1 + 1) + o1), nout1 + 1);
          for (int i=0; i<=nout1; ++i) uf(m,k,j,i) = row[i];
        }
      }
    }
  }
  Kokkos::deep_copy(pmbp->phydro->u0, uh);
  Kokkos::deep_copy(pvet->u0, ur);
  if (hf0) Kokkos::deep_copy(pvet->f0x1, uf);
  if (global_variable::my_rank == 0) {
    std::cout << "star_vet: star_ic_remap_file: u0, VET u0" << (hf0 ? ", f0x1" : "")
              << " loaded (" << pmbp->nmb_thispack << " MeshBlock(s) on rank 0); "
              << "no ic_rad, no star_ic_equil, no seed" << std::endl;
  }
}

void StarFinal(ParameterInput *pin, Mesh *pm) {
  (void) pin; (void) pm;
  star_rho_ = DvceArray1D<Real>();
  star_eint_ = DvceArray1D<Real>();
  star_bd_ = DvceArray1D<Real>();
  star_fm_ = DvceArray1D<Real>();
  star_fsub_ = DvceArray1D<Real>();
  star_tgp_ = DvceArray1D<Real>();
  star_shs_ = DvceArray2D<Real>();
  star_shp_ = DvceArray2D<Real>();
  star_be_ = DvceArray1D<Real>();
  star_phig_ = DvceArray1D<Real>();
  star_gmr_ = DvceArray1D<Real>();
  star_fb_ = DvceArray1D<Real>();
  star_min_ = DvceArray3D<Real>();
  star_sfac_ = DvceArray4D<Real>();
  stardrive::Final();
}

//----------------------------------------------------------------------------------------
//! \fn void StarN2Compute
//! \brief the IC column's convective growth rate on the radial cells of the
//! (single, full-length) x1 MeshBlock column, from the initial column as the IC kernel
//! lays it down (star_ic_balance cells, else the fine column), once (star_n2_):
//!   T from the code EOS; P = p_gas + a T^4/3; Gamma_1 of gas + radiation from the first
//!   law, chi_rho = p_g chi_rho,g/P, chi_T = (p_g chi_T,g + 4 p_r)/P, rho T c_v = rho T
//!   c_v,g + 4 a T^4, Gamma_1 = chi_rho + P chi_T^2/(rho T c_v);  g = G m(r)/r^2 (the
//!   star_gm_column mass, else star_gm);  N^2 = g (dlnP/dr/Gamma_1 - dln rho/dr) by
//!   centred differences;  growth rate sigma = sqrt(max(-N^2, 0)), H_p = P/(rho g).
//! Domain: active cells with r <= r_int (tau = star_tau_int, the code's photosphere).
void StarN2Compute(MeshBlockPack *pmbp, const EOS_Data &eos, const Real ar) {
  StarN2Col &c2 = star_n2_;
  if (c2.done) return;
  Mesh *pm = pmbp->pmesh;
  auto &indcs = pm->mb_indcs;
  if (indcs.nx1 != pm->mesh_indcs.nx1) {
    StarFatal("star_seed_band = n2 / mlt_ramp auto need one MeshBlock across x1 (or "
              "star_seed_band = old and explicit mlt_ramp_start/_time)", __LINE__);
  }
  const int is = indcs.is, ie = indcs.ie;
  auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
  auto hrho = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_rho_);
  auto hein = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_eint_);
  const int nf = star_nf_;
  auto lint = [&](const decltype(hrho) &a, const Real r, const bool lg) {
    Real x = (r - star_rlo_)/star_dr_;
    int i = static_cast<int>(floor(x));
    i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
    Real w = x - i;
    w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
    return lg ? exp((1.0 - w)*log(a(i)) + w*log(a(i+1))) : (1.0 - w)*a(i) + w*a(i+1);
  };
  decltype(hrho) hbd, hbe, hgm;
  if (star_bal_) {
    hbd = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_bd_);
    hbe = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_be_);
  }
  if (star_gmc_) hgm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), star_gmr_);
  const int n1 = ie + 2;
  std::vector<Real> &r = c2.r, &g1 = c2.g1, &gg = c2.gg, &pt = c2.pt, &dd = c2.dd,
                    &tt = c2.tt, &sg = c2.sg, &dm = c2.dm;
  r.assign(n1, 0.0); g1.assign(n1, 0.0); gg.assign(n1, 0.0); pt.assign(n1, 0.0);
  dd.assign(n1, 0.0); tt.assign(n1, 0.0); sg.assign(n1, 0.0); dm.assign(n1, 0.0);
  std::vector<Real> lp(n1), ld(n1);
  for (int i = is - 1; i <= ie + 1; ++i) {
    r[i] = hx1v(0,i);
    const Real d = star_bal_ ? hbd(i) : lint(hrho, r[i], true);
    const Real e = star_bal_ ? hbe(i) : lint(hein, r[i], true);
    const Real t = eos.Temperature(d, e);
    Real eg, pg, crg, ctg, cvg;
    eos.ThermoAt(d, t, eg, pg, crg, ctg, cvg);
    const Real pr = ar*SQR(SQR(t))/3.0, p = pg + pr;
    const Real chr = pg*crg/p, cht = (pg*ctg + 4.0*pr)/p;
    g1[i] = chr + p*cht*cht/(d*t*cvg + 12.0*pr);
    gg[i] = (star_gmc_ ? lint(hgm, r[i], false) : star_gm_)/SQR(r[i]);
    lp[i] = log(p);
    ld[i] = log(d);
    pt[i] = p;
    dd[i] = d;
    tt[i] = t;
    dm[i] = d*SQR(r[i])*(hx1f(0,i+1) - hx1f(0,i));
  }
  c2.is = is;
  c2.iex = is - 1;   // last active cell inside r_int
  c2.smax = 0.0;
  c2.nu = 0;
  for (int i = is; i <= ie; ++i) {
    if (r[i] > star_rint_) break;
    c2.iex = i;
    const Real h = r[i+1] - r[i-1];
    const Real n2 = gg[i]*((lp[i+1] - lp[i-1])/h/g1[i] - (ld[i+1] - ld[i-1])/h);
    sg[i] = (n2 < 0.0) ? sqrt(-n2) : 0.0;
    if (sg[i] > c2.smax) {
      c2.smax = sg[i];
      c2.ipk = i;
    }
    if (sg[i] > 0.0) {
      if (c2.nu == 0) c2.iu0 = i;
      c2.iu1 = i;
      ++c2.nu;
    }
  }
  c2.done = true;
}

//----------------------------------------------------------------------------------------
//! \fn Real StarN2Median
//! \brief the mass-weighted (rho r^2 dr) median of 1/sigma over the cells
//! below r_int with N^2 < 0 (star_n2_), for mlt_ramp_time = auto.
Real StarN2Median(MeshBlockPack *pmbp, const EOS_Data &eos, const Real ar) {
  StarN2Compute(pmbp, eos, ar);
  StarN2Col &c2 = star_n2_;
  std::vector<std::pair<Real, Real>> q;
  Real mt = 0.0;
  for (int i = c2.is; i <= c2.iex; ++i) {
    if (c2.sg[i] > 0.0) {
      q.emplace_back(1.0/c2.sg[i], c2.dm[i]);
      mt += c2.dm[i];
    }
  }
  if (q.empty()) return 0.0;
  std::sort(q.begin(), q.end());
  Real acc = 0.0;
  c2.med = q.back().first;
  for (auto &p : q) {
    acc += p.second;
    if (acc >= 0.5*mt) {
      c2.med = p.first;
      break;
    }
  }
  return c2.med;
}

//----------------------------------------------------------------------------------------
//! \fn bool StarN2Band
//! \brief the radial band and Y_lm l range of the star_seed_band = n2 seed,
//! from the growth rates of StarN2Compute.
//! Band rule: cells with sigma >= frac sigma_max; gaps of <= ngap cells between such
//! cells are closed; connected runs shorter than nmin cells (isolated IC-column ripples)
//! are discarded; the band is the surviving run with the largest sigma (ties: the outer
//! one).  rlo/rhi = the centres of the cells just outside the run (the sin envelope then
//! vanishes there).  l: l_Hp = pi r/H_p at the run's peak cell; the grid limit
//! l_grid = (cells around a great circle)/ppw (cubed sphere 4 nx2; spherical polar 2 pi
//! over the smaller of the theta, phi spacings); lmax = min(l_Hp, l_grid, lcap) (>= 2),
//! lmin = max(1, lmax/2).  Every rank computes the same band (no communication).
//! Returns false (no seed) when no cell is unstable.
bool StarN2Band(MeshBlockPack *pmbp, const EOS_Data &eos, const Real ar, const Real tk,
                const bool hcs, const Real frac, const int nmin, const int ngap,
                const Real ppw, const int lcap, Real &rlo, Real &rhi, int &lmin,
                int &lmax) {
  StarN2Compute(pmbp, eos, ar);
  Mesh *pm = pmbp->pmesh;
  const int is = pm->mb_indcs.is, ie = pm->mb_indcs.ie;
  auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->x1v);
  const StarN2Col &c2 = star_n2_;
  const std::vector<Real> &r = c2.r, &g1 = c2.g1, &gg = c2.gg, &pt = c2.pt, &dd = c2.dd,
                            &tt = c2.tt, &sg = c2.sg;
  const int n1 = ie + 2, iex = c2.iex;
  const Real smax = c2.smax;
  const bool root = (global_variable::my_rank == 0);
  if (!(smax > 0.0)) {
    if (root) {
      std::cout << "star_vet: star_seed_band = n2: N^2 >= 0 in every cell below r_int "
                << "= " << star_rint_ << ": NO SEED" << std::endl;
    }
    return false;
  }
  std::vector<int> c(n1, 0);
  for (int i = is; i <= iex; ++i) c[i] = (sg[i] >= frac*smax) ? 1 : 0;
  // close gaps of <= ngap cells between candidate cells
  for (int i = is; i <= iex; ++i) {
    if (c[i] != 1) continue;
    int j = i + 1;
    while (j <= iex && c[j] == 0) ++j;
    if (j <= iex && j - i - 1 <= ngap) {
      for (int q = i + 1; q < j; ++q) c[q] = 2;
    }
  }
  // connected runs: keep those of >= nmin cells; the band = the one with the largest peak
  int bs = -1, be = -1, bp = -1;
  Real bpk = -1.0;
  if (root) {
    std::cout << "star_vet: star_seed_band = n2: sigma_max (all cells below r_int) = "
              << smax << " /s (1/sigma " << 1.0/smax/86400.0 << " d); runs with sigma >= "
              << frac << " sigma_max (gaps <= " << ngap << " closed):" << std::endl;
  }
  for (int i = is; i <= iex; ++i) {
    if (c[i] == 0 || (i > is && c[i-1] != 0)) continue;
    int j = i;
    while (j + 1 <= iex && c[j+1] != 0) ++j;
    int ip = i;
    for (int q = i; q <= j; ++q) {
      if (sg[q] > sg[ip]) ip = q;
    }
    const bool keep = (j - i + 1 >= nmin);
    if (root) {
      std::cout << "  run r/r_int " << r[i]/star_rint_ << " - " << r[j]/star_rint_ << " ("
                << j - i + 1 << " cells), peak 1/sigma " << 1.0/sg[ip]/86400.0
                << " d at r/r_int " << r[ip]/star_rint_
                << (keep ? "" : "  [< star_seed_n2_nmin cells: discarded]") << std::endl;
    }
    if (keep && sg[ip] >= bpk) {
      bs = i;
      be = j;
      bp = ip;
      bpk = sg[ip];
    }
  }
  if (bs < 0) {
    if (root) {
      std::cout << "star_vet: star_seed_band = n2: no run of >= " << nmin
                << " cells: NO SEED" << std::endl;
    }
    return false;
  }
  rlo = hx1v(0,bs-1);   // envelope zero in the cells just outside the band
  rhi = hx1v(0,be+1);
  const Real hp = pt[bp]/(dd[bp]*gg[bp]);
  const Real lhp = M_PI*r[bp]/hp;
  // l window: l_hi = min(pi r/H_p, l_8), l_8 = 2 pi r/(ppw dx_lat) (ppw = 8:
  // pi r/(4 dx_lat)), dx_lat = r dangle with dangle the SMALLEST angular distance between
  // neighbouring cell centres (cubed sphere: over a whole equiangular panel; spherical
  // polar: min(dtheta, sin(theta) dphi) over the theta centres); l_lo = max(1, round(l_hi
  // /4)) (two octaves); l_hi above lcap is clipped to lcap with a warning.
  Real dang = 1.0e30;
  if (hcs) {
    const int nn = pm->mesh_indcs.nx2;
    const Real dxi = 0.5*M_PI/nn;
    auto uv = [&](const int a, const int b, Real *u) {
      const Real x = tan(-0.25*M_PI + (a + 0.5)*dxi), y = tan(-0.25*M_PI + (b + 0.5)*dxi);
      const Real n = sqrt(1.0 + x*x + y*y);
      u[0] = 1.0/n; u[1] = x/n; u[2] = y/n;
    };
    for (int a = 0; a < nn; ++a) {
      for (int b = 0; b + 1 < nn; ++b) {
        Real u[3], v[3], w[3];
        uv(a, b, u);
        uv(a, b + 1, v);
        uv(b, a, w);
        const Real c1 = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
        dang = std::min(dang, acos(fmin(c1, 1.0)));
        uv(b + 1, a, v);
        const Real c2 = w[0]*v[0] + w[1]*v[1] + w[2]*v[2];
        dang = std::min(dang, acos(fmin(c2, 1.0)));
      }
    }
  } else {
    const int n2 = pm->mesh_indcs.nx2;
    const Real dth = (pm->mesh_size.x2max - pm->mesh_size.x2min)/n2;
    const Real dph = (pm->mesh_size.x3max - pm->mesh_size.x3min)/pm->mesh_indcs.nx3;
    dang = dth;
    for (int j = 0; j < n2; ++j) {
      const Real th = pm->mesh_size.x2min + (j + 0.5)*dth;
      dang = std::min(dang, sin(th)*dph);
    }
  }
  const Real l8 = 2.0*M_PI/(ppw*dang);
  Real lhi = std::floor(std::min(lhp, l8));
  bool clip = false;
  if (lhi > lcap) {
    lhi = lcap;
    clip = true;
  }
  lmax = std::max(static_cast<int>(lhi), 1);
  lmin = std::max(1, static_cast<int>(std::lround(lmax/4.0)));
  if (root) {
    std::cout << "star_vet: star_seed_band = n2 BAND r " << rlo << " .. " << rhi
              << " (r/r_int " << rlo/star_rint_ << " - " << rhi/star_rint_ << ", r_int "
              << star_rint_ << "), " << be - bs + 1 << " cells " << bs - is << ".."
              << be - is
              << "; peak at r/r_int " << r[bp]/star_rint_ << ": sigma " << bpk
              << " /s, 1/sigma " << 1.0/bpk/86400.0 << " d, T " << tt[bp]*tk
              << " K, Gamma1 " << g1[bp] << ", H_p/r " << hp/r[bp]
              << "; pi r/H_p " << lhp << ", l_8 " << l8 << " (" << ppw
              << " cells per wavelength, min dx_lat/r " << dang << ") -> l_lo " << lmin
              << ", l_hi " << lmax << std::endl;
    if (clip) {
      std::cout << "### WARNING star_vet: star_seed_band = n2: l_hi clipped to the cap "
                << lcap << std::endl;
    }
  }
  return true;
}
}  // namespace
