//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file he_star_m1.cpp
//! \brief Problem generator for the 4.0 Msun presupernova He star (Woosley 2019) on a
//! spherical-polar WEDGE with the implicit grey M1 module (<rad_m1>), used with
//! `cmake -D PROBLEM=he_star_m1`.  Input: inputs/radiation/he_presn_m1_wedge.athinput,
//! initial column: docs/handover/scripts/make_ic_he_presn_m1.py.
//!
//! LEAN by design: a point mass, the well-balanced x1 gravity, the <rad_m1> hookup and
//! nothing else.  No two-stream, no MLT closure, no rt_* keys, no conduction operator,
//! no MHD (refused).  Derived from the sph_wedge test pgen (tests/rad_m1_wedge.cpp,
//! which still serves grey and box-sized tests) with the features this star does not
//! need removed and the ones it does need added:
//!   * the gas-only tabulated EOS (<hydro>/eos_radiation = false, enforced by the module)
//!     or an ideal gas (eos = ideal; <units>/mu sets its kelvin scale, which
//!     <rad_m1>/temp_unit_kelvin must equal; no EOS range check then)
//!     and the Rosseland + Planck tables (problem/he_opac_table, he_planck_table);
//!   * the initial column he_ic_file, columns  r rho eint F_r  (cgs, ascending r):
//!     F_r is the DIFFUSIVE part of the luminosity flux, F_r = L/(4 pi r^2) - F_MLT, so
//!     the radiation force kappa_R F_r/c that the gas is in balance with is known
//!     everywhere; E = a T(rho,eint)^4 is set here from the run's own EOS;
//!   * force_reference = wb_arad with the reference acceleration a_ref = kappa_t F_r/c
//!     of that column, delivered through the EFFECTIVE potential of the x1
//!     well-balanced pair Phi_eff = Phi - int a_ref dr (<rad_m1> then applies only the
//!     residual).  The pgen also gives the reference WORK (fref_wsplit_ok), which is what
//!     lets <rad_m1>/force_reference_work = split be accepted;
//!   * x1 boundaries (ix1_bc = ox1_bc = user): inner = closed hydro wall carrying the
//!     scaled initial profile, the luminosity entering through the M1 face flux
//!     (<rad_m1>/implicit_bc_x1min = flux, implicit_flux_x1min = L/(4 pi r_in^2));
//!     outer = hydro outflow with NO inflow (mass may leave through the top, none
//!     enters), M1 Marshak top (implicit_bc_x1max = marshak);
//!   * a startup CHECK that the whole initial column, ghost margin included, lies inside
//!     the EOS table (rho, T) and inside the Rosseland/Planck table grid (fatal);
//!   * history columns (user_hist): face luminosities, L_in, E_rad, the total energy of
//!     gas + radiation, the mass through both x1 faces, the interior kinetic sums, the
//!     mean v_r^2 of the wall cells and the mean Picard passes per implicit solve since
//!     the previous history output;
//!   * an optional multi-mode entropy/temperature seed (he_seed_nk, he_seed_rad,
//!     he_seed_signs, he_seed_kdist).
//!   * (hegiant-1006, problem/he_gm_column, default off) the gravity of the enclosed
//!     mass of the initial column instead of the point mass: G m(r) = he_gm + G int 4 pi
//!     r^2 rho_IC dr from r_in (he_gm = G m(r_in)), a static monopole frozen at t = 0.
//!
//! RESTARTS: the pgen carries NO state that is not recomputed here.  It is not skipped on
//! a restart: the column, the tables, the potentials, the reference acceleration and the
//! hooks are rebuilt from the input file (the evolved arrays come from the restart file).
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
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/coordinates.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "utils/wb_background.hpp"
#include "coordinates/cubed_sphere.hpp"
#include "coordinates/cell_locations.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_implicit.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"

// defaults-1002: <hydro>/fofc defaults to TRUE for this problem generator (every
// stellar/box input ran it explicitly); see hydro::fofc_pgen_default (hydro.hpp).
namespace {
[[maybe_unused]] const bool kFofcPgenDefault = (hydro::fofc_pgen_default = true);
}  // namespace

namespace {
// the initial column on a uniform fine radial grid (device), read by the x1 BCs
DvceArray1D<Real> hs_rho_, hs_eint_;
Real hs_rlo_ = 0.0, hs_dr_ = 1.0;
int hs_nf_ = 0;
Real hs_gm_ = 0.0, hs_rin_ = 1.0, hs_rint_ = 0.0, hs_fin_ = 0.0;
// hegiant-1006: problem/he_gm_column (read only when named, default false = the point
// mass GM = he_gm, bit for bit): the gravity of the ENCLOSED MASS of the initial column,
// G m(r) = he_gm + G int_{r_in}^{r} 4 pi r'^2 rho_IC dr' (he_gm = G m(r_in), the mass
// inside the wall), FROZEN at t = 0 (a static monopole, rebuilt from he_ic_file on every
// start).  hs_phig_ = Phi(r) = int_{r_in}^{r} G m/r'^2 dr' and hs_gmr_ = G m(r) on the
// uniform fine grid of the column (hs_rlo_, hs_dr_, hs_nf_).
bool hs_gmc_ = false;
DvceArray1D<Real> hs_phig_, hs_gmr_;
Real hs_sp_rate_ = 0.0, hs_sp_r0_ = 0.0, hs_rtop_ = 1.0;
bool hs_zflux_ = true;
// he-wind-bc: problem/he_bc_inner = inflow (hs_binf_) and he_bc_outer = outflow
// (hs_bout_); defaults wall / noinflow = the he-presn-m1 behaviour
bool hs_binf_ = false, hs_bout_ = false;
// bsg-arm2: problem/he_bc_outer = hse (hs_bhse_): the top ghosts are the HYDROSTATIC
// continuation of the last active cell (HeStarBC); he_bc_hse_gmax clamps the Eddington
// factor used there, he_bc_hse_tgrad continues T along the interior's dT/dPhi_eff
// (polytropic WB walk) instead of holding it (isothermal).  problem/he_sponge_mode =
// radial (hs_sp_rad_) makes the top sponge damp v_r only (default all = v).
// problem/he_sponge_dmax (default 0 = off): the sponge also acts, at its full rate, on
// every cell with rho < he_sponge_dmax (the floor gas above the photosphere, which can
// not be hydrostatic and falls freely)
bool hs_bhse_ = false, hs_bhse_tg_ = false, hs_sp_rad_ = false;
// problem/he_bc_hse_flux = face (default on a fresh start) | cell: the radial flux in the
// hse top's Gamma.  face = the transported comoving face flux (f0x1, mean of the top
// cell's two x1 faces); cell = the old clipped lab-frame cell F1.  A restart whose
// embedded input lacks the key keeps cell (no silent change mid-run).
bool hs_bhse_face_ = true;
Real hs_bhse_gmax_ = 0.9, hs_sp_dmax_ = 0.0;
// problem/he_ic_balance: the discretely balanced initial cells (all x1 cells incl. ghosts
// of the ONE MeshBlock along x1; every block has the same radial grid)
bool hs_bal_ = false;
DvceArray1D<Real> hs_bd_, hs_be_;
// problem/mlt_flux_frozen: the frozen MLT flux of the IC on the x1 faces (index = face)
bool hs_mlt_ = false;
DvceArray1D<Real> hs_fm_;
// problem/mlt_ramp_start, mlt_ramp_time: the time-only weight w(t) on the frozen MLT
// deposit (1 before the start, cosine to 0 over the ramp time, 0 after; start < 0 = no
// ramp, w = 1 always); hs_mw_ = the w now in <rad_m1>/esrc; hs_lmx_ = max over the x1
// faces of the UNSCALED L_MLT = F_MLT r^2 Omega (the wedge's solid angle)
Real hs_mrs_ = -1.0, hs_mrt_ = 0.0, hs_mw_ = 1.0, hs_lmx_ = 0.0;
// problem/he_esrc_const, he_esrc_rmax (the esrc test, kept when w is re-applied)
Real hs_esc_ = 0.0, hs_esrmx_ = 0.0;
// problem/mlt_closure = adaptive (default frozen): the ADAPTIVE SHELL-MEAN DEFICIT
// closure.  The applied subgrid flux F_sub (x1 faces) relaxes toward
// s(r) max(0, L - L_rad - L_conv,res)/A over problem/mlt_relax_time, updated every
// problem/mlt_closure_every cycles at the START of the cycle (before any stage), from
// shell sums over all ranks: L_rad = sum A F_0 (the M1 comoving face flux), L_conv,res =
// the FLUCTUATION flux A (<v_r'(e+p)_gas'> + (4/3) <v_r'E'>) of the cell shells averaged
// to the face.  s = clamp(fmlt_IC/problem/mlt_closure_taper, 0, 1): the FeCZ of the IC
// (F_MLT > 0) with its edges tapered.  F_sub(t=0) = F_MLT of the IC.  hs_fsub_ and the
// update clock hs_ad_tl_ are RESTART STATE (the pgen block of pgen.hpp).
bool hs_ad_ = false;
Real hs_ad_tau_ = 1500.0, hs_ad_taper_ = 0.01, hs_ad_tl_ = 0.0;
Real hs_ad_dout_ = 0.0, hs_ad_tout_ = 0.0, hs_lsub_ = 0.0;
int hs_ad_every_ = 1, hs_ad_cyc_ = -1;
DvceArray1D<Real> hs_fsub_;
DvceArray2D<Real> hs_shs_;
std::vector<Real> hs_tap_, hs_r2o_;
char hs_ad_file_[1024] = {0};
// Picard counters at the previous history output (HeStarHist)
Real hs_pic_n0_ = 0.0, hs_pic_s0_ = 0.0;

//! the frozen-MLT weight w(t): a function of the time only (restart-safe)
Real HsMltW(const Real t) {
  if (hs_mrs_ < 0.0 || t < hs_mrs_) return 1.0;
  if (hs_mrt_ <= 0.0 || t >= hs_mrs_ + hs_mrt_) return 0.0;
  return 0.5*(1.0 + cos(M_PI*(t - hs_mrs_)/hs_mrt_));
}

//! log-linear interpolation on the fine grid, clamped to its end nodes
KOKKOS_INLINE_FUNCTION
Real HsLogInterp(const DvceArray1D<Real> &a, const Real rlo, const Real dr, const int nf,
                 const Real r) {
  Real x = (r - rlo)/dr;
  int i = static_cast<int>(floor(x));
  i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
  Real w = x - i;
  w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
  return exp((1.0 - w)*log(a(i)) + w*log(a(i+1)));
}

//! linear interpolation on the fine grid, clamped
KOKKOS_INLINE_FUNCTION
Real HsLinInterp(const DvceArray1D<Real> &a, const Real rlo, const Real dr, const int nf,
                 const Real r) {
  Real x = (r - rlo)/dr;
  int i = static_cast<int>(floor(x));
  i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
  Real w = x - i;
  w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
  return (1.0 - w)*a(i) + w*a(i+1);
}

[[noreturn]] void HsFatal(const std::string &msg, const int line) {
  std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << line << std::endl
            << "he_star_m1: " << msg << std::endl;
  std::exit(EXIT_FAILURE);
}

//! the merged opacity table format of box_convection.cpp's ReadOpacityTable (comment
//! lines, one of them "# nT nD lTmin dlT lDmin dlD", then nT*nD values of log10 kappa,
//! T slowest), filled on the host and deep-copied to the device; the grid ranges are
//! returned for the startup range check
void HsReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                        DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD,
                        Real &lt_lo, Real &lt_hi, Real &ld_lo, Real &ld_hi,
                        const Real ld_ext = 1.0e300, const bool ld_hold = false) {
  std::ifstream f(fname);
  if (!f.good()) HsFatal("cannot open opacity table '" + fname + "'", __LINE__);
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
    HsFatal("opacity table '" + fname + "' has no grid line or the wrong value count",
            __LINE__);
  }
  // problem/he_opac_logd_min (he-wind-bc, default none): the grid is EXTENDED to lower
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
      // problem/he_opac_extend_hold (gd_physfix_1005): no extrapolation in rho, the
      // extended columns hold the table's edge value (kappa constant below its rho edge)
      if (ld_hold) sl = 0.0;
      for (int j=0; j<nadd; ++j) vn[i*nDn + j] = v0 - sl*(nadd - j)*dld;
      for (int j=0; j<nD; ++j) vn[i*nDn + nadd + j] = vals[i*nD + j];
    }
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: opacity table '" << fname << "' extended from log10 rho "
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
}  // namespace

namespace {
void HeStarGravity(Mesh *pm, const Real bdt);
void HeStarBC(Mesh *pm);
void HeStarHist(HistoryData *pdata, Mesh *pm);
void HsApplyMltW(Mesh *pm, const Real w);
void HsAdaptiveUpdate(Mesh *pm);
std::vector<char> HsAdaptiveRstWrite();
void HeStarFinal(ParameterInput *pin, Mesh *pm);
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief see the file header.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto *pm1 = pmbp->pradm1;
  // STAGE CS3 (m1-cs-implicit): the cubed sphere as well (x1 = r on both meshes; the
  // column, the potentials, the WB pair and the x1 BCs are radial).  On cs the solid
  // angle is 4 pi and the multi-mode seed is laid down in the global (theta, phi) of the
  // cell centre (see the seed below).
  const bool hcs = pmy_mesh_->use_cubed_sphere;
  if (pm1 == nullptr || pmbp->phydro == nullptr ||
      !(pmy_mesh_->use_spherical_polar || hcs)) {
    HsFatal("needs <hydro>, <rad_m1> and mesh/use_spherical_polar or use_cubed_sphere",
            __LINE__);
  }
  if (pmbp->pmhd != nullptr) {
    HsFatal("<mhd> is not supported by this pgen (hydro only)", __LINE__);
  }
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto *ph = &fl;
  const Real cl = pm1->c_light, efl = pm1->e_floor, ar = pm1->arad;
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
  const Real hs_tk = eos.tbl.active ? eos.temp_cgs :
                     ((pmbp->punit != nullptr) ? pmbp->punit->temperature_cgs() : 1.0);
  const bool etg = ph->use_etotgrav;
  const bool wbdyn = ph->use_wellbalance_dynamic;
  if (!(wbdyn && ph->use_wb_x1 && etg)) {
    HsFatal("needs <hydro>/etotgrav = true, wellbalance_dynamic = true and wb_x1 = true "
            "(the gravity source reads the x1 well-balanced cache)", __LINE__);
  }
  if (pm1->opacity_type != radm1::M1_OPAC_TABLE) {
    HsFatal("needs <rad_m1>/opacity = table", __LINE__);
  }
  if (pm1->force_ref != radm1::M1_FREF_WB_ARAD) {
    HsFatal("needs <rad_m1>/force_reference = wb_arad (the well-balanced pair carries "
            "the reference force of the initial column)", __LINE__);
  }

  // ---- parameters
  hs_gm_ = pin->GetReal("problem","he_gm");
  hs_gmc_ = pin->DoesParameterExist("problem","he_gm_column") &&
            pin->GetBoolean("problem","he_gm_column");
  hs_rin_ = pmy_mesh_->mesh_size.x1min;
  const Real rtop = pmy_mesh_->mesh_size.x1max;
  hs_fin_ = pin->GetReal("rad_m1","implicit_flux_x1min");
  const std::string fn = pin->GetString("problem","he_ic_file");
  const int nf = pin->GetOrAddInteger("problem","he_nfine",16384);
  const Real tau_int = pin->GetOrAddReal("problem","he_tau_int",1.0);
  // he_ic_cols = 5: the file also carries the radiation energy density E (a state that
  // was relaxed by a run); 4 (default): E = a T(rho,eint)^4
  const int ncols = pin->GetOrAddInteger("problem","he_ic_cols",4);
  if (ncols != 4 && ncols != 5) HsFatal("problem/he_ic_cols must be 4 or 5", __LINE__);
  // problem/mlt_flux_frozen (default false; in 3-D only as the start-up scaffold that is
  // ramped off by problem/mlt_ramp_start/_time while convection grows, or relaxed by
  // problem/mlt_closure = adaptive): the MLT flux of the IC, F_MLT = F_r fmlt/(1 - fmlt)
  // with fmlt = F_MLT/F the file's column 7 (make_ic mlt: r rho eint F_r E T fmlt), so
  // that F_r + F_MLT = L/(4 pi r^2) exactly, is deposited as the conservative energy
  // source -(A F_MLT|_{i+1/2} - A F_MLT|_{i-1/2})/V into the RADIATION energy through
  // <rad_m1>/esrc (the implicit solve's old vector, every stage with its weight).
  // Rebuilt from the file on every start (restarts).
  hs_mlt_ = pin->GetOrAddBoolean("problem","mlt_flux_frozen",false);
  if (hs_mlt_ && ncols != 5) {
    HsFatal("problem/mlt_flux_frozen = true needs he_ic_cols = 5 and the 7-column "
            "file of make_ic_he_presn_m1.py mlt", __LINE__);
  }
  // ONE MeshBlock along x1 (the margins below, the history's face picks and the balanced
  // IC read the radial grid of block 0)
  if (pmy_mesh_->mesh_indcs.nx1 != indcs.nx1) {
    HsFatal("needs ONE MeshBlock along x1 (meshblock/nx1 = mesh/nx1)", __LINE__);
  }
  // the fine grid spans the OUTERMOST GHOST FACES (the code's own x1f, stretched or
  // not) plus two end-cell widths on each side: on a uniform grid this is the old
  // (ng + 2) dx margin exactly; on a coarse stretched grid the old "at least the mean
  // width" margin reached rho < 1e-14 above the top (v2/runs/P03 fatal)
  Real rlo_g = hs_rin_, rhi_g = rtop, dxlo = 0.0, dxhi = 0.0;
  {
    auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pmbp->pcoord->xx1f);
    rlo_g = hx1f(0,0);
    rhi_g = hx1f(0,n1m1+1);
    dxlo = hx1f(0,indcs.is+1) - hx1f(0,indcs.is);
    dxhi = hx1f(0,indcs.ie+1) - hx1f(0,indcs.ie);
    // problem/he_grid_dump (default false): rank 0 writes the code's own radial grid of
    // block 0, "i x1f(i) x1v(i) x1f(i+1)-x1f(i)" for every cell incl. ghosts (i = 0 is
    // the first ghost; active cells is..ie), to <basename>.x1grid.txt
    if (pin->DoesParameterExist("problem","he_grid_dump") &&
        pin->GetBoolean("problem","he_grid_dump") &&
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
  hs_rlo_ = rlo_g - 2.0*dxlo;
  const Real rhi = rhi_g + 2.0*dxhi;
  hs_nf_ = nf;
  hs_dr_ = (rhi - hs_rlo_)/(nf - 1);
  std::vector<Real> hr(nf), hd(nf), he(nf), hF(nf), hE(nf, -1.0), hM(nf, 0.0);
  for (int n=0; n<nf; ++n) hr[n] = hs_rlo_ + n*hs_dr_;

  // problem/he_wind_ic (he-wind-bc, default false): above the radius r_c where the
  // column's density falls below that of a beta-law wind (he_wind_mdot [g/s, full
  // sphere], he_wind_vinf, he_wind_beta, he_wind_v0 [cm/s], he_wind_r0 [cm]: v = max(v0,
  // vinf (1 - r0/r)^beta), rho = Mdot/(4 pi r^2 v)), the IC is that wind: F_r = F_c (r_c/
  // r)^2, E = (F_r/c) (1 + (c E_c/F_c - 1)(r_c/r)^2) (free streaming far out), T_gas =
  // T_c (E/E_c)^(1/4) (continuous at r_c), v_r = the beta law.  The file then need not
  // cover the mesh top (it is extrapolated in log up to r_c, which must lie inside it).
  const bool wic = pin->GetOrAddBoolean("problem","he_wind_ic",false);
  std::vector<Real> hV(nf, 0.0), hdc, hec;
  Real wrc = 0.0;
  // ---- the column file: r rho eint F_r, ascending r, resampled in log
  {
    std::ifstream f(fn);
    if (!f.good()) HsFatal("cannot open problem/he_ic_file '" + fn + "'", __LINE__);
    std::vector<Real> fr, fd, fe, fF, fE, fM;
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      Real a, b, c, d, e = -1.0;
      if (!(ss >> a >> b >> c >> d)) continue;
      if (ncols == 5 && !(ss >> e)) continue;
      Real tcol = 0.0, fm = 0.0;
      if (hs_mlt_ && !(ss >> tcol >> fm)) {
        HsFatal("problem/mlt_flux_frozen: he_ic_file line without column 7 (fmlt)",
                __LINE__);
      }
      fM.push_back(fm);
      fr.push_back(a); fd.push_back(b); fe.push_back(c); fF.push_back(d);
      fE.push_back(e);
    }
    if (fr.size() < 2 || fr.front() > hr.front() || (!wic && fr.back() < hr.back())) {
      std::ostringstream os;
      os << "problem/he_ic_file '" << fn << "' must cover r = " << hr.front() << " .. "
         << hr.back() << " (the mesh plus the ghost margin), it has "
         << (fr.empty() ? 0.0 : fr.front()) << " .. " << (fr.empty() ? 0.0 : fr.back());
      HsFatal(os.str(), __LINE__);
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
      if (hs_mlt_) hM[n] = (1.0 - w)*fM[kk] + w*fM[kk+1];
    }
    if (wic) {
      const Real mdot = pin->GetReal("problem","he_wind_mdot");
      const Real vinf = pin->GetReal("problem","he_wind_vinf");
      const Real wbeta = pin->GetOrAddReal("problem","he_wind_beta",1.0);
      const Real wv0 = pin->GetOrAddReal("problem","he_wind_v0",1.0e6);
      const Real wr0 = pin->GetReal("problem","he_wind_r0");
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
        HsFatal("problem/he_wind_ic: the wind density never exceeds the column's inside "
                "the file's range above he_wind_r0", __LINE__);
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
        std::cout << "he_star_m1: he_wind_ic: Mdot = " << mdot << " g/s, vinf = " << vinf
                  << ", beta = " << wbeta << ", r0 = " << wr0 << "; the wind starts at "
                  << "r_c = " << wrc << " (rho " << hd[nc] << ", T [K] "
                  << tc*hs_tk << ", v "
                  << hV[nc] << ")" << std::endl;
      }
    }
  }

  // ---- the Rosseland + Planck tables (one shared grid), handed to the module
  Real tab_lt_lo = 0.0, tab_lt_hi = 0.0, tab_ld_lo = 0.0, tab_ld_hi = 0.0;
  {
    const std::string rt = pin->GetString("problem","he_opac_table");
    const std::string pt = pin->GetString("problem","he_planck_table");
    DvceArray2D<Real> krt, kpt;
    DvceArray1D<Real> mlT, mlD, plT, plD;
    int mnT = 0, mnD = 0, pnT = 0, pnD = 0;
    Real a1, a2, a3, a4, b1, b2, b3, b4;
    // problem/he_opac_logd_min (default none): extend both tables to this log10 rho
    const Real ldx = pin->DoesParameterExist("problem","he_opac_logd_min") ?
                     pin->GetReal("problem","he_opac_logd_min") : 1.0e300;
    // problem/he_opac_extend_hold (gd_physfix_1005, read only when named, default false):
    // the he_opac_logd_min extension holds both tables at their rho edge (slope 0)
    // instead of continuing log kappa in log rho; rho >= the table edge is unchanged
    const bool ldh = pin->DoesParameterExist("problem","he_opac_extend_hold") &&
                     pin->GetBoolean("problem","he_opac_extend_hold");
    HsReadOpacityTable(rt, krt, mlT, mlD, mnT, mnD, tab_lt_lo, tab_lt_hi, tab_ld_lo,
                       tab_ld_hi, ldx, ldh);
    HsReadOpacityTable(pt, kpt, plT, plD, pnT, pnD, a1, a2, a3, a4, ldx, ldh);
    b1 = tab_lt_lo; b2 = tab_lt_hi; b3 = tab_ld_lo; b4 = tab_ld_hi;
    if (mnT != pnT || mnD != pnD || fabs(a1 - b1) + fabs(a2 - b2) + fabs(a3 - b3) +
        fabs(a4 - b4) > 1.0e-9) {
      HsFatal("the Rosseland and Planck tables are not on the same grid", __LINE__);
    }
    pm1->SetOpacityTables(krt, kpt, mlT, mlD, mnT, mnD);
    // the lookup's units against the EOS's own code temperature and <units> density
    const Real tcgs = hs_tk;
    const Real dcgs = (pmbp->punit != nullptr) ? pmbp->punit->density_cgs() : 1.0;
    if (std::fabs(pm1->otab.tunit/tcgs - 1.0) > 1.0e-5 ||
        std::fabs(pm1->otab.dunit/dcgs - 1.0) > 1.0e-5) {
      std::ostringstream os;
      os << "<rad_m1>/temp_unit_kelvin = " << pm1->otab.tunit << ", rho_unit_cgs = "
         << pm1->otab.dunit << " but this run's units are (" << tcgs << ", " << dcgs
         << ")";
      HsFatal(os.str(), __LINE__);
    }
  }

  // ---- the column on the device: T, E = a T^4 and kappa_t
  Kokkos::realloc(hs_rho_, nf);
  Kokkos::realloc(hs_eint_, nf);
  DvceArray1D<Real> cE("hs_E", nf), cT("hs_T", nf), ckt("hs_kt", nf), cF("hs_F", nf);
  {
    auto h1 = Kokkos::create_mirror_view(hs_rho_);
    auto h2 = Kokkos::create_mirror_view(hs_eint_);
    auto h4 = Kokkos::create_mirror_view(cF);
    auto h5 = Kokkos::create_mirror_view(cE);
    for (int n=0; n<nf; ++n) {
      h1(n) = hd[n];
      h2(n) = he[n];
      h4(n) = hF[n];
      h5(n) = hE[n];
    }
    Kokkos::deep_copy(cE, h5);
    Kokkos::deep_copy(hs_rho_, h1);
    Kokkos::deep_copy(hs_eint_, h2);
    Kokkos::deep_copy(cF, h4);
  }
  {
    auto crho = hs_rho_, ceint = hs_eint_;
    radm1::M1OpacTab ot = pm1->otab;
    par_for("hs_col", DevExeSpace(), 0, nf-1, KOKKOS_LAMBDA(const int n) {
      const Real d = crho(n);
      const Real t = eos.Temperature(d, ceint(n));
      Real op, oe, of, os;
      radm1::M1TableOpacities(ot, d, t, op, oe, of, os);
      cT(n) = t;
      ckt(n) = of + os;
      if (!(cE(n) > 0.0)) cE(n) = ar*t*t*t*t;
    });
  }
  auto hT = Kokkos::create_mirror_view(cT);
  auto hkt = Kokkos::create_mirror_view(ckt);
  Kokkos::deep_copy(hT, cT);
  Kokkos::deep_copy(hkt, ckt);

  // ---- startup check (fatal): the column, ghost margin included, is inside the EOS
  // table (tabulated EOS only; an ideal gas has no range) and the opacity table grid,
  // in (rho [g/cm^3], T [K])
  {
    const Real tk = pm1->otab.tunit, dk = pm1->otab.dunit;
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
    if (bad) HsFatal(os.str(), __LINE__);
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: column inside the tables: log10 rho " << ldmin << " .. "
                << ldmax << ", log10 T " << ltmin << " .. " << ltmax << std::endl;
    }
  }
  // the imposed bottom flux must be the column's own F_r at r_in
  {
    const Real x = (hs_rin_ - hs_rlo_)/hs_dr_;
    const int i = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w = x - i;
    const Real fcol = (1.0 - w)*hF[i] + w*hF[i+1];
    if (fabs(hs_fin_/fcol - 1.0) > 1.0e-3) {
      std::ostringstream os;
      os << "<rad_m1>/implicit_flux_x1min = " << hs_fin_ << " but the column has F_r("
         << hs_rin_ << ") = " << fcol;
      HsFatal(os.str(), __LINE__);
    }
  }
  // tau from the top face, r_int (tau = he_tau_int), Phi_eff = Phi - int a_ref dr
  std::vector<Real> tau(nf, 0.0), aref(nf), cum(nf, 0.0);
  for (int n=0; n<nf; ++n) aref[n] = hkt(n)*hF[n]/cl;
  for (int n=nf-2; n>=0; --n) {
    const Real a0 = (hr[n] < rtop) ? hd[n]*hkt(n) : 0.0;
    const Real a1 = (hr[n+1] < rtop) ? hd[n+1]*hkt(n+1) : 0.0;
    tau[n] = tau[n+1] + 0.5*hs_dr_*(a0 + a1);
  }
  hs_rint_ = hs_rin_;
  for (int n=nf-1; n>=0; --n) {
    if (tau[n] >= tau_int) {
      hs_rint_ = std::max(hr[n], hs_rin_);
      break;
    }
  }
  for (int n=1; n<nf; ++n) cum[n] = cum[n-1] + 0.5*hs_dr_*(aref[n-1] + aref[n]);
  Real cin = 0.0;
  {
    const Real x = (hs_rin_ - hs_rlo_)/hs_dr_;
    int i = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w = x - i;
    cin = (1.0 - w)*cum[i] + w*cum[i+1];
  }
  // he_gm_column: G m(r) and Phi(r) of the column on the fine grid (trapezoid sums,
  // both zero-referenced at r_in like the point-mass Phi = GM (1/r_in - 1/r))
  std::vector<Real> hphig;
  if (hs_gmc_) {
    if (pmbp->punit == nullptr) HsFatal("problem/he_gm_column needs <units>", __LINE__);
    const Real lu = pmbp->punit->length_cgs(), tu = pmbp->punit->time_cgs();
    const Real gcode = Units::grav_constant_cgs*pmbp->punit->mass_cgs()*tu*tu/
                       (lu*lu*lu);
    std::vector<Real> cm(nf, 0.0), hgm(nf), cp(nf, 0.0);
    for (int n=1; n<nf; ++n) {
      cm[n] = cm[n-1] + 0.5*hs_dr_*4.0*M_PI*(hr[n-1]*hr[n-1]*hd[n-1] + hr[n]*hr[n]*hd[n]);
    }
    const Real x = (hs_rin_ - hs_rlo_)/hs_dr_;
    const int i0 = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w0 = x - i0;
    const Real cm0 = (1.0 - w0)*cm[i0] + w0*cm[i0+1];
    for (int n=0; n<nf; ++n) hgm[n] = hs_gm_ + gcode*(cm[n] - cm0);
    for (int n=1; n<nf; ++n) {
      cp[n] = cp[n-1] + 0.5*hs_dr_*(hgm[n-1]/(hr[n-1]*hr[n-1]) + hgm[n]/(hr[n]*hr[n]));
    }
    const Real cp0 = (1.0 - w0)*cp[i0] + w0*cp[i0+1];
    hphig.resize(nf);
    for (int n=0; n<nf; ++n) hphig[n] = cp[n] - cp0;
    Kokkos::realloc(hs_phig_, nf);
    Kokkos::realloc(hs_gmr_, nf);
    auto h1 = Kokkos::create_mirror_view(hs_phig_);
    auto h2 = Kokkos::create_mirror_view(hs_gmr_);
    for (int n=0; n<nf; ++n) {
      h1(n) = hphig[n];
      h2(n) = hgm[n];
    }
    Kokkos::deep_copy(hs_phig_, h1);
    Kokkos::deep_copy(hs_gmr_, h2);
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: he_gm_column: G m(r_in) = " << hs_gm_
                << ", G m(r_top) = " << hgm[nf-1] << " (G = " << gcode
                << " code), static monopole of the initial column" << std::endl;
    }
  }
  DvceArray1D<Real> cphi("hs_phieff", nf);
  {
    auto hp = Kokkos::create_mirror_view(cphi);
    for (int n=0; n<nf; ++n) {
      if (hs_gmc_) {
        hp(n) = hphig[n] - (cum[n] - cin);
      } else {
        hp(n) = hs_gm_*(1.0/hs_rin_ - 1.0/hr[n]) - (cum[n] - cin);
      }
    }
    Kokkos::deep_copy(cphi, hp);
  }

  // ---- potentials: the TRUE one (etotgrav, the WB default) ...
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  const Real gm = hs_gm_, rin = hs_rin_;
  if (hs_gmc_) {
    auto phicc = ph->phicc0;
    auto ph1 = ph->phi0->x1f, ph2 = ph->phi0->x2f, ph3 = ph->phi0->x3f;
    auto pg = hs_phig_;
    const Real rlo0 = hs_rlo_, dr0 = hs_dr_;
    par_for("hs_phi_col", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real pc = HsLinInterp(pg, rlo0, dr0, nf, x1v(m,i));
      phicc(m,k,j,i) = pc;
      ph1(m,k,j,i) = HsLinInterp(pg, rlo0, dr0, nf, x1f(m,i));
      if (i == n1m1) ph1(m,k,j,i+1) = HsLinInterp(pg, rlo0, dr0, nf, x1f(m,i+1));
      ph2(m,k,j,i) = pc;
      ph3(m,k,j,i) = pc;
      if (j == n2m1) ph2(m,k,j+1,i) = pc;
      if (k == n3m1) ph3(m,k+1,j,i) = pc;
    });
  } else {
    auto phicc = ph->phicc0;
    auto ph1 = ph->phi0->x1f, ph2 = ph->phi0->x2f, ph3 = ph->phi0->x3f;
    par_for("hs_phi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
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
  const Real rlo = hs_rlo_, dr = hs_dr_;
  {
    ph->EnableWBEffectivePotential();
    auto pwc = ph->phicc_wb, pwf = ph->phi_wb_x1f;
    par_for("hs_phieff", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pwc(m,k,j,i) = HsLinInterp(cphi, rlo, dr, nf, x1v(m,i));
      pwf(m,k,j,i) = HsLinInterp(cphi, rlo, dr, nf, x1f(m,i));
      if (i == n1m1) pwf(m,k,j,i+1) = HsLinInterp(cphi, rlo, dr, nf, x1f(m,i+1));
    });
  }
  // the reference acceleration, face form: kt(cell) (F(r_l) + F(r_r))/(2c), from the
  // INITIAL column (identical on a restart)
  {
    DvceArray4D<Real> aref_d("hs_aref", nmb1+1, n3m1+1, n2m1+1, n1m1+1);
    par_for("hs_aref", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real kt = HsLinInterp(ckt, rlo, dr, nf, x1v(m,i));
      const Real rl = x1f(m,i), rr = x1f(m,i+1);
      aref_d(m,k,j,i) = 0.5*kt*(HsLinInterp(cF, rlo, dr, nf, rl)
                                + HsLinInterp(cF, rlo, dr, nf, rr))/cl;
    });
    // HeStarGravity gives the reference work (split); set before SetForceReference,
    // which resolves force_reference_work = auto
    pm1->fref_wsplit_ok = true;
    pm1->SetForceReference(aref_d);
  }

  // ---- DISCRETE hydrostatic balance of the initial cells (problem/he_ic_balance,
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
  // The inner ghosts are the column scaled to the edge cell, as HeStarBC makes them;
  // their T enters the stencil of the edge cell, which is iterated.  The inner wall face
  // (wall_closed_ix1, mirror state) is then balanced by construction.  E and F_r keep the
  // column's values.  Rebuilt on a restart (the BCs read the balanced ghost profiles).
  hs_bal_ = pin->GetOrAddBoolean("problem","he_ic_balance",false);
  if (hs_bal_) {
    if (ph->wb_option != WBOption::polytropic) {
      HsFatal("problem/he_ic_balance = true needs <hydro>/wb_option = polytropic",
              __LINE__);
    }
    const int n1 = n1m1 + 1;
    const int is = indcs.is, ie = indcs.ie, k0 = indcs.ks, j0 = indcs.js;
    auto hx1v = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1v);
    auto hpc = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ph->phicc_wb);
    auto hpf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), ph->phi_wb_x1f);
    std::vector<Real> dc(n1), ec(n1), tc(n1), bd(n1), be(n1);
    auto logint = [&](const std::vector<Real> &a, const Real r) {
      Real x = (r - hs_rlo_)/hs_dr_;
      int i = static_cast<int>(floor(x));
      i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
      Real w = x - i;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return exp((1.0 - w)*log(a[i]) + w*log(a[i+1]));
    };
    for (int i=0; i<n1; ++i) {
      dc[i] = logint(hd, hx1v(0,i));
      ec[i] = logint(he, hx1v(0,i));
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
    // THE WALL CELL IN THE MODULE'S OWN DISCRETE FORM.  In the IC state <rad_m1> gives
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
      Real x = (r - hs_rlo_)/hs_dr_;
      int n = static_cast<int>(floor(x));
      n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
      Real w = x - n;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return (1.0 - w)*hF[n] + w*hF[n+1];
    };
    auto phig = [&](const Real r) {   // he_gm_column: the column's Phi, linear
      Real x = (r - hs_rlo_)/hs_dr_;
      int n = static_cast<int>(floor(x));
      n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
      Real w = x - n;
      w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
      return (1.0 - w)*hphig[n] + w*hphig[n+1];
    };
    auto seg = [&](const Real a_, const Real ra, const Real rb) {   // Phi_eff(rb) - (ra)
      if (hs_gmc_) return phig(rb) - phig(ra) - a_*(rb - ra);
      return hs_gm_*(1.0/ra - 1.0/rb) - a_*(rb - ra);
    };
    // he_wind_ic: the march starts below r_c (the wind cells keep their IC)
    int ianc = ie;
    if (wic) {
      while (ianc > is && !(hx1v(0,ianc) < wrc)) --ianc;
    }
    // problem/he_ic_balance_rmax [cm] (default none = the top ghost): the march starts
    // below this radius and the cells above keep the column state.  For a column that
    // ends in a constant density floor (make_ic_mlt_star.py --atm iso): the floor is not
    // hydrostatic, and anchored there the march would rescale the whole column
    if (pin->DoesParameterExist("problem","he_ic_balance_rmax")) {
      const Real rbm = pin->GetReal("problem","he_ic_balance_rmax");
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
    const int NW = std::min(4, ie - is);
    DvceArray1D<Real> kd_d("hs_kd", NW + 1), kt_d("hs_kt", NW + 1), kk_d("hs_kk", NW + 1);
    auto hkd = Kokkos::create_mirror_view(kd_d);
    auto hkt = Kokkos::create_mirror_view(kt_d);
    radm1::M1OpacTab ot1 = pm1->otab;
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
      par_for("hs_kap1", DevExeSpace(), 0, NW, KOKKOS_LAMBDA(const int n) {
        Real op, oe, of, os;
        radm1::M1TableOpacities(ot1, kd_d(n), kt_d(n), op, oe, of, os);
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
      auto har = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pm1->arad_ref);
      for (int m=0; m<=nmb1; ++m) {
        for (int k=0; k<=n3m1; ++k) {
          for (int j=0; j<=n2m1; ++j) {
            for (int n=0; n<NW; ++n) har(m,k,j,is+n) = aw[n];
          }
        }
      }
      Kokkos::deep_copy(ph->phicc_wb, hpc);
      Kokkos::deep_copy(ph->phi_wb_x1f, hpf);
      Kokkos::deep_copy(pm1->arad_ref, har);
    }
    Real dmax = 0.0;
    int imax = is;
    for (int i=is; i<=ie; ++i) {
      if (fabs(bd[i]/dc[i] - 1.0) > dmax) {
        dmax = fabs(bd[i]/dc[i] - 1.0);
        imax = i;
      }
    }
    if (wic && ianc < ie) {   // the anchor cell is a wind cell again
      bd[ianc+1] = logint(hd, hx1v(0,ianc+1));
      be[ianc+1] = logint(he, hx1v(0,ianc+1));
    }
    Kokkos::realloc(hs_bd_, n1);
    Kokkos::realloc(hs_be_, n1);
    auto h1 = Kokkos::create_mirror_view(hs_bd_);
    auto h2 = Kokkos::create_mirror_view(hs_be_);
    for (int i=0; i<n1; ++i) {
      h1(i) = bd[i];
      h2(i) = be[i];
    }
    Kokkos::deep_copy(hs_bd_, h1);
    Kokkos::deep_copy(hs_be_, h2);
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: he_ic_balance: max |rho/rho_col - 1| = " << dmax
                << " at r = " << hx1v(0,imax) << ", max face residual |PR/PL - 1| = "
                << rmax << ", wall cells a_ref (module form), last change " << wchg
                << std::endl;
    }
  }

  // ---- the frozen MLT flux on the x1 faces (problem/mlt_flux_frozen)
  if (hs_mlt_) {
    auto hx1f = Kokkos::create_mirror_view_and_copy(HostMemSpace(), x1f);
    const int nfc = n1m1 + 2;
    Kokkos::realloc(hs_fm_, nfc);
    auto hfm = Kokkos::create_mirror_view(hs_fm_);
    Real fmax = 0.0, rmx = 0.0;
    for (int i=0; i<nfc; ++i) {
      const Real r = hx1f(0,i);
      Real x = (r - hs_rlo_)/hs_dr_;
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
    Kokkos::deep_copy(hs_fm_, hfm);
    {
      const auto &ms = pmy_mesh_->mesh_size;
      const Real omg = hcs ? (4.0*M_PI)
                           : (cos(ms.x2min) - cos(ms.x2max))*(ms.x3max - ms.x3min);
      hs_lmx_ = 0.0;
      hs_r2o_.assign(nfc, 0.0);
      for (int i=0; i<nfc; ++i) {
        hs_r2o_[i] = SQR(hx1f(0,i))*omg;
        hs_lmx_ = std::max(hs_lmx_, hfm(i)*hs_r2o_[i]);
      }
    }
    // the adaptive closure's radial taper, from the IC's F_MLT/F on the faces
    hs_tap_.assign(nfc, 0.0);
    {
      const Real tp = pin->GetOrAddReal("problem","mlt_closure_taper",0.01);
      for (int i=0; i<nfc; ++i) {
        const Real r = hx1f(0,i);
        Real x = (r - hs_rlo_)/hs_dr_;
        int n = static_cast<int>(floor(x));
        n = (n < 0) ? 0 : ((n > nf - 2) ? (nf - 2) : n);
        Real w = x - n;
        w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
        const Real fm = (1.0 - w)*hM[n] + w*hM[n+1];
        Real s = (tp > 0.0) ? fm/tp : ((fm > 0.0) ? 1.0 : 0.0);
        s = (s < 0.0) ? 0.0 : ((s > 1.0) ? 1.0 : s);
        hs_tap_[i] = (hfm(i) > 0.0) ? s : 0.0;
      }
    }
    // the RADIATION takes the deposit (<rad_m1>/esrc, added to the implicit solve's old
    // vector with the stage weights): beta = Pg/P ~ 0.01, so the gas holds ~1 % of the
    // heat capacity, and an explicit gas-side deposit was -0.5 e_gas per step at cfl 0.3
    // (v2/runs/mF03: force 0.94 rho g).  The implicit exchange gives the gas its share.
    if (pm1->transport == radm1::M1_TRANSPORT_EXPLICIT) {
      HsFatal("problem/mlt_flux_frozen needs <rad_m1>/transport = implicit (esrc)",
              __LINE__);
    }
    Kokkos::realloc(pm1->esrc, nmb1+1, n3m1+1, n2m1+1, n1m1+1);
    Kokkos::deep_copy(pm1->esrc, 0.0);
    auto es = pm1->esrc;
    auto fmd = hs_fm_;
    auto &area1 = pmbp->pcoord->area.x1f;
    auto &volume = pmbp->pcoord->volume;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    par_for("hs_esrc", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      es(m,k,j,i) = -(area1(m,k,j,i+1)*fmd(i+1) - area1(m,k,j,i)*fmd(i))/volume(m,k,j,i);
    });
    pm1->esrc_on = true;
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: mlt_flux_frozen: max F_MLT/F on the faces = " << fmax
                << " at r = " << rmx << "; F_MLT at the bottom/top faces = "
                << hfm(indcs.is) << " / " << hfm(indcs.ie + 1) << std::endl;
    }
  }

  // TEST of <rad_m1>/esrc (default off): a constant source he_esrc_const (erg/cm^3/s) in
  // the cells with r < he_esrc_rmax, ON TOP of the frozen MLT deposit if that is on
  const Real esc = pin->GetOrAddReal("problem","he_esrc_const",0.0);
  hs_esc_ = esc;
  if (esc != 0.0) {
    const Real rmx = pin->GetOrAddReal("problem","he_esrc_rmax",hs_rin_);
    hs_esrmx_ = rmx;
    if (!pm1->esrc_on) {
      Kokkos::realloc(pm1->esrc, nmb1+1, n3m1+1, n2m1+1, n1m1+1);
      Kokkos::deep_copy(pm1->esrc, 0.0);
      pm1->esrc_on = true;
    }
    auto es = pm1->esrc;
    const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
    const int ks = indcs.ks, ke = indcs.ke;
    par_for("hs_esrc_t", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (x1v(m,i) < rmx) es(m,k,j,i) += esc;
    });
  }

  // problem/mlt_ramp_start (s, default -1 = no ramp) and mlt_ramp_time (s): the weight
  // w(t) on the frozen MLT deposit, w = 1 for t < start, 0.5 (1 + cos(pi (t - start)/
  // time)) during the ramp, 0 after.  A function of t only, re-applied on a restart (both
  // keys may be overridden on the restart's command line).  Unused: esrc bitwise.
  hs_mrs_ = pin->GetOrAddReal("problem","mlt_ramp_start",-1.0);
  hs_mrt_ = pin->GetOrAddReal("problem","mlt_ramp_time",0.0);
  hs_mw_ = 1.0;
  if (hs_mrs_ >= 0.0 && !hs_mlt_) {
    HsFatal("problem/mlt_ramp_start >= 0 needs problem/mlt_flux_frozen = true", __LINE__);
  }
  // problem/mlt_closure: frozen (default: the IC's F_MLT, times w(t)) or adaptive (see
  // hs_ad_ above; w(t) caps it as an overall factor)
  {
    const std::string mc = pin->GetOrAddString("problem","mlt_closure","frozen");
    if (mc != "frozen" && mc != "adaptive") {
      HsFatal("problem/mlt_closure must be frozen or adaptive", __LINE__);
    }
    hs_ad_ = (mc == "adaptive");
  }
  if (hs_ad_) {
    if (!hs_mlt_) {
      HsFatal("problem/mlt_closure = adaptive needs problem/mlt_flux_frozen = true (the "
              "IC's F_MLT is its start profile and its FeCZ mask)", __LINE__);
    }
    if (pm1->f0x1.extent_int(0) <= 0) {
      HsFatal("problem/mlt_closure = adaptive needs the implicit M1 face flux f0x1",
              __LINE__);
    }
    hs_ad_tau_ = pin->GetOrAddReal("problem","mlt_relax_time",1500.0);
    hs_ad_every_ = pin->GetOrAddInteger("problem","mlt_closure_every",1);
    hs_ad_dout_ = pin->GetOrAddReal("problem","mlt_closure_dout",0.0);
    if (!(hs_ad_tau_ > 0.0) || hs_ad_every_ < 1) {
      HsFatal("problem/mlt_relax_time must be > 0 and mlt_closure_every >= 1", __LINE__);
    }
    const int nfc = n1m1 + 2;
    Kokkos::realloc(hs_fsub_, nfc);
    Kokkos::deep_copy(hs_fsub_, hs_fm_);
    Kokkos::realloc(hs_shs_, nfc, 8);
    hs_ad_tl_ = pmy_mesh_->time;
    hs_ad_cyc_ = -1;
    if (restart) {
      const auto &blk = pgen_rststate;
      if (blk.size() < 2*sizeof(std::int32_t) + sizeof(Real)) {
        if (global_variable::my_rank == 0) {
          std::cout << "### he_star_m1: this restart file carries no adaptive MLT state; "
                    << "F_sub re-seeded with the IC's F_MLT" << std::endl;
        }
      } else {
        std::int32_t hdr[2];
        std::memcpy(&(hdr[0]), blk.data(), sizeof(hdr));
        if (hdr[0] != nfc ||
            blk.size() != sizeof(hdr) + sizeof(Real) + nfc*sizeof(Real)) {
          HsFatal("the restart file's adaptive MLT profile does not match this radial "
                  "grid", __LINE__);
        }
        std::vector<Real> tmp(nfc + 1);
        std::memcpy(tmp.data(), blk.data() + sizeof(hdr), (nfc + 1)*sizeof(Real));
        hs_ad_tl_ = tmp[0];
        auto hf = Kokkos::create_mirror_view(hs_fsub_);
        for (int i=0; i<nfc; ++i) hf(i) = tmp[i+1];
        Kokkos::deep_copy(hs_fsub_, hf);
      }
    }
    {
      auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), hs_fsub_);
      hs_lsub_ = 0.0;
      for (int i=0; i<nfc; ++i) hs_lsub_ = std::max(hs_lsub_, hf(i)*hs_r2o_[i]);
    }
    if (hs_ad_dout_ > 0.0) {
      snprintf(hs_ad_file_, sizeof(hs_ad_file_), "%s.fsub.txt",
               pin->GetString("job","basename").c_str());
      hs_ad_tout_ = hs_ad_dout_*floor(pmy_mesh_->time/hs_ad_dout_ + 1.0);
      if (!restart) hs_ad_tout_ = 0.0;
    }
    pgen_rst_write_func = HsAdaptiveRstWrite;
    // esrc from F_sub (bitwise the start-up esrc when F_sub = F_MLT and w = 1)
    HsApplyMltW(pmy_mesh_, HsMltW(pmy_mesh_->time));
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: mlt_closure = adaptive, relax time " << hs_ad_tau_
                << " s, every " << hs_ad_every_ << " cycle(s), taper "
                << pin->GetReal("problem","mlt_closure_taper") << ", update clock "
                << hs_ad_tl_ << ", max L_sub = " << hs_lsub_ << std::endl;
    }
  }
  if (hs_mlt_) {
    const Real w = HsMltW(pmy_mesh_->time);
    if (w != hs_mw_) HsApplyMltW(pmy_mesh_, w);
    if (global_variable::my_rank == 0) {
      std::cout << "he_star_m1: mlt ramp start = " << hs_mrs_ << ", time = " << hs_mrt_
                << ", w(t = " << pmy_mesh_->time << ") = " << hs_mw_
                << ", max L_MLT (unscaled) = " << hs_lmx_ << std::endl;
    }
  }

  // top sponge (default off) and hooks
  hs_zflux_ = pin->GetOrAddBoolean("problem","he_wall_zero_flux",true);
  // the inner wall face sees the exact mirror of the interior-side state
  // (hydro_fluxes.cpp, wall_closed_ix1, as deep_hot_jupiter_rt problem/wall_closed): zero
  // mass flux to round-off.  Measured: it does not change the first-cell force
  // (-2.1e-3 rho g at 640 zones, first order in dx, gate 1a), which is not the wall flux.
  // problem/he_bc_inner (he-wind-bc): wall (default: the closed wall above) | inflow:
  // the ghosts hold the (balanced) initial column at fixed rho and eint, UNSCALED, with
  // v_r copied from the first active cell (zero gradient), so mass may enter or leave;
  // no wall-face correction.  problem/he_bc_outer: noinflow (default) | outflow: the
  // ghosts continue the edge cell as a constant-velocity wind, rho ~ r^-2, the edge's
  // specific internal energy (the no-inflow face correction of he_wall_zero_flux stays)
  // | hse (bsg-arm2): the ghosts are the edge cell walked outward in hydrostatic balance
  // with the effective gravity of the code's own forces (see HeStarBC), T held (or
  // continued along the interior's dT/dPhi_eff with he_bc_hse_tgrad), v_r copied where
  // it points out and zero where it points in, v_t copied (box_convection bc_mode_top 4).
  {
    const std::string bi = pin->GetOrAddString("problem","he_bc_inner","wall");
    const std::string bo = pin->GetOrAddString("problem","he_bc_outer","noinflow");
    if ((bi != "wall" && bi != "inflow") ||
        (bo != "noinflow" && bo != "outflow" && bo != "hse")) {
      HsFatal("problem/he_bc_inner must be wall|inflow, he_bc_outer noinflow|outflow|hse",
              __LINE__);
    }
    hs_binf_ = (bi == "inflow");
    hs_bout_ = (bo == "outflow");
    hs_bhse_ = (bo == "hse");
    if (hs_bhse_) {
      hs_bhse_gmax_ = pin->GetOrAddReal("problem","he_bc_hse_gmax",0.9);
      hs_bhse_tg_ = pin->GetOrAddBoolean("problem","he_bc_hse_tgrad",false);
      std::string hf;
      if (restart && !pin->DoesParameterExist("problem","he_bc_hse_flux")) {
        hf = "cell";
        pin->SetString("problem","he_bc_hse_flux",hf);
        if (global_variable::my_rank == 0) {
          std::cout << "### he_star_m1: restart input has no problem/he_bc_hse_flux; "
                    << "keeping the old cell-flux Gamma (he_bc_hse_flux = cell)"
                    << std::endl;
        }
      } else {
        hf = pin->GetOrAddString("problem","he_bc_hse_flux","face");
      }
      if (hf != "face" && hf != "cell") {
        HsFatal("problem/he_bc_hse_flux must be face|cell", __LINE__);
      }
      hs_bhse_face_ = (hf == "face");
      if (hs_bhse_face_ && (pmbp->pradm1 == nullptr ||
                            pmbp->pradm1->f0x1.extent_int(0) == 0)) {
        HsFatal("problem/he_bc_hse_flux = face needs the implicit rad_m1 transport "
                "(face fluxes f0x1); use cell", __LINE__);
      }
      if (!(hs_bhse_gmax_ >= 0.0 && hs_bhse_gmax_ < 1.0)) {
        HsFatal("problem/he_bc_hse_gmax must be in [0,1)", __LINE__);
      }
    }
  }
  pmbp->phydro->wall_closed_ix1 = pin->GetOrAddBoolean("problem","he_wall_closed",
                                                       !hs_binf_);
  if (hs_binf_ && pmbp->phydro->wall_closed_ix1) {
    HsFatal("problem/he_bc_inner = inflow needs he_wall_closed = false", __LINE__);
  }
  hs_sp_rate_ = pin->GetOrAddReal("problem","he_sponge_rate",0.0);
  hs_sp_r0_ = pin->GetOrAddReal("problem","he_sponge_r0",hs_rint_);
  if (pin->DoesParameterExist("problem","he_sponge_mode")) {
    const std::string sm = pin->GetString("problem","he_sponge_mode");
    if (sm != "all" && sm != "radial") {
      HsFatal("problem/he_sponge_mode must be all|radial", __LINE__);
    }
    hs_sp_rad_ = (sm == "radial");
  }
  if (pin->DoesParameterExist("problem","he_sponge_dmax")) {
    hs_sp_dmax_ = pin->GetReal("problem","he_sponge_dmax");
  }
  hs_rtop_ = rtop;
  user_srcs_func = HeStarGravity;
  user_bcs_func = HeStarBC;
  user_hist_func = HeStarHist;
  pgen_final_func = HeStarFinal;
  if (global_variable::my_rank == 0) {
    const int i0 = static_cast<int>((hs_rin_ - hs_rlo_)/hs_dr_);
    std::cout << "he_star_m1: GM = " << hs_gm_ << ", F_in = " << hs_fin_
              << ", tau(r_in) = " << tau[i0] << ", r_int (tau = " << tau_int
              << ") = " << hs_rint_ << ", r_top = " << rtop
              << ", sponge rate = " << hs_sp_rate_ << " above r = " << hs_sp_r0_
              << std::endl;
  }
  if (restart) return;

  // ---- the initial state (active + ghost cells)
  const Real seed = pin->GetOrAddReal("problem","he_seed",0.0);
  const Real seedk = pin->GetOrAddReal("problem","he_seed_k",4.0);
  const Real srlo = pin->GetOrAddReal("problem","he_seed_rlo",hs_rin_);
  const Real srhi = pin->GetOrAddReal("problem","he_seed_rhi",hs_rint_);
  // multi-mode seed (he_seed_nk > 0, default 0 = the single mode above): nk random
  // lateral modes sin(2 pi (k2 (theta - theta_a)/L_theta + k3 (phi - phi_a)/L_phi) + ph),
  // integer k2, k3 drawn uniformly in [he_seed_kmin, he_seed_kmax] (periodic in the
  // wedge), amplitudes normalised to sum a_n^2 = 1 (field rms 1/sqrt 2), under the radial
  // envelope sin(pi (r - rlo)/(rhi - rlo)); box_convection's vpert_var = eint recipe.
  // he_seed_rad = true scales E by (1 + d)^4 with eint by (1 + d): a TEMPERATURE seed at
  // fixed rho, gas and radiation in equilibrium (a gas-only eint seed in a radiation-
  // dominated layer is removed by the stiff exchange in the first step).
  const int snk = pin->GetOrAddInteger("problem","he_seed_nk",0);
  const bool srad = pin->GetOrAddBoolean("problem","he_seed_rad",false);
  const int skmin = pin->GetOrAddInteger("problem","he_seed_kmin",1);
  const int skmax = pin->GetOrAddInteger("problem","he_seed_kmax",2);
  const int srng = pin->GetOrAddInteger("problem","he_seed_rng",1234);
  // he_seed_signs = positive (default, the historical draw: k2, k3 >= 0, so every
  // wavefront tilts the same diagonal way) | random (an independent random sign for k2
  // and for k3 of each mode).  he_seed_kdist = box (default: k2, k3 independently
  // uniform in [kmin, kmax]) | shell (|k| = sqrt(k2^2 + k3^2) uniform in [kmin, kmax],
  // direction uniform, both rounded to integers, k = 0 rejected; with signs = positive
  // the direction is limited to the first quadrant).  Both extra draws come AFTER the
  // historical per-mode draws (a separate pass on the same RNG stream), so the default
  // keys reproduce the old seed bitwise.
  const std::string ssign = pin->GetOrAddString("problem","he_seed_signs","positive");
  const std::string skdist = pin->GetOrAddString("problem","he_seed_kdist","box");
  if (ssign != "positive" && ssign != "random") {
    HsFatal("he_seed_signs must be positive or random", __LINE__);
  }
  if (skdist != "box" && skdist != "shell") {
    HsFatal("he_seed_kdist must be box or shell", __LINE__);
  }
  const bool srsign = (ssign == "random"), sshell = (skdist == "shell");
  if (snk > 0 && (skmin < 0 || skmax < skmin || !(srhi > srlo))) {
    HsFatal("need 0 <= he_seed_kmin <= he_seed_kmax and he_seed_rhi > he_seed_rlo",
            __LINE__);
  }
  if (snk > 0 && sshell && skmax < 1) {
    HsFatal("he_seed_kdist = shell needs he_seed_kmax >= 1", __LINE__);
  }
  DualArray2D<Real> smd("hs_seed_modes", std::max(snk, 1), 4);   // k2, k3, amp, phase
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
      std::cout << "he_star_m1: seed " << snk << " modes, amplitude " << seed
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
  auto ur = pm1->u0;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto crho = hs_rho_, ceint = hs_eint_;
  const bool gmc = hs_gmc_;
  auto pgc = hs_phig_;
  const bool bal = hs_bal_;
  auto cbd = hs_bd_, cbe = hs_be_;
  DvceArray1D<Real> cwv("hs_wv", wic ? nf : 1);
  if (wic) {
    auto hw = Kokkos::create_mirror_view(cwv);
    for (int n=0; n<nf; ++n) hw(n) = hV[n];
    Kokkos::deep_copy(cwv, hw);
  }
  par_for("hs_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const Real d = bal ? cbd(i) : HsLogInterp(crho, rlo, dr, nf, r);
    Real e = bal ? cbe(i) : HsLogInterp(ceint, rlo, dr, nf, r);
    // seed (default off): one smooth lateral mode in eint, wave number he_seed_k
    Real er = HsLogInterp(cE, rlo, dr, nf, r);
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
      if (snk > 0) {
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
      e *= fac;
      if (srad) er *= SQR(SQR(fac));
    }
    const Real vr = wic ? HsLinInterp(cwv, rlo, dr, nf, r) : 0.0;
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = d*vr;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    if (gmc) {
      uh(m,IEN,k,j,i) = e + d*HsLinInterp(pgc, rlo, dr, nf, r);
    } else {
      uh(m,IEN,k,j,i) = e + d*gm*(1.0/rin - 1.0/r);
    }
    if (wic) uh(m,IEN,k,j,i) += 0.5*d*vr*vr;
    ur(m,radm1::M1_E,k,j,i) = fmax(er, efl);
    ur(m,radm1::M1_F1,k,j,i) = HsLinInterp(cF, rlo, dr, nf, r);
    ur(m,radm1::M1_F2,k,j,i) = 0.0;
    ur(m,radm1::M1_F3,k,j,i) = 0.0;
  });
  // the comoving x1 face fluxes of the implicit transport (allocated only there)
  auto ff0 = pm1->f0x1;
  if (ff0.extent_int(0) >= nmb1 + 1) {
    const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
    const int e1 = ff0.extent_int(3) - 1;
    par_for("hs_icf", DevExeSpace(), 0, nmb1, 0, e3, 0, e2, 0, e1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      ff0(m,k,j,i) = HsLinInterp(cF, rlo, dr, nf, x1f(m,i));
    });
  }
  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn void HeStarGravity()
//! \brief the point-mass source with the x1 well-balanced pressure form (etotgrav), the
//! reference work of the WB kick (force_reference_work = split), the wall/top mass-flux
//! correction and the optional top sponge.

void HeStarGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto *ph = &fl;
  auto u0 = ph->u0;
  auto w0 = ph->w0;
  auto eos = ph->peos->eos_data;
  auto wbq0 = ph->wbq0;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  // the adaptive closure: once per cycle, at its first stage (before any M1 solve of
  // the cycle), every mlt_closure_every cycles; it re-applies w(t) itself
  // the closure is skipped once the ramp has ended (w(t) = 0 AND w = 0 already applied
  // to esrc): F_sub then enters nothing (esrc = w F_sub deposit, history slot 21 =
  // w L_sub), so F_sub and its clock stay FROZEN at their last values (also in the
  // restart files), except when a problem/mlt_closure_dout profile is due (diagnostic
  // cadence kept: the update runs and writes its lines)
  if (hs_ad_ && pm->ncycle != hs_ad_cyc_) {
    hs_ad_cyc_ = pm->ncycle;
    if (pm->ncycle % hs_ad_every_ == 0) {
      const bool off = (hs_mrs_ >= 0.0 && HsMltW(pm->time) == 0.0 && hs_mw_ == 0.0);
      const bool dout = (hs_ad_dout_ > 0.0 && pm->time >= hs_ad_tout_);
      if (!off || dout) HsAdaptiveUpdate(pm);
    }
  }
  // the frozen-MLT ramp: w at the step's start time (the same on every stage)
  if (hs_mlt_ && hs_mrs_ >= 0.0) {
    const Real w = HsMltW(pm->time);
    if (w != hs_mw_) HsApplyMltW(pm, w);
  }
  par_for("hs_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real d = w0(m,IDN,k,j,i);
    Real pl, pr, d1, d2, d3;
    WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
    const Real e_ = w0(m,IEN,k,j,i);
    const Real p = (e_ > 0.0) ? eos.Pressure(d, e_) : 0.5*(pl + pr);
    u0(m,IM1,k,j,i) += bdt*(area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))
                       /volume(m,k,j,i);
  });
  // <rad_m1>/force_reference_work = split: the work of the rho a_ref part of the WB kick
  // at the stage-start velocity; the radiation pays it in its solve
  auto *pm1 = pmbp->pradm1;
  if (pm1 != nullptr && pm1->fref_wsplit) {
    auto aref = pm1->arad_ref;
    if (pm1->FrefWaccOn()) {
      // fref-split-cons-1006 (time_scheme = be): the same increment also goes to the
      // accumulator the be solve makes E pay (acc = gam0 acc + increment)
      auto wacc = pm1->FrefWacc();
      const Real g0 = pm1->FrefWaccGam0(bdt, pm->dt);
      par_for("hs_grav_fws_x", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        const Real w = bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
        u0(m,IEN,k,j,i) += w;
        wacc(m,k,j,i) = g0*wacc(m,k,j,i) + w;
      });
    } else {
      par_for("hs_grav_fws", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
      KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
        u0(m,IEN,k,j,i) += bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
      });
    }
  }
  // x1 walls: the scaled-profile ghosts are not a mirror image of the edge cell, so the
  // Riemann flux through a "closed" wall carries mass.  The inner wall loses its mass,
  // transverse momentum and energy fluxes (u0 was just updated with them); the pressure
  // flux stays.  The outer face keeps OUTFLOW and loses only inflow (no mass enters).
  if (hs_zflux_) {
    auto flx1 = ph->uflx->x1f;
    auto &mbbcs = pmbp->pmb->mb_bcs;
    const bool binf = hs_binf_;
    par_for("hs_zflux", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (!binf && mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
        const Real f = bdt*area1(m,k,j,is)/volume(m,k,j,is);
        u0(m,IDN,k,j,is) -= f*flx1(m,IDN,k,j,is);
        u0(m,IM2,k,j,is) -= f*flx1(m,IM2,k,j,is);
        u0(m,IM3,k,j,is) -= f*flx1(m,IM3,k,j,is);
        u0(m,IEN,k,j,is) -= f*flx1(m,IEN,k,j,is);
      }
      if (mbbcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user &&
          !(flx1(m,IDN,k,j,ie+1) > 0.0)) {
        const Real f = bdt*area1(m,k,j,ie+1)/volume(m,k,j,ie);
        u0(m,IDN,k,j,ie) += f*flx1(m,IDN,k,j,ie+1);
        u0(m,IM2,k,j,ie) += f*flx1(m,IM2,k,j,ie+1);
        u0(m,IM3,k,j,ie) += f*flx1(m,IM3,k,j,ie+1);
        u0(m,IEN,k,j,ie) += f*flx1(m,IEN,k,j,ie+1);
      }
    });
  }
  if (hs_sp_rate_ <= 0.0) return;
  // top sponge: the velocity relaxes to zero at the rate he_sponge_rate*w,
  // w = ((r - r0)/(r_top - r0))^2 above r0; the kinetic energy removed leaves the total
  // energy (the heat is not kept, so that the sponge adds no buoyancy).
  // he_sponge_mode = radial: only v_r (the transverse momenta are untouched)
  auto &x1v = pmbp->pcoord->x1v;
  const Real rate = hs_sp_rate_, r0 = hs_sp_r0_, rt = hs_rtop_, spd = hs_sp_dmax_;
  if (hs_sp_rad_) {
    par_for("hs_sponge_r", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real r = x1v(m,i);
      const bool low = (u0(m,IDN,k,j,i) < spd);
      if (r <= r0 && !low) return;
      const Real w = low ? 1.0 : SQR((r - r0)/(rt - r0));
      const Real f = exp(-bdt*rate*w);
      const Real ke1 = 0.5*SQR(u0(m,IM1,k,j,i))/u0(m,IDN,k,j,i);
      u0(m,IM1,k,j,i) *= f;
      u0(m,IEN,k,j,i) -= (1.0 - f*f)*ke1;
    });
    return;
  }
  par_for("hs_sponge", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const bool low = (u0(m,IDN,k,j,i) < spd);
    if (r <= r0 && !low) return;
    const Real w = low ? 1.0 : SQR((r - r0)/(rt - r0));
    const Real f = exp(-bdt*rate*w);
    const Real d = u0(m,IDN,k,j,i);
    const Real ke0 = 0.5*(SQR(u0(m,IM1,k,j,i)) + SQR(u0(m,IM2,k,j,i))
                          + SQR(u0(m,IM3,k,j,i)))/d;
    u0(m,IM1,k,j,i) *= f;
    u0(m,IM2,k,j,i) *= f;
    u0(m,IM3,k,j,i) *= f;
    u0(m,IEN,k,j,i) -= (1.0 - f*f)*ke0;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void HeStarBC()
//! \brief x1 ghosts carrying the scaled initial profile (rho and eint scaled by the ratio
//! of the adjacent active cell to its own initial value), the radial velocity mirrored
//! (wall) or, at the outer edge where it points outward, copied (outflow).  M1 ghosts by
//! M1FillGhost (the implicit solve imposes its own face BCs).

void HeStarBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto &x1v = pmbp->pcoord->x1v;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto *ph = &fl;
  auto uh = ph->u0;
  auto phicc = ph->phicc0;
  auto crho = hs_rho_, ceint = hs_eint_;
  const Real rlo = hs_rlo_, dr = hs_dr_;
  const int nf = hs_nf_;
  const bool bal = hs_bal_;
  auto cbd = hs_bd_, cbe = hs_be_;
  const bool binf = hs_binf_, bout = hs_bout_, bhse = hs_bhse_, bhtg = hs_bhse_tg_;
  auto ur = pmbp->pradm1->u0;
  const Real cl = pmbp->pradm1->c_light;
  const Real efl = pmbp->pradm1->e_floor;
  // he_bc_outer = hse
  auto eos = ph->peos->eos_data;
  auto pwc = ph->phicc_wb;
  auto aref = pmbp->pradm1->arad_ref;
  auto opac = pmbp->pradm1->opac;
  const bool haveop = (opac.extent_int(0) > 0);
  const Real gm = hs_gm_, gmax = hs_bhse_gmax_;
  const bool gmc = hs_gmc_;
  auto cgm = hs_gmr_;
  const Real dfl = eos.dfloor;
  const bool bhfc = bhse && hs_bhse_face_;
  auto ff1 = pmbp->pradm1->f0x1;
  par_for("hs_bc", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    for (int side=0; side<2; ++side) {
      const bool lo = (side == 0);
      if (mbbcs.d_view(m, lo ? BoundaryFace::inner_x1 : BoundaryFace::outer_x1)
          != BoundaryFlag::user) continue;
      const int ia = lo ? is : ie;
      const Real da = uh(m,IDN,k,j,ia);
      const Real kea = 0.5*(SQR(uh(m,IM1,k,j,ia)) + SQR(uh(m,IM2,k,j,ia))
                            + SQR(uh(m,IM3,k,j,ia)))/da;
      const Real ea = uh(m,IEN,k,j,ia) - kea - da*phicc(m,k,j,ia);
      const bool inflo = lo && binf, outfl = !lo && bout;
      if (!lo && bhse) {
        // THE HYDROSTATIC TOP (he_bc_outer = hse).  Every ghost layer continues the LAST
        // ACTIVE cell a = ie (box_convection bc_mode_top 4, ported): with the WB pair's
        // effective potential Phi_eff (gravity minus the reference radiation force
        // a_ref of the column, what the x1 well-balanced stencil walks along) plus the
        // RESIDUAL radiation force the M1 module applies on top of it,
        //   a_rad = opac_T F_1/(c rho) of cell a  (the code's own force; the reference
        //           a_ref when the opacity is not yet filled, i.e. on the first step),
        //           F_1 per he_bc_hse_flux: face (default) = 0.5 (f0x1(ie) + f0x1(ie+1)),
        //           the comoving flux the transport carries and the gas force uses;
        //           cell = the clipped lab-frame cell F1 (old; noisy, ~0.76 of the face
        //           flux in the 3-D BSG top where |F| = cE clips)
        //   Gamma = a_rad/g(r_a) clamped to [0, he_bc_hse_gmax], scaled ~ r^-2,
        // the ghost is the cell walked by dPhi = Phi_eff(r_g) - Phi_eff(r_a)
        //   - (Gamma g_a - a_ref(a)) r_a^2 (1/r_a - 1/r_g)
        // with the WB closure (utils/wb_background.hpp WBAdvance): ISOTHERMAL (T_g = T_a,
        // exact for the ideal gas: rho_g = rho_a exp(-dPhi rho_a/p_a)) or, with
        // he_bc_hse_tgrad, polytropic along the interior's one-sided dT/dPhi_eff.  So
        // the total effective gravity over the ghosts is g (1 - Gamma), and with
        // a_rad = a_ref the ghost IS the WB background of the edge cell: the top-face
        // Riemann problem then carries no hydrostatic residual (no one-signed mass flux).
        // Not finite / not positive -> copy of the edge cell; floors after.
        const int ia1 = ie - 1;
        const Real ta = eos.Temperature(da, ea);
        const Real ra = x1v(m,ia);
        const Real ga = gmc ? HsLinInterp(cgm, rlo, dr, nf, ra)/SQR(ra) : gm/SQR(ra);
        const Real ar0 = aref(m,k,j,ia);
        Real arad = ar0;
        if (haveop && opac(m,radm1::M1_OP_T,k,j,ia) > 0.0) {
          // face (he_bc_hse_flux): the comoving radial flux the transport carries, the
          // mean of the cell's two x1 faces = the unclipped cell F1 minus its lab-frame
          // (v + v.D) E; cell: the clipped lab-frame cell F1 (old)
          const Real f1 = bhfc ? 0.5*(ff1(m,k,j,ia) + ff1(m,k,j,ia+1))
                               : ur(m,radm1::M1_F1,k,j,ia);
          arad = opac(m,radm1::M1_OP_T,k,j,ia)*f1/(cl*da);
        }
        Real gam = arad/ga;
        gam = (gam < 0.0) ? 0.0 : ((gam > gmax) ? gmax : gam);
        const Real dres = gam*ga - ar0;
        Real dltdphi = 0.0;
        if (bhtg) {
          const Real d1 = uh(m,IDN,k,j,ia1);
          const Real e1 = uh(m,IEN,k,j,ia1) - 0.5*(SQR(uh(m,IM1,k,j,ia1))
                          + SQR(uh(m,IM2,k,j,ia1)) + SQR(uh(m,IM3,k,j,ia1)))/d1
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
          radm1::M1FillGhost(ur, m, k, j, ig, k, j, ia, 1, 2, 1.0, cl, efl);
        }
        continue;
      }
      const Real sd = inflo ? 1.0 :
                      da/(bal ? cbd(ia) : HsLogInterp(crho, rlo, dr, nf, x1v(m,ia)));
      const Real se = inflo ? 1.0 :
                      ea/(bal ? cbe(ia) : HsLogInterp(ceint, rlo, dr, nf, x1v(m,ia)));
      for (int g=0; g<ng; ++g) {
        const int ig = lo ? (is - 1 - g) : (ie + 1 + g);
        const int im = lo ? (is + g) : (ie - g);        // the mirror cell
        const Real rg = x1v(m,ig);
        const Real dg = outfl ? da*SQR(x1v(m,ia)/rg) :
                        sd*(bal ? cbd(ig) : HsLogInterp(crho, rlo, dr, nf, rg));
        const Real eg = outfl ? ea*(dg/da) :
                        se*(bal ? cbe(ig) : HsLogInterp(ceint, rlo, dr, nf, rg));
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
        radm1::M1FillGhost(ur, m, k, j, ig, k, j, ia, 1, lo ? 0 : 2,
                           lo ? -1.0 : 1.0, cl, efl);
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void HeStarHist()
//! \brief history columns, all plain sums (see the file header).  Etot = gas total energy
//! (conserved IEN incl. kinetic and the static potential) + radiation energy: with closed
//! walls its rate is L_in - L_top(lab) (the comoving L_top differs from it by O(v/c)).
//! Mdot_top is the mass really leaving through the top (the inflow part is removed by
//! the gravity-source correction), Min_top the inflow that was removed.

void HeStarHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  // w_mlt, L_MLT (= w max_faces F_MLT r^2 Omega) only with mlt_flux_frozen
  pdata->nhist = hs_mlt_ ? 22 : 20;
  const char *lab[22] = {"L_bot", "L_mid", "L_int", "L_top", "L_in", "E_rad", "e_gas",
                         "M_int", "KE_int", "KEr_int", "Mr_int", "PV_int", "V_int",
                         "M_tot", "Mdot_top", "Mdot_bot", "Etot", "Min_top",
                         "v1sq_wall", "Picard", "w_mlt", "L_MLT"};
  for (int n=0; n<pdata->nhist; ++n) pdata->label[n] = lab[n];
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2, ks = indcs.ks, nx3 = indcs.nx3;
  const int imid = is + nx1/2;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto *ph = &fl;
  auto w0 = ph->w0;
  auto u0 = ph->u0;
  auto fl1 = ph->uflx->x1f;
  auto eos = ph->peos->eos_data;
  auto *pm1 = pmbp->pradm1;
  auto ur = pm1->u0;
  auto ff0 = pm1->f0x1;
  const bool haveff = (ff0.extent_int(0) > 0);
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  const Real rint = hs_rint_, fin = hs_fin_;
  // v1sq_wall: mean v_r^2 over the wall (first active) cells of the whole mesh
  const Real inwall = 1.0/(static_cast<Real>(pm->mesh_indcs.nx2)*pm->mesh_indcs.nx3);
  array_sum::GlobalSum sum_this;
  Kokkos::parallel_reduce("hs_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
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
      if (i == is) h.the_array[0] = ff0(m,k,j,i+1)*area1(m,k,j,i+1);
      if (i == imid) h.the_array[1] = ff0(m,k,j,i)*area1(m,k,j,i);
      if (x1f(m,i) <= rint && x1f(m,i+1) > rint) {
        h.the_array[2] = ff0(m,k,j,i)*area1(m,k,j,i);
      }
      if (i == ie) h.the_array[3] = ff0(m,k,j,i+1)*area1(m,k,j,i+1);
    }
    if (i == is) h.the_array[4] = fin*area1(m,k,j,i);
    h.the_array[5] = ur(m,radm1::M1_E,k,j,i)*vol;
    h.the_array[6] = w0(m,IEN,k,j,i)*vol;
    h.the_array[13] = w0(m,IDN,k,j,i)*vol;
    if (i == ie) {
      const Real fm = fl1(m,IDN,k,j,i+1)*area1(m,k,j,i+1);
      h.the_array[14] = (fm > 0.0) ? fm : 0.0;
      h.the_array[17] = (fm > 0.0) ? 0.0 : fm;
    }
    if (i == is) {
      h.the_array[15] = fl1(m,IDN,k,j,i)*area1(m,k,j,i);
      h.the_array[18] = SQR(w0(m,IVX,k,j,i))*inwall;
    }
    h.the_array[16] = (u0(m,IEN,k,j,i) + ur(m,radm1::M1_E,k,j,i))*vol;
    if (x1v(m,i) <= rint) {
      const Real d = w0(m,IDN,k,j,i);
      const Real v1 = w0(m,IVX,k,j,i), v2 = w0(m,IVY,k,j,i), v3 = w0(m,IVZ,k,j,i);
      h.the_array[7] = d*vol;
      h.the_array[8] = 0.5*d*(v1*v1 + v2*v2 + v3*v3)*vol;
      h.the_array[9] = 0.5*d*v1*v1*vol;
      h.the_array[10] = d*v1*vol;
      h.the_array[11] = eos.Pressure(d, w0(m,IEN,k,j,i))*vol;
      h.the_array[12] = vol;
    }
    msum += h;
  }, Kokkos::Sum<array_sum::GlobalSum>(sum_this));
  for (int n=0; n<20; ++n) pdata->hdata[n] = sum_this.the_array[n];
  if (hs_mlt_) {
    // host values, summed over ranks by the history output: rank 0 contributes
    const bool r0 = (global_variable::my_rank == 0);
    pdata->hdata[20] = r0 ? hs_mw_ : 0.0;
    pdata->hdata[21] = r0 ? hs_mw_*(hs_ad_ ? hs_lsub_ : hs_lmx_) : 0.0;
  }
  // Picard: mean passes per implicit solve since the previous history output (the
  // counters are MPI_MAX-reduced, identical on every rank: rank 0 contributes)
  pdata->hdata[19] = 0.0;
  if (global_variable::my_rank == 0) {
    const Real dn = pm1->impl_nstep - hs_pic_n0_, ds = pm1->impl_itsum - hs_pic_s0_;
    pdata->hdata[19] = (dn > 0.0) ? ds/dn : 0.0;
    hs_pic_n0_ = pm1->impl_nstep;
    hs_pic_s0_ = pm1->impl_itsum;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void HsApplyMltW()
//! \brief rebuild <rad_m1>/esrc = w (frozen MLT deposit) + the he_esrc_const test source,
//! with exactly the operations of the start-up kernels (w = 1: bitwise the same esrc).

void HsApplyMltW(Mesh *pm, const Real w) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto es = pmbp->pradm1->esrc;
  auto fmd = hs_ad_ ? hs_fsub_ : hs_fm_;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  auto &x1v = pmbp->pcoord->x1v;
  const Real esc = hs_esc_, rmx = hs_esrmx_;
  par_for("hs_esrc_w", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real e = -(area1(m,k,j,i+1)*fmd(i+1) - area1(m,k,j,i)*fmd(i))/volume(m,k,j,i);
    e = w*e;
    if (esc != 0.0 && x1v(m,i) < rmx) e += esc;
    es(m,k,j,i) = e;
  });
  hs_mw_ = w;
}

//----------------------------------------------------------------------------------------
//! \fn void HsAdaptiveUpdate()
//! \brief the adaptive shell-mean deficit closure (hs_ad_): shell sums (deterministic:
//! one thread per radial index, serial over the columns; MPI sum), the target, the
//! relaxation F_sub += (1 - exp(-(t - t_last)/tau)) (target - F_sub), clip >= 0, esrc.

void HsAdaptiveUpdate(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  radm1::FluidRef fl = radm1::FluidRef::Get(pmbp);
  auto w0 = fl.w0;
  auto eos = fl.peos->eos_data;
  auto *pm1 = pmbp->pradm1;
  auto ur = pm1->u0;
  auto f0 = pm1->f0x1;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  auto sh = hs_shs_;
  par_for("hs_ad_shell", DevExeSpace(), is, ie+1, KOKKOS_LAMBDA(const int i) {
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
            const Real v = w0(m,IVX,k,j,i), er = ur(m,radm1::M1_E,k,j,i);
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
  auto hsh = Kokkos::create_mirror_view_and_copy(HostMemSpace(), sh);
  const int nfc = hsh.extent_int(0);
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &hsh(is,0), 8*(ie + 2 - is), MPI_ATHENA_REAL, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  // luminosities of the wedge (the units of the history's L_in = F_in A(r_in))
  std::vector<Real> lc(nfc, 0.0), lrad(nfc, 0.0), lcf(nfc, 0.0), tgt(nfc, 0.0);
  for (int i=is; i<=ie; ++i) {
    const Real sv = hsh(i,2);
    const Real vm = hsh(i,3)/sv;
    const Real fc = (hsh(i,6)/sv - vm*hsh(i,4)/sv)
                    + (4.0/3.0)*(hsh(i,7)/sv - vm*hsh(i,5)/sv);
    lc[i] = fc*0.5*(hsh(i,1) + hsh(i+1,1));
  }
  const Real lin = hs_fin_*hsh(is,1);
  for (int i=is; i<=ie+1; ++i) {
    lrad[i] = hsh(i,0);
    lcf[i] = (i > is && i <= ie) ? 0.5*(lc[i-1] + lc[i]) : 0.0;
    const Real def = lin - lrad[i] - lcf[i];
    tgt[i] = hs_tap_[i]*((def > 0.0) ? def : 0.0)/hsh(i,1);
  }
  const Real t = pm->time;
  const Real a = 1.0 - exp(-(t - hs_ad_tl_)/hs_ad_tau_);
  hs_ad_tl_ = t;
  auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), hs_fsub_);
  hs_lsub_ = 0.0;
  for (int i=0; i<nfc; ++i) {
    Real f = hf(i) + a*(tgt[i] - hf(i));
    hf(i) = (f > 0.0) ? f : 0.0;
    hs_lsub_ = std::max(hs_lsub_, hf(i)*hs_r2o_[i]);
  }
  Kokkos::deep_copy(hs_fsub_, hf);
  HsApplyMltW(pm, HsMltW(t));
  // profiles (problem/mlt_closure_dout > 0): rank 0 appends, per output, three lines
  // "t <kind> values over the faces is..ie+1" in units of L_in: fsub (w A F_sub),
  // lrad, lconv; and one "t fmlt" line (A F_MLT) at the first output
  if (hs_ad_dout_ > 0.0 && t >= hs_ad_tout_ && global_variable::my_rank == 0) {
    auto hm = Kokkos::create_mirror_view_and_copy(HostMemSpace(), hs_fm_);
    FILE *fp = fopen(hs_ad_file_, "a");
    if (fp != nullptr) {
      const char *kind[4] = {"fsub", "lrad", "lconv", "fmlt"};
      for (int q=0; q<4; ++q) {
        if (q == 3 && hs_ad_tout_ > 0.0) break;
        fprintf(fp, "%.10e %s", t, kind[q]);
        for (int i=is; i<=ie+1; ++i) {
          Real v = 0.0;
          if (q == 0) v = hs_mw_*hf(i)*hsh(i,1);
          if (q == 1) v = lrad[i];
          if (q == 2) v = lcf[i];
          if (q == 3) v = hm(i)*hsh(i,1);
          fprintf(fp, " %.6e", v/lin);
        }
        fprintf(fp, "\n");
      }
      // STAGE CS3 (gate (c)): the shell-mean radiation energy density per CELL is..ie
      // (code units, not normalised; a trailing 0 keeps the face count)
      fprintf(fp, "%.10e emean", t);
      for (int i=is; i<=ie; ++i) {fprintf(fp, " %.6e", hsh(i,5)/hsh(i,2));}
      fprintf(fp, " 0\n");
      fclose(fp);
    }
    hs_ad_tout_ = hs_ad_dout_*floor(t/hs_ad_dout_ + 1.0);
  } else if (hs_ad_dout_ > 0.0 && t >= hs_ad_tout_) {
    hs_ad_tout_ = hs_ad_dout_*floor(t/hs_ad_dout_ + 1.0);
  }
}

//----------------------------------------------------------------------------------------
//! \fn std::vector<char> HsAdaptiveRstWrite()
//! \brief the adaptive closure's restart state: int32 nface, int32 0, Real t_last (the
//! update clock), nface Reals F_sub.

std::vector<char> HsAdaptiveRstWrite() {
  std::vector<char> out;
  if (!hs_ad_ || hs_fsub_.extent_int(0) <= 0) return out;
  const int nfc = hs_fsub_.extent_int(0);
  auto hf = Kokkos::create_mirror_view_and_copy(HostMemSpace(), hs_fsub_);
  const std::int32_t hdr[2] = {static_cast<std::int32_t>(nfc), 0};
  std::vector<Real> tmp(nfc + 1);
  tmp[0] = hs_ad_tl_;
  for (int i=0; i<nfc; ++i) tmp[i+1] = hf(i);
  out.resize(sizeof(hdr) + (nfc + 1)*sizeof(Real));
  std::memcpy(out.data(), &(hdr[0]), sizeof(hdr));
  std::memcpy(out.data() + sizeof(hdr), tmp.data(), (nfc + 1)*sizeof(Real));
  return out;
}

//----------------------------------------------------------------------------------------
//! \fn void HeStarFinal()
//! \brief release the namespace-scope Views before Kokkos::finalize.

void HeStarFinal(ParameterInput *pin, Mesh *pm) {
  (void) pin; (void) pm;
  hs_rho_ = DvceArray1D<Real>();
  hs_eint_ = DvceArray1D<Real>();
  hs_bd_ = DvceArray1D<Real>();
  hs_fm_ = DvceArray1D<Real>();
  hs_fsub_ = DvceArray1D<Real>();
  hs_shs_ = DvceArray2D<Real>();
  hs_be_ = DvceArray1D<Real>();
}
}  // namespace
