//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file red_giant.cpp
//! \brief A self-luminous stellar envelope: the outer shell of a red giant (or any cool
//! giant) between an inner wall that injects the stellar luminosity and an outer wall
//! above the photosphere.  Point-mass gravity, stellar Rosseland opacities, and the
//! radiative-conduction operator carrying the flux; the initial column is the
//! radiative-convective hydrostatic equilibrium of that same operator.
//!
//! Runs on the cubed sphere and on spherical polar (x1 = r), and as a Cartesian x1
//! column for one-dimensional tests (x1 is then r - rin + mesh x1min).
//!
//!   problem/mstar        mass inside the inner wall [g] (point mass; g = G M / r^2)
//!   problem/lstar        luminosity injected through the inner wall [erg/s]
//!   problem/teff         effective temperature [K]: the top of the column has
//!                        T^4 = Teff^4/2 (Eddington grey at tau = 0)
//!   problem/ptop         pressure at the outer wall [dyn/cm^2]
//!   problem/opac_table   log10 kappa_R on (log10 T, log10 rho): the file written by
//!                        tools/stellar_opac/merge_rosseland.py
//!   problem/mu           mean molecular weight for an IDEAL gas (ignored by a general
//!                        EOS, which carries its own composition)
//!   problem/vpert        velocity seed, in units of the local sound speed (0 = rest)
//!   problem/kappa_fac    scales kappa_R everywhere (the column here AND the conduction
//!                        operator, which is told).  A TEST knob: a giant envelope is
//!                        convective from tau ~ 10 down, so radiation carries only
//!                        grad_ad/grad_rad of L there and a 1-D column, which cannot
//!                        convect, cannot be steady.  kappa_fac ~ 1e-3 makes the
//!                        whole envelope radiative, and then the column must sit still.
//!   problem/kappa_const  > 0: a CONSTANT kappa_R [cm^2/g] replaces the table (column and
//!                        operator).  The grey Eddington atmosphere with constant
//!                        opacity has kappa p = g tau exactly, so grad_rad = tau/(4 tau +
//!                        8/3) < 0.25 < grad_ad: the whole ideal-gas column is RADIATIVE
//!                        and its equilibrium is known in closed form,
//!                        T^4 = (3/4) Teff^4 (tau + 2/3), p = g tau/kappa -- the 1-D
//!                        validation of the operator, the walls and the relaxation.
//!                        (kappa_fac cannot do this: scaling a T-dependent kappa leaves
//!                        kappa p / (g tau), and with it the convection zone, unchanged.)
//!   problem/mlt_alpha    > 0: a MIXING-LENGTH convective flux on the radial faces, with
//!                        l = alpha H_p (Boehm-Vitense; Kippenhahn et al. eq. 7.6),
//!                        wherever grad > grad_ad.  A giant envelope is convective from
//!                        just below the photosphere down, so a column that cannot
//!                        convect cannot be steady; this carries L - F_rad in 1-D, and
//!                        is a sub-grid option in 3-D (default 0 = off: resolve it).
//!   problem/rt_ck        true: the CORRELATED-K two-stream carries the optically thin
//!                        layers instead of the grey relaxation below.  The tables are
//!                        the Exo-FMS ones (problem/ck_table, ck_data_dir), the stellar
//!                        sweep is off (Teq = 0) and the internal temperature is teff,
//!                        so the band solver sees a self-luminous atmosphere.  The tau
//!                        blend hands over to radiative diffusion underneath exactly as
//!                        it does for the hot Jupiter.
//!   problem/inner_bc     wall (default) or open.  THE INNER WALL SITS INSIDE THE
//!                        CONVECTION ZONE -- for this star the radiative-convective
//!                        boundary is at 0.97 of the outer edge, so 94 % of the domain
//!                        convects and the wall is 12 density scale heights inside it.
//!                        A reflecting wall there bounces plumes, forces the convective
//!                        flux to zero and pins the envelope entropy at its initial
//!                        value.  `open` replaces it with the Stein-Nordlund / CO5BOLD
//!                        (Freytag+ 2012) treatment that solar_convection.cpp already
//!                        uses: the ghosts continue hydrostatically at the interior
//!                        temperature with the velocity copied, so plumes pass through,
//!                        and the lowest active layer is relaxed each step -- upflows
//!                        toward a prescribed deep adiabat (this is the energy input and
//!                        it sets the emergent flux), pressure toward the shell mean,
//!                        the mean density restored, and the net mass flux driven to
//!                        zero.  With `open` the luminosity is an OUTPUT, so
//!                        rad_flux_inner must be 0: the inflow entropy carries the
//!                        energy, not a diffusive flux through a wall.
//!   problem/s_relax_cs, problem/s_relax_cp   the two relaxation rates (0.1, 0.3)
//!   problem/column_dump  if set, rank 0 writes the initial column to this file
//!   problem/user_srcs    must be true (the gravity source lives here)
//!
//! <hydro|mhd>/isotropic_conduction = radiative with rad_kappa_src = table_rho makes the
//! conduction operator read the same opacity table; rad_flux_inner < 0 lets this file
//! set the inner flux to L/(4 pi rin^2).  Use ix1_bc = ox1_bc = user: both walls are
//! reflecting, with the ghost column continued hydrostatically from the initial state.
//!
//! THE INITIAL COLUMN.  From the outer wall inward (and outward through the ghosts):
//!     d ln p / dr = -rho g / p,     d ln T / d ln p = min(grad_rad, grad_ad),
//!     grad_rad = 3 kappa_R p L / (16 pi a c G M T^4),
//! with rho(p, T) and grad_ad from the run's own EOS and kappa_R from the table.  The
//! diffusion form of grad_rad is exactly the Eddington grey slope, so the only thing the
//! optically thin top adds is the boundary value T^4(tau = 0) = Teff^4/2.  Where the
//! radiative gradient exceeds the adiabatic one the column follows the adiabat -- the
//! Schwarzschild criterion, mixing length omitted, which is what a resolved convection
//! zone is supposed to supply.  The result is the steady state of the conduction operator
//! in the radiative layers and of adiabatic convection below, so the run starts from
//! its own equilibrium rather than relaxing to it.

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif
#include "globals.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "diffusion/conduction.hpp"
#include "utils/wb_background.hpp"
#include "units/units.hpp"
#include "pgen.hpp"
#include "pgen_eos_utils.hpp"
#include "utils/correlated_k.hpp"
#include "utils/two_stream_rt.hpp"

using pgen_eos::DensFromPT;
using pgen_eos::EintFromDensT;
using pgen_eos::GradAd;
using pgen_eos::TempKelvin;

void RedGiantGravity(Mesh *pm, Real bdt);
void RedGiantBC(Mesh *pm);
void RedGiantFinal(ParameterInput *pin, Mesh *pm);

namespace {
// physical constants, cgs
constexpr Real kGrav = 6.674e-8;
constexpr Real kBoltz = 1.380649e-16;
constexpr Real kMH = 1.6726e-24;
constexpr Real kArad = 7.5657e-15;
constexpr Real kClight = 2.99792458e10;
constexpr Real kSigmaSB = 5.670374419e-5;

// the star and the grid, code units unless said otherwise
Real gm_ = 0.0;          // G M, code units
Real rin_ = 1.0;         // radius of the inner wall
Real x1min_ = 0.0;       // mesh x1min (the Cartesian column maps x1 -> r)
bool curv_ = false;      // spherical polar or cubed sphere: x1 IS r
Real rgas_ = 1.0;        // k/(mu m_H) for an ideal gas, code units
Real gm1_ = 0.4;
bool etotgrav_ = false;

// the initial column on a fine uniform grid in r: ln p [code] and T [K]
DvceArray1D<Real> lnp_d_, tk_d_;
// The MIXING-LENGTH convective velocity on the same fine column, for the velocity seed.
// problem/vpert_mlt (default true when the column is built with MLT): the seed amplitude
// is v_c(r) rather than a Mach number.  A red giant's surface convection is transonic --
// v_c ~ (F/rho)^(1/3) is 1.4e6 cm/s against c_s = 5.4e5 at the photosphere -- so a Mach
// 1e-3 seed starts three orders of magnitude below the amplitude the flow has to reach,
// and the layers cool faster than that many e-foldings take.  problem/vpert then
// multiplies v_c instead of c_s.
DvceArray1D<Real> vc_d_;
bool vpert_mlt_ = true;
// problem/ic_tau_rad: the optical depth above which the INITIAL COLUMN follows grad_rad
// rather than the convective gradient.  It used to be <hydro>/rad_tau_hi, which made the
// initial star depend on where the RT/diffusion handover was put: moving the handover to
// tau = 300 for stability also built the column radiative down to tau = 300, where
// grad_rad/grad_ad ~ 1e3, i.e. a different and far hotter star.  Defaults to rad_tau_hi,
// so nothing changes unless it is set.
Real ic_tau_rad_ = -1.0;
// problem/vpert_mach_max (default 0.3): cap the seed at this fraction of the local sound
// speed.  v_c peaks just under the photosphere, where convection is least efficient and
// the density is lowest, and there it is TRANSONIC -- seeding a single smooth mode at
// that amplitude drove the top 50 radial cells non-finite within 3400 cycles.  Deep down
// v_c/c_s ~ 1e-4 and the cap never binds, so this only tames the surface.
Real vpert_mach_max_ = 0.3;
Real rlo_ = 0.0, drf_ = 1.0;
int nfine_ = 0;
// the opacity table (for the thin-layer relaxation) and the star's Teff
DvceArray2D<Real> ktab_;
DvceArray1D<Real> klT_, klD_;
int knT_ = 0, knD_ = 0;
Real teff_ = 0.0;
Real kfac_ = 1.0;        // problem/kappa_fac, applied to every table lookup
// The window in logR = log10 rho - 3 log10 T + 18 over which the opacity table is DATA.
// Both of its sources are tabulated in that variable over a finite range; outside it the
// merge fills with the edge value, which is constant in density and is not physics.  The
// table's header records the window and the initial column is checked against it below.
Real opac_lR_lo_ = -1.0e30, opac_lR_hi_ = 1.0e30;
// the open inner boundary (problem/inner_bc = open) and the deep adiabat it relaxes to
bool open_inner_ = false;
Real p_base_ = 0.0, t_base_ = 0.0;      // the initial column AT the inner wall, cgs
Real cs_change_ = 0.1, cp_change_ = 0.3;
// problem/open_budget (cycles, 0 = off): the open boundary is meant to be the star's
// ENERGY SOURCE, but each of its three passes rewrites the base shell's total energy and
// none of them is conservative by construction.  Accumulate what each one adds, in erg,
// and report it against L so a drain of the size of the luminosity cannot hide.
int open_budget_ = 0;
// problem/open_conserve (default true): make the density-restoring pass energy neutral.
// That pass adds the SAME drho to every cell of the base shell at the cell's own
// specific energy, so it moves mass between hot and cold columns and the total energy
// changes by drho * sum_i e_i V_i -- a sink of thousands of L once the shell has any
// entropy contrast (measured: -4e3 L, against the +1 L the boundary is meant to supply).
// With this on, whatever it adds is removed again uniformly over the shell.
bool open_conserve_ = true;
// problem/face_budget (cycles, 0 = off): the energy and mass the RADIAL FACES of the
// domain carry, integrated over each shell.  The open boundary's own passes are only
// half of its budget -- the other half is what the Riemann solver advects across the
// same face, and nothing else in the code reports it.  Sampled once per cycle, at the
// first stage, times dt: an O(dt) estimate of the cycle's transfer, which is all that
// is needed to compare a suspected sink against L.
// problem/open_nomassflux (default true with inner_bc = open): hold the NET mass flux
// through the inner FACE at zero.  The boundary's fourth pass drives the shell-mean
// radial velocity at the lowest cell CENTRE to zero, but the mass leaves across the
// face, and measured (problem/face_budget) that face drained 3.0e30 g and 2.0e44 erg --
// the whole of the run's energy loss -- in 2.3e4 s.  Each stage this removes the net
// mass the face just delivered, uniformly over the base shell, together with the energy
// and momentum that mass carried at the shell-mean specific enthalpy and velocity.  A
// zero-net-mass convective enthalpy flux, which is the luminosity the boundary is meant
// to supply, passes through untouched.
bool open_nomassflux_ = true;
// problem/wall_noflux (default true, inner_bc = wall): make the reflecting wall exactly
// impermeable.  The ghost carries the INITIAL COLUMN's density and pressure with only
// v_r mirrored -- which is what keeps it in balance with the well-balanced background --
// so once the interior drifts from the column the Riemann problem at the wall is no
// longer symmetric, HLLC's contact speed is not zero, and the face advects mass.
// Measured on the cubed sphere it reached -100 L and was still growing.  Each stage this
// removes the mass the face just delivered, cell by cell (a wall is impermeable
// locally, not on average), together with the energy that mass carried at the cell's own
// specific total enthalpy.  It does NOT touch the momentum flux, whose pressure term is
// what holds the star up, and it does NOT touch the conduction flux: the wall's imposed
// luminosity is added to the same flx1(IEN) and must survive.
bool wall_noflux_ = true;
// problem/sponge (default true), sponge_zbot (0.96, as a fraction of the radial domain),
// sponge_c (0.1, in units of cs/dr): an absorbing layer under the outer wall.  Both
// radial walls are reflecting, so the acoustic flux convection drives up through the
// photosphere has nowhere to go; a wave amplitude grows as rho^-1/2 and this envelope is
// far more stratified than the solar box that needed the same layer.  The default zbot
// sits ABOVE tau = 2/3 (0.952 of the domain for the standard input) so the layer never
// damps the convection zone itself.  The kinetic energy removed is DISCARDED, not
// thermalised -- the layer stands in for an atmosphere that carries the flux away, and
// depositing it locally would keep inflating the top, which is the problem being fixed.
bool sponge_on_ = true;
Real sponge_zbot_ = 0.96, sponge_c_ = 0.1;
int open_debug_ = 0;          // problem/open_debug: boundary calls to trace, 0 = off
int open_dbg_calls_ = 0;
int face_budget_ = 0;
int face_cycle_ = -1;
Real face_E_in_ = 0.0, face_E_out_ = 0.0, face_M_in_ = 0.0, face_M_out_ = 0.0;
Real face_E_in_l_ = 0.0, face_E_out_l_ = 0.0, face_t_last_ = 0.0;
Real open_dE_ent_ = 0.0, open_dE_prs_ = 0.0, open_dE_den_ = 0.0;
Real open_lstar_ = 0.0, open_t_last_ = 0.0;
Real open_dE_ent_l_ = 0.0, open_dE_prs_l_ = 0.0, open_dE_den_l_ = 0.0;
bool relax_ = false;     // the optically thin relaxation is on (radiative + tau blend)
bool rt_ck_ = false;     // problem/rt_ck: the band solver replaces that relaxation
// problem/rt_grey: the GREY two-stream on the same machinery, sharing the conduction
// module's opacity table.  Mutually exclusive with rt_ck; it replaces the same
// relaxation.  See two_stream_rt::rt_grey for why a grey control exists at all.
bool rt_grey_ = false;
// problem/ck_dump_t2, ck_dump_file2: re-arm the solver's one-shot column dump once the
// run reaches t2, so the emergent flux can be compared BEFORE and AFTER the atmosphere
// has adjusted to the band opacities.  The global energy budget cannot do this: the
// envelope holds ~1e48 erg while the flux imbalance integrates to ~1e40 over a short run.
Real ck_dump_t2_ = -1.0;
std::string ck_dump_file2_;
bool ck_dumped2_ = false;
Real mlt_alpha_ = 0.0;   // problem/mlt_alpha; > 0 turns the convective flux on
std::string mlt_dump_ = "";  // problem/mlt_dump: write the faces of one column once
bool mlt_dumped_ = false;
DvceArray2D<Real> fdiag_;   // (i, 8): grad, grad_ad, H_p, c_p, v, F, F_cap, F_used
DvceArray4D<Real> fconv_; // its radial face flux, code units, (m,k,j,i) on x1 faces

KOKKOS_INLINE_FUNCTION Real GravAt(const Real gm, const Real r) {
  return gm/(r*r);
}
// increasing outward, zero at the inner wall (the code's sign convention for phi)
KOKKOS_INLINE_FUNCTION Real PotAt(const Real gm, const Real rin, const Real r) {
  return gm*(1.0/rin - 1.0/r);
}
KOKKOS_INLINE_FUNCTION Real RadiusOf(const bool curv, const Real x1, const Real rin,
                                     const Real x1min) {
  return curv ? x1 : (rin + x1 - x1min);
}
// the column at radius r: ln p (code) and T (K), linear in r between fine nodes
KOKKOS_INLINE_FUNCTION void ColumnAt(const DvceArray1D<Real> &lnp,
                                     const DvceArray1D<Real> &tk, const int nfine,
                                     const Real rlo, const Real drf, const Real r,
                                     Real &lp, Real &t) {
  Real s = (r - rlo)/drf;
  int ii = static_cast<int>(s);
  ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
  const Real f = s - ii;
  lp = lnp(ii)*(1.0 - f) + lnp(ii+1)*f;
  t = tk(ii)*(1.0 - f) + tk(ii+1)*f;
}
// the same as ColumnAt, on host mirrors, for the start-up checks
template <typename V1>
void ColumnAtHost(const V1 &lnp, const V1 &tk, const int nfine, const Real rlo,
                  const Real drf, const Real r, Real &lp, Real &t) {
  Real s = (r - rlo)/drf;
  int ii = static_cast<int>(s);
  ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
  const Real f = s - ii;
  lp = lnp(ii)*(1.0 - f) + lnp(ii+1)*f;
  t = tk(ii)*(1.0 - f) + tk(ii+1)*f;
}

// log-bilinear kappa_R(T, rho) on the (log10 T, log10 rho) table, clamped
KOKKOS_INLINE_FUNCTION Real KappaTab(const DvceArray2D<Real> &tab,
                                     const DvceArray1D<Real> &lT,
                                     const DvceArray1D<Real> &lD, const int nT,
                                     const int nD, const Real tk, const Real rho) {
  const Real x = log10(tk), y = log10(rho);
  int i = 0, j = 0;
  Real fx = 0.0, fy = 0.0;
  if (x > lT(0)) {
    if (x >= lT(nT-1)) {
      i = nT-2; fx = 1.0;
    } else {
      i = static_cast<int>((x - lT(0))/(lT(1) - lT(0)));
      i = (i < 0) ? 0 : ((i > nT-2) ? nT-2 : i);
      fx = (x - lT(i))/(lT(i+1) - lT(i));
    }
  }
  if (y > lD(0)) {
    if (y >= lD(nD-1)) {
      j = nD-2; fy = 1.0;
    } else {
      j = static_cast<int>((y - lD(0))/(lD(1) - lD(0)));
      j = (j < 0) ? 0 : ((j > nD-2) ? nD-2 : j);
      fy = (y - lD(j))/(lD(j+1) - lD(j));
    }
  }
  const Real lk = (1.0-fx)*((1.0-fy)*tab(i,j) + fy*tab(i,j+1))
                +      fx *((1.0-fy)*tab(i+1,j) + fy*tab(i+1,j+1));
  return pow(10.0, lk);
}

//----------------------------------------------------------------------------------------
//! \fn ReadOpacityTable
//! \brief the merged Rosseland table: comment lines, one of them "# nT nD lTmin dlT
//! lDmin dlD", then nT*nD values of log10 kappa_R with T slowest

void ReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                      DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD) {
  std::ifstream f(fname);
  if (!f.good()) {
    std::cout << "### FATAL ERROR in red_giant: cannot open problem/opac_table '"
              << fname << "'" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string line;
  Real lt0 = 0.0, dlt = 0.0, ld0 = 0.0, dld = 0.0;
  bool have_grid = false;
  std::vector<Real> vals;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      // "# valid_logR lo hi": outside that window the table is edge-filled, not data
      const std::size_t vp = line.find("valid_logR");
      if (vp != std::string::npos) {
        std::istringstream vs(line.substr(vp + 10));
        Real lo, hi;
        if (vs >> lo >> hi) { opac_lR_lo_ = lo; opac_lR_hi_ = hi; }
      }
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
    std::cout << "### FATAL ERROR in red_giant: opacity table '" << fname
              << "' has no grid line or " << vals.size() << " values for "
              << nT << " x " << nD << std::endl;
    std::exit(EXIT_FAILURE);
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
  std::cout << "red_giant: opacity table '" << fname << "', " << nT << " x " << nD
            << " nodes, log10 T " << lt0 << ".." << lt0 + (nT-1)*dlt
            << ", log10 rho " << ld0 << ".." << ld0 + (nD-1)*dld << std::endl;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_srcs_func = RedGiantGravity;
  user_bcs_func = RedGiantBC;
  pgen_final_func = RedGiantFinal;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  const Real gamma = eos.gamma;
  const Real gm1 = gamma - 1.0, igm1 = 1.0/gm1;
  const bool etotgrav = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool wbdyn = is_mhd ? pmbp->pmhd->use_wellbalance_dynamic
                            : pmbp->phydro->use_wellbalance_dynamic;
  if (pmbp->punit == nullptr) {
    std::cout << "### FATAL ERROR in red_giant: a <units> block is required" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const Real lunit = pmbp->punit->length_cgs();
  const Real dunit = pmbp->punit->density_cgs();
  const Real punit_ = pmbp->punit->pressure_cgs();
  const Real vunit = pmbp->punit->velocity_cgs();
  const bool curv = pmy_mesh_->use_spherical_polar || pmy_mesh_->use_cubed_sphere;
  if (!curv && !(pmy_mesh_->one_d)) {
    std::cout << "### FATAL ERROR in red_giant: on a Cartesian mesh only a 1D x1 column "
              << "is supported (or use spherical polar / the cubed sphere)" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // --- the star
  const Real mstar = pin->GetReal("problem", "mstar");
  const Real lstar = pin->GetReal("problem", "lstar");
  const Real teff = pin->GetReal("problem", "teff");
  const Real ptop_cgs = pin->GetReal("problem", "ptop");
  const Real mu = pin->GetOrAddReal("problem", "mu", 0.62);
  const Real vpert = pin->GetOrAddReal("problem", "vpert", 0.0);
  vpert_mlt_ = pin->GetOrAddBoolean("problem", "vpert_mlt", true);
  ic_tau_rad_ = pin->GetOrAddReal("problem", "ic_tau_rad", -1.0);
  vpert_mach_max_ = pin->GetOrAddReal("problem", "vpert_mach_max", 0.3);
  const Real kfac = pin->GetOrAddReal("problem", "kappa_fac", 1.0);
  const Real kconst = pin->GetOrAddReal("problem", "kappa_const", 0.0);
  mlt_alpha_ = pin->GetOrAddReal("problem", "mlt_alpha", 0.0);
  // problem/mlt_alpha_ic: the mixing length used to build the INITIAL COLUMN, which is
  // a separate question from whether the MLT flux runs as a source term.  Defaults to
  // mlt_alpha, so setting one knob does both; set it alone (with mlt_alpha = 0) to start
  // a RESOLVED-convection run from a stratification that already carries L.
  const Real mlt_alpha_ic = pin->GetOrAddReal("problem", "mlt_alpha_ic", mlt_alpha_);
  mlt_dump_ = pin->GetOrAddString("problem", "mlt_dump", "");
  sponge_on_ = pin->GetOrAddBoolean("problem", "sponge", true);
  sponge_zbot_ = pin->GetOrAddReal("problem", "sponge_zbot", 0.96);
  sponge_c_ = pin->GetOrAddReal("problem", "sponge_c", 0.1);
  const std::string opac = pin->GetString("problem", "opac_table");
  const std::string dump = pin->GetOrAddString("problem", "column_dump", "");
  x1min_ = pmy_mesh_->mesh_size.x1min;
  const Real x1max = pmy_mesh_->mesh_size.x1max;
  rin_ = curv ? x1min_ : pin->GetReal("problem", "rin")/lunit;
  const Real rout = curv ? x1max : rin_ + (x1max - x1min_);
  gm_ = kGrav*mstar/(lunit*lunit*lunit)*(pmbp->punit->time_cgs()
                                         *pmbp->punit->time_cgs());
  // ideal gas only: p = rho Rgas T with T in KELVIN (pgen_eos_utils convention), so
  // Rgas carries velocity^2 per kelvin in code units
  rgas_ = kBoltz/(mu*kMH)/(vunit*vunit);
  curv_ = curv; gm1_ = gm1; etotgrav_ = etotgrav;
  const Real gm = gm_, rin = rin_, x1min = x1min_, rgas = rgas_;

  // --- the opacity table, for the column here and for the conduction operator
  DvceArray2D<Real> ktab;
  DvceArray1D<Real> klT, klD;
  int knT = 0, knD = 0;
  ReadOpacityTable(opac, ktab, klT, klD, knT, knD);
  if (kconst > 0.0) {
    auto h = Kokkos::create_mirror_view(ktab);
    for (int i=0; i<knT; ++i) {
      for (int j=0; j<knD; ++j) h(i,j) = log10(kconst);
    }
    Kokkos::deep_copy(ktab, h);
    std::cout << "red_giant: kappa_const = " << kconst << " cm^2/g REPLACES the table"
              << std::endl;
  }
  Conduction *pc = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
  if (pc != nullptr && pc->iso_cond_type.compare("radiative") == 0) {
    if (pc->rad_kappa_tab && !pc->rad_kappa_rho) {
      std::cout << "### FATAL ERROR in red_giant: use rad_kappa_src = table_rho, the "
                << "stellar table is on (T, rho)" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    pc->rad_kappa_fac = kfac;
    if (pc->rad_kappa_rho) {
      pc->rad_kr_tab = ktab;
      pc->rad_kr_lT = klT;
      pc->rad_kr_lP = klD;
      pc->rad_kr_nT = knT;
      pc->rad_kr_nP = knD;
    }
    if (pc->rad_flux_inner < 0.0) {
      const Real fin = lstar/(4.0*M_PI*SQR(rin*lunit));       // erg/cm^2/s
      pc->rad_flux_inner = fin/(punit_*vunit);
      std::cout << "red_giant: inner flux L/(4 pi rin^2) = " << fin
                << " erg/cm^2/s = sigma (" << pow(fin/kSigmaSB, 0.25) << " K)^4"
                << std::endl;
    }
    // THE OPTICALLY THIN LAYERS.  Diffusion is the wrong physics above tau ~ 1 and,
    // worse, explicitly stiff there: kappa_rad = 16 sigma T^3/(3 kappa_R rho) diverges
    // as kappa_R rho -> 0 while the layer's heat capacity vanishes, so the explicit
    // diffusion step collapses to microseconds even with the flux limiter (which does
    // nothing where the gradient is flat).  So the tau blend is REQUIRED: diffusion
    // carries the flux with weight w(tau) and, with weight 1 - w, each thin cell relaxes
    // toward the grey Eddington temperature of its own optical depth,
    //     T_eq^4 = (3/4) Teff^4 (tau + 2/3),
    // on its radiative time  t_rad = (rho c_v) / (4 kappa_P rho sigma T^3), applied as
    // the exact exponential decay so it is stable at any timestep.  In radiative
    // equilibrium that profile IS what the diffusion delivers below, so the two hand
    // over consistently for the luminosity Teff encodes.  See RedGiantGravity.
    if (!pc->rad_tau_mode) {
      std::cout << "### FATAL ERROR in red_giant: radiative conduction needs the tau "
                << "blend here (set rad_tau_lo/rad_tau_hi, e.g. 1 and 10): the optically "
                << "thin layers are relaxed toward the Eddington profile, not diffused"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    relax_ = true;
  }

  // --- the CORRELATED-K two-stream for the optically thin layers
  // (utils/two_stream_rt.hpp)
  rt_ck_ = pin->GetOrAddBoolean("problem", "rt_ck", false);
  rt_grey_ = pin->GetOrAddBoolean("problem", "rt_grey", false);
  if (rt_ck_ && rt_grey_) {
    std::cout << "### FATAL ERROR in red_giant: problem/rt_ck and problem/rt_grey are "
              << "two solvers for the same layers. Pick one." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (rt_grey_) {
    namespace ts = two_stream_rt;
    ts::rt_grey = true;
    ts::rt_split = true;             // implied: the grey kernel lives on the split path
    ts::rt_de_max = pin->GetOrAddReal("problem", "rt_de_max", 0.5);
    ts::rt_int_at_cut = false;       // the same handover argument as rt_ck below
    ts::rt_tint_override = teff;
    ts::rt_star_teff = 0.0;
    ts::rt_dump_file = pin->GetOrAddString("problem", "ck_dump_file", "");
    ts::rt_dump_m = pin->GetOrAddInteger("problem", "ck_dump_m", 0);
    ts::rt_dump_j = pin->GetOrAddInteger("problem", "ck_dump_j", -1);
    ts::rt_dump_k = pin->GetOrAddInteger("problem", "ck_dump_k", -1);
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    ck_dump_t2_ = pin->GetOrAddReal("problem", "ck_dump_t2", -1.0);
    ck_dump_file2_ = pin->GetOrAddString("problem", "ck_dump_file2", "");
    // ck_nquad picks the angular quadrature for the grey sweep too: 1 = hemispheric
    // mean, 2 = two-point Gauss
    correlated_k::ck_nq = pin->GetOrAddInteger("problem", "ck_nquad", 1);
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: GREY two-stream (problem/rt_grey), opacity from the "
                << "conduction module's table, " << pin->GetOrAddInteger("problem",
                   "ck_nquad", 1) << "-point angular quadrature" << std::endl;
    }
  }
  if (rt_ck_) {
    namespace ck = correlated_k;
    namespace ts = two_stream_rt;
    ts::rt_ck = true;
    ts::rt_split = true;             // implied by rt_ck
    ts::rt_ck_pcut = pin->GetOrAddReal("problem", "ck_pcut_bar", 1.0e30);
    ts::rt_de_max = pin->GetOrAddReal("problem", "rt_de_max", 0.5);
    // THE SOLVER INJECTS NO ENERGY, and cannot: its upward intensity at the cut is the
    // local band Planck function, which is the right handover to an optically thick,
    // diffusive interior.  The internal-flux term is off twice over -- set false here,
    // and forced false by the tau blend, which this problem generator requires above.
    // So the only energy input is the inner boundary: rad_flux_inner through a wall, or
    // the entropy of the inflow through an open one.
    ts::rt_int_at_cut = false;
    // T_int reaches only the grey picket-fence path and the albedo, neither of which runs
    // with rt_ck and the blend on, and the albedo multiplies a stellar flux that is zero.
    // It is set to a sane value rather than left at 0 so those helpers, which are
    // evaluated unconditionally, are not handed a degenerate temperature.
    ts::rt_tint_override = teff;
    ts::rt_dump_file = pin->GetOrAddString("problem", "ck_dump_file", "");
    ts::rt_dump_m = pin->GetOrAddInteger("problem", "ck_dump_m", 0);
    ts::rt_dump_j = pin->GetOrAddInteger("problem", "ck_dump_j", -1);
    ts::rt_dump_k = pin->GetOrAddInteger("problem", "ck_dump_k", -1);
    ts::rt_apply_debug = pin->GetOrAddInteger("problem", "rt_apply_debug", 0);
    ts::rt_apply_debug_n = pin->GetOrAddInteger("problem", "rt_apply_debug_n", 8);
    ck_dump_t2_ = pin->GetOrAddReal("problem", "ck_dump_t2", -1.0);
    ck_dump_file2_ = pin->GetOrAddString("problem", "ck_dump_file2", "");
    ts::rt_star_teff = 0.0;          // no host: the stellar band fractions are unused
    ts::rt_nchain = ck::CK_NB*ck::CK_NG*pin->GetOrAddInteger("problem", "ck_nquad", 1);
    ck::ck_nq = pin->GetOrAddInteger("problem", "ck_nquad", 1);
    ck::read_ck_table(pin->GetString("problem", "ck_table"), ts::rt_ck_pcut);
    ck::build_planck_fractions(ts::rt_ck_pcut);
    ck::read_ck_continuum(pin->GetString("problem", "ck_data_dir"),
                          pin->GetOrAddString("problem", "ck_swflux",
                                              "sw_band_flux_W121_11.txt"), 0.0);
    ck::ck_selftest();
    ck::ck_rt_selftest();
    // problem/opac_compare: write the Rosseland mean DERIVED FROM THE CK TABLE on its own
    // (T, p) grid, so it can be compared against the stellar table the diffusion uses.
    // The blend hands over between the two, and they are different data for the same gas:
    // if they disagree where they overlap, the handover is not conservative.  Built into
    // the conduction object and then undone, since the run must keep the stellar table.
    {
      const std::string ocmp = pin->GetOrAddString("problem", "opac_compare", "");
      if (!ocmp.empty() && pc != nullptr && pc->rad_kappa_tab) {
        auto sv_tab = pc->rad_kr_tab;
        auto sv_lT = pc->rad_kr_lT;
        auto sv_lP = pc->rad_kr_lP;
        const int sv_nT = pc->rad_kr_nT, sv_nP = pc->rad_kr_nP;
        const bool sv_rho = pc->rad_kappa_rho;
        pc->rad_kappa_rho = false;
        ck::ck_build_rosseland_table(pc);
        auto hk = Kokkos::create_mirror_view(pc->rad_kr_tab);
        auto hT = Kokkos::create_mirror_view(pc->rad_kr_lT);
        auto hP = Kokkos::create_mirror_view(pc->rad_kr_lP);
        Kokkos::deep_copy(hk, pc->rad_kr_tab);
        Kokkos::deep_copy(hT, pc->rad_kr_lT);
        Kokkos::deep_copy(hP, pc->rad_kr_lP);
        if (global_variable::my_rank == 0) {
          std::ofstream f(ocmp);
          f.precision(10);
          f << "# Rosseland mean from the CORRELATED-K table + continuum\n"
            << "# nT nP, then nT values of log10 T[K], nP of log10 p[dyn/cm^2], then\n"
            << "# nT*nP values of log10 kappa_R [cm^2/g], T slowest\n"
            << pc->rad_kr_nT << " " << pc->rad_kr_nP << "\n";
          for (int i = 0; i < pc->rad_kr_nT; ++i) f << hT(i) << "\n";
          for (int j = 0; j < pc->rad_kr_nP; ++j) f << hP(j) << "\n";
          for (int i = 0; i < pc->rad_kr_nT; ++i) {
            for (int j = 0; j < pc->rad_kr_nP; ++j) f << hk(i, j) << "\n";
          }
          std::cout << "red_giant: ck-derived Rosseland table written to '" << ocmp
                    << "'" << std::endl;
        }
        pc->rad_kr_tab = sv_tab;
        pc->rad_kr_lT = sv_lT;
        pc->rad_kr_lP = sv_lP;
        pc->rad_kr_nT = sv_nT;
        pc->rad_kr_nP = sv_nP;
        pc->rad_kappa_rho = sv_rho;
      }
    }
    // the star and grid numbers the solver reads; Teq = 0 switches the stellar sweep off
    hot_jupiter_param.Teq = 0.0;
    hot_jupiter_param.omega = 0.0;
    hot_jupiter_param.grav = kGrav*mstar/SQR(rin*lunit);
    hot_jupiter_param.ap = rin;
    hot_jupiter_param.Rgas = rgas_;
    hot_jupiter_param.met = pin->GetOrAddReal("problem", "met", 0.0);
    hot_jupiter_param.grav_point_mass = true;
    hot_jupiter_param.stellar_tide = false;
    hot_jupiter_param.rot_potential = false;
    std::cout << "red_giant: correlated-k two-stream ON, T_int = " << teff
              << " K, no irradiation" << std::endl;
  }
  ktab_ = ktab; klT_ = klT; klD_ = klD; knT_ = knT; knD_ = knD; teff_ = teff;
  kfac_ = kfac;
  if (mlt_alpha_ > 0.0) {
    Kokkos::realloc(fconv_, pmbp->nmb_thispack, n3m1+1, n2m1+1, n1m1+2);
    Kokkos::realloc(fdiag_, n1m1+2, 8);
  }

  // --- the initial column: fine grid in r from below the inner ghosts to above the
  // outer ones (10 % margins cover any stretch), integrated from the outer wall
  const int nfine = 40*pmy_mesh_->mesh_indcs.nx1 + 2*ng*40;
  const Real rlo = rin - 0.1*(rout - rin), rhi = rout + 0.1*(rout - rin);
  const Real drf = (rhi - rlo)/(nfine - 1);
  const int itop = static_cast<int>((rout - rlo)/drf + 0.5);
  DvceArray1D<Real> lnp("rg_lnp", nfine), tk("rg_tk", nfine);
  DvceArray1D<Real> kap("rg_kap", nfine), grad("rg_grad", nfine), tau("rg_tau", nfine);
  DvceArray1D<Real> vcc("rg_vc", nfine);
  {
    const Real ttop = teff*pow(0.5, 0.25);
    const Real ptop = ptop_cgs/punit_;
    const Real lum = lstar, mass = mstar, lun = lunit, dun = dunit, pun = punit_;
    const Real alpha_ic = mlt_alpha_ic;
    const Real vunit_c = vunit;
    const Real tun = pmbp->punit->temperature_cgs();
    // one serial sweep on the device (the EOS and the table live there); RK2 in r
    // where the thin-layer relaxation acts (tau < rad_tau_hi) the column is RADIATIVE
    // whatever the Schwarzschild criterion says: the relaxation targets the Eddington
    // profile there, and an adiabatic start would be cooled toward it from the first
    // step (in a real star convection is inefficient that far out anyway)
    const Real tauhi_blend = (pc != nullptr && pc->rad_tau_mode) ? pc->rad_tau_hi : 0.0;
    const Real tauhi = (ic_tau_rad_ >= 0.0) ? ic_tau_rad_ : tauhi_blend;
    par_for("rg_column", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(const int dummy) {
      // gradient at (p [code], T [K], tau): d ln T / d ln p and the state
      auto nabla = [&](const Real r, const Real p, const Real t, const Real ta,
                       Real &rho, Real &kr, Real &gr, Real &vc) {
        vc = 0.0;
        rho = DensFromPT(eos, rgas, p, t);
        kr = kfac*KappaTab(ktab, klT, klD, knT, knD, t, rho*dun);
        const Real grad_rad = 3.0*kr*(p*pun)*lum
                              /(16.0*M_PI*kArad*kClight*kGrav*mass*t*t*t*t);
        const Real grad_ad = GradAd(eos, gamma, rgas, p, t);
        if (ta < tauhi || grad_rad < grad_ad) {
          gr = grad_rad;                       // radiative, or the two-stream's domain
          return gr;
        }
        if (!(alpha_ic > 0.0)) {
          gr = grad_ad;             // perfectly efficient convection: no driving at all
          return gr;
        }
        // MIXING LENGTH.  grad_ad is the limit of infinitely efficient convection: it
        // carries the flux at zero superadiabaticity, so a column built on it is
        // marginally stable, has no convective flux, and cannot start convecting.  The
        // gradient that actually carries F = L/(4 pi r^2) splits it between radiation
        // and the same MLT closure the source term uses (see the rg_mlt_flux kernel),
        //   F (1 - grad/grad_rad) = rho cp T sqrt(g delta) l^2 x^(3/2)
        //                           / (4 sqrt2 Hp^(3/2)),   x = grad - grad_ad,
        // whose left side falls and right side rises in grad, so the root in
        // [grad_ad, grad_rad] is unique; bisect for it.  Deep down the convection is
        // efficient and it sits a hair above grad_ad, which is why this changes nothing
        // there and everything in the superadiabatic layer below the photosphere.
        const Real rcm = r*lun;
        const Real g = kGrav*mass/(rcm*rcm);                    // cgs
        const Real dcgs = rho*dun, pcgs = p*pun;
        const Real ftot = lum/(4.0*M_PI*rcm*rcm);               // erg/cm^2/s
        const Real hp = pcgs/(dcgs*g);
        const Real ell = alpha_ic*hp;
        const Real e = EintFromDensT(eos, rgas, igm1, rho, t);
        Real cv = (EintFromDensT(eos, rgas, igm1, rho, 1.01*t) - e)/(0.01*t)*pun/dcgs;
        if (!(cv > 0.0)) cv = 1.5*pcgs/(dcgs*t);
        Real delta = 1.0, cp = gamma*cv;
        if (eos.IsGeneral()) {
          const Real chit = eos.ChiT(rho, e), chir = eos.ChiRho(rho, e);
          delta = chit/chir;
          cp = cv*eos.Gamma1(rho, e)/chir;
        }
        const Real cmlt = dcgs*cp*t*sqrt(g*delta)*ell*ell/(4.0*sqrt(2.0)*hp*sqrt(hp));
        if (!(cmlt > 0.0) || !(ftot > 0.0)) {
          gr = grad_ad;
          return gr;
        }
        Real lo = grad_ad, hi = grad_rad;
        for (int it = 0; it < 60; ++it) {
          const Real mid = 0.5*(lo + hi);
          const Real x = mid - grad_ad;
          if (ftot*(1.0 - mid/grad_rad) > cmlt*x*sqrt(x)) {
            lo = mid;
          } else {
            hi = mid;
          }
        }
        gr = 0.5*(lo + hi);
        // the convective velocity that carries the rest of the flux, the same expression
        // the rg_mlt_flux source term uses.  This is what the seed is scaled to.
        vc = sqrt(g*delta*ell*ell*(gr - grad_ad)/(8.0*hp))/vunit_c;
        return gr;
      };
      // a hydrostatic RK2 step from (r, lnp, T, tau) over dr (signed; tau grows inward)
      auto step = [&](const Real r, const Real lp0, const Real t0, const Real dr,
                      Real &lp1, Real &t1, Real &ta) {
        Real rho, kr, gr, vcd;
        const Real p0 = exp(lp0);
        nabla(r, p0, t0, ta, rho, kr, gr, vcd);
        const Real dlnp_a = -rho*GravAt(gm, r)/p0;
        const Real lpm = lp0 + 0.5*dr*dlnp_a;
        const Real tm = t0*exp(0.5*dr*dlnp_a*gr);
        const Real rm = r + 0.5*dr;
        nabla(rm, exp(lpm), tm, ta, rho, kr, gr, vcd);
        const Real dlnp_m = -rho*GravAt(gm, rm)/exp(lpm);
        lp1 = lp0 + dr*dlnp_m;
        t1 = t0*exp(dr*dlnp_m*gr);
        ta -= kr*rho*dun*dr*lun;        // dr < 0 inward: tau increases
        if (ta < 0.0) ta = 0.0;
      };
      lnp(itop) = log(ptop);
      tk(itop) = ttop;
      Real ta = 0.0;
      for (int i = itop+1; i < nfine; ++i) {
        step(rlo + (i-1)*drf, lnp(i-1), tk(i-1), drf, lnp(i), tk(i), ta);
      }
      ta = 0.0;
      for (int i = itop-1; i >= 0; --i) {
        step(rlo + (i+1)*drf, lnp(i+1), tk(i+1), -drf, lnp(i), tk(i), ta);
      }
      // diagnostics: tau from the top down, kappa, and the gradient actually used
      Real tsum = 0.0;
      for (int i = nfine-1; i >= 0; --i) {
        Real rho, kr, gr, vcd;
        nabla(rlo + i*drf, exp(lnp(i)), tk(i), tsum, rho, kr, gr, vcd);
        kap(i) = kr; grad(i) = gr; vcc(i) = vcd;
        tau(i) = tsum;
        tsum += kr*rho*dun*drf*lun;
      }
      (void) tun;
    });
  }
  lnp_d_ = lnp; tk_d_ = tk; vc_d_ = vcc; rlo_ = rlo; drf_ = drf; nfine_ = nfine;

  // --- report the column, and dump it if asked
  {
    auto hlnp = Kokkos::create_mirror_view(lnp);
    auto htk = Kokkos::create_mirror_view(tk);
    auto hkap = Kokkos::create_mirror_view(kap);
    auto hgrad = Kokkos::create_mirror_view(grad);
    auto htau = Kokkos::create_mirror_view(tau);
    Kokkos::deep_copy(hlnp, lnp);
    Kokkos::deep_copy(htk, tk);
    Kokkos::deep_copy(hkap, kap);
    Kokkos::deep_copy(hgrad, grad);
    Kokkos::deep_copy(htau, tau);
    int iphot = nfine-1, iconv = -1;
    for (int i = nfine-1; i >= 0; --i) {
      if (htau(i) >= 2.0/3.0 && iphot == nfine-1) iphot = i;
    }
    // the outermost point (below the wall) where the column is on the adiabat, i.e.
    // where the radiative gradient EXCEEDS the one used (which is then grad_ad)
    for (int i = itop; i >= 0; --i) {
      Real p = exp(hlnp(i)), t = htk(i);
      Real grad_rad = 3.0*hkap(i)*(p*punit_)*lstar
                      /(16.0*M_PI*kArad*kClight*kGrav*mstar*t*t*t*t);
      if (grad_rad > hgrad(i)*(1.0 + 1.0e-6) && iconv < 0) {
        iconv = i;
      }
    }
    const int ibot = static_cast<int>((rin - rlo)/drf + 0.5);
    // THE OPACITY TABLE MUST BE DATA WHERE THIS STAR LIVES.  Outside the recorded logR
    // window the merge filled with the edge value, constant in density; a run there is
    // not using an opacity at all.  Checked on the ACTUAL radial grid including the
    // ghosts, since the boundary cells feed the wall flux -- not merely on the domain.
    // A first line of defence, not a guarantee: cells can still wander out as the
    // atmosphere relaxes.
    if (opac_lR_lo_ > -1.0e29) {
      auto &x1v_chk = pmbp->pcoord->x1v;
      auto hx1v = Kokkos::create_mirror_view(x1v_chk);
      Kokkos::deep_copy(hx1v, x1v_chk);
      int nbad = 0;
      Real wl = 0.0, wh = 0.0, wr = 0.0;
      for (int m = 0; m <= nmb1; ++m) {
        for (int i = 0; i <= n1m1; ++i) {
          const Real rr_ = curv ? hx1v(m,i) : rin;
          if (!(rr_ > 0.0)) continue;
          Real lp, tt;
          ColumnAtHost(hlnp, htk, nfine, rlo, drf, rr_, lp, tt);
          const Real dd = DensFromPT(eos, rgas, exp(lp), tt)*dunit;
          if (!(dd > 0.0)) continue;
          const Real lR = log10(dd) - 3.0*log10(tt) + 18.0;
          if (lR < opac_lR_lo_ || lR > opac_lR_hi_) {
            if (nbad == 0) { wl = lR; wh = lR; wr = rr_*lunit; }
            wl = fmin(wl, lR);
            wh = fmax(wh, lR);
            ++nbad;
          }
        }
      }
      if (nbad > 0) {
        std::cout << "### FATAL ERROR in red_giant: " << nbad << " cells (ghosts "
                  << "included) lie OUTSIDE the opacity table's valid window logR = ["
                  << opac_lR_lo_ << ", " << opac_lR_hi_ << "]: they reach logR " << wl
                  << " .. " << wh << ", first at r = " << wr << " cm. There the table is "
                  << "edge-filled, not data. Move the domain, or extend the tables."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      std::cout << "red_giant: every cell (ghosts included) is inside the opacity "
                << "table's valid window logR = [" << opac_lR_lo_ << ", "
                << opac_lR_hi_ << "]" << std::endl;
    }

    // --- the inner boundary, and the deep adiabat an OPEN one relaxes toward
    {
      const std::string ibc = pin->GetOrAddString("problem", "inner_bc", "wall");
      if (ibc.compare("open") == 0) {
        open_inner_ = true;
      } else if (ibc.compare("wall") != 0) {
        std::cout << "### FATAL ERROR in red_giant: problem/inner_bc must be wall or open"
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
      cs_change_ = pin->GetOrAddReal("problem", "s_relax_cs", 0.1);
      cp_change_ = pin->GetOrAddReal("problem", "s_relax_cp", 0.3);
      open_budget_ = pin->GetOrAddInteger("problem", "open_budget", 0);
      face_budget_ = pin->GetOrAddInteger("problem", "face_budget", 0);
      open_debug_ = pin->GetOrAddInteger("problem", "open_debug", 0);
      wall_noflux_ = pin->GetOrAddBoolean("problem", "wall_noflux", true);
      open_nomassflux_ = pin->GetOrAddBoolean("problem", "open_nomassflux", true);
      open_conserve_ = pin->GetOrAddBoolean("problem", "open_conserve", true);
      open_lstar_ = lstar;
      p_base_ = exp(hlnp(ibot))*punit_;
      t_base_ = htk(ibot);
      if (open_inner_) {
        if (pc != nullptr && pc->rad_flux_inner != 0.0) {
          std::cout << "### FATAL ERROR in red_giant: inner_bc = open carries the energy "
                    << "in as the ENTROPY of the inflow, so <block>/rad_flux_inner must "
                    << "be 0. It is " << pc->rad_flux_inner << ". With both, the "
                    << "luminosity is double counted." << std::endl;
          std::exit(EXIT_FAILURE);
        }
        std::cout << "red_giant: OPEN inner boundary (Stein-Nordlund / CO5BOLD). The "
                  << "deep adiabat is anchored at p = " << p_base_ << " dyn/cm^2, T = "
                  << t_base_ << " K; the emergent luminosity is now an OUTPUT."
                  << std::endl;
      }
    }
    // the top state must have a gas solution: under a general EOS with radiation the
    // total pressure cannot fall below a T^4/3, and SolveDensity then hands back the
    // table floor -- a column of floor density that looks like a run and is not one
    {
      const Real p_t = exp(hlnp(itop)), t_t = htk(itop);
      const Real rho_t = DensFromPT(eos, rgas, p_t, t_t)*dunit;
      const Real prad = kArad*t_t*t_t*t_t*t_t/3.0;
      const bool eos_rad = eos.IsGeneral() && eos.tbl.radiation;
      if (!(rho_t > 1.0e-30) || !std::isfinite(rho_t)
          || (eos_rad && prad > 0.9*ptop_cgs)) {
        std::cout << "### FATAL ERROR in red_giant: no gas solution at the outer wall: "
                  << "problem/ptop = " << ptop_cgs << " dyn/cm^2 against a radiation "
                  << "pressure a T^4/3 = " << prad << " at T = " << t_t << " K (rho = "
                  << rho_t << " g/cm^3, eos_radiation "
                  << (eos_rad ? "on" : "off")
                  << "). Raise ptop above a T^4/3, or switch eos_radiation off."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    if (global_variable::my_rank == 0) {
      std::cout << "red_giant: M = " << mstar << " g, L = " << lstar << " erg/s, Teff = "
                << teff << " K, rin = " << rin*lunit << " cm, rout = " << rout*lunit
                << " cm" << std::endl
                << "           wall p/T: outer " << ptop_cgs << " / " << htk(itop)
                << " K, inner " << exp(hlnp(ibot))*punit_ << " / " << htk(ibot) << " K"
                << std::endl
                << "           tau = 2/3 at r = " << (rlo + iphot*drf)*lunit
                << " cm, T = " << htk(iphot)
                << " K; radiative-convective boundary at r = "
                << ((iconv >= 0) ? (rlo + iconv*drf)*lunit : -1.0) << " cm"
                << std::endl;
      if (!dump.empty()) {
        std::ofstream df(dump);
        df.precision(10);
        df << std::scientific;
        auto hvc = Kokkos::create_mirror_view(vcc);
        Kokkos::deep_copy(hvc, vcc);
        df << "# red_giant initial column\n# r[cm] p[dyn/cm2] T[K] rho[g/cm3] kappa_R "
           << "grad tau v_mlt[cm/s]\n";
        for (int i = 0; i < nfine; ++i) {
          Real p = exp(hlnp(i)), t = htk(i);
          df << (rlo + i*drf)*lunit << " " << p*punit_ << " " << t << " "
             << DensFromPT(eos, rgas, p, t)*dunit << " " << hkap(i) << " " << hgrad(i)
             << " " << htau(i) << " " << hvc(i)*vunit << "\n";
        }
      }
    }
  }

  // --- the potential, needed on restarts too (rebuilt here, not stored)
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = is_mhd ? pmbp->pmhd->phi0.x1f : pmbp->phydro->phi0.x1f;
  DvceArray4D<Real> ph2 = is_mhd ? pmbp->pmhd->phi0.x2f : pmbp->phydro->phi0.x2f;
  DvceArray4D<Real> ph3 = is_mhd ? pmbp->pmhd->phi0.x3f : pmbp->phydro->phi0.x3f;
  const bool have_phi = (etotgrav || wbdyn);
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &x1f_ = pmbp->pcoord->xx1f;
  if (have_phi) {
    par_for("rg_phi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
      const Real xl = curv ? x1f_(m,i) : LeftEdgeX(i-is, indcs.nx1, x1lo, x1hi);
      const Real xr = curv ? x1f_(m,i+1) : LeftEdgeX(i+1-is, indcs.nx1, x1lo, x1hi);
      const Real phi_c = PotAt(gm, rin, RadiusOf(curv, xc, rin, x1min));
      phicc(m,k,j,i) = phi_c;
      ph1(m,k,j,i) = PotAt(gm, rin, RadiusOf(curv, xl, rin, x1min));
      if (i == n1m1) ph1(m,k,j,i+1) = PotAt(gm, rin, RadiusOf(curv, xr, rin, x1min));
      ph2(m,k,j,i) = phi_c;
      ph3(m,k,j,i) = phi_c;
      if (j == n2m1) ph2(m,k,j+1,i) = phi_c;
      if (k == n3m1) ph3(m,k+1,j,i) = phi_c;
    });
  }
  if (restart) return;

  // --- the initial state
  const int nf = nfine;
  auto vcd = vc_d_;
  const bool vpert_mlt = vpert_mlt_ && (mlt_alpha_ic > 0.0);
  const Real vmachmax = vpert_mach_max_;
  par_for("rg_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real x2lo = size.d_view(m).x2min, x2hi = size.d_view(m).x2max;
    const Real x3lo = size.d_view(m).x3min, x3hi = size.d_view(m).x3max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    Real lp, t;
    ColumnAt(lnp, tk, nf, rlo, drf, r, lp, t);
    const Real p = exp(lp);
    const Real d = DensFromPT(eos, rgas, p, t);
    const Real e = EintFromDensT(eos, rgas, igm1, d, t);
    Real v1 = 0.0, v2 = 0.0, v3 = 0.0;
    if (vpert > 0.0) {
      const Real x2v = CellCenterX(j-js, indcs.nx2, x2lo, x2hi);
      const Real x3v = CellCenterX(k-ks, indcs.nx3, x3lo, x3hi);
      // the amplitude: the MLT convective velocity where the column has one, else the
      // sound speed (the historical Mach-number seed)
      const Real csl = sqrt(gamma*p/d);
      // With the MLT seed, a cell the column gives no convective velocity -- the
      // radiative interior below the RCB, and the optically thin top -- gets NO seed.
      // Falling back to the sound speed there means seeding the densest material in the
      // star at Mach vpert, which with vpert = 1 is a transonic slug in the radiative
      // core: 8 orders of magnitude more kinetic energy than the rest of the seed put
      // together, and the upper domain went non-finite within 3400 cycles.
      Real amp = vpert_mlt ? 0.0 : csl;
      if (vpert_mlt) {
        const Real xf = (r - rlo)/drf;
        int i0 = static_cast<int>(xf);
        i0 = (i0 < 0) ? 0 : ((i0 > nf-2) ? nf-2 : i0);
        const Real f = xf - static_cast<Real>(i0);
        const Real vcl = vcd(i0)*(1.0 - f) + vcd(i0+1)*f;
        if (vcl > 0.0) amp = fmin(vcl, vmachmax*csl);
      }
      const Real cs = amp;
      const Real tp = 2.0*M_PI;
      v1 = vpert*cs*sin(3.0*tp*(xc - x1lo)/(x1hi - x1lo))
           *cos(2.0*tp*(x2v - x2lo)/(x2hi - x2lo))*cos(tp*(x3v - x3lo)/(x3hi - x3lo));
      v2 = vpert*cs*sin(2.0*tp*(x2v - x2lo)/(x2hi - x2lo))
           *cos(tp*(x3v - x3lo)/(x3hi - x3lo));
      v3 = vpert*cs*sin(tp*(x3v - x3lo)/(x3hi - x3lo))
           *cos(2.0*tp*(x2v - x2lo)/(x2hi - x2lo));
    }
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*v2;
    u0(m,IM3,k,j,i) = d*v3;
    u0(m,IEN,k,j,i) = e + 0.5*d*(v1*v1 + v2*v2 + v3*v3);
    if (etotgrav) u0(m,IEN,k,j,i) += d*PotAt(gm, rin, r);
  });
  if (is_mhd) {
    auto &b = pmbp->pmhd->b0;
    par_for("rg_b", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      b.x1f(m,k,j,i) = 0.0;
      b.x2f(m,k,j,i) = 0.0;
      b.x3f(m,k,j,i) = 0.0;
      if (i == n1m1) b.x1f(m,k,j,i+1) = 0.0;
      if (j == n2m1) b.x2f(m,k,j+1,i) = 0.0;
      if (k == n3m1) b.x3f(m,k+1,j,i) = 0.0;
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantFaceBudget
//! \brief what the two radial faces of the domain carry, in erg/s and g/s against L.
//! The update adds +area(is) F(is) at the bottom and -area(ie+1) F(ie+1) at the top, so
//! both numbers below are signed as a GAIN by the domain.  Called from the source term,
//! which runs after the fluxes and before the next stage overwrites them.

void RedGiantFaceBudget(Mesh *pm) {
  if (face_budget_ <= 0 || pm->ncycle == face_cycle_) return;
  face_cycle_ = pm->ncycle;
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &flx1 = is_mhd ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const Real dt = pm->dt;
  Real sEi = 0.0, sEo = 0.0, sMi = 0.0, sMo = 0.0;
  Kokkos::parallel_reduce("rg_facebud", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
  KOKKOS_LAMBDA(const int m, Real &aEi, Real &aEo, Real &aMi, Real &aMo) {
    const bool inb = (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user);
    const bool oub = (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        if (inb) {
          aEi += area1(m,k,j,is)*flx1(m,IEN,k,j,is);
          aMi += area1(m,k,j,is)*flx1(m,IDN,k,j,is);
        }
        if (oub) {
          aEo -= area1(m,k,j,ie+1)*flx1(m,IEN,k,j,ie+1);
          aMo -= area1(m,k,j,ie+1)*flx1(m,IDN,k,j,ie+1);
        }
      }
    }
  }, sEi, sEo, sMi, sMo);
#if MPI_PARALLEL_ENABLED
  {
    Real g[4] = {sEi, sEo, sMi, sMo};
    MPI_Allreduce(MPI_IN_PLACE, g, 4, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    sEi = g[0]; sEo = g[1]; sMi = g[2]; sMo = g[3];
  }
#endif
  const Real eunit = pmbp->punit->pressure_cgs()*SQR(pmbp->punit->length_cgs())
                     *pmbp->punit->length_cgs();
  const Real munit = pmbp->punit->density_cgs()*SQR(pmbp->punit->length_cgs())
                     *pmbp->punit->length_cgs();
  face_E_in_  += sEi*dt*eunit;
  face_E_out_ += sEo*dt*eunit;
  face_M_in_  += sMi*dt*munit;
  face_M_out_ += sMo*dt*munit;
  if (pm->ncycle % face_budget_ == 0 && global_variable::my_rank == 0) {
    const Real tnow = pm->time*pmbp->punit->time_cgs();
    const Real dtw = tnow - face_t_last_;
    if (dtw > 0.0) {
      const Real iL = 1.0/open_lstar_;
      std::cout << "face budget: gain/L  inner=" << (face_E_in_ - face_E_in_l_)/dtw*iL
                << " outer=" << (face_E_out_ - face_E_out_l_)/dtw*iL
                << " | erg in=" << face_E_in_ << " out=" << face_E_out_
                << " | g in=" << face_M_in_ << " out=" << face_M_out_
                << " (t = " << tnow << " s)" << std::endl;
    }
    face_t_last_ = tnow;
    face_E_in_l_ = face_E_in_;
    face_E_out_l_ = face_E_out_;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantGravity
//! \brief the radial gravity source, -rho G M / r^2, or under wellbalance_dynamic the
//! background's own pressure difference across the cell (face-sum form on the
//! curvilinear grids, which carries the geometric term; plain difference on the column)

void RedGiantGravity(Mesh *pm, Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  RedGiantFaceBudget(pm);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  const bool etotgrav = is_mhd ? pmbp->pmhd->use_etotgrav : pmbp->phydro->use_etotgrav;
  const bool wbdyn = is_mhd ? pmbp->pmhd->use_wellbalance_dynamic
                            : pmbp->phydro->use_wellbalance_dynamic;
  const bool wbx1 = is_mhd ? pmbp->pmhd->use_wb_x1 : pmbp->phydro->use_wb_x1;
  const WBOption wbo = is_mhd ? pmbp->pmhd->wb_option : pmbp->phydro->wb_option;
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = is_mhd ? pmbp->pmhd->phi0.x1f : pmbp->phydro->phi0.x1f;
  DvceArray5D<Real> wbq0 = is_mhd ? pmbp->pmhd->wbq0 : pmbp->phydro->wbq0;
  auto &x1v_ = pmbp->pcoord->x1v;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  const bool curv = curv_;
  const Real gm = gm_, rin = rin_, x1min = x1min_;
  par_for("rg_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    const Real d = w0(m,IDN,k,j,i);
    const Real g = GravAt(gm, r);
    Real src = -bdt*g*d;
    if (!etotgrav) u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
    if (wbdyn) {
      Real pl, pr, d1, d2, d3;
      if (wbx1) {
        WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      } else {
        hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
            w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i), w0(m,IDN,k,j,i+1),
            w0(m,IEN,k,j,i-1), w0(m,IEN,k,j,i), w0(m,IEN,k,j,i+1),
            phicc(m,k,j,i-1), ph1(m,k,j,i), phicc(m,k,j,i), ph1(m,k,j,i+1),
            phicc(m,k,j,i+1), d1, pl, d2, pr, d3);
      }
      if (curv) {
        const Real p = eos.Pressure(d, w0(m,IEN,k,j,i));
        src = bdt*(area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))/volume(m,k,j,i);
      } else {
        src = bdt*(pr - pl)/((x1hi - x1lo)/indcs.nx1);
      }
    }
    u0(m,IM1,k,j,i) += src;
  });

  // --- the mixing-length convective flux (see the header): on each interior radial
  // face, where the face's d ln T / d ln p exceeds grad_ad,
  //     F_conv = rho c_p T sqrt(g delta) l^2 (grad - grad_ad)^(3/2)
  //              / (4 sqrt2 H_p^(3/2)),
  // l = alpha H_p, H_p = p/(rho g), delta = chi_T/chi_rho, c_p = c_v Gamma_1/chi_rho, all
  // in cgs then converted; zero on both walls.  Capped so one step moves at most a tenth
  // of the smaller neighbour's internal energy.
  if (mlt_alpha_ > 0.0) {
    auto fconv = fconv_;
    auto fdiag = fdiag_;
    const Real alpha = mlt_alpha_, rgas = rgas_, igm1 = 1.0/gm1_, gamma = gm1_ + 1.0;
    const Real dunit = pmbp->punit->density_cgs();
    const Real punit_ = pmbp->punit->pressure_cgs();
    const Real vunit = pmbp->punit->velocity_cgs();
    const Real lunit = pmbp->punit->length_cgs();
    const Real tunit = pmbp->punit->time_cgs();
    auto &x1f_ = pmbp->pcoord->xx1f;
    par_for("rg_mlt_flux", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      fconv(m,k,j,i) = 0.0;
      if (m == 0 && k == ks && j == js) {
        for (int q = 0; q < 8; ++q) fdiag(i,q) = 0.0;
      }
      if (i == is || i == ie+1) return;
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xf = curv ? x1f_(m,i) : LeftEdgeX(i-is, indcs.nx1, x1lo, x1hi);
      const Real rf = RadiusOf(curv, xf, rin, x1min);
      const Real dl = w0(m,IDN,k,j,i-1), dr_ = w0(m,IDN,k,j,i);
      const Real el = w0(m,IEN,k,j,i-1), er = w0(m,IEN,k,j,i);
      const Real pl = eos.Pressure(dl, el), pr = eos.Pressure(dr_, er);
      const Real tl = TempKelvin(eos, rgas, dl, el, pl);
      const Real tr = TempKelvin(eos, rgas, dr_, er, pr);
      if (!(pl > pr) || !(tl > tr)) return;             // not stratified the right way
      const Real grad = log(tl/tr)/log(pl/pr);
      const Real pf = 0.5*(pl + pr), tf = 0.5*(tl + tr), df = 0.5*(dl + dr_);
      const Real grad_ad = GradAd(eos, gamma, rgas, pf, tf);
      if (!(grad > grad_ad)) return;
      const Real ef = 0.5*(el + er);
      // thermodynamics at the face, cgs
      const Real g = GravAt(gm, rf)*lunit/(tunit*tunit);
      const Real hp = pf*punit_/(df*dunit*g);
      const Real ell = alpha*hp;
      Real cv = (EintFromDensT(eos, rgas, igm1, df, 1.01*tf)
                 - EintFromDensT(eos, rgas, igm1, df, tf))/(0.01*tf)*punit_/(df*dunit);
      if (!(cv > 0.0)) cv = 1.5*pf*punit_/(df*dunit*tf);
      Real delta = 1.0, cp = gamma*cv;
      if (eos.IsGeneral()) {
        const Real chit = eos.ChiT(df, ef), chir = eos.ChiRho(df, ef);
        delta = chit/chir;
        cp = cv*eos.Gamma1(df, ef)/chir;
      }
      const Real x = grad - grad_ad;
      const Real vc = sqrt(g*delta*ell*ell*x/(8.0*hp));
      Real f = df*dunit*cp*tf*sqrt(g*delta)*ell*ell*x*sqrt(x)
               /(4.0*sqrt(2.0)*hp*sqrt(hp));
      f /= (punit_*vunit);                                // code flux
      // STABILITY.  This is an explicit nonlinear diffusion of T with diffusivity
      //   chi = (dF/d(dT/dr)) / (rho c_p) = (3/2) (F/x) (H_p/T) / (rho c_p),
      // which is unbounded as the superadiabaticity x grows and rho falls (the
      // tabulated column blew up 1e6-fold this way).  An explicit step is stable for
      // chi <= dx^2/(2 dt); the flux is scaled down to hold chi at half that, so a
      // strongly superadiabatic layer under-transports rather than explodes.  A tenth
      // of the smaller internal energy per step caps it as well.
      const Real dxl = curv ? (x1f_(m,i) - x1f_(m,i-1)) : (x1hi - x1lo)/indcs.nx1;
      const Real dxc = dxl*lunit;
      const Real chi = 1.5*(f*punit_*vunit/x)*(hp/tf)/(df*dunit*cp);   // cm^2/s
      const Real chimax = 0.25*dxc*dxc/(bdt*tunit);
      if (chi > chimax) f *= chimax/chi;
      const Real fmax_ = 0.1*fmin(el, er)*dxl/bdt;
      fconv(m,k,j,i) = fmin(f, fmax_);
      if (m == 0 && k == ks && j == js) {
        fdiag(i,0) = grad; fdiag(i,1) = grad_ad; fdiag(i,2) = hp; fdiag(i,3) = cp;
        fdiag(i,4) = vc; fdiag(i,5) = f*punit_*vunit; fdiag(i,6) = fmax_*punit_*vunit;
        fdiag(i,7) = fconv(m,k,j,i)*punit_*vunit;
      }
    });
    if (!mlt_dumped_ && !mlt_dump_.empty() && global_variable::my_rank == 0) {
      mlt_dumped_ = true;
      auto hd = Kokkos::create_mirror_view(fdiag);
      Kokkos::deep_copy(hd, fdiag);
      std::ofstream df(mlt_dump_);
      df.precision(6);
      df << std::scientific;
      df << "# red_giant MLT faces, block 0, first column, first source call (cgs)\n"
         << "# i grad grad_ad H_p c_p v_conv F_mlt F_cap F_used\n";
      for (int i = is; i <= ie+1; ++i) {
        df << i;
        for (int q = 0; q < 8; ++q) df << " " << hd(i,q);
        df << "\n";
      }
    }
    par_for("rg_mlt_div", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      Real div;
      if (curv) {
        div = (area1(m,k,j,i)*fconv(m,k,j,i) - area1(m,k,j,i+1)*fconv(m,k,j,i+1))
              /volume(m,k,j,i);
      } else {
        div = (fconv(m,k,j,i) - fconv(m,k,j,i+1))/((x1hi - x1lo)/indcs.nx1);
      }
      u0(m,IEN,k,j,i) += bdt*div;
    });
  }

  // --- OPEN INNER BOUNDARY: the CO5BOLD (Freytag+ 2012) relaxation of the lowest active
  // layer, the same treatment solar_convection.cpp uses, on a spherical shell.  Four
  // steps: relax UPFLOWING gas toward the deep adiabat (this is the energy input and it
  // sets the emergent luminosity), damp pressure toward the shell mean, restore the mean
  // density, and remove the net radial mass flux.  The means are global over the whole
  // shell at i = is, so they need a reduction across ranks.
  if (open_inner_) {
    auto &mb_bcs = pmbp->pmb->mb_bcs;
    const Real rgas = rgas_, igm1 = 1.0/gm1_, gamma = gm1_ + 1.0;
    const Real punit_ = pmbp->punit->pressure_cgs();
    const Real p_base = p_base_/punit_, t_base = t_base_;
    const Real ccs = cs_change_, ccp = cp_change_;
    auto &x1f_r = pmbp->pcoord->xx1f;
    // Pass 0: cancel the NET mass the inner face delivered this stage.  The fluxes the
    // RKUpdate just applied are still in uflx, so this is exact per stage rather than a
    // correction chasing a drift.
    if (open_nomassflux_) {
      auto &flx1 = is_mhd ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
      auto &area1_ = pmbp->pcoord->area.x1f;
      auto &vol0_ = pmbp->pcoord->volume;
      Real sM = 0.0, sV = 0.0, sH = 0.0, sD = 0.0, sPv = 0.0;
      Kokkos::parallel_reduce("rg_co50", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
      KOKKOS_LAMBDA(const int m, Real &aM, Real &aV, Real &aH, Real &aD, Real &aPv) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real vol = vol0_(m,k,j,is);
            const Real d = u0(m,IDN,k,j,is);
            Real ei = u0(m,IEN,k,j,is)
                      - 0.5*(SQR(u0(m,IM1,k,j,is)) + SQR(u0(m,IM2,k,j,is))
                             + SQR(u0(m,IM3,k,j,is)))/d;
            if (etotgrav) ei -= d*phicc(m,k,j,is);
            aM += area1_(m,k,j,is)*flx1(m,IDN,k,j,is);
            aV += vol;
            aH += (u0(m,IEN,k,j,is) + eos.Pressure(d, ei))*vol;   // total enthalpy
            aD += d*vol;
            aPv += u0(m,IM1,k,j,is)*vol;
          }
        }
      }, sM, sV, sH, sD, sPv);
#if MPI_PARALLEL_ENABLED
      {
        Real g[5] = {sM, sV, sH, sD, sPv};
        MPI_Allreduce(MPI_IN_PLACE, g, 5, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
        sM = g[0]; sV = g[1]; sH = g[2]; sD = g[3]; sPv = g[4];
      }
#endif
      if (sV > 0.0 && sD > 0.0) {
        const Real dd = -bdt*sM/sV;      // code density, the same in every base cell
        const Real hbar = sH/sD;         // shell-mean specific total enthalpy
        const Real vbar = sPv/sD;        // shell-mean radial velocity
        par_for("rg_co50b", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
          if (u0(m,IDN,k,j,is) + dd <= 0.0) return;
          u0(m,IDN,k,j,is) += dd;
          u0(m,IM1,k,j,is) += dd*vbar;
          u0(m,IEN,k,j,is) += dd*hbar;
        });
      }
    }
    // Pass A: shell means of density and pressure
    Real srho0 = 0.0, sumP = 0.0, sN = 0.0;
    Kokkos::parallel_reduce("rg_co5A", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &a, Real &b, Real &c) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real r = w0(m,IDN,k,j,is);
            a += r;
            b += eos.Pressure(r, w0(m,IEN,k,j,is));
            c += 1.0;
          }
        }
      }
    }, srho0, sumP, sN);
#if MPI_PARALLEL_ENABLED
    {
      Real g[3] = {srho0, sumP, sN};
      MPI_Allreduce(MPI_IN_PLACE, g, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      srho0 = g[0]; sumP = g[1]; sN = g[2];
    }
#endif
    const Real meanP = (sN > 0.0) ? sumP/sN : 0.0;
    // Pass B: the inflow entropy relaxation, then the pressure damping.  Both rewrite
    // u0(IEN); the budget (problem/open_budget) accumulates each one separately, which
    // is the only way to see which of them is the star's energy source and which is not.
    auto &vol_ = pmbp->pcoord->volume;
    const Real eunit = pmbp->punit->pressure_cgs()
                       *SQR(pmbp->punit->length_cgs())*pmbp->punit->length_cgs();
    Real dEent = 0.0, dEprs = 0.0;
    Kokkos::parallel_reduce("rg_co5B", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &sEent, Real &sEprs) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      const Real dr_ = x1f_r(m,is+1) - x1f_r(m,is);
      Real r = u0(m,IDN,k,j,is);
      const Real v1 = u0(m,IM1,k,j,is)/r, v2 = u0(m,IM2,k,j,is)/r;
      const Real v3 = u0(m,IM3,k,j,is)/r;
      const Real Ein = u0(m,IEN,k,j,is);
      Real ei = u0(m,IEN,k,j,is) - 0.5*r*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) ei -= r*phicc(m,k,j,is);
      Real es = ei/r;
      Real P = eos.Pressure(r, ei);
      Real g1 = eos.Gamma1(r, ei);
      const Real cs = sqrt(g1*P/r);
      if (v1 > 0.0) {
        // the point on the adiabat through the base state AT THIS CELL'S PRESSURE; a
        // fractional step toward it at constant pressure IS the entropy relaxation
        const Real gad = GradAd(eos, gamma, rgas, p_base*punit_, t_base);
        const Real T_ad = t_base*pow(P/p_base, gad);
        const Real r_ad = DensFromPT(eos, rgas, P, T_ad);
        if (r_ad > 0.0) {
          const Real es_ad = EintFromDensT(eos, rgas, igm1, r_ad, T_ad)/r_ad;
          const Real rlx = ccs*bdt*cs/dr_;
          r  += rlx*(r_ad - r);
          es += rlx*(es_ad - es);
        }
      }
      {
        Real Em = r*es + 0.5*r*(v1*v1 + v2*v2 + v3*v3);
        if (etotgrav) Em += r*phicc(m,k,j,is);
        sEent += (Em - Ein)*vol_(m,k,j,is)*eunit;
        sEprs -= Em*vol_(m,k,j,is)*eunit;    // completed after the pressure damping
      }
      // damp pressure toward the shell mean, adiabatically
      Real P1 = eos.Pressure(r, r*es);
      Real g1p = eos.Gamma1(r, r*es);
      const Real cs2 = g1p*P1/r;
      const Real rlxp = ccp*bdt*sqrt(cs2)/dr_;
      r  += rlxp*(1.0/cs2)*(meanP - P1);
      es += rlxp*(1.0/(g1p*r))*(meanP - P1);
      u0(m,IDN,k,j,is) = r;
      u0(m,IM1,k,j,is) = r*v1;
      u0(m,IM2,k,j,is) = r*v2;
      u0(m,IM3,k,j,is) = r*v3;
      Real E = r*es + 0.5*r*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) E += r*phicc(m,k,j,is);
      u0(m,IEN,k,j,is) = E;
      sEprs += E*vol_(m,k,j,is)*eunit;
      }}
    }, dEent, dEprs);
    // Pass C: the shell means again, after those two steps
    Real srho2 = 0.0, srho2v = 0.0, sv = 0.0;
    Kokkos::parallel_reduce("rg_co5C", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &a, Real &b, Real &c) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            const Real r = u0(m,IDN,k,j,is);
            const Real v1 = u0(m,IM1,k,j,is)/r;
            a += r; b += r*v1; c += v1;
          }
        }
      }
    }, srho2, srho2v, sv);
#if MPI_PARALLEL_ENABLED
    {
      Real g[3] = {srho2, srho2v, sv};
      MPI_Allreduce(MPI_IN_PLACE, g, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      srho2 = g[0]; srho2v = g[1]; sv = g[2];
    }
#endif
    const Real drho4 = (sN > 0.0) ? (srho0 - srho2)/sN : 0.0;
    const Real cvel = (srho0 > 0.0) ? (srho2v + drho4*sv)/srho0 : 0.0;
    // Pass D: restore the mean density, and drive the net radial mass flux to zero
    Real dEden = 0.0;
    Kokkos::parallel_reduce("rg_co5D", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb1+1),
    KOKKOS_LAMBDA(const int m, Real &sEden) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      const Real Ein = u0(m,IEN,k,j,is);
      const Real r0 = u0(m,IDN,k,j,is);
      const Real v1 = u0(m,IM1,k,j,is)/r0 - cvel;
      const Real v2 = u0(m,IM2,k,j,is)/r0, v3 = u0(m,IM3,k,j,is)/r0;
      Real ei = u0(m,IEN,k,j,is) - 0.5*r0*(SQR(u0(m,IM1,k,j,is)/r0) + v2*v2 + v3*v3);
      if (etotgrav) ei -= r0*phicc(m,k,j,is);
      const Real es = ei/r0;
      const Real r = r0 + drho4;
      if (!(r > 0.0)) continue;
      u0(m,IDN,k,j,is) = r;
      u0(m,IM1,k,j,is) = r*v1;
      u0(m,IM2,k,j,is) = r*v2;
      u0(m,IM3,k,j,is) = r*v3;
      Real E = r*es + 0.5*r*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) E += r*phicc(m,k,j,is);
      u0(m,IEN,k,j,is) = E;
      sEden += (E - Ein)*vol_(m,k,j,is)*eunit;
      }}
    }, dEden);
    // The report.  Rates over the interval since the last one, and the running totals,
    // both in units of L: the boundary is supposed to feed the star, so anything that
    // is not O(L) is a bug in one of the three passes above.
#if MPI_PARALLEL_ENABLED
    {
      Real g[3] = {dEent, dEprs, dEden};
      MPI_Allreduce(MPI_IN_PLACE, g, 3, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      dEent = g[0]; dEprs = g[1]; dEden = g[2];
    }
#endif
    // Pass E: give back what Pass D took.  Uniformly in volume, so the shell's own
    // structure is untouched and only the spurious global term goes away.
    if (open_conserve_ && dEden != 0.0) {
      Real svol = 0.0;
      Kokkos::parallel_reduce("rg_co5Evol", Kokkos::RangePolicy<>(DevExeSpace(),0,nmb1+1),
      KOKKOS_LAMBDA(const int m, Real &a) {
        if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
          for (int k=ks; k<=ke; ++k) {
            for (int j=js; j<=je; ++j) { a += vol_(m,k,j,is); }
          }
        }
      }, svol);
#if MPI_PARALLEL_ENABLED
      MPI_Allreduce(MPI_IN_PLACE, &svol, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
      if (svol > 0.0) {
        const Real de = -dEden/(svol*eunit);      // code energy density
        par_for("rg_co5E", DevExeSpace(), 0, nmb1, ks, ke, js, je,
        KOKKOS_LAMBDA(const int m, const int k, const int j) {
          if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
          u0(m,IEN,k,j,is) += de;
        });
        dEden = 0.0;
      }
    }
    if (open_budget_ > 0) {
      open_dE_ent_ += dEent;
      open_dE_prs_ += dEprs;
      open_dE_den_ += dEden;
      if (pm->ncycle % open_budget_ == 0 && global_variable::my_rank == 0) {
        const Real tnow = pm->time*pmbp->punit->time_cgs();
        const Real dtw = tnow - open_t_last_;
        const Real iL = 1.0/open_lstar_;
        if (dtw > 0.0) {
          std::cout << "open BC budget: rate/L  entropy="
                    << (open_dE_ent_ - open_dE_ent_l_)/dtw*iL << " pressure="
                    << (open_dE_prs_ - open_dE_prs_l_)/dtw*iL << " density="
                    << (open_dE_den_ - open_dE_den_l_)/dtw*iL
                    << " | total erg = " << open_dE_ent_ + open_dE_prs_ + open_dE_den_
                    << " (t = " << tnow << " s)" << std::endl;
        }
        open_t_last_ = tnow;
        open_dE_ent_l_ = open_dE_ent_;
        open_dE_prs_l_ = open_dE_prs_;
        open_dE_den_l_ = open_dE_den_;
      }
    }
  }

  // --- THE WALL IS IMPERMEABLE.  See wall_noflux_ above: cancel the mass the inner face
  // advected this stage and the energy it carried, cell by cell, leaving the momentum
  // flux (its pressure term holds the star up) and the conduction flux (it carries the
  // imposed luminosity through the same face) alone.
  if (!open_inner_ && wall_noflux_) {
    auto &flx1w = is_mhd ? pmbp->pmhd->uflx.x1f : pmbp->phydro->uflx.x1f;
    auto &area1w = pmbp->pcoord->area.x1f;
    auto &volw = pmbp->pcoord->volume;
    auto &mb_bcs = pmbp->pmb->mb_bcs;
    par_for("rg_wallflux", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) != BoundaryFlag::user) return;
      const Real dm = bdt*area1w(m,k,j,is)*flx1w(m,IDN,k,j,is)/volw(m,k,j,is);
      if (dm == 0.0) return;
      const Real d = u0(m,IDN,k,j,is);
      if (!(d - dm > 0.0)) return;
      Real ei = u0(m,IEN,k,j,is)
                - 0.5*(SQR(u0(m,IM1,k,j,is)) + SQR(u0(m,IM2,k,j,is))
                       + SQR(u0(m,IM3,k,j,is)))/d;
      if (etotgrav) ei -= d*phicc(m,k,j,is);
      const Real h = (u0(m,IEN,k,j,is) + eos.Pressure(d, ei))/d;   // specific total
      u0(m,IDN,k,j,is) -= dm;
      u0(m,IEN,k,j,is) -= dm*h;
    });
  }

  // --- THE SPONGE under the outer wall.  See sponge_on_ above.
  if (sponge_on_) {
    const Real r0s = pm->mesh_size.x1min, r1s = pm->mesh_size.x1max;
    const Real zs = r0s + sponge_zbot_*(r1s - r0s);
    const Real izw = (r1s > zs) ? 1.0/(r1s - zs) : 0.0;
    const Real sc = sponge_c_;
    auto &x1vs = pmbp->pcoord->x1v;
    auto &dx1s = pmbp->pcoord->dx1;
    const Real gamma_ = gm1_ + 1.0;
    par_for("rg_sponge", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
      const Real xc = curv ? x1vs(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
      if (!(xc > zs)) return;
      const Real ramp = SQR((xc - zs)*izw);
      const Real d = u0(m,IDN,k,j,i);
      Real v1 = u0(m,IM1,k,j,i)/d, v2 = u0(m,IM2,k,j,i)/d, v3 = u0(m,IM3,k,j,i)/d;
      Real ei = u0(m,IEN,k,j,i) - 0.5*d*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) ei -= d*phicc(m,k,j,i);
      const Real cs2 = eos.IsGeneral() ? eos.Gamma1(d, ei)*eos.Pressure(d, ei)/d
                                       : gamma_*gm1_*ei/d;
      const Real dz = curv ? dx1s(m,k,j,i) : (x1hi - x1lo)/indcs.nx1;
      const Real fac = 1.0/(1.0 + sc*ramp*bdt*sqrt(cs2)/dz);
      v1 *= fac; v2 *= fac; v3 *= fac;
      u0(m,IM1,k,j,i) = d*v1;
      u0(m,IM2,k,j,i) = d*v2;
      u0(m,IM3,k,j,i) = d*v3;
      Real E = ei + 0.5*d*(v1*v1 + v2*v2 + v3*v3);
      if (etotgrav) E += d*phicc(m,k,j,i);
      u0(m,IEN,k,j,i) = E;
    });
  }

  // --- the optically thin layers: the correlated-k or the grey two-stream if one of
  // them is on, else the old grey Eddington relaxation
  if (rt_ck_ || rt_grey_) {
    if (!ck_dumped2_ && ck_dump_t2_ >= 0.0 && !ck_dump_file2_.empty() &&
        pm->time >= ck_dump_t2_) {
      ck_dumped2_ = true;
      two_stream_rt::rt_dump_file = ck_dump_file2_;
      two_stream_rt::rt_dump_done = false;      // re-arm the one-shot dump
    }
    two_stream_rt::picket_fence_two_stream_RT(pm, bdt);
  }
  // --- the grey relaxation (see UserProblem): with weight 1 - w each cell
  // decays toward the Eddington temperature of its optical depth on its radiative time
  if (relax_ && !rt_ck_ && !rt_grey_) {
    Conduction *pc = is_mhd ? pmbp->pmhd->pcond : pmbp->phydro->pcond;
    auto &wf = pc->rad_w;
    auto &tf = pc->rad_tauf;
    auto ktab = ktab_;
    auto klT = klT_;
    auto klD = klD_;
    const int knT = knT_, knD = knD_;
    const Real teff4 = SQR(SQR(teff_));
    const Real kfac = kfac_;
    const Real rgas = rgas_, igm1 = 1.0/gm1_, gm1 = gm1_;
    const Real dunit = pmbp->punit->density_cgs();
    const Real punit_ = pmbp->punit->pressure_cgs();
    const Real dt_cgs = bdt*pmbp->punit->time_cgs();
    par_for("rg_relax", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real wc = 0.5*(wf(m,k,j,i) + wf(m,k,j,i+1));
      if (wc >= 1.0) return;
      const Real tauc = 0.5*(tf(m,k,j,i) + tf(m,k,j,i+1));
      const Real d = w0(m,IDN,k,j,i);
      const Real e = w0(m,IEN,k,j,i);
      const Real p = eos.Pressure(d, e);
      const Real t = TempKelvin(eos, rgas, d, e, p);
      const Real teq = pow(0.75*teff4*(tauc + 2.0/3.0), 0.25);
      const Real eeq = EintFromDensT(eos, rgas, igm1, d, teq);
      // rho c_v [erg/cm^3/K] by a finite difference of the EOS's own e(rho, T)
      const Real e1 = EintFromDensT(eos, rgas, igm1, d, 1.01*t);
      const Real e0 = EintFromDensT(eos, rgas, igm1, d, t);
      // fall back to the ideal-gas value if the table's finite difference misbehaves
      Real rcv = (e1 - e0)/(0.01*t)*punit_;
      if (!(rcv > 0.0)) rcv = 1.5*p*punit_/t;
      const Real kp = kfac*KappaTab(ktab, klT, klD, knT, knD, t, d*dunit);
      const Real trad = rcv/(4.0*kp*(d*dunit)*kSigmaSB*t*t*t);
      if (!(trad > 0.0)) return;
      const Real de = -(1.0 - wc)*(e - eeq)*(1.0 - exp(-dt_cgs/trad));
      u0(m,IEN,k,j,i) += de;
      (void) gm1;
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantBC
//! \brief both radial walls: the ghost cells hold the initial column at their own
//! radius, with the radial velocity mirrored (a reflecting wall whose ghost state is in
//! hydrostatic balance with the interior, so the well-balanced scheme sees no jump)

void RedGiantBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  const bool is_mhd = (pmbp->pmhd != nullptr);
  auto &u0 = is_mhd ? pmbp->pmhd->u0 : pmbp->phydro->u0;
  auto &w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  auto eos = is_mhd ? pmbp->pmhd->peos->eos_data : pmbp->phydro->peos->eos_data;
  auto &x1v_ = pmbp->pcoord->x1v;
  const bool curv = curv_, etotgrav = etotgrav_;
  const Real gm = gm_, rin = rin_, x1min = x1min_, rgas = rgas_, igm1 = 1.0/gm1_;
  DvceArray4D<Real> phicc = is_mhd ? pmbp->pmhd->phicc0 : pmbp->phydro->phicc0;
  const Real rlo = rlo_, drf = drf_;
  const int nfine = nfine_;
  auto lnp = lnp_d_;
  auto tk = tk_d_;
  // OPEN inner boundary: continue the interior hydrostatically at ITS OWN temperature and
  // copy the velocity, so a plume crosses the boundary instead of bouncing off it.  The
  // pressure follows d ln p = -(rho/p) dPhi from the lowest active cell, which keeps the
  // ghost in balance with the interior rather than with the initial column -- tying it to
  // the column would fight the relaxation and pressurise the envelope.
  const bool open_in = open_inner_;
  // problem/open_debug: print the open ghost's inputs and outputs for one column on the
  // first few boundary calls.  A ghost that comes out non-finite says nothing about
  // which of the four EOS calls did it.
  const bool open_dbg = (open_debug_ > 0) && (open_dbg_calls_++ < open_debug_);
  auto fill_open = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                                 const int im) {
    const Real d_i = w0(m,IDN,k,j,im);
    const Real e_i = w0(m,IEN,k,j,im);
    const Real p_i = eos.Pressure(d_i, e_i);
    const Real t_i = TempKelvin(eos, rgas, d_i, e_i, p_i);
    const Real dphi = phicc(m,k,j,i) - phicc(m,k,j,im);
    const Real p_g = p_i*exp(-(d_i/p_i)*dphi);
    const Real d_g = DensFromPT(eos, rgas, p_g, t_i);
    const Real e_g = EintFromDensT(eos, rgas, igm1, d_g, t_i);
    const Real v1 = w0(m,IVX,k,j,im), v2 = w0(m,IVY,k,j,im), v3 = w0(m,IVZ,k,j,im);
    if (open_dbg && m == 0 && k == 2 && j == 2) {
      Kokkos::printf("rg_open_bc i=%d im=%d: d_i=%.6e e_i=%.6e p_i=%.6e t_i=%.6e "
                     "dphi=%.6e p_g=%.6e d_g=%.6e e_g=%.6e v1=%.6e\n",
                     i, im, d_i, e_i, p_i, t_i, dphi, p_g, d_g, e_g, v1);
    }
    w0(m,IDN,k,j,i) = d_g;
    w0(m,IEN,k,j,i) = e_g;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    u0(m,IDN,k,j,i) = d_g;
    u0(m,IM1,k,j,i) = d_g*v1;
    u0(m,IM2,k,j,i) = d_g*v2;
    u0(m,IM3,k,j,i) = d_g*v3;
    Real et = e_g + 0.5*d_g*(v1*v1 + v2*v2 + v3*v3);
    if (etotgrav) et += d_g*phicc(m,k,j,i);
    u0(m,IEN,k,j,i) = et;
  };
  auto fill = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            const int im) {
    const Real x1lo = size.d_view(m).x1min, x1hi = size.d_view(m).x1max;
    const Real xc = curv ? x1v_(m,i) : CellCenterX(i-is, indcs.nx1, x1lo, x1hi);
    const Real r = RadiusOf(curv, xc, rin, x1min);
    Real lp, t;
    ColumnAt(lnp, tk, nfine, rlo, drf, r, lp, t);
    const Real p = exp(lp);
    const Real d = DensFromPT(eos, rgas, p, t);
    const Real e = EintFromDensT(eos, rgas, igm1, d, t);
    const Real v1 = -w0(m,IVX,k,j,im), v2 = w0(m,IVY,k,j,im), v3 = w0(m,IVZ,k,j,im);
    w0(m,IDN,k,j,i) = d;
    w0(m,IEN,k,j,i) = e;
    w0(m,IVX,k,j,i) = v1;
    w0(m,IVY,k,j,i) = v2;
    w0(m,IVZ,k,j,i) = v3;
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = d*v2;
    u0(m,IM3,k,j,i) = d*v3;
    Real et = e + 0.5*d*(v1*v1 + v2*v2 + v3*v3);
    if (etotgrav) et += d*PotAt(gm, rin, r);
    u0(m,IEN,k,j,i) = et;
  };
  par_for("rg_bc_x1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int n) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      if (open_in) {
        fill_open(m, k, j, is-1-n, is);      // always from the lowest ACTIVE cell
      } else {
        fill(m, k, j, is-1-n, is+n);
      }
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      fill(m, k, j, ie+1+n, ie-n);
    }
  });
  if (is_mhd) {
    auto &b = pmbp->pmhd->b0;
    const int n1m1 = indcs.nx1 + 2*ng - 1;
    par_for("rg_bc_b", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (i < is || i > ie) {
        b.x1f(m,k,j,i) = 0.0;
        b.x2f(m,k,j,i) = 0.0;
        b.x3f(m,k,j,i) = 0.0;
        if (i == n1m1) b.x1f(m,k,j,i+1) = 0.0;
        if (j == n2m1) b.x2f(m,k,j+1,i) = 0.0;
        if (k == n3m1) b.x3f(m,k+1,j,i) = 0.0;
      }
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RedGiantFinal

void RedGiantFinal(ParameterInput *pin, Mesh *pm) {
  // every file-scope View must be released before Kokkos::finalize, or it aborts
  lnp_d_ = DvceArray1D<Real>();
  tk_d_ = DvceArray1D<Real>();
  ktab_ = DvceArray2D<Real>();
  fconv_ = DvceArray4D<Real>();
  fdiag_ = DvceArray2D<Real>();
  klT_ = DvceArray1D<Real>();
  klD_ = DvceArray1D<Real>();
  return;
}
