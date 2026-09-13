//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file box_convection.cpp
//! \brief A Cartesian plane-parallel LOCAL BOX of compressible convection, heated by a
//! fixed radiative flux through the bottom wall and cooled by a Newton-relaxation layer
//! at the top.  Written for the base of a red-giant convection zone, but nothing in it is
//! specific to that star: the state is fixed by (g0, rho_base, t_base) plus the EOS.
//!
//! GEOMETRY.  x1 is the vertical, gravity is the CONSTANT -g0 x1-hat, so the potential is
//! Phi(z) = g0 (z - x1min), increasing upward as the code's well-balanced machinery
//! expects.  x2 and x3 are periodic.  2-D (nx3 = 1) and 3-D both work.
//!
//! INITIAL STATE.  An exactly ISENTROPIC hydrostatic column, integrated with the RUN's
//! OWN EOS by marching in (p, T):
//!     dln p/dz = -rho g0/p,     dln T/dz = grad_ad(p,T) dln p/dz,
//!     rho = DensFromPT(p, T),   e = EintFromDensT(rho, T),
//! with grad_ad and DensFromPT taken from the tabulated EOS (pgen_eos_utils).  The march
//! runs on a fine grid covering the mesh plus its ghosts, in a ONE-THREAD DEVICE KERNEL:
//! a tabulated EOS lives in a DvceArray and cannot be evaluated from host code.  Because
//! the column is isentropic, `<hydro>/wb_option = isentropic` is the EXACT closure for
//! the dynamic well-balanced scheme, and the hydrostatic residual is at round-off.
//!
//! NOTE that a constant-g adiabat has a finite height: T falls linearly and reaches zero
//! at z = H_p(base)/grad_ad, i.e. 2.5 base scale heights for grad_ad = 0.4.  The domain
//! MUST end below that; the start-up report prints where it is.
//!
//! HEATING.  The bottom wall carries a fixed inward radiative flux, supplied by the
//! conduction module: set `<hydro>/rad_flux_inner` to the flux in erg/cm^2/s, with
//! ix1_bc = user (conduction.cpp adds it at the i = is face, and that path is
//! geometry-free -- the Cartesian RKUpdate divides by dx1).  The Rosseland opacity table
//! is read here and installed into the Conduction object, exactly as red_giant.cpp does.
//!
//! COOLING.  In the top `cool_depth` of the box the SPECIFIC internal energy is relaxed
//! toward the initial column's on a timescale `cool_tau`, with a smoothstep ramp:
//!     d(u_IEN)/dt = -ramp(z) rho (eps - eps_0(z))/cool_tau.
//! Relaxing the specific energy rather than the energy density keeps the sink from
//! forcing the density.  In steady state the layer removes exactly what the bottom
//! puts in.
//!
//! FLUX BOOST.  At a star's true flux the convective Mach number here is ~1e-3 and a
//! turnover costs 1e6-1e7 sound-limited steps.  Nothing in this file knows about that:
//! the boost is applied by the INPUT FILE, by raising `<hydro>/rad_flux_inner` and
//! lowering `<hydro>/rad_kappa_fac` by the same factor (so the radiative and convective
//! fluxes keep their ratio).  v ~ F^(1/3) and (nabla - nabla_ad) ~ F^(2/3) scale the
//! answer back.
//!
//! <problem> keys
//!   g0            constant gravity [cm/s^2] (positive; points to x1min)
//!   rho_base      density at x1min [g/cm^3]
//!   t_base        temperature at x1min [K]
//!   dgrad         SUPERADIABATIC EXCESS of the initial column: the march uses
//!                 dln T/dln p = grad_ad + dgrad instead of grad_ad.  0 gives the exact
//!                 adiabat, for which `wb_option = isentropic` is an exact closure and
//!                 the hydrostatic residual is at round-off.  A small positive value
//!                 starts the box CLOSE TO the convective steady state: building an
//!                 excess from an adiabat is limited by the imposed flux and takes a
//!                 Kelvin-Helmholtz time (hundreds of turnovers), whereas SHEDDING one is
//!                 done by the convection itself in a few turnovers.  The WB closure
//!                 stays second-order accurate at any dgrad; at 1e-2 the residual is
//!                 ~1e-4 of the gravity term.
//!   cool_depth    thickness of the cooling layer below x1max [cm]; < 0 -> 0.3 H_p(base)
//!   cool_tau      relaxation time of that layer [s]; < 0 -> 0.1 H_p(base)/v*
//!   bc_mode       the x1 walls.  0: the ghost carries the INITIAL column with the normal
//!                 velocity mirrored (what wb_column.cpp does -- correct only while the
//!                 state stays on that column; in a convecting box the ghost and the
//!                 evolved interior drift apart and the wall runs a steady WIND through
//!                 the box).  1: a plain reflecting mirror, which makes the wall-face
//!                 Riemann problem exactly symmetric and the mass flux exactly zero.
//!                 2 (DEFAULT): the mirror RESCALED by the initial column's own ratio
//!                 across the wall, rho_g = rho_m rho_col(z_g)/rho_col(z_m) and likewise
//!                 for e -- impermeable like 1 up to the stratification over one cell,
//!                 and hydrostatic like 0 at t = 0.
//!   ic_profile    if set, a text file "z rho eint" (cgs, one node per line, increasing
//!                 z, '#' comments) REPLACES the isentropic march.  It must already cover
//!                 the ghosts -- analysis/mkprofile.py writes the horizontally and time
//!                 averaged profile of a finished case and pads both ends
//!                 hydrostatically.  This is how a finer case is started from a coarser
//!                 one's RELAXED stratification: an AMR-style upsampled restart is not
//!                 possible across different meshes, but the 1-D mean profile is exactly
//!                 the part that takes a thermal time to establish.  The flow itself
//!                 still has to grow from the seed, which takes a few turnovers, not a
//!                 Kelvin-Helmholtz time.
//!   vpert         velocity seed amplitude, in units of the LOCAL sound speed
//!   vpert_nk      number of random horizontal modes (default 16)
//!   vpert_seed    RNG seed for those modes (default 1234)
//!   nfine         nodes in the column march (default 8192)
//!   opac_table    Rosseland table, "# nT nD lTmin dlT lDmin dlD" then nT*nD log10 kappa
//!   column_dump   if set, write the initial column to this file
//!   mu            ideal-gas branch only: mean molecular weight
//!   user_srcs     must be true (gravity and the cooling layer live in the source term)

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "coordinates/cell_locations.hpp"
#include "mesh/mesh.hpp"
#include "globals.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "utils/wb_background.hpp"
#include "diffusion/conduction.hpp"
#include "units/units.hpp"
#include "pgen_eos_utils.hpp"
#include "pgen.hpp"

void BoxConvSrcs(Mesh *pm, Real bdt);
void BoxConvBC(Mesh *pm);
void BoxConvFinal(ParameterInput *pin, Mesh *pm);

namespace {
// the column, on a uniform fine grid covering the mesh plus its ghosts
DvceArray1D<Real> cd_, ce_, cp_, ct_;   // density, eint, pressure, temperature [K]
Real g0_ = 0.0, zlo_ = 0.0, dzf_ = 1.0, zmin_ = 0.0;
Real zcool_ = 0.0, zmax_ = 0.0, tcool_ = 1.0;
int nfine_ = 0, bc_mode_ = 2;
bool etotgrav_ = false;

//----------------------------------------------------------------------------------------
//! \fn ReadOpacityTable
//! \brief the merged Rosseland table: comment lines, one of them
//! "# nT nD lTmin dlT lDmin dlD", then nT*nD values of log10 kappa_R with T slowest.
//! Same format red_giant.cpp reads; the box never goes near the table edges, so none of
//! that file's edge repair is carried over.

void ReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                      DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD) {
  std::ifstream f(fname);
  if (!f.good()) {
    std::cout << "### FATAL ERROR in box_convection: cannot open problem/opac_table '"
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
    std::cout << "### FATAL ERROR in box_convection: opacity table '" << fname
              << "' has no grid line, or " << vals.size() << " values for "
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
  return;
}
//----------------------------------------------------------------------------------------
//! \fn par_file_tp
//! \brief p and T of a supplied (rho, e) column, for the start-up report and the dump.
//! A tabulated EOS lives in a DvceArray, so this has to happen on the device.

void par_file_tp(const EOS_Data &eos, DvceArray1D<Real> d, DvceArray1D<Real> e,
                 DvceArray1D<Real> p, DvceArray1D<Real> t, const int n, const Real gm1,
                 const Real rgas) {
  par_for("boxconv_filetp", DevExeSpace(), 0, n-1, KOKKOS_LAMBDA(const int i) {
    Real pp, tk;
    pgen_eos::PresTempFromEint(eos, gm1, rgas, d(i), e(i), -1.0, pp, tk);
    p(i) = pp;
    t(i) = tk;
  });
  return;
}
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  user_srcs_func = BoxConvSrcs;
  user_bcs_func = BoxConvBC;
  pgen_final_func = BoxConvFinal;
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (pmbp->phydro == nullptr) {
    std::cout << "### FATAL ERROR in box_convection: <hydro> is required" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (pin->GetOrAddBoolean("mesh", "use_cubed_sphere", false) ||
      pin->GetOrAddBoolean("mesh", "use_spherical_polar", false)) {
    std::cout << "### FATAL ERROR in box_convection: Cartesian meshes only" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  auto &indcs = pmy_mesh_->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, js = indcs.js, ks = indcs.ks;
  const int n1m1 = indcs.nx1 + 2*ng - 1;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;

  const Real g0 = pin->GetReal("problem", "g0");
  const Real rho_b = pin->GetReal("problem", "rho_base");
  const Real t_b = pin->GetReal("problem", "t_base");
  const Real vpert = pin->GetOrAddReal("problem", "vpert", 0.0);
  const int nk = pin->GetOrAddInteger("problem", "vpert_nk", 16);
  const int kseed = pin->GetOrAddInteger("problem", "vpert_seed", 1234);
  const int nfine = pin->GetOrAddInteger("problem", "nfine", 8192);
  const Real mu = pin->GetOrAddReal("problem", "mu", 1.3);
  const Real dgrad = pin->GetOrAddReal("problem", "dgrad", 0.0);
  bc_mode_ = pin->GetOrAddInteger("problem", "bc_mode", 2);
  if (bc_mode_ < 0 || bc_mode_ > 2) {
    std::cout << "### FATAL ERROR in box_convection: problem/bc_mode must be 0, 1 or 2"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  const std::string dump = pin->GetOrAddString("problem", "column_dump", "");
  const std::string icprof = pin->GetOrAddString("problem", "ic_profile", "");
  g0_ = g0;

  auto &eos = pmbp->phydro->peos->eos_data;
  auto &u0 = pmbp->phydro->u0;
  const Real gamma = eos.gamma;
  const Real igm1 = 1.0/(gamma - 1.0);
  const bool etotgrav = pmbp->phydro->use_etotgrav;
  const bool wbdyn = pmbp->phydro->use_wellbalance_dynamic;
  etotgrav_ = etotgrav;
  // ideal-gas branch only; the general branch carries composition in the table
  Real vunit = 1.0, lunit = 1.0, punit = 1.0;
  if (pmbp->punit != nullptr) {
    vunit = pmbp->punit->velocity_cgs();
    lunit = pmbp->punit->length_cgs();
    punit = pmbp->punit->pressure_cgs();
  }
  const Real rgas = 1.380649e-16/(mu*1.66053906660e-24)/(vunit*vunit);

  // --- the column, on a fine grid over the mesh's x1 extent plus its ghosts
  const Real zmin = pmy_mesh_->mesh_size.x1min;
  const Real zmax = pmy_mesh_->mesh_size.x1max;
  const Real dz = (zmax - zmin)/pmy_mesh_->mesh_indcs.nx1;
  const Real zlo = zmin - (ng + 1)*dz, zhi = zmax + (ng + 1)*dz;
  const Real dzf = (zhi - zlo)/(nfine - 1);
  zlo_ = zlo; dzf_ = dzf; zmin_ = zmin; zmax_ = zmax; nfine_ = nfine;
  DualArray1D<Real> cd("cd", nfine), ce("ce", nfine), cp("cp", nfine), ct("ct", nfine);
  auto cd_d = cd.d_view, ce_d = ce.d_view, cp_d = cp.d_view, ct_d = ct.d_view;
  const int i0 = static_cast<int>((zmin - zlo)/dzf + 0.5);
  // ONE THREAD: the tabulated EOS cannot be evaluated on the host, and the march is
  // sequential anyway.  RK2 (midpoint) in (ln p, T), isentropic by construction.
  par_for("boxconv_col", DevExeSpace(), 0, 0, KOKKOS_LAMBDA(int) {
    // the base pressure from (rho_base, t_base), through the internal energy: the general
    // EOS evaluates e(rho,T) and p(rho,e,T) directly, so no root find is needed here
    const Real e_b = pgen_eos::EintFromDensT(eos, rgas, igm1, rho_b, t_b);
    Real p = pgen_eos::PresFromEint(eos, gamma - 1.0, rho_b, e_b);
    Real tk = t_b;
    cp_d(i0) = p;
    ct_d(i0) = tk;
    // upward
    Real pc = p, tc = tk;
    for (int i=i0+1; i<nfine; ++i) {
      Real d = pgen_eos::DensFromPT(eos, rgas, pc, tc);
      Real gad = pgen_eos::GradAd(eos, gamma, rgas, pc, tc) + dgrad;
      Real dlnp = -d*g0/pc;
      Real pm = pc*exp(0.5*dzf*dlnp);
      Real tm = tc + 0.5*dzf*tc*gad*dlnp;
      d = pgen_eos::DensFromPT(eos, rgas, pm, tm);
      gad = pgen_eos::GradAd(eos, gamma, rgas, pm, tm) + dgrad;
      dlnp = -d*g0/pm;
      pc = pc*exp(dzf*dlnp);
      tc = tc + dzf*tm*gad*dlnp;
      cp_d(i) = pc;
      ct_d(i) = tc;
    }
    // downward
    pc = p; tc = tk;
    for (int i=i0-1; i>=0; --i) {
      Real d = pgen_eos::DensFromPT(eos, rgas, pc, tc);
      Real gad = pgen_eos::GradAd(eos, gamma, rgas, pc, tc) + dgrad;
      Real dlnp = -d*g0/pc;
      Real pm = pc*exp(-0.5*dzf*dlnp);
      Real tm = tc - 0.5*dzf*tc*gad*dlnp;
      d = pgen_eos::DensFromPT(eos, rgas, pm, tm);
      gad = pgen_eos::GradAd(eos, gamma, rgas, pm, tm) + dgrad;
      dlnp = -d*g0/pm;
      pc = pc*exp(-dzf*dlnp);
      tc = tc - dzf*tm*gad*dlnp;
      cp_d(i) = pc;
      ct_d(i) = tc;
    }
  });
  par_for("boxconv_colde", DevExeSpace(), 0, nfine-1, KOKKOS_LAMBDA(const int i) {
    const Real d = pgen_eos::DensFromPT(eos, rgas, cp_d(i), ct_d(i));
    cd_d(i) = d;
    ce_d(i) = pgen_eos::EintFromDensT(eos, rgas, igm1, d, ct_d(i));
  });
  cd.modify_device();  cd.sync_host();
  ce.modify_device();  ce.sync_host();
  cp.modify_device();  cp.sync_host();
  ct.modify_device();  ct.sync_host();
  // --- an externally supplied stratification REPLACES the march
  if (!icprof.empty()) {
    std::ifstream pf(icprof);
    if (!pf.good()) {
      std::cout << "### FATAL ERROR in box_convection: cannot open problem/ic_profile '"
                << icprof << "'" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::vector<Real> zf, df, ef;
    std::string line;
    while (std::getline(pf, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      Real a, b, c;
      if (!(is >> a >> b >> c)) continue;
      zf.push_back(a);
      df.push_back(b);
      ef.push_back(c);
    }
    if (zf.size() < 2 || zf.front() > zlo || zf.back() < zhi) {
      std::cout << "### FATAL ERROR in box_convection: problem/ic_profile has "
                << zf.size() << " nodes spanning ["
                << (zf.empty() ? 0.0 : zf.front()) << ", "
                << (zf.empty() ? 0.0 : zf.back())
                << "], which does not cover the mesh plus its ghosts ["
                << zlo << ", " << zhi << "].  Pad it with analysis/mkprofile.py."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::size_t kk = 0;
    for (int i=0; i<nfine; ++i) {
      const Real z = zlo + i*dzf;
      while (kk + 2 < zf.size() && zf[kk+1] < z) ++kk;
      const Real w = (z - zf[kk])/(zf[kk+1] - zf[kk]);
      // logarithmic in both: the profile spans two decades in density
      cd.h_view(i) = std::exp(std::log(df[kk])*(1.0 - w) + std::log(df[kk+1])*w);
      ce.h_view(i) = std::exp(std::log(ef[kk])*(1.0 - w) + std::log(ef[kk+1])*w);
    }
    cd.modify_host();  cd.sync_device();
    ce.modify_host();  ce.sync_device();
    // T and p of the supplied state, for the report and the dump only
    {
      auto cd_dv = cd.d_view, ce_dv = ce.d_view, cp_dv = cp.d_view, ct_dv = ct.d_view;
      par_file_tp(eos, cd_dv, ce_dv, cp_dv, ct_dv, nfine, gamma - 1.0, rgas);
    }
    cp.modify_device();  cp.sync_host();
    ct.modify_device();  ct.sync_host();
    if (global_variable::my_rank == 0) {
      std::cout << "box_convection: initial stratification READ FROM " << icprof
                << " (" << zf.size() << " nodes); the isentropic march is overridden"
                << std::endl;
    }
  }
  cd_ = cd.d_view; ce_ = ce.d_view; cp_ = cp.d_view; ct_ = ct.d_view;

  // --- the derived scales of the base state, and the cooling layer
  const Real p_b = cp.h_view(i0);
  const Real hp0 = p_b/(rho_b*g0);
  const Real g1 = pgen_eos::HostGamma1FromP(eos, rho_b, p_b);
  // grad_ad at the base, read back off the column: dlnT/dlnp between the two nodes
  const Real gad_b = std::log(ct.h_view(i0+1)/ct.h_view(i0-1))
                     /std::log(cp.h_view(i0+1)/cp.h_view(i0-1));
  const Real cs0 = std::sqrt(g1*p_b/rho_b);
  Conduction *pc = pmbp->phydro->pcond;
  Real fin = 0.0;
  if (pc != nullptr) fin = pc->rad_flux_inner*punit*vunit;   // back to erg/cm^2/s
  const Real vstar = (fin > 0.0) ? std::cbrt(fin/rho_b) : cs0;
  const Real tturn = hp0/vstar;
  Real cdep = pin->GetOrAddReal("problem", "cool_depth", -1.0);
  if (cdep < 0.0) cdep = 0.3*hp0;
  Real ctau = pin->GetOrAddReal("problem", "cool_tau", -1.0);
  if (ctau < 0.0) ctau = 0.1*tturn;
  zcool_ = zmax - cdep;
  tcool_ = ctau;
  if (zcool_ <= zmin) {
    std::cout << "### FATAL ERROR in box_convection: cooling layer fills the box"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // --- the Rosseland table, for the conduction operator
  const std::string opac = pin->GetOrAddString("problem", "opac_table", "");
  if (pc != nullptr && pc->iso_cond_type.compare("radiative") == 0) {
    if (!opac.empty()) {
      if (!(pc->rad_kappa_tab && pc->rad_kappa_rho)) {
        std::cout << "### FATAL ERROR in box_convection: use <hydro>/rad_kappa_src = "
                  << "table_rho, the stellar table is on (T, rho)" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      DvceArray2D<Real> ktab;
      DvceArray1D<Real> klT, klD;
      int knT = 0, knD = 0;
      ReadOpacityTable(opac, ktab, klT, klD, knT, knD);
      pc->rad_kr_tab = ktab;
      pc->rad_kr_lT = klT;
      pc->rad_kr_lP = klD;
      pc->rad_kr_nT = knT;
      pc->rad_kr_nP = knD;
    }
    if (pc->rad_kappa_rmax > 0.0) {
      std::cout << "### FATAL ERROR in box_convection: <hydro>/rad_kappa_rmax compares "
                << "against x1v, which a Cartesian mesh never allocates" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // --- the start-up report: every number the design rests on
  if (global_variable::my_rank == 0) {
    std::printf("box_convection: base rho = %.5e g/cm^3, T = %.5e K, p = %.5e\n",
                rho_b, t_b, p_b);
    std::printf("  g0 = %.5e cm/s^2, H_p = %.5e cm, Gamma_1 = %.5f, c_s = %.5e cm/s\n",
                g0, hp0, g1, cs0);
    std::printf("  domain x1 = [%.5e, %.5e] = %.4f H_p; a constant-g adiabat reaches "
                "T = 0 at z - x1min ~ %.5e cm = %.3f H_p\n", zmin, zmax, (zmax-zmin)/hp0,
                hp0/gad_b, 1.0/gad_b);
    std::printf("  initial column: dln T/dln p = grad_ad + %.5e (dgrad; 0 = exact "
                "adiabat)\n", dgrad);
    std::printf("  column top:  T = %.5e K, rho = %.5e, p = %.5e, ln(p_b/p_top) = %.4f\n",
                ct.h_view(nfine-1-(ng+1)), cd.h_view(nfine-1-(ng+1)),
                cp.h_view(nfine-1-(ng+1)),
                std::log(p_b/cp.h_view(nfine-1-(ng+1))));
    std::printf("  F_bot = %.5e erg/cm^2/s -> v* = %.5e cm/s, Mach = %.4e, "
                "turnover H_p/v* = %.5e s\n", fin, vstar, vstar/cs0, tturn);
    std::printf("  cooling layer: z > %.5e (top %.4f H_p), tau = %.5e s = %.4f"
                " turnover\n", zcool_, cdep/hp0, ctau, ctau/tturn);
    std::printf("  x1 walls: bc_mode = %d (0 column ghost, 1 mirror, 2 mirror x the"
                " column ratio)\n", bc_mode_);
    if (pc != nullptr) {
      std::printf("  rad_kappa_fac = %.5e (conductivity is 1/rad_kappa_fac x physical)\n",
                  pc->rad_kappa_fac);
    }
    if (ct.h_view(nfine-1) <= 0.0) {
      std::printf("### FATAL: the adiabat runs out of temperature inside the domain;"
                  " shorten x1max\n");
      std::exit(EXIT_FAILURE);
    }
  }
  if (!dump.empty() && global_variable::my_rank == 0) {
    std::ofstream fo(dump);
    fo << "# z rho T p eint\n";
    for (int i=0; i<nfine; ++i) {
      fo << (zlo + i*dzf) << " " << cd.h_view(i) << " " << ct.h_view(i) << " "
         << cp.h_view(i) << " " << ce.h_view(i) << "\n";
    }
  }
  if (restart) return;

  // --- the random horizontal modes of the velocity seed
  DualArray2D<Real> md("md", std::max(nk,1), 4);   // k2, k3, amplitude, phase
  {
    std::mt19937 rng(kseed);
    std::uniform_real_distribution<double> u01(0.0, 1.0);
    const bool three_d = (indcs.nx3 > 1);
    Real norm = 0.0;
    for (int n=0; n<nk; ++n) {
      md.h_view(n,0) = 1.0 + std::floor(4.0*u01(rng));
      md.h_view(n,1) = three_d ? (1.0 + std::floor(4.0*u01(rng))) : 0.0;
      md.h_view(n,2) = u01(rng) + 0.25;
      md.h_view(n,3) = 2.0*M_PI*u01(rng);
      norm += md.h_view(n,2)*md.h_view(n,2);
    }
    norm = (norm > 0.0) ? 1.0/std::sqrt(norm) : 1.0;
    for (int n=0; n<nk; ++n) md.h_view(n,2) *= norm;
  }
  md.modify_host();
  md.sync_device();
  auto md_d = md.d_view;

  // --- the state and the potentials, every cell and every face, ghosts included
  DvceArray4D<Real> phicc = pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = pmbp->phydro->phi0.x1f;
  DvceArray4D<Real> ph2 = pmbp->phydro->phi0.x2f;
  DvceArray4D<Real> ph3 = pmbp->phydro->phi0.x3f;
  const bool have_phi = (etotgrav || wbdyn);
  const Real x2min_m = pmy_mesh_->mesh_size.x2min, x2max_m = pmy_mesh_->mesh_size.x2max;
  const Real x3min_m = pmy_mesh_->mesh_size.x3min, x3max_m = pmy_mesh_->mesh_size.x3max;
  const Real lz = zmax - zmin;
  par_for("boxconv_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real x2min = size.d_view(m).x2min, x2max = size.d_view(m).x2max;
    const Real x3min = size.d_view(m).x3min, x3max = size.d_view(m).x3max;
    const Real z = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real x2v = CellCenterX(j-js, indcs.nx2, x2min, x2max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, x3min, x3max);
    const Real x1l = LeftEdgeX(i-is, indcs.nx1, x1min, x1max);
    const Real x1r = LeftEdgeX(i+1-is, indcs.nx1, x1min, x1max);
    Real s = (z - zlo)/dzf;
    int ii = static_cast<int>(s);
    ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
    const Real f = s - ii;
    const Real d = cd_d(ii)*(1.0 - f) + cd_d(ii+1)*f;
    const Real e = ce_d(ii)*(1.0 - f) + ce_d(ii+1)*f;
    const Real p = cp_d(ii)*(1.0 - f) + cp_d(ii+1)*f;
    Real v1 = 0.0;
    if (vpert > 0.0 && z > zmin && z < zmax) {
      const Real cs = sqrt(gamma*p/d);
      Real amp = 0.0;
      for (int n=0; n<nk; ++n) {
        amp += md_d(n,2)*sin(2.0*M_PI*(md_d(n,0)*(x2v - x2min_m)/(x2max_m - x2min_m)
                                     + md_d(n,1)*(x3v - x3min_m)/(x3max_m - x3min_m))
                             + md_d(n,3));
      }
      v1 = vpert*cs*sin(M_PI*(z - zmin)/lz)*amp;
    }
    u0(m,IDN,k,j,i) = d;
    u0(m,IM1,k,j,i) = d*v1;
    u0(m,IM2,k,j,i) = 0.0;
    u0(m,IM3,k,j,i) = 0.0;
    u0(m,IEN,k,j,i) = e + 0.5*d*v1*v1;
    if (have_phi) {
      const Real phi_c = g0*(z - zmin);
      phicc(m,k,j,i) = phi_c;
      if (etotgrav) u0(m,IEN,k,j,i) += d*phi_c;
      ph1(m,k,j,i) = g0*(x1l - zmin);
      if (i == n1m1) ph1(m,k,j,i+1) = g0*(x1r - zmin);
      ph2(m,k,j,i) = phi_c;
      ph3(m,k,j,i) = phi_c;
      if (j == n2m1) ph2(m,k,j+1,i) = phi_c;
      if (k == n3m1) ph3(m,k+1,j,i) = phi_c;
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvSrcs
//! \brief constant gravity along x1 -- in the well-balanced form under
//! wellbalance_dynamic -- plus the top cooling layer.

void BoxConvSrcs(Mesh *pm, Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  auto eos = pmbp->phydro->peos->eos_data;
  const bool etotgrav = pmbp->phydro->use_etotgrav;
  const bool wbdyn = pmbp->phydro->use_wellbalance_dynamic;
  const bool wbx1 = pmbp->phydro->use_wb_x1;
  const WBOption wbo = pmbp->phydro->wb_option;
  DvceArray4D<Real> phicc = pmbp->phydro->phicc0;
  DvceArray4D<Real> ph1 = pmbp->phydro->phi0.x1f;
  DvceArray5D<Real> wbq0 = pmbp->phydro->wbq0;
  const Real g0 = g0_, zlo = zlo_, dzf = dzf_;
  const Real zcool = zcool_, zmax = zmax_, tcool = tcool_;
  const int nfine = nfine_;
  auto cd_d = cd_, ce_d = ce_;
  const bool use_cache = wbx1;
  const Real cwid = (zmax > zcool) ? (zmax - zcool) : 1.0;

  par_for("boxconv_srcs", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real z = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real dzc = (x1max - x1min)/indcs.nx1;
    const Real d = w0(m,IDN,k,j,i);
    Real src = -bdt*g0*d;
    if (!etotgrav) u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
    if (wbdyn) {
      Real pl, pr, d1, d2, d3;
      if (use_cache) {
        WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      } else {
        hydro::Hydro::getWBq0(eos, wbo, WBVar::wb_pres,
            w0(m,IDN,k,j,i-1), w0(m,IDN,k,j,i), w0(m,IDN,k,j,i+1),
            w0(m,IEN,k,j,i-1), w0(m,IEN,k,j,i), w0(m,IEN,k,j,i+1),
            phicc(m,k,j,i-1), ph1(m,k,j,i), phicc(m,k,j,i), ph1(m,k,j,i+1),
            phicc(m,k,j,i+1), d1, pl, d2, pr, d3);
      }
      // Cartesian: equal face areas, so the background's own pressure drop over dz
      src = bdt*(pr - pl)/dzc;
    }
    u0(m,IM1,k,j,i) += src;
    // the cooling layer: relax the SPECIFIC internal energy toward the initial column's
    if (z > zcool) {
      Real s = (z - zcool)/cwid;
      s = (s > 1.0) ? 1.0 : s;
      const Real ramp = s*s*(3.0 - 2.0*s);
      Real t = (z - zlo)/dzf;
      int ii = static_cast<int>(t);
      ii = (ii < 0) ? 0 : ((ii > nfine-2) ? nfine-2 : ii);
      const Real f = t - ii;
      const Real d0 = cd_d(ii)*(1.0 - f) + cd_d(ii+1)*f;
      const Real e0 = ce_d(ii)*(1.0 - f) + ce_d(ii+1)*f;
      u0(m,IEN,k,j,i) -= bdt*ramp*d*(w0(m,IEN,k,j,i)/d - e0/d0)/tcool;
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvBC
//! \brief the initial column continued into the x1 ghosts with the normal velocity
//! mirrored: a reflecting wall whose ghost state is hydrostatic, so the well-balanced
//! stencil of the wall cell sees a consistent background.  The bottom wall's conductive
//! flux is REPLACED by <hydro>/rad_flux_inner, so the ghost temperature there is inert.

void BoxConvBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &u0 = pmbp->phydro->u0;
  auto &w0 = pmbp->phydro->w0;
  const Real g0 = g0_, zlo = zlo_, dzf = dzf_, zmin = zmin_;
  const int nfine = nfine_;
  const bool etotgrav = etotgrav_;
  auto cd_d = cd_, ce_d = ce_;
  const int bcm = bc_mode_;
  auto fill = KOKKOS_LAMBDA(const int m, const int k, const int j, const int i,
                            const int km, const int jm, const int im) {
    // (k,j,i) the ghost cell, (km,jm,im) the active cell it mirrors
    const Real x1min = size.d_view(m).x1min, x1max = size.d_view(m).x1max;
    const Real zg = CellCenterX(i-is, indcs.nx1, x1min, x1max);
    const Real zm = CellCenterX(im-is, indcs.nx1, x1min, x1max);
    Real d, e;
    if (bcm == 1) {
      d = w0(m,IDN,km,jm,im);
      e = w0(m,IEN,km,jm,im);
    } else {
      Real sg = (zg - zlo)/dzf;
      int ig = static_cast<int>(sg);
      ig = (ig < 0) ? 0 : ((ig > nfine-2) ? nfine-2 : ig);
      const Real fg = sg - ig;
      const Real dg = cd_d(ig)*(1.0 - fg) + cd_d(ig+1)*fg;
      const Real eg = ce_d(ig)*(1.0 - fg) + ce_d(ig+1)*fg;
      if (bcm == 0) {
        d = dg;
        e = eg;
      } else {
        Real sm = (zm - zlo)/dzf;
        int im2 = static_cast<int>(sm);
        im2 = (im2 < 0) ? 0 : ((im2 > nfine-2) ? nfine-2 : im2);
        const Real fm = sm - im2;
        const Real dm = cd_d(im2)*(1.0 - fm) + cd_d(im2+1)*fm;
        const Real em = ce_d(im2)*(1.0 - fm) + ce_d(im2+1)*fm;
        d = w0(m,IDN,km,jm,im)*(dg/dm);
        e = w0(m,IEN,km,jm,im)*(eg/em);
      }
    }
    const Real v1 = -w0(m,IVX,km,jm,im);
    const Real v2 = w0(m,IVY,km,jm,im);
    const Real v3 = w0(m,IVZ,km,jm,im);
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
    if (etotgrav) et += d*g0*(zg - zmin);
    u0(m,IEN,k,j,i) = et;
  };
  par_for("boxconv_bc_x1", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, ng-1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int n) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      fill(m, k, j, is-1-n, k, j, is+n);
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      fill(m, k, j, ie+1+n, k, j, ie-n);
    }
  });
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoxConvFinal
//! \brief releases the file-scope device views before Kokkos::finalize()

void BoxConvFinal(ParameterInput *pin, Mesh *pm) {
  cd_ = DvceArray1D<Real>();
  ce_ = DvceArray1D<Real>();
  cp_ = DvceArray1D<Real>();
  ct_ = DvceArray1D<Real>();
  return;
}
