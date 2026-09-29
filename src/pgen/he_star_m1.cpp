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
//!     gas + radiation, the mass through both x1 faces, and the interior kinetic sums.
//!
//! RESTARTS: the pgen carries NO state that is not recomputed here.  It is not skipped on
//! a restart: the column, the tables, the potentials, the reference acceleration and the
//! hooks are rebuilt from the input file (the evolved arrays come from the restart file).
//!
//! CUDA-safe: no lambdas inside kernels, no host reads of device Views (the fine column
//! is built on the host and copied), no class members inside kernels.

#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
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
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "rad_m1/m1_fluid.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"

namespace {
// the initial column on a uniform fine radial grid (device), read by the x1 BCs
DvceArray1D<Real> hs_rho_, hs_eint_;
Real hs_rlo_ = 0.0, hs_dr_ = 1.0;
int hs_nf_ = 0;
Real hs_gm_ = 0.0, hs_rin_ = 1.0, hs_rint_ = 0.0, hs_fin_ = 0.0;
Real hs_sp_rate_ = 0.0, hs_sp_r0_ = 0.0, hs_rtop_ = 1.0;
bool hs_zflux_ = true;

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
                        Real &lt_lo, Real &lt_hi, Real &ld_lo, Real &ld_hi) {
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
void HeStarFinal(ParameterInput *pin, Mesh *pm);
}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief see the file header.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto *pm1 = pmbp->pradm1;
  if (pm1 == nullptr || pmbp->phydro == nullptr || !pmy_mesh_->use_spherical_polar) {
    HsFatal("needs <hydro>, <rad_m1> and mesh/use_spherical_polar = true", __LINE__);
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
  hs_rin_ = pmy_mesh_->mesh_size.x1min;
  const Real rtop = pmy_mesh_->mesh_size.x1max;
  hs_fin_ = pin->GetReal("rad_m1","implicit_flux_x1min");
  const std::string fn = pin->GetString("problem","he_ic_file");
  const int nf = pin->GetOrAddInteger("problem","he_nfine",16384);
  const Real tau_int = pin->GetOrAddReal("problem","he_tau_int",1.0);
  // the fine grid spans the mesh plus a margin of (ng + 2) mean radial widths
  const Real dxm = (rtop - hs_rin_)/indcs.nx1;
  const Real marg = (ng + 2)*dxm;
  hs_rlo_ = hs_rin_ - marg;
  const Real rhi = rtop + marg;
  hs_nf_ = nf;
  hs_dr_ = (rhi - hs_rlo_)/(nf - 1);
  std::vector<Real> hr(nf), hd(nf), he(nf), hF(nf);
  for (int n=0; n<nf; ++n) hr[n] = hs_rlo_ + n*hs_dr_;

  // ---- the column file: r rho eint F_r, ascending r, resampled in log
  {
    std::ifstream f(fn);
    if (!f.good()) HsFatal("cannot open problem/he_ic_file '" + fn + "'", __LINE__);
    std::vector<Real> fr, fd, fe, fF;
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      Real a, b, c, d;
      if (!(ss >> a >> b >> c >> d)) continue;
      fr.push_back(a); fd.push_back(b); fe.push_back(c); fF.push_back(d);
    }
    if (fr.size() < 2 || fr.front() > hr.front() || fr.back() < hr.back()) {
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
    HsReadOpacityTable(rt, krt, mlT, mlD, mnT, mnD, tab_lt_lo, tab_lt_hi, tab_ld_lo,
                       tab_ld_hi);
    HsReadOpacityTable(pt, kpt, plT, plD, pnT, pnD, a1, a2, a3, a4);
    b1 = tab_lt_lo; b2 = tab_lt_hi; b3 = tab_ld_lo; b4 = tab_ld_hi;
    if (mnT != pnT || mnD != pnD || fabs(a1 - b1) + fabs(a2 - b2) + fabs(a3 - b3) +
        fabs(a4 - b4) > 1.0e-9) {
      HsFatal("the Rosseland and Planck tables are not on the same grid", __LINE__);
    }
    pm1->SetOpacityTables(krt, kpt, mlT, mlD, mnT, mnD);
    // the lookup's units against the EOS's own code temperature and <units> density
    const Real tcgs = eos.temp_cgs;
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
    for (int n=0; n<nf; ++n) {
      h1(n) = hd[n];
      h2(n) = he[n];
      h4(n) = hF[n];
    }
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
      cE(n) = ar*t*t*t*t;
    });
  }
  auto hT = Kokkos::create_mirror_view(cT);
  auto hkt = Kokkos::create_mirror_view(ckt);
  Kokkos::deep_copy(hT, cT);
  Kokkos::deep_copy(hkt, ckt);

  // ---- startup check (fatal): the column, ghost margin included, is inside the EOS
  // table and the opacity table grid, in (rho [g/cm^3], T [K])
  {
    const Real tk = pm1->otab.tunit, dk = pm1->otab.dunit;
    const Real elo_d = pin->GetReal("hydro","eos_logd_min");
    const Real ehi_d = pin->GetReal("hydro","eos_logd_max");
    const Real elo_t = pin->GetReal("hydro","eos_logt_min");
    const Real ehi_t = pin->GetReal("hydro","eos_logt_max");
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
  DvceArray1D<Real> cphi("hs_phieff", nf);
  {
    auto hp = Kokkos::create_mirror_view(cphi);
    for (int n=0; n<nf; ++n) {
      hp(n) = hs_gm_*(1.0/hs_rin_ - 1.0/hr[n]) - (cum[n] - cin);
    }
    Kokkos::deep_copy(cphi, hp);
  }

  // ---- potentials: the TRUE one (etotgrav, the WB default) ...
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  const Real gm = hs_gm_, rin = hs_rin_;
  {
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

  // top sponge (default off) and hooks
  hs_zflux_ = pin->GetOrAddBoolean("problem","he_wall_zero_flux",true);
  hs_sp_rate_ = pin->GetOrAddReal("problem","he_sponge_rate",0.0);
  hs_sp_r0_ = pin->GetOrAddReal("problem","he_sponge_r0",hs_rint_);
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
  const Real x2a = pmy_mesh_->mesh_size.x2min, x3a = pmy_mesh_->mesh_size.x3min;
  const Real lth = pmy_mesh_->mesh_size.x2max - x2a;
  const Real lph = pmy_mesh_->mesh_size.x3max - x3a;
  auto uh = ph->u0;
  auto ur = pm1->u0;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto crho = hs_rho_, ceint = hs_eint_;
  par_for("hs_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const Real d = HsLogInterp(crho, rlo, dr, nf, r);
    Real e = HsLogInterp(ceint, rlo, dr, nf, r);
    // seed (default off): one smooth lateral mode in eint, wave number he_seed_k
    if (seed != 0.0 && r >= srlo && r <= srhi) {
      e *= 1.0 + seed*sin(2.0*M_PI*seedk*(x2v(m,j) - x2a)/lth + 0.3)
                     *sin(2.0*M_PI*seedk*(x3v(m,k) - x3a)/lph + 1.1);
    }
    const Real t = eos.Temperature(d, HsLogInterp(ceint, rlo, dr, nf, r));
    const Real er = ar*t*t*t*t;
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = e + d*gm*(1.0/rin - 1.0/r);
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
    par_for("hs_grav_fws", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      u0(m,IEN,k,j,i) += bdt*w0(m,IDN,k,j,i)*aref(m,k,j,i)*w0(m,IVX,k,j,i);
    });
  }
  // x1 walls: the scaled-profile ghosts are not a mirror image of the edge cell, so the
  // Riemann flux through a "closed" wall carries mass.  The inner wall loses its mass,
  // transverse momentum and energy fluxes (u0 was just updated with them); the pressure
  // flux stays.  The outer face keeps OUTFLOW and loses only inflow (no mass enters).
  if (hs_zflux_) {
    auto flx1 = ph->uflx->x1f;
    auto &mbbcs = pmbp->pmb->mb_bcs;
    par_for("hs_zflux", DevExeSpace(), 0, nmb1, ks, ke, js, je,
    KOKKOS_LAMBDA(const int m, const int k, const int j) {
      if (mbbcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user) {
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
  // energy (the heat is not kept, so that the sponge adds no buoyancy)
  auto &x1v = pmbp->pcoord->x1v;
  const Real rate = hs_sp_rate_, r0 = hs_sp_r0_, rt = hs_rtop_;
  par_for("hs_sponge", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    if (r <= r0) return;
    const Real w = SQR((r - r0)/(rt - r0));
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
  auto ur = pmbp->pradm1->u0;
  const Real cl = pmbp->pradm1->c_light;
  const Real efl = pmbp->pradm1->e_floor;
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
      const Real sd = da/HsLogInterp(crho, rlo, dr, nf, x1v(m,ia));
      const Real se = ea/HsLogInterp(ceint, rlo, dr, nf, x1v(m,ia));
      for (int g=0; g<ng; ++g) {
        const int ig = lo ? (is - 1 - g) : (ie + 1 + g);
        const int im = lo ? (is + g) : (ie - g);        // the mirror cell
        const Real rg = x1v(m,ig);
        const Real dg = sd*HsLogInterp(crho, rlo, dr, nf, rg);
        const Real eg = se*HsLogInterp(ceint, rlo, dr, nf, rg);
        const Real dm = uh(m,IDN,k,j,im);
        // wall: v1 mirrored; outer edge: v1 of the edge cell where it flows out
        const Real v1e = uh(m,IM1,k,j,ia)/da;
        const Real v1 = (!lo && v1e > 0.0) ? v1e : -uh(m,IM1,k,j,im)/dm;
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
  pdata->nhist = 18;
  const char *lab[18] = {"L_bot", "L_mid", "L_int", "L_top", "L_in", "E_rad", "e_gas",
                         "M_int", "KE_int", "KEr_int", "Mr_int", "PV_int", "V_int",
                         "M_tot", "Mdot_top", "Mdot_bot", "Etot", "Min_top"};
  for (int n=0; n<18; ++n) pdata->label[n] = lab[n];
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
    if (i == is) h.the_array[15] = fl1(m,IDN,k,j,i)*area1(m,k,j,i);
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
  for (int n=0; n<pdata->nhist; ++n) pdata->hdata[n] = sum_this.the_array[n];
}

//----------------------------------------------------------------------------------------
//! \fn void HeStarFinal()
//! \brief release the namespace-scope Views before Kokkos::finalize.

void HeStarFinal(ParameterInput *pin, Mesh *pm) {
  (void) pin; (void) pm;
  hs_rho_ = DvceArray1D<Real>();
  hs_eint_ = DvceArray1D<Real>();
}
}  // namespace
