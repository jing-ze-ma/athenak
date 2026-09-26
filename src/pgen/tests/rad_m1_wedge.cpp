//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file rad_m1_wedge.cpp
//! \brief <problem>/m1_test = sph_wedge (branch m1-wedge, docs/dev/m1_wedge_0926.md):
//! a GRAVITY-BEARING radiation-hydrodynamics wedge on the spherical-polar mesh with the
//! implicit M1 module.  Reached from RadiationM1Tests2(); every other m1_test value is
//! untouched.
//!
//! Physics.  A point mass G M = <problem>/wg_gm at the origin, the gas moving freely,
//! the luminosity entering through the implicit solver's inner face
//! (<rad_m1>/implicit_bc_x1min = flux, implicit_flux_x1min = F_in at r_in) and leaving
//! through the Marshak outer face.  F(r) = F_in (r_in/r)^2 is the radiative-equilibrium
//! flux.  The initial state is either
//!   wg_ic = grey : the grey (opacity = const) atmosphere in hydrostatic AND radiative
//!                  equilibrium, integrated by RK4 on a fine radial grid from the top
//!                  face inward (Eddington closure, E_top = F_top/(c marshak_q)):
//!                    dE/dr   = -3 rho kappa_t F/c
//!                    dp/dr   = -rho (GM/r^2 - kappa_t F/c),   p = rho T,
//!                    T = (E/a)^(1/4),  rho_top = <problem>/wg_rho_top;
//!                  an IDEAL gas, e = p/(gamma-1).  This is gate (A).
//!   wg_ic = file : <problem>/wg_ic_file, columns  r rho eint [E]  in code units,
//!                  ascending r (# comments), resampled in log onto the fine grid; any
//!                  EOS (gas only).  Without the E column E = a T(rho,eint)^4.
//!
//! Gravity (RadM1WedgeGravity, user_srcs).  The form of red_giant.cpp's RedGiantGravity:
//! plain -rho g (plus the work term unless <hydro>/etotgrav), or, with
//! <hydro>/wellbalance_dynamic + wb_x1, the background's own area-weighted pressure
//! drop.  With <problem>/wg_phi_eff = true (default when <rad_m1>/force_reference =
//! wb_arad) the x1 well-balanced pair is built with the EFFECTIVE potential
//! Phi_eff = Phi - int a_ref dr, a_ref = kappa_t F/c of the initial state, and the M1
//! coupling hands the gas only the RESIDUAL force (SetForceReference), the spherical
//! analogue of box_convection's wb_phi_eff + force_reference = wb_arad.
//!
//! x1 boundaries (RadM1WedgeBC, mesh ix1_bc = ox1_bc = user).  Closed walls: the ghost
//! cells carry the INITIAL profile at their radius scaled by the ratio of the adjacent
//! active cell to its own initial value (rho and eint separately), with the radial
//! velocity mirrored and the transverse velocity copied; so the initial state is its
//! own ghost state, and a drifting edge cell drags its ghosts with it.  The M1 ghosts
//! are a copy (inner) and dark (outer); the implicit solve imposes its own face BCs.
//!
//! <problem>/wg_spot_amp > 0 adds an isochoric Gaussian temperature perturbation
//! T -> T(1 + A exp(-d^2/(2 w^2))), E -> E (1 + ...)^4 (test C).
//!
//! History (user_hist): face luminosities of the comoving face flux f0x1 at the first
//! interior face, the middle face, the face at r_int (tau = wg_tau_int of the initial
//! state) and the Marshak face; L_in; E_rad; e_gas; and, over the interior r <= r_int,
//! mass, KE, radial KE, radial momentum and sum p dV (Mach number = sqrt(2 KE/(Gamma
//! PV))).
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
#include "outputs/outputs.hpp"
#include "utils/wb_background.hpp"
#include "rad_m1/rad_m1.hpp"
#include "rad_m1/rad_m1_closure.hpp"
#include "rad_m1/rad_m1_opacity.hpp"
#include "units/units.hpp"
#include "pgen/pgen.hpp"

namespace {
// the initial column on a uniform fine radial grid (device), read by the x1 BCs
DvceArray1D<Real> wg_rho_, wg_eint_;
Real wg_rlo_ = 0.0, wg_dr_ = 1.0;
int wg_nf_ = 0;
Real wg_gm_ = 0.0, wg_rin_ = 1.0, wg_rint_ = 0.0, wg_fin_ = 0.0;

//! log-linear interpolation on the fine grid, clamped to its end nodes
KOKKOS_INLINE_FUNCTION
Real WgLogInterp(const DvceArray1D<Real> &a, const Real rlo, const Real dr, const int nf,
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
Real WgLinInterp(const DvceArray1D<Real> &a, const Real rlo, const Real dr, const int nf,
                 const Real r) {
  Real x = (r - rlo)/dr;
  int i = static_cast<int>(floor(x));
  i = (i < 0) ? 0 : ((i > nf - 2) ? (nf - 2) : i);
  Real w = x - i;
  w = (w < 0.0) ? 0.0 : ((w > 1.0) ? 1.0 : w);
  return (1.0 - w)*a(i) + w*a(i+1);
}

//! the grey ODE right-hand side, y = (E, p)
void WgGreyRHS(const Real r, const Real *y, const Real gm, const Real kt, const Real cl,
               const Real ar, const Real fr2, Real *dy) {
  const Real f = fr2/(r*r);
  const Real t = pow(fmax(y[0], 1.0e-300)/ar, 0.25);
  const Real d = y[1]/t;
  dy[0] = -3.0*d*kt*f/cl;
  dy[1] = -d*(gm/(r*r) - kt*f/cl);
}

//! one RK4 step of the grey ODE from r to r + h
void WgGreyStep(const Real r, const Real h, Real *y, const Real gm, const Real kt,
                const Real cl, const Real ar, const Real fr2) {
  Real k1[2], k2[2], k3[2], k4[2], yt[2];
  WgGreyRHS(r, y, gm, kt, cl, ar, fr2, k1);
  for (int n=0; n<2; ++n) yt[n] = y[n] + 0.5*h*k1[n];
  WgGreyRHS(r + 0.5*h, yt, gm, kt, cl, ar, fr2, k2);
  for (int n=0; n<2; ++n) yt[n] = y[n] + 0.5*h*k2[n];
  WgGreyRHS(r + 0.5*h, yt, gm, kt, cl, ar, fr2, k3);
  for (int n=0; n<2; ++n) yt[n] = y[n] + h*k3[n];
  WgGreyRHS(r + h, yt, gm, kt, cl, ar, fr2, k4);
  for (int n=0; n<2; ++n) y[n] += h*(k1[n] + 2.0*k2[n] + 2.0*k3[n] + k4[n])/6.0;
}
//! the merged opacity table format of box_convection.cpp's ReadOpacityTable (copied:
//! comment lines, one of them "# nT nD lTmin dlT lDmin dlD", then nT*nD values of
//! log10 kappa with T slowest), filled on the host and deep-copied to the device
void WgReadOpacityTable(const std::string &fname, DvceArray2D<Real> &tab,
                        DvceArray1D<Real> &lT, DvceArray1D<Real> &lD, int &nT, int &nD) {
  std::ifstream f(fname);
  if (!f.good()) {
    std::cout << "### FATAL ERROR in sph_wedge: cannot open opacity table '" << fname
              << "'" << std::endl;
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
    std::cout << "### FATAL ERROR in sph_wedge: opacity table '" << fname
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
}
}  // namespace

void RadM1WedgeGravity(Mesh *pm, const Real bdt);
void RadM1WedgeBC(Mesh *pm);
void RadM1WedgeHist(HistoryData *pdata, Mesh *pm);
void RadM1WedgeFinal(ParameterInput *pin, Mesh *pm);

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::RadiationM1Wedge()
//! \brief m1_test = sph_wedge, see the file header.  Runs on a restart as well (the
//! hooks, the potentials and the force reference are not in the restart file) and
//! returns before writing the evolved arrays.

void ProblemGenerator::RadiationM1Wedge(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto *pm1 = pmbp->pradm1;
  auto *ph = pmbp->phydro;
  if (ph == nullptr || !pmy_mesh_->use_spherical_polar) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "<problem>/m1_test = sph_wedge needs a <hydro> block and "
      << "mesh/use_spherical_polar = true" << std::endl;
    std::exit(EXIT_FAILURE);
  }
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
  if (wbdyn && !ph->use_wb_x1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "sph_wedge: <hydro>/wellbalance_dynamic needs wb_x1 = true (the "
      << "gravity source reads the x1 well-balanced cache)" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ---- parameters
  wg_gm_ = pin->GetReal("problem","wg_gm");
  wg_rin_ = pmy_mesh_->mesh_size.x1min;
  const Real rtop = pmy_mesh_->mesh_size.x1max;
  wg_fin_ = pin->GetOrAddReal("rad_m1","implicit_flux_x1min",0.0);
  const Real fr2 = wg_fin_*wg_rin_*wg_rin_;          // F(r) = fr2/r^2
  const std::string icm = pin->GetOrAddString("problem","wg_ic","grey");
  const int nf = pin->GetOrAddInteger("problem","wg_nfine",8192);
  const Real tau_int = pin->GetOrAddReal("problem","wg_tau_int",1.0);
  const bool fref = (pm1->force_ref == radm1::M1_FREF_WB_ARAD);
  const bool phieff = pin->GetOrAddBoolean("problem","wg_phi_eff",fref);
  if (fref && !phieff) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "sph_wedge: force_reference = wb_arad needs problem/wg_phi_eff = "
      << "true (the well-balanced pair must carry the reference)" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (phieff && !wbdyn) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "sph_wedge: problem/wg_phi_eff needs <hydro>/wellbalance_dynamic "
      << "and wb_x1" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // the fine grid spans the mesh plus a margin for the ghost cells (stretched grids:
  // allow 4x the mean radial width per ghost)
  const Real dxm = (rtop - wg_rin_)/indcs.nx1;
  const Real marg = 4.0*(ng + 1)*dxm;
  wg_rlo_ = wg_rin_ - marg;
  const Real rhi = rtop + marg;
  wg_nf_ = nf;
  wg_dr_ = (rhi - wg_rlo_)/(nf - 1);
  std::vector<Real> hr(nf), hd(nf), he(nf), hE(nf, -1.0);
  for (int n=0; n<nf; ++n) hr[n] = wg_rlo_ + n*wg_dr_;

  if (icm.compare("grey") == 0) {
    if (!eos.is_ideal || pm1->opacity_type != radm1::M1_OPAC_CONST) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "sph_wedge: wg_ic = grey needs an ideal EOS and <rad_m1>/"
        << "opacity = const" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    const Real kt = pm1->kappa_f + pm1->kappa_s;
    const Real rhot = pin->GetReal("problem","wg_rho_top");
    const Real gm1 = eos.gamma - 1.0;
    const Real qm = pm1->marshak_q;
    Real y0[2];
    y0[0] = fr2/(rtop*rtop)/(cl*qm);
    y0[1] = rhot*pow(y0[0]/ar, 0.25);
    // node nt: the last node at or below the top face
    int nt = static_cast<int>(floor((rtop - wg_rlo_)/wg_dr_));
    nt = std::min(std::max(nt, 0), nf - 2);
    // inward: a partial step to node nt, then full steps
    Real y[2] = {y0[0], y0[1]};
    Real r = rtop;
    for (int n=nt; n>=0; --n) {
      const Real h = hr[n] - r;
      const int nsub = 4;
      for (int s=0; s<nsub; ++s) {
        WgGreyStep(r, h/nsub, y, wg_gm_, kt, cl, ar, fr2);
        r += h/nsub;
      }
      r = hr[n];
      const Real t = pow(y[0]/ar, 0.25);
      hd[n] = y[1]/t; he[n] = y[1]/gm1; hE[n] = y[0];
    }
    // outward from the top face
    y[0] = y0[0]; y[1] = y0[1]; r = rtop;
    for (int n=nt+1; n<nf; ++n) {
      const Real h = hr[n] - r;
      const int nsub = 4;
      for (int s=0; s<nsub; ++s) {
        WgGreyStep(r, h/nsub, y, wg_gm_, kt, cl, ar, fr2);
        r += h/nsub;
      }
      r = hr[n];
      const Real t = pow(y[0]/ar, 0.25);
      hd[n] = y[1]/t; he[n] = y[1]/gm1; hE[n] = y[0];
    }
  } else if (icm.compare("file") == 0) {
    const std::string fn = pin->GetString("problem","wg_ic_file");
    std::ifstream f(fn);
    if (!f.good()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "sph_wedge: cannot open problem/wg_ic_file '" << fn << "'"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::vector<Real> fr, fd, fe, fE;
    std::string line;
    bool hasE = true;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      Real a, b, c, d;
      if (!(ss >> a >> b >> c)) continue;
      if (!(ss >> d)) {
        hasE = false;
        d = -1.0;
      }
      fr.push_back(a); fd.push_back(b); fe.push_back(c); fE.push_back(d);
    }
    if (fr.size() < 2 || fr.front() > hr.front() || fr.back() < hr.back()) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "sph_wedge: problem/wg_ic_file '" << fn << "' must cover r = "
        << hr.front() << " .. " << hr.back() << " (the mesh plus the ghost margin)"
        << std::endl;
      std::exit(EXIT_FAILURE);
    }
    std::size_t kk = 0;
    for (int n=0; n<nf; ++n) {
      const Real r = hr[n];
      while (kk + 2 < fr.size() && fr[kk+1] < r) ++kk;
      const Real w = (r - fr[kk])/(fr[kk+1] - fr[kk]);
      hd[n] = exp((1.0 - w)*log(fd[kk]) + w*log(fd[kk+1]));
      he[n] = exp((1.0 - w)*log(fe[kk]) + w*log(fe[kk+1]));
      hE[n] = hasE ? exp((1.0 - w)*log(fE[kk]) + w*log(fE[kk+1])) : -1.0;
    }
  } else {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
      << std::endl << "sph_wedge: problem/wg_ic = '" << icm << "' (grey | file)"
      << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ---- <rad_m1>/opacity = table: the pgen hands the Rosseland + Planck tables over
  // (problem/wg_opac_table, wg_planck_table; one shared grid), as box_convection does
  if (pm1->opacity_type == radm1::M1_OPAC_TABLE) {
    const std::string rt = pin->GetString("problem","wg_opac_table");
    const std::string pt = pin->GetString("problem","wg_planck_table");
    DvceArray2D<Real> krt, kpt;
    DvceArray1D<Real> mlT, mlD, plT, plD;
    int mnT = 0, mnD = 0, pnT = 0, pnD = 0;
    WgReadOpacityTable(rt, krt, mlT, mlD, mnT, mnD);
    WgReadOpacityTable(pt, kpt, plT, plD, pnT, pnD);
    if (mnT != pnT || mnD != pnD) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "sph_wedge: the Rosseland and Planck tables are not on the "
        << "same grid" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    pm1->SetOpacityTables(krt, kpt, mlT, mlD, mnT, mnD);
    // the lookup's units against the EOS's own code temperature and <units> density
    const Real tcgs = eos.temp_cgs;
    const Real dcgs = (pmbp->punit != nullptr) ? pmbp->punit->density_cgs() : 1.0;
    if (std::fabs(pm1->otab.tunit/tcgs - 1.0) > 1.0e-5 ||
        std::fabs(pm1->otab.dunit/dcgs - 1.0) > 1.0e-5) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "sph_wedge: <rad_m1>/temp_unit_kelvin = " << pm1->otab.tunit
        << ", rho_unit_cgs = " << pm1->otab.dunit << " but this run's units are ("
        << tcgs << ", " << dcgs << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // ---- the column on the device: T, E (if not given), kappa_t and a_ref = kt F/c
  Kokkos::realloc(wg_rho_, nf);
  Kokkos::realloc(wg_eint_, nf);
  DvceArray1D<Real> cE("wg_E", nf), cT("wg_T", nf), ckt("wg_kt", nf);
  {
    auto h1 = Kokkos::create_mirror_view(wg_rho_);
    auto h2 = Kokkos::create_mirror_view(wg_eint_);
    auto h3 = Kokkos::create_mirror_view(cE);
    for (int n=0; n<nf; ++n) {
      h1(n) = hd[n];
      h2(n) = he[n];
      h3(n) = hE[n];
    }
    Kokkos::deep_copy(wg_rho_, h1);
    Kokkos::deep_copy(wg_eint_, h2);
    Kokkos::deep_copy(cE, h3);
  }
  {
    auto crho = wg_rho_, ceint = wg_eint_;
    const int otype = pm1->opacity_type;
    const Real kp = pm1->kappa_p, ke = pm1->kappa_e, kf = pm1->kappa_f;
    const Real ks = pm1->kappa_s;
    const Real rref = pm1->opac_rho_ref, tref = pm1->opac_t_ref;
    const Real aa = pm1->opac_a, bb = pm1->opac_b;
    radm1::M1OpacTab ot = pm1->otab;
    par_for("wg_col", DevExeSpace(), 0, nf-1, KOKKOS_LAMBDA(const int n) {
      const Real d = crho(n);
      const Real t = eos.Temperature(d, ceint(n));
      Real op, oe, of, os;
      if (otype == radm1::M1_OPAC_TABLE) {
        radm1::M1TableOpacities(ot, d, t, op, oe, of, os);
      } else {
        radm1::M1Opacities(otype, d, t, kp, ke, kf, ks, rref, tref, aa, bb,
                           op, oe, of, os);
      }
      cT(n) = t;
      ckt(n) = of + os;
      if (cE(n) < 0.0) cE(n) = ar*t*t*t*t;
    });
  }
  auto hT = Kokkos::create_mirror_view(cT);
  auto hkt = Kokkos::create_mirror_view(ckt);
  auto hEd = Kokkos::create_mirror_view(cE);
  Kokkos::deep_copy(hT, cT);
  Kokkos::deep_copy(hkt, ckt);
  Kokkos::deep_copy(hEd, cE);
  // tau from the top face, r_int, Phi_eff = GM(1/r_in - 1/r) - int_{r_in}^r a_ref dr
  std::vector<Real> tau(nf, 0.0), aref(nf), cum(nf, 0.0);
  for (int n=0; n<nf; ++n) aref[n] = hkt(n)*fr2/(hr[n]*hr[n])/cl;
  for (int n=nf-2; n>=0; --n) {
    const Real a0 = (hr[n] < rtop) ? hd[n]*hkt(n) : 0.0;
    const Real a1 = (hr[n+1] < rtop) ? hd[n+1]*hkt(n+1) : 0.0;
    tau[n] = tau[n+1] + 0.5*wg_dr_*(a0 + a1);
  }
  wg_rint_ = wg_rin_;
  for (int n=nf-1; n>=0; --n) {
    if (tau[n] >= tau_int) {
      wg_rint_ = std::max(hr[n], wg_rin_);
      break;
    }
  }
  for (int n=1; n<nf; ++n) cum[n] = cum[n-1] + 0.5*wg_dr_*(aref[n-1] + aref[n]);
  // shift the running integral so that it is zero at r_in
  Real cin = 0.0;
  {
    const Real x = (wg_rin_ - wg_rlo_)/wg_dr_;
    int i = std::min(std::max(static_cast<int>(floor(x)), 0), nf - 2);
    const Real w = x - i;
    cin = (1.0 - w)*cum[i] + w*cum[i+1];
  }
  DvceArray1D<Real> cphi("wg_phieff", nf);
  {
    auto hp = Kokkos::create_mirror_view(cphi);
    for (int n=0; n<nf; ++n) {
      hp(n) = wg_gm_*(1.0/wg_rin_ - 1.0/hr[n]) - (cum[n] - cin);
    }
    Kokkos::deep_copy(cphi, hp);
  }

  // ---- potentials: the TRUE one (etotgrav, the WB default) ...
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  const Real gm = wg_gm_, rin = wg_rin_;
  if (etg || wbdyn) {
    auto phicc = ph->phicc0;
    auto ph1 = ph->phi0.x1f, ph2 = ph->phi0.x2f, ph3 = ph->phi0.x3f;
    par_for("wg_phi", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
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
  const Real rlo = wg_rlo_, dr = wg_dr_;
  if (phieff) {
    ph->EnableWBEffectivePotential();
    auto pwc = ph->phicc_wb, pwf = ph->phi_wb_x1f;
    par_for("wg_phieff", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      pwc(m,k,j,i) = WgLinInterp(cphi, rlo, dr, nf, x1v(m,i));
      pwf(m,k,j,i) = WgLinInterp(cphi, rlo, dr, nf, x1f(m,i));
      if (i == n1m1) pwf(m,k,j,i+1) = WgLinInterp(cphi, rlo, dr, nf, x1f(m,i+1));
    });
  }
  // the reference acceleration, face form: kt(cell) (F(r_l) + F(r_r))/(2c), from the
  // INITIAL column (identical on a restart)
  if (fref) {
    DvceArray4D<Real> aref_d("wg_aref", nmb1+1, n3m1+1, n2m1+1, n1m1+1);
    par_for("wg_aref", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real kt = WgLinInterp(ckt, rlo, dr, nf, x1v(m,i));
      const Real rl = x1f(m,i), rr = x1f(m,i+1);
      aref_d(m,k,j,i) = 0.5*kt*fr2*(1.0/(rl*rl) + 1.0/(rr*rr))/cl;
    });
    pm1->SetForceReference(aref_d);
  }

  user_srcs_func = RadM1WedgeGravity;
  user_bcs_func = RadM1WedgeBC;
  user_hist_func = RadM1WedgeHist;
  pgen_final_func = RadM1WedgeFinal;
  if (global_variable::my_rank == 0) {
    Real a0 = aref[0], g0 = wg_gm_/(hr[0]*hr[0]);
    std::cout << "sph_wedge: wg_ic = " << icm << ", GM = " << wg_gm_ << ", F_in = "
              << wg_fin_ << ", Gamma(r_lo) = kt F/(c g) = " << a0/g0
              << ", tau(r_in) = " << tau[static_cast<int>((wg_rin_ - wg_rlo_)/wg_dr_)]
              << ", r_int (tau = " << tau_int << ") = " << wg_rint_
              << ", phi_eff = " << (phieff ? "on" : "off")
              << ", force_reference = " << (fref ? "wb_arad" : "none") << std::endl;
  }
  if (restart) return;

  // ---- the initial state (active + ghost cells)
  const Real samp = pin->GetOrAddReal("problem","wg_spot_amp",0.0);
  const Real sr = pin->GetOrAddReal("problem","wg_spot_r",0.5*(wg_rin_ + rtop));
  const Real sth = pin->GetOrAddReal("problem","wg_spot_th",0.5*M_PI);
  const Real sph = pin->GetOrAddReal("problem","wg_spot_ph",
                     0.5*(pmy_mesh_->mesh_size.x3min + pmy_mesh_->mesh_size.x3max));
  const Real sw = pin->GetOrAddReal("problem","wg_spot_w",0.05*(rtop - wg_rin_));
  // wg_seed: a small deterministic cell-to-cell perturbation of the gas internal energy
  // (a seed for convection, as sph_atm's atm_seed), applied where tau(IC) >= wg_tau_int;
  // a function of the GLOBAL cell indices, so independent of the rank decomposition
  const Real seed = pin->GetOrAddReal("problem","wg_seed",0.0);
  const Real seed_rmax = (seed != 0.0) ? wg_rint_ : 0.0;
  const Real x2a = pmy_mesh_->mesh_size.x2min, x3a = pmy_mesh_->mesh_size.x3min;
  const Real idx2 = pmy_mesh_->mesh_indcs.nx2/(pmy_mesh_->mesh_size.x2max - x2a);
  const Real idx3 = pmy_mesh_->mesh_indcs.nx3/(pmy_mesh_->mesh_size.x3max - x3a);
  const Real sx = sr*sin(sth)*cos(sph), sy = sr*sin(sth)*sin(sph), sz = sr*cos(sth);
  auto uh = ph->u0;
  auto ur = pm1->u0;
  auto x2v = pmbp->pcoord->x2v;
  auto x3v = pmbp->pcoord->x3v;
  auto crho = wg_rho_, ceint = wg_eint_;
  const bool ideal = eos.is_ideal;
  par_for("wg_ic", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1, 0, n1m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const Real d = WgLogInterp(crho, rlo, dr, nf, r);
    Real e = WgLogInterp(ceint, rlo, dr, nf, r);
    Real er = WgLogInterp(cE, rlo, dr, nf, r);
    if (samp != 0.0) {
      const Real th = x2v(m,j), ph_ = x3v(m,k);
      const Real x = r*sin(th)*cos(ph_), y = r*sin(th)*sin(ph_), z = r*cos(th);
      const Real d2 = SQR(x - sx) + SQR(y - sy) + SQR(z - sz);
      const Real del = samp*exp(-0.5*d2/(sw*sw));
      // isochoric: T -> T (1 + del); for a general EOS e(rho,T) is not linear in T,
      // so the gas perturbation is the ideal-gas one (a test knob only)
      e *= (1.0 + del);
      er *= SQR(SQR(1.0 + del));
      (void) ideal;
    }
    if (seed != 0.0 && r <= seed_rmax) {
      const int jg = static_cast<int>(floor((x2v(m,j) - x2a)*idx2));
      const int kg = static_cast<int>(floor((x3v(m,k) - x3a)*idx3));
      e *= 1.0 + seed*cos(2.3*jg + 1.7*kg + 0.9*i);
    }
    uh(m,IDN,k,j,i) = d;
    uh(m,IM1,k,j,i) = 0.0;
    uh(m,IM2,k,j,i) = 0.0;
    uh(m,IM3,k,j,i) = 0.0;
    uh(m,IEN,k,j,i) = e + (etg ? d*gm*(1.0/rin - 1.0/r) : 0.0);
    ur(m,radm1::M1_E,k,j,i) = fmax(er, efl);
    ur(m,radm1::M1_F1,k,j,i) = fr2/(r*r);
    ur(m,radm1::M1_F2,k,j,i) = 0.0;
    ur(m,radm1::M1_F3,k,j,i) = 0.0;
  });
  // the comoving x1 face fluxes of the implicit transport (allocated only there)
  auto ff0 = pm1->f0x1;
  if (ff0.extent_int(0) >= nmb1 + 1) {
    const int e3 = ff0.extent_int(1) - 1, e2 = ff0.extent_int(2) - 1;
    const int e1 = ff0.extent_int(3) - 1;
    par_for("wg_icf", DevExeSpace(), 0, nmb1, 0, e3, 0, e2, 0, e1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      const Real r = x1f(m,i);
      ff0(m,k,j,i) = fr2/(r*r);
    });
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1WedgeGravity()
//! \brief the point-mass source, red_giant.cpp's RedGiantGravity without its extras.

void RadM1WedgeGravity(Mesh *pm, const Real bdt) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, js = indcs.js, je = indcs.je;
  const int ks = indcs.ks, ke = indcs.ke;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto *ph = pmbp->phydro;
  auto u0 = ph->u0;
  auto w0 = ph->w0;
  auto eos = ph->peos->eos_data;
  const bool etg = ph->use_etotgrav;
  const bool wbdyn = ph->use_wellbalance_dynamic;
  auto wbq0 = ph->wbq0;
  auto &x1v = pmbp->pcoord->x1v;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  const Real gm = wg_gm_;
  par_for("wg_grav", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    const Real r = x1v(m,i);
    const Real d = w0(m,IDN,k,j,i);
    Real src = -bdt*d*gm/(r*r);
    if (!etg) u0(m,IEN,k,j,i) += src*w0(m,IVX,k,j,i);
    if (wbdyn) {
      Real pl, pr, d1, d2, d3;
      WBReadCache(wbq0, WBVar::wb_pres, m, k, j, i, d1, pl, d2, pr, d3);
      const Real e_ = w0(m,IEN,k,j,i);
      const Real p = (e_ > 0.0) ? eos.Pressure(d, e_) : 0.5*(pl + pr);
      src = bdt*(area1(m,k,j,i+1)*(pr - p) + area1(m,k,j,i)*(p - pl))/volume(m,k,j,i);
    }
    u0(m,IM1,k,j,i) += src;
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1WedgeBC()
//! \brief closed x1 walls carrying the scaled initial profile (see the file header), for
//! the hydro conserved variables and the M1 ghosts.  Called from both modules' BC tasks;
//! both fills are functions of the active cells only.

void RadM1WedgeBC(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is, ie = indcs.ie;
  const int n2m1 = (indcs.nx2 > 1) ? (indcs.nx2 + 2*ng - 1) : 0;
  const int n3m1 = (indcs.nx3 > 1) ? (indcs.nx3 + 2*ng - 1) : 0;
  const int nmb1 = pmbp->nmb_thispack - 1;
  auto &mbbcs = pmbp->pmb->mb_bcs;
  auto &x1v = pmbp->pcoord->x1v;
  auto *ph = pmbp->phydro;
  const bool have_phi = ph->use_etotgrav || ph->use_wellbalance_dynamic;
  const bool etg = ph->use_etotgrav;
  auto uh = ph->u0;
  auto phicc = ph->phicc0;
  auto crho = wg_rho_, ceint = wg_eint_;
  const Real rlo = wg_rlo_, dr = wg_dr_;
  const int nf = wg_nf_;
  const bool rad = (pmbp->pradm1 != nullptr);
  auto ur = rad ? pmbp->pradm1->u0 : DvceArray5D<Real>();
  const Real cl = rad ? pmbp->pradm1->c_light : 1.0;
  const Real efl = rad ? pmbp->pradm1->e_floor : 0.0;
  (void) have_phi;
  par_for("wg_bc", DevExeSpace(), 0, nmb1, 0, n3m1, 0, n2m1,
  KOKKOS_LAMBDA(const int m, const int k, const int j) {
    for (int side=0; side<2; ++side) {
      const bool lo = (side == 0);
      if (mbbcs.d_view(m, lo ? BoundaryFace::inner_x1 : BoundaryFace::outer_x1)
          != BoundaryFlag::user) continue;
      const int ia = lo ? is : ie;
      const Real da = uh(m,IDN,k,j,ia);
      const Real kea = 0.5*(SQR(uh(m,IM1,k,j,ia)) + SQR(uh(m,IM2,k,j,ia))
                            + SQR(uh(m,IM3,k,j,ia)))/da;
      const Real ea = uh(m,IEN,k,j,ia) - kea - (etg ? da*phicc(m,k,j,ia) : 0.0);
      const Real sd = da/WgLogInterp(crho, rlo, dr, nf, x1v(m,ia));
      const Real se = ea/WgLogInterp(ceint, rlo, dr, nf, x1v(m,ia));
      for (int g=0; g<ng; ++g) {
        const int ig = lo ? (is - 1 - g) : (ie + 1 + g);
        const int im = lo ? (is + g) : (ie - g);        // the mirror cell
        const Real rg = x1v(m,ig);
        const Real dg = sd*WgLogInterp(crho, rlo, dr, nf, rg);
        const Real eg = se*WgLogInterp(ceint, rlo, dr, nf, rg);
        const Real dm = uh(m,IDN,k,j,im);
        const Real v1 = -uh(m,IM1,k,j,im)/dm;
        const Real v2 = uh(m,IM2,k,j,im)/dm;
        const Real v3 = uh(m,IM3,k,j,im)/dm;
        uh(m,IDN,k,j,ig) = dg;
        uh(m,IM1,k,j,ig) = dg*v1;
        uh(m,IM2,k,j,ig) = dg*v2;
        uh(m,IM3,k,j,ig) = dg*v3;
        uh(m,IEN,k,j,ig) = eg + 0.5*dg*(v1*v1 + v2*v2 + v3*v3)
                           + (etg ? dg*phicc(m,k,j,ig) : 0.0);
        if (rad) {
          radm1::M1FillGhost(ur, m, k, j, ig, k, j, ia, 1, lo ? 0 : 2,
                             lo ? -1.0 : 1.0, cl, efl);
        }
      }
    }
  });
}

//----------------------------------------------------------------------------------------
//! \fn void RadM1WedgeHist()
//! \brief the sph_wedge history columns (see the file header), all plain sums.

void RadM1WedgeHist(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  pdata->nhist = 13;
  pdata->label[0] = "L_bot";
  pdata->label[1] = "L_mid";
  pdata->label[2] = "L_int";
  pdata->label[3] = "L_top";
  pdata->label[4] = "L_in";
  pdata->label[5] = "E_rad";
  pdata->label[6] = "e_gas";
  pdata->label[7] = "M_int";
  pdata->label[8] = "KE_int";
  pdata->label[9] = "KEr_int";
  pdata->label[10] = "Mr_int";
  pdata->label[11] = "PV_int";
  pdata->label[12] = "V_int";
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is, ie = indcs.ie, nx1 = indcs.nx1;
  const int js = indcs.js, nx2 = indcs.nx2, ks = indcs.ks, nx3 = indcs.nx3;
  const int imid = is + nx1/2;
  const int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1, nji = nx2*nx1;
  auto *ph = pmbp->phydro;
  auto w0 = ph->w0;
  auto eos = ph->peos->eos_data;
  auto *pm1 = pmbp->pradm1;
  auto ur = pm1->u0;
  auto ff0 = pm1->f0x1;
  const bool haveff = (ff0.extent_int(0) > 0);
  auto &x1v = pmbp->pcoord->x1v;
  auto &x1f = pmbp->pcoord->xx1f;
  auto &area1 = pmbp->pcoord->area.x1f;
  auto &volume = pmbp->pcoord->volume;
  const Real rint = wg_rint_, fin = wg_fin_;
  array_sum::GlobalSum sum_this;
  Kokkos::parallel_reduce("wg_hist", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
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
//! \fn void RadM1WedgeFinal()
//! \brief release the namespace-scope Views before Kokkos::finalize.

void RadM1WedgeFinal(ParameterInput *pin, Mesh *pm) {
  (void) pin; (void) pm;
  wg_rho_ = DvceArray1D<Real>();
  wg_eint_ = DvceArray1D<Real>();
}
