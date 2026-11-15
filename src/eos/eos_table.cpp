//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_table.cpp
//! \brief builds the (log10 rho, log10 T) table of the general EOS on the host by
//! sampling the analytic composition model (eos_composition.hpp).  Node derivatives come
//! from central differences of the model, not of the table, so the bicubic Hermite
//! interpolant is fourth-order accurate.

#include <math.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "eos/eos_composition.hpp"
#include "eos/eos_table.hpp"
#include "utils/rad_taper.hpp"

namespace {

//----------------------------------------------------------------------------------------
//! \fn void SampleModel
//! \brief the three tabulated surfaces at one (x,y) = (log10 rho, log10 T), in cgs

void SampleModel(const EOSCompositionModel &m, const double x, const double y,
                 double &fe, double &fp, double &fmu, double &fxe) {
  EOSCompositionState s = m.Evaluate(pow(10.0, x), pow(10.0, y));
  fe = log10(s.e_spec);
  fp = log10(s.p_spec);
  fmu = s.mu;
  // log10 of the electron fraction: it spans ~20 decades, and the floor keeps the
  // logarithm finite where ionization is switched off entirely
  fxe = log10(s.xe > 1.0e-300 ? s.xe : 1.0e-300);
}

//----------------------------------------------------------------------------------------
//! \fn double EnergyAtPressure
//! \brief internal energy density (gas + radiation) at a density and total pressure in
//! cgs, by bisection; used for the pressure-floor bound

double EnergyAtPressure(const EOSCompositionModel &m, const bool rad,
                        const double wrad, const double rho, const double ptarget) {
  double ylo = -2.0, yhi = 12.0;
  for (int it=0; it<200; ++it) {
    double y = 0.5*(ylo + yhi);
    double t = pow(10.0, y);
    EOSCompositionState s = m.Evaluate(rho, t);
    double p = rho*s.p_spec + (rad ? wrad*eos_cgs::a_rad*t*t*t*t/3.0 : 0.0);
    if (p > ptarget) {
      yhi = y;
    } else {
      ylo = y;
    }
    if ((yhi - ylo) < 1.0e-13) break;
  }
  double t = pow(10.0, 0.5*(ylo + yhi));
  EOSCompositionState s = m.Evaluate(rho, t);
  return rho*s.e_spec + (rad ? wrad*eos_cgs::a_rad*t*t*t*t : 0.0);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void BuildEOSTable
//! \brief reads the table and composition parameters, samples the model onto the grid and
//! copies the table to the device

void BuildEOSTable(EOSTable &tbl, ParameterInput *pin, const std::string &block,
                   const Real dens_cgs, const Real pres_cgs, const Real temp_cgs,
                   const Real pfloor) {
  // ------------------------------------------------------------------ composition model
  EOSCompositionModel model;
  model.xhyd = pin->GetOrAddReal(block, "eos_xh", 0.7381);
  model.yhel = pin->GetOrAddReal(block, "eos_yhe", 0.2485);
  model.a_metal = pin->GetOrAddReal(block, "eos_a_metal", 16.0);
  model.include_h2 = pin->GetOrAddBoolean(block, "eos_h2", true);
  model.include_ion = pin->GetOrAddBoolean(block, "eos_ionization", true);
  model.include_metal_ion = pin->GetOrAddBoolean(block, "eos_metal_ionization", false);
  // [M/H] of the metal donors, default <problem>/met when present so the EOS and the
  // opacity use the same metallicity
  Real mh_default = 0.0;
  if (pin->DoesParameterExist("problem", "met")) {
    mh_default = pin->GetReal("problem", "met");
  }
  model.metal_mh = pin->GetOrAddReal(block, "eos_metal_mh", mh_default);
  model.include_metal_cond = pin->GetOrAddBoolean(block, "eos_metal_condensation", false);
  tbl.radiation = pin->GetOrAddBoolean(block, "eos_radiation", false);

  // utils/rad_taper.hpp; eos_rad_rho_hi = 0 (default) = no taper
  const Real rad_rho_hi = pin->GetOrAddReal(block, "eos_rad_rho_hi", 0.0);
  const Real rad_rho_lo = pin->GetOrAddReal(block, "eos_rad_rho_lo", 0.0);
  tbl.rad_taper = (rad_rho_hi > 0.0);
  if (tbl.rad_taper) {
    if (!(rad_rho_lo > 0.0) || !(rad_rho_lo < rad_rho_hi)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/eos_rad_rho_lo must be positive and strictly less "
                << "than eos_rad_rho_hi; got rho_lo = " << rad_rho_lo
                << ", rho_hi = " << rad_rho_hi << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!tbl.radiation) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/eos_rad_rho_hi tapers a radiation term that is not "
                << "there: set eos_radiation = true, or drop the taper." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    tbl.rad_lrho_lo = log10(static_cast<double>(rad_rho_lo));
    tbl.rad_lrho_hi = log10(static_cast<double>(rad_rho_hi));
  }

  // w = max(w_rho, w_T): cells hotter than eos_rad_t_hi keep their radiation
  // (rad_taper::WeightGated); eos_rad_t_hi = 0 (default) = gate off
  const Real rad_t_hi = pin->GetOrAddReal(block, "eos_rad_t_hi", 0.0);
  const Real rad_t_lo = pin->GetOrAddReal(block, "eos_rad_t_lo", 0.0);
  tbl.rad_tgate = (rad_t_hi > 0.0);
  if (tbl.rad_tgate) {
    if (!tbl.rad_taper) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/eos_rad_t_hi gates a taper that is not there: the "
                << "gate raises w = max(w_rho, w_T), and without eos_rad_rho_hi w is 1 "
                << "everywhere already.  Set the density window, or drop the gate."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(rad_t_lo > 0.0) || !(rad_t_lo < rad_t_hi)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "<" << block << ">/eos_rad_t_lo must be positive and strictly less "
                << "than eos_rad_t_hi; got T_lo = " << rad_t_lo
                << ", T_hi = " << rad_t_hi << std::endl;
      std::exit(EXIT_FAILURE);
    }
    tbl.rad_lt_lo = log10(static_cast<double>(rad_t_lo));
    tbl.rad_lt_hi = log10(static_cast<double>(rad_t_hi));
  }

  if (model.xhyd < 0.0 || model.yhel < 0.0 || (model.xhyd + model.yhel) > 1.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/eos_xh and eos_yhe must be non-negative and sum "
              << "to at most one (they are mass fractions); got X = " << model.xhyd
              << ", Y = " << model.yhel << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // defaults span irradiated atmospheres to stellar envelopes, molecular to fully
  // ionized
  Real xlo = pin->GetOrAddReal(block, "eos_logd_min", -14.0);
  Real xhi = pin->GetOrAddReal(block, "eos_logd_max", 2.0);
  Real ylo = pin->GetOrAddReal(block, "eos_logt_min", 1.5);
  Real yhi = pin->GetOrAddReal(block, "eos_logt_max", 8.0);
  // default resolution is 5x finer in log T: ionization and dissociation fronts are
  // ~kT/chi wide in ln T, while nothing varies fast in density
  Real dlog = pin->GetOrAddReal(block, "eos_dlog", 0.05);
  Real dlogd = pin->GetOrAddReal(block, "eos_dlogd", dlog);
  Real dlogt = pin->GetOrAddReal(block, "eos_dlogt", 0.2*dlog);

  if (!(xhi > xlo) || !(yhi > ylo) || !(dlogd > 0.0) || !(dlogt > 0.0)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << "> table bounds must satisfy eos_logd_max > "
              << "eos_logd_min, eos_logt_max > eos_logt_min, and positive eos_dlogd "
              << "and eos_dlogt"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  int nx = static_cast<int>(std::round((xhi - xlo)/dlogd)) + 1;
  int ny = static_cast<int>(std::round((yhi - ylo)/dlogt)) + 1;
  nx = std::max(nx, 4);
  ny = std::max(ny, 4);

  tbl.nx = nx;
  tbl.ny = ny;
  tbl.xmin = xlo;
  tbl.ymin = ylo;
  tbl.dx = (xhi - xlo)/static_cast<Real>(nx - 1);
  tbl.dy = (yhi - ylo)/static_cast<Real>(ny - 1);
  tbl.xmax = xlo + (nx - 1)*tbl.dx;
  tbl.ymax = ylo + (ny - 1)*tbl.dy;
  tbl.dxi = 1.0/tbl.dx;
  tbl.dyi = 1.0/tbl.dy;
  tbl.dens_cgs = dens_cgs;
  tbl.pres_cgs = pres_cgs;
  tbl.temp_cgs = temp_cgs;
  tbl.arad = eos_cgs::a_rad;
  // recorded for metallicity consistency checks in problem generators
  tbl.metal_ion = model.include_metal_ion;
  tbl.metal_mh = model.metal_mh;

  // ------------------------------------------------------------------------ sample grid
  Kokkos::realloc(tbl.tbl, ny, nx, ITNVAR);
  Kokkos::realloc(tbl.efbnd, nx-1);
  auto h_tbl = Kokkos::create_mirror_view(tbl.tbl);
  auto h_efb = Kokkos::create_mirror_view(tbl.efbnd);

  // central-difference steps of the node derivatives; the cross derivative needs the
  // larger step (round-off ~ 1/h^2)
  const double h1 = 1.0e-5;
  const double h2 = 1.0e-3;

  for (int j=0; j<ny; ++j) {
    double y = ylo + j*tbl.dy;
    for (int i=0; i<nx; ++i) {
      double x = xlo + i*tbl.dx;
      double fe, fp, fmu, fxe;
      double exp_, exm, eyp, eym, pxp, pxm, pyp, pym, mxp, mxm, myp, mym;
      double epp, epm, emp, emm, ppp, ppm, pmp, pmm, mpp, mpm, mmp, mmm;
      double xxp, xxm, xyp, xym, xpp, xpm, xmp, xmm;

      SampleModel(model, x, y, fe, fp, fmu, fxe);
      SampleModel(model, x+h1, y, exp_, pxp, mxp, xxp);
      SampleModel(model, x-h1, y, exm, pxm, mxm, xxm);
      SampleModel(model, x, y+h1, eyp, pyp, myp, xyp);
      SampleModel(model, x, y-h1, eym, pym, mym, xym);
      SampleModel(model, x+h2, y+h2, epp, ppp, mpp, xpp);
      SampleModel(model, x+h2, y-h2, epm, ppm, mpm, xpm);
      SampleModel(model, x-h2, y+h2, emp, pmp, mmp, xmp);
      SampleModel(model, x-h2, y-h2, emm, pmm, mmm, xmm);

      const double i2h1 = 0.5/h1;
      const double i4h2 = 0.25/(h2*h2);

      h_tbl(j,i, ITE) = fe;
      h_tbl(j,i, ITE+1) = (exp_ - exm)*i2h1;
      h_tbl(j,i, ITE+2) = (eyp - eym)*i2h1;
      h_tbl(j,i, ITE+3) = (epp - epm - emp + emm)*i4h2;

      h_tbl(j,i, ITP) = fp;
      h_tbl(j,i, ITP+1) = (pxp - pxm)*i2h1;
      h_tbl(j,i, ITP+2) = (pyp - pym)*i2h1;
      h_tbl(j,i, ITP+3) = (ppp - ppm - pmp + pmm)*i4h2;

      h_tbl(j,i, ITMU) = fmu;
      h_tbl(j,i, ITMU+1) = (mxp - mxm)*i2h1;
      h_tbl(j,i, ITMU+2) = (myp - mym)*i2h1;
      h_tbl(j,i, ITMU+3) = (mpp - mpm - mmp + mmm)*i4h2;

      h_tbl(j,i, ITXE) = fxe;
      h_tbl(j,i, ITXE+1) = (xxp - xxm)*i2h1;
      h_tbl(j,i, ITXE+2) = (xyp - xym)*i2h1;
      h_tbl(j,i, ITXE+3) = (xpp - xpm - xmp + xmm)*i4h2;
    }
  }

  // per x cell, the maximum of the floor energy over the cell from the analytic model,
  // slightly inflated so it also bounds the interpolant (EOSTable::EnergyFloorBound).
  // The applied floor is clamped to the lowest tabulated temperature
  // (EOSTable::SolveTemperatureFromP), so the bound covers that value too.
  const double pfl_cgs = pfloor*pres_cgs;
  const double tmin_cgs = pow(10.0, ylo);
  double efmax = 0.0;
  for (int i=0; i<nx-1; ++i) {
    double emax = 0.0;
    for (int k=0; k<=4; ++k) {
      double x = xlo + (i + 0.25*k)*tbl.dx;
      const double rho = pow(10.0, x);
      const double wr0 = tbl.rad_taper
          ? rad_taper::WeightOnly(x, tbl.rad_lrho_lo, tbl.rad_lrho_hi) : 1.0;
      // with the temperature gate w is unknown here: bound over both w = w_rho and w = 1
      const double wrs[2] = {wr0, tbl.rad_tgate ? 1.0 : wr0};
      for (int q=0; q<(tbl.rad_tgate ? 2 : 1); ++q) {
        const double wr = wrs[q];
        emax = std::max(emax, EnergyAtPressure(model, tbl.radiation, wr, rho, pfl_cgs));
        EOSCompositionState st = model.Evaluate(rho, tmin_cgs);
        const double etmin = rho*st.e_spec
                           + (tbl.radiation ? wr*eos_cgs::a_rad*tmin_cgs*tmin_cgs
                                              *tmin_cgs*tmin_cgs : 0.0);
        emax = std::max(emax, etmin);
      }
    }
    h_efb(i) = (emax/pres_cgs)*(1.0 + 1.0e-6);
    efmax = std::max(efmax, static_cast<double>(h_efb(i)));
  }
  tbl.efmax = efmax;

  Kokkos::deep_copy(tbl.tbl, h_tbl);
  Kokkos::deep_copy(tbl.efbnd, h_efb);
  tbl.tbl_h = h_tbl;       // host-side copies for host callers (see eos_table.hpp)
  tbl.efbnd_h = h_efb;
  tbl.active = true;

  // <block>/eos_table_dump: rank 0 writes the e and p node values to an ASCII file for
  // post-processing with the run's own EOS
  std::string dumpfile = pin->GetOrAddString(block, "eos_table_dump", "");
  if (!dumpfile.empty() && global_variable::my_rank == 0) {
    std::ofstream f(dumpfile);
    // every metadata line is a comment, so the body loads with a bare np.loadtxt
    f << "# AthenaK general-EOS table dump\n"
      << "# grid: nx ny xmin dx ymin dy   (x = log10 rho[cgs], y = log10 T[K])\n"
      << "# " << nx << " " << ny << " " << xlo << " " << tbl.dx << " "
      << ylo << " " << tbl.dy << "\n"
      << "# comp: X Y Z H2 ion metal_ion mh\n"
      << "# " << model.xhyd << " " << model.yhel << " " << model.MetalFraction() << " "
      << model.include_h2 << " " << model.include_ion << " "
      << model.include_metal_ion << " " << model.metal_mh << "\n"
      << "# then ny*nx rows, y (T) slowest: log10(e/rho)[cgs] log10(p/rho)[cgs]\n";
    f.precision(12);
    for (int j=0; j<ny; ++j) {
      for (int i=0; i<nx; ++i) {
        f << h_tbl(j,i,ITE) << " " << h_tbl(j,i,ITP) << "\n";
      }
    }
    f.close();
    std::cout << "General EOS: table dumped to '" << dumpfile << "'" << std::endl;
  }

  // ---------------------------------------------------------------------------- report
  EOSCompositionState smid = model.Evaluate(pow(10.0, 0.5*(xlo+xhi)),
                                            pow(10.0, 0.5*(ylo+yhi)));
  std::cout << "General EOS: tabulated on " << nx << " x " << ny
            << " nodes (" << (ITNVAR*sizeof(Real)*nx*ny)/(1024*1024) << " MB), "
            << "log10 rho[cgs] in [" << xlo << ", " << xhi
            << "], log10 T[K] in [" << ylo << ", " << yhi << "]" << std::endl
            << "             X = " << model.xhyd << ", Y = " << model.yhel
            << ", Z = " << model.MetalFraction()
            << ", H2 " << (model.include_h2 ? "on" : "off")
            << ", ionization " << (model.include_ion ? "on" : "off")
            << ", radiation " << (tbl.radiation ? "on" : "off") << std::endl;
  if (tbl.rad_taper) {
    std::cout << "             radiation TAPERED off between rho = " << rad_rho_hi
              << " (w = 1, LTE) and rho = " << rad_rho_lo
              << " g/cm^3 (w = 0, the two-stream owns it)" << std::endl;
  }
  if (tbl.rad_tgate) {
    std::cout << "             ...GATED on temperature: w = max(w_rho, w_T), w_T = 1 at "
              << "and above T = " << rad_t_hi << " K, 0 at and below " << rad_t_lo
              << " K -- a hot cell keeps LTE radiation whatever its density"
              << std::endl;
  }
  std::cout << "             metal ionization "
            << (model.include_metal_ion ? "on" : "off");
  if (model.include_metal_ion) {
    std::cout << " ([M/H] = " << model.metal_mh;
    if (model.include_metal_cond) {
      std::cout << ", condensation on";
    }
    std::cout << ")";
  }
  std::cout << std::endl
            << "             mu at the grid centre = " << smid.mu
            << ", n_e/n_tot = " << smid.xe << std::endl;
  return;
}
