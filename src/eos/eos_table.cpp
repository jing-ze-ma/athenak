//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_table.cpp
//! \brief builds the (log10 rho, log10 T) interpolation table for the general EOS by
//! sampling the analytic composition model of eos_composition.hpp.
//!
//! Everything here runs ONCE, on the host, at problem setup. The node derivatives that
//! the bicubic Hermite interpolation needs are obtained by central differencing the
//! analytic model itself, not by differencing the table: differencing the model gives
//! node derivatives whose accuracy is independent of the grid spacing, which is what
//! makes the interpolant fourth-order accurate rather than second.

#include <math.h>

#include <algorithm>
#include <iostream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "eos/eos_composition.hpp"
#include "eos/eos_table.hpp"

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
//! \brief total internal energy density (gas plus radiation if enabled) at a given
//! density and TOTAL pressure, in cgs. Used only to precompute the pressure-floor bound,
//! so a plain bisection is fast enough and avoids depending on the table being built.

double EnergyAtPressure(const EOSCompositionModel &m, const bool rad,
                        const double rho, const double ptarget) {
  double ylo = -2.0, yhi = 12.0;
  for (int it=0; it<200; ++it) {
    double y = 0.5*(ylo + yhi);
    double t = pow(10.0, y);
    EOSCompositionState s = m.Evaluate(rho, t);
    double p = rho*s.p_spec + (rad ? eos_cgs::a_rad*t*t*t*t/3.0 : 0.0);
    if (p > ptarget) {
      yhi = y;
    } else {
      ylo = y;
    }
    if ((yhi - ylo) < 1.0e-13) break;
  }
  double t = pow(10.0, 0.5*(ylo + yhi));
  EOSCompositionState s = m.Evaluate(rho, t);
  return rho*s.e_spec + (rad ? eos_cgs::a_rad*t*t*t*t : 0.0);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void BuildEOSTable
//! \brief read the table and composition parameters, sample the analytic model onto the
//! grid, and copy the result to the device.

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
  // [M/H] in dex for the metal donors. It DEFAULTS to <problem>/met when that exists,
  // because a problem generator that has a metallicity uses it for its opacity, and an
  // atmosphere that is opaque at one metallicity and conducting at another is not a
  // physical configuration -- it is a parameter that was updated in one place only. Set
  // eos_metal_mh explicitly to override, which a pgen may then reject.
  Real mh_default = 0.0;
  if (pin->DoesParameterExist("problem", "met")) {
    mh_default = pin->GetReal("problem", "met");
  }
  model.metal_mh = pin->GetOrAddReal(block, "eos_metal_mh", mh_default);
  model.include_metal_cond = pin->GetOrAddBoolean(block, "eos_metal_condensation", false);
  model.metal_tcond = pin->GetOrAddReal(block, "eos_metal_tcond", 0.0);
  tbl.radiation = pin->GetOrAddBoolean(block, "eos_radiation", false);

  if (model.xhyd < 0.0 || model.yhel < 0.0 || (model.xhyd + model.yhel) > 1.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "<" << block << ">/eos_xh and eos_yhe must be non-negative and sum "
              << "to at most one (they are mass fractions); got X = " << model.xhyd
              << ", Y = " << model.yhel << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // ------------------------------------------------------------------------- table grid
  // Defaults span the conditions this EOS is meant for: from a tenuous irradiated
  // atmosphere to a gas-giant or stellar envelope interior, and from molecular
  // temperatures up to fully ionized.
  Real xlo = pin->GetOrAddReal(block, "eos_logd_min", -14.0);
  Real xhi = pin->GetOrAddReal(block, "eos_logd_max", 2.0);
  Real ylo = pin->GetOrAddReal(block, "eos_logt_min", 1.5);
  Real yhi = pin->GetOrAddReal(block, "eos_logt_max", 8.0);
  // Resolution is specified separately in the two directions, and by default is five
  // times finer in temperature. That is not a tuning knob but a property of the physics:
  // an ionization or dissociation front is a factor exp(-chi/kT), so its width in ln T is
  // kT/chi -- about a hundredth of a decade at the low densities where hydrogen ionizes
  // near 3000 K -- while nothing in the EOS varies faster than smoothly in density. A
  // grid that is uniform in the two logs therefore wastes most of its nodes and still
  // fails to resolve the only feature that is hard.
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
  // recorded so that a problem generator can verify its own metallicity matches the one
  // the EOS was built with; a mismatch is silent otherwise
  tbl.metal_ion = model.include_metal_ion;
  tbl.metal_mh = model.metal_mh;

  // ------------------------------------------------------------------------ sample grid
  Kokkos::realloc(tbl.tbl, ny, nx, ITNVAR);
  Kokkos::realloc(tbl.efbnd, nx-1);
  auto h_tbl = Kokkos::create_mirror_view(tbl.tbl);
  auto h_efb = Kokkos::create_mirror_view(tbl.efbnd);

  // Central-difference steps for the node derivatives. The first derivatives use the
  // smaller step, where truncation and round off are both near 1e-10; the cross
  // derivative needs the larger one because its round-off error scales as 1/h^2.
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

  // ------------------------------------------------- pressure-floor bound, in CODE units
  // For each x cell, the maximum of e(rho, pfloor) over that cell, evaluated from the
  // ANALYTIC model rather than from the interpolant and then inflated slightly, so that
  // it bounds the interpolated value too. See EOSTable::EnergyFloorBound().
  const double pfl_cgs = pfloor*pres_cgs;
  double efmax = 0.0;
  for (int i=0; i<nx-1; ++i) {
    double emax = 0.0;
    for (int k=0; k<=4; ++k) {
      double x = xlo + (i + 0.25*k)*tbl.dx;
      emax = std::max(emax, EnergyAtPressure(model, tbl.radiation, pow(10.0,x), pfl_cgs));
    }
    h_efb(i) = (emax/pres_cgs)*(1.0 + 1.0e-6);
    efmax = std::max(efmax, static_cast<double>(h_efb(i)));
  }
  tbl.efmax = efmax;

  Kokkos::deep_copy(tbl.tbl, h_tbl);
  Kokkos::deep_copy(tbl.efbnd, h_efb);
  tbl.active = true;

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
            << ", radiation " << (tbl.radiation ? "on" : "off") << std::endl
            << "             metal ionization "
            << (model.include_metal_ion ? "on" : "off");
  if (model.include_metal_ion) {
    std::cout << " ([M/H] = " << model.metal_mh;
    if (model.metal_tcond > 0.0) {
      std::cout << ", rainout below " << model.metal_tcond << " K";
    } else if (model.include_metal_cond) {
      std::cout << ", condensation on";
    }
    std::cout << ")";
  }
  std::cout << std::endl
            << "             mu at the grid centre = " << smid.mu
            << ", n_e/n_tot = " << smid.xe << std::endl;
  return;
}
