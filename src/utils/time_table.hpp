#ifndef UTILS_TIME_TABLE_HPP_
#define UTILS_TIME_TABLE_HPP_
//========================================================================================
// AthenaK astrophysical fluid dynamics code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file time_table.hpp
//! \brief hg262-drive-1113: a host-side table of columns against time, read from a text
//! file (whitespace separated, '#' starts a comment, first column = time, strictly
//! increasing), evaluated piecewise-linearly in t and held at the end values outside its
//! range.  Integral() is the EXACT integral of that piecewise-linear function, so a
//! budget that books Integral(t_n, t_n+1) per step sums to Integral(t_0, t_N) up to
//! round-off.  A pure function of t: restart-safe without any stored state.

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"

namespace timetable {

struct Table {
  std::vector<Real> t;                 // the time nodes
  std::vector<std::vector<Real>> c;    // c[q][n]: column q+1 of the file at node n
  int NCol() const {return static_cast<int>(c.size());}
  bool Empty() const {return t.empty();}

  //! read the file; ncol = number of data columns wanted after t (the rest ignored).
  //! Returns an empty string on success, else the reason.
  std::string Read(const std::string &fn, const int ncol) {
    t.clear();
    c.assign(ncol, std::vector<Real>());
    std::ifstream f(fn);
    if (!f) {return "cannot open " + fn;}
    std::string line;
    int ln = 0;
    while (std::getline(f, line)) {
      ++ln;
      const std::size_t h = line.find('#');
      if (h != std::string::npos) {line.resize(h);}
      std::istringstream ss(line);
      double x;
      std::vector<double> v;
      while (ss >> x) {v.push_back(x);}
      if (v.empty()) {continue;}
      if (static_cast<int>(v.size()) < ncol + 1) {
        return fn + " line " + std::to_string(ln) + ": fewer than " +
               std::to_string(ncol + 1) + " columns";
      }
      if (!t.empty() && !(static_cast<Real>(v[0]) > t.back())) {
        return fn + " line " + std::to_string(ln) + ": time not strictly increasing";
      }
      t.push_back(static_cast<Real>(v[0]));
      for (int q = 0; q < ncol; ++q) {c[q].push_back(static_cast<Real>(v[q+1]));}
    }
    if (t.empty()) {return fn + ": no data rows";}
    return std::string();
  }

  //! the segment n with t[n] <= x < t[n+1] (clamped to 0 .. N-2); N = 1: 0
  int Seg(const Real x) const {
    const int n = static_cast<int>(t.size());
    if (n < 2 || x <= t[0]) {return 0;}
    if (x >= t[n-1]) {return n - 2;}
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const int mid = (lo + hi)/2;
      if (t[mid] <= x) {lo = mid;} else {hi = mid;}
    }
    return lo;
  }

  //! column q at time x: linear inside, held at the end values outside
  Real At(const int q, const Real x) const {
    const int n = static_cast<int>(t.size());
    const std::vector<Real> &a = c[q];
    if (n == 1 || x <= t[0]) {return a[0];}
    if (x >= t[n-1]) {return a[n-1];}
    const int s = Seg(x);
    const Real w = (x - t[s])/(t[s+1] - t[s]);
    return (1.0 - w)*a[s] + w*a[s+1];
  }

  //! exact integral of column q over [x0, x1] (x1 >= x0) of the function At()
  Real Integral(const int q, const Real x0, const Real x1) const {
    if (!(x1 > x0)) {return 0.0;}
    const int n = static_cast<int>(t.size());
    const std::vector<Real> &a = c[q];
    if (n == 1) {return a[0]*(x1 - x0);}
    Real s = 0.0;
    // below the first node and above the last: the held end values
    if (x0 < t[0]) {s += a[0]*(std::fmin(x1, t[0]) - x0);}
    if (x1 > t[n-1]) {s += a[n-1]*(x1 - std::fmax(x0, t[n-1]));}
    const Real lo = std::fmax(x0, t[0]), hi = std::fmin(x1, t[n-1]);
    if (hi > lo) {
      for (int k = Seg(lo); k < n - 1 && t[k] < hi; ++k) {
        const Real a0 = std::fmax(lo, t[k]), a1 = std::fmin(hi, t[k+1]);
        if (a1 > a0) {s += 0.5*(At(q, a0) + At(q, a1))*(a1 - a0);}
      }
    }
    return s;
  }
};

}  // namespace timetable

#endif  // UTILS_TIME_TABLE_HPP_
