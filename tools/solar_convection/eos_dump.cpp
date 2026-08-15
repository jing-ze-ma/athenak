// Dump the AthenaK general-EOS composition model onto a (log10 rho, log10 T) grid so that
// offline analysis can recover T, mu, Gamma_1 and grad_ad from the (rho, e) pairs in the
// .bin dumps.  The solver itself cannot output temperature for Newtonian hydro, and the
// tabulated EOS lives in a device array, so this is the only faithful route: it includes
// the SAME header the code builds its table from, so the physics cannot drift.
//
// Build:  g++ -O2 -I<athenak>/src -o eos_dump eos_dump.cpp
// Run:    ./eos_dump <out.bin> [xh] [yhe] [a_metal] [h2] [ion]
//
// Output: little-endian float64. Header: nd nt logd_min dlogd logt_min dlogt (6 doubles),
// then nd*nt records of (log10 e_spec, log10 p_spec, mu, xh2, xhii), T fastest.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "eos/eos_composition.hpp"

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s out.bin [xh yhe a_metal h2 ion]\n", argv[0]);
    return 1;
  }
  EOSCompositionModel m;
  if (argc > 2) m.xhyd    = atof(argv[2]);
  if (argc > 3) m.yhel    = atof(argv[3]);
  if (argc > 4) m.a_metal = atof(argv[4]);
  if (argc > 5) m.include_h2  = atoi(argv[5]) != 0;
  if (argc > 6) m.include_ion = atoi(argv[6]) != 0;

  // Finer than the run's own table: this is offline, cost does not matter.
  const double logd_min = -14.0, logd_max = 0.0, dlogd = 0.02;
  const double logt_min = 1.5,   logt_max = 7.0, dlogt = 0.004;
  const int nd = static_cast<int>((logd_max - logd_min)/dlogd) + 1;
  const int nt = static_cast<int>((logt_max - logt_min)/dlogt) + 1;

  FILE *f = fopen(argv[1], "wb");
  if (f == nullptr) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
  double hdr[6] = {static_cast<double>(nd), static_cast<double>(nt),
                   logd_min, dlogd, logt_min, dlogt};
  fwrite(hdr, sizeof(double), 6, f);

  std::vector<double> row(5*nt);
  for (int i = 0; i < nd; ++i) {
    double rho = pow(10.0, logd_min + i*dlogd);
    for (int j = 0; j < nt; ++j) {
      double t = pow(10.0, logt_min + j*dlogt);
      EOSCompositionState s = m.Evaluate(rho, t);
      row[5*j+0] = log10(s.e_spec);
      row[5*j+1] = log10(s.p_spec);
      row[5*j+2] = s.mu;
      row[5*j+3] = s.xh2;
      row[5*j+4] = s.xhii;
    }
    fwrite(row.data(), sizeof(double), row.size(), f);
  }
  fclose(f);
  fprintf(stderr, "wrote %s: %d x %d grid\n", argv[1], nd, nt);
  return 0;
}
