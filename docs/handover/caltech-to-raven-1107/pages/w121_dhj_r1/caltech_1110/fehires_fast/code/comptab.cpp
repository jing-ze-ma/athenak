// comptab: composition on a (log rho, log T) grid with the AthenaK EOSCompositionModel
// (H2 + H/He ionisation + Na/K/Ca/Al/Mg/Fe donor ionisation + metal condensation), the
// model the dhj runs use.  The keys must match the run input: eos_xh, eos_yhe, eos_metal_mh.
//   WASP-121b 3x MHD runs (s10_e15 / b100_e13 / w121fresh_3x_f6): 0.7188 0.2420 0.4771
//   Caltech 1x hydro runs (w121r1_1101, r1planets_1101):          0.7381 0.2485 0.0
// Build (header-only, no Kokkos):  g++ -O2 -I <athenak>/src comptab.cpp -o comptab
// Run:  ./comptab <out.bin> [xh yhe metal_mh]        (defaults: the WASP-121b 3x values)
// out: int nx ny; double x0 dx y0 dy; then ny*nx*5 doubles: log10 xe, xh2, xhii, mu,
// log10 p_spec
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "eos/eos_composition.hpp"
int main(int argc, char **argv) {
  if (argc != 2 && argc != 5) {
    fprintf(stderr, "usage: comptab out.bin [xh yhe metal_mh]\n");
    return 1;
  }
  EOSCompositionModel m;
  m.xhyd = 0.7188; m.yhel = 0.2420; m.a_metal = 16.0; m.include_h2 = true;
  m.include_ion = true; m.include_metal_ion = true; m.metal_mh = 0.4771;
  m.include_metal_cond = true; m.metal_tcond = 0.0;
  if (argc == 5) {
    m.xhyd = atof(argv[2]); m.yhel = atof(argv[3]); m.metal_mh = atof(argv[4]);
  }
  double x0 = -19.0, dx = 0.05, y0 = 2.5, dy = 0.005;
  int nx = 381, ny = 441;
  FILE *f = fopen(argv[1], "wb");
  if (f == nullptr) { perror(argv[1]); return 1; }
  int h[2] = {nx, ny}; fwrite(h, sizeof(int), 2, f);
  double g[4] = {x0, dx, y0, dy}; fwrite(g, sizeof(double), 4, f);
  for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
    EOSCompositionState s = m.Evaluate(pow(10.0, x0 + i*dx), pow(10.0, y0 + j*dy));
    double o[5] = {log10(s.xe > 1e-300 ? s.xe : 1e-300), s.xh2, s.xhii, s.mu, log10(s.p_spec)};
    fwrite(o, sizeof(double), 5, f);
  }
  fclose(f);
  printf("comptab: %s xh %g yhe %g metal_mh %g\n", argv[1], m.xhyd, m.yhel, m.metal_mh);
  return 0;
}
