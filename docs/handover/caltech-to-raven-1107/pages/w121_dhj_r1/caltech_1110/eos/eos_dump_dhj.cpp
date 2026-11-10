// Dump the AthenaK general-EOS composition model (same header the run's table is built from,
// eos_table.cpp BuildEOSTable) in the docs/handover/scripts/dhjcs.py text format.
// usage: eos_dump_dhj out.txt xh yhe metal_mh   (h2, ion, metal_ion, metal_cond ON, a_metal 16)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "eos/eos_composition.hpp"
int main(int argc, char **argv) {
  EOSCompositionModel m;
  m.xhyd = atof(argv[2]); m.yhel = atof(argv[3]); m.metal_mh = atof(argv[4]);
  m.include_h2 = true; m.include_ion = true; m.include_metal_ion = true;
  m.include_metal_cond = true; m.a_metal = 16.0;
  const double xlo=-18.0, xhi=0.0, ylo=1.5, yhi=6.0, dx=0.02, dy=0.004;
  const int nx = (int)std::lround((xhi-xlo)/dx)+1, ny=(int)std::lround((yhi-ylo)/dy)+1;
  FILE *f=fopen(argv[1],"w");
  fprintf(f,"# AthenaK general-EOS model dump (offline, eos_dump_dhj.cpp)\n");
  fprintf(f,"# grid: nx ny xmin dx ymin dy   (x = log10 rho[cgs], y = log10 T[K])\n");
  fprintf(f,"# %d %d %.10g %.10g %.10g %.10g\n",nx,ny,xlo,dx,ylo,dy);
  fprintf(f,"# comp: X Y mh = %g %g %g\n",m.xhyd,m.yhel,m.metal_mh);
  fprintf(f,"# \n# then ny*nx rows, y (T) slowest: log10(e/rho) log10(p/rho) mu\n");
  for (int j=0;j<ny;++j) for (int i=0;i<nx;++i) {
    EOSCompositionState s=m.Evaluate(pow(10.0,xlo+i*dx),pow(10.0,ylo+j*dy));
    fprintf(f,"%.12g %.12g %.8g\n",log10(s.e_spec),log10(s.p_spec),s.mu);
  }
  fclose(f); return 0;
}
