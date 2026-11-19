# rp16k: reduced RSG cold-start repro (viper -> Orion, 11-19)

This is a spherical-polar wedge cut from rsg.01658:
- 16x16 cells at the cs cell size, about the centre of cs panel 1;
- the full 548-cell radial grid;
- tangential momenta kept;
- E, F_r = f0x1 + 4/3 E v_r.

It is a fresh start via problem/he_ic_remap_file (the remap file is next to the input). The RSG prod keys are used, with
vet_gd_gl_nphi=32 and vet_gd_allow_clamp=true.
make_repro.py builds it from a cs restart.

To run: 1 Raven CPU node, ~80 s, about 3 cycles at dt 1000. The paths in rp16k.athinput point to /raven/ptmp/jinma; override
problem/ remap path keys as needed.

What clean-cc 6f3aa7d07 gives:
- at dt 1000, pass-0 dE 2.0 (the cap) and 2.6k E<0 cells in the inner bins, not converged;
- at dt 63, it converges.
