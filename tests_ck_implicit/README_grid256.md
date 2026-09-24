# dhj-grid256: team-scratch level fallback + 8-coefficient radial stretch (09-25)

Code: commit 1be3eb7a on dhj-grid256 (base 82d7ceb1).

1. `TeamScratchLevel` (src/athena.hpp) asks Kokkos whether the level-0 team scratch of a
   kernel can launch with any team size; only if not (MHD fluxes at nx1 > ~264 on MI300A,
   64 KB LDS) the kernel uses level 1. Applied to hydro/MHD/dyn-GRMHD fluxes, hydro/MHD and
   resistive updates, viscosity and ohmic EMF kernels. Every size that launches today stays
   on level 0; host backends unchanged.
2. `NSTRETCH_R_POLY` 4 -> 8 (src/coordinates/grid_stretch.hpp); c5..c8 are read only when
   the input gives them (a command-line override of c5..c8 therefore needs the lines in the
   input file, see below).

Work dir: /viper/ptmp2/jinma/grid256_0925 (cpu.sh, cpu2.sh, gpu.sub, gridcheck.py, log/).

## Gates

- CPU, base vs new, every output file compared with cmp:
  cs MHD from bench/cs_mhd_prod4/rst/dhj.00159.rst (rot 79.50): 4/4 files identical;
  cs hydro from bench/cs_hyd4_prod/rst/dhj.00195.rst (rot 97.50): 4/4 identical;
  sp (inputs/tests/dhj_ck_spherical.athinput, 3 cycles): 7/7 identical (log/cpu.out,
  log/cpu2.out; the first cs "new" arms died with a bus error in the 00:44 login-node crash
  and were rerun by cpu2.sh).
- tst: test_hydro_fofc_sod_cpu (3 passed), test_mhd_balsara_vortex_cpu (1),
  test_hydro_fofc_sp_cpu (1) (log/tst.out).
- GPU nx1 = 128, prod4 rst 159, 10 cycles, base vs new: restart and .hst bitwise
  (job 11967857).
- GPU nx1 = 256 / 320, 8-coefficient grid, 64x64 per face, MHD + T4 + c2 + ck_impl_every=4,
  22 cycles from the IC, apudev, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1 (job 11967951,
  runs/n*_r*): no team-size abort.

| nx1 | ms/cycle (cycles 2-22), r1 / r2 | dt cycle 0 / 22 [s] |
| --- | --- | --- |
| 256 | 182.2 / 182.4 | 11.45 / 10.87 |
| 320 | 223.8 / 218.4 | 9.16 / 8.62 |

Job 11967857 failed at 256/320 because `mesh/f_stretch_r_c5..c8=` on the command line is a
FATAL when the input file has no such parameter; gpu.sub now uses dhj_c8.athinput (the
dhj_weak input plus c5..c8 = 0 lines, overridden on the command line).

## Grid (gridcheck.py -> gridcheck.txt, gridcheck.png)

nx1 = 256 with c = 0.418017, 4.818585, -84.625150, 397.754052, -955.667782, 1270.732651,
-887.750373, 254.765233: x1f[0] = 9.44e9 and x1f[256] = 2.0556e10 cm (= prod4 x1min/x1max),
monotonic, dr 1.91e7 (i = 61) .. 7.17e7 cm (i = 232), max |dr ratio - 1| = 0.165. Cells per
scale height (p10/median, day | night) match design256.txt: top 16.4/17.2 | 2.6/5.3, upper
14.9/15.5 | 4.7/7.8, photosphere/jet base 7.5/15.0 | 4.8/5.6, jet 4.7/5.2 | 4.6/5.1, deep
4.6/4.9 | 4.6/4.9 (today's 128 grid: 7.5/7.6, 6.6/7.0, 3.4/5.2, 3.7/4.6, 4.8/5.0 on the day
side).
