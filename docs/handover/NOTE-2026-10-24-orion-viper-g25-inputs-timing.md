# NOTE 2026-10-24 orion -> viper/raven: g25 inputs bundle + g25 IDORT timing (reply to TASK-2026-10-24-viper-orion-g25-inputs)

## 1. g25 inputs bundle (staged on Orion, the USER carries it to Raven)
`orion02.bc.mpcdf.mpg.de:/orion/ptmp/jinma/idort_dev/bundles/g25_raven/` (778 MB, `MD5SUMS` inside), already laid out as
the Raven target `/raven/ptmp/jinma/arepo2_gpu/g25/`:

| bundle path | md5 |
| --- | --- |
| ics/g25_red.hdf5 | 999c75f387dfc36d3db1040b3836d2c3 |
| ics/g25_full.hdf5 (the IC of the full long_bd run) | 40c5ce81bc13bc465ee9afc485e93269 |
| tables/EOS5_data | 7bdf6fc262375540bbf72d37ab657913 |
| tables/gs98_z0.02_x0.7_ma_ross.data | 352f7b44fc25b53b518cafc665501e6e |
| tables/gs98_z0.02_x0.7_ma_planck.data | 8c5453f63f02fe933261d92a0dd6cece |
| tables/coreB_side_g25_c1s.hdf5 | 7cc8a743fbb2920ebedd41623f5bb796 |
| param_long_bd_raven.txt (long_bd param.txt, the four table paths rewritten to .../g25/tables/) | edeff4d355a8a9d2556b1a0f8691834d |

The long_bd param.txt reads no other files (OutputListOn 0, SphericalSourceGravityTable none, InitCondFile g25_full).
Copy (from the laptop):

    rsync -av orion02.bc.mpcdf.mpg.de:/orion/ptmp/jinma/idort_dev/bundles/g25_raven/ ~/g25_raven/ && \
    rsync -av ~/g25_raven/ raven:/raven/ptmp/jinma/arepo2_gpu/g25/
    # on Raven: cd /raven/ptmp/jinma/arepo2_gpu/g25 && md5sum -c MD5SUMS

## 2. g25 IDORT timing (full g25 run long_bd, Arepo2-IDORT-dev bd36de23f3, -O3, 560 CPU ranks on Orion, 5.07 h wall)
Files in `orion_g25_timing/`: `idort_solver_rows.txt` (header, init solve, one full sync, one partial step),
`idort_timing_lines.txt` (the matching "IDORT: timing" lines), `cpu_step15575.txt` / `cpu_step16458.txt` (cpu.txt blocks),
`param_long_bd.txt`.

Whole run: 16537 IDORT steps = 314 full syncs (17-18.8 M cells) + 16223 partial steps (median 4.6 k cells).
Full syncs: mean 21.5 s, 7.8 Newton outers, 76 FGMRES its, 162 sweeps, 1.0 s positivity re-solves; total 1.87 h.
Partial steps: mean 0.23 s (setup + couple dominate, the solve is ~10 %); total 1.05 h. IDORT is under `misc` in cpu.txt:
78 % of the full-sync top step, 65.5 % of the whole run.

- Full sync, Sync-Point 15575 (17.46 M cells): step 21.02 s = setup 1.10 + solve 19.34 + couple 0.40 + diag 0.14;
  solve: sweeps 6.36, halo 3.71, orth 3.21, transport 1.91, lin 1.34, block 0.22, norms 0.28, other 2.31; reductions 1.56;
  8 outers, 85 matvecs, 77 FGMRES its, 163 sweeps; operator work max/mean 1.32. Top step 27.6 s (voronoi 2.8, treegrav 1.4, misc 21.6).
- Partial step, Sync-Point 16458 (4563 cells): step 0.090 s = setup 0.053 + solve 0.009 + couple 0.025;
  operator work max/mean 20 (the active cells sit on a few ranks).
- Worst step: the last one (Sync-Point 16537), 55.2 s, 21 outers / 235 FGMRES its, 35.5 s positivity re-solves.

## 3. GMRES basis high-water mark
Not logged in long_bd. IdortKrylovRestart = IdortKrylovMaxIter = 100, so each linear solve is at most ONE FGMRES cycle and
the basis (V: j+2 vectors, Z: j+1, float with IdortKrylovBasisF 1) grows lazily to the largest inner count of any single solve.
From the per-step totals: FGMRES its per Newton outer at full syncs median 9.2, max 14.2; max ceil(its/outers) over the run 15.
The absolute cap is 100 (101 V + 100 Z vectors); no solve is shown to reach it, but the per-solve maximum is not recorded.
An exact number needs a short run with the verbose solver print (one fgmres-cycle line per solve); say if you want it.

## 4. Open questions -> USER (not decided on Orion)
- IdortDsa 0 / IdortPrecSweepMode 4: current g25 settings. DSA 0 because FGMRES+DSA stagnated at full syncs (STATUS_perf,
  2026-10-15); sweep mode 4 (+ GmresOrth 1 + KrylovBasisF 1) chosen for star runs (-25 %), not the general default.
  Not declared final for production.
- GPU acceptance criteria: open, the user decides.

FYI: 926 steps report outcome 2 with 0 outers; all are 1-9 cell active sets already below IdortTol at entry (res_cell <= 4.9e-4; harmless label).
