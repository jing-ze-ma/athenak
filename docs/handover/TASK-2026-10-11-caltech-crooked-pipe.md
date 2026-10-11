# TASK raven -> Caltech: crooked pipe (JSD12 Sect 5.2) on the merged prodcand code, REC and REC-be (CPU, ~2-3 h)

This comes from the Raven literature-suite session (suite-1011). The same 4 runs are queued on Raven as job 31068508,
but the account's job cap holds it with no start estimate. **Whichever site finishes first wins.** Post your NOTE as
soon as you have results. If Raven's job lands first, Raven will post a NOTE and this task is closed.

## Code (fork github.com/jing-ze-ma/athenak)
- Branch `prodcand-1011`, commit **4c033dcf** ("Merge moveq-1011 ... into prodcand-1011"). This is the merged
  production candidate the Raven suite uses (Raven binary athena_cpu_4c033dcf, md5 1f79726459c6e3484c680b578bdeae07).
  If the branch head has moved, check out 4c033dcf itself.
- Build: **CPU, built-in pgens (no -D PROBLEM), MPI, OpenMP**, Release, e.g.
  `cmake -S . -B build_cpu -D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx`.
  Report the md5 of your binary.

## Files: docs/handover/suite-pipe-1011/ on this branch
- `pipe.athinput` (md5 in the commit): the JSD12 tophat in Cartesian.
  - (0,7)x(-2,2) at 512x256, one MeshBlock.
  - 7 opaque walls (rho 10, sigma_a 2000) via `m1_test = cyl` rectangles; the pipe is rho 0.01, sigma_a 0.2 (powerlaw kappa ∝ rho^(1/3)).
  - c = 300; units cm, 1e-8 s, 0.5 keV; gamma - 1 = 1/58.32 (c_v = 1e15 erg/g/keV); T0 = 0.1 (0.05 keV) in equilibrium.
  - Isotropic T = 1 (0.5 keV) bath on the whole x1min face (`implicit_ebath_x1min = 1`). Walls cover x1min for |y| > 0.5, so this acts as the |y| < 0.5 source.
  - x2 periodic (vet_sc needs it; the top and bottom rows are opaque wall, so this is the same as outflow).
  - Static gas: `evolution = kinematic`, `rsolver = advect`, `dbg_gas_force(_trans) = false`, so v stays 0 and rho is constant while T evolves (JSD12: "velocities are always kept zero").
- `run_pipe.sh`: runs the 4 cases in parallel. It contains the exact arm keys.
- `fixkeys.py`: adds command-line keys that are missing from the deck, since AthenaK makes such keys fatal.
- `pipe.py`: probe histories. Set `ATHENAK_VIS=<repo>/vis/python` for bin_convert.

## Runs (one node; 4 processes x NT OpenMP threads, e.g. NT = cores/4)
```
BIN=<your athena> NT=16 PY=python3 bash docs/handover/suite-pipe-1011/run_pipe.sh <outdir>
```
| name | arm | implicit_cfl (c dt/dx) | dt (1e-8 s) | tlim | steps |
|---|---|---|---|---|---|
| REC_c220 | REC | 220 | 1.0e-2 | 10 | 1000 |
| RECbe_c220 | REC-be | 220 | 1.0e-2 | 10 | 1000 |
| REC_c2200 | REC | 2200 | 0.1 | 10 | 100 |
| REC_c22 | REC | 22 | 1.0e-3 | 1 | 1000 |

- **REC**: the recommended vet_sc keys (hr keys + `implicit_blend_xthin=0 implicit_thin_switch=fse implicit_thin_corr=hcs implicit_thin_recon=jse implicit_thin_bq=lin`) with `time_scheme=hesdirk2 time2_vet_sc=predict time2_vstage=true`, and `implicit_vimp=false` because the medium is static.
- **REC-be**: the same keys with `time_scheme=be`.
- Raven timing: about 5 s per step on 8 threads, so the 1000-step runs take about 1-1.5 h with 16+ threads.

## Report: NOTE-2026-10-11-caltech-crooked-pipe.md
- Every run: rc, FATAL, NON-CONVERGED, Picard mean/max, the hesdirk2 "NOT ADMISSIBLE"/fallback counts, wall time.
- Run `python3 pipe.py <outdir>/*/` and give its summary lines (T at the end of each probe and the time each probe reaches T = 0.5).
- Attach `probes.txt` of every run (t, the 5 probe T, max|vx|), or paste it if short.
- Check `max|vx| = 0` (static gas).
- Name the dt-dependence: REC_c220 vs REC_c2200 at t = 1, 3 and 10, and REC_c22 vs REC_c220 at t = 0.5 and 1.0, as the maximum relative difference over the probes.
- Ledger line as usual.
