# NOTE 2026-11-09 Orion -> viper: raycheap2 restart-fix gates (rstabe-1109) + convergence reconciliation (answer to TASK-2026-11-09-orion-raycheap2)

Short answer:
- **R1 PASS:** the FIX restart is bitwise across the back-off.
- **R2 DIFFERS, as expected:** on BASE, m1 0/4 bitwise.
- **R3 PASS** on H and on P.
- **R4 PASS.**
- All 9 gate runs and all 8 convergence runs ended rc=0, 0 FATAL.
- **C1 (Q64) reproduces Raven to 3-4 digits:** d 8.680e-5, 3.029e-5, 7.497e-6; p 1.52, 2.01.
- **C2 (Q128T8, implicit_tol 1e-13, lin_tol 1e-14) gives nearly the same numbers as the 1e-10 / 1e-12 run:**
  - d 6.975e-6, 5.135e-6, 2.850e-6; p 0.44, 0.85.
  - The 1e-10 run gave d 7.197e-6, 5.116e-6, 2.850e-6; p 0.49, 0.84.
  - So the d ~ 3e-6 floor of the Q128T8 arm is **not** the solver tolerance.

## Builds
- Three commits, from `git archive` (the fork branches are on remote `origin` = jing-ze-ma/athenak), kokkos d8e9af03:
  - FIX `4e9f602c8` (rstabe-1109)
  - BASE `f34ada236` (prod-1112-cand)
  - MOD `9fbd2d3c6` (raycheap2-1109)
- modules gcc/13 (g++ 13.1.0), openmpi/4.1 (4.1.8), cmake/3.30. Kokkos SERIAL, PROBLEM none, Release, `CMAKE_CXX_COMPILER=mpicxx`,
  `-ffp-contract=off`. All three builds: rc 0, 0 warnings, 241 s with the three running side by side.
- Inputs for all runs: `tests_raycheap/` of MOD.

## Restart-fix gates
- Slurm job 227924 ran on one p.exclusive node (`--exclusive --mem=0`).
- Each run: `srun -n <np> --exact --mem-per-cpu=24G --cpu-bind=cores`. Runs that do not depend on each other ran at the same time.
- The restart legs use `-r <cycle-20 rst> time/nlim=40`.

```
R1 restart FIX:  m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
R2 restart BASE: m1: cycles 40/40 bitwise 0/4 vars, max |diff| 1.191e-03 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
R3 H FIX/BASE:   m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
R3 P FIX/BASE:   m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
R4 FIX 4 vs 1:   m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
```

### Fallback lists
Stage fallbacks are `NOT ADMISSIBLE` -> BE, given as stage@cycle. "Skipped" is the `time2_auto_be` skipped-stage BE steps.

| run | fallbacks | stage steps / BE steps / fallbacks | skipped / retries (succeeded) | wall |
|---|---|---|---|---|
| FIX H 40 (fix_H) | 1@1 1@3 2@6 1@11 2@20 | 3 / 37 / 5 | 31 / 5 (1) | 1357 s |
| FIX H 20 (fix_H_p1) | 1@1 1@3 2@6 1@11 | 0 / 20 / 4 | 15 / 3 (0) | 953 s |
| FIX H rst 20->40 (fix_H_p2) | 2@20 | 3 / 17 / 1 | 16 / 2 (1) | 410 s |
| FIX H 4 ranks (fix_H_np4) | 1@1 1@3 2@6 1@11 2@20 | 3 / 37 / 5 | 31 / 5 (1) | 766 s |
| BASE H 40 (base_H) | 1@1 1@3 2@6 1@11 2@20 | 3 / 37 / 5 | 31 / 5 (1) | 1347 s |
| BASE H 20 (base_H_p1) | 1@1 1@3 2@6 1@11 | 0 / 20 / 4 | 15 / 3 (0) | 946 s |
| BASE H rst 20->40 (base_H_p2) | 2@20 | 18 / 2 / 1 | 1 / 1 (1) | 493 s |
| FIX P / BASE P | none (be_only deck) | - | - | 835 / 834 s |

- **R1 log line:** fix_H_p2 prints `<rad_m1> rstabe-1109: time2_auto_be skip run restored (len 8, left 0)`. After the
  cycle-20 fallback, the restarted leg takes the same 16-step back-off as the straight run, and it ends with the same
  3 / 37 totals over the two legs.
- **The predictor warning** (`restart file has no <rad_m1> implicit_predictor increment ... not bitwise`) still prints on
  both restarted legs, but R1 is bitwise. This confirms that no predictor record is needed when `pred_ok` is false at the dump.
- **R2 (BASE):** the restarted leg has its first fallback at the restart step (2@20), then 1 BE step. It runs in stage mode
  from cycle 22 to 40 (18 stage steps) and has **no fallback after cycle 20**. The straight run is BE for cycles 21-36.
  The first differing variable is `m1_e` (cmpbin order). The m1 difference is 1.2e-3, and m1_lat is bitwise.
  - This differs from the earlier a39c9ec49 mod-deck G3, where the restarted leg failed again at 22, 25, 30 and 39.
  - That was a different deck (the mod keys) and a different base, a39c9ec49 / 3330bc0dd.
  - So BASE f34ada236 on deck H restarts into clean stage steps instead.

## Convergence reconciliation (MOD binary, 24 ranks per run)
- Two Slurm jobs, each on one p.exclusive node with the four levels running together (96 ranks per node):
  - 227926: C1
  - 227927: C2
- Each run: `srun -n 24 --exact --mem-per-cpu=2G --cpu-bind=cores`.
- The C2 keys were accepted on the command line, which is fresh-start safe because both keys exist in the deck. The log shows `lin_tol=1e-14`.

C1 `Z_Q64_c*.athinput` as they are:
```
t_end [400.0, 400.0, 400.0, 400.0]
d ['8.680e-05', '3.029e-05', '7.497e-06']
p ['1.52', '2.01']
```
C2 `Z_Q128T8_c*.athinput rad_m1/implicit_tol=1.0e-13 rad_m1/implicit_lin_tol=1.0e-14`:
```
t_end [400.0, 400.0, 400.0, 400.0]
d ['6.975e-06', '5.135e-06', '2.850e-06']
p ['0.44', '0.85']
```

| run | stage steps | BE steps | stage fallbacks | wall |
|---|---|---|---|---|
| C1 c0.15 | 293 | 6 | 2 | 621 s |
| C1 c0.075 | 597 | 1 | 0 | 1085 s |
| C1 c0.0375 | 1194 | 1 | 0 | 1851 s |
| C1 c0.01875 | 2388 | 1 | 0 | 3369 s |
| C2 c0.15 | 291 | 8 | 3 | 622 s |
| C2 c0.075 | 597 | 1 | 0 | 1073 s |
| C2 c0.0375 | 1194 | 1 | 0 | 1876 s |
| C2 c0.01875 | 2388 | 1 | 0 | 3688 s |

Reading:
- Q64 (EMA64, T_rot 10.72) matches viper/Raven's numbers on Orion: d 8.7e-5, 3.0e-5, 7.5e-6; p 1.5 / 2.0.
- Q128T8 (EMA128, T_rot 21.44) has d 10-30x smaller and p < 1. Tightening the tolerances by 1000x changes d by < 3 %
  (the finest d is the same to 4 digits).
- Inference: the Q128T8 differences are dominated by something that is not dt^2 truncation and not solver
  tolerance, e.g. an O(dt) or dt-independent term from the physical-clock rotation/EMA events at that cadence. Its
  small signal (d ~ 5e-6) is buried under it.
- In both sets, the c0.15 level has 2-3 stage fallbacks, which also spoils the first p.

## Wall times
- Builds: 241 s.
- Gates job: 31 min.
- C1 job: 56 min; C2 job: 62 min (C1 and C2 ran at the same time).

## Files
- `/orion/ptmp/jinma/raycheap2_1109`: logs (`gates.out`, `conv_rc2_C1.out`, `conv_rc2_C2.out`, run.log/wall.txt per run) and sbatch scripts.
- Builds, sources, rst and bin outputs were deleted after the comparisons (ptmp quota).
