# NOTE 2026-11-09 Orion -> viper: prod-1113-cand gates (answer to TASK-2026-11-09-orion-r136cap13)

Short answer: gates 1, 2 and 3 are bitwise PASS, and gate 4 PASS.
- No time2_auto_be back-off is active anywhere in the cap run: 0 stage fallbacks, and no M1AUTOBE block in any rst. So the cycle-20 split is the only split; no repeat at another cycle C was possible.
- NEW + cap gives the same results as the f34ada236 + cap run of NOTE-2026-11-09-orion-viper-r136fast-gates:
  - the hst files are byte-identical (16-digit %e);
  - the floorcap counters are identical;
  - |dEtot|/(dKE+dEint) is identical.

## Code and build
- NEW = `prod-1113-cand` bc35bfca8: merge of fastfloor-1109 (a24871d81, 263ef0393, 138e6ecb5, e5cccab04) into prod-1112-cand.
- BASE = `prod-1112-cand` 4e9f602c8.
- Both branches are on remote `origin`.
- Detached worktrees `/orion/ptmp/jinma/r136cap13_1109/{wt_new,wt_base}`, kokkos d8e9af03.
- Modules gcc/13 openmpi/4.1 cmake/3.30, Kokkos SERIAL. The two builds ran at the same time: 232 s (NEW) and 228 s (BASE), 0 warnings:
  `cmake -B build -D Athena_ENABLE_MPI=ON -D PROBLEM=he_star_m1 -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off && make -C build -j32`

## Deck
- The deck is the same as in the r136fast NOTE (same sed lines; the he_presn_m1_wedge.athinput file is identical in bc35bfca8, 4e9f602c8 and 138e6ecb5). Only the data-file paths point to the new dir.
- Cap keys are the previous values: `floorcap_rho = 9.818878e-13` (the TASK's "9.8e-13" rounded; this is 30 x the t = 0 minimum rho, as before), `floorcap_v = 1.0e6`, `floorcap_report = 1`.
- One Slurm job, 228016: p.shared, `--qos=serial`, -N 1 -n 17. Each arm ran as `srun -n <r> --exact --mem=15G --cpu-bind=cores`.
  - All arms started at once.
  - The restart leg S40 (`-r S20/rst/hepresn.00001.rst time/nlim=40`) started right after S20.

## Gates (final hydro_w + m1 double bins, `bin_convert.read_binary`, np.array_equal over all 163840 cells, 9 arrays)
| gate | arms | result |
|---|---|---|
| 1 keys absent, NEW vs BASE, 4 ranks, nlim 40 | G1_new vs G1_base | **PASS** bitwise (max rel 0) |
| 2 cap on, 1 rank vs 4 ranks, nlim 40 | C1 vs C4 | **PASS** bitwise; counters identical |
| 3 cap on, straight 40 vs 20 + restart to 40 | S40 vs C4 | **PASS** bitwise |
| 4 cap active | C4 | **PASS**: n_v > 0 from cycle 8; energy ratio as in f34ada236 |

Gate 3, back-off state:
- C4 (40 cycles): `stage steps=39 backward-Euler steps=1 stage fallbacks=0`. The one BE step is the first step.
  `time2_auto_be=16: skipped-stage BE steps=0 stage retries=0`. Wall split: `failed stage attempts 0`, `auto_be steps 0`. `backward-Euler redo 0`.
- S20: stage steps 19, BE 1 (first step), fallbacks 0, skipped-stage BE 0.
- S40: stage steps 20, BE 0, fallbacks 0. The restart continued with the stage slope, not a restart BE step.
- rst files (grep -a): M1AUTOBE = 0 in S20/hepresn.00001.rst, S40/hepresn.00002.rst and C4/hepresn.00001.rst; M1TIME2A = 1 in each.
  So no skip run was active at cycle 20. The back-off block is written only while one is active.
- No cycle of the straight run had an active back-off (0 stage fallbacks means t2_abe_len stays 0). So a repeat split at another C was not possible.
- Caveat: on this deck the rstabe restore path is therefore not exercised together with the cap. What is shown is that the merge keeps restarts bitwise with the cap on.

Gate 4, counter lines (cumulative, all calls):
```
C4  cycle 39: n_v 6.293504e+06 n_t 1.181082e+07 dKE 4.088578e+32 dEint 9.728892e+29 (erg, to E_rad)
C1  cycle 39: identical
S20 cycle 19: n_v 1.149696e+06 n_t 5.758976e+06 dKE 7.773454e+31 dEint 3.407634e+29
first n_v > 0: cycle 8, n_v 7.424e+03
```
- These are the same to all printed digits as f34ada236 + cap.
- e5cccab04 (T ceiling >= the EOS tfloor) is inert on this deck: T_rad > tfloor wherever the cap acts.
- Energy check (NEW + cap vs NEW, nlim 40, user hst Etot col 19; method of the previous NOTE):
  - dEtot = -1.0956213265e32 erg, so raw |dEtot|/(dKE+dEint) = 0.267.
  - The ledger residual difference is +1.4117e31, so |dResid|/(dKE+dEint) = 0.034.
  - Both are identical to the f34ada236 numbers (0.267 / 0.034). The trajectory difference between the two code versions is 0: C4/hepresn.user.hst and hydro.hst are byte-identical to r136fast_1109/C4, and so are G1_new's.

## Wall times (srun step, all arms concurrent on one node)
| arm | ranks | cycles | wall |
|---|---|---|---|
| G1_new / G1_base (no cap) | 4 | 40 | 205 / 217 s |
| C4 (cap) | 4 | 40 | 201 s |
| C1 (cap) | 1 | 40 | 406 s |
| S20 + S40 | 4 | 20 + 20 | 107 + 104 s |

The cap and no-cap walls are close this time: 201 vs 205 s, against 194 vs 263 s in r136fast. Node load differed between the two runs, so this is not a clean timing comparison.

## Files
- `/orion/ptmp/jinma/r136cap13_1109/<arm>/`: deck, log, hst, wall.txt; scripts/ (cmpbin.py, etot.py, percyc.py).
- All rst and bin files were deleted after the comparisons.
- The build worktrees wt_new / wt_base (0.55 GB each) are kept.
