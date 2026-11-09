# NOTE 2026-11-09 Orion -> viper: raycheap2-1109 CPU gates + convergence check (answer to TASK-2026-11-09-orion-raycheap)

Short answer:
- **G1 PASS, G2 PASS, G4 PASS. G3 (restart) FAILS.** The first differing variable is `m1_e`: 0/4 m1 vars bitwise, max |diff| 3.2e-1.
- Cause, from the logs (no code was changed): the restarted run does not carry the hesdirk2 `time2_auto_be` back-off
  state, and the restart file holds no `implicit_predictor` increment. The restarted leg therefore takes a different
  sequence of stage attempts and backward-Euler (BE) steps from cycle 21 on (details below).
  The same fallback sequence occurs on the deck without the mod keys (`G_cs_lo_shd_ema16_h2`, NEW and BASE), so this
  is probably not caused by the new key. I did not run that restart test, because the TASK says to stop at a failing gate.
- **Convergence check: FAIL against the PASS criterion.** The finest p is 0.84 (criterion >= 1.8). d is ~7e-6 to 3e-6,
  i.e. 10-30x below the expected ~1e-4 / 3e-5 / 8e-6. All 4 runs ended rc=0, 0 FATAL, t_end = 400.
  The four decks had already run, at the same time as the gates, before G3 failed.

## Build
- NEW `a39c9ec49` (origin/raycheap2-1109; the fork branches are on `origin` = jing-ze-ma/athenak here) and BASE
  `3330bc0dd` (origin/prod-1112-cand). Built from `git archive <sha>`, with kokkos from the gitlink
  d8e9af0316d77cff5059748958a81a8961db2a0b (4.6.2) via `git -C kokkos archive`.
- modules gcc/13 (g++ 13.1.0) and openmpi/4.1 (4.1.8), cmake/3.30. Kokkos SERIAL, PROBLEM none:
  `cmake -S src_<sha> -B b_<sha> -D CMAKE_BUILD_TYPE=Release -D Athena_ENABLE_MPI=ON -D CMAKE_CXX_COMPILER=mpicxx -D CMAKE_CXX_FLAGS=-ffp-contract=off; make -j24`
  (the C compiler that cmake picked is the system gcc 7.5, but there are no C sources). Both builds: 0 warnings, rc 0,
  about 7 min with the two builds running side by side.

## Gates (`tests_raycheap/run_gates.sh`, launcher adapted)
- Slurm job 227906 ran on one p.exclusive node (orion055, `--exclusive --mem=0`).
- The only change: `mpirun -np $np $b` became `srun -n $np --exact --mem-per-cpu=24G --cpu-bind=cores $b`, in a copy of
  the script (`run_gates_srun.sh`). A first try with 3G per rank was killed by OOM (a 1-rank cs-lo run peaks at about
  5 GB RSS). It was cancelled and its outputs deleted.
- All 10 runs: rc=0, FATAL=0.

```
G1 prod:   m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
G1 ema16:  m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
G1 h2:     m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
G2 mod acts (must DIFFER): m1: cycles 40/40 bitwise 0/4 vars, max |diff| 3.813e-01 m1_lat: cycles 40/40 bitwise 5/12 vars, max |diff| 3.364e-01
G3 restart: m1: cycles 40/40 bitwise 0/4 vars, max |diff| 3.222e-01 m1_lat: cycles 40/40 bitwise 7/12 vars, max |diff| 3.364e-01
G4 1 vs 4:  m1: cycles 40/40 bitwise 4/4 vars, max |diff| 0.000e+00 m1_lat: cycles 40/40 bitwise 12/12 vars, max |diff| 0.000e+00
```

### G3 detail (m_p2 = m_p1 cycle-20 rst + 20 cycles, compared with m_full = 40 cycles)
- Differing vars: m1_e, m1_f1, m1_f2, m1_f3 (every cell). In m1_lat, tt0, lat0, lat1, lat2 and lat3 differ;
  lat4, lat5, rx, ry, rz, jg and hg are bitwise.
- m_p1 matches m_full through cycle 20: the same fallback events, at cycles 1, 3, 6 and 11.
- After the restart, the m_p2 log has
  `### WARNING: restart file has no <rad_m1> implicit_predictor increment ...; the first implicit step cold starts and this restart is not bitwise.`
  At cycle 20, `pred_ok` was false, so `outputs/restart.cpp` wrote no predictor block.
- Stage fallbacks (`NOT ADMISSIBLE` -> BE):
  - m_full: cycles 1, 3, 6, 11, 20. After the 5th fallback, `time2_auto_be=16` makes cycles 21-36 BE.
    Totals: stage steps 3, BE steps 37, fallbacks 5.
  - m_p2: cycles 20, 22, 25, 30, 39. The back-off restarts at 1, 2, 4, 8 BE steps. So the 1, 2, 4, ... 16 escalation
    counter is not in the restart file, and m_p2 retries the stage solve at cycles 22, 25, 30 and 39, where m_full takes BE.
- The m_np4 and h_new (no mod keys) runs show the identical event list 1, 3, 6, 11, 20 and 3/37/5 counts. So the
  restart mismatch should also appear on the h2 deck without the mod keys (and on BASE). This is an inference; it was not run.
- If you want it: the discriminating test is the G3 pair on `G_cs_lo_shd_ema16_h2` with BASE (about 40 min at 1 rank on Orion).

## Convergence check
- Slurm job 227910 ran on one p.exclusive node (orion056). The 4 decks ran at the same time, 24 ranks each:
  `srun -n 24 --exact --mem-per-cpu=4G --cpu-bind=cores <NEW athena> -i Z_Q128T8_c<c>.athinput -d <dir>`.
  All runs: rc=0, 0 FATAL, "Terminating on time limit".

```
t_end [400.0, 400.0, 400.0, 400.0]
d ['7.197e-06', '5.116e-06', '2.850e-06']
p ['0.49', '0.84']
```

| cfl (dt) | wall | hesdirk2 stage steps | BE steps | stage fallbacks |
|---|---|---|---|---|
| 0.15 (1.34) | 575 s | 291 | 8 | 3 |
| 0.075 (0.67) | 1029 s | 597 | 1 | 0 |
| 0.0375 (0.335) | 1742 s | 1194 | 1 | 0 |
| 0.01875 (0.1675) | 3240 s | 2388 | 1 | 0 |

- The coarsest run had 3 stage fallbacks and 8 BE steps (12 `NOT POSITIVE|NOT ADMISSIBLE` log lines). The other three runs
  had only the first-step BE and no such lines. All runs: bicgstab breakdowns 3, line_jacobi fallbacks 1.
- Reading of p (inference): d is small (3e-6) and falls only 1.8x per halving. That looks like an error floor of the
  same size as the d values, not the O(dt^2) term, e.g. the solver tolerance (`implicit_tol = 1e-10` in the deck) or
  the physical-clock rotation/EMA events. The coarsest pair is also contaminated by the c0.15 fallbacks.
  I did not investigate further.

## Wall times
- Builds: about 7 min, both together.
- Gates job: 70 min in total. Per run, 1 rank unless noted:
  - p_new / p_base 837 / 831 s
  - e_new / e_base 785 / 788 s
  - h_new / h_base 1338 / 1340 s
  - m_full 1733 s; m_p1 1134 s; m_p2 1131 s cpu; m_np4 (4 ranks) 977 s
- Convergence job: 54 min (the finest deck sets it).

## Files
- Work dir `/orion/ptmp/jinma/raycheap_1109`: run logs, `gates.out`, `conv.out`, the sbatch scripts and the adapted `run_gates_srun.sh`.
- Restart dumps, bin outputs, source trees and build dirs were deleted after the comparisons (ptmp quota).
