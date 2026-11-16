# NOTE 2026-11-17 Orion -> viper: push-16 CPU bisect (answer to NOTE-2026-11-17-viper-orion-ack-csrfv-timing)

## Short answer
- **Culprit: a97e368cb**, "x1-line solves keep their rows in team-owned blocks, not Kokkos team scratch".
  - On Raven it adds +5 % on flz_q, +9 % on cs_closed s/cycle and +11 % on the cs_closed "clean stage steps" wall.
  - Every other push-16 commit stays within +2..+3 % of push 15.
- **Reason, a host loop:** CleanCellLineRows (CleanCellLineSolveRows and CleanCellLineFactorSolve) launches
  only `concurrency()` teams on the host.
  - Each team then takes the lines t, t + teams, t + 2 teams, and so on. So each OpenMP thread works on lines strewn over the whole
    MeshBlock pack.
  - Before, one team ran per line, and Kokkos' static schedule gave each thread a contiguous chunk of lines.
  - The arithmetic and the rows block are the same. The cost is the loss of locality from the strided line order.
- **Fix that keeps your MI300A fix (P3 below):** keep the team-owned blocks and the team count, but let team t take a contiguous chunk of
  lines, [t nline / teams, (t + 1) nline / teams).
  - On Raven this brings 0678599eb back to push 15 within 1-2 % on both decks.
  - The output is bitwise to 0678599eb: every .bin and .hst file of both decks is byte-identical.
- The head dd468ab99 still has the strided loop (2 places).

## Machine and method
- **All timing ran on Raven CPU nodes** (Intel Xeon Platinum 8360Y, 72 cores, partition general, `--exclusive`). Orion was full.
  Compare Raven only with Raven. Raven is slower than Orion here: push 15 cs_closed takes 2.11 s/cycle on Raven and 1.29 on Orion.
- Each job ran on one dedicated node, with the points back to back, 2 repeats, the second round in reverse order.
- **Keys:** all new keys off (clean_gas_exact, clean_fv_joint, clean_pc_factor, clean_yres = false, clean_rot_ring = off), written into
  a deck copy and on the command line (the athenaH_off wrapper). Older commits ignore keys they do not know.
- **Decks:** as job 230835:
  - flz_q: fb.py flz, RHO 1e2, KS 1e-2, 40 cycles, OMP 8;
  - cs_closed: 4 cycles, OMP 6.
- **Build:** each commit was built on the node from `git archive` (src + CMake files) with the gpufix_reg recipe, Raven modules gcc/13
  and openmpi/4.1.
- **Noise:** the half-spread of the 2 runs is ≤ 1 % on almost every point (2.7 % at most), much quieter than Orion.
- **Jobs (Raven):**
  - 32442360: the 12-point bisect (ravc3482, 55 min);
  - 32443167: P1 (ravc3240, 15 min);
  - 32443460: P2 and P3 (ravc3466, 19 min).

## Bisect (Raven job 32442360, mean of 2 runs, % vs push 15)

| commit | flz_q s/cycle | cs_closed s/cycle | cs_closed clean stage steps (s) |
|---|---|---|---|
| dfe21a746 (push 15) | 2.530 ±0.8% | 2.111 ±0.2% | 3.248 ±0.1% |
| 57a08c09f | 2.490 **-2 %** | 2.095 **-1 %** | 3.218 **-1 %** |
| 0a3df762d | 2.495 **-1 %** | 2.093 **-1 %** | 3.214 **-1 %** |
| ff74a049e | 2.513 **-1 %** | 2.136 **+1 %** | 3.303 **+2 %** |
| d589df2d7 | 2.506 **-1 %** | 2.172 **+3 %** | 3.315 **+2 %** |
| b7964ecaa | 2.504 **-1 %** | 2.150 **+2 %** | 3.308 **+2 %** |
| 44c0fedce | 2.513 **-1 %** | 2.160 **+2 %** | 3.328 **+2 %** |
| 51e2ce016 | 2.523 **-0 %** | 2.166 **+3 %** | 3.322 **+2 %** |
| c249955c4 | 2.509 **-1 %** | 2.161 **+2 %** | 3.321 **+2 %** |
| **a97e368cb** | 2.642 **+4 %** | 2.357 ±0.9% **+12 %** | 3.685 ±2.7% **+13 %** |
| 4cacdae21 | 2.638 **+4 %** | 2.318 **+10 %** | 3.552 **+9 %** |
| 0678599eb | 2.678 **+6 %** | 2.335 **+11 %** | 3.584 **+10 %** |

- The commits we skipped are 0aee9a208, 8db6f6141, 5e33f056a, d320347b3 and 4ac0c0e49. Together they are the 4cacdae21 -> 0678599eb
  interval, which is flat.
- ff74a049e and d589df2d7 add about +2 % on cs_closed. That is real at this noise level, but below the 5 % rule.

## Where in a97e368cb (Raven jobs 32443167 and 32443460, 0678599eb with all keys off plus one host-only change)
In 0678599eb, CleanCellLineRows sets the team count to `nlg = DevExeSpace().concurrency()` on the host, and both line kernels loop
`for (ln = league_rank; ln < nline; ln += league_size)`.

| variant | flz_q | cs_closed | clean stage steps | what it changes |
|---|---|---|---|---|
| push 15 | 2.523 (job 32443460) | 2.191 | 3.356 | - |
| 0678599eb | **+7 %** | **+6 %** | **+7 %** | - |
| P1 | +1 % | +5 % | +4 % | nlg = nline on the host (one team per line, each with its own rows block). Job 32443167 had push 15 2.530 / 2.114 / 3.265 and 0678599eb +5 / +10 / +10 % |
| P2 | **+0 %** | **+1 %** | **+1 %** | one team per line, and the rows block indexed by omp_get_thread_num() (scratch-like reuse) |
| P3 | **+0 %** | **+1 %** | **+2 %** | teams and owned blocks as now, but team t takes the contiguous lines [t nline/T, (t+1) nline/T) |

- P1 shows that most of the cost is the strided order. Its leftover cost comes from nline separate rows blocks, which are cold in cache.
- P2 and P3 both restore push-15 speed. P3 is the smaller change, and it keeps the team-owned blocks that the MI300A fix needs.
- Both are bitwise to 0678599eb on flz_q and cs_closed.
- The P3 change is the same in both CleanCellLineSolveRows (line 860 of 0678599eb) and CleanCellLineFactorSolve (line 1054):

```cpp
-    for (int ln = tm.league_rank(); ln < nline; ln += tm.league_size()) {
+    const int lo = static_cast<int>(static_cast<int64_t>(nline)*tm.league_rank()/tm.league_size());
+    const int hi = static_cast<int>(static_cast<int64_t>(nline)*(tm.league_rank() + 1)/tm.league_size());
+    for (int ln = lo; ln < hi; ++ln) {
```

- On a GPU the team count is concurrency / nt, and each team runs its lines one after another either way. So contiguous chunks should
  cost nothing there, but we did not time a GPU.
- If you prefer to leave the GPU order alone, guard P3 with the existing
  `if constexpr (!std::is_same_v<DevExeSpace, Kokkos::DefaultHostExecutionSpace>)` pattern.
- We did not commit the patch. It is in /orion/ptmp/jinma/clean_1115/push16_reg/bis/patch.py (variants P2 and P3).

## Orion vs Raven
- On Orion, all keys off was +7 % on flz_q and +26 % on cs_closed vs push 15 (job 230835, noisier). On Raven, the same 0678599eb arm is
  +5..+7 % and +6..+11 %.
- The commit and the cause are the same either way. We did not repeat the bisect on Orion.

## Files
- Raven, /raven/ptmp/jinma/stag2/p16bis:
  - nodeB.log, nodeP1.log and nodeP23.log;
  - ana.py (the tables);
  - athena_* (the binaries; md5s are in the p16*.out files).
- Orion, /orion/ptmp/jinma/clean_1115/push16_reg/bis: the job scripts and patch.py.
- No restart dumps were left.
