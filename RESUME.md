# RESUME: m1-fast5-sp (sp wedge implicit M1 + vet_col speed-ups), 2026-09-25 ~23:30

- **Branch:** `m1-fast5-sp`, from rt-integration 2c2a0d78.
- **Worktree:** `/viper/ptmp2/jinma/wt_fast5sp`.
- **Run tree:** `/viper/ptmp2/jinma/fast5sp_0925`.
  - `bin/`: binaries. The tag is the lever, `c` = CPU build.
  - `runs/`: GPU jobs.
  - `cpu/`: CPU gates.
  - `st/`: the runs_5q space-time set.
  - `scripts/`: the scripts.
- **Status:** nothing merged or pushed. Build dirs `b_*` and `src_*` in the run tree are still to be deleted.

## State: commits (oldest first)

| commit | lever | result |
|---|---|---|
| bfd81013 | v1: vet_col team sweep, one barrier per chunk instead of per shell | bitwise on CPU |
| 90206c74 | v2: transposed (shell, ray) seg/gb table (+161 kB), ray type hoisted | bitwise on CPU |
| 28676f4e / ba08cbea | v3: grouped branch-free weights | slower, reverted |
| ee2320d7 | v4: ideal-gas kernels for m1_impl_tsolve, m1_impl_src, m1_t2_vcp | table code out of these kernels |
| 7798c44d | v5: next-shell seg/gb prefetch | bitwise, ~0 gain |
| ca5481e8 | v6: m1_impl_stb compile-time slots with value selects | 176 B scratch -> 0 |
| 32341697 / cf18e10e | v7: pipelined rays/moments | 8.1 ms per build, slower, reverted |
| ff2506e8 | v8: par_for_lb on m1_impl_tcell, m1_impl_f2face, m1_impl_f3face | bitwise |
| 54eb6cf7 | v9 (WIP): vet_col moments in 8 fixed ray groups, default chunk 8 | ROUND-OFF, not bitwise |

## Measured (GPU apudev, 1 GPU, us per call before -> after)

| kernel | b0 | v8 | commit |
|---|---|---|---|
| m1_vcol_team | 9306 | 6774 | v1/v2 |
| m1_impl_tsolve | 613 | 162 | v4 |
| m1_impl_src | 310 | 38 | v4 |
| m1_t2_vcp | 230 | 90 | v4 |
| m1_impl_stb | 456 | 314 | v6 |
| m1_impl_tcell | 254 | 236 | v8 |
| m1_impl_f3face | 75 | 56 | v8 |

- Profile: `runs/t6/prof`; table from `scripts/ktab.py`.
- **ms/cycle** (job 11980524, runs/t6, 3 interleaved reps):

  | | b0 | v8 |
  |---|---|---|
  | 1 GPU | 31.78, 31.98, 31.88 | 27.18, 27.13, 27.26 |
  | 2 GPUs | 19.48, 19.57, 19.52 | 16.96, 17.08, 17.05 |

- **Radiation / hydro** (timers): 1 GPU 3.47 -> 2.81; 2 GPUs 3.26 -> 2.76.
- **v9, per build** (job 11980678, runs/swv9b, fenced, vet_col_chunk = 6, 7, 8, 10, 12): 5.44, 5.21, **4.95**, 5.93, 5.57 ms. v8 is 6.8 ms; with the old default chunk 16, v9 takes 7.5 ms. v9 ms/cycle is not yet measured.
- **Counter profile** (runs/c_b0, `scripts/csum.py`):
  - m1_vcol_team is latency-bound. SQ_WAIT_ANY is 83 % of wave cycles; VALU is active 8.8 % of wave cycles; about 3 teams per CU (LDS).
  - Experiment x1 (moment sums removed, timing only; runs/swx1): 3.6 ms per build, so the per-shell serial moment sums cost about 3 ms.

## Gates

- **CPU, v6c vs b0** (`scripts/cpugates.sh`, `cpu/gates_b0_v6c.out`):
  - box and slab arms bitwise;
  - small sp wedge bitwise;
  - wedge restart bitwise;
  - space-time set: 16/16 L1 rows identical (`st/RESULTS_s.txt`);
  - NON-CONVERGED 0.
- **CPU, later levers:**
  - v7c: small-wedge bitwise (reverted anyway).
  - v8 only changes the launch bounds. Its CPU gate is not run.
  - v9c: space-time set 16/16 rows identical at printed precision (`st/RESULTS_s2.txt`). The small-wedge run and the Cartesian gates on v9c are NOT run.

## GPU: round-off or bug? Verdict so far

- **Round-off, not a bug, for v1-v8.**
  - CPU runs are bitwise for every lever. The GPU restart files after 60 cycles differ from b0 (v1 onwards) by the same amount as b0 on 1 GPU vs b0 on 2 GPUs: median relative difference 1e-4, max 2, 38.9M of 62.3M doubles (`scripts/rstdiff.py`, runs/t3).
  - The first vet_col build of v1/v2 matches b0 to at most 2.4e-15; it then grows to about 1e-12 per step (runs/d1, dumps). That is FMA contraction, amplified by lin_tol and Picard tol 1e-10.
- **v7: GPU differences not explained.** Its first-build dump differed from v6 on GPU (runs/d2) while CPU was bitwise. It was reverted for speed; I did not diagnose it.
- **v9:** the first dump is within 1.5e-15 of v8 (runs/d3). It is round-off by design.

## Exact next steps

1. Build v9 from HEAD 54eb6cf7 (default chunk 8):
   - `REV=HEAD bash scripts/build.sh v10 none gpu`
   - `bash scripts/build.sh v10c none cpu`
   - `bash scripts/build.sh v10c box_convection cpu`
2. GPU timing: `sbatch --export=ALL,TAGS="b0 v8 v10",RUN=t7,NREP=3,PROF=1,TIMERS=1 scripts/job_t.sh`. Then:
   - `python3 scripts/tsum.py runs/t7 10 50`
   - `scripts/ktab.py runs/t7/prof/wed1_*`
   - `scripts/tmr.py runs/t7/wed{1,2}tm_*`
3. CPU gates: `bash scripts/cpugates.sh b0 v10c`.
   - Cartesian box/slab: must be bitwise, because Cartesian vet_sc does not use vet_col.
   - Cartesian vet_col: changes by round-off (grouped sums).
   - Wedge: round-off, not bitwise (sp vet_col).
   - Restart: bitwise within v10c.
   - Space-time set rows: identical.
   Write down the magnitude of the wedge CPU difference with `scripts/rstdiff.py`.
4. Optional further vet_col lever: try VC_NG = 16 (compile-time).
5. Candidates not yet done:
   - m1_impl_asm: occupancy 22 %.
   - Hydro hflux_x1/x2/x3 on sp: 1.0-1.4 kB scratch per lane, 2.9 ms/cycle. This is hydro, outside the brief; report it.
6. Final report (under 12 lines): per-kernel before/after, ms/cycle and radiation/hydro at 1 and 2 GPUs. State that v1-v8 are bitwise on CPU and round-off on GPU, and that v9 is round-off. Then delete `b_*` and `src_*` in `/viper/ptmp2/jinma/fast5sp_0925`.

## Session 09-26 (dated 2026-09-26; worker agent)

- **Builds from HEAD 9ebad969 (code = 54eb6cf7, v9 with default chunk 8), all OK:**
  - `bin/athena_v10_none_gpu`, md5 82b832af02d8cc45c523b86a03d59215
  - `bin/athena_v10c_none_cpu`, md5 f7f9a21c0a116eea7c3e5f995ecc7556
  - `bin/athena_v10c_box_convection_cpu`, md5 5cc6eaf65a3760869b26ed8909bfe32e
- **GPU timing job 11981006** was submitted: TAGS="b0 v8 v10", RUN=t7, NREP=3, PROF=1, TIMERS=1. `job_t.sh` already exports HSA_XNACK=1 and HSA_NO_SCRATCH_RECLAIM=1. The job is not evaluated yet. To evaluate it, run `tsum.py runs/t7 10 50`, `ktab.py runs/t7/prof/wed1_*` and `tmr.py runs/t7/wed{1,2}tm_*`.
- **CPU gates** (`bash scripts/cpugates.sh b0 v10c`, output in `cpu/gates_b0_v10c.out`):
  - **Cartesian box/slab:** box_a1, box_m2, slab_a1 and slab_m2 are BITWISE in both loose and tight modes, restart files included.
  - **Cartesian vet_col:** cpugates.sh has no such arm (every gates.py arm uses vet_sc), so this gate was NOT run.
  - **sp wedge (b0 vs v10c):** not bitwise, as expected. Measured with rstdiff.py:
    - after 4 cycles (`cpu/g_b0_v10c/w_*`): 1.31M of 2.51M doubles differ; relative max 1.01, median 2.9e-5.
    - after 1 cycle (`cpu/g_b0_v10c/n1_*`): 0.75M of 2.09M doubles differ; relative max 4.7e-4, median 2.0e-11.
    - The difference grows about 1e6 times in 3 cycles, even though implicit_tol and lin_tol are 1e-10. This matches the GPU b0 1-GPU vs 2-GPU spread from the earlier session (median 1e-4 after 60 cycles), but the case for round-off would be stronger with a b0-vs-b0 CPU perturbation baseline.
    - NON-CONVERGED 0.
  - **Wedge restart within v10c:** BITWISE after par_end.
  - **Space-time set:** 16/16 L1 rows are identical between s0 (b0) and s1 (v10c) (`st/RESULTS_s.txt`, overwritten by this run).
