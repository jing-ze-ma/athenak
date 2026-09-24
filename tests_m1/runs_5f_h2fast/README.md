# runs_5f_h2fast: a cheaper `time_scheme = hesdirk2` step

- **Branch.** `m1-h2fast` from rt-integration 9c1a12d7; worktree `/viper/ptmp2/jinma/wt_h2fast`.
- **Build and run tree.** `/viper/ptmp2/jinma/h2fast_0924`:
  - binaries in `bin/`: `base` = 9c1a12d7; `v1` = lever defaults only; `v3` = + stage
    one-pass safety; `v4` = final code (f557ebf6; the later commits change comments only).
    GPU v4 md5 b7b9aebe, CPU box v4 67e2b38f, CPU mpi v4 6d936c55;
  - `gpu*/` GPU jobs, `cpu/` CPU gates, `rw/<arm>/` radwave runs.
- **Environment of every GPU job:** apudev, 1 x gfx942 (2 for the bitwise gate),
  `HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1`.
- **Box used for timing:** the 3-D He box (84x104x104, 4 blocks, plm + vimp), 120 cycles,
  ms/cycle over cycles 20-120. The arms are interleaved and repeat 2 runs them in reverse
  order (`scripts/sweep.sh`).

## Result

hesdirk2 / be per step, on 1 GPU, same job, interleaved, 2 repeats (job 11960025):

| closure | be (ms/cycle) | hesdirk2 before (9c1a12d7) | ratio before | hesdirk2 after (v4 defaults) | ratio after |
|---|---|---|---|---|---|
| Eddington | 40.08, 39.54 -> 39.8 | 74.02, 73.27 -> 73.6 | **1.85** | 55.70, 51.83 -> 53.8 | **1.35** |
| vet_sc (full) | 45.78, 47.23 -> 46.5 | 83.71, 86.55 -> 85.1 | **1.83** | 62.36, 62.49 -> 62.4 | **1.34** |

The same final set, measured one job earlier (job 11959729, arm `h2lt9`: v3 +
`time2_lin_tol = 1e-9`, which is what the v4 defaults give), came out at 50.5 / 52.0 ms
(Eddington) and 59.2 / 58.9 ms (vet_sc). Against be at 40.8 / 46.9 ms in that job, both
ratios are **1.26**. The apudev node was shared with other agents' jobs, so single runs
scatter by up to about 4 ms; one run in job 11959401 took 256 ms/cycle.

**Why the brief's "1.48x / 1.57x" became "before 1.85x / 1.83x".** The brief's ratios
(runs_3x) compare against be without the accel levers. be has had the levers by default
since m1-accmerge (40 ms), while hesdirk2 kept none of them (74 ms).

## Levers (all `<rad_m1>`; a be run reads none of the new keys)

| lever | status | Eddington ms/cycle | vet_sc ms/cycle |
|---|---|---|---|
| hesdirk2 before (9c1a12d7) | | 73.6 | 85.1 |
| (1) the runs_4a_accel levers default on for hesdirk2 too: `implicit_fast_kernels`, `implicit_vimp_fold`, `implicit_one_pass = 8`, `implicit_predictor_order = 2` | KEPT | 57.9 (h2EL) | 68.8 (h2VL) |
| (1b) `time2_one_pass_safety = 30`: the one-pass safety factor of the stage solves only | KEPT, needed for G1 | 58.5 (job 11959729) | 69.2 |
| (2) `time2_lin_tol_fac = 10`: stage-solve Krylov tolerance = 10 x `implicit_lin_tol` | KEPT | 53.8 | 62.4 |
| stage-2 predictor from stage 1 (`d2 = d1/2` + history of `d2 - d1/2`) | DROPPED: Krylov iterations per solve 21.2 -> 21.8, Newton fallbacks 26 -> 147 (job 11959226, h2sEn) | 56.4 vs 56.8 | 66.5 vs 69.4 (noise) |
| `implicit_one_pass = 16` | DROPPED: 57.2 vs 57.4 (job 11959401) | | |
| `implicit_one_pass = 4` | not taken: G1 passes (1.86), no speed gain | 57.9 (job 11959729) | |

**(1) The accel levers.** runs_4a_accel had already validated them under hesdirk2. The
per-solve-kind one_pass q and the stage-2 predictor `ipred2` exist already. Only the
default rule changed, in `ImplicitInit`: `ldef` accepts `time_scheme = be | hesdirk2`. It
follows the same restart rule: a restart whose file lacks a key keeps the old value.

**(1b) Why the stage safety is needed.** With one_pass 8 at safety 3, G1 at the gate
tolerance (tol 1e-11) drops to a median order of 1.74 at (P, tau) = (100, 1e3), for both
closures (`RESULTS_g1_arms.txt`, arm n1).
- **The cause.** An accepted one-pass solve leaves a Picard error of about
  tol / safety. The two-pass test leaves about q tol, with q ~ 1e-4 to 2e-3. The radwave
  (amplitude 1e-5) sees that difference at nt >= 1024.
- **Isolated arms at the two stiff corners** (`rw/iso`, lowest single order). Each of
  predictor_order 2 and one_pass alone perturbs the order: one_pass 0 gives 1.52 and
  predictor_order 1 gives 1.30. With both off it is 1.96. So this is a solver-tolerance
  effect, as runs_4a_accel found.
- **The fix.** `time2_one_pass_safety = 30` applies to the stage solves only, so be keeps
  its 3. It brings every G1 median back to >= 1.94 (arm sf30). It costs about 0 ms/cycle
  in the box: Picard 1.20 -> 1.35 passes per solve, but the ms/cycle is unchanged within
  noise.

**(2) The stage-solve linear tolerance.** The inputs set `implicit_lin_tol = implicit_tol/100`
(1e-10 / 1e-8). The stage solves now take 10x that: 1e-9 in the box, 1e-12 in the radwave.
- **Box.** Krylov iterations per linear solve fall from 21.2 to 14.0 (Eddington) and from
  22.1 to 14.8 (vet_sc). The cost falls by 4-7 ms/cycle.
- **G1.** It passes together with safety 30 (arm sf30lt12 and the final v4 defaults,
  below). Alone, with safety 3, it does not: arm lt12 gives 1.76.
- **Implementation.** `impl_lin_tol` is swapped at the top of `ImplicitSolve` for a stage
  solve and put back at its single exit.

**Not a lever: the vet_sc opacity-caching fix.** It had landed before this branch:
`Time2VetStart` has the fused fast path, which saves and copies the stage opacities.

### Profile (rocprofv3, job 11959228, 40 cycles, `gpu/P*/split.txt`; v1 binary: levers on)

be at 42.2 ms/cycle against hesdirk2 at 59.8. The +17.6 ms is:
- **Krylov, +9.1 ms.** The Krylov part is 26.6 ms against 17.4, as 51 against 35 inner
  iterations per step:
  - ImplicitStencilOp +3.1;
  - M1PCRX +2.7;
  - vector/reduction +1.3;
  - halo +0.7.

  Lever (2) targets this.
- **Per-pass Picard work, +6.5 ms.** The second solve pays its own opacity, assembly, vimp
  build, stencil build and transverse terms:
  - Opacity 0.9 ms per solve;
  - ImplicitSolve#10 0.37 ms per pass;
  - VimpBuild 0.47;
  - StencilBuild 0.23.

  Every one of them depends on the stage state: gas T, opacities, velocity. Nothing is
  unchanged between the stages that could be reused. The PCR line preconditioner is rebuilt
  from the stencil on every apply, so it keeps no factorisation to reuse.
- **Hydro side, +2.0 ms.** One more HydroConToPrim (0.7 ms) and BoxConvBC x2 (0.5 ms)
  after the stage-1 solve. The next Heun stage needs these.

## Gates

| gate | result | evidence |
|---|---|---|
| be bitwise, CPU: `gates/bitwise_off.sh` base vs v4 box (box 1/2o/4o ranks, slab 1/2o, slab vet_sc 2, slab_nd 2 with no keys) | PASS 7/7: hst, rst and every row identical | `RESULTS_bitwise_off.txt` |
| be bitwise, GPU (job 11960026, 2 GPUs): He slab 1 rank, box3d_nd 2 ranks, box3d_nd + vet_sc 2 ranks, 9c1a12d7 vs v4 | PASS 3/3 BITWISE (10/10 files each) | `RESULTS_gpu_gate.txt` |
| `tst/test_suite/rad_m1` (ATHENAK_M1_DATA = faces_0924/m1data, snapshot of v4) | 3 passed (420 s) | `cpu/pytest.log` |
| `gates/gates.py` run + eval (v4 box) | GATES: PASS, 14/14 pairs | `RESULTS_gates_py.txt` |
| hesdirk2 G1 radwave, v4 defaults, 12 cases x 2 closures, tol 1e-11 | PASS: median order >= 1.94 (Eddington), >= 1.87 (vet_sc, at (100,10); base 1.95). Base: 1.89 / 1.92. Lowest single order 1.58, at (100,1e3) and (10,1e3). e1024 is 0.66x-1.02x of base. | `RESULTS_g1_final.txt` |
| stiff P=100, tau=1e5, 40 periods (N32 nt24, N64 nt32/48) | same as base to 3 digits (end 0.719 / 0.800 / 0.801, exact 0.799; NC 1 / 0 / 0 in both) | `RESULTS_stiff.txt` |
| He slab 2-D hesdirk2, v4 defaults: Eddington 1 rank and vet_sc 2 ranks to 200 s and 1000 s | NON-CONVERGED 0, stage fallbacks 0 (1 BE step, the first) | `cpu/Zh2*` |
| forced fallback (`time2_dbg_fail = 50`) | 1 fallback, 2 BE steps, runs through, NC 0 | `cpu/Zh2Efail` |
| restart hesdirk2 (t=100 -> 200): Eddington 1 rank, vet_sc 2 ranks | PASS: bin and rst identical, 100/100 hst rows identical (rt_profile.bin differs, as always across a restart: runs_4j) | `cpu/RZh2*` |
| 3-D box, GPU, 120 cycles, every arm | NON-CONVERGED 0 | `gpu5/*/run.log` |
| implicit_op_check under hesdirk2 (v4 defaults, vimp_fold on), solves 1-5 (BE step + 2 stage steps): box 1 rank, box 2 ranks halo_mpi, box vet_sc, slab 1, slab 2 halo_mpi, slab vet_sc 2 | PASS 5/5 solves in all 6 runs, worst dy/row 9.7e-16 | `RESULTS_opcheck_h2.txt` |

**Flag.** The per-cell gas Newton falls back to the bracketed root more often under
hesdirk2 with the levers:
- in 120 box cycles: 44 (Eddington) and 15 (vet_sc) in about 2.6e8 cell-passes;
- before: 8 (Eddington) and 0 (vet_sc).

be with its defaults has 0-2 (runs_4j reported the same kind of count). No implicit-solve
fallback fired.

## Files

- `run_radwave.py`, `time_table.py`, `med_table.py` (`EXCL=vetf|edd_` splits the closures),
  `stiff_table.py`, `rw.sh`, `lists/`: the radwave runner (runs_4a_accel copy, paths to
  h2fast_0924).
- `scripts/`: `sweep.sh` (GPU timing / rocprof), `gate_gpu.sh`, `cpu_h2.sh` (slab gates +
  restart), `opchk_h2.sh`, `pytest_m1.sh`, `gates_run.sh`, `build.sh`, `gbuild.sh`, `gsum.py`,
  `gcmp.py`.
- `RESULTS_*.txt`: the tables above.

## HANDOVER

- **Merge candidate.** m1-h2fast is a candidate for rt-integration.
- **What changes, and for which runs.**
  - hesdirk2 runs change: new defaults, round-off to solver-tolerance level.
  - be runs are bitwise, on the CPU and the GPU.
  - A restart of an older hesdirk2 file keeps the old behaviour. It gets no levers, safety
    falls back to `implicit_one_pass_safety`, and lin_tol_fac is 1.
- **Cost now.** hesdirk2 costs 1.26-1.35x be per step on 1 GPU, for both closures.
- **What is left: the second solve's fixed Picard work.** About 6.5 ms/cycle of assembly,
  vimp, stencil and opacity per pass, plus the extra HydroConToPrim. The rest is Krylov
  iterations at the looser tolerance.
- **Next ideas, not tried:**
  - fuse the Opacity task into the solve's start-T kernel. It would help be too, and must
    stay bitwise;
  - the admissibility reduction could share the write-back kernel (about 0.1 ms).
- **G1 at tol 1e-11 sits at the solver-noise floor.** At (P, tau) = (10, 1e3) and
  (100, 1e3), any change of the starting iterate or of the pass count moves single orders
  between 1.3 and 2.6. Judge by the median, or run at tol 1e-12 (`lists/g1_t12.txt`).
