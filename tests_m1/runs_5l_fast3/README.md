# runs_5l_fast3: cheap levers on the implicit M1 + VET cycle (m1-fast3)

- **Date:** 2026-09-24/25, viper. **Branch:** `m1-fast3` from rt-integration d2572b4f
  (= `base`), worktree `/viper/ptmp2/jinma/wt_fast3`. **Run tree:** `/viper/ptmp2/jinma/fast3`
  (sources `src_base`, `src_n3` = f3ee5a0b + the 2-rank switch fix, `src_n4` = this commit; binaries in `bin/`,
  md5 in `build_*.log`; job logs `runs/log.out.<job>`; build dirs deleted).
- **Cases** (as in /viper/ptmp2/jinma/prof_m1_0924): **box** = He box `inp/box_vsc.athinput`,
  closure vet_sc, 84x104x104 on 1 GPU, the w2 box 84x208x104 on 2 GPUs (weak);
  **wedge** = sp He wedge `inp/hewedge.athinput` 96x128x128, closure vet_col (2 GPUs strong).
  Both start from their input ICs (no restart involved).
- **Timing:** apudev MI300A, `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1`, 60 cycles,
  ms/cycle over cycles 10-50 (`tsum.py`), 3 interleaved repeats per arm, same binaries.
  Final jobs: **11968166** (1 GPU, `job_t5.sh`, n3), **11968167** (2 GPUs, `job_t6.sh`,
  n3) and **11969142** (2 GPUs, `job_t7.sh`, n4 = team_red also in the 2-rank parts);
  earlier jobs 11967747 / 11967748 (n2 binary, split_red, krylov_dev, pipe, lin_tol arms)
  agree with them to 0.3 ms. Repeat spread <= 1.1 ms (1 GPU), <= 0.3 ms (2 GPUs).

## What changed (all switches default off unless stated)

| lever | change | default | results |
|---|---|---|---|
| 2 (inert conduction) | `AddIsotropicHeatFluxRadiative` returns before its 3 flux kernels when `BuildRadWeights` proved every weight 0, no wall flux, plain explicit path | **on** (no switch) | bitwise |
| 5 (Picard/opacity) | `RosselandTable`: bisection instead of a linear scan for (i, j) | **on** | bitwise (same index for any non-decreasing table) |
| 1 (fused stencil + dots) | `implicit_op_team_red`: stencil one cell per thread in teams of 256, `team_reduce` into per-team partials, then a small reduction | off | round-off (sum order) |
| 1 | `implicit_op_split_red`: stencil par_for + separate read-only reduction (also at 2+ ranks) | off | bitwise at 1 rank; round-off at 2+ |
| 5 | `implicit_eos_cache_check_every = N`: the eos-cache accuracy check (a measurement, nothing reads it but the final report) on every N-th cycle | 1 (= old) | bitwise state |

`implicit_op_team_red` covers both operators: `ImplicitStencilOp` (1 rank) and the
overlapped interior + shell `ImplicitStencilOpPart` (2+ ranks, per-team partials, then a
small asynchronous reduction into the pinned host sums). `implicit_op_split_red` no longer
turns on with team_red at 2+ ranks (it measured 2-3 % slower there). On a host backend the
team size is 1 (Serial allows no more), so CPU gates exercise the same code.

## Timings (ms/cycle, mean of 3; jobs 11968166 / 11968167, 2-GPU team rows 11969142)

| arm | box 1 GPU | box 2 GPUs | wedge 1 GPU | wedge 2 GPUs |
|---|---|---|---|---|
| base d2572b4f | 55.26 | 62.70 | 41.96 | 23.87 |
| n3 defaults (levers 2 + bisection) | **50.82** (-8.0 %) | **58.41** (-6.8 %) | 41.72 | 23.79 |
| + op_team_red | **49.70** | **55.80** (n4; base 62.53, n0 58.14 in the same job) | **36.65** (-12.7 %) | **23.11** (n4; n0 23.70) |
| + eos_cache_check_every=10 alone | 49.86 | - | n/a (no eos cache) | - |
| + team_red + krylov_dev=2 | 49.92 | - | 42.01 (loses the team gain) | 23.76 |
| + team_red + krylov_dev + check_every=10 | **48.50** (-12.2 % vs base) | **54.16** (n4, -13.4 %) | - | - |

Earlier arms (jobs 11967747/48, n2 binary, same base within 0.3 ms):
split_red box 51.10 / wedge 41.00 (1 GPU), box 57.83 / wedge 24.50 (2 GPUs);
krylov_dev box 50.03 / wedge 42.10 (1 GPU); krylov_pipe (2 GPUs) box 58.02 / wedge 24.70.

**Profiles** (rocprofv3, cycles 10-40, `prof.py`; runs/t3/p*_ns, runs/t5/p*_nt):
- wedge 1 GPU, matvec+dot per cycle: split 7.21 ms (stencil 36.5 us + read-only
  4-value `parallel_reduce` **134 us**) -> team 2.64 ms (**47.4 us** per fused call).
  The slow part was never the stencil: it is the flat Kokkos multi-value reduction on
  this shape (the profile's fused kernel 189 us, the split read-only one still 134 us).
- box 1 GPU: split 65 + 39 us -> team 77 us per call; matvec+dot 7.97 -> 6.73 ms/cycle.
- box 2 GPUs (rank 0): matvec+dot 9.42 (split, job 11967748) -> 8.95 ms/cycle (team parts,
  11969142); the prof_m1_0924 fused parts were 11.13 ms (71 us per interior part).
- Idle after a reducing op (lever 3): box 2.8-3.0 ms, wedge 1.7-2.1 ms per cycle with the
  default host-synchronous BiCGStab. `implicit_krylov_dev=2` (device-resident
  convergence test) saves ~1 ms on the box only when combined with the other levers and
  nothing on the wedge; `implicit_krylov_pipe` at 2 GPUs: box -0.4 ms, wedge +0.9 ms.
  Not worth a default change.
- Remaining big items (not touched): ConsToPrim 5.7 ms (box, general EOS),
  VetColBuildTeam 6.7 ms (wedge), M1PCRX preconditioner 4.5-5.2 ms (m1-precond agent).

## lin_tol (lever 4): keep 1e-10

| test | 1e-10 | 1e-9 | 1e-8 |
|---|---|---|---|
| box 1 GPU ms/cycle (split_red arm, job 11967747) | 51.10 | 45.82 | 41.45 |
| wedge 1 GPU ms/cycle | 41.00 | 35.89 | 32.40 |
| CPU He box (4 ranks, 30 cyc) vs tight ref, max rel hst dyn | 3.8e-10 | 1.1e-8 | 8.5e-8 |
| CPU He slab (100 cyc) vs tight ref, max rel hst dyn | 8.0e-10 | 3.9e-9 | 7.8e-7 |
| radwave (hesdirk2, vet, Picard tol 1e-8) L1 at nt=128, P/tau = 1/10 | 3.96e-4 | 3.94e-4 | **1.01e-2** |
| radwave P/tau = 10/1000 | 2.58e-3 | **5.31e-2** | **9.78e-1** |
| radwave P/tau = 100/1000 | 2.26e-3 | **1.65e-2** | **8.91e-1** |
| T-S4 wedge atmosphere (vet_col, n=64, 800 cyc, CPU) L1 vs exact | 2.473e-4 | 2.473e-4 | 2.473e-4 |

All 0 non-converged. The steady sp atmosphere does not care, but the stiff radiation waves
lose accuracy by 20x at 1e-9 and are wrong at 1e-8 (the order-2 convergence is gone). The
default stays 1e-10; a production that is steady-state dominated could use 1e-9 per input
(-10 % box, -14 % wedge), knowing that transients are resolved worse.
Files: `fast3/cpu/lt.txt` (`ltol_cpu.sh`), `fast3/rw/RESULTS.txt` (`rw.sh`, runs_4a_accel
driver), `fast3/ts4/result.txt` (`ts4_run.sh`, runs_5d_vetcol/ts4v.py).

## Gates

- CPU regression set, base vs n2 (= f3ee5a0b; the default path is unchanged since) defaults (`fast3/cpu/g0.txt`; box1, box2o, box4o, slab1,
  slab2o, slabnd2, slabv2): hst, rows and restarts **bitwise**; restart bitwise True.
  n3 differs from n2 only in the team_red paths (default path identical).
- `tst/test_suite/rad_m1` with ATHENAK_M1_DATA (n2 CPU): 3 passed (`fast3/cpu/pytest.txt`).
- implicit_op_check matrix (`tests_m1/gates/opcheck_matrix.sh`, n3 CPU, 144 runs):
  n3 and n4 defaults: 0 FAIL; with implicit_op_team_red (n3: 1-rank op; n4: also the
  overlap parts): 0 FAIL, reductions within 1.25e-14 of the reference
  (`fast3/cpu/opc{,4}_{def,team}.txt`).
- GPU restarts (data after the parameter block), base vs n3 defaults: box and wedge,
  1 and 2 GPUs **bitwise**; runs are deterministic (repeat 1 vs 2 bitwise).
  check_every=10 vs default: bitwise. n4 team_red at 1 GPU = n3 team_red, bitwise.
- team_red at 1 GPU vs default after 60 cycles: round-off grown by the flow; the same
  size as the accepted krylov_dev difference (wedge: 66 % of doubles differ > 1e-6 for
  both; box: 33 %; median relative difference 8e-8 box, 7e-5 wedge, the wedge's
  near-zero velocities dominate). Inner iteration means 17.10 -> 17.23 (box 1 GPU),
  17.69 -> 17.64 (box 2 GPUs), 9.38 -> 9.41 (wedge 1 GPU), 9.43 -> 9.39 (wedge 2 GPUs);
  0 non-converged everywhere.

## Recommendation

- Keep on (already default): inert conduction skip, bisection (bitwise, -4.4 ms box).
- `implicit_op_team_red = true` 1 GPU: wedge -5.1 ms (-12 %), box -1.1 ms; 2 GPUs: box -2.3 ms,
  wedge -0.6 ms. Round-off only.
  Candidate for default-on (user decision: it changes bits).
- `implicit_eos_cache_check_every = 10` in box inputs (-1 ms, bitwise state).
- lin_tol: keep 1e-10.
