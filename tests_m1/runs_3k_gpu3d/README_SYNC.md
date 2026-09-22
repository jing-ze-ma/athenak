# Speed-up #3: fewer host syncs and launches in the implicit-M1 BiCGStab

Date 2026-09-22/23, viper. Base: HEAD 1d22a118 (rt-integration, includes PCR 1715190b and VET).
NOT COMMITTED. Patch: `bench/m1_sync_0922/sync.patch` (`git apply --check` clean on 1d22a118;
files `src/rad_m1/rad_m1.{hpp,cpp}`, `src/rad_m1/rad_m1_implicit.cpp`). Snapshot and builds:
`bench/m1_sync_0922/{snap,athena_cpu,athena_gpu}` (GPU md5 15bb8dcc...).

## What changed

A new input key, `<rad_m1>/implicit_bcg_sync`, with levels 0, 1 and 2. The default is 0, the
original loop, which stays unchanged.

- **Level 1: `ImplicitBiCGStabFused`.** The recurrence, breakdown and restart rules,
  true-residual check and fallback are all the same as the original.
  - The kernel that updates x and r also returns max|r| and (rhat, r_{k+1}). The original
    used three kernels and two syncs for this.
  - (t,s) and (t,t) come from one kernel.
  - r0 and the restart kernels return (r,r) and max|r|.
  - The p and s updates also stage `M1_IW_TR`. `ImplicitPrecond(-1, z)` then skips its
    copy kernel.
  - Reductions use one flattened RangePolicy and a 3-sum + 1-max reducer. Across ranks
    this is one MPI_Allreduce, with a user op on a 4-Real type.
- **Level 2 = level 1, plus rhat.v reduced into a device scalar (1 rank only).** The
  s-kernel reads alpha on the device. The host gets rhat.v, exactly, from the (t,s)
  reduction and only then applies the breakdown test. If that test fails, s, z and t are
  discarded; x has not been touched. With more than 1 rank, level 2 acts as level 1.
- **Shared end of solve.** The fallback, output and statistics tail is shared by both
  loops (`ImplicitBiCGStabEnd`, a pure move).
- **Style.** cpplint shows 8 errors in the 3 files both before and after the patch, so the
  patch adds none. No line is over 90 columns.
- **Other options.** Checking convergence less often was not needed: max|r| rides in the
  same reduction as rho, so it costs no sync. A pipelined BiCGStab was not tried (see below).

## Gate 1: CPU, 2-D seeded slab, 200 steps (`bench/m1_sync_0922/cpu/`)

Setup: thomas, 1 rank, plus a 2-rank run with meshblock nx2=16.

- **Levels 1 and 2 are BITWISE identical to level 0** on the CPU, for both 1 and 2 ranks.
  - The bin files differ only in the echoed input line, and the hst files are identical.
  - Counts: 805 Picard passes, 6347 inner iterations, 0 breakdowns. On 2 ranks all three
    levels give 6358 inner iterations.
  - The +/-1 ulp controls (rad_flux_inner) give 6345 and 6335 inner iterations, with hst
    differences of 3.4e-7 and 5.1e-7.
- **Blocking reductions per inner iteration:** 4.19 at level 0, 3.19 at level 1, 2.19 at
  level 2.

## Gate 2: GPU, apudev, 3-D box with `implicit_line_solver = pcr`, 200 steps

- **Jobs:** 11943323 (1 GPU, timing), 11943324 (2 GPUs, and rocprofv3 with nlim = 20),
  11943419 (the +/-1 ulp controls).
- **Method:** same binary for every arm. Arms were interleaved and repeated, with the
  second repeat in reverse order.
- **Files:** logs in `bench/m1_sync_0922/runs/<arm>/run.log`. Tables come from
  `summarize.py` (ms/cycle = elapsed from cycle 50 to cycle 200), `gapattr.py` and
  `gapsum.py`.

| box | GPUs | level | ms/cycle (rep a, b) | zone-cyc/s | inner its/step | blocking red./step | red./its |
|---|---|---|---|---|---|---|---|
| 84x104x104 | 1 | 0 | 267.1, 264.5 | 3.30e6, 3.38e6 | 110.54 | 449 | 4.06 |
| 84x104x104 | 1 | 1 | 244.4, 243.1 (**1.09x**) | 3.65e6, 3.67e6 | 110.59 | 339 | 3.06 |
| 84x104x104 | 1 | 2 | 241.8, 243.7 (1.09x) | 3.70e6, 3.67e6 | 110.59 | 228 | 2.06 |
| 84x52x52 | 1 | 0 | 175.3, 175.1 | 1.28e6 | 99.44 | 405 | 4.07 |
| 84x52x52 | 1 | 1 | 149.5, 146.1 (**1.18x**) | 1.49e6, 1.53e6 | 99.34 | 305 | 3.07 |
| 84x52x52 | 1 | 2 | 147.0, 144.0 (1.20x) | 1.53e6, 1.56e6 | 99.34 | 206 | 2.07 |
| 84x104x104 | 2 | 0 | 243.1, 244.7 | 3.60e6, 3.65e6 | 110.54 | 449 | 4.06 |
| 84x104x104 | 2 | 1 | 214.7, 217.9 (**1.13x**) | 4.15e6, 4.10e6 | 110.53 | 339 | 3.06 |
| 84x104x104 | 2 | 2 (= 1) | 213.7, 213.6 | 4.17e6 | 110.53 | 339 | 3.06 |

- **Stability:** NON-CONVERGED is 0 in every arm. The Picard mean is 4.33 to 4.39 in all arms.
- **Timing noise:** the level-0 +/-1 ulp arms ran at 263.9 and 274.0 ms/cycle, so level-0
  timings vary by about 4 % between runs.
- **Round-off (GPU runs are deterministic: rep a == rep b bitwise).** Measured as the
  largest hst column difference (normalised by the column max), level 1 or 2 against
  level 0, next to the +/-1 ulp controls of level 0:

| box | level 1/2 vs level 0 | +1 ulp control | -1 ulp control |
|---|---|---|---|
| 84x104x104 | 3.35e-4 | 3.74e-4 | 5.67e-4 |
| 84x52x52 | 7.66e-4 | 4.70e-4 | 4.10e-4 |

  - These largest differences are all in the 2-mom and 3-mom columns, net transverse
    momenta that are about 0. KE differs by about 1e-9, and user.hst by 1.2e-8 to 1.9e-8
    (controls 0.9e-8 to 2.0e-8).
  - Inner iterations per step on the half box: 99.34 at levels 1/2 and 99.44 at level 0,
    against 99.08 and 99.17 for the controls.
  - On the GPU the change is therefore round-off, not bitwise.
- **Trace (rocprofv3, nlim = 20, window after ProblemGenerator):**

| arm | launches/cycle | BiCGStab + precond-copy kernels | idle ms/cycle (%) | idle before BiCGStab kernels | idle before bvals pack/unpack |
|---|---|---|---|---|---|
| full box, 1 GPU, level 0 | 3890 | 48.2 ms | 117.9 (33.4 %) | 20.6 | 60.9 |
| full box, 1 GPU, level 1 | 3211 | 26.6 ms | 128.5 (37.4 %) | 17.8 | 65.4 |
| full box, 1 GPU, level 2 | 3211 | 27.0 ms | 126.2 (37.0 %) | 11.9 | 66.3 |
| half box, 1 GPU, level 0 | 3317 | 34.4 ms | 101.3 (45.1 %) | 24.0 | 54.5 |
| half box, 1 GPU, level 2 | 2724 | 10.5 ms | 87.4 (47.0 %) | 9.7 | 55.0 |
| full box, 2 GPUs, level 0 (rank 0) | 3891 | - | 149.4 (44.9 %) | 30.6 | 86.0 |
| full box, 2 GPUs, level 1 (rank 0) | 3221 | - | 139.8 (47.0 %) | 20.9 | 89.1 |

- **Where the gain comes from.** Most of it is kernel time: the BiCGStab kernels fall
  from 48 to 27 ms/cycle, and the old MDRange dot products were slow. Launches fall by
  17 %. Idle time before the BiCGStab kernels falls by 3 to 14 ms/cycle.
- **Why the idle fraction rises.** The busy time shrank faster than the idle time. The
  "other" idle in the profiled full-box level 1/2 arms is also 9 to 12 ms/cycle higher. It
  comes from the gap before the hydro-side kernel `RemoveGravEtot` (13 to 21 ms/cycle across
  arms, with no BiCGStab kernel in between), so it is not caused by the patch.
- **Level 2 adds almost nothing over level 1** (1 % or less, inside the noise). The
  blocking reduction it removes is followed straight away by a halo exchange, and that
  exchange syncs anyway (below).
- **Recommendation:** `implicit_bcg_sync = 1` for every GPU M1 input. The gain is the same
  with level 2, which only helps on 1 rank.

## What the remaining idle time is (next target, NOT in this patch)

About half of all idle time comes just before `MeshBoundaryValuesCC::PackAndSendCC` and
`RecvAndUnpackCC`:
- 61 to 66 ms/cycle on 1 GPU and 86 to 89 ms/cycle on 2 GPUs;
- 162 to 171 us per pack call;
- 2 halo exchanges per inner iteration.

**Cause.** Both functors are larger than 32 kB. Kokkos HIP therefore launches them through
`hip_parallel_launch_global_memory`, and that path calls
`HIPInternal::stage_functor_for_execution`. That function starts with `hipStreamSynchronize`
(`kokkos/core/src/HIP/Kokkos_HIP_Instance.cpp` l. 322). So every pack and every unpack launch
is a full host sync.

On top of that, PackAndSendCC calls `Kokkos::fence()` before its MPI sends.

**Possible fixes.** Shrink the captured buffer structs below 32 kB (or below 4 kB), or give
the 1-variable Krylov exchange its own small pack/unpack kernels. Estimated gain: tens of
ms/cycle, not measured.

**Smaller item.** ImplicitOffDiagOp captures three DualViews, so it launches through
constant memory. The gap before it is 12 us x 290 per cycle, about 3.5 ms/cycle. Capturing
the `.d_view`s would bring it under the 512-byte local-memory limit. Its kernel time is now
54.5 ms/cycle at 1d22a118 (it was 40.6 at 886bd677).

A pipelined BiCGStab (1 sync per iteration) was not tried. The syncs left in the BiCGStab
loop cost at most 12 to 21 ms/cycle, less than the halo-exchange syncs.
