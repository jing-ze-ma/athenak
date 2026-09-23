# Speed-up #4: the halo-exchange host stalls (PackAndSendCC / RecvAndUnpackCC)

Date 2026-09-23, viper. Base: HEAD dccdf506. NOT COMMITTED.
- Patch: `bench/m1_halo_0923/halo.patch` (md5 050ee0a5). It passes `git apply --check`
  on c02ddd34 as well.
- Files changed: `src/bvals/{bvals.hpp, bvals.cpp, bvals_cc.cpp, bvals_fc.cpp,
  prolongation.cpp, prolong_prims.cpp, flux_correct_cc.cpp, flux_correct_fc.cpp,
  flux_seam_cc.cpp}` and `src/rad_m1/rad_m1_implicit.cpp`.
- Builds and gate runs are in `bench/m1_halo_0923/`: binaries in `bin/`, scripts and runs in
  `gate/`, trace tools in `tools/`.
  - `base` = dccdf506.
  - `v1` / `new` = the first patch (no split pack).
  - `new2` = the final patch.

## Mechanism (Kokkos 4.6.02, HIP)

Kokkos HIP chooses how to launch a kernel from `sizeof(driver)`. The logic is in
`DeduceHIPLaunchMechanism`, `kokkos/core/src/HIP/Kokkos_HIP_KernelLaunch.hpp` l. 124-182. The
thresholds are in `Kokkos_HIP_Instance.hpp` l. 52-54.

| driver size | launch path | host cost per launch |
|---|---|---|
| < 512 B (`ConstantMemoryUseThreshold`) | kernel argument (local memory) | none |
| 512 B to 32 kB (`ConstantMemoryUsage`) | constant memory | `hipEventSynchronize` on the previous constant-memory launch (l. 520-522) |
| >= 32 kB | global memory | `stage_functor_for_execution`, which begins with `hipStreamSynchronize` (`Kokkos_HIP_Instance.cpp` l. 322). This is a full drain of the stream. |

With `HintLightWeight`, a kernel-argument launch is allowed up to 4 kB (`KernelArgumentLimit`).

**Why the halo kernels went through the global path.** Every bvals kernel captured
`sendbuf[56]` / `recvbuf[56]` by value, and each array is about 31.5 kB.

**Trace at base** (rocprofv3, M1 3-D box, 1 GPU):
- Global-memory launches per cycle: `PackAndSendCC` 276.4 and `RecvAndUnpackCC` 276.4, and
  nothing else in the M1 step.
- Constant-memory launches: `ImplicitOffDiagOp` 267/cycle, plus 86 others.

**Trace of prod4 at base:** 20 global-memory launches per cycle. They are
`PackAndSend/RecvAndUnpack{CC,FC,FluxSeamCC}`, `PackAndSendFluxFC`, `SumBoundaryFluxes` and
`AverageBoundaryFluxes`.

**A second mechanism, found while fixing the first.** Once the syncs were gone, a gap of
about 142 us was still left before every pack, with no host call in it. The launch was issued
4.5 us after the preceding `ImplicitHaloCopy` launch and started 302 us after it (median).
- The pack is the only kernel in the halo sequence with a private segment: 704 B at base and
  496 B after v1. The seam/pole branch keeps per-thread arrays (stencil indices, 3
  SeamXform).
- A scratch dispatch that follows scratch-free dispatches stalls the queue while ROCr
  re-attaches scratch.
- This is confirmed by `HSA_NO_SCRATCH_RECLAIM=1`, which removes the gap. See the timing
  tables below.

## The patch

1. **`MeshBufferDv` / `SendBufDv()` / `RecvBufDv()` (bvals.hpp/.cpp).**
   - This is a device View of the 56 buffers' index ranges, sizes and (pointer, row-stride)
     handles on `vars` and `flux`.
   - It is rebuilt by `SyncBufferDv()` at the end of `InitializeBuffers()`. It is also
     rebuilt whenever an accessor finds, by host-side pointer compare, that a buffer View
     was reallocated.
   - All 26 bvals kernels now capture this View instead of the host arrays, with identical
     indexing. `MPI_Test` on the host now uses `recvbuf[n]`. The kernel text is otherwise
     unchanged.
2. **`BvalsTeamFor()`.** Every bvals team kernel is launched with `HintLightWeight`, i.e.
   the functor is passed as a kernel argument.
   - A `static_assert(sizeof(functor) <= 3072)` stops a functor that grows from falling back
     silently to the global path.
3. **Lazy fence in the MPI send loops.** The `Kokkos::fence()` before `MPI_Isend` now runs
   only before the first off-rank send, in all 5 senders. On 1 rank there is nothing to wait
   for, because the unpack runs on the same stream.
4. **`PackAndSendCC` compiled twice.** The body is a generic lambda with a compile-time seam
   flag.
   - Cubed sphere / polar runs launch the full body.
   - Every other run launches the body without the seam branch. There `do_cs` and `do_pole`
     are false anyway, so the arithmetic is the same.
   - Still 1 launch per call. The non-seam pack has private segment 0 in the trace (5986 of
     5986 launches). The cs pack keeps 496 B.
5. **`ImplicitOffDiagOp` captures `.d_view`s instead of whole DualViews.** It still lands on
   the constant-memory path (the functor is still > 512 B), but the host syncs it adds fell
   from 353 to 86 per cycle.

**cpplint** on the 10 files reports 126 errors before and 126 after the patch, so the patch
adds none. No line is over 90 columns.

## Gate 1: bitwise (payload md5 after `<par_end>` of every bin/rst/tab, plus hst `cmp`)

GPU: jobs 11944415 (v1) and 11944727 (new2); CPU: `gate/cpu_gate.sh` (v1) and
`gate/cpu_gate2.sh` (new2). Each run is compared with base.

- **prod4 and hyd4 dhj** (`deep_hot_jupiter_cs_{prod4,hyd4}.athinput`, from scratch):
  identical on 1 and 2 ranks.
  - GPU: 300 cycles, both v1 and new2.
  - CPU: v1 at 300 cycles on 2 ranks; new2 at 20 cycles on 1 and 2 ranks (see
    `gate/cpu2/`). The 300-cycle 1-rank CPU runs were stopped at cycle 100-200 on request.
- **prod4 from the rot-283 restart** (300 cycles, 2 GPUs): bin, rst and hst are identical in
  all 16 timing arms (base, v1, new2, each with and without the env var).
- **Cartesian hydro linear wave** (uniform grid and SMR): identical on 1 and 2 ranks, CPU and
  GPU.
- **MHD SMR linear wave** (the FC prolongation and flux-correction paths): identical on 1 and
  2 ranks, CPU and GPU. The CPU runs were repeated 2x.
- **M1 slab** (`he_slab_m1_2d_V3edd`, 200 cycles; the 2-rank run uses meshblock nx2 = 16):
  identical on 1 and 2 ranks, CPU and GPU.
- **M1 3-D box timing arms:** hst identical between base, v1 and new2 on 1 GPU and on
  2 GPUs.
- **Not usable as a gate: `linear_wave_mhd_amr`.** HEAD itself is not reproducible on it:
  - The 2-rank run segfaults at cycle 4 in both base and new.
  - The 1-rank run gives run-to-run different rst files, even at cycle 0, in base as well as
    in new (`gate/cpu_rep/`).
  - This comes from HEAD, not from the patch. It needs its own look.

## Gate 2: GPU timing (apudev; same binary per arm; interleaved, repeated)

**M1 3-D box.** Setup:
- Input: 84x104x104 (4 MeshBlocks), `implicit_line_solver = pcr`, `implicit_bcg_sync = 1`
  (`gate/inp/m1_3d.athinput`).
- 200 cycles; ms/cycle is measured from cycle 50 to cycle 200.
- Jobs 11944655 (v2) and 11944401 (v1), node vipa1001.

| GPUs | arm | ms/cycle (rep a, b) | zone-cyc/s | speed-up |
|---|---|---|---|---|
| 1 | base | 246.1, 248.2 | 3.62e6 | - |
| 1 | v1 | 233.9, 235.5 | 3.79e6 | 1.05x |
| 1 | **new2** | **191.8, 193.1** | **4.62e6** | **1.28x** |
| 1 | new2 + `HSA_NO_SCRATCH_RECLAIM=1` | 187.6, 187.3 | 4.78e6 | 1.32x |
| 2 | base | 213.5, 215.0 | 4.16e6 | - |
| 2 | v1 | 205.1, 203.6 | 4.36e6 | 1.04x |
| 2 | **new2** | **160.7, 160.7** | **5.55e6** | **1.33x** |
| 2 | new2 + `HSA_NO_SCRATCH_RECLAIM=1` | 154.4, 155.2 | 5.76e6 | 1.38x |

In every arm: Picard mean 4.385 (1 GPU) and 4.375 (2 GPUs), NON-CONVERGED 0, and the same
inner-iteration totals.

**Traces** (rocprofv3, nlim = 20, cycles 2-20, rank 0; `tools/trace.py` on
`gate/prof/*`, `gate/scr/pnew_nsr`):

| arm | ms/cycle in the trace | idle ms/cycle (%) | launches/cycle, global / constant / arg | idle before pack, unpack (ms/cycle) | hipStreamSync / hipEventSync / hipDeviceSync per cycle |
|---|---|---|---|---|---|
| 1 GPU base | 297.2 | 95.5 (32.1 %) | 553 / 353 / 2049 | 49.1, 15.5 | 1003 / 353 / 321 |
| 1 GPU v1 | 276.9 | 72.5 (26.2 %) | 0 / 86 / 2869 | 39.3, 5.6 | 451 / 86 / 45 |
| 1 GPU new2 | 225.7 | 26.7 (11.8 %) | 0 / 86 / 2869 | < 0.6, < 0.6 | 451 / 86 / 45 |
| 1 GPU v1 + no-reclaim | 222.4 | 21.8 (9.8 %) | 0 / 86 / 2869 | < 0.3, < 0.3 | 451 / 86 / 45 |
| 2 GPUs base | 254.8 | 108.6 (42.6 %) | 553 / 354 / 2049 | 47.5, 27.4 | 1004 / 354 / 321 |
| 2 GPUs new2 | 204.3 | 60.8 (29.7 %) | 0 / 87 / 2869 | < 0.6, 18.2 | 451 / 87 / 321 |

- **What is left on 1 GPU** is the BiCGStab reductions: about 15 ms/cycle, 3 x 128 gaps of
  about 40 us.
- **On 2 GPUs** the `RecvAndUnpackCC` gap (18 ms/cycle) is the MPI wait plus the
  still-needed fence. The 321 `hipDeviceSynchronize` per cycle are those fences.

**prod4 (cubed-sphere dhj MHD) from `bench/cs_mhd_prod3/rst/dhj.00567.rst`, 300 cycles,
2 GPUs** (jobs 11944519 (v1) and 11944726 (v2), `cpu time used`, s):

| arm | runs | mean |
|---|---|---|
| base | 22.29, 22.14, 21.85, 22.25; 22.51, 22.61 | 22.28 |
| v1 | 21.83, 22.06, 21.67, 21.82 | 21.85 |
| new2 | 23.14, 22.37 | 22.76 |
| base + `HSA_NO_SCRATCH_RECLAIM=1` | 16.94, 17.09 | 17.02 (**1.31x**) |
| new2 + `HSA_NO_SCRATCH_RECLAIM=1` | 17.12, 16.97 | 17.05 |

**The patch alone is neutral on prod4, within about 2 % noise.**
- Its 20 global-memory launches per cycle go to 0, and `hipStreamSynchronize` falls from 99
  to 79 per cycle.
- But prod4's idle is MPI waits and the 150 `hipDeviceSynchronize` per cycle, and its cs
  pack keeps its scratch.

**The environment variable is worth 24 % of prod4 cpu time by itself**, with bitwise
identical output (same md5 in all arms).

## Recommendations

- **Apply halo.patch.** It gives M1 1.28x on 1 GPU and 1.33x on 2 GPUs, bitwise.
- **Export `HSA_NO_SCRATCH_RECLAIM=1` in every GPU submit script.** This covers the prod4 and
  hyd4 run dirs at a link boundary. The job scripts in `bench/cs_mhd_prod4` and
  `bench/cs_hyd4_prod` were NOT touched here.
- **Next targets:**
  - The BiCGStab reductions: about 15 ms/cycle on 1 GPU.
  - The MPI-path fences: 321 `hipDeviceSynchronize` per cycle on 2 GPUs. `Kokkos::fence()`
    could become a stream fence.
  - The 86 constant-memory launches per cycle, from `par_for` wrappers over large functors.
