# vet_sc sweep on several GPUs: profile, fixes, measurements, TO-DO

Date 2026-09-23, viper. NOT COMMITTED.

- **Patch:** `bench/m1_scscale_0923/scscale.patch`. It touches only `src/rad_m1/rad_m1_vet.cpp`.
  - Made against d14f96de.
  - `git apply --check` passes on d14f96de, a6202226 and 66d5c59a. `rad_m1_vet.cpp` is the
    same file in all three.
- **Where things are:**
  - Snapshots, CPU gate runs and scripts moved to `/viper/ptmp2/jinma/m1_scscale_0923/`
    (`ref`, `new`, `ref2`, `new2`, `new3`, `cpu`, `runs`). Build dirs deleted.
  - Binaries: `/viper/ptmp2/jinma/m1_scscale_0923/bin_*`.
    - `bin_ref2_gpu` = a6202226 (md5 e4b6e039).
    - `bin_new3_gpu` = a6202226 + patch (md5 7eb1d89f).
    - `bin_ref_cpu` = d14f96de (md5 5ef48abc).
    - `bin_final_cpu` = d14f96de + patch (md5 e8fe6ae2).
  - Scripts (in `bench/m1_scscale_0923`): `time.sh`, `prof.sh`, `stab.sh`,
    `tools/{tsum,vetprof}.py`. `runs/...` paths below are under ptmp2.
  - Inputs: `he3d_vet_a62.athinput`. This is the 3-D box of RESULTS.txt gate 3/4, plus the
    a6202226 recommended settings `lres_test=false conv_est=true predictor=step
    lin_ew_max=1e-2 bcg_sync=1`.

## 1. Profile of HEAD (rocprofv3 kernel trace)

Setup: 3-D box 84x104x104, 4 MeshBlocks of 84x52x52, vet_tensor = full, 64 rays. Ten calls
per run.

- Per-call split from `runs/prof_ref_{1,2}`, analysed with `tools/vetprof.py`:

| per call | 1 GPU | 2 GPUs |
|---|---|---|
| span of the vet kernels | 18.8 ms | 24.3 ms |
| sweep kernel, 84 launches | 14.45 ms (172 us each) | 14.2 ms (169 us each) |
| column source kernel | 1.21 ms | 1.02 ms |
| band copy/pack/unpack | 0.55 ms | 1.14 ms |
| gaps (host/MPI) | 2.2 ms | 7.7 ms |

- **Diagnosis:** the cost is not the exchange. It is the **sweep kernel itself**.
  - It runs one thread per (k, j) column, i.e. 10816 threads, and each thread loops over
    64 rays. That is latency bound at under 1 wave per CU.
  - On 2 GPUs each GPU gets half the threads and **the same time per launch**.
  - It also went through Kokkos' constant-memory launch path, which puts a host wait
    between launches.
- The band-exchange functors never reached the > 32 kB global path, so d14f96de was not
  their problem. Their DualView captures did put them on the constant-memory path.

## 2. What the patch does (all in rad_m1_vet.cpp)

Each item is exact (same arithmetic, same numbers) unless it says otherwise.

1. **Ray-parallel sweep** (`VetRayLaunch`): one thread per (block, ray, cell) per layer.
2. **Moments in a separate launch** (`VetMomLaunch`).
   - One launch sums the moments of B = 4 layers (`vet_mb_mom_batch`, default 4). It keeps
     a ring of max(B, 2) intensity planes.
   - There is one thread per (layer, hemisphere, cell). The rays are summed in the old
     order, and the intensity loads are batched 8 at a time.
   - A cell gets exactly two sums, one up and one down, onto 0. That addition is exact in
     either order.
3. **Per-ray band exchange** (`VetPlanePost/Wait`).
   - Each ray sends only the ghost cells it reads: one side, floor(reach)+1 cells. For
     most rays that is 1 cell instead of 4.
4. **Kernel-argument launches for every sweep kernel** (`VetFor`, HintLightWeight, with
   `static_assert` < 3 kB). There are no host waits between launches.
5. **Cell-parallel source kernel** (when Milne is off).
6. **`vet_mb_halo = K`**: a band K reaches wide is exchanged every K layers, and the
   overlap is recomputed redundantly per ray.
7. **`vet_mb_overlap`**: the interior of the next layer is computed while the MPI is in
   flight.
8. **`vet_mb_kernel = cell`**: keeps the old kernel.
9. **PROTOTYPE `vet_mb_angles`**, angle decomposition (`VetSweepAng`).
   - Every rank gathers chi and S of the whole mesh (point-to-point).
   - Each rank sweeps rays r = p, p+P, ... over the global periodic mesh with no per-layer
     exchange.
   - The partial moments go to their owners and are summed in rank order.
   - MPI collectives on device buffers took 115 ms per call (run t2 `A2`). They were
     replaced by Isend/Irecv.

All new keys are read only when they are given, so inputs without them keep the restart
echo.

## 3. Gates

**Gate 1, CPU byte-identity.** Runs are in ptmp2 `cpu/`. The summary is in
`cpu/gate_final.out`, final code.

- Default arms against d14f96de, every file `cmp`-identical:
  - seeded 2-D slab, tlim 3, bin/rst every 1 s: uniaxial and full; 1 block; x2 4 blocks on
    1 and 4 ranks; x2 2 blocks on 2 ranks; x1 4 blocks on 4 ranks; x1 2 x x2 2 on 4 ranks
    (15/15 files each);
  - 3-D box 84x32x32, 4 blocks of 16x16, on 1, 2 and 4 ranks (11/11 files each);
  - closure m1 (default), 15/15, run `q2_def`.
- Option arms against HEAD on the same input:
  - The arms: halo 2/3/4, overlap, batch 1/3/7/8, kernel=cell, 2-D and 3-D.
  - The hst files are identical.
  - bin/rst differ only in 1-6 header bytes, which are the echoed option values.
- `vet_mb_angles`:
  - On 1 rank with 4 blocks it is identical, apart from the 5-byte echo.
  - On 2 and 4 ranks (x2, x1 and x1 x x2 splits):
    - At call 0, K/J differs from HEAD by at most 6.7e-16 and J by 1.1e-15 relative
      (`cmp_planes.py`).
    - Later calls drift to 2e-11. That is the same drift as the decomposed runs of
      RESULTS.txt.

**Gate 3, GPU stability.** 3-D box, 6 min each (`runs/stab`, jobs 11945284 and 11945285,
binary `bin_new2_gpu`).

- `bin_new2_gpu` (since deleted) is the same algorithm with the earlier moment kernel.
- S1 (1 GPU), S2 (2 GPUs), S2H (halo 3 + overlap) and S2A (angles) ran 3180-3700 steps.
- In every run: NON-CONVERGED 0, Picard 2.90/4, min E 1.609e5, F1top/Fin 0.9998-1.0000.
- The four runs are identical in these statistics.

## 4. GPU cost

Setup:
- Job 11945324 on node vipa1001, apudev.
- Binaries `bin_ref2_gpu` / `bin_new3_gpu`, interleaved, 2 repeats.
- nlim 100, recommended settings. `HSA_NO_SCRATCH_RECLAIM=1`.
- Runs are in `runs/t4`; the table was made with `tools/tsum.py`.
- SC = the fenced "SC seconds per call".

| arm | SC ms/call | ms/step (c10-100) |
|---|---|---|
| HEAD, 1 GPU | 18.20, 18.28 | 116.4, 118.5 |
| patch, 1 GPU | **3.96, 3.98** | **103.6, 105.1** |
| patch B=1, 1 GPU | 6.35, 6.35 | 108.3, 107.4 |
| angles, 1 GPU | 3.65, 3.78 | 104.0, 103.6 |
| HEAD, 2 GPUs | 23.47, 23.26 | 107.0, 107.9 |
| patch, 2 GPUs | 8.19, 8.20 | 92.6, 92.3 |
| + halo 2 | 6.16, 6.32 | 90.2, 90.8 |
| + halo 3 | 5.64, 5.61 | 90.0, 89.9 |
| + halo 3 + overlap | 5.71, 5.71 | 89.7, 90.2 |
| **angles, 2 GPUs** | **4.19, 4.19** | **88.8, 89.2** |

- **SC on 1 GPU:** 18.2 → 4.0 ms, i.e. 4.6x.
- **SC on 2 GPUs:**
  - 23.4 → 8.2 ms with the default options.
  - 5.6 ms with halo 3.
  - 4.2 ms with angles.
- **Step time:** −13 ms/step on 1 GPU and −15 to −19 ms/step on 2 GPUs.
- **The 1-GPU cost per cell is not recovered on 2 GPUs.** The best 2-GPU run (angles) takes
  4.2 ms, where 1 GPU takes 4.0 ms. Traces of the patched binary (`runs/prof_n3_*`,
  `prof_n3a_2`), per call:
  - **1 GPU:** 3.8 ms of kernels, with no gaps.
    - Ray kernel: 84 x 17 us.
    - Moments: 21 x 70 us.
  - **2 GPUs, spatial split:** 3.5 ms of kernels in a 7.7 ms span. The remaining 4.2 ms is
    the per-layer MPI.
  - **2 GPUs, angles:** 2.5 ms of kernels in a 4.0 ms span. The rest is the gather and the
    moment exchange.
  - Per-layer kernels barely shrink with half the work (17 → 12 us), because they are
    launch/latency bound.
- **4 GPUs were not measured.** apudev nodes have 2 GPUs and MaxNodes = 1.
- **Memory per GPU at production size:**
  - The production size is this box: 84x104x104, 64 rays.
  - Banded path: 126 MB on 1 rank with 4 blocks, 63 MB per rank on 2 ranks, from the
    start-up line.
  - Angles prototype: 194-212 MB per GPU on top of that, and O(N_global) per GPU:
    - gathered chi, S: 2N;
    - partial moments: 10N;
    - received partials: 10N;
    - planes: max(B, 2) x nray/P x plane.
  - For a 10x larger mesh the angle buffers are about 2 GB per GPU. The moment exchange is
    about 10N doubles per call per rank, i.e. 73 MB here.

## 5. Ranked TO-DO for the next session

1. **Angle decomposition, hybrid version.**
   - Angles within a node (xGMI) and the banded spatial split across nodes.
   - Also drop the unused banded buffers when angles are on (63 MB per rank).
   - Aggregate the moment messages to one per rank pair, and overlap the gather with the
     source kernel.
   - It is already the fastest 2-GPU arm.
2. **Aggregate band messages per neighbour rank.** Today there are up to 8 slots x nmb
   messages per exchange. Combined with halo 3 this cuts the 4.2 ms/call of 2-GPU MPI
   further. This is the cheap fix for the spatial path.
3. **Fuse the moment sum into the ray kernel.**
   - Use team scratch, with the ordered sum done by one lane per moment.
   - Or raise B to 8 and re-measure.
   - Moments are still 1.5 ms of the 3.8 ms on 1 GPU.
4. **Make the default halo 3 for multi-rank runs.**
   - It is exact, and 8.2 → 5.6 ms on 2 GPUs.
   - Overlap gave nothing measurable (5.64 vs 5.71). Drop it, or retest with aggregated
     messages.
5. **Single-block path** (`vet_mbs == nullptr`): it still uses the old per-column kernel.
   Route it through the ray-parallel kernel as well.
6. **4-GPU / 2-node timing on the apu partition.** It was not allowed on apudev.
