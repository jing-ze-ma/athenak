# vet_sc SC sweep on several GPUs, step 1: messages, hybrid rays x space, moments

Date 2026-09-23, viper. Branch `sc-scale` (worktree `/viper/ptmp2/jinma/wt_scscale`), base
5cd0a86d. Only `src/rad_m1/rad_m1_vet.cpp` changes. Follows
`tests_m1/runs_3m_vetmb/README_SCALING.md` sect. 5, items 1-3.

- **Binaries** (`/viper/ptmp2/jinma/scscale_0923/`):
  - GPU `gpu/bin/athena_new_gpu` md5 13bd729d (this branch), `gpu/bin/athena_ref_gpu` md5
    ebe72ffe (HEAD 5cd0a86d);
  - CPU `bin/athena_new_cpu` md5 0f9009c5, `bin/athena_ref_cpu` md5 a7ea154a.
  - Build dirs deleted. The Kokkos source came from the GPFS snapshot
    `/viper/u2/.snapshots/daily-20260923063002/jinma/ATHENAK/bench/wt_rgbox/kokkos`: see
    the WARNING at the end.
- **Scripts** (copies in `scripts/`): `time.sh` (interleaved GPU timing, repeat b in
  reverse order), `prof.sh` + `vetprof.py` (rocprofv3 split), `tsum.py` (tables), `gate*.sh`,
  `compare.sh`, `hstdiff.py`, `cmp_planes.py` (CPU gates), `mkinp.py` (inputs).
- **Runs:** `/viper/ptmp2/jinma/scscale_0923/{cpu,gpu/runs}`.
- **Every GPU job** exports `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.

## 1. What changed

All new keys are read only when given (an input without them writes the restart it always
wrote).

1. **One message per neighbour rank and exchange** (`vet_mb_agg`, default **true**, exact).
   - Before: one MPI message per (block, horizontal slot), i.e. up to 8 x nmb per exchange.
   - Now: the pieces of all (block, slot) pairs that go to one rank are packed into one
     flat buffer, ordered by the receiver's (block, slot); plan P for the intensity planes,
     plan B for the extinction/source band. Pack and unpack are one kernel each, as
     before.
   - `vet_mb_agg = false` keeps the old path (used as the A/B arm).
2. **Hybrid decomposition** (`vet_mb_agroup = G`, default 1 = off, round-off).
   - The G consecutive ranks of a group (G = GPUs per node) all hold the extinction,
     source and bottom E, F of the group's blocks (one gather per call). Each sweeps the
     rays r = p, p + G, ... (p = rank in the group) over all the group's blocks.
   - The per-layer band exchange stays on the rank for group neighbours (device copy) and
     goes to the rank with the same p in the neighbour group otherwise (item 1 applies).
   - The partial moments go to the block owners, which sum them in group-rank order.
   - It generalises the `vet_mb_angles` prototype: G = number of ranks is the pure angle
     split and reproduces the prototype bit for bit (gate C). It needs no whole-mesh copy
     per GPU: memory is O(group blocks), 200 MB/rank of sweep buffers on 2 GPUs (94 MB
     spatial).
   - Excludes `vet_mb_angles`, `vet_mb_lag`, `vet_mb_kernel = cell`, `vet_milne`.
3. **Moments** (exact, always on in the ray-kernel path).
   - **The moment launch was slow for a layout reason, not a latency one.** It added the
     ten sums into `vet_cell(m, n, k, j, i)`, i fastest, from threads with j fastest: every
     addition was an uncoalesced read-modify-write.
   - Now the launches add into `mom(block, n, i, k, j)` (j fastest, coalesced), and one
     tiled transpose per call (`VetMomOut`, team scratch) writes `vet_cell`. The numbers
     are the same.
   - The moment launch went from 71 to 34 us per launch (1 GPU, `runs/prof1` vs `prof2`).
   - **`vet_mb_mom_fuse = true`** (exact, default off) is the fusion the TO-DO asked for:
     one team per (block, row, chunk of <= 512/nray cells), one thread per (ray, cell),
     intensities in team scratch, then the ordered per-hemisphere sums with `VetMomAdd`.
     It is **slower** on MI300A (7.5 vs 3.9 ms/call on 1 GPU). The team layout loses the
     coalescing of the flat ray kernel (8 rays x 8 cells per wavefront). Launch bounds 512
     changed nothing. Kept off, for the record.
   - A moment kernel with one thread per moment (10x the threads) was tried and removed:
     5.6 vs 4.8 ms/call.
4. The extinction/source copy no longer re-zeroes the moments (`VetShortChar` already
   does). The `vet_mb_lag` diagnostic still zeroes them. 212 -> 42 us.
5. **Bug fixed during development (hybrid only):** the x1 edge copy of the full-tensor D
   ghosts in `VetFullTensor` indexed the (sweep-block) `lx1` table with the own-block
   index. It showed up as a 1e-7 K/J drift in the x1 x x2 split with `vet_mb_agroup`.
   The own blocks now start at `sown`. Default runs are not affected (sown = 0).

## 2. Gates (CPU, `cpu/gate_final2.out`, final binary 0f9009c5 vs HEAD a7ea154a)

- **A. Default arms, byte identity of every output file** (bin, rst, hst; 15 or 11 files):
  - seeded 2-D slab, tlim 3: U1, F1 (1 block); x2 4 blocks on 4 and 1 ranks, uniaxial
    and full; x2 2 blocks on 2 ranks; x1 4 blocks on 4 ranks; x1 x x2 full on 4 ranks;
  - **odd global nx1 = 81**, 3 x1 blocks x 2, 3 ranks (the lu == ld layer);
  - 3-D box 84x32x32, 4 blocks, on 1, 2 and 4 ranks; 3-D with nx1 = 83 on 4 ranks.
  - **13/13 identical.** `vet_mb_agg` is on by default, so every multi-rank arm exercises
    the new message path. The same holds on the inputs that spell out every key (A2).
- **B. Exact options against HEAD:** agg off (x2r4, T3r4); fuse on x2b2r2, x2r4, x1r4,
  x12r4F, odd r3, T3r4, T3r1, T3odd r4; fuse + halo 3; fuse + halo 2 + overlap; halo 3;
  agg off + halo 3 + overlap on 2 ranks.
  - **Both hst files are identical in all 14 arms.**
  - bin/rst differ only in 1-11 header bytes (offset < 14 kB), which are the echoed keys.
- **C. Hybrid (round-off):**
  - Planes against HEAD, 2-D x2r4 and x1 x x2 full, G = 2 and 4: at call 0, K/J
    <= 6.7e-16 and J <= 1.1e-15 relative. Calls 5-18: K/J <= 2.6e-11, J <= 5.7e-10. This
    is the drift the prototype and the decomposed runs show.
  - G = 4 on 4 ranks equals the `vet_mb_angles` prototype bit for bit (e_hy4 = e_ang, all
    19 calls).
  - 3-D T3r4 G = 2 / 4 and T3r2 G = 2: user.hst <= 1.6e-8 relative after 6 steps.
  - fuse + halo 3 and agg off on top of G = 2 are hst-identical to G = 2.
- **GPU:**
  - The exact options are hst-identical on the GPU as well: fuse, halo 1/3 and agg off vs
    on, at 1, 2, 4 and 8 GPUs, 100 cycles (`runs/t3`, `t4`, `g4`, `g8`).
  - The transposed moment layout is hst-identical to the old one (`runs/t4/g1_base` vs
    `t3/g1_base`).

## 3. GPU cost

**Setup:**
- 3-D He box, 84 x 104 x 104, **8 MeshBlocks of 84 x 52 x 26**, the same decomposition for
  every GPU count.
- vet_sc full tensor, 64 rays.
- FAST switches spelled out in the input: `implicit_halo_direct`, `implicit_od_cache`,
  `implicit_krylov_fuse = 3`, `implicit_op_stencil`, `implicit_precond = rbgs_fwd`.
- `vet_mb_halo = 3`. Input `gpu/he3d_fast.athinput` (`scripts/mkinp.py`).
- 100 cycles; ms/cycle is measured over cycles 10-100.
- Repeats a and b, interleaved, the same two binaries in each job.
- "hydro" = the same box without `<rad_m1>` (`gpu/he3d_hydro.athinput`), same job.

- **Jobs:**
  - 11947695 (`runs/f12`, 1 and 2 GPUs, apudev);
  - 11947696 (`runs/fw2`, weak 2 GPUs, apudev);
  - 11947723 (`runs/f4`, 2 nodes, apu);
  - 11947724 (`runs/f8`, 4 nodes, apu).
- **Binaries:** new 13bd729d, HEAD ebe72ffe.
- **Arms:**
  - "HEAD" = the unmodified binary;
  - "new" = this branch with its defaults (agg on);
  - "no agg" = new with `vet_mb_agg = false`;
  - "hybrid" = new with `vet_mb_agroup = 2`.
- SC is the fenced "SC seconds per call". Every entry gives repeat a, then repeat b.

**Strong scaling** (84 x 104 x 104 = 908 k cells):

| GPUs (nodes) | SC ms/call HEAD | SC new | SC no agg | SC hybrid | ms/cycle HEAD | ms/cycle new | ms/cycle hybrid | **hydro ms/cycle** |
|---|---|---|---|---|---|---|---|---|
| 1 (1) | 4.83, 4.83 | **3.92, 3.92** | - | - | 57.7, 55.3 | 56.1, 55.1 | - | 27.2, 22.7 |
| 2 (1) | 5.91, 6.01 | 4.71, 4.70 | 5.65, 5.69 | **4.15, 4.15** | 54.7, 54.6 | 54.5, 53.4 | 54.4, 53.2 | 16.8, 16.6 |
| 4 (2) | 5.09, 5.24 | **4.10, 4.35** | 4.87, 5.04 | 4.45, 4.71 | 44.3, 44.7 | 43.4, 47.3 | 44.6, 48.3 | 14.2, 14.1 |
| 8 (4) | 4.91, 4.73 | **4.23, 4.16** | 4.50, 4.35 | 4.29, 4.29 | 39.1, 38.7 | 40.4, 39.2 | 39.9, 39.2 | 12.6, 12.5 |

**Weak scaling** (8 blocks = 908 k cells per GPU; the box doubles with the GPU count:
84 x 208 x 104 on 2, 84 x 208 x 208 on 4, 84 x 416 x 208 on 8):

| GPUs | SC ms/call HEAD | SC new | SC hybrid | ms/cycle HEAD | ms/cycle new | ms/cycle hybrid | **hydro ms/cycle** |
|---|---|---|---|---|---|---|---|
| 1 | 4.83, 4.83 | 3.92, 3.92 | - | 57.7, 55.3 | 56.1, 55.1 | - | 27.2, 22.7 |
| 2 | 9.98, 9.91 | **6.68, 6.62** | 6.74, 6.72 | 83.2, 83.4 | 83.5, 81.6 | 80.9, 86.5 | 24.2, 24.9 |
| 4 | 9.70, 9.87 | **6.88, 7.01** | 8.92, 8.94 | 84.0, 84.3 | 83.4, 83.5 | 86.0, 82.8 | 23.4, 23.5 |
| 8 | 10.61, 10.57 | **7.52, 7.54** | 9.77, 9.75 | 99.7, 93.5 | 90.1, 92.8 | 93.7, 91.5 | 25.7, 25.9 |

Other measurements:
- **Fused moments** (`vet_mb_mom_fuse`): 7.47, 7.47 ms/call on 1 GPU (`runs/f12`); 6.67 and
  6.71 on 2 GPUs (`runs/t4`).
- **Angle prototype on 2 GPUs:** 4.18, 4.21 (`runs/g2c`, earlier binary).
- Picard 2.700 per step (2.71-2.76 on the larger boxes) and NON-CONVERGED 0 in every
  arm.
- **GPU exactness:** the "new" and "no agg" user/hydro hst are **identical to HEAD** at 1,
  2, 4 and 8 GPUs and on every weak box.
- **Hybrid differences:** the hybrid differs by 4.2e-9 (user.hst) from HEAD on 2 GPUs.
  HEAD's own 1-GPU vs 2-GPU difference is 4.8e-9.
- **Sweep memory** on 2 GPUs: 94 MB/rank spatial, 200 MB/rank hybrid.

**What this says:**
- **SC per call, new vs HEAD:** -19 % on 1 GPU, -21 % on 2 GPUs, -14 to -20 % on 4-8 GPUs.
  - Agg alone gives 5.65 -> 4.71 (2 GPUs) and 4.87 -> 4.10 (4 GPUs).
  - The moment layout gives the 1-GPU gain.
  - The weak boxes gain 29-33 %: 10.6 -> 7.5 ms at 8 GPUs, from many more messages per
    exchange.
- **The hybrid is the best 2-GPU arm** (4.15 ms, as fast as the angle prototype,
  without its whole-mesh copy). Across nodes it loses to the plain spatial split with
  aggregated messages: 4 GPUs 4.45-4.71 vs 4.10-4.35; weak 8 GPUs 9.8 vs 7.5. Recommend
  `vet_mb_agroup = 2` only for single-node runs.
- **SC still does not strong-scale:** 3.9 ms on 1 GPU, 4.2 ms on 8. The weak-scaling SC
  cost doubles from 1 GPU to 2+ GPUs (3.9 -> 6.6-7.5 ms at a fixed 0.9 M cells/GPU), which
  is the per-layer exchange latency (sect. 4).
- **SC share of the cycle:** 7-8 % of the cycle and 12-13 % of the radiation part.
- **Against hydro:**
  - Radiation (cycle - hydro) is 29-33 ms on 1 GPU, 38 on 2, 30-33 on 4, 27 on 8.
  - Hydro scales 23-27 -> 12.5 ms.
  - At 0.9 M cells/GPU (weak), the whole M1 cycle is 3.3-3.6x hydro (e.g. 8 GPUs: 90-93
    vs 25.8 ms).
  - The SC sweep alone (7.5 ms) is now 0.3x hydro. The Krylov side is the remaining
    60 ms: see sect. 6.

## 4. Where the SC time goes now (rocprofv3, `runs/prof2`, per call)

Final binary, 8 blocks, halo 3, 10 calls (`python3 vetprof.py runs/prof2/<arm>/prof/rank_0_kernel_trace.csv`):

| per call | 1 GPU | 2 GPUs, spatial (agg) | 2 GPUs, hybrid G = 2 |
|---|---|---|---|
| span of the vet kernels | 3.85 ms | 4.64 ms | 3.93 ms |
| kernel busy | 3.70 ms | 3.02 ms | 2.38 ms |
| ray launch, 84 per call | 2.17 ms (25.8 us each) | 1.37 ms (16.3 us) | 1.38 ms (16.5 us) |
| moment launch, 21 per call | 0.71 ms (33.6 us) | 0.64 ms (30.5 us) | 0.47 ms (22.3 us) |
| local band copy, 28 per call | 0.42 ms (14.9 us) | 0.23 ms (8.1 us) | 0.25 ms (9.0 us) |
| pack + unpack, 28 each | - | 0.50 ms | - |
| VetMomOut (transpose, once) | 0.10 ms | 0.06 ms | 0.06 ms |
| gaps (MPI, host) | 0.15 ms | 1.62 ms | 1.55 ms |

- **Before** (`runs/prof1`, the HEAD message path and moment layout, 1 GPU):
  - ray launches 2.16 ms;
  - moment launches 1.49 ms (71 us each);
  - csw copy 0.21 ms;
  - span 4.70 ms.
- **The hybrid gaps** are the group gather and the moment exchange: two host-blocking
  point-to-point rounds per call.

## 5. Ranked TO-DO: SC sweep

**Diagnosis.** The sweep is now a chain of short launches. Per call there are:
- 84 ray launches;
- 21 moment launches;
- 28 local band copies;
- on several ranks, 28 x (pack, fence, MPI, unpack).

**Evidence** (`runs/prof2`):
- On 1 GPU the chain is busy 3.70 of 3.85 ms.
- On 2 GPUs, half the work per launch only cuts the ray launch from 25.8 to 16.3 us.
  The MPI waits add 1.6 ms (span 4.64, busy 3.02 ms).
- More GPUs therefore cannot bring the call below about 84 x (ray launch + exchange
  latency).

**Ranked items:**

1. **Several layers per launch (temporal blocking on the device).**
   - A team takes one ray (or a few rays with similar reach) and a tile of the plane, and
     sweeps K layers in scratch with a team barrier per layer.
   - The upwind halo it needs is its own reach per layer, `rdep`: 1 cell for most rays, 4
     for the most oblique, times K. That is the per-ray band of `vet_mb_halo = K`, which is
     already exchanged every K layers.
   - It cuts the ray launches from 84 to 84/K and removes the per-layer global round trip
     of the intensity plane.
   - The moments can stay in the separate launch: one per K layers, from a ring of K
     planes.
   - Expected: 84 -> 28 launches at K = 3. This is the only item that changes the
     scaling slope.
2. **Hide the exchange.**
   - On 2 GPUs, 1.6 ms of the 4.6 ms call is host-side MPI wait: fence, Waitall, and 2 x
     28 pack/unpack launches.
   - Retest `vet_mb_overlap` with the aggregated messages. Before aggregation it gave
     nothing.
   - Merge the local-copy kernel into the pack kernel: one launch fewer per exchange,
     0.42 ms/call of 14.9 us launches on 1 GPU.
3. **Choose the halo per decomposition.**
   - On 1 GPU, halo 3 costs +0.6 ms against halo 1 (4.78 vs 4.19 ms, `runs/g1c`): the
     redundant overlap on 26-cell-wide blocks.
   - Halo 3 pays only when the exchange goes through MPI.
   - Default suggestion: `vet_mb_halo = 1` on one rank, 3 on several.
4. **Hybrid (`vet_mb_agroup`).**
   - Worth it within a node: 2 GPUs, 4.15 vs 4.71 ms/call, sect. 3.
   - Across nodes it is no better than the spatial split at this size.
   - Overlap the group gather with the source kernel, and send the moment partials as
     soon as the last launch that touches a block is done.
5. **Single-block path** (`vet_mbs == nullptr`, one MeshBlock on one rank).
   - Still the old per-column kernel, as in README_SCALING item 5.
   - Route it through the banded ray kernel (`vet_mb_force` already does; the cost is the
     one-block band copies).

## 6. Ranked TO-DO: the Krylov side, and multi-group

The M1 part of the cycle (ms/cycle - hydro, sect. 3) is almost flat from 1 to 8 GPUs. So is
the SC share. README_FAST sect. 4 has the cause: 588 launches per cycle, about 41 Krylov
iterations per step, 2 blocking reductions per iteration. Every launch and every
reduction costs the same on 1 or 8 GPUs, and the multi-rank run adds MPI latency on top.
Ranked by expected gain for multi-GPU:

1. **Latency hiding of the two reductions per iteration.**
   - On several ranks, every reduction is a device -> host copy plus an `MPI_Allreduce`,
     about 2 x (20-50) us per iteration, i.e. 2-4 ms/cycle at 41 iterations.
   - A pipelined BiCGStab (Ghysels-Vanroose p-BiCGStab: one fused reduction per
     iteration) with `MPI_Iallreduce` overlapped with the operator application removes most
     of it.
   - Keep the convergence test lagged by one iteration, as `implicit_krylov_fuse = 3`
     already does.
2. **Direct halo exchange on several ranks.**
   - `implicit_halo_direct` (-34 ms/cycle on 1 GPU) only works when every neighbour is on
     the rank. On 2+ GPUs the implicit halos fall back to the general boundary machinery:
     per-variable buffers, one message per block and neighbour.
   - Port the pattern used here instead: a precomputed neighbour table, one pack kernel,
     one message per neighbour rank, one unpack kernel. That is about 2 + 2 launches and
     one message pair per neighbour per halo, instead of the full machinery, about 100
     times per cycle.
3. **Device-resident scalars and a convergence check every N iterations.**
   - On 1 rank this removes the host sync per iteration (README_FAST TO-DO 6). On several
     ranks it makes item 1 cheaper.
   - It is also the precondition for HIP graphs: replay one Krylov iteration (operator,
     preconditioner half-sweeps, vector updates) as one graph, i.e. one launch per
     iteration instead of ~14.
4. **Load at least 1 M cells per GPU in production.**
   - At 0.9 M cells/GPU the weak rows of sect. 3 stay at 81-93 ms/cycle from 2 to 8 GPUs.
   - Strong scaling below ~0.5 M cells/GPU buys almost nothing for M1 while the launch
     count is fixed.
5. **Multi-group: make the cost sublinear by batching the groups in every launch.**
   - Put the group index g in the pack dimension, next to m. Every kernel (operator,
     preconditioner line solves, vector updates, SC ray launch) then covers all groups in
     one launch.
   - One halo message per neighbour rank carries all groups (item 2, and item 1 of this
     README for the SC planes).
   - One reduction carries the per-group dot products: the multi-value reducer, with
     NREDUCTION_VARIABLES raised, or a 2-stage team reduction.
   - Then the launch count and the reduction count do not grow with the number of groups,
     only the work per launch does. Today a launch carries 450-900 k cells per GPU and runs
     at 10-90 us, most of it latency: the SC ray launch takes 25.8 us for 8 blocks and
     16.3 us for 4 (`runs/prof2`), so half the work saves only 37 %.
   - So the first several groups ride in the idle part of each launch. The estimate from
     these numbers is about 1.5-2x the grey cost for 4 groups, until the launches become
     bandwidth-bound (roughly 4-8 M cell-group updates per launch on MI300A).
   - The group coupling (the absorption/emission exchange between groups) belongs in the
     cell-local Picard/Newton level, not in the Krylov operator. The transport operator
     is then block-diagonal in g and the Krylov vectors simply get longer.
   - The SC sweep batches the same way: rays x groups per thread index, one launch per
     layer for all groups. Per-group opacity only changes chi and S, which are per-cell
     loads.

## WARNING: the kokkos symlink is broken for every checkout

- **The broken link:** the repository tracks `kokkos` as a SYMLINK to
  `/viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos`. `bench/wt_rgbox` was deleted on
  2026-09-23 between 07:00 and 07:15, so the link now dangles in:
  - the main checkout;
  - `bench/wt_he4`;
  - every new worktree.
- **Consequence:** no build can configure or compile.
- **How this work built:** from snapshots whose `kokkos` points at the read-only GPFS
  snapshot `/viper/u2/.snapshots/daily-20260923063002/jinma/ATHENAK/bench/wt_rgbox/kokkos`
  (Kokkos 4.6.2, d8e9af031, the same tree the earlier builds used).
- **Not fixed here:** the user should decide how to repair it, e.g. restore that directory
  from the snapshot, or make `kokkos` a real submodule checkout.
