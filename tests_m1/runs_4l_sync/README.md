# runs_4l_sync: where the implicit M1 step loses time between 1 and 2 GPUs

- **Date:** 2026-09-24, viper.
- **Branch:** `m1-sync`, based on rt-integration e0074273 (it already contains the m1-launch
  merge). Worktree: `/viper/ptmp2/jinma/wt_m1sync`.
- **Run tree:** `/viper/ptmp2/jinma/sync_0924`:
  - `base/` = git archive e0074273;
  - `new/` = e0074273 plus this work (code-identical to the commit; the commit
    only rewraps one comment line);
  - `bin/`, `inp/`, `cpu/`, `rstgate/`, `runs/<tag>/<arm>_<rep>/`, `scripts/`.
  - The scripts are copied from launch_0923; `gaps.py` and `rgate.sh` are new. Copies
    of `gaps.py` and the arm files are in `scripts/` here.
- **GPU jobs:** apudev (4 GPUs: apu), `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`. Same
  binaries, arms interleaved, even repeats in reverse order.
- **Metric:** ms/cycle over cycles 20-110 of 120 (`scripts/tsum.py <tag>`).
- **Boxes:** the w2 setup of runs_4k_launch.
  - 1 GPU: 84x104x104 in 4 blocks of 84x52x52.
  - 2 GPUs weak: 84x208x104, 0.9 M cells per GPU.
  - 4 GPUs weak: 84x208x208.
- **Arms:** L = implicit_vimp_fold + fast_kernels + one_pass=4 + predictor_order=2.
  - E = Eddington, be.
  - V = vet_sc with the full tensor, be.
  - Vh = vet_sc, hesdirk2.
  - 2+ GPUs add implicit_halo_mpi = true and implicit_halo_overlap = true. These are
    the defaults since m1-accmerge, where they are valid.

## 1. Attribution (no code changes): the premise does not hold

Profile: `launch_0923/runs/pw2/w2_V_1/prof` (job 11955223, V, 2 GPUs, cycles 10-30).
Compared with `launch_0923/runs/p1/V_off_1/prof` (1 GPU).

Tools:
- `scripts/split.py <prof> <rank> 30 10`: kernel time and idle per category.
- `scripts/gaps.py <prof> <rank> 30 10`, new:
  - every GPU-idle gap, keyed by the kernel before, the kernel after and the host sync
    call that covers most of the gap;
  - for each HIP sync or memcpy call, its host wait and the part of that wait during
    which the GPU is idle.

**The host syncs are almost free. They wait while the GPU is still busy.**

| HIP call (2 GPUs, rank 0, per cycle) | calls | host wait ms | of which GPU idle ms |
|---|---|---|---|
| hipStreamSynchronize | 146.0 | 16.79 | 0.39 |
| hipEventSynchronize | 157.6 | 13.54 | 0.45 |
| hipMemcpyAsync | 93.5 | 9.07 | 0.72 |
| hipDeviceSynchronize | 82.7 | 7.64 | 0.35 |
| **all** | 480 | **47.0** | **1.9** |

**Where the 1 -> 2 GPU loss is (profiled ms/cycle, rank 0; rank 1 the same within 1 ms):**

| item | 1 GPU | 2 GPUs | loss | source |
|---|---|---|---|---|
| total (profiled) | 49.7 | 62.4 | **+12.7** | |
| GPU busy | 42.9 | 52.3 | **+9.4** | |
| - Krylov operator kernel | 9.14 (85 x 107 us) | 16.96 (172 x 99 us) | **+7.8** | `ImplicitStencilOpPart` (rad_m1_implicit.cpp) under implicit_halo_overlap. See the bug below. |
| - SC sweep, Picard, rest | | | +1.6 | more planes and bands (VetPlane*, VetBand*) |
| GPU idle | 6.8 | 10.1 | **+3.3** | |
| - SC-sweep plane MPI wait | 0 | 1.25 (28 gaps) | +1.25 | VetPlanePost -> VetPlaneWait, hipDeviceSync (rad_m1_vet.cpp, m1-sctb) |
| - hydro boundary MPI wait | ~0 | 1.3 + 0.3 | +1.6 | PackAndSendCC -> RecvAndUnpackCC (bvals) |
| - halo pack -> on-rank copy | 0 | 0.39 (94 gaps) | +0.4 | M1EvtB().record() between the two kernels (ImplicitHaloMPIPost) |
| - Krylov reduction round trip | 3.06 (84 gaps, 36 us) | 3.09 | **0** | the fence after each operator reduction, then the host scalars (and the MPI_Allreduce on 2 GPUs), then the next launch |
| unprofiled (runs_4k, w2 / s1c) | 43.3 | 53.4 | +10.1 | |

**What the table says:**
- The Krylov dot/norm reductions cost 3.1 ms of GPU idle per cycle on 1 GPU too. The
  2-rank MPI_Allreduce adds nothing measurable. Fusing the dots or using
  MPI_Iallreduce would not recover the 1 -> 2 GPU loss.
- There is no per-iteration convergence check or timer fence on this path worth
  cutting. Together they are well under 0.5 ms of GPU idle.
- **77 % of the loss is GPU busy time in the overlapped operator.**

**The operator loses time twice:**
1. **A bug.** Under implicit_vimp_fold, `ImplicitStencilOpPart` added the unfolded
   `M1VimpRow` to the folded stencil.
   - The x2/x3 +-1 vimp terms were counted twice.
   - It read 10 extra vimp slots per cell.
   - vfold + overlap is the default for be runs on 2+ ranks with halo_mpi, i.e. every
     multi-GPU He box since m1-accmerge. It is not used on the cubed sphere (dhj),
     where halo_mpi is off.
   - The effect on the results is at the solver-tolerance level, because the Picard
     residual is checked with the same operator (CPU gate below). base, vfold,
     overlap vs no overlap, 2 ranks: KE2 5.2e-8. With the fix: 6.9e-11. The vfold-off
     pair gives 3e-10.
2. **The shell split.** Part 1 (the interior) is 108 us for 81 % of the cells. Part 2
   (the w=2 shell on all six faces) is 69 us for 19 %.
   - The shell's x1-end cells are 4 of every 84 in a row, which makes the loads
     uncoalesced.
   - Yet x1 never needs a shell here: the x1 faces are physical (coefficient 0) or
     on-rank.

## 2. Changes

1. **Bug fix (unconditional).** `ImplicitStencilOpPart` computes the row exactly as
   `ImplicitStencilOp`: the folded slots 19-24 under vfold, and `M1VimpRow` only when
   the vimp is not folded.
   - Every path without vfold is unchanged, bitwise.
2. **Switch `<rad_m1>/implicit_halo_ovl_faces`** (bool). It is read only when named,
   default false.
   - `ImplicitHaloMPIInit` records which faces of the pack (x1-, x1+, x2-, x2+, x3-,
     x3+) have any off-rank ghost, from the receive regions, edges and corners
     included. This is `hm_face`, the union over the pack.
   - Under the switch the interior box reaches every other face. The ghosts of those
     faces are in place before the interior kernel: the on-rank copy
     (ImplicitHaloDirect) is queued ahead of it on the same stream, and a physical
     ghost has coefficient 0.
   - The shell is then only the MPI faces, in full contiguous rows.
   - With no MPI face in the pack, the shell launch is skipped.
   - With the switch off, the shell enumeration is the old one (same cell order).

Files: `src/rad_m1/rad_m1_implicit.cpp` (ImplicitStencilOpPart, the switch),
`src/rad_m1/rad_m1_krylov.cpp` (hm_face), `src/rad_m1/rad_m1.hpp`.

## 3. Gates (CPU, login node)

Setup:
- 3-D box `cpu/box3d_be_x`, 84x32x32 in blocks 84x16x16, 30 cycles, L, halo_mpi.
- Lists: `scripts/cpu_g0.txt` (base) and `cpu_g1.txt` (new). Compared with
  `scripts/cmp.py`.

| gate | result |
|---|---|
| new vs base, 1 rank | **BITWISE** (hst) |
| new vs base, 2 ranks, overlap off, vfold on | **BITWISE** |
| new vs base, 2 ranks, overlap off, vfold off | **BITWISE** |
| new vs base, 2 ranks, overlap on, vfold off (switch off) | **BITWISE** |
| new vs base, 2 ranks, overlap on, vfold on | differs (the bug fix): KE2 5.2e-8 |
| fix: new overlap on vs off, vfold on, 2 ranks | KE2 6.9e-11, KE1 3.2e-10, totE 1.7e-12 (base: 5.2e-8, 4.2e-9) |
| faces on vs off, 2 ranks, vfold on / off | KE1 7.5e-11 / 2.0e-10, KE2 1.6e-10 / 6.2e-11, totE < 4e-12 |
| faces on vs off, 4 ranks (blocks 84x8x16) | KE2 1.7e-9, KE1 1.1e-10, totE 2.7e-12 |
| NON-CONVERGED | **0** in all 14 CPU runs |
| **restart, faces on, 2 ranks** (N = 30 vs 15 + restart + 15, `scripts/rgate.sh`) | **BITWISE**: rst data after `<par_end>` identical at the end; the hst rows common to both runs are identical (B lacks the t = 3.07 row and has one at the leg end, as for any restart) |

- The faces switch changes only the summation split of the reductions. It is round-off,
  like the overlap itself (runs_3y).
- The inner-iteration totals are bimodal: 1134-1157 vs 1289-1303. One solve hits the
  200-iteration cap in some runs, for example 1 rank, overlap off with vfold off, and
  faces on with vfold on. The cause is round-off sensitivity of that one solve, not
  the switch.

## 4. GPU timing (ms/cycle, 3 repeats, interleaved, same binaries)

Jobs:
- t12 = 11955903, 1 and 2 GPUs.
- t13 = 11956317, 2 GPUs: faces vs overlap off.
- t4a = 11956311, hydro.

Binaries: `bin/athena_base_gpu` (e0074273, md5 98ed0111) and `bin/athena_new_gpu`
(md5 ce8ae35d).

| arm | 1 GPU | 2 GPUs: base (overlap, bug) | overlap, fixed | **faces on** | overlap off |
|---|---|---|---|---|---|
| E | 36.0, 36.2, 36.2 (**36.1**) | 43.2 (runs_4k w2, same code) | 41.8, 41.6, 41.8 (41.7) | t12: 40.5, 40.3, 40.3 (40.4); t13: 39.9, 40.6, 40.3 (**40.3**) | t13: 40.2, 40.4, 40.4 (**40.4**) |
| V | 42.7, 43.2, 42.5 (**42.8**) | 53.4, 53.3, 53.4 (**53.4**) | 52.1, 52.3, 52.5 (52.3) | t12: 50.0, 49.7, 49.7 (49.8); t13: 49.8, 49.6, 49.8 (**49.7**) | t12: 49.3, 49.5, 48.3 (49.1); t13: 49.3, 49.3, 49.4 (**49.3**) |
| Vh | 59.4, 59.5, 59.8 (**59.6**) | 73.3 (runs_4k w2, same code) | 71.7, 71.7, 71.5 (71.6) | t12: 68.2, 68.1, 68.0 (68.1); t13: 68.2, 68.2, 68.0 (**68.2**) | t13: 67.5, 67.2, 67.3 (**67.3**) |
| hydro | 21.0, 21.0, 21.0 (21.0) | | | | 21.6, 21.6, 21.7 (21.7) |

- NC = 0 in every arm.
- Picard counts are identical across the 2-GPU columns. Inner iterations per step are
  within 1.5 %.

**The 1 -> 2 GPU loss:**

| arm | before (base) | faces on | overlap off |
|---|---|---|---|
| E | +7.1 | +4.2 | +4.3 |
| V | +10.6 | +6.9 | +6.5 |
| Vh | +13.7 | +8.6 | +7.7 |

- The bug fix alone gives -1.1 to -1.7 ms.

**The halo overlap does not pay on 2 GPUs.**
- Even with the faces-only shell, **overlap off** is as fast (E) or faster (V -0.4,
  Vh -0.9) than faces on.
- Profile `runs/pf/w2_V_fac_1` (job 11955905): interior 109 us plus shell 32 us per
  operator, against 107 us for the whole operator on 1 GPU. The exchange it hides is
  shorter than the 30+ us that the second launch and its reduction cost.
- The w2 Krylov halo message is one 84 x 52 x 2 face per block. It is small.

**4 GPUs:** t4b = job 11956312, 2 apu nodes, E/V/Vh with overlap vs faces, plus hydro.
It was **PENDING** (apu full; apudev had one node down) when this was written.
Analyse: `python3 /viper/ptmp2/jinma/sync_0924/scripts/tsum.py t4b`. It has no
overlap-off arm. If 4 GPUs matter, add one (copy an `_ovl` line of
`scripts/arms_t4b.txt`, set `implicit_halo_overlap=false`).

## HANDOVER

**Done (commit on m1-sync, not merged):**
- The ImplicitStencilOpPart vfold bug fix. It is unconditional; paths without vfold
  are bitwise.
- `implicit_halo_ovl_faces` (default off, read only when named).
- Gates: see section 3. GPU timings: see section 4.

**Recommendation (not done, it is a default change for the user to decide):**
- Set the `implicit_halo_overlap` default to **false**. It is the fastest 2-GPU
  configuration measured (V 49.3, E 40.4, Vh 67.3 ms/cycle).
- The overlap code path, with the fix and the faces switch, stays available for
  larger faces or slower networks.

**What is left of the 1 -> 2 GPU loss (V, about 6.5 ms/cycle):**
- about 1.6 ms of extra SC-sweep kernels (more planes and bands);
- 1.25 ms of SC-sweep plane MPI waits (VetPlanePost -> VetPlaneWait, rad_m1_vet.cpp,
  the m1-sctb fences);
- about 1.6 ms of hydro boundary MPI waits (rank skew arriving at the hydro exchange;
  hydro alone loses only 0.7);
- the operator + halo on 2 ranks against the fused 1-rank operator.

**Not a lever:** the Krylov reduction round trip. It costs 3.1 ms/cycle of GPU idle on
1 GPU and on 2 GPUs alike. Cutting it (implicit_krylov_dev on several ranks, or the
pipe) helps both, not the scaling.

**Pending:** job 11956312 (t4b, 4 GPUs on apu).

## Files

| file | content |
|---|---|
| `scripts/gaps.py` | GPU-idle gaps by (kernel before -> after, host call), and the host wait split into GPU-busy and GPU-idle |
| `scripts/arms_t12.txt`, `arms_t13.txt`, `arms_t4a.txt`, `arms_t4b.txt`, `arms_pf.txt` | GPU arms (run with launch_0923 `job.sh`, paths in sync_0924) |
| `scripts/cpu_g0.txt`, `cpu_g1.txt`, `rgate.sh` | CPU gates and the restart gate |
