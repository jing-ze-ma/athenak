# vet_sc SC sweep: temporal blocking (`vet_mb_tblock`)

Date 2026-09-24, viper. Branch `m1-sctb` (worktree `/viper/ptmp2/jinma/wt_sctb`), base
rt-integration 5c9e4432. Only `src/rad_m1/rad_m1_vet.cpp` changes. This is TO-DO #1 of
`tests_m1/runs_3q_scscale/README.md`.

- **Binaries** (`/viper/ptmp2/jinma/sctb_0923/bin/`):
  - GPU: `athena_tb3_gpu` md5 61d19b9e (= the commit; also copied to `athena_new_gpu`
    before job 11954621) and `athena_ref_gpu` md5 cd17bce0 (HEAD 5c9e4432).
    `athena_v1_gpu` 5d65cae9 and `athena_tb2_gpu` c5c44aa0 are the earlier
    development states (no direct reads; tb2 adds the empty-pass skip).
  - CPU: `athena_tb3_cpu` md5 edee29d7 (the commit) and `athena_ref_cpu` md5 84b1c39d
    (HEAD). `athena_new_cpu` 2b91c4d8 is the first development state.
- **Scripts** (copies in `scripts/`): `build.sh`, `gate.sh` / `gate3.sh` (CPU gates),
  `rstgate*.sh`, `compare.sh`, `run.sh`, `time.sh` / `timex.sh` (interleaved GPU timing:
  repeat a in order, repeat b in reverse order), `profx.sh` + `vetprof.py` (rocprofv3),
  `tsum.py`, `arms_*.txt`.
- **Runs:** `/viper/ptmp2/jinma/sctb_0923/{cpu,gpu/runs}`.
- **Every GPU job** exports `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.

## 1. What changed

New keys, all read only when given (an input without them writes the restart it always
wrote). Everything is exact.

1. **`vet_mb_tblock = B`** (default 1 = off).
   - The launches of the sweep are cut into chunks of at most B layers. A chunk lies in
     one band-exchange group of `vet_mb_halo` (it ends at the exchange) and in one x1
     block. So B > `vet_mb_halo` acts as B = `vet_mb_halo`.
   - Each chunk is ONE kernel launch (`VetTBLaunch`):
     - one team per (sweep block, local ray, plane tile);
     - the team sweeps the chunk's layers in turn, with a team barrier between layers;
     - the intensity planes stay in the ring `ipl` (global memory, L2).
   - A ray only reads its own plane, so the teams never depend on each other.
   - The tile core is a piece of the region the chunk's last layer must fill. Layer l
     covers the core widened by (lb - l) reaches of the ray (`rdep`), clipped to the ray's
     valid overlap. That is exactly the upwind cone of the core. Neighbouring tiles
     recompute the cone cells they share with the same arithmetic
     (`VetRayK::Cell`, term for term), so every value written is the value of the
     unblocked sweep.
   - The moments stay in `VetMomLaunch`. A batch closes once `vet_mb_mom_batch` layers
     have accumulated (exact in any batching; see `VetMomLaunch`). The ring grows to
     max(nbat - 1 + B, B + 1) planes.
   - Without `vet_mb_halo` in the input, the default halo becomes max(3, B), narrowed
     to fit a block.
   - Excludes `vet_mb_lag`, `vet_mb_mom_fuse`, `vet_mb_angles` and
     `vet_mb_kernel = cell`. Composes with `vet_mb_agg`, `vet_mb_halo` and
     `vet_mb_agroup`.
2. **`vet_mb_tblock_tj`, `vet_mb_tblock_tk`**: the tile in x2 and x3 cells (0 = the whole
   plane). **`vet_mb_tblock_team`**: threads per team (default 256).
3. **`vet_mb_tblock_direct`** (default true).
   - The first launch after a band exchange reads the ghost band of a rank-local
     neighbour from that neighbour's plane (`VetRayK::Up`). The local band copy
     (`m1_vet_plane_loc`, 28 launches per call) is then skipped.
   - The copy is still made before an x1 block face, because the x1 move copies whole
     planes.
   - Remote slots read the band as received.
4. **`vet_mb_tblock_ovl`** (default true, only with MPI neighbours).
   - After the band is posted, the moment launch runs, and so do the next chunk's
     INTERIOR tiles: those whose cone reads no ghost cell. Then comes the wait and
     unpack, then the remaining tiles.
   - The host skips the interior pass when no (ray, tile) is interior. That is always
     the case with tiles only in x3 (they span the x2 ghosts).
   - Measured: the interior pass does not pay (sect. 3).

## 2. Gates

**CPU** (`cpu/gate3.out`, final binary edee29d7 against HEAD 84b1c39d; the direct-off arms
were rerun after adding the key to the inputs, same result; login node,
<= 12 ranks). The earlier binary gave the same results (`cpu/gate.out`, 2b91c4d8).
- **A. Off = bitwise.** Default arms, every output file (bin, rst, hst) byte-identical
  to HEAD, **10/10**:
  - 2-D slab, tlim 3: U1; Fx2r4; Ux2r1; Ux2b2r2; Ux1r4; Fx12r4; odd nx1 = 81 on 3 ranks;
  - 3-D box 84x32x32: T3 on 1, 2 and 4 ranks.
- **B. On = bitwise, 25 arms.** Each arm is compared with HEAD (hst) and with the
  tblock = 1 run on the same input.
  - Against HEAD, **both hst files are identical** in every arm.
  - Against the base run, bin/rst differ only in 1-6 bytes of the echoed keys (offset
    <= 9.8 kB).
  - NON-CONVERGED = 0 in every run.
  - The arms:
    - B = 3 on every arm of A (it acts as B = 2 where the default halo narrows to 2);
    - B = 2 with halo 3 (chunks shorter than the group), 3-D r4 and 2-D r2;
    - B = 4 with halo 4, and B = 4 with the default-halo rule (halo 4): 3-D r4 and x1 x
      x2 r4;
    - tiles tj x tk = 3 x 5 (T3r4), 4 x 3 (T3r1), tj = 2 (2-D), tk = 5 (T3r2);
    - team size 1;
    - halo 1 (chunks of one layer);
    - ovl off; agg off; direct off (T3r4 and Fx12r4).
- **C. Hybrid.** `vet_mb_agroup = 2` + B = 3 has hst identical to HEAD's
  `vet_mb_agroup = 2` (T3r4 with tj = 5, Fx12r4).
- **Restart** (`cpu/rstgate3.out`): the B = 3 Fx2r4 run restarted from `rst` 00001
  reproduces the continuous run byte for byte: bin 2, 3; rst 2, 3; the hst tails.
- **GPU (hst identity).** user.hst and hydro.hst are identical to the HEAD binary's in
  every arm of every timing job: 1, 2, 4 and 8 GPUs, all B, tiles, team sizes, halo 4
  and 6, direct and ovl on/off.
  - x1: 20/20; x2: 18/18; x3: 30/30; x4: 24/24; x48: 24/24.
  - The hybrid arms are identical to HEAD's hybrid.

## 3. GPU cost

**Setup** (as in runs_3q_scscale):
- 3-D He box 84 x 104 x 104, 8 MeshBlocks of 84 x 52 x 26, vet_sc full tensor, 64 rays.
- FAST switches spelled out in `gpu/he3d_fast.athinput`, `vet_mb_halo = 3`.
- 100 cycles; ms/cycle is measured over cycles 10-100.
- Repeat a, then repeat b in reverse order, both binaries in the same job.
- SC = the fenced "SC seconds per call". Every entry gives repeat a, then b.

**Jobs:**
- 1-2 GPUs: 11955089 (`runs/x3`) and 11955134 (`runs/x4`), apudev, binary 61d19b9e.
- 4-8 GPUs: 11954621 (`runs/x48`), apu, 4 nodes, binary 61d19b9e.
- Exploration with earlier binaries: 11954688 (`runs/x1`), 11954689 (`runs/x2`),
  binary 5d65cae9.

**Chosen setting:** B = 3 (= the default halo), tiles in x3 only: `vet_mb_tblock_tk = 13`
on 1 GPU and 9 on 2-8 GPUs.

**Strong scaling:**

| GPUs (nodes) | SC ms/call HEAD | SC off | **SC B=3** (tk) | ms/cycle HEAD | ms/cycle off | ms/cycle B=3 |
|---|---|---|---|---|---|---|
| 1 (1) | 3.92, 3.92 | 3.92, 3.92 | **3.08, 3.08** (13) | 54.5, 53.9 | 54.6, 53.9 | 55.1, 53.2 |
| 2 (1) | 4.73, 4.72 | 4.69, 4.68 | **4.04, 4.03** (9) | 44.0, 44.0 | 43.2, 43.7 | 42.4, 44.1 |
| 4 (2) | 4.67, 4.69 | 4.75, 4.69 | **4.26, 4.19** (9) | 37.0, 132.4* | 36.9, 36.6 | 37.0, 38.3 |
| 8 (4) | 4.27, 4.23 | 4.23, 4.24 | **4.02, 4.01** (9) | 33.6, 33.1 | 34.1, 34.8 | 33.9, 32.7 |

\* The 132.4 is a one-off stall of that run; the SC time of the same run is normal.

**Other arms** (SC ms/call, repeat a and b):

| arm | 1 GPU | 2 GPUs | 4 GPUs | 8 GPUs |
|---|---|---|---|---|
| B=3, whole plane (no tiles) | 3.61, 3.58 (x1) | 4.68, 4.66 (x2) | - | - |
| B=3, tk=13 | **3.08** | 4.22, 4.19 | 4.50, 4.44 | 4.19, 4.16 |
| B=3, tk=9 | 3.54, 3.54 | **4.04, 4.03** | **4.26, 4.19** | 4.02, 4.01 |
| B=3, tk=13, direct off | 3.24, 3.23 | - | - | - |
| B=3, tk=13, 512 threads / 128 threads | 3.58, 3.62 / 3.88, 3.88 | - | - | - |
| B=3, tj=26 tk=13 | 3.30, 3.29 | - | - | - |
| B=3, tj=18 tk=9 (interior overlap) | - | 4.84, 4.87 | 4.31, 4.27 | 4.10, 4.05 |
| same, ovl off | - | 4.21, 4.17 | - | - |
| B=3, tk=9, ovl off | - | 4.01, 4.04 | - | - |
| B=4 (halo 4), tk=13 / tk=9 | 3.31, 3.32 | 3.91, 3.92 (tk=9) | - | - |
| B=6 (halo 6, batch 6), tk=13 | 3.48, 3.47 | 4.09, 4.08 (x3) | 4.29, 4.21 | **3.83, 3.82** |
| hybrid G=2: HEAD | - | 4.13, 4.14 | - | - |
| hybrid G=2 + B=3, tk=13 / tk=9 | - | 3.81, 3.86 / **3.65, 3.65** | - | - |

Picard 2.700 per step and NON-CONVERGED 0 in every arm.

**Launches per call** (rocprofv3 kernel trace, job 11955135, `runs/prof3`, 10 calls,
`vetprof.py`):

| per call | 1 GPU HEAD path | 1 GPU B=3 tk=13 | 2 GPUs HEAD path | 2 GPUs B=3 tk=9 |
|---|---|---|---|---|
| ray launches | 84 x 25.6 us = 2.15 ms | **29** x 64.9 us = 1.88 ms | 84 x 16.6 us = 1.39 ms | **29** x 47.3 us = 1.37 ms |
| moment launches | 21 (0.70 ms) | 15 (0.57 ms) | 21 (0.64 ms) | 15 (0.49 ms) |
| local band copy | 28 (0.42 ms) | **0** | 28 (0.23 ms) | **0** |
| pack + unpack | - | - | 28 + 28 | 28 + 28 |
| all vet launches | 144 | **55** | 202 | **113** |
| span / busy | 3.84 / 3.69 ms | 3.02 / 2.86 ms | 4.64 / 3.04 ms | 4.03 / 2.63 ms |

The 29 launches are layer 0 alone, 27 groups of 3, and the last group {82, 83}.

**What this says:**
- **SC per call, B = 3 against HEAD:** -21 % on 1 GPU, -15 % on 2, -10 % on 4, -6 % on
  8 GPUs.
  - Hybrid G = 2 + B = 3 on 2 GPUs: 3.65 ms, -23 % against HEAD's spatial split.
  - B = 6 is the best on 8 GPUs (3.83, -10 %): it halves the exchanges.
  - ms/cycle does not move beyond the noise: the SC sweep is 6-12 % of the cycle.
- **Launches go down 3x, but the time only goes down 1.2x.** The blocked launch takes
  about as long as the three launches it replaces minus their gaps: 64.9 us against
  3 x 25.6 us on 1 GPU.
  - The per-layer launches were throughput-bound more than latency-bound. The barrier
    loop runs at 1 team per (block, ray, tile), 256 threads, a few cells per thread per
    layer: low occupancy.
  - Tiles help up to the point where the redundant cones cost more than the occupancy
    gains (tk = 13 best with 8 blocks per GPU, tk = 9 with 1-4).
  - Wider teams (512) and smaller ones (128) are slower.
- **Direct reads** remove the 28 local band copies per call: -0.16 ms on 1 GPU (3.24 ->
  3.08).
- **Overlap:**
  - The moment launch under the MPI band gains nothing measurable (tk = 9, ovl on/off:
    4.03 / 4.02).
  - The interior-tile pass (tiles in x2 and x3) is **slower**: 4.84 against 4.19 ms on
    2 GPUs. The two passes cost more than the MPI they hide.
  - The MPI gap is still 1.4 ms per call on 2 GPUs (span - busy), 28 exchanges. Only
    fewer exchanges (larger halo / B) reduce it: B = 6 is best at 8 GPUs.
- **SC still does not strong-scale:** 3.08 ms on 1 GPU, 4.0 on 2-8. The cause is the
  28 x (pack, fence, MPI, unpack) round trips, not the ray launches any more.

## 4. Recommendation and next steps

- **Recommended:**
  - `vet_mb_tblock = 3` with `vet_mb_tblock_tk = 13` for 8 blocks per GPU, and 9 for
    1-4 blocks per GPU.
  - On 8+ GPUs with 1 block per GPU: `vet_mb_tblock = 6`, `vet_mb_halo = 6`,
    `vet_mb_mom_batch = 6`, tk = 13.
  - Within one node: add `vet_mb_agroup = 2` (3.65 ms on 2 GPUs).
  - Default left off, as the brief asked.
- **Next steps:**
  1. **Cut the exchange count.** A band exchange every hk layers costs about
     50 us/exchange of host round trip on 2+ GPUs.
     - A GPU-aware one-sided or persistent MPI (`MPI_Send_init`, device-initiated) would
       remove the fence.
     - Alternatively, a larger hk only on the ranks with remote neighbours.
  2. **More occupancy in the blocked kernel.**
     - Split the rays of a team into sub-groups: two rays per team, 128 threads each,
       with the same barrier.
     - Or hold the upwind cone in LDS for the most oblique rays only (rd > 1). Their cone
       dominates the redundant work.
  3. **An automatic tile rule** (teams ~ 1000 per GPU), after checking it on more block
     counts.
