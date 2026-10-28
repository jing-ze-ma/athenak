# NOTE 2026-10-27 Delta -> viper + deltaai (USER): hsmlfix-1026 PASSES on Delta at the failing 256-rank layout; g25_full 3-arm timing done (g2m 55.8 / c1m 62.9 / c0 67.2 s per top step)

Full reply to TASK 385edc70 (TASK-2026-10-27-viper-delta-hsmlfix.md). Follows the interim NOTE cdebc9cf.

**1. hsmlfix-1026 test, job 23238294** (4 gpuA100x4 nodes, interactive, 256 ranks, 64/node, 16/island, `AREPO_HSML_DIAG=2`, no `MPICH_SHARED_MEM_COLL_OPT`).

| arm | build | result |
|---|---|---|
| g2m | 2d3967b136 (bundle md5 f547eaaf, verify OK) | **PASS**: first Sync-Point reached after 112 s, 0 TERMINATE lines |
| c1m | 2d3967b136 | not run separately (the USER dropped it for the queue); it PASSes in the 3-arm job below |

- `HSMLDIAG` summary: **16909 cells checked, 0 with Hsml > 0.2 box; 0 with tree ngb != brute force** (max rel 9.0e-16).
- Density loop: **17 iterations**, the same count as Raven's gates. d9c0 at this layout took 26-27 iterations and converged to bad Hsml on ~29-30k lattice cells.
- The 14 IDs from d6aa24a9 were not printed one by one: this run printed only the two summary lines, with no per-cell `HSMLDIAG` lines. So no single converged value for them. All checked cells are below 0.2 box. The d9c0 runs had those 14 at 0.44-0.48 box.
- **Reading:** on Delta, the same layout, nodes and environment fail on d9c0 (V1, V4) and pass on hsmlfix. So the shared-buffer race in the top-leaf exchange (9307d33d49) is the cause of the abort here.
- An A40 cross-check (23239103, fix + d9c0 control) gave no data: both arms hit a 10 s cap. Ignore it.

**2. g25_full 3-arm timing, job 23248921** (fixed builds builds_hf g2m/c1m/c0 from 2d3967b136, 4 gpuA100x4 nodes, 256 ranks, interactive, 10 top steps).
"Top step" = 7275.58 s sim time. Steps 0 and 256 are the all-active steps. s/top step = (cumulative at step 256 - cumulative at step 0) / 10, taken from cpu.txt.

| arm | rc | total | step 0 (init + first all-active) | **s / top step** |
|---|---|---|---|---|
| g2m (GPU) | 0 | 660.8 s | 102.3 s | **55.8** |
| c1m (CPU node-local) | 0 | 753.3 s | 124.1 s | **62.9** |
| c0 (CPU classic) | 0 | 803.4 s | 131.3 s | **67.2** |

- GPU vs CPU on the same nodes: only **1.13x** (vs c1m) and **1.20x** (vs c0).
- **Delta is slower than both references:** Orion CPU old Arepo is 35.6 s/top step on 560 cores; Raven 4 GPU nodes is 20.9 s. The USER's target (~10x Orion) is far off.
- GPU memory high-water mark is 1.54 GB per task out of 32 GB.
- Delta now looks at the per-section breakdown (cpu.txt, g2m vs c0). Results come in a follow-up NOTE.

**3. V1-V4 (diag 23219733):** as reported in cdebc9cf. V1/V2/V4 FAIL, V3 TIMEOUT. V4 = d9c0 + `MPICH_SHARED_MEM_COLL_OPT=0` FAILs, so the env setting is not a workaround.

**4. DeltaAI:** Delta's PASS answers the "does the fix cure it" question at the failing layout. DeltaAI's combined job 3466654 is now optional. The USER decides whether DeltaAI cancels it.

**SU:**
- hsmlfix 23238294: 2.0 min x 16 GPUs x2 (interactive) = 1.08 GPU-h.
- A40 23239103: 0.26 GPU-h.
- 3-arm 23248921: 38.4 min x 16 x2 = 20.48 GPU-h.
- bsg-files-1009 total on Delta: 42.08 GPU-h. bivj-delta-gpu has 952 of 997 GPU-h left.

All changes stay local. No Arepo source or diffs in git.
