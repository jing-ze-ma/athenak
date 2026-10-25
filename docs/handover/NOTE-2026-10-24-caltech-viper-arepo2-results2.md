# NOTE caltech -> viper: Arepo2-IDORT on the Caltech H200 node, results 2 (2026-10-24 late)

Results only. No Arepo source, diffs or patches are in this repo; all Arepo changes are on LOCAL branches on Caltech
(/resnick/groups/carnegie_poc/jingze/arepo2_gpu): nltree-g25-caltechfix 31a5b39ff, idort-gpu-caltech-opt 0ded16763 / 30993e695 / c663a54c5.
Nothing pushed to Bitbucket. Full write-up on Caltech: arepo2_gpu/caltech/RESULTS_CALTECH.md.

## nltree_mns_3d_cuda CHECKS failure: fixed (test bound, not forces)
- The forces were correct; the CHECKS cross-check bound (0.05 x gross) does not hold for the relative opening criterion under
  HIERARCHICAL_GRAVITY: a fine-bin call sees few sources, so each accepted node may deviate by up to ErrTolForceAcc |a_old|
  (6.9e-5 at the failing call vs |diff| 8.2e-6 vs old bound 1.8e-6).
- Fix (CHECKS-only TU): under the relative criterion the bound adds ErrTolForceAcc |a_old| per interaction. Forces untouched;
  snapshots bitwise unchanged. nltree_mns_3d and nltree_mns_3d_cuda both PASS.
- nltree g25_red (32 ranks): gravity 2.72 -> 0.88 s per top step; physics agrees with the CPU tree (L(r) < 1 %).

## IDORT suite (examples/idort, 48 cases, GPU build, IdortGpu 1): no GPU-specific failure
- 45 PASS (incl. reruns). I2a / I3c: 1-rank sub-runs saw 4 GPUs -> pass with one visible GPU (as in the gates).
- C3 fails identically on the CPU path (check.py needs the Orion paper data). B1_dynamic_diffusion_levels exceeded the 30-min
  suite timeout without error (rerun projected > 4 h, stopped). DSA and ONEDIMS cases run the CPU solve by design.

## g25_full (6.1 -> 8.4 M cells), 1 node, 4 H200, 32 ranks, 4 islands; all levers gated 10/10 ALL SAME
| per top step | pure CPU | GPU baseline | GPU now |
|---|---|---|---|
| whole step | 234.6 s | 93.5 s | 83.4 s (68.3 s mean over 8 top steps) |
| IDORT total | 176.2 s | 56.4 s | 46.8 s |
| IDORT full sync | 143.8 s | 23.4 s | 19.5 s |
- Same-node speed-ups: full sync 7.4x, all IDORT 3.8x, whole step 2.8x.
- Helped (exact): direct writes into the island staging window (gpu io 24.6 -> ~12 s/run); lazy radiation time-step criteria
  for active cells + neighbours (env switch, default off; -5.5 s/top step). Run parameters (no defaults changed):
  IdortGpuMinCells 2e4, MeshDeviceMinCells 2000 (small gains).
- No gain: sparse state upload, per-task GS schedule, IdortDomainCostFine 4 (new parameter, default 1), cells/block 8 or 16.
- vs Orion (35.6 s/top step on 560 cores): 2.3x slower on 32 cores; ~55 s/top step is CPU-bound (partial steps 19 s,
  time-step criteria 7 s, Voronoi host stages + inter-island halo 15 s, other 14 s). 10x Orion is out of reach on this layout.
  Estimate with a full 128-core node: ~35 s/top step (not measured; job pending).
- Existing GPU paths: FORCETREE_NODELOCAL_GPU tree already resident (85 % sinks on device); VORONOI_CELL_GPU works but half
  the steps have < 2000 active cells (host builds) and the inter-island halo is host-staged (4.4 s/top step).
- Partial steps cost O(NumGas) per step regardless of the active set (<= 10 active cells still ~0.3 s): the next lever is an
  O(active) partial step (per-step workspaces, active-set halo plan, cheaper couple) = multi-day CPU refactor.
- Device memory: peak 47-51 GB/GPU at 8.4 M cells (~75-80 GB at 13 M); tree replica 0.8/1.1/1.7 GB at 6.1/8.4/13 M;
  full-mesh replication per GPU ~300 GB (does not fit).
