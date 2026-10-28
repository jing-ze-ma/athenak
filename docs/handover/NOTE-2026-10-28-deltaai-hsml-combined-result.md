# NOTE deltaai -> Delta + viper (USER): result of hsmlfix + probe job 3466654 (for the record; moot after Delta PASS b0301516)

- Job 3466654 ran 10-27 04:33-04:51 CDT on gh[104,116], 2 GH200 nodes, 1071 s, about 4 SU.
- All runs: `MPICH_SHARED_MEM_COLL_OPT` unset, `AREPO_HSML_DIAG=2`, `GpuMaxMemSize 50000`. Each stops at the first Sync-Point or TERMINATE.
- **1) fix 2d3967b136, 576 ranks: FAIL. This is GPU setup, not Hsml.**
  - Density converged in 17 iterations (72 left on the last one).
  - HSMLDIAG: 16909 cells checked, **0 with Hsml > 0.2 box**, tree ngb == brute force (max rel 9.0e-16).
  - Then 8 tasks TERMINATE in `reserve_backing_memory` (gpu_memory_manager.cu:297), on the 'vcell_respw' allocation at voronoi_cell_gpu.cu:1639.
  - Cause: 50000 MB requested, only 48052 MB free per GH200 at 288 ranks/node.
  - So 576 ranks needs `GpuMaxMemSize` <= ~45000.
- **2) probe d9c0 + HSMLPROBE, 576 ranks: TIMEOUT at the 520 s cap. Not reproduced.**
  - HSMLPROBE: 0 cells with Hsml > 0.2 box at the initial guess.
  - The density loop was killed at iteration 10, with 4.57 M of 6.09 M still unconverged.
  - Old-code iterations took 5-14 s each; the fix's took 0.06-8.5 s.
  - The converged Hsml was never reached, so the probe neither reproduces nor rules out the wall-cell blow-up here.
- **3) fix 2d3967b136, 256 ranks: PASS (148 s), reaching the first Sync-Point.**
  - Density took 17 iterations; HSMLDIAG gave 0 cells with Hsml > 0.2 box and brute == tree (max rel 9.0e-16).
  - VCELL mesh built with no abort.
- **Verdict:** the fix runs clean on DeltaAI at 256 ranks. The probe did not reproduce the bug here, so this PASS is not independent proof of the fix. Delta's PASS (b0301516) stands as the real test.
- DeltaAI: nothing running or queued; job3 stays cancelled. Logs: /work/nvme/bivj/jma20/arepo2_gpu_1025/logs/a2hcomb.3466654.out and runs/hcomb_3466654/.
