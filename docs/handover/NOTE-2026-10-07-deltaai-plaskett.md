# NOTE 2026-10-07 (DeltaAI): Plaskett q12_s1 / q12f_s1 per TASK-2026-10-07-plaskett-remote -- BOTH STARTED

- **Binary:** accretor-1006 ea04cd19 (clean worktree), ry_per_accretor, CUDA Hopper90 + Grace, md5
  f6552f3c69c55b67b87b085b98915133. Inputs md5 = plaskett-1007/MD5SUMS.
- **Layout:** 1 node ghx4, 4 GH200, 4 ranks (2 blocks per rank).
- **q12_s1:** job 3333330 (ghx4-interactive). Smoke PASS: dt cycle 0 1.956427e-06, time cycle 10 1.942194e-05,
  rc 0, 0 FATAL/nan, 1.6e8 zone-cycles/s. **Production started 2026-10-07 14:57:11 CDT** (21:57 CEST), ~67 cycles/s
  at dt 1.92e-6.
- **q12f_s1:** job 3333331 (ghx4). Smoke PASS: dt cycle 0 1.956427e-06, time cycle 10 1.942194e-05, rc 0,
  0 FATAL/nan. **Production started 2026-10-07 15:03:40 CDT** (22:03 CEST).
- Results follow in this NOTE when each arm ends. Data: /work/nvme/bivj/jma20/plaskett_1007/run/<arm>/.
- **q12_s1 DONE:** reached tlim 0.2291 (0.50 orbit) at 2026-10-07 15:29:24 CDT, cycle 127041, rc 0, 0 FATAL/nan, final dt
  1.70e-6, min dt 1.429097e-06 (no collapse), wall 32.1 min (2.70e8 zone-cycles/s on 4 GH200). Data: /work/nvme/bivj/jma20/plaskett_1007/run/q12_s1/
  (run.log, ryper.hydro.hst, ryper.user.hst, bin/ 21 hydro_w dumps, rst/ 12 files).
- **q12f_s1 DONE:** reached tlim 0.2291 at 2026-10-07 15:40:02 CDT, cycle 126357, rc 0, 0 FATAL/nan, final dt 1.84e-6,
  min dt 1.21e-6 (no collapse), wall 36.3 min (2.38e8 zone-cycles/s). Data: /work/nvme/bivj/jma20/plaskett_1007/run/q12f_s1/ (slurm logs as slurm-<job>.txt).
- 2026-10-07 DeltaAI: Plaskett page https://claude.ai/artifact/6LeTufDkTVZqbNXC6WALLj ; arms q12_s1, q12f_s1 t_final
  0.500 orbits. Through R_acc: M 8.29 / 8.65 of 17.85 rho_s Rsun^3 in; j (last 0.1 orbit) 0.486 / 0.471 j_K, cumulative
  0.670 / 0.649 j_K. Arms nearly identical (FOFC changes little). numbers_*.json in plaskett-1007/results/.
  Caveat: stream width = default c_s/Omega (2.8x wider than Ryu+2025 L1 width), see TASK-2026-10-07-viper-accretor-handover.
