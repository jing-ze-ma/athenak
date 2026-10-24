# NOTE 2026-10-24 caltech -> viper: run12 STARTED on Caltech (start rule: please cancel the pending viper + Raven copies)

- **STARTED 10-24 06:53:44 PDT**, static phase from t=0 (problem/stream=false time/tlim=0.05), then stream in the same link.
  1 node x 4 H200, 4 ranks (2 MeshBlocks/GPU), inside held allocation 5274909 (~23 h left; continues in a 7-d 1-node holder).
- Binary: 958cac780 ry_per_accretor, CUDA 12.9, HOPPER90, -ffp-contract=off, MALLOC_ASYNC off (incremental tree builds_accr):
  /resnick/groups/carnegie_poc/jingze/builds_accr/bin/athena_gpu_ry_per_accretor_958cac78_nofma, md5 af41c336d8a9aa3bcd8db870f6491b2b
- Deck: acc3d_prod_D_B1v6_extrap2_caltech_abs.athinput (md5 a93b8f05, from dff29a35); tables = acc3d_prod_1018/files (md5 = bundle).
- Link script: acc3d_run12_move/link_caltech.sh (= link_viper.sbatch phases + stop rules, srun --overlap, keep newest 3 stream rst).
- Run dir: /resnick/groups/carnegie_poc/jingze/acc3d_run12_move/run12/{static,stream}, chain.txt.

## Port check PASS (static nlim 10, 4 ranks, check/)
| | Caltech H200 | Raven A100 31453560 (4 ranks) |
|---|---|---|
| t=0 dt | 1.416096577200850e-05 | 1.416096577200850e-05 |
| t=0 mass | 1.142252572723363e+06 | 1.142252572723363e+06 |
| t=0 tot-E | -6.402933515344165e+11 | -6.402933515344166e+11 |
| cycle 10 | time 1.416097e-04, Picard 3.379 / 4, NC 0, hesdirk2 BE steps 1 | same |
