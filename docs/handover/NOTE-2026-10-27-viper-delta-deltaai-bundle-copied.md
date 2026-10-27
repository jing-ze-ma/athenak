# NOTE viper -> Delta + deltaai (USER): hsmlfix bundle COPIED to both machines; TASK 385edc70 is GO

- The user copied `arepo_hsmlfix_1026.bundle` (md5 **f547eaaf285bd307c72a32c5483801e2**) to Delta and to DeltaAI.
  - DeltaAI: `/work/nvme/bivj/jma20/arepo2_gpu_1025/bundle/`.
  - Delta: your bundle directory. If it is not where you expect, find it by md5.
- Run `git bundle verify` first (head 2d3967b136, prerequisite d9c0152bce).
- Then follow TASK-2026-10-27-viper-delta-hsmlfix.md / TASK-2026-10-27-viper-deltaai-hsmlfix.md:
  - Delta: one job, g2m + c1m, 256 ranks, `AREPO_HSML_DIAG=2`.
  - DeltaAI: only after probe 3466043, and only if it reproduces the bad Hsml.
