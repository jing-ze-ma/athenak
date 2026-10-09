# TASK for the Delta runner: mem batch2 = vet_gd_twin_lowmem gate on A100 (mem-1009 6c5d8fb2, user 10-09)

New opt-in key `rad_m1/vet_gd_twin_lowmem = true` (with vet_gd_twin_fuse): the twin is swept first in the main
vet_gd intensity array instead of a second array (the 2nd copy found by DeltaAI) -> ~half the band memory, expected
BITWISE identical. The same gate runs on viper (apudev) and DeltaAI; first answer wins (viper may withdraw this).
**Budget: 1 node x 4 A100, <= 30 min, <= 2 GPU-h; ledger line "mem".** No production.

## Builds (`build_delta.sh he_gpu <sha>`)
BO = **98835d99a** (built for batch1/batch3; reuse), BN = **6c5d8fb2c** (`mem-1009` on the fork). md5s.

## Files
Inputs: `/work/nvme/bivj/jma20/delta_1008/bsg_hrdet/files` (bsg_hr_dc5.athinput). Scripts:
`git fetch origin bsg-files-1009 && git archive FETCH_HEAD docs/handover/mem-1009 | tar -x -C /work/nvme/bivj/jma20/delta_1008/bsg_mem`
(use **mem_cuda2.sh**; read its header). KT = the kp_memory_events.so of batch1 if you built it (optional).

## Job (1 node, gpuA100x4, 4 ranks, `--time=00:30:00`, account bivj-delta-gpu)
```
ARMS="lm_r4:N:4:red:5:1 new_r4:N:4:red:5:0 base_r4:O:4:red:5:0 lm_r2:N:2:red:5:1 new_r2:N:2:red:5:0 base_r2:O:2:red:5:0" \
  bash .../mem-1009/mem_cuda2.sh <BO> <BN> <files dir> <out dir> $KT
```
(r4 = 1 block per A100, r2 = 2 blocks per A100; the fused r2 arms may OOM on 40 GB: that is a result.)
Smoke rule: if lm_r4 dies with anything but OOM in the first minute, cancel and report. AG Car is gated on viper
(not needed here).

## NOTE (NOTE-2026-10-09-delta-mem-batch2.md on this branch)
Verdict line first: `bitwise lowmem vs base/fused: PASS|FAIL, keys-off vs base: PASS|FAIL`. Then the per-arm lines
(rc, oom, s/cycle(2-4), peak per GPU, the "ragged per-shell band" log line), all `bitwise ... vs ...` lines, md5s,
job id, GPU-h. Keep run dirs.
