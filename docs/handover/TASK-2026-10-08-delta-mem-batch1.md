# TASK for the Delta runner: mem batch1 = BSG GPU memory per block (mem-1009, user 10-09)

The BSG true-repro run is to move to Caltech (2 x H200, 8 blocks per GPU; target <= ~15 GB per block). viper apudev
is jammed; the same test is also sent to DeltaAI (first result wins). **Budget: 1 node x 4 A100, <= 30 min,
<= 2 GPU-h; ledger line "mem".** No production.

## Builds (`build_delta.sh he_gpu <sha>`, same build dir/flags as the withdrawn batch3)
- BO = **98835d99a** (`hrup-bsg-1009`), BN = **af5147cc2** (`mem-1009`, fork). Record md5s.

## Files
Inputs: `/work/nvme/bivj/jma20/delta_1008/bsg_hrdet/files` (SETUP.sh of batch3; if absent run step 2 of
TASK-2026-10-09-delta-beam-batch3.md). Scripts: `git fetch origin bsg-files-1009` then
`git archive FETCH_HEAD docs/handover/mem-1009 | tar -x -C /work/nvme/bivj/jma20/delta_1008/bsg_mem` (mem_cuda.sh,
ana_mem.py; read the header of mem_cuda.sh).
Kokkos Tools memory-events (wanted; host-only build):
`cd /work/nvme/bivj/jma20/delta_1008/bsg_mem && git clone --depth 1 https://github.com/kokkos/kokkos-tools kt && cd kt/profiling/memory-events && make CXX=g++`
-> KT = .../kt/profiling/memory-events/kp_memory_events.so (if it fails, run without and say so).

## Job (1 node, gpuA100x4, 4 ranks, `--time=00:30:00`, account bivj-delta-gpu)
```
ARMS="mem_r4 base_r4 mem_c4" bash .../mem-1009/mem_cuda.sh <BO> <BN> <files dir> <out dir> $KT
```
(A100 40 GB: the 8-blocks-per-GPU arm does not fit; mem_c4 = full 16-block mesh, 4 blocks per A100, may OOM:
that is a result.) Smoke rule: if mem_r4 dies with anything but OOM in the first minute, cancel and report.

## NOTE (NOTE-2026-10-09-delta-mem-batch1.md on this branch)
Per-arm `rc= fatal= oom= wall= s/cycle= peak:` lines, the `bitwise base_r4 vs mem_r4` line, md5s, job id, GPU-h; for
each arm with KT the `ana_mem.py <one rank's .mem_events> 1283800 <blocks per rank>` table (top 40 lines). Keep run dirs.
