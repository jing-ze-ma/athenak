# NOTE 2026-11-07 DeltaAI -> viper (+Delta): facts for the default-candidate tests (956be9b35) and holders/prod (c0da11b5)

Facts only. Nothing started, nothing queued (only `sbatch --test-only`); waiting for a TASK.

## 1. Build (both NOTEs, item 3/4): prod-1107v BUILDS on DeltaAI
- 698b8241329e (origin/prod-1107v), PROBLEM=he_star_m1, CUDA **12.9** (NVHPC 25.5), Kokkos HOPPER90 (sm_90), MPI ON (cray-mpich 9.0.1),
  host `-ffp-contract=off` + device `--fmad false` (target `hesnf_gpu`, same as the Oct RSG E / HG builds).
- **No nvcc "unexpected number of captures" error**; only unused-variable warnings. Fully compiled and linked (not just compile-only):
  `/u/jma20/ATHENAK/builds/bin/athena_hesnf_gpu_698b8241329e`, md5 `136cf6db3ce6e79950c0e0ca7177c09a` (11-07 13:44 CST, login node, 502 s).
- Not run yet (no GPU job).

## 2. Queue (`sbatch --test-only`, 11-07 13:45 CST; ghx4 has ~2280 pending / ~140 running)
| request | est. start |
|---|---|
| ghx4-interactive 1 node / 1 h | 13:53 CST (~10 min) |
| ghx4 1 node / 1 h | 18:06 CST (~4.5 h) |
| ghx4 1 node / 24 h | 18:06 CST (~4.5 h) |
| ghx4 2 nodes / 24 h | 18:53 CST (~5 h) |
| ghx4 3 nodes / 24 h | 19:58 CST (~6 h) |
| ghx4 6 nodes / 24 h | 22:13 CST (~8.5 h) |
These estimates move a lot (in Oct the same requests were +2.8 d).

## 3. Budget and billing
- bivj-dtai-gh: **811 SU left** of 1001.
- ghx4: 1 SU per GPU-h = **4 SU per node-h** (whole node, 4 GH200). ghx4-interactive bills x2 = 8 SU per node-h.
- 811 SU = ~200 node-h on ghx4. Examples: RSG E 2 nodes x 24 h = 192 SU; HG 3 nodes x 24 h = 288 SU; 1-node test job 1 h = 4 SU.
  So the budget covers roughly ONE production node-day set (e.g. RSG E ~4 days, or HG ~2.8 days), or many short tests, not both.

## 4. Holder allocations (c0da11b5 Q1)
- Limits: ghx4 MaxTime **48 h**, MaxNodes per job unlimited (partition), no per-user running/pending job limit set in QoS qos_ghx4.
  ghx4-interactive: MaxTime 2 h, max 2 nodes, 1 running + 1 submitted job per user.
- An allocation bills for its whole elapsed time whether or not GPUs are busy (billing = allocated TRES x elapsed). An idle holder that
  exits after ~60 min idle costs 4 SU per node for that hour.
- Technically we have run this pattern here (Oct: chained links, srun steps inside the allocation, lock file, self-trimming
  TimeLimit to an SU budget). I found no NCSA rule against a job that runs queued srun steps; I have not asked NCSA.
- Given 811 SU, holders larger than 1 node x 24 h use up the budget quickly.

## 5. Scratch (shared quota with Delta, same project bivj)
- /work/nvme/bivj: 749 GB used of 1000 GB (~250 GB free). /work/hdd/bivj: 207 GB used of 1000 GB (~790 GB free). /projects/bivj: ~950 GB free.
- Tests (HG 6.7 + RSG 4.0 + R136a1 12.1 + deps 1.2 GB) fit on nvme. Production rsts (RSG 45.4, HG 97.5 GB) + outputs belong on hdd or /projects
  (Oct RSG E wrote ~29 GB/h of output + ~150 GB rst; nvme would fill in < 1 d).
- Inbound copy speed from Raven/Caltech: not measured yet (copies go through the USER).

## 6. Fit with D1-D4 and productions
- GH200 96 GB per GPU: a whole small state (hgb/rsb) fits on 1 GPU, so D1-D3 incl. the 1-vs-2-GPU layout gate fit in 1 node / 1 h
  (interactive starts soonest). D4 here = GH200/Grace (aarch64) cross-check vs Raven.
- Measured in Oct on DeltaAI: **RSG E, 2 nodes, 8 GPU x 3 blocks: 52.9 GB/GPU, 2.25 s/cycle** (Caltech 1.35 s/cycle on 8 H200), clean.
  Older note: without the memcut keys, 3 blocks/GPU ran out of memory on 96 GB (10-20).
- HG 8 x 3 (~75 GB/GPU) untested here; 12 x 2 on 3 nodes fits.

## 7. Data route
- USER copies; either source is fine. Raven `/raven/ptmp/jinma/caltech_1107/` is preferred if it already holds hgb/rsb + deps.
- Target dirs once a TASK comes: tests on /work/nvme/bivj/jma20/, production on /work/hdd/bivj/jma20/ or /projects/bivj/jma20/.
