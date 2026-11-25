# NOTE 2026-11-25z caltech -> orion: merged next line e8d3e64d9, CUDA gates (your 25z)

## Builds (worktrees wt_orion1121/next2, skip; CUDA 12.9 + g++13; 02-16, 10 cores each)
- BUILD next2 rc=0 wall=1455 s TUs=250
- BUILD_OK next2 md5=f3c50eae1b2ebb30a02dd3f498722bf6 prov=e8d3e64d9fd77ad94de952fb906a1aa449072df1
- GATE PASS next2: athena_next2_e8d3e64d9 md5 f3c50eae node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 05:22:14 end: Terminating on cycle limit
- BUILD suite rc=0 wall=1376 s
- BUILD_OK suite md5=f9b272e9ec4c57183b716609a6fbc2eb prov=e8d3e64d9fd77ad94de952fb906a1aa449072df1
- BUILD skip rc=0 wall=1494 s TUs=250
- BUILD_OK skip md5=8c11b73936a0f71a90df93725f930aba prov=ab99949076adfbead8df203083ba385a5739b443
- GATE PASS skip: athena_skip_ab9994907 md5 8c11b739 node hpc-sm-02-16 pass1=3.686e-11 log date 11-25 05:22:57 end: Terminating on cycle limit
- Tiny gate deck as before (clean_face_hctau = false in <rad_m1>).

## HG (hegiant.00205, T6; ABBA 79c1cabef BASE / e8d3e64d9 NEW, 8 ranks 2 nodes, 26 cycles; MT_SKIP = ab9994907 with libmpit.so; ONE4 = e8d3e64d9 on 02-03, 4 ranks)
| arm | cycles timed | s/cyc mean | median | min/max | dt first/last | its (sum/cycle) | CW/FATAL/inf/DT_RETRY | end |
|---|---|---|---|---|---|---|---|---|
| BASE1 | 24 | 2.8703 | 2.7974 | 2.692/3.535 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 |
| NEW1 | 24 | 2.4765 | 2.4011 | 2.290/3.155 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0 FATAL 0 |
| NEW2 | 24 | 2.4322 | 2.3598 | 2.252/3.116 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0 |
| BASE2 | 24 | 2.8495 | 2.7739 | 2.667/3.568 | 90.000/90.000 | 31.9 | 0/0/0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 |
- base (both arms) mean 2.8599 median 2.7861; new mean 2.4544 median 2.3923; new / base mean 0.8582, median 0.8587; base arms 2.8703 vs 2.8495, new arms 2.4765 vs 2.4322 (mean)
- NEW1 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW1 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- NEW2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 user.hst against BASE1: 26 / 26 rows, printed-bitwise
- BASE2 hydro.hst against BASE1: 26 / 26 rows, printed-bitwise

| arm | s/cyc mean / median | peak MiB per GPU (min / max) | user.hst / hydro.hst vs BASE1 | CW/FATAL | end |
|---|---|---|---|---|---|
| ABBA/BASE1 | 2.8703 / 2.7974 | 5 / 76139 (10 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| ABBA/NEW1 | 2.4765 / 2.4011 | 53409 / 53839 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0  |
| ABBA/NEW2 | 2.4322 / 2.3598 | 1465 / 53633 (9 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ABBA/BASE2 | 2.8495 / 2.7739 | 60295 / 76139 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| MT_SKIP | 2.3880 / 2.3137 | 14591 / 53549 (9 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ONE4 | 4.0572 / 3.9414 | 104835 / 105343 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=156 s CW 0  |
- MT_SKIP: Waitall cln_hd_copy 234.4 ms per cycle (215.5/249.8 min/max rank, 535 calls), all Waitall 357.4 ms per cycle; the ab9994907 and the base timer: compare with T6r (Waitall cln_hd_copy 584 ms, 644 calls per cycle, 79c1cabef)
- MT_SKIP: 2026-11-25 05:32:41 rc=0 ended: Terminating on cycle limit wall=106 s CW 0 FATAL 0
- ONE4: 2026-11-25 05:35:17 rc=0 ended: Terminating on cycle limit wall=156 s CW 0 FATAL 0

## HG: memory table corrected, and the suite
- The memory minima printed in the HG table above are wrong (two nodes append to one mem.log on NFS, samples get lost or garbled). The peak (max) is right. Parsed again with a strict host / index pattern; minima below are still lower bounds where samples were lost:
| arm | s/cyc mean / median | peak MiB per GPU (min / max) | user.hst / hydro.hst vs BASE1 | CW/FATAL | end |
|---|---|---|---|---|---|
| ABBA/BASE1 | 2.8703 / 2.7974 | 639 / 76139 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| ABBA/NEW1 | 2.4765 / 2.4011 | 53409 / 53839 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=111 s CW 0  |
| ABBA/NEW2 | 2.4322 / 2.3598 | 1465 / 53633 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ABBA/BASE2 | 2.8495 / 2.7739 | 60295 / 76139 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=116 s CW 0  |
| MT_SKIP | 2.3880 / 2.3137 | 52633 / 53549 (8 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=106 s CW 0  |
| ONE4 | 4.0572 / 3.9414 | 104835 / 105343 (4 GPUs) | bitwise(26) / bitwise(26) | 0/0 | rc=0 ended: Terminating on cycle limit wall=156 s CW 0  |
- Read: peak per GPU 79c1cabef 76.1 GB (1000 MiB), e8d3e64d9 and ab9994907 53.8 GB, one node 105.3 GB. Speed (T6, ABBA): e8d3e64d9 / 79c1cabef 0.858 (2.454 vs 2.860 s/cycle mean; base arms 2.870, 2.850, new 2.477, 2.432). The ab9994907 arm with the MPI timer: 2.388 s/cycle, Waitall cln_hd_copy 234 ms per cycle in 535 calls (79c1cabef T6r: 584 ms, 644 calls), so the fill's exchange counts fell by 17 %. ONE4: 4.057 s/cycle (median 3.94), hst bitwise with the 8-GPU base.
- The e8d3e64d9 suite: my first run used the 79c1cabef tests (wrong worktree variable), so it had no star_cs_ring / ray_cut_slab; stopped after 25 of the old tests had passed, restarted 05:36 with the e8d3e64d9 tests. The logs of the wrong run are kept under b1_1124/suiten2_oldtests and runs_suiten2*_oldtests.
