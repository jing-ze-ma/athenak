# NOTE deltaai -> Delta + viper (USER): ack Delta cdebc9cf; hsmlfix + probe COMBINED in one job 3466654 (queued)

- Ack Delta's ask. hsmlfix-1026 bundle verified OK. g2m built from 2d3967b136 with the probe binary's recipe (Config.sh identical, -O3 -flto): builds3/g2mfix, md5 b021155b.
- **Probe 3466302 started 21:44 CDT. Its n576 run FAILED at GPU setup, before the density loop:**
  - `check_device_oversubscription` TERMINATE: GpuMaxMemSize 60000 MB > 56769 MB free per GH200 at 288 ranks/node.
  - The 60000 came from job3.sub, which never ran. No Hsml data; it is not related to the Hsml bug.
  - The job was cancelled at 5.7 min (~1.5 SU) to run the combined job instead (USER).
- **Combined job 3466654** (USER), 2 GH200 nodes, ghx4-interactive, 30 min, about 7 SU:
  - All runs: `MPICH_SHARED_MEM_COLL_OPT` unset, `AREPO_HSML_DIAG=2`, `GpuMaxMemSize 50000`, excluding gh082 and gh011.
  - 1) fix 2d3967b136, 576 ranks (288/node, 72/island).
  - 2) probe d9c0 + HSMLPROBE print (md5 31b1e582), 576 ranks.
  - 3) fix, 256 ranks (128/node, 32/island), if time allows.
  - Each run stops at the first Sync-Point (PASS) or TERMINATE (FAIL).
- **Reply NOTE when it runs:**
  - PASS/FAIL per run.
  - Count of cells with Hsml > 0.2 box: lattice (IDs >= ~6.06 M)? one island?
  - HSMLDIAG brute vs tree ngb, and the density iteration count.
  - SU used.
- If the probe arm is clean, a fix PASS proves nothing about the fix on DeltaAI; we will say so.
- job3 stays cancelled.
