# NOTE 2026-11-20bf orion -> caltech (USER): HIGH PRIORITY layout-invariance check, 2 RSG arms (B_v20n), plus a compare script

The batch-4 LS arms gave different hst md5s for different rank counts. These two arms tell whether only the history reduction differs or
the state itself. Arms are in `arms_li.txt`:
- LI4_rsg: NR=4 on 1 node.
- LI8_rsg: NR=8 on 2 nodes, with UCX_IB_GPU_DIRECT_RDMA=no as usual for multi-node.
Both run 3 cycles from rsg.01295.
- Keep the bins and the final rst on your disk. armc deletes rst, so skip that step for these two arms. Check that the bins exist at the last
  cycle; outputs 2/3/6 are set to dt=950, and they must include hydro_w, m1 and m1_vet. If the names differ, set the matching outputN.
- Then run `python3 docs/handover/caltech_runs_1121/licmp.py $OUT/LI4_rsg $OUT/LI8_rsg` and put its output in
  `caltech_results/LI_cmp.txt`. It reports: hst differing columns, the first differing row and the max relative difference; then, per
  final bin and variable, identical or the count of differing cells with the max relative difference (MeshBlocks matched by logical location).
Start them next.

-- orion (for the USER)
