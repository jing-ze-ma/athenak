# NOTE 2026-10-21 caltech -> raven + viper: FINAL setup validated on BOTH decks; RSG E final production running

- Final setup: f55be371 (origin/rst-detsum-1021, md5 50772674) + HCOLL on + UCX_IB_GPU_DIRECT_RDMA=yes + per-GPU NIC pinning.
- RSG E deck (rsg_E_det N=32, 2 nodes x 2 H200): C2==C (201 rows + bins), R==C (100 rows + bins): PASS.
- HG global deck (hgglob_fresh_clip at N=32, 1 x 4 H200): C2==C (201 rows + final hydro_w/m1 bins), R (from cycle-97 rst) == C (103 rows + bins): PASS.
- RSG E final production 5240541 started 12:05 PDT (rsg_E_1020/prodF), clean at cycle 450 / 5.2 d. HG global final 5240544-46 (hgglob_1021/prodF) queued.
- TASK 11e4e0f3 (+ addendum 59d41b7a) for Viper/Raven stands: please run the same triplets on your machines.
