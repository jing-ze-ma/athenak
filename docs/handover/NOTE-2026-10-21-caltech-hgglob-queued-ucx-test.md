# NOTE 2026-10-21 caltech -> raven: HG global queued (held on the UCX test); TASK 1ba37054 test running

## TASK 1ba37054 (RSG UCX rdiff)
- Job 5224962 RUNNING (1 x 4 H200, hpc-sm-02-10; normal QoS, as the debug QoS rejects the walltime). Dir rsg_E_1020/rdiff_ucx, deck A (rsg_E_det N=32), 9a56b4a1.
  Variants a UCX_RNDV_THRESH=inf, b UCX_IB_GPU_DIRECT_RDMA=no, c control; each C, C2, R, R2; compare C==C2, R==C, R2==R + wall/zone-cycles.
- RSG E production ALREADY RUNNING since ~06:40 PDT (5223807 -> 08 -> 09, 2 x 4 H200) with UCX_IB_GPU_DIRECT_RDMA=no (NOTE 324ea964), ~2.9 s/cycle.
  It could not be held: link 1 started minutes after the requeue. If (a) is bitwise and (b) is not, I switch link.sh for link 2 (or restart early from the 5-d rst, USER call).

## TASK 1d352abb (HG global)
- Dir /resnick/groups/carnegie_poc/jingze/hgglob_1021. BIN 4de8168d 9a56b4a1. Deck template 4e6ea902 -> hgglob_fresh_clip.athinput md5 06001321
  (@FILES@ -> hgglob_1021/files; ic f9678c2a, rosseland e2868b75, planck 23ee543f, from docs/handover/hg145-prod-1018). run.cfg md5 a61047c8.
- link.sh = RSG E link.sh (same STOP rules, 2 x 4 H200, 3 blocks/GPU, 7 d) with UCX env read from run.cfg UCX_ENV at job start (now provisional
  UCX_IB_GPU_DIRECT_RDMA=no; set from the 5224962 result).
- Jobs: 5224982 (afterany:5224962, so it can't start before the test is done) -> 5224983 -> 5224984 afterok.
