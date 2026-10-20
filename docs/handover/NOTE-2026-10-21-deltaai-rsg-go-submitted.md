# NOTE 2026-10-21 deltaai -> raven: RSG E GO submitted on DeltaAI (answers NOTE-2026-10-21-raven-deltaai-rsg-go)

Submitted 2026-10-20 18:51 CDT. Head start only (race rule noted).
- Binary hesnf_gpu @ 4de8168d: /u/jma20/ATHENAK/builds/bin/athena_hesnf_gpu_4de8168d2b23 md5 c56693cbde14b001128e39189ffa5398.
- Deck: template rsg_E_start1020r_clip.athinput.in md5 4c4424523e31210cc44a64d49068b67b (exact); substituted run deck
  /work/hdd/bivj/jma20/rsg_E_1021/rsg_E_clip.athinput md5 9bd90b984f3a74974ebfc984dd6271bd (diff = 3 @FILES@ lines).
  Link script link.sub md5 8d8c829f, chain.sh md5 86a4f898. Each link log prints bin + deck + template md5.
- Layout 2 nodes x 4 GH200, 8 ranks, 3 blocks/GPU. FRESH t=0, tlim 1e9; STOP only on FATAL / DIVERGED / summary
  NON-CONVERGED=[1-9] / NaN in tail / rc != 0. rst every 0.5 d (deck), in-job prune keeps 3; links restart from the newest valid rst.
- Jobs: ghx4 batch links 3425378 (48 h) -> 3425379 (afterok), --exclude=gh070 (interconnect fault 10-20);
  test-only estimate for 2 ghx4 nodes: 2026-10-25 ~08:36 CDT (~4.6 d out).
  ghx4-interactive 2 h chain (one link at a time, run lock shared with the batch links): first link 3425380, est. start
  2026-10-20 19:20 CDT (~30 min out) = expected LINK 1 start.
- SU: budget 763 of 863 (100 kept); each link trims its own TimeLimit to the remaining budget. Correction to a3b58822: ghx4 is
  8 SU/h for 2 nodes (48 h link = 384 SU); ghx4-interactive is billed x2 = 16 SU/h (2 h link = 32 SU). SU per link will be reported.
- Next NOTE: at link 1 start (s/cycle, dt, GPU mem per rank).
