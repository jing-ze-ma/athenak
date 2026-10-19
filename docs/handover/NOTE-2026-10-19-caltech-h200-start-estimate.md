# NOTE 2026-10-19 caltech -> raven: H200 start estimates (reply to ce6272a4)

`sbatch --test-only` at 10-19 13:11 PDT (22:11 CEST), partition gpu, 24 h, 4 H200/node, 120 GB host mem/GPU, --exclude=hpc-sm-01-09:
- 2 H200 (1 node):  start 10-19 16:01 PDT (~3 h from now)
- 4 H200 (1 node):  start 10-20 22:32 PDT (~33 h)
- 8 H200 (2 nodes): start 10-20 22:31 PDT (~33 h)
Backfill estimates; they can move either way. Our queued acc3d trc (2 H200 check + 6 links) and the rdiff job count against our fair share.
Total GPUs: no per-user/account GPU cap set on carnegie_poc (assoc GrpTRES empty); partition gpu has 4 H200/node, ~18 H200 nodes,
roughly 16 fully allocated + 20 mixed right now. So 8 H200 is feasible, at the same ~33 h as 4.
