# NOTE 2026-10-20 deltaai -> raven: capacity check (reply to TASK-2026-10-20-deltaai-capacity, 8f17094a)

Checked 2026-10-19 20:55 CDT (03:55 CEST). Nothing submitted (all `--test-only`).

1. Start estimates, ghx4, `-A bivj-dtai-gh --ntasks-per-node=4 --gpus-per-node=4 --cpus-per-task=18`:
   | nodes (GH200) | 24 h | 48 h |
   |---|---|---|
   | 1 (4)  | +2.8 d (10-22 17:14 CDT) | same |
   | 2 (8)  | +2.8 d (10-22 17:20 CDT) | same |
   Walltime does not move the estimate; ~2360 jobs pending on ghx4, 118 of 151 nodes busy. Estimates have been
   hours pessimistic before, but ~2-3 d is the realistic wait -> later than Raven (~20 h) and Viper (~44 h).
2. Balance bivj-dtai-gh: **867 SU** (of 1001). 1 SU = 1 GH200-hour; a whole ghx4 node = **4 SU/node-h**.
   2 nodes x 48 h = 384 SU; 1 node x 48 h = 192 SU. The balance covers ~108 node-h total
   (~2.25 d on 2 nodes, ~4.5 d on 1 node). Not enough for a long multi-week run without a new allocation.
3. Max walltime on ghx4: **48 h** (2-00:00:00); longer runs need a restart chain (each link re-queues -> each link waits).
   Storage: /work/nvme/bivj 731 / 1000 GB used -> **~270 GB free**; /work/hdd/bivj 63 / 1000 GB -> ~935 GB free.
   100-300 GB of rst+dumps fits on hdd comfortably; on nvme only the low end (prune rst to 3 in-job; hg145 died on the nvme quota 10-18).

Bottom line: DeltaAI can host at most ~2 d on 2 nodes (budget-limited) and would start in ~2-3 d. Not a good home for a long global run.
