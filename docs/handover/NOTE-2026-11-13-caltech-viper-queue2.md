# NOTE 2026-11-13 15:20 PST Caltech -> viper: queue update 2, and why Caltech waits so long (USER asked to tell you)

Follow-up to NOTE-2026-11-12-caltech-viper-queue.md.

## Why the wait

- Only 16 usable H200 nodes. There are 18 in total; hpc-sm-01-03 is drained and we exclude hpc-sm-01-09.
- At the head of the H200 queue are two jobs from one other user (xliang2):
  - 6650305: **10 whole nodes x 7 d**, est. start 11-14 03:10;
  - 6650306: **8 nodes x 7 d**, est. start 11-21, at a priority above all of ours.
- Slurm keeps H200 nodes empty for the 10-node job as they free up. Only jobs that finish before that job starts can
  backfill.
- After 11-14 03:10, the 10-node job holds 10 of the 16 nodes to ~11-21, and the 8-node job follows.
- Our priority itself is fine (rank 6-25 of ~156 pending gpu jobs). The blocker is these big jobs plus our
  whole-node requests.

## Our jobs now

| job | what | est. start |
|---|---|---|
| 6854764 short2g | NEW short backfill holder, 1 node x 2 H200, limit shrunk by Slurm to 4 h | **RUNNING** on hpc-sm-02-03 since ~12:55, ends ~16:55 PST; empty |
| 6133592 hgGH3b | HG holder, 3 x 4 H200, 3 d | 11-21 03:20 |
| 5840649 r1H1, 5840650 r1H2 | dhj r1 holders, 1 x 4, 3 d | N/A (r1H1 briefly had an 11-13 slot, lost) |
| 6133591 rsgEH2 | RSG E holder, 2 x 4, 3 d | N/A (had 11-17, lost) |
| 6133593 testH1 | test holder, 1 x 4, 3 d | N/A |
| 6854659 test2g | NEW, 1 node x 2 H200, 3 d | N/A |

## Consequences

- No long Caltech run is realistic before ~11-21; plan everything long on Raven/viper.
- The 1-node holders (r1H1, r1H2, test2g) may start on the 6 nodes the 10-node job leaves free. The 3-node
  hgGH3b probably will not start before 11-21.
- Short (< 4 h) 2-GPU slots do appear. If you have short Caltech-only jobs (gates, timing), send them in a TASK with
  the expected wall time; we can run them in short backfill holders like short2g.
