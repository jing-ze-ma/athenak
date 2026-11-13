# NOTE 2026-11-12 18:15 PST Caltech -> viper: Caltech GPU queue status (USER asked to tell you)

## Status

Nothing is running on Caltech. All 5 holders are pending (reason Priority).

On 11-12 at 00:40 the USER cut all of them from 7 days to a **3-day limit**. The change was made with
`scontrol update`, so they keep their queue age.

They are sleep holders: Slurm ends each one after 3 days. A run inside a holder therefore gets at most 3 days of
wall time; anything longer needs a new holder or a chained job.

## Estimated starts

| job | holder for | nodes x H200 | est. start (PST) |
|---|---|---|---|
| 5840649 r1H1 | dhj round-1 runs | 1 x 4 | 11-13 15:30 (in ~21 h; probably a backfill slot, may move) |
| 6133591 rsgEH2 | RSG E continuation | 2 x 4 | 11-17 12:30 |
| 5840650 r1H2 | dhj round-1 runs | 1 x 4 | 11-21 03:20 |
| 6133592 hgGH3b | HG global continuation | 3 x 4 | 11-21 03:20 |
| 6133593 testH1 | test / agents | 1 x 4 | 11-24 03:20 |

## Why it is slow

Our fairshare factor is 0.178. The carnegie_poc account is using 2.1x its share, and usage decays with a 14-day
half-life.

Even so, our priority is already in the top ~10 of the ~200 pending gpu jobs. The real blocker is node count x
walltime, which is why the 3-day cut was made. r1H1's estimate moved from 11-21 to 11-13 after it.

## Plan

HG, RSG and R136 stay on Raven as planned. Tell us (NOTE) if you want a particular run moved into r1H1 when it starts.
