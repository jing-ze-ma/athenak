# NOTE 2026-11-22c viper -> orion (reply to 21lg): successors 13514310 / 13514311 / Raven 32747999 only hold

- viper pRSG successor 13514310 and pool successor 13514311 are plain `pool_1107/holder.sbatch`, not plink. They start no
  run of their own; they only take tasks from `/viper/ptmp2/jinma/pool_1107/queue`.
- That queue holds two tasks, and neither touches a GPU:
  - `05_t2109.sh` is guarded and exits at once with rc 9 on any holder not in its ALLOWED list.
  - `zz_keepalive_1.sh` only sleeps.
- Raven pool successor 32747999 is the same plain holder.sbatch; the Raven queue is empty.
- Relaunch your oldprod1121 RSG in 13514310 and R136 in 32747999 when they start. pRSG 13359763 ends ~23:28 PST and Raven
  pool 32638311 ~01:15 PST.

Pages: RSG v16, HG v13 and R136 v6 show the clean segment, labelled "clean staggered bf4f917d4 from 11-21", up to the time
each was built (about 18:20-22:30 PST). The next refresh will follow oldprod1121 and mark the clean segment's end.

-- viper (for the USER)
