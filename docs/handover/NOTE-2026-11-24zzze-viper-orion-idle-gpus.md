# NOTE 2026-11-24zzze viper -> orion: USER: the idle viper GPUs are yours (3 nodes)

The USER (11-24 ~18:55 PST): give the idle viper GPUs to Orion. Two viper allocations started 11-24 17:17-17:35 PST and run
only keepalive scripts. Orion may use them for anything (tests, arms, a production link):

| holder | name | nodes | GPUs | ends (PST) | pool queue |
|---|---|---|---|---|---|
| 13610529 | pr136f | 2 (vipa1310, vipa1324) | 8 MI300A | 11-25 ~17:17 | yes, pool_1107 holder |
| 13610530 | pool   | 1 (vipa1114)           | 4 MI300A | 11-25 ~17:34 | yes, pool_1107 holder |

- Both run `/viper/ptmp2/jinma/pool_1107` keepalives: drop a task script into `pool_1107/queue/` (the existing format; the holder takes
  it at the next keepalive tick) or use `srun --overlap --jobid=<holder> ...` directly.
- The only queued task, `05_t2109.sh`, rejects both ("wrong holder (nodes 1/2)", rc=9, delayed requeue). Retarget it or remove it
  as you see fit.
- R136 has no viper production link in pr136f (Raven is in maintenance until 11-27 22:55 PST). If you want R136 to continue here, that's
  your call, but keep it within the USER's R136 plan.
- Untouched (production, already Orion's): HG262 in 13587215, HG in 13587214, RSG in 13585068 (old scheme from rsg.01809, see zzzd).

-- viper
