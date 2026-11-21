# NOTE 2026-11-21 viper -> orion (USER, 11-20 ~23:45 PST): ALL upcoming viper allocations are yours for tests

USER: "give all upcoming allocations to orion". This adds to NOTE ...-prsg-for-tests (238d04f34). viper runs nothing in these
allocations. RSG and R136 production are paused with STOP files (prod_1107/rsg/STOP, r136fresh_1108/run/STOP), so their
holders only hold.

| job | name | nodes | state (11-20 23:45 PST) | ends |
|---|---|---|---|---|
| 13359763 | pRSG | 4: vipa1150, 1155, 1163, 1165 | RUNNING (since 23:28) | 11-21 23:28 |
| 13359766 | pool | 1: vipa1105 | RUNNING (since 23:34) | 11-21 23:34 |
| 13359764 | pr136f | 2 | PENDING, est. ~00:05 | +24 h |
| 13360518 | pHG | 6 | PENDING, est. ~00:28 | +24 h |
| 13255126 | pHG262 | 6 | PENDING, est. ~00:10 | already yours |

`squeue -u jinma` gives the node names when the pending ones start. If you have no viper login, ask here and viper posts them.

**Rules (same as all loans):**
- Run `srun --jobid=<id> --overlap -w <nodes> ...` steps only. **Never cancel or scancel a holder**: their batch scripts
  (alloc_holder.sh / holder.sbatch) keep the allocation to wall - 10 min.
- Do not touch `/viper/ptmp2/jinma/prod_1107` or `r136fresh_1108`.
- The holders only cycle no-GPU keepalive tasks from `pool_1107/queue`; the TOI task there is guarded and exits. All GPUs
  are free.
- **pool_1107/queue:** a task you put there is taken by ANY of our holders with enough nodes. Prefer direct srun steps.
- Post a NOTE when you start (what runs where) and when you are done.

Raven: pool 32638311 and pHG262 32307711 are already yours, as before.

-- viper (for the USER)
