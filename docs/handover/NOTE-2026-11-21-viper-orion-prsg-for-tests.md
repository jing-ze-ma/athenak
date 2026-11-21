# NOTE 2026-11-21 viper -> orion (USER): viper pRSG 13359763 (4 nodes, 8 MI300A) is yours for tests

The USER asked viper to point pRSG at you for tests. The RSG production that had just relinked there is paused (killed
before its first cycle; newest rst rsg.01707 kept; `prod_1107/rsg/STOP` blocks relinks).

- **Allocation:** viper pRSG **13359763**, nodes **vipa1150, vipa1155, vipa1163, vipa1165** (2 MI300A each), started
  2026-11-21 08:28 CET = 11-20 23:28 PST, EndTime=2026-11-22T08:28:23 (CET).
- **Rules (as for the earlier pRSG / pHG loans):**
  - Run `srun --jobid=13359763 --overlap -w <nodes> ...` steps only. **Never cancel the holder** (its batch is
    alloc_holder.sh, which keeps it alive to wall - 10 min).
  - Do not touch `/viper/ptmp2/jinma/prod_1107` (RSG/HG/R136 production dirs).
- The holder itself only cycles a 1-node keepalive pool task (sleep + requeue, no GPU). All 8 GPUs are free.
- Please post a NOTE when you start and what you run, and one when you are done.

-- viper (for the USER)
