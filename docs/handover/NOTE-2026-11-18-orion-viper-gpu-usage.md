# NOTE 2026-11-18 Orion -> viper (USER): which viper GPUs Orion is using

FYI, so the viper session knows what Orion touched on viper MI300A and what is still queued. The USER approved all of it ("you can stop some prod runs"; it chose the HG262 task).

## Stopped production (USER-approved)
- **HG262 on viper (pool holder pHG262 13199032):** Orion stopped the athena step 13199032.0 at 01:33 CET, right after its restart dump hegiant.00046 (01:31, 67.8 GB, complete).
  - The task wrote `hg262_1113/run_viper/STOP` ("srun rc=137 in run.log.13199032"). **It is still there.** HG262 will not resume until it is lifted and 01_hg262.sh is requeued. That is the USER's decision; Orion has NOT touched STOP / CLAIM.d / DONE.
- **13199032 then died:** "CANCELLED DUE TO NODE FAILURE" on vipa1085 at 02:47 CET. Slurm requeued it, and it is PENDING again. When it restarts, it runs holder.sbatch and takes queue tasks as usual.

## What Orion used
- **Pool holder pHG262 13199032 (6 nodes, until it failed at 02:47 CET):** Orion stag2 test tasks:
  - full HG/RSG restart tests;
  - fresh R136 / RSG / HG comparisons, staggered vs old prod, he_ic_equil;
  - diagnostics.
- **Pool holder 13199033 (1 node, vipa1089):** Orion tasks since ~01:33 CET: the MI300A speed rows/profiles, MEMOUT, seam76, and the integrator's MI300A gate round (running now).
- **apudev:** Orion's small 15-min jobs (orion_ordP, orion_itb47; pending).
- **Not used:** viper CPU partitions (USER rule), and the prod link allocations pRSG 13199029, pr136f 13199030, pHG 13199031 (no steps by Orion).
  - NOTE: those holders do pick up Orion's queue tasks and bounce them back with rc 9 (each task checks SLURM_JOB_ID = 13199032/13199033 only). E.g. "holder 13199031 ... START 00a_orion_gate2.sh" in pool.log is such a bounce; it runs nothing.
- **Orion's own holders 13249507/08:** still PENDING; not used.

## Orion tasks still in pool_1107/queue (all named *_orion_*)
- **Run only on 13199033 / 13199032:** 00a_orion_gate2, 00e_orion_mi300a5, 00f_orion_memout3, 00f_orion_rstdiff3, 00g_orion_mi300a4, 00h_orion_solveport1, 00z_orion_itb47, 01b_orion_lp_rsg, 01e_orion_lprst, 02_orion_fresh_hg.
- Each has a STARTED guard and its own run dirs under /viper/ptmp2/jinma/orion_stag2/. Each deletes its restart dumps.
- **If they get in the way of viper's work:** delete them from queue/, or touch /viper/ptmp2/jinma/orion_stag2/HOLD_ALL and tell Orion.

-- Orion
