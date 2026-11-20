# NOTE 2026-11-21e orion -> caltech (USER): amendments to batch 1 and the v20g gate (please read before starting these arms)

1. **v20g gate** (`arms_v20g.txt`, updated; still ahead of batch 2/3):
   - The rt arms now use **NH=2**: leg b runs 2 cycles, writes the rst, restarts and runs 2 more.
   - New arm **V20GR_rsg** is a repeat of V20G_rsg leg a (repro check). Its hst must equal V20G_rsg/a bitwise, and so must V20GP_rsg.
   - Please give s/cycle for each arm (res.txt has it).
   - No R136 arm here (it needs a restart from Raven).
2. **CONSERV ctl build changed.** CS_rsg_ctl and CS_hg_ctl now use **B_conserv-ctl2** = caltech-runs/conserv-ctl2, commit
   5ca370ffcd7e98488f05bb15592cd4329f7196c3 (arms_batch1.txt is updated). If you already ran them with the old build
   conserv-ctl e973b3cbd, still send those results, and rerun with ctl2 when a node is free. CS_*_tie is unchanged.
   Return: run.log (CAUD / CLEAN LZ / PROF / EAUDIT lines) and the `*.user.hst` and `*.hydro.hst` files.
3. **LAG LW1170_hg / LW260_hg**: there is no need to send the bins. Run
   `python3 docs/handover/caltech_runs_1121/real.py $OUT/LW1170_hg $OUT/LW260_hg` (it needs numpy; bin_convert.py is next to
   it) and send its output. That is one line per arm: per m1 dump, counts of |F_r|/cE > 1.05 / 1.2 / 1.5 in the top 120 shells.
   After that the bins can go.
4. **LAG LV_* 300-cycle arms**: the full run.log is needed (gzipped), because the pass counts are read from it
   (`npass.py run.log`; you can also send its output).

-- orion (for the USER)
