# TASK for viper (from Caltech, user 10-07): push the He giant small input files so Caltech can smoke a FRESH START

Caltech cannot reach viper (no ssh) and the user cannot copy `hegiant_deltaai_bundle` (2.2 GB) here. User decision: Caltech runs a
**fresh-start smoke** instead of the N897 restart smoke. **Smoke only, no production on Caltech.** Leave hegiant897 and the DeltaAI copy alone.

## Please push to this branch (`hegiant-opn-1007`, fast-forward, never force), under `docs/handover/hegiant-files-1007/`:

1. `ic_giant_own.txt` (`problem/he_ic_file`, md5 ddc72692...).
2. `rosseland_tops_hegiant_blend.txt` and `planck_tops_hegiant_blend.txt` (`he_opac_table` / `he_planck_table`, md5 c1ef063c... / b699c60d...).
3. The fresh-start input(s): `hegiant_scout128_N445_fresh.athinput` (the viper run.cfg `IN`) and, if one exists, an N897 fresh-start input;
   plus the run.cfg XKEYS / command-line keys used for the scout128 fresh start.
4. `MD5SUMS` for all files.
5. **A fresh-start reference smoke** so Caltech can compare (not just "it runs"): the N445 (or smallest sensible) fresh start, 60 cycles,
   same binary as the bundle (f3a66907, md5 fc7b33a1), with run.log, hst files, job output and run.cfg, plus the table that section 5 of
   TASK-2026-10-07-deltaai-hegiant gives for the restart smoke (dt at cycles 0/1/10/20/59, last hst row t/dt/mass/tot-E, Picard
   solves/mean/max/NON-CONV, `he_ic_balance` line, s/cycle with the layout). If 60 cycles of the fresh start is not representative
   (scaffold ramp, startup transient), say so and pick the cycle count.
6. Anything else the fresh start reads from disk (check the input up to `<par_end>`).

Reply as `docs/handover/NOTE-2026-10-0x-viper-hegiant-files-caltech.md` on this branch. Caltech builds f3a66907 (H200, CUDA) meanwhile.
