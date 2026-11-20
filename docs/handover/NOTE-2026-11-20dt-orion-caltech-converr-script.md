# NOTE 2026-11-20dt orion -> caltech (USER): CONVERR analysis script (run once, commit 4 small text files)

The scripts are in `docs/handover/caltech_runs_1121/`: `converr.py` (numpy only; it imports the bin_convert.py already there) and `converr.sh`. They only read the CONVERR run dirs and write text files.
- **When:** after all CONVERR arms are done. For earlier results, it can already run once CEP / CET3 / CET2 (+ CET1, CEP0) of a star exist; it uses whatever runs exist.
- **How** (inside a holder; ~40 GB RAM, ~10-30 min per star):
  `srun --overlap --jobid=$HOLDER -w $NODE -N1 -n1 -c 8 -t 120 bash docs/handover/caltech_runs_1121/converr.sh /resnick/groups/carnegie_poc/jingze/orion_runs_1121`
- **Output:** `$OUT/CONVERR/converr_{rsg,hg}.txt` and `converr_{rsg,hg}_3c.txt`. Please commit these 4 files to `docs/handover/caltech_results/CONVERR/`.
- **After that commit, the CONVERR bins can be deleted.**

-- orion (for the USER)
