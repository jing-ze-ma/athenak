# NOTE 2026-11-20r orion -> caltech (USER): LW bins, one diagnostic before you delete them

Thank you for the real.py run. Its counts (max ratio ~1e8) look like an artefact of the script, not physics. Before you delete the
LW1170_hg / LW260_hg bins, please run (it needs a few minutes of CPU and numpy):
```
python3 docs/handover/caltech_runs_1121/real_diag.py $OUT/LW1170_hg
python3 docs/handover/caltech_runs_1121/real_diag.py $OUT/LW260_hg
```
Put the output in `docs/handover/caltech_results/LW_real_diag.txt` (one short block per bin: variables, block layout, percentiles).
It reads the first and last m1 bin of each run, and `bin_convert.py` next to it is used. After that, the bins can go.
Low priority: run it whenever a CPU is free; it is not a GPU job.

-- orion (for the USER)
