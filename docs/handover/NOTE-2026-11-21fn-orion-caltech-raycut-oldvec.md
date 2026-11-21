# NOTE 2026-11-21fn orion -> caltech (PRODMON): ray-cut / OLDVEC / cycle-time text extracts -- READ-ONLY, small text

These come from the kept full run.log(.gz) of **t19l64** and **t19l64_corona**, plus the **live t9l62** log for contrast (read-only). Commit into `caltech_results/PROD_<star>/`:
- `raycut.txt`: `grep -E 'ray cut cycle' <log> | sed -E 's/.*cycle ([0-9]+) shell ([0-9]+) -> ([0-9]+).*/\1 \2 \3/'`
- `oldvec.txt`: `grep -E 'CLEAN OLDVEC' <log> | sed -E 's/.*cycle ([0-9]+) t2s ([0-9]+) eint<0 ([0-9]+).*E<0 ([0-9]+).*/\1 \2 \3 \4/'`
- `cyctime.txt`: every 10th cycle line (`grep -E '^elapsed=|cycle=' <log> | awk 'NR%10==1'`), giving cycle, time and dt.
Why: in production T19 the sc ray cut jumps between shell ~27 and ~84 every few cycles; the smoke at cycle 100 had 117-119. The corona saw-tooths 31 -> 95 and back, while T9 holds steady at 182-187. We want to know when the cut started diving and whether the KE growth follows it.
This is small text and can go out with the next PROD monitor refresh.

-- orion (PRODMON)
