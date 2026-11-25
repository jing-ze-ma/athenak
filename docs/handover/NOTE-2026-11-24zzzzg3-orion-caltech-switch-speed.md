# NOTE 2026-11-24zzzzg3 orion -> Caltech: SWITCH t9l62, m35mms, M35ZAMS to c4speed (USER go; gates passed)

Your 1d55c17df (100 cycles identical) and 921ab1d53 (continuation identical) pass; Orion's gates too. Switch, **one star at a time, t9l62 first, then m35mms, then M35ZAMS; bsg is NOT touched.** Binary athena_c4speed_d979d5268 (md5 3627c9a8, checked before srun). For each star:
1. Stop cleanly (the STOP mechanism you use; it writes the final rst: no model time lost), after its newest periodic rst is finished; check the final rst exists, its size is right and it is not being written; keep it and the previous one (hard links in rst_keep). Note the cycle, t and the wall time lost.
2. Restart from that rst on c4speed in the SAME production directory / chain, same holder and node, the same args as now (every key named; no deck or key change; `time/cfl_number=0.3`, `clean_conv_max=1e-5`, the candidate keys, per-star mg levels); if your launcher stores the binary name, keep the chain's restart logic on c4speed afterwards (and the flag 'c4 at the next restart' of bsg stays as is).
3. Check the first 5 rows: t, dt, L_top, KE_int, M_tot, Etot against the last rows before the stop (continuity; the final rst's own row), CW / FATAL / BE / fallback 0, the log's "not in the restart file" lines (none expected), s/cyc over the first 100 cycles against the old value.
4. Post each star as it is done (start line with md5, holder, node, dir; the 5 rows; s/cyc) and only then go to the next. **If anything is wrong (FATAL, CW, a row that does not continue, no final rst): don't relaunch on c4speed; restart the star from the kept rst on its OLD binary at once and post.**
5. Watch (watch.sh) as before on all three. If the t9l62 burst has begun when you get there (dt < 48 s), switch it anyway; the rst is the same.

-- orion (coordinator)
