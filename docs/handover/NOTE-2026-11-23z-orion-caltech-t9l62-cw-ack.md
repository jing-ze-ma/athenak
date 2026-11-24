# NOTE 2026-11-23z orion -> Caltech: thanks for the t9l62 CW; keep running IA11

Thanks. One stage fallback (BE redo) at cycle 6338 with WBFCT left 0, and you're now past IA9's degradation point (t 4.83e5 vs 4.76e5) at dt 7.95. That is the behaviour we hoped IA11 would show. Keep it running.

Please:
1. Post the same short context block for every new CLEAN WARNING (first 5, then a count per hour).
2. Every ~2 h, post one line: cycle, t, dt, s/cycle, and counts per hour of stage fallback, dt_retry, CLEAN FCONE faces and CLEAN XT closure max.
3. Stop and post at once if dt falls below 1 or dt_retry repeats more than 3 times in an hour.

Not needed now: the full log. We may ask for the lines around cycle 6338 (±50 cycles, CLEAN lines only) for the write-back diagnosis.

IA12.1 (with a write-back exchange fix for this phase) is being gated tonight. We will send it with decks when it passes and the USER approves; no action until then.

-- orion (coordinator)
