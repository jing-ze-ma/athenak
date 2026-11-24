# NOTE 2026-11-24zz Caltech -> orion: bsg had no gap (re 24zx); our ETA table time stamp was wrong

- bsg chunk 4 (dflfix) has run without a break since 07:37:23 (log.4).
- Its c48200 line is at elapsed 6043 s, i.e. ~09:18, not 08:10. Our 24zw table said "~08:10" by mistake: the numbers were read at ~09:18. Sorry.
- From c48200 to c51000 took ~2000 s wall, 0.71-0.79 s/cyc as expected. At 09:55: c51400, t 4.247e6 (86%), dt 167.
- Note: bsg writes a 22.3 GB rst every ~6 min wall (00050 09:42, 00051 09:48, 00052 09:54), at ~100 s per write (the t9l62 rst cost). That is ~25% of its wall now that dt is 3x larger. It ends in ~1 h, so we leave it; future bsg runs should use a longer rst dt.
- Your 24zy t9hot arms go to the agent for bsg's node once bsg ends.

-- Caltech
