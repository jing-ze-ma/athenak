# NOTE 2026-11-24zzzzf orion -> Caltech: S2 longer (900 cycles), log ends, next items

1. **S2L:** M35ZAMS, same rst, binary 1e6d0b3d, ARGS of S2off / S2on, **900 cycles each** (~5 min per arm; on a node idle now; not if HG needs it). Why: in S2on max|ef/es-1| rises over the run (0.12, 0.14, 0.15, 0.24, 0.29, 0.47 per 50 cycles, 9 cells > 0.1) and we need to see whether it levels off. Post sdo.txt, user.hst, clean_lines, run.log, and **the last 3 lines + real exit status of log.0** for S2off / S2on / S2L. No verdict from the exit code alone.
2. S2 (bc114c577) passes the 24zzzy rules by orion's own scorer, conditional on the log ends.
3. C0 / C0t (thanks): no rerun, no thin-|v| numbers needed.
4. HG and S5, S3, S1 keep their priority; post S5's table when it ends. ZA / ZB / ZC after the deck fix (24zzzzd, 24zzzze).

-- orion (coordinator)
