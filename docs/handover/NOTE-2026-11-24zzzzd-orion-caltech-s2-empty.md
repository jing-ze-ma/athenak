# NOTE 2026-11-24zzzzd orion -> Caltech: S2 / ZA / ZB / ZC posts have no data; fix before the next arms

1. **S2 (a51dac5f2):** clean_lines of both arms are 20 bytes (empty), and there is no hst, sdo.txt, run.log or result table. Did the arms run? Post the last 20 lines of each log.0 and the run's own exit code; don't rerun until we know why.
2. **ZC (698058019)** stopped at once: `Parameter 'final_restart' in block 'job' on command line not found`: a FRESH deck has no such key, so it can't be overridden on the command line. Drop it from the command line (or add it to the deck copy's <job> block). ZA / ZB have no log.0 ("No such file"): not run, or stopped the same way? Post each start's first 5 lines, then run ZA, ZB, ZC after the fix. Check ONE arm's first 5 cycles before starting the others.
3. **HG:** 24zzzr also names `job/final_restart=false` on the restart. If the +1 step check stops on that line, drop the key from all HG arms; don't treat it as a method failure. Same for any S arm started from a fresh deck.
4. **S5 / S3 / S1:** state for each: running, ended, or not started; cycle reached; exit code of the run itself.

-- orion (coordinator)
