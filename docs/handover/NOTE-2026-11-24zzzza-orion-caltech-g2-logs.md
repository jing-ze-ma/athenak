# NOTE 2026-11-24zzzza orion -> Caltech: ray-set G2 produced no cycle; logs to ship; RSG panel-edge number

1. **t9l62 G2 (13x36) did not run a cycle:** 62 s, no hst row, no dump, the log stops before "ragged per-shell band intensity array"; rc=0 is the wrapper's. Please post the last 30 lines of its stdout AND stderr, and the same for the +1 preflight (57 s: did it complete a cycle?). Then run the fallback named in 24zzzt: **7x36 (252 directions)**, else identical to G1, same 8 cycles, after checking that the first hst row appears. If it dies the same way, stop and post the logs.
2. **RSG thinf pair (thanks):** we also need the quantity thinf was set false for: the panel-edge maximum and the fraction of faces with |F| > cE (24zzzs), both arms, at the end of the 100 cycles (and at +1). If no dump was kept, say so; don't rerun yet. Which source commit is md5 09f50abe?
3. **C0t:** state its binary md5 and C0's. They must be the same binary to be a pair.
4. **sdirk2 star pairs (24zzzy):** ship the FULL clean log lines including every `SDO` line, the CW / DT_RETRY / BE / fallback counts, and user.hst of both arms. A run.log with only start / rc lines can't be scored. Both arms on the TEST binary a840d7e92 with SDO_PRINT=1.

-- orion (coordinator)
