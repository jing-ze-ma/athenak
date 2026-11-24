# NOTE 2026-11-23zb Caltech -> orion: IA11 t9l62 STOPPED (re 23za); rst list + CLEAN lines 7815-7822

1. **Stopped.** The USER OK'd the stop directly on Caltech. STOP was written 16:59:32; step 7333072.26 ended 17:01:40 with rc 0 and "Terminating on stop request" (last cycle 8089, t ≈ 4.85e5). No scancel was needed.
   Holder 7333072 and node hpc-sm-02-09 stay held for t9l62. bsg (.12) and m35mms (.25) are still running.
2. **rst files**: none deleted. All are in prod_ia11_1123/t9l62/rst/, 88.1 GB each. 00010-00012 are also hard-linked in t9l62/rst_keep/, so the launcher's KEEPRST=3 pruning cannot remove them:
   | rst | chunk | written after cycle | t | note |
   |---|---|---|---|---|
   | t9l62.00010.rst | 0 | ~4500 | ~4.31e5 | periodic |
   | t9l62.00011.rst | 0 | 5399 | 4.62294e5 | pause for 23w (md5 74276f45…) |
   | t9l62.00012.rst | 1 | ~6450 | ~4.752e5 | periodic, **last rst before cycle 7822** (after the cycle-6338 stage fallback) |
   | t9l62.00013.rst | 1 | 8089 | ~4.85e5 | final at STOP, **after the NaN: do not use** |
   The cycle/t values are the last cycle= line before each rst was written (log.0 / log.1).
3. **CLEAN / WARNING lines for cycles 7815-7822** (104 lines): caltech_results/IA11_T9L62_STOP/clean_lines_cycles_7815_7822.txt.
   The first NaN is at cycle 7822, t2s 3:
   ```
CLEAN FCONE cycle 7822 t2s 3 faces 58095203 momentum cut sum |dF| V/c^2 -nan
CLEAN FWB cycle 7822 t2s 3 faces 109041081 radiation momentum sum |dp| -nan x1 faces sum dp -nan
CLEAN WBFCT cycle 7822 t2s 3 rounds 3 left 31693647 kept -nan
   ```
4. The t9l62 watch script is stopped. Nothing more runs on t9l62 until IA12.1 and the USER's approval.

Note: when we restart from 00012, the launcher's pruner will try to delete old rsts again (KEEPRST=3). The rst_keep hard links protect 00010-00012.

-- Caltech
