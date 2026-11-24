# NOTE 2026-11-23za orion -> Caltech: STOP IA11 t9l62 now (dt collapse + NaN); keep the rsts

Thanks for the context blocks. This meets the alarm:
- NaN in the stage-3 FWB / WBFCT sums from cycle 7822;
- unpaid floor energy 1.8e41 at 7906;
- dt 3.45 -> 1.5e-3 at 7907, then write-back closure inf.

## Do now
1. Stop the IA11 t9l62 run cleanly (STOP file; scancel if it doesn't respond within a few minutes).
2. Do NOT delete any t9l62 rst. Post the list of rst files with cycle and time, so we can pick the last one before cycle 7822.
3. Post the CLEAN lines for cycles 7815-7822 (all CLEAN / WARNING lines, no hydro lines), so we can see the first NaN.
4. Keep the GPUs for t9l62. bsg and m35mms keep running.

## Next
IA12.1 is being gated tonight. It contains the write-back exchange fix for exactly this phase (on HG's local reproducer: 0 floor raises, 0 column take). We will send the branch + decks to restart t9l62 from the last good rst once the gates pass and the USER approves. No action until then.

-- orion (coordinator)
