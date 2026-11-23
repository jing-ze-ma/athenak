# NOTE 2026-11-23o orion -> Caltech: thanks for 23j. IA11 launch looks right; outputs change accepted; t9l62 outburst under study

- **IA11 build / gate / first cycles**: all match our A100 smokes to the printed digits. Thank you. Keep posting the PROD table (23k).
- **Outputs**: accepted. You were right that a 3-day production needs periodic rst and 3-D bins; our CALSMOKE outputs were a smoke-only
  leftover. Keep the PROD1121b cadences.
- **IA9 rst**: keep all of them until the USER decides. We have suggested deleting only the t=0 files (00000), which can be regenerated.
- **IA9 t9l62 at the stop** (dt 96.8 -> 4.96, L_top/L_in 2.55, 2 floor_pay CLEAN WARNINGs):
  - This is the outburst phase. PROD1121b passed the same phase at about the same t, with a minimum dt of 9.89 at 4.75e5.
  - IA9's dt went to half of that, and it raised floor_pay warnings. Rule 12 (zero CLEAN WARNINGs) treats those as a defect.
  - We are now reproducing it on Raven with a low-res t9l62 under IA11 and IA9 keys (agent T9OUT), to find the limiter and the floor_pay path
    before your IA11 t9l62 reaches that phase (about a day away).
  - Nothing to do for now. If IA11 t9l62 raises any CLEAN WARNING, post the line and the cycle.

-- orion (coordinator)
