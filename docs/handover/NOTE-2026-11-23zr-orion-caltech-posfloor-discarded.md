# NOTE 2026-11-23zr orion -> Caltech: yes, the 5.5e47 was in a DISCARDED stage; keep running

Your reading is right.
- The "implicit_pos_floor created" warning at 7605 was printed by the stage-2 attempt that then failed. A stage fallback redoes the whole stage with backward Euler from the stage start, so that attempt's state, including its floor raise, is thrown away.
- Your hydro.hst column 7 staying at 1.33501e47 confirms nothing was kept.
- It still matters: it shows the write-back's work term in a few low-density outburst cells reaching 1e16 against a radiation energy of 3e8. That is the bug we're fixing now (T9WK). We'll also make the warning say "(attempt discarded)" when it is.

Keep running per the USER. Stop rules as before: dt < 0.01, NaN / inf, FATAL. Also post at once if a pos_floor warning is NOT followed by a fallback / redo in the same cycle, since that one would be kept.

-- orion (coordinator)
