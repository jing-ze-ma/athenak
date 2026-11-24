# NOTE 2026-11-23zz25 orion -> Caltech: watch production logs for "implicit_pos_floor created inf"

T9FLOOR found the cause of T0da's c30385 warning `implicit_pos_floor created inf of E`:
- a diverged stage iterate (e = -inf) was floored to e_floor;
- that hid the non-finite value from the rejection check;
- T0da's fallback discarded that attempt, but where no fallback follows, the step can be ACCEPTED. A CPU test reproduces it: dE 0.028.
- The fix is in the next candidate (rc2).

Until then, please add to the PROD-table poster for m35mms production (IA12.2 B): **grep the log for "created inf". If it ever appears, stop the run at once (STOP file / scancel), keep the last rst, and post the lines.** The production m35mms has had CW 0 so far, so it hasn't happened there.

-- orion (coordinator)
