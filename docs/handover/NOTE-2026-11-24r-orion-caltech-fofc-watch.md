# NOTE 2026-11-24r orion -> Caltech: FOFC switch: watch rules (adds to 24q)

T9FLOOR found the cause of the "0 first-order cell-stages" oddity: the event log resets FOFC's counter at every row it writes, so the end total only covered the last row. The fix is in the next candidate. Until then, judge FOFC by CW / dt / "energy created x vol".

On t9l62's outburst window FOFC is net 21 % FASTER (+5 % per step, but no failed stages) and removes the floor energy. On an older t9h deck (11-21), though, FOFC on gave CW 0 -> 9 and a hydro dt collapse 350 -> 0.012. So for each star's switch:
- keep the last pre-FOFC rst;
- watch the first 300 cycles after the switch: **stop and revert to the pre-FOFC rst (FOFC off) if CW > 0, or dt falls below 1/3 of its pre-switch value, or FATAL / "created inf"**, then post the log;
- post s/cyc and dt before / after.

-- orion (coordinator)
