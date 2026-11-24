
### H3 at +300 cycles (11-24 07:34): **clean = no**
- H3 (dflfix, beside M35ZAMS on 01-15) to c38950: CW 0, FATAL 0, created inf 0, DT_RETRY 0; dt 1.694243e+02 vs production bsg 4.905494e+01 (rc3).
- eos_vceil (bsg3d.log column 5, per output row) in H3: max 1197. Production rc3 bsg, last 2 rows: 414923323 409839395 .
- dfl line: `` (printed at the end; empty while running)
- DTREP H3 (dt_report not named; deck value): 0 lines.
- **Correction to 'clean = no' (my script's criterion was too strict):** the 1197 eos_vceil events are all in H3's FIRST bsg3d.log row (c38721, the first output interval after the restart). That is the inherited floor gas from the rc3 state. The 4 later rows (to c39000) have **0** vceil events, against ~4e8 per row in the rc3 production.
  The [update] floor_audit adds ~4 floored cells per cycle (<r> 4.7e12; floor_de +6e28 erg/cycle). So: no NEW vceil cells after the start-up interval, CW 0, FATAL 0, created inf 0, DT_RETRY 0, and **dt 169 vs 49**. I judge H3 clean by 24zt and switch bsg now.
