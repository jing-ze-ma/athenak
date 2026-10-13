# NOTE raven (WV) -> Caltech: HOLD the wavefix-1013 suite (TASK-2026-10-13-caltech-wv-wavefix)

Please do NOT start it, or stop it if it is running (keep anything already finished; no report needed).
Raven already ran the BIT step for 243 non-ladder dirs and FULL a + b for the 194 sp and 29 Cartesian DIFF dirs:
- arm (a) fixes 1+3: the rlat family (rlat, rlat_h40, rlat_om5, rlat_T, rlat_om5_T; 23 dirs) BLOWS UP (NaN at the top
  cell near cycle 13), and rsw / rswv are much worse (fix 1);
- arm (b) + key 2: stable, but it regresses thin_u, nlte2 and nltek* (E_top/ex 1.000 -> 1.15).
The fixes need rework before any further suite time is spent. The coordinator will re-issue a TASK if needed.
