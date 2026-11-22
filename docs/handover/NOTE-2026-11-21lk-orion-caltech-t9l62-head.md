# NOTE 2026-11-21lk orion -> Caltech: USER — move t9l62 to the current head from the last good B_j3 restart (cfl 0.9, fwb on)

USER (11-21 ~23:05 PST), after your 22:46-22:51 CW notes (t9l62 log.8 stuck at cycle 33461: stage fallbacks 1-2 + backward-Euler not converged, dt_retry 3, last change 0.52 -> 0.41):
"resume from the last good restart" on the current head, cfl 0.9, fwb on.

## Do
1. Stop the B_j3 run on testH1 cleanly (STOP file), if it has not already died.
2. Restart from the NEWEST t9l62 rst written by B_j3 BEFORE cycle 33458 (the first CW).
3. Binary: branch `caltech-prod-462210d51` on GitHub jing-ze-ma/athenak (= Orion integration head 462210d51).
   - Build with CUDA 12.9 + g++ 13 (NOT CUDA 12.2.1 + g++ 11.5: that miscompiles the ray closure on H200).
   - Gate first: the tiny BSG smoke (pass 1 = 4.020e-06), as before.
4. Keys: the head's defaults + the rst's own physics keys (no physics change vs the B_j3 segment), plus on the command line:
   - `time/cfl_number=0.9`
   - `rad_m1/clean_conv_max=1e-5` (the head's new default; name it explicitly anyway)
   - fwb stays ON (head default); fcone on (head default); flim off (head default).
5. Watch the first few hundred cycles vs the B_j3 segment, and post as usual:
   - CLEAN WARNING / stage fallbacks / dt_retry,
   - dt, s/cycle,
   - your jh noise metric on the dumps (rms 2nd difference of log T, phi / r, 0.3-1.5 R_ph).
   If CW/dt_retry recur or the jh metric rises above the B_j3 level, stop and post; we revert to B_j3.

## Why this should work (Orion evidence on head 07e1e5543 / 462210d51)
- clean_conv_max 1e-5: t9h9 noise 0.0003 vs 0.18-0.20 at 1e-4 (A100 and MI300A; 0 T-floor cells; 15-26 % faster);
  T19-class wedge (prod T19 deck cut to T19HOT) ran 150 ks at 1e-5 with 0 CW / 0 fallbacks / 0 dt_retry (your T1 stall was the
  old flim + xdonor switch limit cycle on B_damp3; the head has flim off and shows no residual floor).
- Not yet tested on Orion: cfl 0.9 together with conv 1e-5 — so please report early.

-- orion (coordinator)
