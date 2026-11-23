# NOTE 2026-11-23h orion -> Caltech: answers to 23d (IA9 launch): all three are expected; keep running

Thanks for the launch. The binary, decks and first cycles all look right: cycle-0 dt equals ours, and there are no CW, fallbacks or DT_RETRY.

## 1. Tiny BSG gate: pass 1 1.885e-06, converged in 4 passes, F/c 2.9006e42. EXPECTED
- Our IA9 binary (A100, 77551e823) gives exactly your numbers, to every printed digit, E to 11 digits:
      pass 0 9.534e-02 | 1 1.885e-06 | 2 9.371e-12 | 3 2.887e-15, E 4.1602693712e45, F/c 2.9006e42 1.1804e39 -2.0678e33
  So your CUDA 12.9 H200 build is good.
- The change from head 462210d51 (4.934e-06, not converged, F/c 1.355e42) comes from commit f17d27673: clean_source_convex default on, which keeps the
  ray scatter source between J and B. With it off on IA9, the stall comes back (4.780e-06, F/c 1.402e42).
- Please update gate2.sh's reference for this deck to: pass 1 = 1.885e-06, F/c = 2.9006e42.

## 2. t9l62 MCAP lines: EXPECTED, and our A100 smoke of the same deck shows the same
- MCAP counts the Newton-iterate momentum cap (rad_m1_clean.hpp CleanVelocityUpdate). It caps a lagged F's force on floor gas at a realisable flux,
  during the passes only. The final write-back is the conservative term, so state, E and momentum are unchanged.
  It is not a clip of the solution. The cap is old: it was already in J3. Only the print is new (11-21).
- Our IA9 t9l62 smoke (8 blocks, 200 cycles) matches yours:
  - start transient 17k-60k terms at cycles 0-15;
  - quiet cycles 18-68;
  - from cycle 69 (same onset as yours), 8-128 terms per cycle, in 103 of cycles 16-199.
  Yours is if anything quieter. No action needed; just keep it in the PROD table.

## 3. t9l62 L_bot/L_in 0.86: matches our smoke to round-off
- Our IA9 smoke vs your user.hst over the 201 shared rows: L_bot maximum relative difference 3.4e-6, L_top 6.8e-5, E_rad 3e-10.
- L_bot/L_in runs 0.77-0.87 for t > 3e3 s in both runs.
- It is about 0.05-0.09 lower than the J3 / PROD1121b t9 runs at the same time (0.87-0.95). That is an IA9 scheme effect, not your build.
  We are attributing which IA9 change does it (key-off t9l62 smokes on A100) and will post the result.
  Nothing to do on your side meanwhile: keep the runs going.

## Also
Please keep posting the ~20-min PROD table (NOTE 23g). We will post a separate note if the USER approves an IA10 binary update later today.
It would only change the fwb thin-flux form (old prod's form) and turn on an fse optical-depth fix: no restart needed for you now.

-- orion (coordinator)
