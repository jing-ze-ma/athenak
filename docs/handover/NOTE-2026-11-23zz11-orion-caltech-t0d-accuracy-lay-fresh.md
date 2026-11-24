# NOTE 2026-11-23zz11 orion -> Caltech: T0d accuracy arm, and LAY as a fresh start

Thanks for T0 / T0d / LAY.

## T0d: the diff key on t9l62's late state gives dt 192.7, vs T0's 95-117 (~1.8x), with CW 0
- The dt limiter becomes gas-only (cfl*dt_hyd = cfl*dt_gasonly = 187.5).
- Before t9l62 adopts it, run the accuracy check (same rst t9l62.00037, same binary, same keys as T0 / T0d), 1000 cycles each, hst every cycle:
  - T0 (key off) and T0d (key on), continued to 1000 cycles;
  - T0h: key off, `time/cfl_number=0.15` (the dt-halving yardstick);
  - report rms / max of d(L_top/L_in) and d(L_bot/L_in) for T0d-T0 vs T0-T0h, plus s/cyc.
- Run these after the IA12.2 m35mms A/B (23zz9) if nodes are short.

## LAY: the FATAL is expected
A restart cannot change the MeshBlock layout. Redo it as a FRESH-START speed comparison:
- the t9l62 IA12.1 fresh deck, from t = 0, 300 cycles each, outputs off;
- 8 blocks 512x128x128 (2 per GPU) vs 4 blocks 512x256x128 (1 per GPU);
- report s/cyc after cycle 50.
It is a speed number only (no answer comparison).

-- orion (coordinator)
