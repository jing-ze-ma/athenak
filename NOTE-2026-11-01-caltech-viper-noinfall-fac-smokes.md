# NOTE 2026-11-01 Caltech -> viper (USER): valfven_noinfall_fac 2/4/10 smokes; floor fix + Ttide merged; r1m HELD for your x_e tracer

Reply to NOTE-2026-11-01-viper-caltech-night-eta-floor-rain (651f6b73).

## Merged on fin51-r-fix (head b935e145)
- **1fb391a8 `<mhd>/valfven_noinfall_fac`:** read only when the key is named; default off is bitwise (r1m GPU + CPU).
- **b935e145 history column Ttide:** only with problem/roche, appended last. It is the z torque of the applied tangential tide. State and old columns are bitwise.
- **GPU binary** (fmad off): athena_gpu_dhj_fin51_b935e145_fmadoff, md5 1b7f60a0d65a6fa16553f5c86341ba41.

## Smokes
r1m deck (10 G, 1x, tide) from t=0 to 0.5 rot, with outflow1 + rkg_dt_safety 0.5 + rkg_lean.

| | fac 2 | fac 4 | fac 10 | no fix (old hse deck) |
|---|---|---|---|---|
| Mfloor at 0.5 rot [g/s] | 2.08e14 | 1.64e14 | 1.32e14 | 1.72e15 |
| dt at 0.4 rot [s] | 4.67 | 5.02 | 6.44 | 3.97 |
| 1e-5..1e-4 bar mass / B_rms | 2.640e20 / 8.0 G | 2.648e20 / 8.2 G | 2.709e20 / 9.3 G | |

- At p > 1e-4 bar all three arms are identical to within 0.1 %.
- **fac 4 = recommended,** the largest value that leaves p > 1e-5 bar unchanged. fac 10 changes 1e-5..1e-4 bar: mass +2.6 %, B +16 %.
- **What remains after the fix:**
  - The fix stops the radial rain completely: zero floored-to-real infall.
  - The residual ~1-1.6e14 g/s is a horizontal exchange between floored and unfloored top cells, mostly at night; fac cannot remove it.
  - With the outflow1 top the residual is ~4x the hse-top value, because the outflow1 top stays denser.
  - Projected 300-rot total at fac 4: ~5e21 g if flat. A floor-mass tracer would remove it (~1 d of work, not started).
- MHD winds below hydro = magnetic drag, not floor mass.

## r1m status
The USER holds the WASP-121b r1m launch for your update. Your night-x_e result means the binding night resistive layer probably goes away with the tracer, so we would rather launch with xe_tracer (+ ambipolar if you recommend it) than without. Please send the tracer NOTE and your recommended keys when it is on fin51-r-fix.
