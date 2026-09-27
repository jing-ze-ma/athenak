---
name: dhj-ideal-vs-general-cost
description: "ANSWERED 2026-08-18 from the four long runs: general EOS is 3.09x ideal per SIMULATED SECOND on dhj (not 1.77x -- the dt advantage evaporates), and max_eta=1e13 beats 1e14 for BOTH EOSs."
metadata:
  node_type: memory
  type: project
---

Four long dhj runs in `/orion/u/jinma/ATHENAK/idflr/L_{id13,id14,gen13,gen14}` (ideal|general
x max_eta 1e13|1e14, all `use_rkg_sts=false`, whole node each, binary `idflr/athena_wb` =
`d1c289dd`). Submitted 2026-08-17 18:15, 23:30 wall -> they end 2026-08-18 ~17:45. **None
reaches tlim=1e7**; at 16 h the ideal pair is at t~2.4e6 and the general pair at t~9.5e5.
That is fine -- cost is a rate over a matched window, and the window is covered.

## THE RESULT (matched window t = 1e5..9e5, the largest all four cover)

| run | wall s / sim s | s/cycle | mean dt | vs ideal-1e13 |
|---|---|---|---|---|
| L_id13  (ideal, 1e13)   | 0.01971 | 0.0189 | 0.959 | 1.00x |
| L_id14  (ideal, 1e14)   | 0.02216 | 0.0192 | 0.865 | 1.12x |
| L_gen13 (general, 1e13) | 0.06085 | 0.0557 | 0.916 | **3.09x** |
| L_gen14 (general, 1e14) | 0.06573 | 0.0565 | 0.859 | 3.34x |

Stable across sub-windows (t=5e5..9e5 gives 2.99x / 3.22x). Extract with
`grep 'elapsed=' <dir>/stdout.txt` -- **the diag line carries wall `elapsed=`**, which is the
only wall-clock trace in the run (the .hst has none); bin-dump mtimes are a coarser backup.

## THE 10 G RUNS ARE DONE (cancelled 2026-08-18 11:31 / 12:13, "we've learnt enough")

Final reach: L_id13 t=2.590e6, L_id14 t=2.480e6, L_gen13 t=1.048e6, L_gen14 t=9.9e5.
**All four have a bin dump at matched t = 1e6**, so the ideal/general x 1e13/1e14 2x2 can be
compared structurally at one time. (L_gen14 was deliberately held ~40 min past the others to
let that dump land -- recreating it would have cost 16 h.) No restart files were configured.

## THE 3 G RUNS -- ALL FOUR CANCELLED 2026-08-18 ~12:40, SUPERSEDED

They used the TRACKED ideal input (x1max = 12.54e9, dfloor 5.44e-12), which is NOT the
configuration the user wanted -- see [[dhj-viper-dt-discrepancy]] for the correct base. The
findings below still stand as measurements of that input, but the live 3 G runs are the
`V_*` set built on `idflr/base_viper.athinput`.

### The cancelled set

All ideal, `deep_hot_jupiter_rt_ideal_xe.athinput`, `problem/bbot=3.0e0`, whole node each.

| job | dir | resistivity | max_eta | floors |
|---|---|---|---|---|
| 190484 | `L_id14_b3` | eos | 1e14 | default |
| 190485 | `L_id13_b3` | eos | 1e13 | default |
| 190486 | `L_id12_b3_p1t10` | eos | 1e12 | pfloor 1, tfloor 4.593e8 (=10 K) |
| 190495 | `L_id12_b3_p1t10_perna` | **perna** | 1e12 | same |

**`tfloor_kelvin` CANNOT be used in the ideal input** -- it converts through `<units>/mu = 1`
while the run's real mu ~ 1.81, so it lands at ~18 K, and it is a FATAL error alongside the
file's existing `tfloor`. Use `mhd/tfloor = 4.593e8`, which is exactly 10 K because the ideal
branch of the resistivity uses `T = p/Rgas/rho` (`resistivity.cpp:176`) -- the same convention
as the dynamics. Confirmed.

### Findings at t <= 1.93e5 (EARLY -- this is the window that gave the 1.77x error)

- **max_eta = 1e13 is the recommendation at 3 G.** 1e13 vs 1e14 is **2.9x** in cost at 3 G
  against only 1.12x at 10 G: the weaker field lifts the hydro dt, giving the cap more room to
  bind. 1e13 -> 1e12 buys ~0% (see below) while capping 9.5% of cells instead of 6.1%.
- **The 1e14 run is measuring its own cap and nothing else**: dt = **1.018 at t = 1e5, 1.5e5
  AND 1.9e5, identical to four digits.** A hard plateau = dt_diff. 3 G buys it only 1.06-1.09x.
- **At 1e13 the weaker field is worth 1.5x in dt and growing** (1.29x at t=5e4 -> 1.56x at
  1.9e5) as the atmosphere magnetizes and the radial Alfven limit takes over from the
  B-independent phi limit. Cost 0.00779 vs 0.01130 -> 3 G is 1.45x cheaper at 1e13.
- **The pfloor=1/tfloor=10 K pair is NOT worth it.** Its early 3.6% edge evaporated: dt ratio
  vs L_id13_b3 went 1.03 (t=5e4) -> 1.11 (1e5) -> 0.99 (1.5e5) -> **0.95 (1.9e5)**, cost now
  identical within 0.8%. Static snapshot analysis had already shown dt unchanged (0.7498 ->
  0.7498) because the limiting cell is Alfven-limited at v_A/c_s = 27. See
  [[dhj-dt-limited-by-alfven-floor]].
- **The floor swap is confirmed in the event log.** Per cycle, `tfloor` firings went 567 -> **0**
  (the 10 K floor is completely inert) while `efloor`/pfloor went 32,756 -> 61,766, i.e. pfloor
  now fires in essentially every cell that also hits dfloor. Exactly the predicted role swap.
- If a less artificial T floor is genuinely wanted, lower `tfloor` ALONE and keep pfloor = 0.1.
  Safe at 1e13 (eta is capped so dt_diff stays ~9.5 s), NOT safe at 1e14 (dt_diff ~0.95 already
  binds).

`perna` skips the electron-fraction table build entirely -- cycle 0 at elapsed 0.19 s vs 4.9-6.6 s
for the `eos` runs, a useful check that it is on the intended branch.

## Why 1.77x was wrong -- the dt advantage is transient

The short runs measured t = 2000..5400, where general had **dt = 15.2 vs ideal 8.3** and so
paid only 1.77x despite being ~3x per cycle. By t > 1e5 both dts have collapsed to ~0.9 and
general's is now marginally **smaller** (0.916 vs 0.959). So per-cycle and per-simulated-second
now coincide: **~3x**, matching the 3.07x/cycle of [[general-eos-table-cost]]. Never quote an
EOS cost ratio from t < 1e5 on this problem.

## max_eta: 1e13 wins for BOTH EOSs -- the old "general wants 1e14" is DEAD

Per-cycle cost is identical between the two caps (0.0189 vs 0.0192; 0.0557 vs 0.0565), so the
entire 8-12% is dt. At 1e14 `dt_diff` has decayed to **~0.95 s**, comparable to dt_hyd, and
clips: neither 1e14 run ever exceeds dt = 0.97 (L_gen14 sits at a 0.9087 plateau for 26
consecutive diag lines), while the 1e13 runs reach 1.32-1.51. At 1e13, dt_diff ~9.5 s never
binds. Yesterday's "dt_diff = 1.9e4, three orders from binding" was an EARLY-TIME reading and
does not hold at late times.

Consequence for `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`: the fix is **`max_eta = 1e13`**,
and then `use_rkg_sts = false` follows for free (at 1e13 diffusion never binds, so STS costs
31% for nothing). 1e14+STS is worse than 1e13 no-STS by estimate (~0.081 vs 0.061 wall s/sim s).
**NOT YET CHANGED -- ask the user.**

Related: [[general-eos-table-cost]], [[dhj-ideal-xe-floor-relaxation]],
[[resistivity-xe-table-for-ideal]], [[session-state-2026-08-17]].
