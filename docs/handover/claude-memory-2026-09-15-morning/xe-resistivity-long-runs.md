---
name: xe-resistivity-long-runs
description: "Ideal-gas + EOS-x_e resistive hot Jupiter: the apudev max_eta tests and the 24 h apu1 runs. THOSE RUNS CARRY THE STELLAR-HEATING BUG fixed in b4e0953c -- they absorbed ~76% of the insolation. Do not restart them with a post-fix binary."
metadata: 
  node_type: memory
  type: project
  originSessionId: be8c40b1-14bc-4240-8f72-1fdd200e6113
  modified: 2026-08-19T11:38:38.957Z
---

## THE COMPLETED RUNS CARRY A KNOWN BUG (flagged 2026-08-22)

Jobs 10971428 (`xe_b3_e13`) and 10971429 (`xe_b10_e13`) COMPLETED normally at
2026-08-21T23:55 after 23:50 of wall. Nothing was cancelled; they simply hit their limit.

**They ran the pre-`b4e0953c` binary throughout, so they absorbed only about 76 % of the
incident stellar flux** -- see [[exofms-cross-validation]] for the deposition bug. Treat
their thermal structure accordingly.

**Do NOT restart them from their restart files with a post-fix binary.** The insolation
would jump about 31 % at the restart, giving a discontinuity in the forcing partway through
the run, which is worse than either the buggy or the correct run alone. If these results are
needed, rerun from scratch with a current binary.

**Launched 2026-08-19:** two 24 h `apu1` runs of the deep hot Jupiter with `eos = ideal` +
`ohmic_resistivity = eos`, `max_eta = 1e13`, in `/viper/u2/jinma/ATHENAK/bench/xe_long/{b3_e13,b10_e13}`
(jobs 10955000 / 10955001, `bbot` = 3 G and 10 G). Physics and floors are
`run/2500_test/ohm_128/deep_hot_jupiter.athinput`'s, unchanged. Each directory holds its OWN
snapshot of `athena` (commit `4cc3690c`) plus `BUILD_COMMIT.txt`, and a `resubmit_viper.sh`
that continues from the newest `bin/*.rst` after the 24 h timeout.

**Continued 2026-08-21 00:05** as jobs 10971428 (b3_e13) / 10971429 (b10_e13), each a fresh
24 h `apu1` slot picking up the newest `rst/*.rst` via `resubmit_viper.sh`. The first 24 h
slot ended clean at the wall on 2026-08-20 17:58 (exit 0, 23:50 elapsed), reaching t = 2.53e7
(b3) and t = 8.27e6 (b10). Both still run the `c39b794c` binary and the 32-block decomposition
of [[meshblock-decomposition-gpu]]; resubmit the same way each time they time out.

**Production partition is `apu1`, not `apu`** — that is where the user's own `sohm_dhj2500`
jobs ran (10884546, 10860275). It is the MaxNodes=1 partition, which suits `ntasks=1
--gres=gpu:1`. `apudev` (15 min, 2 nodes) is for tests. See [[viper-hip-build-recipe]].

**Measured on apudev, 12.5 min each, one MI300A** (jobs 10954628-31), 64x64x128, ideal gas +
EOS x_e:

| run | bbot | max_eta | t reached | mean dt | wall/sim-s | zone-cycles/s |
|---|---|---|---|---|---|---|
| b3_e13 | 3 G | 1e13 | 1.443e5 | 4.08 | 5.16e-3 | 2.49e7 |
| b3_e14 | 3 G | 1e14 | 4.512e4 | 1.27 | 1.65e-2 | 2.51e7 |
| b10_e13 | 10 G | 1e13 | 1.125e5 | 3.19 | 6.62e-3 | 2.49e7 |
| b10_e14 | 10 G | 1e14 | 4.583e4 | 1.28 | 1.62e-2 | 2.52e7 |

All four ran clean. `max_eta = 1e13` is 3.2x faster at 3 G and 2.5x at 10 G at IDENTICAL cost
per cycle — the whole difference is timestep, and at matched t = 4.40e4 every global agrees to
<= 3.4e-4 (3 G) / 2.0e-4 (10 G) between the two caps. At 1e13 dt decays smoothly with the flow;
at 1e14 it locks to 1.582031 then steps to 0.78, which is the diffusive limit binding.

**The `tfloor` worry did NOT materialise.** `docs/ideal_gas_resistive.md` recommends
`tfloor = 9.186e9` (200 K) and predicts a stall or NaN without it, but with the user's
`dfloor = 7.26e-12` (1.33x the doc's) the temperature floor fired EXACTLY ZERO times in all
four runs, as did `vceil` and C2P failures. dfloor fires at 8-11% of C2P calls and efloor at
8-9%; mass and total energy drift up ~0.15% / ~0.26% per 1.4e5 s from the density floor.

**Gotcha that cost a submission round:** do not copy a trailing `</output>` into an athinput.
That tag is the Read tool's result wrapper, not file content; the parser reads it as a block
named `/output` and every rank aborts at `parameter_input.cpp:121`.
