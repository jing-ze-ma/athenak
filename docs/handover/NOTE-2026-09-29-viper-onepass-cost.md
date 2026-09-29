# NOTE 2026-09-29 (viper): cost of `implicit_one_pass = 8` on top of Newton, He box, MI300A

Gate (c) of `TASK-2026-09-29-viper-onepass-newton.md` only (gates (a)/(b) not asked, user 09-29).

## Setup

- Binary `/viper/ptmp2/jinma/builds/bin/athena_box_gpu72_3058c6bd` (src = rt-integration d26b7364,
  ROCm 7.2, gcc/16, openmpi_gpu/5.0), `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1`, apudev, 1 node.
- Input `docs/handover/bench-2026-09-29-hebox/inputs/hebox_bench.athinput` (84x104x104, fresh start,
  hesdirk2, Newton default), `time/nlim=300`, `time/cfl_number` = 0.3 or 0.9 on the command line.
- Arms (keys edited in input copies): **N** = input as is (`implicit_one_pass = 0`);
  **O30** = `implicit_one_pass = 8` (auto on by default), `time2_one_pass_safety = 30`;
  **O10** = the same with `time2_one_pass_safety = 10`.
- Jobs 12021097 (1 GPU) and 12021098 (2 GPUs), each r1: {cfl 0.3: N O30 O10; cfl 0.9: N O30 O10}, then r2.
  Run tree `/viper/ptmp2/jinma/onepass_cost_0929` (runs/, log_n1.12021097, log_n2.12021098).
- `ana_bench_hebox.py --c0 100 --c1 300`: ms/cyc = median of 10-cycle windows; ms/cyc_w = whole window;
  wall/sim-s = median; Picard/solve over the whole run.

## Results (mean of the 2 repeats; repeats agree within 0.3 % at 1 GPU, 1.4 % at 2 GPUs)

| GPUs | cfl | arm | ms/cyc | ms/cyc_w | wall/sim-s | dt | Picard/solve | one_pass: accepted / measurements / switch-offs / re-probes | vs N (median) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.3 | N   | 49.35 | 49.20 | 0.3061 | 0.16121 | 2.020 | - | - |
| 1 | 0.3 | O30 | 46.35 | 46.56 | 0.2875 | 0.16121 | 1.329 | 410 / 141 / 1 / 0 (end on,off,on) | -6.1 % |
| 1 | 0.3 | O10 | 45.18 | 45.25 | 0.2803 | 0.16121 | 1.179 | 500 / 96 / 0 / 0 (end on,on,on) | -8.4 % |
| 1 | 0.9 | N   | 57.19 | 55.32 | 0.1183 | 0.48339 | 2.030 | - | - |
| 1 | 0.9 | O30 | 57.33 | 56.48 | 0.1186 | 0.48339 | 2.015 | 0 / 84 / 2 / 8 (end on,probe,probe) | +0.2 % (_w +2.1 %) |
| 1 | 0.9 | O10 | 57.12 | 56.33 | 0.1182 | 0.48339 | 2.015 | 0 / 84 / 2 / 8 (end on,probe,probe) | -0.1 % (_w +1.8 %) |
| 2 | 0.3 | N   | 40.06 | 39.98 | 0.2485 | 0.16121 | 2.020 | - | - |
| 2 | 0.3 | O30 | 37.52 | 37.72 | 0.2328 | 0.16121 | 1.329 | 410 / 141 / 1 / 0 | -6.3 % |
| 2 | 0.3 | O10 | 36.64 | 36.71 | 0.2273 | 0.16121 | 1.179 | 500 / 96 / 0 / 0 | -8.5 % |
| 2 | 0.9 | N   | 47.23 | 45.55 | 0.0977 | 0.48339 | 2.030 | - | - |
| 2 | 0.9 | O30 | 47.54 | 46.73 | 0.0983 | 0.48339 | 2.015 | 0 / 84 / 2 / 8 | +0.7 % (_w +2.6 %) |
| 2 | 0.9 | O10 | 47.13 | 46.55 | 0.0975 | 0.48339 | 2.015 | 0 / 84 / 2 / 8 | -0.2 % (_w +2.2 %) |

All 24 runs: NON-CONVERGED = 0, fallbacks = 0, no FATAL/NaN, dt identical across arms.
"accepted" = solves accepted after one pass; "measurements" = contraction measurements; switch-offs /
re-probes from the `implicit_one_pass_auto` line; "end" = state (be, stage 1, stage 2) at cycle 300.

## Verdict

- **cfl 0.3:** the DeltaAI gain reproduces on MI300A: O30 -6.1 % (1 GPU) / -6.3 % (2 GPUs), O10 -8.4 % /
  -8.5 % (DeltaAI -7 % / -10 %). Picard/solve 2.02 -> 1.33 (O30) / 1.18 (O10), the same as DeltaAI.
  At safety 30 the auto logic switches the stage-1 solves off once (still off at cycle 300); at safety 10
  nothing is switched off.
- **cfl 0.9:** nothing is accepted (0 of 84 measurements, both safeties). Auto switches stage 1 and stage 2
  off (2 switch-offs), then re-probes every 64 solves (8 re-probes in 300 cycles: 6 failed, 0 confirmed, 2 still probing at the end). Cost:
  median of the 10-cycle windows +0.2..+0.7 % (O30) / -0.1..-0.2 % (O10), i.e. noise, but the whole-window
  mean is +1.8..+2.6 % (the re-probe cycles). A small loss, not a gain, in the steady-state phase.
- Gain >= 3 % on MI300A only at cfl 0.3; at cfl 0.9 (the decided M1 run setting, hesdirk2 cfl 0.9) there is
  no gain. A default
  `implicit_one_pass = 8` wherever Newton is on would cost ~2 % at cfl 0.9 unless the re-probe period is
  raised; gates (a)/(b) (ke-dt, order) are still to be run before any default change. The user decides.

**USER DECISION 09-29 ~09:25: `implicit_one_pass` stays OPT-IN** (gain < 10 % at cfl 0.3, below the 20 % speed-work rule; nothing accepted at cfl 0.9, +2 % from re-probes). Gates (a)/(b) not run.
