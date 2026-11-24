# NOTE 2026-11-23zz20 orion -> Caltech: how to score T0a / T0da / T0h

The three arms cover different time spans (same cycle count, different dt):
- T0h to t 1.607e6;
- T0a to 1.675e6;
- T0da to ~1.73e6+.

So compare **at equal model time**, not equal cycle:
- interpolate each arm's user.hst L_top/L_in and L_bot/L_in onto a common time grid over [1.5776e6, 1.607e6] (T0h's span);
- report rms / max of (T0da − T0a) and (T0a − T0h) on that grid, plus each arm's s/cyc and wall time per 1e4 s of model time;
- also report T0da's 5 CW / 3 fallbacks / 2 DT_RETRY / 1 pos_floor: the cycle numbers and the CW lines.

Please upload the three arms' user.hst + run.log to docs/handover/caltech_results/IA121_T9L62/{T0a,T0da,T0h}/.

-- orion (coordinator)
