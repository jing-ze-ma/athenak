# NOTE 2026-11-21ks orion -> caltech: t9r 26485 diagnostic — 3 short arms on TEST nodes (not production)

The 26485 stall on B_fwbs is a non-contracting pass cycle (finite, ~1 Krylov it/pass, 200 passes, BE last change 1e-5), not the
26418 NaN divergence. To tell whether it is the xdonor donor/cap flicker or the gas E<->T loop, on three idle TEST nodes in parallel:
from a COPY of t9r's restart (t9l62.00033, the same rst kp used), B_fwbs, the kp keys (incl. fwb_stage on, thinf_held off, cfl 0.9),
+ rad_m1/clean_nl_print=true, run to ~5 cycles past 26485, rst off:
  D_a: as is
  D_b: + rad_m1/clean_stag_xdonor_freeze=20
  D_c: + rad_m1/clean_gas_newton_switch=4 rad_m1/clean_gas_newton_switch_min=6
Post per arm: CW/dt_retry lines, and grep -B40 of the cycle-26485 stages (per-pass changes, CleanStagXCount donor/cap counts, worst-cell
positions if printed) -> caltech_results/V214/<arm>/. Diagnostics only; production t9r stays on cfl 0.6 (NOTE kr).

-- orion (for the USER)
