# NOTE 2026-11-21kr orion -> caltech: t9r still fails at cfl 0.9 with fwb_stage — apply the kp fallback (cfl 0.6) now

t9r log.3 cycle 26485 (~85 cycles after the kp restart): stage 2/3 solves hit 200 passes with failed linear solves (it -2), BE not
converged, dt_retry x2 - finite values this time (no NaN), but within the first 100 cycles, so per NOTE kp rule 3: restart t9r at
cfl 0.6 with the SAME kp keys (fwb_stage on, thinf_held off) from the same backed-up rst. If it recurs at 0.6, flim at 0.6.
Please post: t9r log.3 (gz) to caltech_results/PROD1121b/t9r/k9b_log.gz, and the first 50 cycles of m35mms / i1_m35mms_f1k / t9l62
on the kp settings (are they clean?).

-- orion (for the USER)
