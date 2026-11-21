# NOTE 2026-11-21kh orion -> caltech: W4g / W4h = RSG oh3 with the cs lateral-sum fix (± gas Newton switch)

Two new test branches (CUDA 12.9 + g++13, smoke gate each):
  caltech-test-opaccs      = 58517aae5 (772ee5ac6 + cs fix: the opacity refresh no longer leaves stale lateral sums for xdonor)
  caltech-test-opaccs-gnsw = 3ecd98bc6 (= 58517aae5 + clean_gas_newton_switch key)
Arms, EXACTLY W4f's setup (rsg.01295, 150 cycles, damp 0, W4f keys = harmonic + opac_update, cap 3 default), on two test nodes in parallel:
  W4g: caltech-test-opaccs binary, W4f keys.
  W4h: caltech-test-opaccs-gnsw binary, W4f keys + rad_m1/clean_gas_newton_switch=10.
Report as W4f (CW, FATAL, passes/stage mean/max, dt range, s/cycle, series) -> V210/W4g, W4h. Priority: next on free test nodes.

-- orion (for the USER)
