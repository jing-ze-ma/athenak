# NOTE 2026-11-21kd orion -> caltech: W4 follow-up — split h from o on RSG (o was uncapped)

W4b RSG oh FATAL (19 CW, 166 passes) and W4d HG oh FATAL; W4c HG off itself has 20 CW (dt down to 0.43). W4 ran branch 00ad333b7,
whose o (opac_update) refreshes EVERY pass. Branch caltech-test-opac is now 772ee5ac6 = + rad_m1/clean_stag_opac_passes (default 3,
old prod's opac_lag = pred N = 3). Rebuild B_opac3 (CUDA 12.9 + g++13, smoke gate), then on the next free test nodes, in parallel with W5:
  W4e RSG h   : W4a's setup + rad_m1/clean_face_opac=harmonic only
  W4f RSG oh3 : W4a's setup + rad_m1/clean_face_opac=harmonic + rad_m1/clean_stag_opac_update=true (cap 3 is the default)
150 cycles, rst off, no damping -> V210. Skip further HG arms (its off baseline already stalls). Also post W4c's first CW context (cw_context).

-- orion (for the USER)
