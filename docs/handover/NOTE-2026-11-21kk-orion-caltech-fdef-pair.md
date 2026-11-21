# NOTE 2026-11-21kk orion -> caltech: FDEF pair — flim vs fwb+fcone on t9l62 (on the NEW production settings)

Branch caltech-test-fdef = 03b027502 (integration head 7184d4997 + fwb counter fix df788f2e9 + fresh-run defaults fwb/fcone on, flim off;
restarts keep their recorded values). Build B_fdef (CUDA 12.9 + g++13, smoke gate).
From a COPY of the current t9l62 production rst, 150 cycles, rst off, two test nodes in parallel. BOTH arms name the new production
settings: rad_m1/clean_conv_max=1.0e-4 rad_m1/clean_face_opac=harmonic rad_m1/clean_stag_opac_update=true time/cfl_number=0.9
(no damp, no gate):
  F_A: + rad_m1/clean_stag_flim=true rad_m1/clean_stag_fwb=false rad_m1/clean_stag_fcone=false
  F_B: + rad_m1/clean_stag_flim=false rad_m1/clean_stag_fwb=true rad_m1/clean_stag_fcone=true
Report (V213/<arm>/res.txt): CW, FATAL, passes/stage, dt range, s/cycle, wall per simulated second, max cell |F|/cE and cells > 1 (the
'CLEAN' realisability lines if printed), E-floor raises, series of KE, L_top, L_in, L_bot. Priority: after W4g/W4h.

-- orion (for the USER)
