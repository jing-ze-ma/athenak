# NOTE 2026-11-21kg orion -> caltech: CFL test on t9l62 (USER: can T9 run much faster?) — 3 test nodes in parallel

Binary B_opac3 (caltech-test-opac 772ee5ac6). Restart: a COPY of the latest settled t9l62 production rst (as V8a/W5). Keys: keys from the rst
+ rad_m1/clean_conv_max=1.0e-4 + rad_m1/clean_face_opac=harmonic + rad_m1/clean_stag_opac_update=true (pass cap 3 = default).
NO clean_pass_damp, do NOT add thinf_held_gate (production t9l62 does not name it). rst off.
  C1: time/cfl_number = 0.3   (the production value)
  C2: time/cfl_number = 0.6
  C3: time/cfl_number = 0.9
Run every arm to the SAME end time: t_end = t_rst + 16500 s (~300 cycles at cfl 0.3). Same node type, one GPU job per node.
Report per arm (caltech_results/V212/<arm>/res.txt): wall per SIMULATED second (the key number), s/cycle, cycles, dt range, CW, FATAL,
passes/stage mean/max, and series every 10 cycles of KE_int, L_top, L_in, L_bot; plus the differences C2-C1 and C3-C1 at t_end.
If any arm hits CW/dt_retry, report where and keep going. Priority: next (after GPUCMP's single H200 arm, which can share).

-- orion (for the USER)
