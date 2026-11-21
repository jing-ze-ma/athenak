# NOTE 2026-11-21jz orion -> caltech: OPAC arms (USER top priority) on the 4 test nodes

Branch caltech-test-opac = 00ad333b7 (af9e8d420 + 154acbabe clean_stag_opac_update + 00ad333b7 clean_face_opac=harmonic; both default
off = v20.7 bitwise). Build with CUDA 12.9 + g++13 (B_opac), tiny BSG smoke gate (pass 1 4.020e-06) first.

Keys for every arm: keys from the rst + rad_m1/clean_conv_max = 1.0e-4 + rad_m1/clean_stag_thinf_held_gate = true. NO clean_pass_damp.
Four variants per restart:  off | o: rad_m1/clean_stag_opac_update = true | h: rad_m1/clean_face_opac = harmonic | oh: both.
300 cycles per arm (as V8), user/hydro hst every 10 cycles, final bins, rst off (delete any).
Waves (one GPU job per node, test nodes only, never production):
  W1: t9l62 (copy of its latest settled rst, as V8a) off / o / h / oh
  W2: T19 (t19r t19l64.00001, N0 1070, as V7a) off / o / h / oh
  W3: t9r off / oh, m35mms off / oh
Report per arm (caltech_results/V210/<arm>/res.txt + the hst files): series every 10 cycles of KE, L_top, L_in, L_bot; CW, FATAL,
solves/stage-pass mean/max, dt range, wall/cycle. Pre-authorised: all three waves.

-- orion (for the USER)
