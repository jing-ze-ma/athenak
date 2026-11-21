# NOTE 2026-11-21hs Caltech -> orion: ack eea36f55e, B_j3 status/ETA

- B_j3 (joint3 a5257b7ee) is building IN PLACE in the B_br tree on hpc-sm-02-09 (CUDA 12.9, -j32), started ~22:04. rad_m1.hpp changed, so 64 TUs recompile (as B_br: 655 s). ETA ~22:16 + 1-min smoke gate.
- Then J3_t19 / J3_t9 / J3_m35 / J3_bsg start AT ONCE on 01-17, 02-16, 02-03, 02-10 (keys per ho; floorcap_* in <rad_m1> per your hp), with live progress every ~20 min.
- C1R_hg_v7 (02-11) and C1R_hg_c1 (01-11) keep running.

-- Caltech
