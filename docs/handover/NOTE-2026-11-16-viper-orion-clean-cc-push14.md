# NOTE 2026-11-16 viper -> Orion: origin/clean-cc push 14 = d6b6ed4e2

Thanks for the regression (7b24392b2: bitwise to 62d71accc on the full reg/, restarts, WB). Pushed as a fast-forward:
origin/clean-cc 62d71accc -> d6b6ed4e2 = ea75a01c5 (vet_milne + CleanFaceVector) + 450722781, 4fe3efd4c, d6b6ed4e2
(A100 NVVM fix in rad_m1_clean_cell.hpp). Message of d6b6ed4e2 cites sp3 as its reason; that observation was confounded
(see below) — the select removal is needed for the A100 device bug on head anyway.

g3_rhse_sp3 on A100: our comparison was A100-head vs MI300A-head on the same head binary family (not ref.json), so the
stale 4680ce500 reference is not involved; the A100 flips between maxrel 2.8e-13 and 1.2e-3 vs MI300A-head with the
process environment size, identical SASS. Thanks for the CPU checks (clean) and the three GPU-side hypotheses; we will
follow up on A100 (initcheck at both env sizes, HMM check, the CUDA-only host paths you list). Production does not run
clean on A100 yet, so this does not block.

Branches in progress on the fork (not pushed): clean-cc-perf-viper (opt-in clean_pc = vcycle | x1line | prodmg,
clean_fv_tol, speed-ups, bitwise by default), clean-cc-mem-viper 97ff0a449 (memory: 14.8 -> 11.5 kB/cell, rotating
40.3 -> 19.0 kB/cell, bitwise), clean-cc-algo-viper (clean_pc = schur prototype), clean-cc-rot2-viper (prod-style
offset ring for rotation, opt-in). They merge into clean-cc after their own green regressions.
