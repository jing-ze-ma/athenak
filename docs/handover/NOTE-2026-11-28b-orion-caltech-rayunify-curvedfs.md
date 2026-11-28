# NOTE 2026-11-28b orion -> Caltech: rayunify, r-face face sums in the sweep on sp/cs (cost gate <= 3 %/cycle on bsg)

From RAYUNIFY (orion). Queue this behind 28a. Idle nodes only; 4 single-GPU arms per node where possible.

The USER takes this change only if it costs <= 3 % per cycle on bsg.

## Binaries
Build all three with CUDA 12.9 + g++13, HOPPER90, star_vet (+ built_in_pgens if cheap). Tiny gate **tiny_bsg_rel6** (pass1 3.720e-11 for REL and CV0; CV may move: post it).
- **REL** = branch `caltech-rayunify-e449d772c` (the vet-next5 e02e092e8 base + the rayunify chain).
- **CV** = branch `caltech-rayunify-8b48dd3ed` = REL + 2 commits: the r-face closure_consistency sums are made in the sweep on sp and cs. The lateral faces are unchanged.
- **CV0** = CV built with `-DVET_SP_FS=0` added to CMAKE_CXX_FLAGS. It must give the same bits as REL.

## Arms
Use the 27zy bsg deck (template of the branch, layout_tile 16, nlim 120, clean_verbose 1, outputs off). S4 = 4 GPUs, meshblock 256 128 128. Cycles 61-120. Order REL, CV, CV0.
- **G1**: s/cycle REL vs CV, the vet_gd fenced line, and builds per cycle. **Pass if CV <= 1.03 x REL.**
- **G2**: kernel timer (libkp_ktime.so) REL vs CV, top 15. Name vet_vgd_shell / _int / _edge, vet_vgd_spfs_sum, cln_gd_fs.
- **G3**: ncu, one launch of vet_vgd_shell_int, REL vs CV: registers, spills, occupancy, duration.
- **G4**: CV0 vs REL: are the CLEAN lines identical cycle by cycle? yes/no.
- **G5**: CV vs REL physics at cycle 120: luminosity (user.hst), Krylov K/T.
- **G6** (optional): RSG star3d_cs (template), 20 cycles from fresh, REL vs CV s/cycle.

Files to caltech_results/RAYUNIFY_CV/<arm>/. Delete dumps; keep no rst.

-- RAYUNIFY (orion)
