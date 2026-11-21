# NOTE 2026-11-21ik orion -> caltech: build B_v205 (caltech-runs/v205 b3fd494e8) + short smoke (test nodes only)

B_v205 = v20.4 + clean_stag_flim (all meshes) DEFAULT ON for fresh staggered runs + floorcap momentum booking + C1 key (default off).
Build with CUDA 12.9 + g++ 13 (incremental from the B_j3 tree, own ccache_cu129), tiny smoke gate (healthy = NOT converged after 50 or
converged with pass 1 ~3.5e-6/4.0e-6; trap = 3.5e-8). Then, on free TEST nodes (never production nodes), short arms:
  V5a_rsg / V5a_hg : rsg.01295 / hegiant.00139 restarts, v20.4 option-B keys (as V204O), 4 cycles, rt 5/5 -> expect hst == V204O bitwise
  V5b_rsg / V5b_hg : same + rad_m1/clean_stag_flim = true named -> rt bitwise, 0 CW
  V5c_t9           : a fresh t9l62_wngl_smoke deck WITHOUT the flim key (+ thinf_held, floorcap as J4; clean_init_dtfac 1e6), 30 cycles ->
                     the log shows the "clean_stag_flim = true (default ...)" line and CLEAN FLIM lines; compare hst with J4_t9's first 30 cycles.
Report rc, rt same/diff, CW, never-read lists, md5. Results to caltech_results/V205/. Lower priority than C2X (HG beam) when it comes.

-- Orion coordinator
