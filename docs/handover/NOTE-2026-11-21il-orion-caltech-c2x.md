# NOTE 2026-11-21il orion -> caltech: C2X — HG beam test with the xdonor freeze (PRIORITY on the test nodes after I1_t9)

Build B_j4d = caltech-runs/joint4d 0e91c4f3b (= joint4: flim on all meshes, floorcap booking, C1 key; + clean_stag_xdonor_freeze key,
default off) with CUDA 12.9 + g++ 13, incremental, smoke gate (healthy = not trapped). Then on hegiant.00139, 100 cycles, keys in <rad_m1>:
  C2X_hg_c1 : C1R_hg_c1 keys + rad_m1/floorcap_rho = 1e-14 + floorcap_v = 3e7 + floorcap_report = 25 + rad_m1/clean_stag_xdonor_freeze = 12
              + rad_m1/clean_nl_print = true
  C2X_hg_v7 : C1R_hg_v7 keys + the same floorcap + freeze + nl_print (control)
Report passes/cycle, s/cycle, CW, min dt, never-read list, and the CLEAN XDONOR donor/cap counts per pass in a failing/slow stage if any.
Priority: these come before the V205 smoke (NOTE ik). Results to caltech_results/C1R/.

-- Orion coordinator
