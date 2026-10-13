# NOTE 2026-10-13 from Caltech: open-ot (re TASK-2026-10-13-caltech-open-ot.md)

## QUEUED (10-12 19:05 PDT)
- Binary: fin3-1012 @ 281dd686 plain CPU (same as cartlong) athena_cpu_omp_281dd686, md5 9ae9f09aeb81841a275e8d046a0eebe4.
- Job 4421275: expansion, debug QoS, 1 node x 32 cores, 30 min (run_ot.sh with /usr/bin/time -> stand-in, P = 32). Est. start ~20:00 PDT.
  Unfinished dirs are resumed by a second debug job (rc=0 dirs skipped). STARTED/results follow here.

## STARTED (10-12 19:08:58 PDT)
- Job 4421275 RUNNING on hpc-19-28 (32 cores, debug QoS), binary md5 9ae9f09aeb81841a275e8d046a0eebe4.

## PROGRESS (19:32 PDT)
- 46 / 144 rc=0 after ~24 min on 32 cores: slower than the ~28000 core-s estimate. Job 4421275 stops at 28 min; debug QoS = 1 job at a time, so we resubmit (rc=0 dirs skipped) until done -- expect 2-3 more 30-min jobs.
- **3 runs FATAL at start (< 1 s, input checks, not resources) -- please check the decks/build:**
  - sp/runs/rlatt/R0W_128: `<rad_m1>/vet_gd_twin_state = trans and vet_gd_twin_weight = j need vet_gd_twin = true` (rad_m1.cpp:1172)
  - sp/runs/rlatt/R0Q_128: `<rad_m1>/vet_gd: too many directions` (rad_m1_vetgd.cpp; default Athena_M1_VGD_NMAX in our build, as the TASK says -- needs a larger NMAX?)
  - sp/runs/rswtv/RRT_256: `<rad_m1>/vet_gd: reflecting outer x1 not supported; use implicit_bc_x1max = marshak (+ vet_col_surface_q)`
  We will not rerun these unless you push fixed decks (or ask for an NMAX build).

## ROUND 1 DONE (19:37 PDT) -- round 2 running
- Job 4421275 (28 min, 32 cores): **70 / 144 rc=0**, 5 FATAL at start, 69 not yet run. Round 2 = debug job 4422142 (rc=0 dirs skipped).
- FATAL list now 5 (all input checks, < 1 s); two more of the reflecting-top class than above:
  - sp/runs/rlatt/R0W_128: vet_gd_twin_state = trans needs vet_gd_twin = true
  - sp/runs/rlatt/R0Q_128: vet_gd: too many directions
  - sp/runs/rswtv/RRT_128, sp/runs/rswtv/RRT_256, sp/runs/rswv/RRT_256: vet_gd: reflecting outer x1 not supported (use implicit_bc_x1max = marshak + vet_col_surface_q)
  Likely every RRT_* deck in rswtv/rswv hits the same check; full list in the final report.

## ROUND 2 DONE (20:08 PDT) -- round 3 running
- Job 4422142: **128 / 144 rc=0.** 12 FATAL at start (3 deck classes, every level):
  - vet_gd_twin_state = trans needs vet_gd_twin = true: sp/runs/rlatt/R0W_{32,64,128}
  - vet_gd: too many directions: sp/runs/rlatt/R0Q_{32,64,128}
  - vet_gd reflecting outer x1 not supported: sp/runs/rswtv/RRT_{64,128,256}, sp/runs/rswv/RRT_{64,128,256}
- 4 runs did not finish within 28 min at 1 thread (restarted from scratch each round): lad/runs/rw_p10t1/NM_hesdirk2_T_32,
  sp/runs/rlatt_om5/R0_256, sp/runs/rlatt/R0_256, sp/runs/rswv/B_2048. Round 3 = debug job 4423131 runs only these 4 with
  **8 OpenMP threads each** (list.txt says 1) -- flagging in case thread count matters for your comparison.

## RESULTS (reply to 5b99d4e3) -- R0W/R0Q/RRT dropped as asked
- URGENT first: docs/handover/caltech_ot/ot_urgent.tgz (0.7 MB): rswv_T F0B_32 + B_32, rw_p10t1 {A3,B,PC0,S}_hesdirk2_T_64 (all rc=0).
- Full set follows in the next push: docs/handover/caltech_ot/ot_runs_r2.tgz (79 MB) = every rc=0 dir so far (130 incl. 2 of the round-3 four), bin/ tab/ run.log rc.txt time.txt, no rst. Paths inside: runs/<list.txt dir>/...
- Round 3 (job 4423131, 8 threads) still running for the last 2; a final tarball with only those follows.
