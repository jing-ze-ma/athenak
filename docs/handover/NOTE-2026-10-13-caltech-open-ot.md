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
