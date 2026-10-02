# NOTE 2026-10-02 (viper): QUEUE the BSG arm-2 production NOW, fresh start, final code later

For TASK-2026-10-01-{caltech,deltaai}-bsg.md and NOTE-2026-10-01-bsg-arm2-in-rt-integration.md.
User decision 10-02: queue the BSG runs on DeltaAI and Caltech now, because the queue wait (DeltaAI estimate ~10 d for a
24 h 1-node job) dominates; the code for the fresh start is being finalised on viper today. Viper also restarts both arms
fresh (duplicates across sites are intended; whichever starts first; they are a cross-machine check).

## Why fresh, and what changes in the code (lands in rt-integration within hours; watch for NOTE-2026-10-02-bsg-code-ready)
The viper arm-2 run stalled from ~4 d: a round-off floor in the gas-eliminated implicit-M1 source row at the inner wall
(K = c dt rho kappa_P ~ 2e8; Picard stuck at resid ~1.2e-7 > tol 1e-8, up to 200 passes, cost x3). Fixed by an exact rewrite,
<rad_m1>/implicit_src_stable (default true for fresh runs). Also new defaults for fresh runs: <rad_m1>/implicit_face_weight
= distance (sp), <problem>/he_bc_hse_flux = face (top-boundary Gamma from the transported face flux), restart guards.
DO NOT run production on bsg-arm2 d0abe21e or rt-integration before the code-ready NOTE.

## Recipe: queue now, run only when ready
1. Submit now a 1-node x 4-GPU, 24 h job (one go: arm 2 at ~0.9-1 s/cycle on GH200 after convection develops is ~15-17 h
   for ~56,000 cycles to tlim 4.96e6 s = 57.4 d) whose script, AT START:
   - reads the binary and input from FIXED paths you fill in later, e.g. $W/bin/athena_bsg and $W/bsg3d_arm2.athinput;
   - exits 0 immediately unless a marker file $W/READY exists (log "not ready, resubmit");
   - otherwise runs a fresh start (no rst in the run dir) or restarts from the newest rst (with the all-output last_time
     reset, see below) and writes STOP on NaN/FATAL.
2. When NOTE-2026-10-02-bsg-code-ready appears: build he_star_m1 at that commit (CUDA), rerun the port gate (column CPU vs GPU,
   numbers in that NOTE) and a 20-cycle 3-D smoke with the exact binary/input/keys, check the log prints
   implicit_src_stable = true, implicit_face_weight = distance, he_bc_hse_flux = face; then put the binary/input in the
   fixed paths and `touch $W/READY`.
3. Input changes vs the bundle bsg3d_arm2.athinput: <output1> dt = 1000.0 (light-curve cadence <= 1158 s, the paper's).
   Nothing else; the new defaults come from the code.
4. If the job starts before READY it exits; just resubmit (it lost only its place in the queue).
5. Restart links (only if a link ends at the wall limit): reset every <outputN>/last_time = floor(t/dt)*dt from the rst header
   (the wall-limit Finalize dump advances last_time by one dt; without the reset each restart skips one output interval).
   viper's helper: bsg_1001/backfill/rst_info.py (prints the keys; t is 232 bytes after <par_end> in the rst).

Report back with NOTE-2026-10-0x-<site>-bsg.md: job id, start time, s/cycle, any FATAL.
