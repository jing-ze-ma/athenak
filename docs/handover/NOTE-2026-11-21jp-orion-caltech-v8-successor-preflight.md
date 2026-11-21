# NOTE 2026-11-21jp orion -> caltech: V8 = successor-key preflight on the idle non-production nodes (replaces the held V8)

Purpose: before the USER decides the successor restarts (holders end 17:46 PST), check the successor key set on each running star.
Use the idle NON-production nodes freed by bsg / i1_bsg / t19l64 / t19r that NOTE jo's V7a/V7c reruns do not take (and any
V7 node as it frees). One GPU job per node; never a production node.

Binary B_v207 (md5 86d31021). Restart from each run's LATEST rst (copy, do not touch the production run dir), keys from the rst
plus: rad_m1/clean_conv_max = 1.0e-4, rad_m1/clean_stag_thinf_held_gate = true, rad_m1/clean_pass_damp = 1. 300 cycles, rst off.
Order (as nodes free): V8a t9l62, V8b t9r, V8c m35mms, V8d i1_m35mms_f1k.
Report per arm (caltech_results/V208/<arm>/res.txt as for V7): CW, FATAL, solves/stage-pass mean/max, dt range vs the production
run's recent dt, KE and Ltop/Lin vs the production run at the same time, wall per cycle vs production.
Pre-authorised: all four arms; nothing else.

-- orion (for the USER)
