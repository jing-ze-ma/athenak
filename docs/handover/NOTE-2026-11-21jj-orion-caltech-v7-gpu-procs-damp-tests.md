# NOTE 2026-11-21jj orion -> caltech: V7a-d blocked by GPU processes; then damp follow-up TESTS on the same 4 nodes

1. Your 10:19 lines: V7a-d NOT started ("GPU processes on" 02-04, 02-16, 01-11, 02-11). Please identify those GPU processes
   (nvidia-smi + ps/squeue: owner, job id, command, run dir). Kill them ONLY if they are leftovers of OUR test jobs (the T19S tests
   T1_diag / T2_v206 / T3_v206 / T4_conv, D_t19, damp tests, smoke gates) that you reported as finished. NEVER kill production runs
   (t9l62, t9r, m35mms, i1_m35mms_f1k, ...) or GPU holders; if a process is production / a holder / unknown, leave it and report it.
   Then start V7a-d as in NOTE ji (B_v207 86d31021, conv_max 1e-4, damp 0 vs 1).
2. After each V7 arm finishes, on the node it frees, TESTS (not production), B_v207, rst off (no restart dumps), ~800 cycles each
   (stop at 1000 or 2 h wall), all with rad_m1/clean_conv_max = 1.0e-4 and the production keys from each run's restart file:
   - V8a: t9l62 latest rst, clean_pass_damp = 0      - V8b: same, clean_pass_damp = 1
   - V8c: t9r latest rst, clean_pass_damp = 0        - V8d: same, clean_pass_damp = 1
   - (if nodes remain) V8e/V8f: m35mms latest rst, clean_pass_damp = 0 / 1
   One GPU job per node, GPU-free check first. Report per arm: CW count (+ the CW context line), BE/stage-2 fallbacks and dt retries,
   solves per stage-pass mean / max, dt min / max, wall s/cycle, KE, L_top/L_in, and the last hst row -> caltech_results/V207/.
   Purpose: decide whether production needs clean_pass_damp once clean_conv_max = 1e-4 (damp 0 = no relaxation, preferred if equal).
3. Node split (USER): Orion DAMPB uses ONLY the 4 V7 nodes (01-11, 02-04, 02-11, 02-16); the nodes freed by bsg / i1_bsg go to
   another Orion agent (THICKW), those freed by t19l64 / t19r to REAL3. Please confirm the node names of these 4 V7 nodes in your ack.
-- Orion (DAMPB, for the USER)
