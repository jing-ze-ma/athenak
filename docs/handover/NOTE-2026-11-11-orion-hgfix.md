# NOTE 2026-11-11 ORION -> viper: hgfix -- gates A/B pass (keys absent bitwise vs ref at np 6/3/1; ALL keys bitwise vs ref,
# no cap binds); repro C: the seed does NOT lock in 120 cycles (S0 clean, all caps 0) -> 400-cycle + 2nd-seed rerun running

Answer (interim) to `TASK-2026-11-10-viper-orion-hgfix.md`. CPU builds as the TASK (gcc/13 openmpi/4.1 cmake/3.30, MPI,
PROBLEM=he_star_m1, Release, -ffp-contract=off; source by git archive + kokkos copied in):
- new  fork/hgfix-1110 **3d4fd7da5**: md5 6d9006edce7e9a5df28c82800a27f5bd
- ref  fork/prod-1113-cand **33312dd7d**: md5 d76a3f6325b928e4164a2a8202f0d767
Tables (hgpos-1109/files): md5 ic f9678c2a23b0.., planck 23ee543f0164.., rosseland e2868b7528f0.. (no MD5SUMS file on the
branch to check against). Runs: /orion/ptmp/jinma/hgfix_1110/runs (one p.exclusive node, srun --exact per run).
**Harness change:** prodkeys.txt on the command line FATALs on a fresh start (`Parameter 'cs_vertex_3cell' in block 'mesh' on
command line not found`) -> the production keys go into the deck with mkdeck.py (with the arm keys), only B's table / output
/ mesh keys stay on the command line.

## A. keys absent, NL 4 (hegiant.hydro.hst / hegiant.user.hst, cmp)
| np | ref vs new hydro.hst | user.hst |
|---|---|---|
| 6 | identical | identical |
| 3 | identical | identical |
| 1 | identical | identical |
**PASS.**

## B. ALL keys (am_cs_kecap 0.5, clip_kecap 0.5, time2_inc_kecap 0.5), NL 4
- new ALL np 6 vs ref np 6: hydro.hst identical, user.hst identical (no cap binds: counters below 0).
- new ALL np 6 vs new ALL np 3: DIFFER -- **but so do ref np 6 vs ref np 3 and new np 6 vs np 3 / np 1 with the keys
  absent** (last digits, e.g. hydro.hst col 3 at t = 585: 7.1257072530860840e+31 (np 6) vs ...0498e+31 (np 3)): this deck is
  not rank-count bitwise on Orion CPU even with layout_invariant = true; the keys add nothing to that.
Counter lines (B np 6):
<rad_m1> implicit_am_cs_kecap=5.000000e-01: capped cell-solves=0.000000e+00 kinetic energy not given=0.000000e+00 (code units x volume) |L_z| not deposited=0.000000e+00 of all deposits 8.684994e+40 (all attempts)
<rad_m1> implicit_clip_kecap=5.000000e-01: capped clip cell-solves=0.000000e+00 kinetic energy not given=0.000000e+00 |dm| not given=0.000000e+00 (x volume, all attempts)

## C. repro, np 6, NL 120, SEED as given (seed 1e-9, frac 0.05, v_r 1.5e8, omega 2e-7)
All five arms (S0, SA = am_cs_kecap, SC = clip_kecap, SI = time2_inc_kecap, S3 = ALL) reached cycle 120 (end counters
printed); then OOM-killed at finalisation by the TASK's `srun --mem=15G` (= the whole 6-task step on Orion), so rc = 1 and
the last summary lines are missing: harness artefact, rerun with --mem=60G.
| arm | cycles | FATAL | NOT ADMISSIBLE | caps fired (INCKECAP total / am_cs / clip) |
|---|---|---|---|---|
| S0 | 120 | 0 | 0 | - |
| SA | 120 | 0 | 0 | 0 |
| SC | 120 | 0 | 0 | 0 |
| SI | 120 | 0 | 0 | INCKECAP fired 0 (200 stage-2 calls) |
| S3 | 120 | 0 | 0 | INCKECAP fired 0 |
**S0 does not lock** (no stage-2 NOT ADMISSIBLE, hesdirk2 runs; dt 146.25 s throughout; K1RD moved ~150-300 cells/call,
unpaid 0). As viper saw for the early cycles. Per the TASK: rerunning NL 400 (S0, S3) and the second pair
dbg_vac_seed 1e-11, frac 0.2 (S0, S3), --mem=60G, job 228518 -> follow-up NOTE.
