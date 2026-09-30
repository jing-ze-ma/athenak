# TASK 2026-09-30 (viper -> orion): He-star wedge + M1 positivity arms on orion CPU

Written on viper 09-30 ~22:50 CEST. Orion is EXTRA CAPACITY for the He presn M1-positivity work; the decisive gates
run on viper GPUs (m1pos_0930, jobs 12046524-30). Read `docs/handover/SESSION-2026-09-30-viper.md` section 6 first.

## Standing rules that apply here
- **ALWAYS smoke-test** (a few cycles, same binary + input + command-line keys, separate dir) before any real job.
- Orion: an orion node has **112 physical cores** (SLURM shows 224 hyperthreads). p.shared: set --cpus-per-task with OMP.
- CPU timings are not cost results (costs are measured on viper GPUs); orion answers pass/fail and physics only.
- Commit only on the branches named here; do not merge; push to the fork only when asked.

## Why
The 3-D He-star wedge (4 Msun presn, pgen he_star_m1) develops a porous Fe-bump outflow from ~3.5 turnovers
(1 turnover = 4700 s) and the implicit M1 solver then fails (E < 0 in fast thin channel cells -> NaN):
64x64 hllc NaN at t 31,626 s (6.7 tt; bad from ~29,200 s), 128x128 NaN at 26,825 s (5.7 tt). Branch
`m1-positivity` adds default-off keys: `implicit_g0_exchange`, `implicit_g0_limit`, `implicit_pos_gas`,
`implicit_pos_floor`, and `implicit_opac_newton_guard_mode` bit 4 (value 6). Open physics question: after the solver
survives, does the outflow settle (steady porous state + wind) or keep draining? (64x64 mass loss had slowed:
-6.6 % at 5.5 tt, -6.9 % at 6.0 tt.)

## Code
Fork github.com:jing-ze-ma/athenak, branch `m1-positivity` (tip at writing **33150622**, on top of `he-presn-m1`
abf547b5). The viper agent may add commits: use the newest tip named in a later viper NOTE if there is one.
CPU build: `cmake -D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Athena_ENABLE_OPENMP=ON` (Freya/orion toolchain
recipe in the memory notes).

## Data (staged on viper, 283 MB): `/viper/ptmp2/jinma/orion_handover_0930/`
- `ad3d_hllc_rst00012_t28200.rst`: 64x64 hllc restart, t = 28,200 s (last clean one; byte-identical to ad3d_hllc rst 00012)
- `inputs/`: he3d_M1.athinput (the base input), r3.sh (the 64x64 job script, KEYS line), arms.sh (the viper arm
  definitions), on.athinput / off.athinput (the restart key files)
- `data/`: ic_he_presn_m1_mlt.txt, rosseland_he_x0.0_z0.02.txt, planck_he_x0.0_z0.02_ferg+tops.txt
The user copies it to orion. **The restart and the input carry viper paths**: override on the command line
`problem/he_ic_file=<orion>/ic_he_presn_m1_mlt.txt problem/he_opac_table=<...>/rosseland_... problem/he_planck_table=<...>/planck_...`
(and fix them in the input for fresh starts).

## Arms (in this order)
A. **Fresh 64x64 from t = 0 with the positivity keys ON** (on.athinput keys + guard_mode 6 + `rad_m1/implicit_vimp=false`
   as in the viper arm onnh64), same KEYS as r3.sh, time/tlim = 56400 (12 turnovers). A fresh start allows small angular
   blocks (meshblock/nx2 = nx3 = 8 or 16; keep ONE block along x1, the sp M1 requirement) so all 112 cores work; a
   restart cannot re-block (the rst has 4 blocks of 184x32x32 -> at most 4 ranks x 28 threads). Estimated ~1 s/cycle
   at perfect scaling (viper column throughput 7e3 zone-cycles/core/s); report the real rate.
B. Control: the same fresh start with the keys OFF, to t >= 33,000 s (should fail near 29,000-31,600 s like viper).
C. Optional, if A passes 7 tt: continue A to 12 tt and report whether the outflow settles.

## Report (NOTE-2026-09-30-orion-he-m1pos.md on the branch you work on, or a note to the user)
Per arm: t reached, NaN yes/no, first-failure time, positivity-limiter counters from the end-of-run M1 summary
(implicit_vimp fallbacks, floor clips, eint <= 0 write-backs, Picard NON-CONVERGED, bicgstab breakdowns), dt(t),
L_top/L_in(t), M_tot(t)/M0 from hepresn.user.hst (cols: 1 t, 2 dt, 6 L_top, 7 L_in, 16 M_tot, 22 Picard), and for
A the mass-loss rate per turnover after 6 tt. Compare with the viper numbers above.
