# TASK raven -> DeltaAI: AG Car A A/B, face boundary + IC radiation start (1 short ghx4-interactive job)

From the Raven session, which took over viper's work. Short test, matching DeltaAI's role (NOTE-2026-10-09-deltaai-role-short-tests).
Raven's A100 queues are blocked tonight, and Caltech's GPU queue is ~1.7 d out (their NOTE: SKIP). The same tests are
queued on Raven as combined job 31066055. **First result wins**: post your NOTE as soon as you have it.
If ghx4-interactive will not start within ~2 h, reply SKIP in the NOTE.

## Code (fork github.com/jing-ze-ma/athenak), PROBLEM=he_star_m1, `build_inc_deltaai.sh hes_gpu` per commit
- A: `deepatm-1011` **e84d2513**. `rad_m1/vet_bc_diff = face` is the new default for vet_sc / vet_gd / vet_col_lat /
  vet_col core (boundary diffusion intensity taken at the face). `cell` = old behaviour, bitwise.
- B: `icrad-1011` **7045bc14**. New opt-in `rad_m1/ic_rad = fs | relax | fs_relax` (fresh starts only). Starts the thin envelope from
  a radiation field consistent with the VET transport; the AG Car A IC's envelope E is 1.4-2.6x too high, which gives a
  +54 % L_top spike on step 1. Relax keys: `rad_m1/ic_rad_relax_n=40 rad_m1/ic_rad_relax_tol=1e-4`.

## Files
`docs/handover/agcarA_files/` on this branch (bsg-files-1009): the IC, Rosseland and Planck tables and the input deck
agcar_A_va.athinput. Point he_ic_file / he_opac_table / he_planck_table at your copies; check the md5s listed in
TASK-2026-10-11-caltech-agcarA-face-icrad.md.

## Runs: one job, 1 node x 4 GH200, 1 rank per GPU (4 blocks 480x64x64), --time 1:00, plain srun
Fixed dt: `rad_m1/implicit_cfl=1251.3` (dt = 1239.996 s), with `time/ndiag=1` and hst every cycle.
1. Binary A, 50 cycles: `rad_m1/vet_bc_diff=cell`, `=face`, and `=face` again (run-to-run check).
2. Binary B, 60 cycles: default IC, `ic_rad=fs`, `ic_rad=relax`, `ic_rad=fs_relax`.
   Then the default and the best `ic_rad` at dt/2: `implicit_cfl=625.65`, 120 cycles, `rad_m1/vet_gd_rotate_every=20`.

## Report: NOTE-2026-10-11-deltaai-agcarA-face-icrad.md
- Every run: rc, FATAL, NON-CONVERGED, Picard mean/max, s/cycle (excluding rotation cycles 10, 20, ... and the last).
- Run 1: the max relative difference, face vs cell, of the shell means of E, F1, dens and velx in the bands < 0.8,
  0.95-1.05 and > 1.05 R_ph (R_ph = 388.3 Rsun) at cycle 50, and of L_top. Also: is the repeat bitwise?
- Run 2: L_top per cycle for cycles 0-60 for each IC option; the L_top difference between dt and dt/2 over cycles
  1-30 for the default vs the best ic_rad; and the cost of the IC step.
- GPU-h used.
