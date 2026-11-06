# TASK 2026-11-06 viper -> Orion (USER): run the 4 x 4 Athena++ column for Caltech (TASK-2026-11-06p)

The USER asks you to run this. Raven's `small` queue is full. **Don't take nodes from your 3-D run;** this needs only 32 cores on one node.

## What it is
`his1d_g` on Caltech's `Lr` column shape: pi/3 wedge, 4 x 4 lateral cells, his 256 radial faces, his 3-D angles, no seed. Caltech wants to know if the 1-D line vs 3-D mesh and angles explains their v_r and KE differences. Full request: `docs/handover/TASK-2026-11-06p-caltech-viper-4x4-column.md`.

## Files, in `docs/handover/athenapp-bsg/repro3d/col4x4/`
- `athinput`: tested deck, LF line endings. `problem_id` writes to `output/`, so `mkdir output` in the run dir. tlim 373.5 = 25 d.
- `bsg_star_col4x4.cpp`: his pgen with the 4-line output edit, plus the two changes in `col4x4.diff`:
  - Radial faces come from a user mesh generator, enrolled when `x1rat < 0`.
  - `problem/seed_amp` is read from the deck; default is his 0.05, this deck sets 0.
- `link.sh`: the Raven job script, for reference only.

## Steps
1. Build a separate binary from your e47fabb7 tree, using this pgen. Don't overwrite the binary your 3-D run uses.
2. Link his tables and converted profile into a new run dir, as in your 3-D run.
3. Run 32 ranks, one rank per block (32 blocks). Estimated ~1.5 h to 25 d on Raven-speed cores, ~30 min to 8 d.

## Raven smoke test, 20 cycles on 4 ranks
No FATAL. dt 1.8332e-3 code = 10.6 s, first cell x1v 20.157, 2.9e4 zone-cycles/s. Please check that you get the same dt and first cell. The Warning about cell sizes differing by > 10 % is not expected except possibly at 35 and 55 Rsun, where spacing changes by 2x on purpose.

## Report to Caltech, on this branch
**First post a one-line NOTE when the job starts.** Job 31980425 is still pending on Raven as a duplicate, writing nothing; I will cancel it when I see your NOTE. If Raven starts first I will say so here, and you can stop.

Then, to at least 8 d, and to 25 d if it stays cheap, the columns in TASK-2026-11-06p:
- `*_ke.txt`: KE at hst cadence, with KE_r and KE_lat separately.
- `*_wall_profiles.txt`: shell-mean first 24 cells, 0.5 d steps to 5 d, then 1 d to 10 d: r, rho, Pgas, Er(lab), T, vr, Fr1(lab), Fr01, N^2.
- dt at those times, Jacobi stats, and whether lateral motion appears above round-off.
- Against `his1d_g`: KE at 2-4 d and v_r at 22.66 Rsun at 1 d. `his1d_g_ke.txt` and `his1d_g_wall_profiles.txt` are in `docs/handover/athenapp-bsg/cflrad_1106/`.

One difference to state in the NOTE: `his1d_g` got this grid from SMR with prolongation at level edges, and this run uses explicit faces with no level boundaries. The wall region (r < 35) is uniform in both.
