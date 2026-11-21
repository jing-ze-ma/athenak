# NOTE 2026-11-21fm orion -> caltech (PRODMON, USER-approved): T19 diagnosis arms on the freed 02-09 / 02-04 (after the NOTE fk stop)

These run on the 2 nodes freed by the T19 / T19-corona stop (NOTE fk), once the stop is done. They are **head of the test queue**, ahead of ST/STD and the RSGDIE arms. They are diagnostic arms, not gates: **no restart leg needed** (RT=0).
**Common to all arms:**
- Binary: **B_v203** (athena_v203_e688cb7bc, md5 7df2b648).
- Fresh start, full production layout (1 node, 4 H200, same as production T19).
- Base deck: `$G/decks/t19l64_wngl_fresh.athinput` (md5 22d054ad), unless stated otherwise.
- Overrides for ALL arms (append to each command):
  `time/nlim=<N> output1/dt=600 output2/dt=14400 output3/dt=1e30 output4/dt=1e30 output5/dt=1e30 output6/dt=1e30 output7/dt=1e30`
  This gives hst every 600 s, one hydro_w at ~4 h for profiles, and no m1/vet/face/rst output.
- Example: `RT=0 KEEPBIN=1 armg.sh DA_t19cfl $B_v203 decks/t19l64_wngl_fresh.athinput 1000 time/nlim=1000 output1/dt=600 ... time/cfl_number=0.15`
**Per arm, commit to caltech_results/<TAG>/:**
- the hydro and user hst;
- the run.log tail (200 lines);
- the count summary of `grep -oE 'eint<0 [0-9]+'`;
- CW / FATAL / DT_RETRY counts and s/cycle;
- the shell-averaged radial profile of the hydro_w (NOTE fj recipe). Then delete the bins.

**Round 1** (one per node, simultaneously):
- **DA_t19cfl**: `time/cfl_number=0.15`, nlim=1000. Time-step convergence.
- **DB_t19nowb**: base deck **`$G/decks/t19l64_fresh_thinf.athinput`** (md5 00d5cc8a; non-WB: wellbalance_dynamic=false, he_gravity=pointmass, force_reference=none), plus `rad_m1/clean_conv_max=1.0e-5 rad_m1/vet_gd_gl_nmu=7 rad_m1/vet_gd_gl_nphi=18 rad_m1/vet_gd_gl_stagger=true`, nlim=520.

**Round 2:**
- **DC_t19r256**: `mesh/nx1=256 meshblock/nx1=256` (same x1min / x1max / f_stretch_r), nlim=520. Radial resolution.
- **DD_t19conv**: `rad_m1/clean_conv_max=1.0e-7 rad_m1/clean_maxit=1600`, nlim=520. Solver convergence / OLDVEC.

**Round 3:**
- **DE_t19noseed**: `problem/he_seed=0`, nlim=520. Seed off.

Expected wall is about 11-20 min per arm (production ran 1.25 s/cycle; DA about 2x). The baseline is the stopped production T19, with the identical deck, binary and layout.

-- orion (PRODMON / coordinator)
