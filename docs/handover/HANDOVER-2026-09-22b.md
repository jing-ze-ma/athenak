# Handover 2026-09-22b (viper), written late on 2026-09-21

Continues HANDOVER-2026-09-22.md (which was written at midday on 09-21).  Everything is on
`rt-integration`, local HEAD **86214be9**; the fork is still at **cde8940e** (53 commits
behind; NOTHING has been pushed, the user has not asked for it).  `rt-integration` contains
every local branch (the "unmerged" ones are cherry-picked or ablations); it is 72 commits
behind upstream `origin/main` (last merge 2025-12-19).  Memory notes as of now are in
`docs/handover/claude-memory-2026-09-22b/` (the index was restructured into sub-indexes).

## USER DECISION: M1 is paused, next project = deep hot Jupiter (dhj)

The M1 module is a stepping stone to a VET scheme whose Eddington tensor will come from
short characteristics, refreshed every hydro step.  Analytic chi(f) fixes do not serve that.
Do not start new M1 work (closures, (E,F) Newton, multigrid, more optimisation) unless the
user reopens it.

## 1. What the M1 session established (all committed; RESULTS.txt in each run dir)

| commit | what |
| --- | --- |
| 4e56dc97 | the seeded 2-D He slab is destroyed by the TRANSVERSE implicit transport in the optically thin top (tau < 0.1), not by the frozen reference acceleration (`tests_m1/runs_3b7`); `implicit_trans_limit = lp` (default off) only caps it (`runs_3b8`) |
| 56c97566 | `closure = minerbo \| kershaw` (implicit only).  Levermore, Minerbo and Kershaw ALL fragment the thin top; `closure = eddington` is clean: the instability is the per-step lagged (chi, n) computed from the cell's own F at c dt/dx ~ 7e3 |
| 3ed12a80 | Eddington-consistent column `V3edd` (`bench/m1_stage2/ic/build_ic.py --closure eddington`): seeded slab clean over 2400 s = 5 turnovers (14894 solves, 0 failures, F1top/Fin 1.0000); convection grows at d ln KE_2/dt = 1.2e-3 /s (radiatively damped).  Anderson acceleration of the closure iteration (`implicit_accel`, default off): NEGATIVE result, does not converge chi(f) in the thin top (`runs_3e_newton`) |
| 8d09252b | `implicit_gas_newton`, `implicit_eos_cache` (default off, same converged answer): 2-D slab 1.37x, 1-D column 5.3x; pass count unchanged (~4) |
| de7a52bb | halo copies restricted to the ghost shell, per-pass halo of what changes only: bitwise on 8 arms, halo share 12 % -> 0.6 %; 2-D implicit M1 now ~5x hydro on CPU (`runs_3f_prof/PROFILE.txt`, `runs_3h_halo`) |
| 86214be9 | VET scaffolding `dbg_tensor = frozen \| tilt \| tau` (`runs_3i_tensor`): a tensor NOT taken from the local flux -- frozen chi(z), prescribed off-diagonals, or rebuilt every step from the column optical depth -- is clean in every arm.  The thin-top instability does not carry over to VET |

Closure survey: `docs/dev/rad_m1_closure_survey.md` (exact Hopf surface chi 0.410 at f 0.577;
every chi(f) closure overshoots; user: do not invent closures).
Upstream comparison (their GR well-balanced scheme and PrimitiveSolver EOS vs ours):
`docs/dev/upstream_wb_eos_comparison.md`.

Open M1 items, parked: the ~1 v_MLT coherent 1-D drift of the He column (also in 1-D); GPU and
MPI runs of the 3-D solver; a tensor with grid-scale noise; 3-D He box.

## 2. dhj: where it stands (from the notes + run dirs; nothing new run yet)

`bench/cs_mhd_prod3` (cs MHD, general EOS, ck RT, semi-implicit source, no WB, 2x2 blocks per
panel) ran CLEAN to rot 283 (stopped by a STOP file 09-14).  cs vs sp agree at rot 20/50;
seams/vertices were found innocent.  Open: (1) `src/bvals/bvals_fc.cpp` seam halo with > 1
block per panel (the cc fix 94c7165d is in HEAD; an agent was working on the fc path at the
time of writing: `tests_seam_fc/`, build dir `build_cs_fc`, UNCOMMITTED, check its README
and diff before trusting it); (2) unexplained non-reproducible deaths at rot 11-12 of the WB
and FOFC arms; (3) `dt_min = 1e-2` too tight; (4) lessons from the massive-star work to
check on dhj: r^2 geometry of the column RT (8bca3dfa) in the dhj pgen, the split
ck / conduction handover at 10 bar (half-CFL test), general-EOS MHD restart caches (`wder`),
`wb_cache_every = 1` if WB is used, the implicit transverse conduction on cs
(`docs/dev/cs_implicit_transverse.md`).

## 3. Housekeeping

09-21: the session crashed once (per-user 108 GB cgroup cap on the login node); 47 stale
build dirs, 13 core files and logs were deleted (inode quota).  Build dirs now: `build_cpu_box`
(box_convection), `build_cpu_m1`, `build_cs_fc` (agent).  The `origin` remote URL still embeds a
GitHub token (user action).
