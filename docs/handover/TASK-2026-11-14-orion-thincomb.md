# TASK 2026-11-14 viper -> Orion: thin-switch combine arms (base / min / aniso) on the CPU gate set (USER 11-13 ~23:00 PST)

Context: the RSG stall is a cold opaque cell whose own C1 = <|I-S|>/J ~ 1 (heating front, not streaming) pins the fse face weight to
streaming through max(g_L, g_R). Two tau-free alternatives (user rule: the blend must never use cell tau) are on
**fork/rsgposlin-1113 = 48558c0d** (record the hash):
- `rad_m1/implicit_thin_combine=min`   : face takes min(g_L, g_R) instead of max.
- `rad_m1/implicit_thin_combine=aniso` : each cell's measure = half-range asymmetry |h+ + h-|/(h+ - h-) of its own rays, max-combined.
- no key = base = prod-1118 radiation (bitwise to max).

## Build (new worktree; never touch main; no destructive git in existing worktrees)
`git fetch fork && git worktree add --detach wt_tc fork/rsgposlin-1113 && git -C wt_tc submodule update --init --recursive`;
build_none with exactly the CMake options of your wt_p18/build_none (and build_omp if you run the OMP rows). Record md5s.

## Arms (3 per row/case): base | TC=min | TC=aniso, i.e. the deck/keys as in your idortp1118 q arm (decks_d6 + EXTRA_HG, splitkeys
insertion as you did) plus nothing / `rad_m1/implicit_thin_combine=min` / `rad_m1/implicit_thin_combine=aniso`.
1. The Raven IDORT set (main 41 + rsx 12 + spxx 100; OMP rows optional): crossing beams cb/cbg/cbh, ws66, shadow shdg, hcyl,
   homogeneous sphere, non-LTE nltek, core/curv. Runner as one_q.sh with an arm tag (base/min/ani).
2. Grey atmosphere gate G1 (vs Hopf, pass max|dE/E| <= 2e-2): deck + script in this commit, `docs/handover/task-thincomb-1114/grey/`
   (from blendall_1009, 10-09; add the EXTRA_HG keys and the thin-switch keys of the IDORT decks if the deck lacks them; state what you added).
3. Stiff pulse_cs time ladder dt = T/1..T/16: decks `task-thincomb-1114/pulse/pulse_T_{1,2,4,8,16}.athinput` (prod-1118 radiation keys,
   blend_dtc unset), order from the final m1 dumps with `pulse/order.py <dir>` (as on viper). Viper base reference: E dL1 2.57e-5,
   1.61e-5, 4.50e-6, 1.50e-6, orders 0.68 / 1.84 / 1.58.

## Pass rules
rc=0, 0 FATAL, NON-CONVERGED=0; every metric <= 1.05 x base (and vs FIN as before); ladder orders kept (>= base - 0.1);
G1 <= 2e-2. List every row where min or aniso != base (bitwise) with the size of the change, and where a metric moves > 1.05x.

## Return: NOTE-2026-11-14-orion-viper-thincomb.md on bsg-files-1009, part by part. Times in PST.
