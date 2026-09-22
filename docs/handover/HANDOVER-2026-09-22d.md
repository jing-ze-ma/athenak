# Handover 2026-09-22d (viper), written 20:15 -- END OF A 27-HOUR SESSION

Continues HANDOVER-2026-09-22c.md.  Branch `rt-integration`, local HEAD = the commit carrying this
file; pushed to the FORK (jing-ze-ma/athenak) at the same commit; `origin` untouched.  Memory as of
now: `docs/handover/claude-memory-2026-09-22d/` -- **read `master-todo-2026-09-22.md` first**; it is
the running log of the whole day with every decision.

## STATE: prod4 LAUNCHED 2026-09-22 ~20:45 (bench/cs_mhd_prod4, jobs 11941995-8, md5 f95130b2, from scratch); see memory cs-mhd-prod4-run.md
`inputs/production/deep_hot_jupiter_cs_prod4.athinput` (DRAFT header still says so; remove the
DRAFT banner when launching).  FROM SCRATCH (user's choice), 2 GPUs / 1 node (measured the efficient
layout), `apu`, chained as bench/cs_mhd_prod3/{submit.sh,chain.sh}.  Binary: HIP build of HEAD from a
clean `git archive` snapshot (recipe bench/prof_0922/README.md; check `strings` for ck_sweep_form and
hlld_bx_zero_tol).  Settings and why (all in the file's header): ck_spherical + ck_beam_sph +
ck_sweep_form = tm (exact spherical sweep), polytropic WB with wb_cache_every 0 and
cs_wellbalanced_src, rad_angular false, dt_min 1e-3, max_eta 5e12 (Ohmic dt cap), pfloor 1e-5
barye (correct top-cell emission), hlld_bx_zero_tol default 1e-4 (the original; 1e-8 dissipated
the field 6x).

## What the session established (commits since 22c, newest last; RESULTS in each tests_* dir)
restart caches bitwise for MHD and cs (0eb7e2c4, 544cb36b) | fc seam halo has no defect (2268e633)
| cc corner slots resampled (ec0a9724) | ck_spherical + corrected beam (3c5d4846, 2db7095c) |
ck_beam_sph pseudo-spherical beam, exact photon budget (13a7cf30, 916dc953) | hipcc build fix
(b4c6a2f2) | HEAD hybrid repaired + ck_implicit (ed188f69) | rad_angular switch (213095e4) |
angular momentum cs = sp (131eabe8) | per-region gate (76359975) | switch inventory (dd317c25) |
GPU A/B + isolation (8e41808b, e31c029f) | CLEAN-UP -4779 src lines (47eb060e, 912ef43c, ac42334d)
| default flips (0e11d533) | HLLD tol parameter (a544a761, 0aaf7f39) | ck_implicit phases 2-4
(4a8daf88, 90c33ce5, 5f348164): 4.2x on GPU, parked | vertex fills cc (212857c1) | Mignone radial
reconstruction (145edc00) | vertex-cell CONSISTENT, order-1.6 retracted (446b7b30) | field-loop
growth = truncation (1a83c06e) | Ohmic cap diagnosed (35498e08) | top slab = correct emission,
pfloor (90c33ce5, 9e5b697f) | probe-free tm/sd sweeps (5f348164) | tm default (0172f3bc) |
MPI gate passed (0f1baed6) | FOP kernel tag, 1.30x regression fixed (cfe92ac7) | profile + cold-
start split-lag gate passed + cache arms (4381ece1).

## Cost, GPU, rot-283 restart, 300 cycles (cycles/s): old production binary 12.84; HEAD defaults
15.26 (tm + FOP tag).  RT 82 % of GPU time, the ck sweep 72 %.  Next levers: sweep rank imbalance
20 % (new with tm), template tags IMP/LEG (probe build halves the kernel), occupancy (1.16
waves/SIMD).  ck_implicit: 4.2x, parked.

## Open small items
eos_efloor nonzero from cold start with pfloor 1e-5 (harmless; check first prod dumps); per-rank
restart files abort on read (build_tree.cpp:429); grey picket-fence rt_split not bitwise vs serial
(7 % KE) and nx1 > 72 fatals; the 36 UNSURE switches; prod3 rst archive 46 GB (retention rule).

## Standing rules added today
cost verdicts on GPU only; Sonnet only for really simple tasks; default-off switches inside the
sweep kernel must be template tags; M1 paused; no pushes to origin.
