---
name: m1-precond-mg
description: 09-25 merged 091d1422 -- <rad_m1>/implicit_precond = mg halves BiCGStab iterations (box 21.5 -> 9.8) but saves only ~7 %/cycle on the box and nothing on the sp wedge; use for Cartesian boxes only
metadata:
  type: project
---
tests_m1/runs_5m_precond/README.md. Diagnosis: the implicit M1 system is nearly singular for sideways-smooth error
(k = 0 mode 98-99.7 % of the box error after rbgs_fwd; wedge k = 1-2 at the base); a coarse correction is the lever,
more sweeps / ADI lines are not. mg = rbgs_fwd + 2x2 sideways aggregates per block (x1 never coarsened), keys
implicit_mg_levels (code default 2; recommended 3) and implicit_mg_halo (true). GPU: box 1 GPU 56.9 -> 52.8 ms/cycle
(-7 %, -8 % with implicit_precond_float), 2 GPU -6 %; wedge -3.5 % (1 GPU), +2 % (2 GPU). Each iteration costs
+0.5 ms (box) / +0.7 ms (wedge), mostly coarse-level kernel launch overhead. Offline: a GLOBAL V-cycle would give 3/2
iterations (box/wedge) -> the real prize.

**How to apply:** Cartesian box inputs: implicit_precond = mg, implicit_mg_levels = 3 (+ implicit_precond_float);
sp wedge: keep rbgs_fwd. Next step if M1 speed matters again: fuse the coarse-level kernels, then a global coarse solve.

**09-25 continued (user: make implicit VET faster AND scale, Cartesian + sp):** branch m1-mgfuse: phase 1 fuse coarse-level kernels; phase 2 GLOBAL coarse space = per-x1-layer sideways means (one allreduce of nx1 values, tridiagonal x1 solve on every rank) + lowest sideways modes for the wedge, as a two-level preconditioner with mg; scaling 1-8 GPUs. m1-fast3 (quick wins) still running; after it: multi-GPU scaling profile.

**m1-fast3 merged 3e7659fd (09-25):** defaults (inert-conduction skip, RosselandTable bisection) bitwise, box -8 %; implicit_op_team_red OFF by default (round-off; box -12 % with krylov_dev=2 + implicit_eos_cache_check_every=10, wedge 1 GPU -13 %; the slow part was the flat 4-value Kokkos reduction, 134 us vs 36.5 us stencil); lin_tol stays 1e-10. team_red DEFAULT ON since 09-25 (user).
