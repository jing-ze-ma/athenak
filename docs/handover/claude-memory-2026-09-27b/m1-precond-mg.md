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

**m1-mgfuse phase 1 FAILED 09-25 (docs merged):** fusing restriction into the coarse line sweeps was exact but SLOWER (9.8 -> 29.6 us per coarse sweep); GPU idle before mg launches only 0.16-0.27 ms/cycle -> launch overhead is NOT the cost; the cost is ~100 us of small slow coarse kernels + one extra halo per application. Box per-iteration 1.12 ms (mg) vs 0.69 (rbgs_fwd). Phase 2 (global band-mean coarse space) written up in tests_m1/runs_5n_mgfuse/README.md (+ untested draft patch); PARKED.

**09-25 (user): M1/VET speed + scaling RESUMED:** branch m1-coarse2 (/viper/ptmp2/jinma/wt_coarse2, runs tests_m1/runs_5p_coarse2): part 1 = scaling baseline 1-8 GPUs box vet_sc + wedge vet_col (limiters); part 2 = global two-level coarse correction (phase 2, un-parked).

**m1-sp-order2 merged 09-25:** sp wedge was already 2nd order in the transport; first order were the Marshak end face and the vet_col formal solution. New keys (default off): implicit_marshak_face = linear, vet_col_order2 = true -> Eddington atmosphere order 2, free streaming exact, vet_col tensor 2.1-2.25, T-S4 L1 3.5-12x better (flattens ~5e-6); still first order: vet_col radiating pulse (f_K clamp 1/3), reflecting outer x1 under vet_col, tensor time lag. Cost +3.9 % on the He wedge.

**09-25 (user): sp wedge implicit VET must be 2nd order in space AND time by default.** Branch m1-sp-order2b (/viper/ptmp2/jinma/wt_sporder2b, tests_m1/runs_5q_sporder2b): make implicit_marshak_face=linear + vet_col_order2 default on sp (old restarts keep old behaviour), per-stage vet_col rebuild or 2nd-order extrapolation, f_K clamp fix, reflecting outer x1 under vet_col, combined space-time convergence incl. stiff hesdirk2 order reduction.

**m1-coarse2 merged 09-25:** implicit_precond = mg_gc (global band-mean coarse space + mg). WEDGE: use mg_gc with implicit_mg_levels 1 -> 1.1 it/solve, -18 %/-22 % ms/cycle (1/2 GPU), eff 0.82. BOX: keep mg levels 3 (mg_gc ties; residual k=1-2 modes not captured by band means; full global V-cycle would reach 3 it). 4-8 GPU scaling scripts ready in /viper/ptmp2/jinma/coarse2_0925 (not run).

**m1-sp-order2b interim 09-25 18:00 (WIP 481e6528):** new sp defaults (old restarts keep old): implicit_marshak_face=linear (fixed-tensor closures only), vet_col_order2, vet_col_fk_min=0 (1/3 clamp caused the pulse kink), vet_col_reflect_top, vet_col_surface_face, time2_vet_col=predict (1 build at U^n+dt K1; rebuild same order +18 %), time2_vstage (vimp: stage-velocity F + gas work; two hidden O(dt) lags). Space-time orders now ~2.0 everywhere (pulse 0.8->2.07, moving slab F 1.0->1.99, lateral 2.0); stiff hesdirk2 2.00, no order reduction up to rho kappa c dt 1e6. Gates pass; GPU He wedge +6 %. Remaining: mg_gc sp default + README (~1.5-2 h).
User 09-25: preconditioner defaults stay geometry-specific: sp wedge = mg_gc (levels 1, default via m1-sp-order2b), Cartesian box = mg (levels 3) as the recommended input setting since it is faster there; no unification needed.

**m1-fast4 STARTED 09-25 ~18:05 (user):** implicit VET radiation well below hydro cost (target <= 0.3-0.5x) on box (vet_sc, mg 3) and wedge (vet_col, mg_gc 1), scaling; profile first, then Picard passes, setup fusion/caching, tensor cadence with time2 predictor, multi-rate radiation, Krylov launch/reduction cuts; lin_tol stays 1e-10. Branch m1-fast4, /viper/ptmp2/jinma/wt_fast4, tests_m1/runs_5r_fast4.

m1-fast4 interim 09-25 18:40 (profile): box vet_sc mg3 1 GPU 47.5 ms (rad 33.0, rad/hydro 2.28; Krylov 13.7, tensor 4.2, setup 4.0, end 3.3, bvals+C2P 2.5); wedge vet_col mg_gc1 31.4 ms (rad 24.5, rad/hydro 3.41; VetColBuild 7.2-7.7, per-solve setup/post kernels ~9, Krylov only 5.3). Radiation = per-solve FIXED cost, not iterations. Multi-rate implicit_mr_every=k (Strang-placed SDIRK2 every k steps): k=2 same error, k=4 2.3x, k=8 4.9x; adds first-order damping in mixed regimes -> opt-in only. Next: GPU MR timing, tensor cadence, fuse per-solve setup kernels.
m1-sp-order2b MERGED 09-25 ~19:00; He wedge 2 GPU 31.1 -> 26.0 ms/cycle with the new defaults incl. mg_gc.
m1-fast4 interim 2 (09-25 ~19:30): MR k=4 gives box 46.7 -> 29.2 ms (rad/hydro 2.28 -> ~1.0), wedge 33.4 -> 18.1 (3.69 -> ~1.5) BUT MR is first order (sp pulse 1.04-1.4, radwave damping) -> opt-in only, not for sp; eos_check_every=10 box -3.3 %; vet_sc_every=2 same accuracy (4 = 45x worse), GPU build fix being timed; no single 4x lever with 2nd order kept: per-solve setup kernels (scratch spills 176-830 B) + tensor build dominate, 0.5-1.5 ms each.
m1-fast4 interim 3 (76241f53, bitwise): Kokkos HIP RangePolicy = 1024 threads/block caps kernels at 128 VGPRs -> spills; par_for_lb LaunchBounds<256,1> + always_inline + flat reductions + fused opacity lookup: m1_impl_eck 742->115 us, m1_opacity 538->347 (wedge 249->69), i1 371->243; cycle -5..-6 % box and wedge. vet_col team kernel 9.3 ms (28 % of wedge) not limited by access pattern. vet_sc_every=2 slower on GPU (spoils predictor) -> not recommended. Useful general lesson: LaunchBounds for heavy kernels on MI300A.
m1-fast4 MERGED 09-25 ~20:30: target <=0.5x NOT reached with accuracy kept (box 2.08, wedge 3.47 rad/hydro); left: m1_vcol_team 9.3 ms (unknown limiter, needs HW counters), box Krylov 13.7 ms, stb/src/vsf scratch, general-EOS C2P. 4-8 GPU scripts in runs_5r_fast4.

**09-25 ~20:40 (user): continue VET acceleration, separately:** m1-fast5-sp (wt_fast5sp, runs_5s: vcol_team kernel via HW counters, sp setup kernels) and m1-fast5-box (wt_fast5box, runs_5t: global low-Fourier-mode coarse space / cross-block MG, Krylov warm start, fused BiCGStab, SC sweep). No accuracy sacrifice.
m1-fast5-sp interim 09-25 ~21:20: m1_vcol_team is LATENCY-bound on scattered table reads (waves wait 83 %, VALU active 8.8 %, occupancy 36 %, 5 L2 requests per load: seg_/gb_(r,l) strided by n1 across ray lanes; 2 barriers per shell). Fixes (bitwise): v1 chunk of 16 shells per barrier (bfd81013); v2 transposed (shell,ray) tables +161 kB. Other sp kernels spill (tsolve 324 B, src 828 B, vcp 392 B); stb memory-bound.
m1-fast5-sp interim 2 (09-25 22:00): vcol_team 9.45 -> 6.89 ms/build (v2), wedge 1 GPU 31.9 -> 29.2 ms/cycle (-8.5 %), 2 GPU 19.46 -> 18.1; CPU bitwise, GPU NOT bitwise vs b0 after 60 cycles (suspected FMA contraction change; column-dump job 11980296 checks); v3 (ILP segment weights) timing 11980297.
m1-sp-order2c (docs) merged 09-25 23:00: space-only order 2 confirmed on all transients; TIME order on the pure-scattering atmosphere transient = H-ESDIRK2 stiff order reduction (1.40 at c rho kappa dt 16 -> 1.92 at 0.25); magnitude 1e-7 relative vs spatial 1e-3. Fix would be a stage-order-2 ESDIRK (3 implicit solves).
m1-fast5-box interim 09-25 23:20: implicit_precond = mg_gf (global lowest 2-D Fourier modes |k|<=2 per x1 layer, 25 modes, Galerkin, pentadiagonal x1, multiplicative projection + mg3): box it/solve 8.5 -> 4.5 (1 GPU), 8.7 -> 5.4 (2 GPU); ms/cycle 44.11 -> 41.65 (-5.6 %) 1 GPU, 50.15 -> 49.28 (-1.7 %) 2 GPU (host Allreduce 25x84 per application); same fixed point as mg3 (1.5e-10 level). Gates running; then levers 2/3.
