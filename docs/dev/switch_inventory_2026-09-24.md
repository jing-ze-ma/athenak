# Switch inventory, rt-integration @ HEAD (2026-09-24)

This is a read-only inventory of the runtime switches this project added since about 2026-08. Keys were collected from
`GetOrAdd*` / `Get*` / `DoesParameterExist` calls at HEAD, and the dates come from `git blame`.
Scope: `<rad_m1>` (all keys); `<problem>` ck_* / rt_* in deep_hot_jupiter_rt.cpp,
two_stream_rt.hpp and two_stream_column_ck.hpp; the recent `<problem>` keys in box_convection / red_giant;
`<hydro>`, `<mhd>`, `<mesh>`, `<time>` keys added for the cs/sp work. shearing_box is out of scope.
Predecessor: `docs/dev/switch_inventory_2026-09-22.md/.csv` (432 rows, before the 09-22 clean-up). Rows that
have not changed since then are carried in section 6 in compact form.

Recommendation tags: **DEFAULT-ON** (validated, pays, no downside) / **KEEP-OPTION** (situational) /
**REMOVE** (failed or superseded) / **NEEDS-TEST** / **DIAG** (debug, dump, check) / **PARAM** (a physical or
tuning number, not a switch). The "gate" column gives two things. The first is "bitwise-off": the default path
is unchanged by the switch's existence. The second is whether correctness was tested (yes, partly, no, or FAILED).
All timings are GPU (MI300A) unless marked CPU. "E / V / Vh" are the Eddington, vet_sc and vet_sc+hesdirk2
He-box M1 arms of tests_m1.

## 0. Summary

Counts, over sections 1-5 (one row per switch or per grouped family; PARAM and DIAG rows grouped):

| recommendation | count |
|---|---|
| DEFAULT-ON (including those already on and to be kept) | 36 |
| KEEP-OPTION (and KEEP) | 35 |
| REMOVE | 52 |
| NEEDS-TEST | 17 |
| DIAG | 17 |
| PARAM (grouped rows) | 7 |

Counted by the first tag in the rec column; a row may carry a secondary tag, e.g. "DEFAULT-ON (3); REMOVE levels 1-2".
Of the DEFAULT-ON rows, **10 would flip a default that is off today**: implicit_gas_newton, vet_tensor = full,
implicit_closure_lag = step, f_source = wb (explicit only), ck_impl_nosync, ck_impl_arat = 1e30 under ck_implicit,
ck_sweep_form = tm under ck_implicit, implicit_krylov_dev (1 rank), grav_point_mass, rot_potential. The other 26 are
already on and should stay on. The T4 flags (frozen_op, lin, fuse, jac_lin, cvsec) are also DEFAULT-ON, but only
under ck_implicit and only after decision 1.

### The most important decisions

1. **ck_implicit (T4) + c2 + ck_impl_every = 4 for the next dhj production**: NEEDS-TEST, then DEFAULT-ON
   under ck_implicit. The user decided this on 09-24. Remaining blockers are a restart that is not bitwise
   with every > 1, and the 3-rotation A/B (ck-restart, jobs 11956956/7). Measured on GPU, RT cost relative to
   hydro: T4 7.24x, c2 2.03x, c2 + every 4 0.94x (tests_ck_implicit/README_fast.md §5, README_cadence.md §5).
   On MHD: e4 34.1 vs c2 52.2 ms/cycle (README_mhdsplit.md). Once gated, make the T4 set the defaults when
   ck_implicit = true, so inputs no longer need 13 lines (frozen_op, lin, fuse, jac_lin, cvsec, arat 1e30,
   nosync, and the c2 set once / xstep 2 / jreuse 0.2 / pred).
2. **time_scheme = hesdirk2 (+ implicit_vimp)**: DEFAULT-ON since m1-defaults2 (user, 09-24; tests_m1/runs_5g_defaults2/README.md). Earlier note: KEEP-OPTION, not default. It is 2nd order (G1 Edd 1.89-2.00,
   vet_sc 1.92-2.00) and restarts are bitwise. It costs 1.39x be+levers (60.5 vs 43.5 ms/cycle, 1 GPU;
   runs_4a_accel jobs 11954113/4). On the He box its accuracy edge is gone by t = 500 (runs_3z_validate).
   Side decision: the accel levers (vimp_fold, fast_kernels, one_pass, predictor_order 2) default on only
   under `be`. If hesdirk2 is kept, extend them to hesdirk2 after one G1 check at the stiff corners.
3. **implicit_halo_overlap**: KEEP DEFAULT-ON for now, contrary to the 09-24 handover's proposed flip. The
   loss with it on is <= 0.9 ms on 2 GPUs of one node (runs_4l_sync t13). Across nodes it gained 1.4-2.6 ms
   at 4 and 8 GPUs (3y f4/f8, pre-fix code). NEEDS-TEST: an overlap-off arm on the current code at 4 GPUs
   (job 11956312 / t4b has none).
4. **implicit_krylov_pipe**: NEEDS-TEST or REMOVE. Its gains were at noise level: 39.15 vs 39.2 ms on 2 GPUs,
   31.6 vs 31.85 on 4 strong. The one real gain was -1.4 ms on 8 GPUs weak (61.05 vs 62.45). Re-time it at
   4-8 GPUs on the current defaults, and remove it if the gain is under 1 ms.
5. **implicit_krylov_dev**: DEFAULT-ON for 1-rank runs, with dev_halo folded in as its only mode. It is bitwise,
   restarts are bitwise, and it saves 1.1-1.4 ms (E 34.7 vs 36.1, V 42.2 vs 43.3, Vh 58.4 vs 59.8; job 11955221).
   It does nothing on multi-GPU runs.
6. **vet_mb_tblock = 3**: DEFAULT-ON once there is an automatic tile rule (NEEDS-TEST of that rule). It is exact
   (25 arms bitwise). SC time per call drops 21 / 15 / 10 / 6 % on 1 / 2 / 4 / 8 GPUs, but ms/cycle stays
   within noise (runs_4i_sctb).
7. **implicit_precond_float**: REMOVE. It was neutral on 1 GPU and gave -1 ms on 2 GPUs with the pre-lever
   code; the README says "not recommended" (runs_3k_gpu3d/README_FAST.md).
8. **implicit_gas_newton**: flip to DEFAULT-ON. It is off by default, yet every GPU timing and gate since 3g
   ran it on (runs_3g_newton_T VERDICT PASS; 20 fallbacks per 3.29e6 cell-passes). Then run one GPU A/B of
   implicit_eos_cache, which is largely redundant with Newton (3g: within 0.3 %).
9. **implicit_solver default = line_jacobi** turns off the whole fast path, because `fdef` needs bicgstab.
   NEEDS-TEST: gate bicgstab as the default for transport = implicit on a multi-D mesh.
10. **Bulk REMOVE list** (failed or superseded): the Anderson family, implicit_closure_relax(_thin),
    implicit_trans_limit/fmax, dbg_trans_memory, the implicit_flux/blend family, implicit_recon plm_dc and its
    tuning keys, time2_vet_extrap, the vet_mb_* dead options, and the ck warm_step / cvkeep / jreuse_xc /
    xstep_thr / jac0 / jneg / lw options. Also rt_floor_consistent, rt_layer_legacy, restart_refill_ghosts,
    cs_vertex_fill_cc and reconstruct_x1.

### Default-on switches turned on recently (since 09-20)

| key | module | on since | condition of the default | evidence |
|---|---|---|---|---|
| ck_spherical | problem (ck) | 0e11d533, 09-22 | cs/sp + rt_ck, not rt_layer_legacy; fatal on Cartesian | tests_ck_sph/README.md §0 (plane-parallel: +0.63 budget error) |
| ck_beam_sph | problem (ck) | 0e11d533, 09-22 | same | tests_ck_sph/README_beam.md §6.6 (P/P_exact 1.000); CPU 1.08-1.12x cost |
| ck_sweep_form = tm | problem (ck) | 0172f3bc, 09-22 | ck_spherical, not ck_implicit, sweep_cache 2 | tests_ck_sweep_form/README*.md: exact; GPU 1.36x faster than four-pass |
| ck_sweep_cache = 2 | problem (ck) | 5f348164, 09-22 | ck_spherical (else 1) | tests_ck_sweep_cost/README.md: bitwise; +3.9..+7 % |
| rt_split, rt_use_cons, rt_col3_ex_iter, mlt_split_deposit | problem | 0e11d533, 09-22 | unconditional (dhj / box / RG) | 09-22 csv; rt_use_cons moved sp_prod3 bitwise |
| ck_impl_colskip | problem (ck) | ed188f69, 09-22 | only acts under ck_implicit | README_phase2/3 |
| ck_impl_sub_max 8 -> 32 | problem (ck) | be02c647, 09-23 | only with ck_impl_glob = ls_sub | README_glob (2000 s H2 columns need 32) |
| implicit_line_solver = pcr | rad_m1 | 66d5c59a, 09-23 | all | runs_3k README_PCR: kernel 15x, -152 ms/cycle production box |
| implicit_bcg_sync = 1 | rad_m1 | 66d5c59a, 09-23 | all | README_SYNC: 1.09-1.18x (1 GPU), 1.13x (2) |
| implicit_conv_est, implicit_lin_ew_max = 1e-2, implicit_lres_test off | rad_m1 | 66d5c59a, 09-23 | fixed-tensor closures (eddington/vet_sc/tau) | README_PICARD/DEFAULTS: Picard 4.58 -> 2.03; -14 ms (E) |
| implicit_predictor = step | rad_m1 | 66f1f889, 09-23 | fixed-tensor closures | README_PICARD/PREDRST: inner its -29 %; restart bitwise |
| implicit_halo_direct, od_cache, krylov_fuse 3, op_stencil, precond rbgs_fwd (the "fast path", fdef) | rad_m1 | 51a9adb2, 09-23 | fixed-tensor closure + transport implicit + bicgstab + bcg_sync 1 + pcr + one block along x1 + no lin_cnorm (op_stencil: not x1-periodic) | README_FAST / 3p_fastdefault: 113 -> ~48 ms (E, 1 GPU) |
| implicit_enthalpy = plm | rad_m1 | c27c7adc, 09-23 | implicit + fixed-tensor closure, not an old restart | runs_3s_space2, 4b_defaults: 2nd order, 4 % faster |
| implicit_halo_mpi | rad_m1 | c27c7adc, 09-23 | halo_direct + same-level mesh + not cs / sp / polar | runs_3w_krylov: -10..-17 ms at 2-8 GPUs |
| implicit_bc_advect | rad_m1 | 7066f13f, 09-23 | all, except old restarts | runs_4f_drift, 4h_bcadv: shock drift -1.64 -> -0.46 cells |
| implicit_fast_kernels, vimp_fold, one_pass 8, predictor_order 2 | rad_m1 | c50279ba, 09-24 | implicit + be + fixed-tensor closure, not an old restart (vimp_fold also needs vimp + op_stencil) | runs_4a_accel, 4j_accmerge: E 52.2 -> 40.5, V 58.1 -> 50.4 ms (job 11955024); radwave order borderline 1.86/1.89 |
| implicit_halo_overlap | rad_m1 | c50279ba, 09-24 | halo_mpi + > 1 rank, not an old restart | runs_3y, 4l (see decision 3) |
| implicit_halo_ovl_faces | rad_m1 | ecadc549, 09-24 | wherever the overlap is on | runs_4l_sync: 4 GPUs -1.6..-3.4 ms (job 11956312) |
| always-on fixes (no switch) | rad_m1 / mhd | 09-23..24 | n/a | ImplicitStencilOpPart vimp double count (68ea1fd0); MHD calls user_split_func (47bf0349) |

Note: on the cubed-sphere dhj, halo_mpi / overlap / ovl_faces resolve to off (hmdef), so their timings apply
only to Cartesian boxes.

## 1. `<rad_m1>`: time scheme, closures, transport physics

| key | default | what | evidence | gate | benefit / cost | rec |
|---|---|---|---|---|---|---|
| transport | explicit | explicit / implicit_x1 / implicit | all tests_m1 RESULTS; HANDOVER-09-24 | top selector | implicit is the only mode used on the He slab/box | KEEP-OPTION |
| time_scheme | hesdirk2 since m1-defaults2 where accepted (transport = implicit, Cartesian, integrator rk2, closure not tau / vet_col); be elsewhere and on restarts whose file lacks the key | be / hesdirk2 (H-ESDIRK2) | runs_3x_hesdirk2, 3z_validate, 4a_accel, 5f_h2fast, h2val_0924, 5g_defaults2 | G0 bitwise CPU+GPU; G1 order 1.89-2.00; restart bitwise; be named = bitwise base (5g) | 1.35x be per step (5f); at cfl 0.9 1.57x cheaper than be at cfl 0.3 at equal accuracy (h2val_0924) | DEFAULT-ON (tests_m1/runs_5g_defaults2/README.md) |
| time2_enth_vel | central (named only) | a_f build in hesdirk2 stages: old / start / central | runs_3x | old and start unstable at P=100, tau >= 1e5 | n/a | REMOVE old/start (keep central, hard-wired) |
| time2_vet_extrap | false (named only) | extrapolate the VET tensor in stages | runs_3x | FAILED (G1 order 1.59) | none | REMOVE |
| time2_dbg_fail | -1 | forces a stage-1 failure (fallback test) | runs_3x | n/a | n/a | DIAG |
| implicit_vimp | on where the resolved time_scheme is hesdirk2 and vimp is valid (multi-D, bicgstab, nghost >= 2, hydro + coupling); else off; off on restarts whose file lacks the key | Newton-implicit gas velocity in the enthalpy flux | runs_3v_vimplicit, 3u, 3z, 5g_defaults2 | bitwise-off; restart bitwise; NC 0; fatal on a 1-D mesh | +7.5 % E, +8 % V (job 11949611); under be it over-damps P=100 waves (0.22 vs 3.6e-3) | DEFAULT-ON with hesdirk2 only (tests_m1/runs_5g_defaults2/README.md) |
| implicit_vimp_jscale | 1.0 (named only) | Jacobian scale of vimp | runs_3v (0.5/1/2 agree to 5e-11) | n/a | n/a | DIAG |
| implicit_vimp_fold | on (be + fixed closure + vimp + op_stencil) | fold vimp terms into the stored stencil | runs_4a, 4j; bug fix 68ea1fd0 | bitwise-off; restart bitwise | -1.1 ms E, -2.0 ms V (hesdirk2) | DEFAULT-ON; extend to hesdirk2 (decision 2) |
| closure | m1 | m1 / minerbo / kershaw / eddington / vet_sc / tau / vet_col | runs_4b (m1 slab NC 6), 5b Finding 1 (m1, kershaw unstable in thin cells), 3j, 3l, 5d | per closure | 1 GPU box: E 40.5, V ~47.5 ms/cycle; vet_col +22 % vs E (job 11958342) | KEEP-OPTION (m1 is the explicit default); add a warning or fatal for implicit + chi(f) closures |
| vet_tensor | uniaxial (named only) | full 3x3 SC tensor | runs_3j README_FULLTENSOR; all later vet_sc gates use full | bitwise-off; NC 0; captures 12 % off-axis structure | +0.03 ms per SC call (job 11943125) | DEFAULT-ON for vet_sc (after a flip gate) |
| vet_axis | flux | uniaxial axis: flux / eigen | runs_3j V1_eigen | eigen FAILED (axis flips, KE_2 3.5e4x) | none | REMOVE eigen |
| vet_eig_min | 0 (named only) | eigenvalue floor of the full tensor | runs_3j gate 3 | partly | none | PARAM |
| vet_col_every | 1 | rebuild the vet_col tensor every k steps | runs_5d | restart bitwise only at multiples of k | build = 21 % of the vet_col run (11.9 ms/step) | NEEDS-TEST: accuracy and cost at every = 2..4 on GPU |
| vet_col_axis | radial | tensor axis r-hat / M1 flux | runs_5d thin-cell table | both stable | none | KEEP-OPTION (radial default) |
| vet_col_surface_q | true since m1-defaults2 (false on restarts whose file lacks the key) | outer-x1 Marshak q from the vet_col formal solution | runs_5e_vetcol2 | bitwise-off (sp gates); false named = bitwise base (5g) | T-S4 L1 2.1e-3 -> 5.2e-5 (n = 256); Milne 1.8e-3 -> 1.1e-4 | DEFAULT-ON (tests_m1/runs_5g_defaults2/README.md) |
| implicit_enthalpy | plm (implicit + fixed closure); else upwind | enthalpy face value upwind / central / plm | runs_3s_space2, 4b | bitwise-off; restart bitwise | 2nd order in space; 4 % faster E; Picard 2.27 -> 2.05 | DEFAULT-ON (done) |
| implicit_bc_advect | on (not on old restarts) | Marshak end face carries the enthalpy flux | runs_4f_drift, 4h_bcadv | bitwise-off; restart bitwise; NC 0 | drift -1.64 -> -0.46 cells; T_gas L1 7.7e-4 -> 8.3e-5 | DEFAULT-ON (done) |
| implicit_closure_lag | pass | closure recomputed per Picard pass, or frozen per step | runs_3b5 (recommends step for 2-D/3-D); every modern input sets step | step changes G-static in the 7th digit | pass ~4x slower per sim-sec (09-22 csv) | DEFAULT-ON (step, for implicit multi-D); needs one flip gate |
| implicit_offdiag | auto (= operator under bicgstab, else lagged; rad_m1_implicit.cpp:972) | off-diagonal Eddington terms: lagged / operator / none | runs_3b5; 5b Finding 3 (lagged has no fixed point in thin cells) | gated | lagging diverges at c dt >> dx | KEEP-OPTION (auto is right) |
| implicit_closure_relax / _thin | 1.0 / false | under-relax the closure update | runs_3b5, 5b Finding 1 | FAILED (not a cure) | none | REMOVE |
| implicit_trans_limit / trans_fmax | none / 1.0 | transverse realizability limiter | runs_3b8, 3d, 3h | FAILED ("only caps the instability") | none | REMOVE (fixed-tensor closures are the cure) |
| dbg_trans_memory | 1.0 | weight of the old transverse flux | runs_3b8, 5b Finding 1 | FAILED | none | REMOVE |
| implicit_flux | central | implicit face flux central / ap_hll / blend (transport=implicit fatals on non-central) | runs_3a2, 3b3, 3c | superseded | none | REMOVE non-central (unless implicit_x1 is kept) |
| implicit_blend, _fmode, _mode, _tau0, _flo, _fhi | tau_f / max / flux / 1.0 / 0.6 / 0.9 | tuning of implicit_flux = blend | runs_3c | inert (forbidden in multi-D) | none | REMOVE (family) |
| implicit_recon | dc | implicit face recon dc / plm_dc | runs_3a2, 3c | plm_dc FAILED (limiter limit cycle) | runs to maxit | REMOVE plm_dc (superseded by implicit_enthalpy = plm) |
| implicit_recon_w / _lag / _npass | -1 / picard / -1 | plm_dc tuning | runs_3a2, 3c | inert | none | REMOVE |
| implicit_bmom_half | true | half-step momentum write-back (a bug fix) | runs_3a2 / 3b / 3b3 | adopted | n/a | REMOVE the false option |
| implicit_bc_x1min / x1max | auto | column BC marshak / flux / reflect / periodic / efix | runs_4d_efix, He slab (flux) | gated | n/a | KEEP-OPTION |
| thick_flux | none | explicit AP flux none / ap_hll / scaled | runs_1c, runs_open T3b | ap_hll gated; none fails T3b | ap_hll is the only one that passes AP | KEEP-OPTION (none, ap_hll); REMOVE scaled + scaled_prefactor |
| ap_form | alpha2 | AP-HLL alpha form | runs_1c | unified worse on every 1b gate | none | REMOVE unified |
| advect_split | on iff hydro + ap_hll | explicit: advect E,F with the gas separately | runs_1c | gated | n/a | KEEP-OPTION (explicit only) |
| split_vel | recon | advective-split velocity recon / cell | runs_1cB | cell = control | n/a | REMOVE cell |
| source_form | full | full / ovc (O(v/c) coupling) | code comment: ovc is a failing control | FAILED by design | none | REMOVE ovc |
| f_source | cell | explicit F source at the cell / interface (wb) | runs_open T3b | bitwise-off; cell FAILS case B | flux error 1.35e-3 -> 1.8e-5 (A), 3.47e-3 -> 4.9e-4 (B) | DEFAULT-ON wb (explicit only; low priority) |
| coupling, gas_feedback | on / true | matter coupling; back-reaction | tests_m1 | adopted | n/a | KEEP-OPTION (false = ablation) |
| force_reference | none | residual force against the pgen WB a_rad_ref | runs_2a..3b6, He slab | adopted | n/a | KEEP-OPTION |
| opacity | const | const / powerlaw / table | tests_m1 | adopted | n/a | KEEP-OPTION |
| subcycle | true | explicit: subcycle transport in a hydro stage | 09-22 csv | adopted | no evidence found | KEEP-OPTION |
| opac_freeze, dbg_gas_force, dbg_gas_heat, dbg_gas_force_trans | false / true / true / true | ablations | runs_2a_diag, 3b6 | ablations | none | DIAG (REMOVE candidates) |
| dbg_tensor, dbg_tensor_tilt | none / 0.3 | frozen / tilt / tau tensor probes | runs_3i, 3l | tau mode superseded by closure = tau | none | DIAG; REMOVE the tau mode |
| dbg_opac_patch, dbg_opac_x1lo/x1hi/x2lo/x2hi | off (named only) | opacity multiplier in a box | runs_3j FULLTENSOR gate 3 | bitwise-off | n/a | DIAG (test fixture) |
| vet_milne, vet_dump, vet_dump_every, vet_col_dump, vet_col_dump_every | off | Milne test mode; dumps | runs_3j, 5d | n/a | n/a | DIAG |
| PARAM | | c_light, chat_over_c 1.0, e_floor, kappa_p/e/f/s, rho_ref, t_ref, opac_a/b, temp_unit_kelvin, rho_unit_cgs, kappa_unit, arad, marshak_q 0.5 (runs_5d: feed the formal-solution g into q next), implicit_flux_x1min/max, implicit_ebath_x1min/max, rsla_warn/vmax/force, reconstruct plm, cfl_rad, vet_nmu 4, vet_nphi 8, vet_col_ncore 8, vet_col_nsub 1, vet_col_nmu 4 | | | | PARAM |

## 2. `<rad_m1>`: implicit solver, Krylov, halo, SC sweep performance

| key | default | what | evidence | gate | benefit / cost | rec |
|---|---|---|---|---|---|---|
| implicit_solver | line_jacobi | x1 line Jacobi / matrix-free bicgstab | runs_3b2-3b5 | bicgstab is the gated multi-D path | the default turns off every fast-path default (fdef) | NEEDS-TEST: bicgstab as the default for implicit multi-D |
| implicit_gas_newton | false | Newton/Schur gas-temperature update | runs_3g_newton_T VERDICT PASS; 09-22 csv DEFAULT ON | bitwise-off; 20 fallbacks per 3.29e6 cell-passes | every GPU timing ran it on; no A/B | DEFAULT-ON |
| implicit_eos_cache | false | per-cell e(T) cache over the step | runs_3g, 3k EOS_CACHE_FLAG.md | bitwise-off; PASS; 2-rank diagnostic normalisation fix uncommitted (HANDOVER-09-23) | 3g: redundant with Newton (within 0.3 %) | NEEDS-TEST: GPU A/B with Newton on; REMOVE if < 1 ms |
| implicit_eos_cache_check | true | exact-table check of the cache | EOS_CACHE_FLAG.md | inert without cache | not measured | DIAG |
| implicit_accel + implicit_anderson_m/_beta/_start | none / 5 / 1.0 / 1 | Anderson acceleration of Picard | runs_3e_newton N-c | FAILED ("does not converge the thin top") | static slab 3.5x fewer passes; seeded slab fails | REMOVE (one_pass + predictor give 1.25 passes/step) |
| implicit_line_solver | pcr (09-23) | thomas / pcr | 3k README_PCR | thomas = bitwise reference; pcr round-off | -152 ms/cycle production box | DEFAULT-ON (done) |
| implicit_pcr_check | false | runs pcr and thomas, prints the difference | README_PCR | n/a | n/a | DIAG |
| implicit_bcg_sync | 1 (09-23) | host syncs 0 / 1 fused / 2 alpha on device | 3k README_SYNC | 0 bitwise old; 1/2 round-off | level 1 1.09-1.18x; level 2 within 1 % of level 1 | DEFAULT-ON (1); REMOVE level 2 |
| implicit_lin_cnorm | 0 | per-cell inner error bound | README_PICARD | disables fdef when > 0 | 137 vs 117 ms (slower) | REMOVE |
| implicit_halo_direct | fdef | on-rank copy kernel exchange | README_FAST, 3p | bitwise | -34 ms (113 -> 80, E) | DEFAULT-ON (done) |
| implicit_od_cache | fdef | cached off-diagonal sums | README_FAST | bitwise | -10 ms | DEFAULT-ON (done) |
| implicit_krylov_fuse | fdef ? 3 : 0 | fused preconditioner + vector updates | README_FAST | 1 bitwise; 2/3 round-off | fuse3 -5 ms vs fuse2 | DEFAULT-ON (3); REMOVE levels 1-2 |
| implicit_op_stencil | fdef and not x1-periodic | stored 19-point operator | README_FAST | round-off | 54.8 -> 48.6 ms | DEFAULT-ON (done) |
| implicit_precond | fdef ? rbgs_fwd : line | line / rbgs / rbgs_fwd | README_FAST | round-off | 65.3 -> 54.8 ms (rbgs_fwd) | DEFAULT-ON (rbgs_fwd); REMOVE rbgs |
| implicit_precond_float | false | preconditioner in float | README_FAST ("not recommended") | round-off | neutral on 1 GPU; -1 ms on 2 GPUs (pre-lever code) | REMOVE |
| implicit_op_split_red | false (named only) | stencil apply + separate reduction | runs_4a | round-off | 68.8 vs 69.2 ms | REMOVE |
| implicit_krylov_pipe | false | pipelined BiCGStab (needs fuse 3) | runs_3w, 3y | round-off; restart bitwise | 2 GPUs 39.15 vs 39.2; 4 strong 31.6 vs 31.85; 8 weak 61.05 vs 62.45 | NEEDS-TEST (4-8 GPUs, current defaults); else REMOVE |
| implicit_halo_mpi | hmdef | direct inter-rank implicit exchanges | runs_3w | bitwise vs the general path | -10 ms (2 strong), -16 (2 weak), -14..-17 (4-8 weak) | DEFAULT-ON (done; off on cs/sp) |
| implicit_halo_overlap | halo_mpi and > 1 rank | Krylov interior overlaps the MPI halo | runs_3y, 4l | round-off; restart bitwise | 2 GPUs/1 node: off faster by <= 0.9 ms; 4/8 GPUs across nodes: -1.4..-2.6 ms (pre-fix) | DEFAULT-ON kept; NEEDS-TEST overlap-off at 4 GPUs |
| implicit_halo_ovl_faces | wherever the overlap is on | shell only on MPI faces | runs_4l | round-off; pytest 3/3, gates 14/14 | 4 GPUs: E 42.1 vs 43.7, V 51.6 vs 54.2 (job 11956312) | DEFAULT-ON (done) |
| implicit_krylov_dev | 0 (named only) | BiCGStab scalars on the device (1 rank) | runs_4k_launch | bitwise-off; CPU bitwise on; restart bitwise | E 34.7 vs 36.1, V 42.2 vs 43.3 (job 11955221) | DEFAULT-ON for 1 rank |
| implicit_krylov_dev_halo | 0 (named only) | krylov_dev reads ghosts directly | runs_4k | CPU bitwise | -0.6 V, -1.0 Vh vs K=2 alone | fold into krylov_dev (REMOVE as a separate key) |
| implicit_op_check, implicit_op_check_tol | 0 (named only) | operator-equivalence check | m1-tests 72dbb7fd, tests_m1/gates | bitwise-off | caught the vimp_fold bug | DIAG (keep; used by pytest) |
| implicit_picard_log | 0 | per-pass log | README_PICARD | n/a | n/a | DIAG |
| implicit_lres_test / conv_est / lin_ew_max / predictor | closure-dependent (09-23) | Picard stopping rules and predictor | README_PICARD/DEFAULTS/PREDRST | PASS for fixed closures | see section 0 table | DEFAULT-ON (done) |
| implicit_fast_kernels / one_pass / predictor_order | be + fixed closure (09-24) | acceleration levers | runs_4a, 4j | vet_sc bitwise, E round-off; radwave order 1.86/1.89 borderline | E 52.2 -> 40.5, V 58.1 -> 50.4 | DEFAULT-ON (done) |
| implicit_opac_update | false | recompute opacity inside the iteration | 09-22 csv (UNSURE) | bitwise-off | no evidence found | NEEDS-TEST (T-dependent kappa A/B) or REMOVE |
| implicit_allow_multid | false | allow implicit_x1 on a multi-D mesh | runs_3a | safety override | none | REMOVE (transport = implicit covers multi-D) |
| implicit_partition | none | gather the x1 column across blocks and ranks | HANDOVER-09-24 (sc_x1_0924) | correct | x1 split +10-27 % M1 cycle | KEEP-OPTION (only if x1 must be split) |
| implicit_res_floor | 0 | absolute floor on the Picard residual | runs_3a2 | no-op | none | REMOVE |
| vet_mb_agg | true | one band message per neighbour rank | runs_3q | exact | 5.65 -> 4.71 ms/call (2 GPUs) | DEFAULT-ON (done); REMOVE the false path |
| vet_mb_mom_fuse | false | moments inside the ray launch | runs_3q | exact | slower (7.47 vs 3.92 ms/call) | REMOVE |
| vet_mb_agroup | 1 | rays split over G ranks of a node | runs_3q, 4i | round-off | 2 GPUs 4.15 vs 4.71 ms/call; worse across nodes; +106 MB/rank | KEEP-OPTION (single node) |
| vet_mb_overlap | false | interior during band MPI | 3m README_SCALING | exact | nothing (5.71 vs 5.64) | REMOVE |
| vet_mb_kernel | ray | ray / cell (old) kernel | 3m | exact | ray 18.2 -> 4.0 ms/call | REMOVE cell |
| vet_mb_tblock (+ _tj, _tk, _team) | 1 (off) | temporal blocking of the SC sweep | runs_4i_sctb | exact (25 arms bitwise) | SC -21 / -15 / -10 / -6 % (1 / 2 / 4 / 8 GPUs); ms/cycle within noise | NEEDS-TEST (auto tile rule), then DEFAULT-ON B=3 |
| vet_mb_tblock_ovl | true | overlap of blocked launches with band MPI | runs_4i | exact | no gain; interior pass slower | REMOVE |
| vet_mb_tblock_direct | true | direct neighbour-plane reads | runs_4i | exact | 3.24 -> 3.08 ms/call | DEFAULT-ON (done); REMOVE the off switch |
| vet_mb_angles | false | angle decomposition with global gather | 3m, 3q | drift 2e-11 | 4.19 ms/call; O(N_global) memory | REMOVE (superseded by agroup) |
| vet_mb_force, vet_mb_lag | named only | one-block multi-block sweep; block-Jacobi emulation | 3m / 3q | n/a | n/a | DIAG |
| vet_x1_periodic | false | periodic x1 SC sweep | runs_3r_radwave | tested | n/a | KEEP-OPTION (radwave tests) |
| PARAM | | cfl_rad 0.4, implicit_cfl -1, implicit_tol 1e-8, implicit_maxit, implicit_lin_tol 1e-10, implicit_lin_maxit 200, implicit_lin_ew_gamma 0.9, implicit_eos_cache_nt 2, implicit_pcr_team 0 (auto), implicit_one_pass_safety 3, vet_x1_npass 64, vet_mb_mom_batch 4, vet_mb_halo 3 | | | | PARAM |

## 3. `<problem>` correlated-k / dhj RT (deep_hot_jupiter_rt.cpp, two_stream_*.hpp)

Production prod4 (bench/cs_mhd_prod4/deep_hot_jupiter.athinput) sets: rt_ck, rt_use_cons, ck_spherical,
ck_beam_sph, ck_sweep_form = 1, ck_nquad 1, ck_pcut_bar 10, rt_split = false (the pgen forces true). Its
ck_implicit is off (semi-implicit).

| key | default | what | evidence | gate | benefit / cost | rec |
|---|---|---|---|---|---|---|
| ck_implicit | false | backward-Euler Newton ck column solve (T1-T4) | tests_ck_implicit/README*.md, wellposed/ | bitwise-off; 20 s: ~35x more accurate than semi; only T4 holds equilibrium at 200-2000 s | RT / hydro: T4 7.24x, c2 2.03x, c2+every4 0.94x | NEEDS-TEST (ck-restart + 3-rot A/B), then DEFAULT-ON with T4+c2+every4 |
| ck_impl_frozen_op | false (T4 = true) | frozen exchange operator over the Newton passes | README_phase3, tests_ck_kernel_tags | bitwise-off | part of T4 | DEFAULT-ON under ck_implicit |
| ck_impl_lin | false (T4 = true) | linear re-apply kernel of the tm sweep | README_T3 | bitwise-off; lin vs chain 4e-15 | T3 2.73x vs T1 4.58x (vs semi) | DEFAULT-ON under ck_implicit |
| ck_impl_fuse | false (T4 = true) | residual + tridiagonal in one kernel | README_T4 §3 | bitwise-off; same fixed point | 25.67 -> 22.19 s | DEFAULT-ON under ck_implicit |
| ck_impl_jac_lin | false (T4 = true) | Jacobian from the stored lin factorisation | README_T4 | bitwise-off | 22.19 -> 20.09 s | DEFAULT-ON under ck_implicit |
| ck_impl_cvsec | false (T4 = true) | secant c_v in the Jacobian rows | README_T4 | bitwise-off; passes 6.20 -> 4.50 | 20.09 -> 18.26 s | DEFAULT-ON under ck_implicit |
| ck_impl_arat | 2.0 | two-level split threshold A <= arat*E | README_T1_stall | FAILED at 2 (9 of 10 ranks non-converged) | all T4/c2 runs use 1e30 | DEFAULT-ON 1e30 under ck_implicit (drop the split) |
| ck_impl_nosync | false | no per-pass allocations or host syncs | README_fast §2, §5 | bitwise on CPU and GPU | 126.3 -> 120.9 ms/cycle (-4.3 %) | DEFAULT-ON |
| ck_impl_once | false | whole ck RT once per cycle (user_split_func; MHD since 0f196ca4) | README_fast §4-5, README_mhdsplit | off = old path; day rms 3.0e-4 vs t4 | RT 3.83x vs 7.24x; MHD 78.2 vs 131.1 ms/cycle | KEEP-OPTION (in c2; next production) |
| ck_impl_xstep | 0 (c2 = 2) | keep the frozen operator over k cycles | README_fast | bitwise-off; kinks +5 % at k4 | c2 46.5, c4 42.3, c8 40.6 ms/cycle; c8 FAILED (103 non-converged) | KEEP-OPTION (k = 2) |
| ck_impl_jreuse (+ _act) | 0.0 (c2 = 0.2) / 0.25 | chord Jacobian with a contraction test | README_fast | bitwise-off; within tol | RT 6.2x vs 7.24x | KEEP-OPTION (c2); supersedes reuse_jac |
| ck_impl_pred (+ _fac) | false (c2 = true) / 0.5 | drop the confirming pass | README_fast | bitwise-off; within tol | 4.00 -> 3.05 passes; RT 6.46x | KEEP-OPTION (c2; DEFAULT-ON candidate under ck_implicit) |
| ck_impl_every | 1 | full implicit call every N cycles, linearised step between | README_cadence, README_mhdsplit | bitwise at 1; e4 day rms 2e-5 vs c2, 0 non-converged; restart NOT bitwise | RT 0.94x vs c2 2.16x; MHD 34.1 vs 52.2 ms/cycle | NEEDS-TEST (restart fix + 3-rot A/B), then 4 in production |
| ck_impl_every_thr | 0.0 | per-column refresh guard | README_cadence | accurate at N = 8/16 | does not pay: 39.1 vs 24.7 ms/cycle | REMOVE (unless a cheap partial call is built) |
| ck_impl_glob (+ ls_c, ls_ntry, sub_max 32) | none | line search / per-column sub-stepping | README_glob | bitwise-off; 2000 s passes at sub_max 32 | same cost as T4 | KEEP-OPTION (dt >> production) |
| ck_impl_esc (+ _rho, _extra) | 0 | cheaper ls_sub escalation | README_jac | bitwise-off; never fires at the production dt | 18.1 -> 14.7 passes at 2000 s (CPU) | KEEP-OPTION (with ls_sub) |
| ck_impl_aa (+ _rst) | 0 | Anderson on the chord iteration | README_jac | bitwise-off; A2 200 s 200x better | same cost as T4 | KEEP-OPTION (large dt); NEEDS-TEST with c2 |
| ck_impl_seed | 0 | pass-0 guess (1 semi, 2 per-cell BE) | README_phase4 | bitwise-off | never measured with c2 | NEEDS-TEST (c2 +/- seed 2), else REMOVE |
| ck_impl_colskip | true | skip converged columns | README_phase2/3 | only under ck_implicit | n/a | DEFAULT-ON (done); REMOVE the off branch |
| ck_impl_frozen_cof | true | store half-layer coefficients | README_phase3 | false gives an identical answer | memory only | REMOVE the false branch |
| ck_impl_reuse_jac | 0 | old chord Newton | README_phase4, README_T4 | bitwise-off | superseded by jreuse | REMOVE (after c2) |
| ck_impl_warm_step | false | force the pass-0 step on a warm start | README_fast | FAILED (drift 1.2e-5) | no gain | REMOVE |
| ck_impl_cvkeep | false | reuse the previous call's secant c_v | README_fast | mispredicts with pred | no gain | REMOVE |
| ck_impl_jreuse_xc | false | carry the chord Jacobian across calls | README_fast | FAILED (209 non-converged at c8) | n/a | REMOVE |
| ck_impl_xstep_thr | 0.0 | rank-wide refresh threshold | README_fast | fires every call | 60.8 vs 42.3 ms/cycle | REMOVE |
| ck_impl_jac0 | false | Jacobian at e^n | README_T4 | FAILED (168/260 non-converged) | n/a | REMOVE |
| ck_impl_jneg | false | keep negative off-diagonal parts | README_T4 | worse gap | no gain | REMOVE |
| ck_impl_lw | false | light-weight launches | README_T4 | bitwise-off | no effect | REMOVE |
| ck_impl_lin_thr | 1 | block-threaded lin variant (4) | README_T3 | identical | 4 slower (2.88 vs 1.99 ms/pass) | REMOVE option 4 |
| ck_impl_tau_min | 0.0 | own-dtau thin/thick split | README_phase2 | "the wrong one" | n/a | REMOVE (superseded by arat) |
| ck_impl_refresh_kappa | false | opacity per Newton pass | code note only; fatal under T3+ | no evidence found | n/a | REMOVE |
| ck_impl_warm | false | warm start of the Newton | README_phase2/fast | no gain recorded | no evidence found | NEEDS-TEST or REMOVE |
| ck_sweep_form | tm if ck_spherical and not ck_implicit (and cache 2), else four-pass | 0 four-pass, 1 tm, 2 sd | tests_ck_sweep_form | exact (tm vs sd 3.5e-17) | tm 1.36x faster than four-pass | DEFAULT-ON tm also under ck_implicit (the T3/T4 lin path needs tm); keep sd as a cross-check |
| ck_sweep_cache | 2 (sph) / 1 | cache per-(cell, chain) layer operator | tests_ck_sweep_cost | bitwise | +3.9..+7 % | DEFAULT-ON (done); REMOVE mode 0 |
| ck_spherical, ck_beam_sph | on for cs/sp | see section 0 table | | | | DEFAULT-ON (done) |
| ck_int_at_cut | true | integrate the flux at the pressure cut | 09-22 csv | on since 09-07 | no evidence found | KEEP-OPTION |
| rt_floor_consistent | false | cap RT cooling at floor energy | tests_ck_topslab §3a | measured no-op | none (the fix was pfloor 1e-5) | REMOVE |
| rt_layer_legacy | false | pre-d8f00d49 staggered layers | 09-22 csv | refused by ck_implicit/sph | n/a | REMOVE (user kept it on 09-22; ~305 dead lines) |
| ck_impl_pred_chk, ck_impl_lin_check, ck_impl_verbose, ck_impl_debug | off | extra residual sweep; lin vs chain; reports | README_fast, T3, T4, jac | n/a | n/a | DIAG |
| rt_test_freeze, rt_test_dt/out/every/col_m/col_k/col_j, rt_test_rho/T/dT/nwave, rt_test_mu0, ck_test_kgrey | off | RT-only harness and test states | tests_ck_implicit/wellposed | inert off | n/a | DIAG (the wellposed suite is the ck accuracy gate) |
| ck_dump_*, rt_outer_verbose, rt_apply_debug(_n), rt_cell_report, rt_report_every | off | dumps and reports | 09-22 csv, tests_ck_topslab | n/a | n/a | DIAG |
| PARAM | | ck_impl_tol / dtol 1e-8 (1e-10 costs 148.6 vs 131.1 ms/cycle), ck_impl_norm_eps, ck_impl_maxit 8, ck_impl_dtmax, ck_impl_demax, ck_star_teff, ck_pcut_bar, ck_nquad, ck_table, ck_data_dir, ck_swflux | | | | PARAM |

## 4. `<hydro>` / `<mhd>` / `<mesh>` / `<time>` (cs/sp work)

| key | default | what | evidence | gate | benefit / cost | rec |
|---|---|---|---|---|---|---|
| hydro+mhd/reconstruct_x1 | plm (09-22) | Mignone curvilinear ppm4 / wenoz on the cs/sp radial sweep | tests_cs_radial_recon/README.md | bitwise-off; no vertex gain (loop 18-22 % worse) | CPU only: ppm4 ~10x, wenoz ~5x per zone-cycle | REMOVE |
| hydro/cs_wellbalanced_src | false (09-23) | well-balanced gnomonic geometric source (hydro) | tests_hyd4/README_wbsrc.md, tests_cs_regions | bitwise-off; seams 1.2-2.1x worse on exact tests | GPU 0-2.4 % | KEEP-OPTION off (or REMOVE the hydro read) |
| mhd/cs_wellbalanced_src | false | same, MHD | tests_cs_regions | bitwise-off | no evidence found | KEEP-OPTION (low beta; NEEDS-TEST there) |
| mesh/cs_vertex_fill_cc | false (09-22) | cc cube-vertex ghosts from neighbour panels | tests_cs_vertex, tests_mpi_gate_0922 | bitwise-off; static 11x better, dynamic unchanged | CPU +0.5 % | REMOVE |
| mesh/cs_vertex_fill | true (08-31) | FC buffers widened for the cube vertex | 09-22 csv | adopted | n/a | DEFAULT-ON (done) |
| mhd/hlld_bx_zero_tol | 1e-4 (09-22) | HLLD rotational-state threshold | HANDOVER-09-22d, bench/bisect_me | rank-invariant | 1e-8 caused the ME-decay regression | PARAM; NEEDS-TEST a proper fix (Bx^2 vs total p) |
| mhd/bs_emf | false | plain four-face corner EMF | memory cs-gs07-emf-gate | bitwise-off | smooth gate L1(B) 1.21e-6 vs 2.11e-6 (GS07) | DIAG |
| mhd/fofc_rsolver | llf | FOFC fallback solver llf / hlle | memory fofc-compatibility-plan | bitwise-off | no evidence found | NEEDS-TEST or REMOVE hlle |
| mhd/polar_hlle_rows | "6" (09-13) | polar-row HLLD -> HLLE mask | memory sp-polar-hlle-swap-is-the-excess | yes | x1 swap = the sp energy excess (KE 7.4x) | KEEP-OPTION (sp) |
| mhd/polar_emf_diss, mhd/cs_lowbeta_fallback, mesh/polar_x3_shift, mesh/use_polar_average_eresist | geometry-conditional on | sp/cs stabilisers | 09-22 csv (KEEP) | adopted | n/a | DEFAULT-ON (done) |
| hydro+mhd/sp_wellbalanced_src | false | WB geometric source on sp | tests_cs_regions; 09-22 csv UNSURE | cs bitwise | interior 2.33e-3 -> 1.56e-3 | NEEDS-TEST (on an sp dhj A/B) or REMOVE |
| hydro+mhd/sp_cart_all_momentum | false | Cartesian momentum on every sp cell | memory sp-pole-fixes | converges at the pole (2.6e-3 -> 5.2e-4) | never A/B'd on dhj | NEEDS-TEST (sp only; cs is production) |
| hydro+mhd/sp_cart_polar_momentum | false | Cartesian momentum in polar rows | memory sp-pole-fixes | FAILED (worse) | n/a | REMOVE |
| hydro+mhd/sp_face_avg (+ _terms 3) | false | 4th-order face averages (sp) | no evidence found | n/a | no evidence found | NEEDS-TEST or REMOVE |
| mesh/polar_quadratic_recon | false | quadratic polar-row reconstruction | memory sp-polar-row-reconstruction | shock-safe | polar residual 7.7x smaller | KEEP-OPTION (sp only) |
| mesh/use_grid_stretch_r_poly | false (08-26) | polynomial radial stretch | no evidence found | n/a | n/a | KEEP-OPTION (grid input) |
| hydro+mhd/wb_cache_every | 0 | rebuild the WB cache every N cycles | fmode-wb-cache-culprit, tests_cs_regions | 10 FAILED (spurious f-mode) | 1 is the fast safe value | KEEP 0 (user 09-22); forbid N > 1 |
| hydro+mhd/wb_rmin, wb_rmax | 0 | WB x1 reconstruction radial window | 09-22 csv | n/a | n/a | PARAM |
| hydro/ausm_mcut, ausm_mcut_p, ausm_wall_hllc | 1e-13 / 1.0 / false (09-07) | AUSM tuning; HLLC at walls | no evidence found | n/a | no evidence found | NEEDS-TEST or REMOVE ausm_wall_hllc |
| hydro/scratch_level | 0 (09-12) | team scratch level (1 = global) | memory nx1-ceiling-lds | yes | GPU +1.1 %; lifts the nx1 264 cap | KEEP-OPTION (nx1 > 264) |
| time/restart_refill_ghosts | false (09-12) | old restart ghost refill | memory restart-bc-reads-w0-bug | the old path is not a bitwise restart | n/a | REMOVE |
| time/dt_min, time/nan_check_cycles | 0 / 100 | dt floor abort; NaN check cadence | no evidence found | n/a | n/a | DIAG |
| mhd/use_rkg_sts, min_xe, max_eta | false / required | resistivity STS and eta caps | 09-22 note (max_eta dt cap) | n/a | n/a | PARAM (use_rkg_sts KEEP-OPTION) |

## 5. box_convection / red_giant `<problem>` keys (recent; verdicts carried from the 09-22 csv)

The 09-22 clean-up deleted the ABANDON switches outside rad_m1; the rows still ABANDON or UNSURE at HEAD are:

| key | default | what | evidence (09-22 csv) | rec |
|---|---|---|---|---|
| rt_strang (RG) | false | Strang-split the grey two-stream | HANDOVER-09-14: 2.1x WORSE at cfl 0.3; 09-20 retraction | REMOVE (but mlt_split_deposit needs it; check the He4 users first) |
| rt_once_per_cycle (RG) | false | RT once per cycle | REFUTED (bench/bstar_fecz/rkstage) | REMOVE |
| rt_split_transverse (RG) | true | move ADI with the column | inert in every production | REMOVE |
| rt_top_re, rt_semi_lin | false (col3) | corona RE top; old linearisation | fatal under mode 3; back-compat only | REMOVE |
| rt_relax_sub, rt_relax_xcrit, rt_relax_submax | 1 / 1.0 / 32 | sub-cycled relaxation | tests_ck_sph §192 (named next) | NEEDS-TEST or REMOVE |
| rt_impl_reuse, rt_impl_reuse_rho | 0 / 0.3 | block-factorisation reuse | no result recorded | NEEDS-TEST or REMOVE |
| rt_impl_redpar | false | PCR reduced system | 1.17x He column, 1.04x gate box | KEEP-OPTION |
| rt_col3_sub | 1 | apply the operator N times per stage | f-mode arm, no verdict | REMOVE |
| rt_weights_per_stage | false | rebuild the tau weights per stage | growth unchanged, 1.24x wall | REMOVE |
| rt_force_tau_gate (+ _lo, _hi) | false | optical-depth gate on the radiative force | 09-20 winning arm G7 has it off | REMOVE (or KEEP-OPTION for He4) |
| stellar_tide | false | host-star tide | validated, 0.5 % at the limb | KEEP-OPTION |
| grav_point_mass, rot_potential (dhj) | false | g ~ r^-2; centrifugal potential | both dhj productions set true | DEFAULT-ON |
| wb_phi_eff, wb_arad_force, wb_arad_file, m1_ic_file, m1_seed_consistent, m1_top_bc, planck_table (box, 09-21) | off / vacuum | M1 He-box inputs | 09-22 csv KEEP | KEEP-OPTION |
| mlt_split_sync, mlt_sync_passes | false / 1 | [1 2 1] sync of the MLT deposit | 09-22 csv KEEP | KEEP-OPTION |
| e_ledger, face_budget, rt_profile_dt/file, rt_surface_dt/file, mlt_dump(_dt), rt_force_verbose, rt_budget_verbose, open_debug, diag_gid (+ cycle window) | off | budgets and dumps | 09-22 csv KEEP | DIAG |

## 6. Everything else (unchanged since 09-22, KEEP)

The other ~150 `<problem>` keys of box_convection / red_giant / dhj are initial-condition, boundary, damping and
seed inputs (vpert_*, vdamp_*, bc_*, bg_*, open_*, mlt_* timing, sponge_*, cool_*, ic_*, opacity tables).
They keep the 09-22 verdict KEEP / PARAM; see docs/dev/switch_inventory_2026-09-22.csv for meaning and gates.
