# Default flips 2026-09-27 (branch defaults-0927)

User 09-27: "make sure all the recommended switches are default on".  These keys now default to
the recommended value.  **An input file that does not name a key changes behaviour.**  An explicit
value always wins, so naming the old value restores the old behaviour exactly.  Every value is
recorded with `GetOrAdd`, so a restart keeps whatever its file stored; where a key was read only
when named (so older restart files do not carry it), the new default is suppressed on restart
(`global_variable::restart_run`) and the old behaviour continues.

The ke-dt keys (`force_reference_work`, `implicit_opac_update`, `implicit_one_pass`,
`implicit_tol`) are handled on a separate branch (m1-keydefault-0927) and are not touched here.

| key | old default | new default | where the new default applies (otherwise old) | evidence |
|---|---|---|---|---|
| `problem/ck_nquad` (deep_hot_jupiter_rt) | 1 | 2 | every dhj run (only read with `rt_ck`) | user 09-25 (exact diffusion limit; nquad 1 ~10 % low); red_giant / box_convection unchanged |
| `problem/ck_impl_frozen_op`, `ck_impl_lin`, `ck_impl_fuse`, `ck_impl_jac_lin`, `ck_impl_cvsec` (the T4 set) | false | true | `ck_implicit = true` with the flags' own requirements: spherical ck, `ck_sweep_form = 1`, `ck_sweep_cache = 2`, no `ck_impl_refresh_kappa`, FP64 chain, `ck_impl_debug <= 0`; `jac_lin` also needs `ck_impl_lin_thr = 1` | switch_inventory_2026-09-24 sect. 3 (same fixed point; T4 = production solver since 09-24) |
| `problem/ck_impl_arat` | 2.0 | 1e30 | `ck_implicit = true` | arat 2 failed (README_T1_stall); every T4/c2 run used 1e30 |
| `problem/ck_sweep_form` | 1 only without `ck_implicit` | 1 (tm) also with `ck_implicit` | spherical ck with `ck_sweep_cache = 2` | tm exact; the T4 linear path needs it |
| `problem/ck_impl_nosync` | false | true | `ck_implicit = true` | bitwise on CPU and GPU, -4.3 %/cycle |
| `problem/ck_beam_par` | false | true | wherever its requirements hold (T4 path, sweep_form 1, ck_spherical, ck_beam_sph) | ck-fast2 (-25 % RT; differs only by beam deposits below e^-60) |
| `rad_m1/implicit_gas_newton` | false | true | `transport = implicit`, fresh runs | runs_3g_newton_T PASS; every GPU gate since 3g ran it |
| `rad_m1/implicit_closure_lag` | pass | step | `transport = implicit` on a multi-D mesh, fresh runs | runs_3b5; every modern input names step; pass ~4x slower |
| `rad_m1/implicit_closure_thin_relax` | 0 (off; read only when named) | 1.5 | chi(f) closures (m1, minerbo, kershaw), `transport = implicit`, multi-D, `implicit_closure_lag = step`, fresh runs | runs_5c_thinstab (thin-cell instability cure) |
| `rad_m1/vet_tensor` | uniaxial (read only when named) | full | `closure = vet_sc`, fresh runs | every vet_sc gate since runs_3j ran full |
| `rad_m1/implicit_precond` | rbgs_fwd (Cartesian fast path) | mg | Cartesian fast path on a multi-D, single-level mesh (not cs / sp / polar), `implicit_krylov_dev` not named > 0, fresh runs | runs_5m_precond: He box -7 %/cycle (1 GPU), -6 % (2 GPUs) |
| `rad_m1/implicit_mg_levels` (under mg) | 2 | 3 | every mg run | runs_5m_precond recommended value |

## Second round (branch defaults2-0927, 2026-09-27/28)

| key | old default | new default | where the new default applies (otherwise old) | evidence |
|---|---|---|---|---|
| `problem/ck_implicit` (deep_hot_jupiter_rt) | false | true | FRESH runs with `rt_ck`, `rt_use_cons`, not `rt_layer_legacy`, `ck_sweep_form != 2` (the switch's own refusals); **false on restarts** (a file written by an older binary records `ck_implicit = 0`; one written before the key existed ran the semi-implicit apply) | user decision 09-27 (the production scheme); the T4 set, `ck_impl_arat` 1e30 and `ck_impl_nosync` follow automatically (conditional on it) |
| `problem/ck_impl_floorbound`, `problem/ck_impl_kkt_demax` | false | true | fresh runs on the fused Newton (`ck_implicit`, `ck_impl_fuse`, `ck_impl_debug <= 0`, `ck_impl_glob = none`, `ck_impl_aa = 0`); false on restarts | bitwise vs off where no cell reaches a bound (CPU gates, floor = kkt = 0); where they bind (both WASP-121b productions, every call) the unconstrained Newton does not converge: 1x 466/500 calls NOT-CONVERGED off vs 6/500 on, 10x 500/500 vs 133/500 (/viper/ptmp2/jinma/defaults2_0927/B) |

The prod4 / hyd4 record inputs (inputs/production) now name `ck_implicit = false` (as run).
`tst` `test_rad_dhj_srclim_cpu.py` tests the semi-implicit `rt_de_max` limiter and now runs a copy
of its input that names `ck_implicit = false`.  Gates and numbers: /viper/ptmp2/jinma/defaults2_0927/RESULTS.md (viper).

Not flipped (user rules or open decisions): `ck_implicit` itself (flipped in the second round), the c2 levers (`ck_impl_once`,
`ck_impl_jreuse`, `ck_impl_pred`), `ck_impl_every`, `ck_impl_floorbound` / `ck_impl_kkt_demax` (second round: flipped),
`ck_pcut_bar`, `grav_point_mass` / `rot_potential`, `f_source = wb`; relaxation-only levers,
cfl 0.9, the floor switches and all tolerances.  Inventory with reasons:
/viper/ptmp2/jinma/defaults_0927/INVENTORY.md (viper).

## Third round (2026-09-28, speed-default audit on MI300A; user decisions 09-28)

Rule: a default that changes results (even at round-off) must buy >~ 3 % per simulated second
and stay within the round-off noise; otherwise it is opt-in by name.  All values are recorded in
the restart file, so a restart keeps the value it recorded; explicit input always wins.

| key | old default | new default | where | evidence |
|---|---|---|---|---|
| `problem/ck_impl_xstep` (deep_hot_jupiter_rt) | 8 wherever `ck_impl_every > 1` with the stored operator | 0 everywhere | fresh runs and restarts that record no value (the WASP-121b inputs name 8) | audit_dhj_0928: +3.2 % (1x) / +3.9 % (10x); deviation 1.3-2x the 1e-14-kick spread at 1e-3..1 and > 1 bar in the first rotation; relaxation-phase option only |
| `rad_m1/implicit_op_team_red` | true | false | everywhere (a restart without the key: true -> false) | audit_m1_0928: +2.2..2.6 % (1 GPU), 0.2-2.2 % (2 GPUs), below the bar |
| `rad_m1/implicit_vimp_fold` | true (fast-kernel path) | false | everywhere | audit_m1_0928: +0.9..2.3 % |
| `rad_m1/implicit_halo_ovl_faces` | true where the halo overlap is on | false | everywhere | audit_m1_0928: -1.4..+1.0 %, no gain |
| `rad_m1/implicit_one_pass` | 8 where not `force_reference = wb_arad` | 0 everywhere | everywhere | audit_m1_0928: -2.8..0 %, and 5-22x the round-off spread |
| `rad_m1/time2_lin_tol_fac` | 10 on fresh runs (1 on restarts without the key) | 1 | fresh runs | audit_m1_0928 rw/RESULTS_time.txt: factor 10 makes the stiff radiation-wave time error 1.6-4x larger (order 1.37 -> 0.48 at tau 10); gain -0.4..+7.5 % |

Kept (they pay): `ck_store_split`, `ck_beam_par`, `ck_impl_colskip` (part of `ck_impl_every > 1`),
`implicit_halo_overlap` (+9.2 % wedge on 2 GPUs in one node), `implicit_mg_levels` 3.  Open (user):
`implicit_gas_newton` (+8.5 % box, 1.6-1.7x noise).  Records: /viper/ptmp2/jinma/audit_dhj_0928 and
audit_m1_0928 (viper); merge gates: /viper/ptmp2/jinma/mergecoord_0928/RESULTS.md.
