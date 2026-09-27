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

Not flipped (user rules or open decisions): `ck_implicit` itself, the c2 levers (`ck_impl_once`,
`ck_impl_jreuse`, `ck_impl_pred`), `ck_impl_every`, `ck_impl_floorbound` / `ck_impl_kkt_demax`,
`ck_pcut_bar`, `grav_point_mass` / `rot_potential`, `f_source = wb`; relaxation-only levers,
cfl 0.9, the floor switches and all tolerances.  Inventory with reasons:
/viper/ptmp2/jinma/defaults_0927/INVENTORY.md (viper).
