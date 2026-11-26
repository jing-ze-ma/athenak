# Input keys of the VET radiation scheme

This page lists the keys a user sets for a run: the scheme, the rays, the time step, the
boundaries, opacities and units, the layout mode, and output. The production values are those
of the star templates in `inputs/stars/` (`star_bsg`, `star_m35mms`, `star_t9l62`:
3-D spherical-polar wedges, one MeshBlock along x1, fixed rays, problem generator
`star_vet`). Solver internals and keys without effect are listed only by name in the
appendix. A full sphere with poles or a stretched theta under `vet_gd` are in the
code but no deck or test covers them ([SCHEME.md](SCHEME.md), "In the code but not tested").

Columns: **default** is the value the code uses when the key is absent; **production** is
the value in the templates (`all:` when the three agree); **read in** is the file and
function that read the key (paths relative to `src/`).

Many defaults differ between a **fresh start** and a **restart** whose file does not hold
the key (the restart default keeps an old run unchanged). The start-up log prints every default
it applied together with the reason, e.g. `(default with vet_gd on a fresh start)`.

## The restart rule

A restart is started with `athena -r <file>.rst [-i <deck>] [block/key=value ...]`. The
parameters are put together in this order (`main.cpp`):

1. the input block **stored in the restart file**;
2. the `-i` deck, if given: every key it names replaces the stored value;
3. the command line: every `block/key=value` replaces the value so far. A key that
   is not in the deck (or the stored block) is added, and the code prints a `### NOTE`. On a
   **fresh** start such a key that the code never reads (a typo) stops the run after the
   first cycle (FATAL); a block that does not exist stops it at once.

**On a restart, a key that is not given in the deck or on the command line takes the value
stored in the restart file, not the code default.** When a `-i` deck is given, each such key whose stored value differs from the code default prints one `CLEAN WARNING restart: block/key = value from the restart file (code default d; the deck does not name it)` after the first cycle; the value is still used (name the key in the deck to choose; a restart without `-i` prints none). If the file does not hold the key either,
a key that corrects the scheme takes the fresh default, with a `### NOTE: ... the fix applies`
(name the old value to keep the old behaviour), and a stored old value of `merge_lat_mom`,
`merge_lat_blockend`, `fofc_flag_exchange`, `flux_balance_lat`, `time2_gtime_slope`,
`time2_gtime_pos` or `mlt_src_time` prints a `### NOTE` with the name to switch; the other keys take their restart default. For a
file of an older scheme use `restart_old_files.txt` (below).

Example: a restart file had stored `solver_maxit = 800` from an earlier deck. The fresh default
is 200, but the restart runs with 800 when the new deck does not name the key.

Rules:

- The decks name only what differs from the default, so `-i <deck>` does not replace the
  other stored values: name on the command line each key a file stores differently. For a
  file of older code see below. The grid keys
  (`<mesh>` x1min..x3max, nx1..nx3, nghost, the stretch keys and `<meshblock>` nx1..nx3)
  must stay as the file was written: a change is FATAL.
- `athena -r <file>.rst -i <deck> -n` prints the merged parameter list and quits (keys the
  modules add later are not shown).
- After the first cycle the code prints `### WARNING: N input parameter(s) were never read by
  the code`: a key there is a typo or a key this build does not have.

- The `<rad_vet>` keys `clean_*` lost "clean" (one method is left): the old names in a deck,
  on the command line or in a restart file are read as the new ones with one `### NOTE`; a
  source that names both spellings of a key with different values is FATAL. A restart file of
  this code stores the new names. Removed keys keep their old names (below). Log tags
  (`CLEAN cycle`, `CLEAN WARNING`, `clean_init steady pass`) are unchanged.

| old | new | old | new | old | new |
|---|---|---|---|---|---|
| `clean_bud_maxnorm` | `newton_maxnorm` | `clean_clos1_budget` | `newton_closure_budget` | `clean_closure` | `closure_form` |
| `clean_consistency` | `closure_consistency` | `clean_conv` | `newton_conv` | `clean_conv_kappa` | `newton_conv_kappa` |
| `clean_conv_max` | `newton_conv_max` | `clean_cs_lat_wb` | `flux_balance_lat` | `clean_cs_lz_k1` | `lz_cs_k1` |
| `clean_cs_lz_ledger` | `lz_cs_ledger` | `clean_cs_lz_tel` | `lz_telescope` | `clean_dbg_inf` | `dbg_inf` |
| `clean_init` | `steady_init` | `clean_init_dtfac` | `steady_dtfac` | `clean_init_n` | `steady_n` |
| `clean_init_tol` | `steady_tol` | `clean_lz_ledger` | `lz_ledger` | `clean_lz_ledger_k1` | `lz_ledger_k1` |
| `clean_lz_spin` | `lz_spin` | `clean_maxit` | `solver_maxit` | `clean_maxnl` | `newton_maxit` |
| `clean_pc_sweep` | `solver_sweep` | `clean_ray_cut_count` | `ray_cut_count` | `clean_ray_cut_now` | `ray_cut_now` |
| `clean_ray_cut_r` | `ray_cut_r` | `clean_ray_cut_wait` | `ray_cut_wait` | `clean_ray_min` | `closure_ray_min` |
| `clean_rot_fuse` | `ring_fuse` | `clean_rot_ring` | `ring` | `clean_sp_lz_bnd` | `lz_sp_bnd` |
| `clean_sp_lz_bnd_k1` | `lz_sp_bnd_k1` | `clean_sphere_radius` | `closure_sphere_radius` | `clean_stag_fcone` | `flux_cone` |
| `clean_stag_fcone_emax` | `flux_cone_emax` | `clean_stag_fcone_lz` | `flux_cone_lz` | `clean_stag_mg_levels` | `solver_mg_levels` |
| `clean_stag_opac_passes` | `gas_opac_passes` | `clean_stag_thinf` | `face_thin` | `clean_stag_xdonor_freeze` | `face_xdonor_freeze` |
| `clean_tol` | `solver_tol` | `clean_tol_nl` | `newton_tol` | `clean_verbose` | `verbose` |
| `clean_w_c0` | `closure_w_c0` | `clean_w_c1` | `closure_w_c1` | `clean_wb_exch_local` | `update_exchange_local` |
| `clean_wb_floor_diag` | `update_floor_diag` | `clean_wb_mcap` | `update_mcap` | `clean_stag_fcone_gas` | `flux_cone_gas` |
| `clean_stag_fwb_line` | `fwb_line` |  |  |  |  |

- The star keys `<problem>/star_*` were named `he_*` (one: `hs_*`) in earlier versions. Old
  names in a deck, on the command line or in a restart file are read as `star_*` with one `### NOTE`; a deck
  that names both spellings of a key with different values is FATAL.
  `tools/rename_star_keys.sh FILE ...` renames the keys in decks and scripts. The problem
  generator `-D PROBLEM=star_vet` was `he_star_vet`; CMake accepts the old name.

### Restart files of older code

Files of the old implicit solver, or of early versions of this scheme, store other values of
some keys, or lack keys whose restart default is not the production value.
`inputs/stars/restart_old_files.txt` names them, with the reason (commented out: what depends
on the star); use it as the `-i` file:

```bash
athena -r old.rst -i $ATHENAK/inputs/stars/restart_old_files.txt mesh/layout_tile=64 ...
```

| key | value | why |
|---|---|---|
| `<rad_vet>/scheme`, `clean_layout`, `clean_rays` | `clean`, `staggered`, `sc` | the old solver and its forms were removed: FATAL without them |
| `<rad_vet>/time2_k1_stiff`, `time2_k1_redist` | `0` | removed keys; older files store 1.0, 0.1 (FATAL) |
| `<rad_vet>/clean_stag_fwb_stage`, `clean_stag_fwb_t2` | `false` | removed keys; early files of this scheme store true (FATAL) |
| `<rad_vet>/clean_w_fse_tau`, `time2_gtime_mom` | `true` | removed keys; early files may store false (FATAL) |
| `<rad_vet>/solver_mg_levels` | `2` | early files may store 4 |
| `<rad_vet>/solver_maxit` | `200` | not stored; a restart keeps 800 |
| `<rad_vet>/vet_gd_wall_interp`, `vet_gd_wall_sym`, `vet_gd_theta_seam_sym` | `true` | spherical-polar wedges with fixed rays (commented out in the file); a theta-periodic wedge without them is FATAL |
| `<hydro>/fofc` | `true` | older star files may store false (HG, RSG, R136 store true) |
| `<hydro>/merge_lat_mom`, `merge_lat_blockend` | `radial`, `exchange` | commented out, a choice of physics: a stored `rigid` (RSG) or `flat` stays (NOTE) |
| `<rad_vet>/vet_rotate_avg_ema`, `_change_form`, `ring` | `128`, `linear`, `live` | commented out: stars with rotating rays (HG, RSG, R136 continuations) |
| `<problem>/mlt_src_time` | `stage` | files of the release ccab558ee store `start` (MLT stars) |
| `<mesh>/layout_tile` | the lateral size of the file's MeshBlocks | not stored (block files): the MeshBlock is the tile; name it to keep the answer on a new layout (commented out) |
| `<time>/restart_refill_ghosts` | `false` | some files store true |

A file this code wrote stores every key it read: nothing needs naming (a file of the release
ccab558ee: `mlt_src_time` only).

Rows marked "hydro only" have no `<mhd>` form: the key named under `<mhd>` (`floorcap_rho`
with `<mhd>`) stops the run. `fofc_retest` is hydro only too.

## Scheme selection

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/transport` | Radiation transport; the scheme needs `implicit`. | `implicit` on a fresh run with a `<rad_vet>` block | all: `implicit` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/closure` | Only `sc`: the Eddington tensor from the short-characteristics rays, on Cartesian meshes and by the global sweep on spherical polar and the cubed sphere. An older name is read as `sc` (NOTE) only where it ran the same closure: `vet_sc` on a Cartesian mesh, `vet_col` with `vet_gd = true` on a curved one; every other value stops the run. | `sc` | - | `rad_vet/rad_vet_params_closure.cpp` (constructor) |
| `rad_vet/closure_form` | How the closure follows the solution in time: `extrap1` (one ray sweep per step; recommended), `iterate`, `predict`, `eddington` (D = I/3, no rays; for tests). | `extrap1` | all: `extrap1` | `rad_vet/rad_vet_params.cpp` (constructor) |
| `hydro/reconstruct`, `mhd/reconstruct` | Reconstruction of the fluid sweeps (`dc`, `plm`, `ppm4`, `ppmx`, `wenoz`) on Cartesian meshes and in the cubed sphere's lateral sweeps (its radial sweep is the grid PLM). On spherical polar every sweep is the grid PLM, and with `wb_dirs = all` the well-balanced PLM: there a value other than `plm` stops the run (FATAL). `wb_x1`, `wb_rho`, `wall_wb_onesided` set without `wellbalance_dynamic = true`, `wb_guard_ideal` / `wb_clamp_ideal` with a general EOS, `mhd/fofc_rsolver = hlle` without `fofc`, `mhd/polar_emf_diss` without `use_polar_boundary` stop the run too (they would not act). `time/dt_retry_views = all` stops the run on the cubed sphere (its retried step differs at a ray-set change). | `plm` | not named | `hydro/hydro.cpp`, `mhd/mhd.cpp` |

## Rays

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/vet_gd_gl_nmu`, `vet_gd_gl_nphi`, `vet_gd_gl_stagger` | Product direction set: nmu Gauss-Legendre polar nodes times nphi azimuths, alternate rings shifted with `stagger`. 2 pi/nphi must divide the phi extent (FATAL otherwise). Takes precedence over the other sets. | 7, 18, `true` on a fresh 3-D spherical-polar run; else `vet_gd_geo = 3` (3-D) | all: 7, 18, `true` (126 directions) | `rad_vet/rad_vet_params_closure.cpp` (constructor) |
| `rad_vet/vet_gd_ls` | Level-symmetric set LQ_N, N = 6, 8, 10, 12 (N = 10: 120 directions). Naming it together with the product set (`vet_gd_gl_nmu` > 0) is FATAL: name one set (0 = off). | 0 (off) | not named | `rad_vet/rad_vet_params_closure.cpp` (constructor) |
| `rad_vet/vet_gd_geo` | Geodesic set, 10 n^2 + 2 directions. | 3 when no set is named | not named | `rad_vet/rad_vet_params_closure.cpp` (constructor) |
| `rad_vet/vet_sc_set` | Cartesian meshes: `gd` sweeps the global direction set of the curved meshes (`vet_gd_gl_nmu`, `vet_gd_gl_nphi`, `vet_gd_gl_stagger`, `vet_gd_ls` or `vet_gd_geo`; with none named the geodesic set 4, 162 directions) with its pole on x1, instead of the double-Gauss set about x1. Replaces `vet_ring` (FATAL with it); a mirrored x1 end needs a set with x1 mirror pairs (FATAL otherwise). Read only when named. | double-Gauss set about x1 | not named | `rad_vet/rays/sc_cart_setup.cpp` (constructor) |
| `rad_vet/vet_gd_rotate_dt` | Rotate the direction set every dt_rot of physical time (`auto` = 20 x the first dt; 0 = fixed set). | `0` on spherical polar, else `auto` (fresh runs that name no rotation key) | all: `0` | `rad_vet/rad_vet_params.cpp` (constructor) |
| `rad_vet/vet_rotate_avg` | With rotation: average the closure over the last K sets (0 = off). | 4 with a rotating set on a fresh run | all: `0` | `rad_vet/rad_vet_params_closure.cpp` (`RotAvgParse`) |
| `rad_vet/vet_rotate_avg_ema` | With rotation: the closure is an exponential moving average over the rotated sets, weight 1/N per block (0 = off; then `vet_rotate_avg` is the K-set window). | 128 on a fresh run with a rotating set; read only when named on a restart | RSG, HG: `128`; all: `0` (fixed rays) | `rad_vet/rad_vet_params_closure.cpp` (`RotAvgParse`) |
| `rad_vet/ring` | Rotating ray ring of the scheme: `off`, `live`, `jump`. | read only when named | all: `off` | `rad_vet/rad_vet_params.cpp` (constructor) |
| `rad_vet/vet_scatter` | Split the opacity into absorption and electron scattering in the ray source (Ma et al. 2026, eq. 7). Changes only the rays' source function. | true on a fresh run with table opacity and the rays | all: `true` | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/vet_scatter_kappa_e` | Electron-scattering opacity: `const` or `eos` (from the EOS ionisation). | `eos` with a tabulated EOS; `const` otherwise | all: `const` | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/vet_gd_wall_sym` | Phi-periodic wedges: centre the set on the wedge, mirror-paired phi-wall maps, the seam face from both sides. Needs `vet_gd_wall_interp = true`. Identical results on wedges that are whole azimuthal periods; required on narrower ones. With a set turning about z it needs a wedge of whole periods (exact wall maps for every turn) and does nothing; FATAL otherwise. No effect without a phi wall (cubed sphere, full 2 pi, 2-D). | false (recorded) | all: `true` | `rad_vet/rad_vet_params_closure.cpp` (constructor), `rad_vet/rad_vet_params_closure.cpp` |
| `rad_vet/vet_gd_theta_seam_sym` | Theta-periodic wedges: the theta seam face's other side is its exact equatorial mirror (both frames, as the phi seam). Required on a theta-periodic mesh (FATAL otherwise); the ray set must be mirror-symmetric about the equator (FATAL otherwise; e.g. GL with even nmu and stagger is not; GL with odd nmu is); a set turning about z keeps its mirror at every turn. No effect without a theta seam. | false (recorded) | all: `true` | `rad_vet/rad_vet_params_closure.cpp` (constructor), `rad_vet/rad_vet_params_closure.cpp` |

The wedge rule: on a phi-periodic spherical-polar wedge the width must be a multiple of the
ray set's azimuthal period (Gauss-Legendre: 360 deg / nphi, 20 deg for nphi = 18; LQ10: 90
deg; the default icosahedral set: 180 deg), else the run stops with a FATAL. A narrow wedge
that cannot meet it (e.g. 0.064 rad) names `vet_gd_wall_sym = true` with
`vet_gd_wall_interp = true` (the `wedge_sym` test). A theta-periodic wedge (all the star
templates) stops with a FATAL unless it names `vet_gd_theta_seam_sym = true`; a restart
from a file written before these keys must name both on the command line.

## Time integration, time step and retries

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `time/integrator` | Hydro integrator; must be `rk2` (the explicit part of the IMEX scheme). | `rk2` | all: `rk2` | `driver/driver.cpp` (constructor) |
| `time/cfl_number` | Courant number of the hydro step (per direction); the radiation enters dt by `rad_signal_speed` where it is trapped. | 0.3 on a fresh VET run | all: `0.3` | `mesh/build_tree.cpp` |
| `hydro/merge_lat_mom` | (hydro only.) Group value of the merged lateral cells (`merge_lat_rc`): `radial` keeps M, E, sum V m_r and L_x, L_y, L_z of every group to round-off and a spherically symmetric state unchanged; where uniform density would cost more kinetic energy than the group has, part of the density contrast stays (no momentum is dropped, thermal energy > 0 if it was before). `rigid`: the old Cartesian rigid field, which drops the group's rotation (L lost) where it exceeds half of the thermal energy. `merge_lat_blockend` (`exchange`, also on a restart whose file lacks it, NOTE; `flat`: the old rule, which depends on the MeshBlock layout). `MLATDIAG`: `nfb` (radial: groups without thermal energy before the merge, left to the C2P floors), `amin` (smallest merged share). `radial` records `merge_lat_state` in the restart file. A running `rigid` star restarted as `radial` re-averages its groups in the first step: its kinetic energy changes there once. | `radial`, also on a restart whose file lacks it (a NOTE in the log; name `merge_lat_mom=rigid` to keep the old) | bsg, RSG cs: not set | `hydro/hydro_mergelat.cpp` |
| `hydro/dt_floor_courant` | (hydro only.) Floor gas (rho < `rad_vet/floorcap_rho`) takes dt <= C / sum_d (abs(v_d) + c_d)/h_d (one Courant bound on the sum of the directions) instead of `cfl_number` per direction, so a floor atmosphere does not set dt; FATAL without `floorcap_rho`. `DTREP` prints `floor=1` for such a binder. Time step only. | 0.6 with `floorcap_rho`, else 0 (off) | all: `0.6` | `hydro/hydro.cpp` |
| `hydro/rad_signal_speed_diff` | (hydro only.) With `rad_signal_speed = true`: the radiation term of the signal speed weighted by 1 - exp(-x), x = tau_d c_eq / c the grid-scale trapping parameter; untrapped radiation (x << 1) does not limit dt. Time step only (accuracy; at cfl 0.3 stability does not need the radiation term). | true with `<rad_vet>` and `rad_signal_speed` | all: `true` | `hydro/hydro.cpp` |
| `hydro/rad_signal_speed_phase` | (hydro only.) With `rad_signal_speed_diff = true`: that weight is 1 - exp(-x^2), the in-phase (restoring) part of the radiation-pressure perturbation; the linear weight also counts its damping part, which grows as rho^-1/2 in thin gas. Same weight at x = 1, larger for x > 1. Time step only. Needs `rad_signal_speed` and `rad_signal_speed_diff` (FATAL otherwise). | true with `rad_signal_speed_diff` | all: `true` | `hydro/hydro.cpp` |
| `problem/mlt_src_time` | Time of the MLT source (`mlt_flux_frozen`: the ramp w(t), the `adaptive` / `pressure` closure's F_sub): `stage` takes it at the time of each implicit stage (second order), `start` at t^n for the whole step (first order; under `sdirk2` dt/2 late). A frozen flux without a ramp is the same either way. A restart file that stores `start` keeps it unless `stage` is named; the first step after the switch of an `adaptive` or `pressure` run takes the target of t^n for both ends (NOTE). | `stage` (a restart file without the key too, NOTE) | none of the templates has an MLT flux | `pgen/star_vet.cpp` |
| `rad_vet/time_scheme` | `hesdirk2`: two implicit stages per step with the Heun hydro (second order with `time2_tableau = hesdirk2`; with the production `sdirk2` by `time2_gtime_slope` where the gas-radiation exchange is stiff); `be`: one backward-Euler solve per step (first order). | `hesdirk2` on a fresh run that supports it; `be` otherwise | all: `hesdirk2` | `rad_vet/rad_vet_params_implicit.cpp` (`ImplicitInit`) |
| `rad_vet/time2_tableau` | Tableau of the two stages: `sdirk2` (Alexander's L-stable SDIRK2, recommended; with `<mhd>` the stage-1 solve takes the field B^n + g (B* - B^n) through \|B\|^2/2), `hesdirk2`, `trbdf2`. | `sdirk2` with `<hydro>` and `<mhd>` (also on a restart whose file lacks it, NOTE; a file that stores `hesdirk2` or `trbdf2` keeps it) | all: `sdirk2` | `rad_vet/rad_vet_params_implicit.cpp` (`Time2Init`) |
| `rad_vet/time2_gtime_slope` | sdirk2: the gas eint of the hydro's stage-2 flux state is that of U* + d plus (1 - g) dt minmod(K1, K2), K1 the last step's stage-2 slope of the exchange, K2 = d/(g dt), added to the gas the caps and floors left (Patankar form where negative; its floors act on the primitives only, nothing is booked). Without it the exchange is first order in time where it is stiff (wave p 1.0-1.5 instead of 2, M5 shock rho L1 +27 % over hesdirk2). Helps where gas is trapped or in fast exchange; no effect where the exchange is slow. Limit: the flux state of cold floor gas can move by a factor 3 (toward T_rad). | true, also on a restart whose file lacks it (NOTE) | all: default | `rad_vet/rad_vet_imex.cpp` (`Time2GasTime`) |
| `rad_vet/time2_gtime_pos` | sdirk2: the flux state is formed from the gas the caps and floors left, stays positive (Patankar) and books nothing; `false` = the older form. Must be true with `time2_gtime_slope` (FATAL otherwise). | true, also on a restart whose file lacks it (NOTE) | all: default | `rad_vet/rad_vet_imex.cpp` (`Time2GasTime`) |
| `rad_vet/time2_auto_be` | After a step whose stages fell back to backward Euler, the next 1, 2, 4, ... N steps take backward Euler at once (0 = off). | 16 | all: `16` | `rad_vet/rad_vet_params_implicit.cpp` (`Time2Init`) |
| `rad_vet/steady_init` | `steady`: before step 1 the radiation is relaxed to the scheme's steady solution in the given gas; `given`: start from the given field. | `given` (star_vet: `steady`) | all: `steady` | `rad_vet/rad_vet_params.cpp` (constructor) |
| `rad_vet/steady_dtfac`, `steady_tol`, `steady_n` | Pseudo time step (in units of the first dt), tolerance and pass limit of the steady start. | 1e6, 1e-8, 40 or 50 (problem dependent) | all: `1.0e6`, not named, not named | `rad_vet/rad_vet_params.cpp` (constructor) |
| `time/dt_retry` | A cycle whose implicit solve fails (not converged, non-finite, unpaid floor energy) is undone and taken again with a smaller dt. Bitwise: a retry restores the full cycle. | true | all: `true` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_fac` | dt factor of a retry. | 0.5 | all: `0.5` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_grow`, `dt_retry_ngood` | After a cut, dt grows back by this factor per good cycle, starting after this many good cycles. | 1.1, 20 | all: `1.1`, `20` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_force_cycle`, `dt_retry_force_fac` | Testing: force one retry at this cycle with this factor (the `dt_retry` test). | off | not named | `driver/dt_retry.cpp` (`DtRetry`) |
| `rad_vet/dbg_inf` | Testing: the first stage-1 write-back of this cycle gives one cell a non-finite E iterate. `implicit_pos_floor` never raises a non-finite E; the stage is rejected and redone (the `inf_reject` test). | off | not named | `rad_vet/rad_vet_params.cpp` (constructor) |

## Boundaries

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/implicit_bc_x1min`, `implicit_bc_x1max` | Radiation x1 boundaries: `marshak`, `outflow`, `reflect`, `periodic`, and on x1 an imposed `flux`. The ray sweep needs `x1max = marshak` or `reflect`. | `auto` (from the mesh boundary); star_vet with `implicit_flux_x1min`: x1min `flux` | all: `flux`, `marshak` | `rad_vet/rad_vet_params_implicit.cpp` (`ImplicitInit`) |
| `rad_vet/implicit_flux_x1min` | Flux through the inner face (= L/(4 pi r_in^2) for a star). | 0 | per deck | `rad_vet/rad_vet_params_implicit.cpp` (`ImplicitInit`) |
| `rad_vet/implicit_flux_x1max` | Incoming bath at the outer Marshak face (0 = vacuum). | 0 | all: `0.0` | `rad_vet/rad_vet_params_implicit.cpp` (`ImplicitInit`) |
| `mesh/ix1_bc` ... `ox3_bc` | Hydro boundaries (the templates: `user` in x1, `periodic` in theta and phi). `ix2_bc` = `ox2_bc` = `polar` needs `use_polar_boundary = true`, which needs both on spherical polar (FATAL otherwise). | required | see the decks | `mesh/mesh.cpp` |
| `problem/star_wall_zero_flux` | star_vet (hydro): the inner x1 wall's mass, transverse-momentum, energy and scalar fluxes are 0 before the update (not with `star_bc_inner = inflow`); so are the top's where the mass flux is inflow or exactly 0 at that point. `dfloor_flux_limit` acts later and closes a top face that it makes massless by the same rule. | true | all: `true` | `pgen/star_vet.cpp` (`StarWallFlux`) |

Transverse Marshak faces are not supported by the ray sweep (FATAL).

## Opacity, units and coupling

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/opacity` | `const`, `powerlaw`, `user` or `table` (Rosseland and Planck tables named in `<problem>`). | `table` with `problem/star_opac_table`, else `const` | all: `table` | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/kappa_p`, `kappa_e`, `kappa_f`, `kappa_s` | Constant opacities per unit mass: Planck (emission), energy mean (absorption), flux mean, scattering. | 0, = kappa_p, = kappa_p, 0 | not named | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/c_light`, `arad` | Speed of light and radiation constant in code units. The start-up log compares them with the values `<units>` implies. | required | per deck | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/temp_unit_kelvin`, `rho_unit_cgs`, `kappa_unit` | Unit conversions for the table lookup: K per code temperature, g/cm^3 per code density, code opacity per cm^2/g. | 1.0 | per deck | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/coupling`, `gas_feedback` | Gas-radiation exchange on; exchange acts back on the gas. | true | all: `true` | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/force_reference`, `force_reference_work` | `wb_arad`: well-balanced radiative envelope; the problem generator supplies the reference radiative acceleration, the gas receives only the residual force; its work is paid as `split` (required by the scheme). | `none` (star_vet: `wb_arad`); `auto` (= `split` there) | all: `wb_arad`, `split` | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/e_floor` | Floor of E (code units). | FLT_MIN | all: `1.0e-4` | `rad_vet/rad_vet_geometry.cpp` (constructor) |
| `rad_vet/floorcap_rho`, `floorcap_v` | (hydro only.) Floor-gas cap: in cells with rho below `floorcap_rho`, |v| is capped at `floorcap_v` and T at max(T_rad, floor); the removed energy goes to E (momentum is removed and reported). | off | all: `1e-14`, `3e7` | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/floorcap_lz` | The floor-gas velocity cap scales only the momentum without L_z: spherical polar r / theta, Cartesian the part off phi_hat = (-y, x)/R, the cubed sphere the part off phi_hat in the metric of the covariant momentum (`GnomonicPhiHat`, p_phi = m_a phi^a). The phi part is kept, also where it alone exceeds `floorcap_v`, so L_z is exact. (The cap itself measures \|v\| in the metric norm on the cubed sphere.) | true with `floorcap_rho` (recorded; a restart file without it: false) | all: default | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `rad_vet/floorcap_lz_spin` | With `floorcap_lz`: the floor-gas cap scales the whole velocity (so \|v\| <= `floorcap_v` in floor gas, which then cannot set dt through a kept v_phi), and the L_z X it removes goes to the gas at or above `floorcap_rho` as one rigid rotation dOmega = X/I, I = sum rho V R^2 (per cell, then in a fixed order: layout-invariant in tile mode); E_rad of each receiving cell pays its kinetic-energy change. L_z and E exact, no ledger. FATAL without `floorcap_lz`, or when the cap removes L_z and no gas is at or above `floorcap_rho`. | true with `floorcap_rho` (recorded; a restart file without it: false), and with `floorcap_lz` | all: default | `rad_vet/rad_vet_params_matter.cpp` (constructor), `rad_vet/rad_vet_gas.cpp` (`FloorGasCap`) |
| `hydro/dfloor`, `pfloor`, `tfloor` | Hydro floors. `hydro/floor_pay = rad_column` pays the internal energy they create from the radiation (needs `hydro/eos = general`). | FLT_MIN | per deck | `eos/eos.cpp`, `hydro/dual_energy.cpp` |
| `hydro/floor_pay_skip` | Pay also the energy that the density floor and the velocity ceiling change. Needs `floor_pay`. | true (general EOS with `<rad_vet>`) | all: `true` | `hydro/dual_energy.cpp` |
| `hydro/floor_from_ekin` | (hydro only.) General-EOS hydro: the primitive recovery keeps E; the velocity ceiling's kinetic energy becomes heat, and a kinetic-energy-dominated cell below e(rho, tfloor) pays the deficit from its own r and theta kinetic energy where that covers it. The phi momentum (L_z) is never scaled (where the phi part alone exceeds `vceil`, the r / theta part goes to 0 and |v| = |v_phi| may stay above `vceil`; Cartesian phi_hat = (-y, x)/R); cells it cannot pay take the floors as before (paid by `floor_pay`). Removes r / theta momentum there (booked: `CLEAN FLOORMOM`). FATAL with another EOS or `dual_energy`. | false (recorded) | not named | `eos/eos.cpp` (constructor) |
| `hydro/dfloor_flux_limit` | (hydro only.) Scale back the outflow of a nearly empty cell so its density stays above the floor (both neighbours see the same flux; a face scaled to zero mass flux keeps only its normal momentum flux, as a closed wall); the cell's momentum but its L_z part is scaled by the same factor, E fixed (the kinetic energy becomes heat); the momentum taken along P (Cartesian: all, spherical polar / cubed sphere: m_r) goes by a face floor flux to the cells the outflow went to (a physical face's share leaves the mesh), so P stays exact. | true with `<rad_vet>` (adiabatic, uniform mesh) | all: `true` | `hydro/hydro.cpp` |
| `hydro/p_ledger` | (hydro only.) The momentum ledger: per accepted step the change of gas + radiation momentum (Cartesian p_x, p_y, p_z on a Cartesian mesh, the radial p_r on spherical polar and cubed sphere), summed exactly, booked by operator class (boundary fluxes `hflx`, non-telescoping flux divergence `ntel`, sources `src` with the reference force the radiation withholds, radiation update `rad`: its boundary stress and geometric terms, the cone cut the gas does not take `cut`, the unprojected fwb residual `fwb`, floors and the floor-gas cap `flr`, the merged lateral cells `mrg`, the rest `gap`; the start-up is not booked); prints `CLEAN PLEDGER` and adds the `P*_` columns to `user.hst` (star_vet, `am_hist`). Reads the state only (bitwise the run without it); about 10 exact mesh sums per stage. The withheld slopes go to the restart file (`p_ledger_k1`), so the ledger restarts bitwise; a file without them: a NOTE. | false (read only when named) | not named | `utils/p_ledger.hpp`, `hydro/hydro.cpp` (constructor) |
| `hydro/fofc`, `fofc_retest` | First-order flux correction: a trial update from the fluxes flags every cell where a floor or the velocity ceiling would act, and its faces take the first-order LLF flux (both neighbours see the same flux: E, mass, momentum and L_z exact). `fofc_retest` re-tests the flagged cells with the metric-correct internal energy (needs `eos = general`). Needs `nghost >= 3` with `plm`. The cell-stages replaced are counted (`FOFC first-order cell-stages`, the run total). | false / false | all: `true` / `true` | `hydro/hydro_fofc.cpp` |
| `hydro/fofc_flag_exchange`, `mhd/fofc_flag_exchange` | Every mesh with more than one MeshBlock (cubed sphere: per panel) and a periodic wrap of spherical polar, no refinement, no shearing box: the FOFC flags of the face ghosts (block edge, periodic wrap) are the owner's, so the two sides of a face carry the same replaced flux (mass and energy exact) (one exchange per stage), so the answer does not depend on the MeshBlock layout where a panel seam meets a MeshBlock edge. | true, also on a restart whose file lacks it (a NOTE in the log; name `fofc_flag_exchange=false` to keep the old) | - | `hydro/hydro_fofc.cpp` |
| `mhd/cs_seam_emf_exchange` | Cubed sphere, GS07 corner EMF (`emf_scheme = gs07`, the only one with `fofc`), more than one MeshBlock per panel: the upwind sides of the corner EMF read the shared seam mass flux in the same-panel ghost rows too (one exchange per stage), so the answer does not depend on the MeshBlock layout (round-off instead of 1e-6 ... 2e-3); conservation and div B are exact either way. | true (a restart file without it: true, with the NOTE) | - | `mhd/mhd_cs_seam.cpp` |

## Layout and reproducibility

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `mesh/layout_invariant` | Global and shell sums in exact, order-independent form: results independent of the MPI rank count. | true | all: `true` | `mesh/mesh.cpp` (constructor) |
| `mesh/layout_invariant_mode` | `tile`: fixed N x N lateral tiles as the summation and multigrid units, so results depend neither on the rank count nor on the layout of whole-tile MeshBlocks (with ray walls in theta: one MeshBlock along x2); `exact`: every cell term exact (about +8 % cost). `block` (older files) is read as `tile`: the same sums. | `tile` | `tile` | `mesh/mesh.cpp` (constructor) |
| `rad_vet/vet_gd_band_ring` | Lend the band intensities per shell from a few slabs: same results, less GPU memory (RUN.md section 6). | on where it applies | A100 sp wedges: `true` | `rad_vet/rad_vet_params_closure.cpp` (constructor) |
| `mesh/layout_tile` | Tile size N for `tile` mode. N dividing `<meshblock>` nx2 and nx3 (FATAL otherwise); a tile smaller than the MeshBlock needs one MeshBlock along x1 (rad_vet). Recorded. Without it (restart of an older file) the MeshBlock is the tile. | fresh start: `16`, so every layout of whole 16 x 16 tiles gives the same answer without naming it; the square lateral MeshBlock where such tiles do not fit (lateral size not a multiple of 16, more than one MeshBlock along x1, SMR/AMR, an odd number of MeshBlocks on a periodic line) | `128` | `mesh/mesh.cpp` (constructor) |
| `mesh/nghost` | Ghost cells; the face line solves across MeshBlocks need >= 3. | 3 on a fresh VET run | all: `3` | `mesh/mesh.cpp` (constructor) |
| `meshblock/nx1`, `nx2`, `nx3` | MeshBlock size. On spherical polar and the cubed sphere nx1 must equal `<mesh>/nx1` (FATAL otherwise); tile mode needs it on every mesh. | = mesh | per deck | `mesh/mesh.cpp` (constructor) |

Changing the mode changes the round-off of the sums: switch at a fresh start, or accept a
round-off change at the restart.

## Output and verbosity

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/verbose` | 0: the `CLEAN cycle` lines, `CLEAN WARNING` lines and event lines (`CLEAN WORST`, `GNSW`, `KPOS` without a face slope); 1: also the per-solve diagnostics (`CLEAN MCAP`, `OLDVEC`, `FWB`, `FCONE`, `LZ`/`LZB`, `WBFCT`, `KPOS`, `KREAL`, `XPAT`, the `XT` write-back closure, the ray cut); 2: also a `CLEAN PASS` line per pass of every solve (Krylov its, the changes de / dv, the tolerance they are tested against, the round-off levels fe / fv, `clos` / `opac` where the pass took a closure / opacity refresh). `CLEAN WARNING` lines always print. | 0 (read only when named) | all: `1` | `rad_vet/rad_vet_params.cpp` (constructor) |
| `rad_vet/floorcap_report` | With `floorcap_rho`: print the cumulative floor-gas cap counters every N cycles (0 = never). | 25 (read only when named) | all: default | `rad_vet/rad_vet_params_matter.cpp` (constructor) |
| `time/ndiag` | Print the cycle line every N cycles. | 1 | all: `50` | `driver/driver.cpp` |
| `hydro/floor_audit` | Print the floor events and their payment every N cycles (`floor_audit cycle` lines). Needs `eos = general`. | 0 | all: `1` | `hydro/dual_energy.cpp` |
| `hydro/floor_audit_worst` | With each `floor_audit` report, every rank prints its cell of the largest energy the inversion created ((E_after - E_before) x vol): its place and its density, total / internal energy, kinetic energy by component and potential before and after (`floor_audit_worst` lines). No effect on the state. Needs `floor_audit`. | false (read only when named) | not named (diagnosis) | `hydro/dual_energy.cpp` |
| `hydro/dt_report` | Print every N cycles which cell and term set dt (`DTREP`, `DTREG`; RUN.md section 5). No effect on the state. | 0 (read only when named) | not named | `hydro/hydro.cpp` |
| `time/nan_check_cycles` | Check the state for NaN every N cycles and stop if found (0 = off). | 100 | all: `100` | `driver/driver.cpp` |

## Removed keys

These keys were removed: the code runs the value given here (the production value, one code path).
**Naming a removed key with that value is accepted; naming another value (deck or command line)
stops the run (FATAL), with the name of the value the code runs.** A value only the restart file
holds (no deck or command line names the key; also a default an older code recorded) prints a
NOTE with the stored value and the value the code runs, and the run continues with the latter.
This is also the only protection
of decks this list was not checked against: the proof that nothing used another value covers the
test suite, the star templates, the RSG / HG / HG262 / R136 continuations and the IA123 decks,
not every deck of a user. A restart file of an older code that does not hold such a key
prints one NOTE with the keys and values the code runs.

| removed key = value | | |
|---|---|---|
| `be_legacy` = `false` | `clean_bc_area` = `true` | `clean_bc` = `halfrange` |
| `clean_bicg_stall_fac` = `0.99` | `clean_bicg_stall_its` = `0` | `clean_conv_est_f` = `false` |
| `clean_cs_lz_abs` = `2` | `clean_cs_lz_bnd_free` = `true` | `clean_face_area` = `false` |
| `clean_face_hctau` = `true` | `clean_face_jse_cell` = `false` | `clean_face_jse` = `true` |
| `clean_face_opac` = `harmonic` | `clean_fill_cs_remap` = `false` | `clean_fill_passes` = `8` |
| `clean_gas_newton_switch_min` = `6` | `clean_gas_newton_switch` = `4` | `clean_gas_secant` = `max_stage2` |
| `clean_gas_tfloor` = `true` | `clean_halo_direct` = `true` | `clean_init_cnorm` = `1e-10` |
| `clean_layout` = `staggered` | `clean_lin_cnorm` = `1e-10` | `clean_lin_pass` = `true` |
| `clean_maxpic` = `12` | `clean_nl_norm` = `global` | `clean_pass_damp` = `0` |
| `clean_pred_stage` = `false` | `clean_pred_xch` = `true` | `clean_predict` = `true` |
| `chat_over_c` = `1` | `clean_ray_cut_margin` = `10` | `clean_ray_cut_taumin` = `0` |
| `clean_rays` = `sc` | `clean_retry_budget` = `true` | `clean_ring_jmin` = `true` |
| `clean_rot_ring_fill_rst` = `true` | `clean_rst_mshare` = `true` | `clean_source_convex` = `true` |
| `clean_stag_central` = `true` | `clean_stag_colflux` = `false` | `clean_stag_csvec_radial` = `true` |
| `clean_stag_eicarry` = `true` | `clean_stag_flim` = `false` | `clean_stag_fs_line` = `pcr` |
| `clean_stag_fs_pcr` = `true` | `clean_stag_fs_pcr_cs` = `true` | `clean_stag_fwb_cart` = `false` |
| `clean_stag_fwb_stage` = `false` | `clean_stag_fwb_t2` = `false` | `clean_stag_lat_sums` = `true` |
| `clean_stag_mg_color` = `mirror` | `clean_stag_opac_drag` = `true` | `clean_stag_opac_update` = `true` |
| `clean_stag_pcr` = `true` | `clean_stag_rayratio` = `false` | `clean_stag_tie_gid` = `true` |
| `clean_stag_vimp` = `true` | `clean_stag_vimp_cs` = `true` | `clean_stag_xdonor` = `lag` |
| `clean_stag_xsrc` = `true` | `clean_stall_fac` = `0.5` | `clean_troot_rtol` = `1e-13` |
| `clean_vel_settled` = `false` | `clean_w_face_min` = `false` | `clean_w_fse` = `true` |
| `clean_w_fse_tau` = `true` | `clean_wb_fct` = `true` | `clean_wb_floor_fail` = `false` |
| `clean_xt_cons` = `false` | `clean_xt_exact` = `cons` | `clean_xt_norm` = `false` |
| `dbg_energy_tally` = `false` | `doppler_face` = `donor` | `floorcap_mesh_vol` = `true` |
| `implicit_am_cs` = `true` | `implicit_enthalpy` = `plm` | `implicit_pos_floor_cons` = `true` |
| `implicit_pos_floor_gas_keep` = `1` | `implicit_pos_gas` = `true` | `implicit_pos_gas_cons` = `true` |
| `implicit_pos_gas_frac` = `1e-3` | `implicit_predictor` = `step` | `implicit_predictor_tcap` = `1.1` |
| `polar_ring_avg` = `false` | `restart_pack_per_mb` = `true` | `rsla_force` = `false` |
| `rsla_vmax` = `-1` | `rsla_warn` = `0.1` | `scaled_prefactor` = `20` |
| `scheme` = `clean` | `time2_be_only` = `false` | `time2_cfl_post` = `true` |
| `time2_face_fct` = `true` | `time2_fref_wacc` = `true` | `time2_gtime_mom` = `true` |
| `time2_k1_redist` = `0` | `time2_k1_stiff` = `0` | `vet_axis` = `flux` |
| `vet_col_lat_cut_thin` = `1` | `vet_gd_halo_compact` = `2` | `vet_gd_rotate_align` = `true` |
| `vet_gd_rotate_stretch` = `0.02` | `vet_gd_rotate_t0` = `0` | `vet_gd_rotate_xfade_m` = `0` |
| `vet_gd_twin` = `false` | `vet_interp` = `auto` | `vet_rotate_avg_ema_change` = `0.05` |
| `vet_rotate_avg_ema_change_j` = `true` | `vet_rotate_avg_ema_dema` = `false` | `vet_rotate_avg_ema_diag` = `0` |
| `vet_rotate_avg_ema_lc` = `0` | `vet_rotate_avg_ema_tau` = `0` | `vet_rotate_avg_weight` = `j` |
| `vet_scatter_form` = `ma` | `vet_scatter_j` = `relaxed` | `vet_seg` = `jump` |
| `vet_seg_plateau` = `false` | `vet_source_noesrc_h2be` = `true` | `vet_tensor` = `full` |
| `zalesak_active` = `true` | `zalesak_halo_skip` = `true` |  |

Removed keys whose value follows from the run (FATAL when named otherwise):

| key | the rule |
|---|---|
| `clean_stag_fwb`, `clean_stag_freal` | as `closure_consistency` |
| `clean_stag_fschur` | true where the mesh allows it (nghost >= 3, MeshBlocks long enough), else the cubed sphere iterates its seams and any other multi-block line is FATAL |
| `clean_stag_fsweeps` | 1; 64 on a cubed sphere without fschur |
| `clean_stag_fwb_ring` | true on 3-D spherical polar with a periodic phi and the rays |
| `vet_gd_rotate_axis` | `3d` on the cubed sphere, `z` on spherical polar |
| `vet_gd_rotate_xfade` | true with a mean of sets rotating in time, not with `ring = jump`; on a Cartesian mesh with `ring = live` only |
| `vet_rotate_advance` | true on the curved meshes and with `ring`, else false |
| `vet_rotate_avg_ema_wb` | true on the curved meshes, else false |
| `vet_gd_band_global` | true on a full sphere |
| `vet_gd_rotate_slope`, `vet_gd_rotate_slope_frac` | `correct`, 1 - 1/sqrt 2 for a Cartesian rotating set without the ring (no effect elsewhere) |
| `implicit_predictor_order` | 2; an older restart file's value is read as the layout of its records |

Removed without effect (named: listed at the start, any value): `advect_split`, `implicit_gas_src`,
`implicit_cs2_vimp_npass`, `implicit_vimp*`, `implicit_recon`, `implicit_solver`, `implicit_offdiag`,
`implicit_trans_limit`, `clean_stag_upwind_central`, `time2_enth_be`, `time2_accept_fallback`, `time2_adm_poslim`,
`time2_one_pass_safety`, `time2_pred_cap`, `time2_vet*`, `cs_seam_avg`, and `vet_gd_twin_*`.
Removed choices: `vet_sc_kernel` (FATAL when named).

## Appendix: advanced and internal keys

These keys hold the scheme (the templates run the defaults) or are solver
internals. Do not change them without a dedicated test (energy and momentum conservation, a
bitwise restart check, and a comparison with the templates).

Scheme options (in brackets the templates' values):

- `closure_consistency` (true with the rays), `face_thin` (true), `flux_cone` (true): the flux bound by the rays, the thin-limit correction, the cone limit |F| <= cE.
- `flux_balance_lat` (default true on the cubed sphere, also on a restart whose file lacks it, with a NOTE in the log: name `flux_balance_lat=false` to keep the old; cubed sphere only, else FATAL): a spherically symmetric radiation pressure pushes the gas along r alone, the same in every cell of a shell. Lateral faces: their cells' lateral stress sums are projected perpendicular to each cell's radial direction (without it a constant lateral force of O(dtheta^2) of the radial one, at every radius, with the cube's symmetry: the interior m = 4 pattern). Radial faces: the normal stress takes the hydro's area r^2 Omega instead of the vector area r^2 |S| (without it a smaller panel pattern of the radial force). Test `shell_force_cs`.
- `flux_cone_gas` (false): `true` gives the gas the momentum the cone cut takes from each inner face (its absorption of the unscaled F), so gas + radiation momentum and L_z are exact, and E pays the kinetic energy; where the radiation cannot pay it the run fails (the radiation-poor exchange box `exch0` of test `floor_momentum`: dt collapse at cycle 1). `false`: the cut is withheld from the gas (`hydro/p_ledger` class `cut`). Old name `clean_stag_fcone_gas`.
- `flux_cone_line` (false; FATAL on a refined mesh; off with `flux_cone_gas`): `true` makes the cut exact: the cone cut each face withholds from the gas less its mean along each mesh line of the directions whose faces carry P, weighted by the faces' gas density (in the write-back's solve): every line's gas takes its cut where its mass is, a telescoping face stress, so the cut moves no P (`hydro/p_ledger` class `cut` 0, tiny BSG: all unphysical classes round-off). Each cell's gas internal energy pays the kinetic energy of the cut it takes back (E unchanged by the return), capped at half its internal energy above e_min; a capped cell is counted (`CLEAN FCLINE`: cells, momentum not returned, which `p_ledger` books as `cut`). Not the default: test `floor_momentum`'s `exch` box (arad 1e-12) caps every cell and still takes two backward-Euler redos with it (`CLEAN WARNING`).
- `fwb_line` (true on a uniform mesh, FATAL on a refined one): the well-balancing residual s less its mean along each mesh line of the directions whose faces carry the momentum P (x1 on spherical polar and cubed sphere; x1, x2, x3 on a Cartesian mesh, unless the old `clean_stag_fwb_cart = false` is named), so fwb gives no P (`hydro/p_ledger` class `fwb` 0). Cost: the time error of F in `order_time` +18 % (finest 1.98e-4 instead of 1.68e-4, order 1.79). `false`: the release's unprojected residual. Old name `clean_stag_fwb_line`.
- `flux_cone_lz` (read only when named): fcone takes a cell's excess from the face families without an L_z arm first (spherical polar x1, x2; cubed sphere x1; Cartesian x3); that cut is counted and printed signed per family (`CLEAN FCONE_LZ`). The faces with an arm are scaled only where they alone exceed cE, and the gas takes their cut at the face, so L_z is exact.
- `lz_spin` (read only when named; needs `flux_cone_lz`; a no-op on the cubed sphere, where `lz_telescope` covers it): the L_z X of the arm faces' fcone cut, the `flux_cone_emax` clamp and the `update_mcap` cap goes to the gas as one rigid rotation dOmega = X/I, I = sum rho V R^2 (no per-cell kick; E pays the kinetic energy; `CLEAN LZSPIN`).
- `flux_cone_emax` (read only when named; needs `flux_cone`): every inner face's F within c max(E_L, E_R), an opaque face's too. L_z-exact with `flux_cone_lz` + `lz_spin`.
- `update_mcap` (read only when named): the write-back gives the gas the face momentum terms at the Newton passes' cap h chi_f E; the withheld momentum is printed per family (`CLEAN MCAP write-back`). L_z-exact only with `lz_spin`.
- `update_exchange_local` (true): the write-back's cell exchange is solved exactly from the transported E (E > 0 and gas >= e(rho, tfloor) by construction); needs an EOS `tfloor` > 0.
- `lz_telescope` (true on the cubed sphere): L_z by telescoping of the face momenta.
- `face_xdonor_freeze` (4): the pass from which the donor form of the off-diagonal stress is held.
- `closure_w_c0` (0.3), `closure_w_c1` (0.6): thresholds of the weight of the kinetic share.
- `gas_opac_passes` (3): passes with opacities at the iterate's temperature.
- `implicit_pos_floor` (true): the E floor, paid from the cell's gas and its radial column. A non-finite E (a diverged iterate) is never raised: the write-back rejects the solve.
- `time2_gtime_slope`, `time2_gtime_pos` (true); `vet_rotate_avg_ema_change_form` (`linear`; HG: `step`); `vet_source_noesrc` (false on every mesh; Cartesian runs took true before).
- `vet_bc_inner`, `vet_bc_diff`, `vet_bc_bath`, `vet_gd_wall_interp`, `vet_gd_allow_clamp`, `closure_ray_min`, `ray_cut_wait`, `ring_fuse`: ray-sweep details.
- `ray_cut_r` (0.003): the ray cut's trigger, a column's lowest cell with |grad S|/(chi S) >= this (SCHEME.md, where the rays start; the closure error below the cut is of the order of its square). 0.01 before 2026-11-26: a restart file that holds 0.01 keeps it, with a NOTE (name `ray_cut_r=0.003` to switch); one without the key takes 0.003, with a NOTE.

Solver internals (in brackets the templates' values): `newton_conv` (`budget`), `newton_conv_max`
(1e-5), `newton_conv_kappa`, `solver_tol`, `newton_tol`, `solver_maxit` (200), `newton_maxit` (200),
`solver_sweep` (default `rb`; `sym` adds a post-sweep), `solver_mg_levels` (default 2, the levels of the
staggered multigrid preconditioner), `newton_closure_budget` (true), `newton_maxnorm` (true)
(false: the closure pass solved in the strict norm; the budget test in the L2 norm),
`implicit_resid_fatal`, `time/dt_retry_views`, `time/dt_retry_host`, `time/restart_refill_ghosts`
(false: restarts are bitwise).

Tests and diagnostics: `dbg_inf`, `time2_dbg_fail`, `dbg_gas_force`, `update_floor_diag`,
`implicit_timers`. State kept in restart files: `rst_version`, `ray_cut_now`,
`ray_cut_count`. Number of `<rad_vet>` keys the code reads: 130 (user 71, physics 26,
scheme and solver 20, tests 8, pgen / state 5).

Removed study and debug keys stop the run with a message when named with another value than
their old default (`CheckKeys` in `rad_vet/rad_vet_params.cpp` lists them); keys of the
removed solvers that old decks and restart files carry are accepted without effect (listed at
the start): also `vet_col_*` (`vet_col_axis` only `radial`, `vet_col_every` only 1), `vet_every` (only 1), `vet_source`, `vet_gd_replace`,
`vet_gd_tangential`, `vet_gd_rotate_slope_reuse`, `time2_vet_source*` and `cs_vertex_mono*`.
`vet_gd` is read only where `closure_form` is `eddington` or `sphere` (true: the rays' arrays
are set up).

Keys that are never read are listed by the `never read` warning after the first cycle.
