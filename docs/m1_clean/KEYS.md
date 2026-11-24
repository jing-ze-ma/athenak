# Input keys of the clean M1 radiation scheme

This page lists the keys a user sets for a run: the scheme, the rays, the time step, the
boundaries, opacities and units, the layout mode, and output. The production values are those
of the star templates in `inputs/rad_m1_clean/` (`star_bsg`, `star_m35mms`, `star_t9l62`:
3-D spherical-polar wedges, one MeshBlock along x1, fixed rays, problem generator
`he_star_m1`). Solver internals and keys without effect are listed only by name in the
appendix.

Columns: **default** is the value the code uses when the key is absent; **production** is
the value in the templates (`all:` when the three agree); **read in** is the file and
function that read the key (paths relative to `src/`).

Many defaults differ between a **fresh start** and a **restart** whose file does not hold
the key (the restart default keeps an old run unchanged). The start-up log prints every default
it applied together with the reason, e.g. `(default with vet_gd on a fresh start)`.

## The restart rule: name every production key

A restart is started with `athena -r <file>.rst [-i <deck>] [block/key=value ...]`. The
parameters are put together in this order (`main.cpp`):

1. the input block **stored in the restart file**;
2. the `-i` deck, if given: every key it names replaces the stored value;
3. the command line: every `block/key=value` replaces the value so far. On a restart a key
   that is not in the stored block is added (the code prints a `### NOTE`). On a **fresh**
   start a command-line key that is not in the deck is FATAL (typo guard), so put every key
   of a fresh run in the deck.

**On a restart, a key that is not given in the deck or on the command line takes the value
stored in the restart file, not the code default.** If the file does not hold the key either,
the key takes its restart default.

Example: a restart file had stored `time2_k1_stiff = 1.0` from an earlier deck. The code
default is 0 (off), but the restart ran with 1.0, because the new deck did not name the key.

Rules:

- When you restart from a restart file written with other keys, restart with
  `-i <production deck>` or name every production key on the command line. The grid keys
  (`<mesh>` x1min..x3max, nx1..nx3, nghost, the stretch keys and `<meshblock>` nx1..nx3)
  must stay as the file was written: a change is FATAL.
- `athena -r <file>.rst -i <deck> -n` prints the merged parameter list and quits (keys the
  modules add later are not shown).
- After the first cycle the code prints `### WARNING: N input parameter(s) were never read by
  the code`: a key there is a typo or a key this build does not have.

### Keys that must always be named

Restart files of earlier runs may carry other values for these keys, so leaving one of them
out of a restart deck changes the scheme silently. The templates name them all.

| key | value | why |
|---|---|---|
| `<rad_m1>/time2_k1_stiff` | `0` | read only when named; older files stored 1.0 |
| `<rad_m1>/clean_stag_fwb_stage` | `false` | older files stored true |
| `<rad_m1>/clean_stag_fwb_t2` | `false` | older files stored true |
| `<rad_m1>/clean_w_fse_tau` | `true` | older files may store false |
| `<rad_m1>/time2_gtime_mom` | `true` | older files may store false |
| `<rad_m1>/clean_pred_stage` | `false` | production value |
| `<rad_m1>/clean_stag_central` | `true` | production value (default false) |
| `<rad_m1>/clean_stag_kpat` | `false` | production value |
| `<rad_m1>/clean_stag_eicarry` | `true` | production value (default false) |
| `<mesh>/layout_invariant_mode` | `tile` | see "Layout and reproducibility" |
| `<mesh>/layout_tile` | `128` | the lateral MeshBlock size of the templates |

## Scheme selection

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_m1/transport` | Radiation transport; the clean scheme needs `implicit`. | `implicit` on a fresh run with a `<rad_m1>` block | all: `implicit` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/scheme` | Radiation solver. Only `clean` exists; any other value is FATAL. | `clean` with implicit transport | all: `clean` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/closure` | Closure. `vet_col` with `vet_gd = true`: variable Eddington tensor from the global ray sweep (required on spherical polar and the cubed sphere). | `vet_col` with `vet_gd` on a fresh sp/cs run; `m1` otherwise | all: `vet_col` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_gd` | Use the global ray sweep for the Eddington tensor. | true with `closure = vet_col` on a fresh sp/cs run | all: `true` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_gd_replace` | The global sweep replaces the per-column sweep. | false (read only when named) | all: `true` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/clean_closure` | How the closure follows the solution in time: `extrap1` (one ray sweep per step; recommended), `iterate`, `predict`, `eddington` (D = I/3, no rays; for tests). | `extrap1` | all: `extrap1` | `rad_m1/rad_m1.cpp` (constructor) |

## Rays

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_m1/vet_gd_gl_nmu`, `vet_gd_gl_nphi`, `vet_gd_gl_stagger` | Product direction set: nmu Gauss-Legendre polar nodes times nphi azimuths, alternate rings shifted with `stagger`. 2 pi/nphi must divide the phi extent (FATAL otherwise). Takes precedence over the other sets. | no set named: `vet_gd_geo = 3` (92 directions) on a fresh 3-D run | all: 7, 18, `true` (126 directions) | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_gd_ls` | Level-symmetric set LQ_N, N = 6, 8, 10, 12 (N = 10: 120 directions). | 0 (off) | all: `10` (no effect: the product set is named) | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_gd_geo` | Geodesic set, 10 n^2 + 2 directions. | 3 when no set is named | not named | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_gd_rotate_dt` | Rotate the direction set every dt_rot of physical time (`auto` = 20 x the first dt; 0 = fixed set). | `auto` on a fresh run that names no rotation key | all: `0` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_rotate_avg` | With rotation: average the closure over the last K sets (0 = off). | 4 with a rotating set on a fresh run | all: `0` | `rad_m1/rad_m1.cpp` (`RotAvgParse`) |
| `rad_m1/clean_rot_ring` | Rotating ray ring of the clean scheme: `off`, `live`, `jump`. | read only when named | all: `off` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_scatter`, `vet_scatter_form` | Split the opacity into absorption and electron scattering in the ray source (`ma`: Ma et al. 2026, eq. 7). Changes only the rays' source function. | true on a fresh run with table opacity and a tabulated EOS; `ma` | all: `true`, `ma` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/vet_scatter_kappa_e` | Electron-scattering opacity: `const` or `eos` (from the EOS ionisation). | `eos` with a tabulated EOS; `const` otherwise | all: `const` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/clean_ray_cut_r` | Adaptive ray depth cut on 3-D sp/cs meshes with one MeshBlock along x1: below the cut the closure is the diffusion limit. | 0.01 where the cut applies | all: `0.01` | `rad_m1/rad_m1.cpp` (constructor) |

## Time integration, time step and retries

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `time/integrator` | Hydro integrator; must be `rk2` (the explicit part of the IMEX scheme). | `rk2` | all: `rk2` | `driver/driver.cpp` (constructor) |
| `time/cfl_number` | Courant number of the hydro step; the radiation does not limit dt. | 0.3 on a fresh clean M1 run | all: `0.3` | `mesh/build_tree.cpp` |
| `rad_m1/time_scheme` | `hesdirk2`: two implicit stages per step with the Heun hydro (second order); `be`: one backward-Euler solve per step (first order). | `hesdirk2` on a fresh run that supports it; `be` otherwise | all: `hesdirk2` | `rad_m1/rad_m1_init.cpp` (`ImplicitInit`) |
| `rad_m1/time2_tableau` | Tableau of the two stages: `sdirk2` (Alexander's L-stable SDIRK2, recommended), `hesdirk2`, `trbdf2`. | `hesdirk2` | all: `sdirk2` | `rad_m1/rad_m1_time2.cpp` (`Time2Init`) |
| `rad_m1/time2_auto_be` | After a step whose stages fell back to backward Euler, the next 1, 2, 4, ... N steps take backward Euler at once (0 = off). | 16 | all: `16` | `rad_m1/rad_m1_time2.cpp` (`Time2Init`) |
| `rad_m1/clean_init` | `steady`: before step 1 the radiation is relaxed to the scheme's steady solution in the given gas; `given`: start from the given field. | `given` | all: `steady` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/clean_init_dtfac`, `clean_init_tol`, `clean_init_n` | Pseudo time step (in units of the first dt), tolerance and pass limit of the steady start. | 1e6, 1e-8, 40 or 50 (problem dependent) | all: `1.0e6`, not named, not named | `rad_m1/rad_m1.cpp` (constructor) |
| `time/dt_retry` | A cycle whose implicit solve fails (not converged, non-finite, rejected stage, unpaid floor energy) is undone and taken again with a smaller dt. Bitwise: a retry restores the full cycle. | true | all: `true` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_fac`, `dt_retry_fac_t2` | dt factor of a retry; of a retry caused by a rejected stage. | 0.5; = `dt_retry_fac` | all: `0.5`, `0.9` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_grow`, `dt_retry_ngood` | After a cut, dt grows back by this factor per good cycle, starting after this many good cycles. | 1.25, 10 | all: `1.1`, `20` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_t2_be` | After this many stage-rejection cuts in one cycle, a rejected stage is redone with backward Euler (0 = at once, < 0 = never). | 0 | all: `0` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_force_cycle`, `dt_retry_force_fac` | Testing: force one retry at this cycle with this factor (the `dt_retry` test). | off | not named | `driver/dt_retry.cpp` (`DtRetry`) |

## Boundaries

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_m1/implicit_bc_x1min`, `implicit_bc_x1max` | Radiation x1 boundaries: `marshak`, `outflow`, `reflect`, `periodic`, and on x1 an imposed `flux`. The ray sweep needs `x1max = marshak` or `reflect`. | `auto` (from the mesh boundary) | all: `flux`, `marshak` | `rad_m1/rad_m1_init.cpp` (`ImplicitInit`) |
| `rad_m1/implicit_flux_x1min` | Flux through the inner face (= L/(4 pi r_in^2) for a star). | 0 | per deck | `rad_m1/rad_m1_init.cpp` (`ImplicitInit`) |
| `rad_m1/implicit_flux_x1max` | Incoming bath at the outer Marshak face (0 = vacuum). | 0 | all: `0.0` | `rad_m1/rad_m1_init.cpp` (`ImplicitInit`) |
| `rad_m1/clean_bc` | Open x1 end: `halfrange` (the end cell's ray half-range flux) or `marshak`. | `halfrange` | all: `halfrange` | `rad_m1/rad_m1.cpp` (constructor) |
| `mesh/ix1_bc` ... `ox3_bc` | Hydro boundaries (the templates: `user` in x1, `periodic` in theta and phi). | required | see the decks | `mesh/mesh.cpp` |

Transverse Marshak faces are not supported by the ray sweep (FATAL).

## Opacity, units and coupling

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_m1/opacity` | `const`, `powerlaw`, `user` or `table` (Rosseland and Planck tables named in `<problem>`). | `const` | all: `table` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/kappa_p`, `kappa_e`, `kappa_f`, `kappa_s` | Constant opacities per unit mass: Planck (emission), energy mean (absorption), flux mean, scattering. | 0, = kappa_p, = kappa_p, 0 | not named | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/c_light`, `arad` | Speed of light and radiation constant in code units. The start-up log compares them with the values `<units>` implies. | required | per deck | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/chat_over_c` | Reduced speed of light; the clean scheme needs 1 (FATAL otherwise). | 1.0 | all: `1.0` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/temp_unit_kelvin`, `rho_unit_cgs`, `kappa_unit` | Unit conversions for the table lookup: K per code temperature, g/cm^3 per code density, code opacity per cm^2/g. | 1.0 | per deck | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/coupling`, `gas_feedback` | Gas-radiation exchange on; exchange acts back on the gas. | true | all: `true` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/force_reference`, `force_reference_work` | `wb_arad`: well-balanced radiative envelope; the problem generator supplies the reference radiative acceleration, the gas receives only the residual force; its work is paid as `split` (required by the clean scheme). | `none`; `auto` | all: `wb_arad`, `split` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/e_floor` | Floor of E (code units). | FLT_MIN | all: `1.0e-4` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/floorcap_rho`, `floorcap_v` | Floor-gas cap: in cells with rho below `floorcap_rho`, |v| is capped at `floorcap_v` and T at max(T_rad, floor); the removed energy goes to E (momentum is removed and reported). | off | all: `1e-14`, `3e7` | `rad_m1/rad_m1.cpp` (constructor) |
| `hydro/dfloor`, `pfloor`, `tfloor` | Hydro floors. `hydro/floor_pay = rad_column` pays the internal energy they create from the radiation. | FLT_MIN | per deck | `eos/eos.cpp`, `hydro/dual_energy.cpp` |

## Layout and reproducibility

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `mesh/layout_invariant` | Global and shell sums in exact, order-independent form: results independent of the MPI rank count. | true | all: `true` | `mesh/mesh.cpp` (constructor) |
| `mesh/layout_invariant_mode` | `block`: sums per MeshBlock (independent of the rank count for fixed MeshBlocks); `exact`: every cell term exact (about +8 % cost); `tile`: fixed N x N lateral tiles as the summation and multigrid units, so results do not depend on the MeshBlock layout either. | `block` | `tile` (recommended) | `mesh/mesh.cpp` (constructor) |
| `mesh/layout_tile` | Tile size N for `tile` mode; required with it. Needs one MeshBlock along x1 and N dividing `<meshblock>` nx2 and nx3. | none | `128` | `mesh/mesh.cpp` (constructor) |
| `mesh/nghost` | Ghost cells; the face line solves across MeshBlocks need >= 3. | 2 | all: `3` | `mesh/mesh.cpp` (constructor) |
| `meshblock/nx1`, `nx2`, `nx3` | MeshBlock size. One MeshBlock along x1 is needed for the ray cut and tile mode. | = mesh | per deck | `mesh/mesh.cpp` (constructor) |

Changing the mode changes the round-off of the sums: switch at a fresh start, or accept a
round-off change at the restart.

## Output and verbosity

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_m1/clean_verbose` | 0: the `CLEAN cycle` lines, `CLEAN WARNING` lines and event lines (`CLEAN WORST`, `GNSW`, `BICGSTAB`, `RETRY`, `KPOS` without a face slope); 1: also the per-solve diagnostics (`CLEAN MCAP`, `OLDVEC`, `FWB`, `FLIM`, `FCONE`, `LZ`/`LZB`, `WBFCT`, `KPOS`, `KREAL`, `XPAT`, the `XT` write-back closure, the ray cut). `CLEAN WARNING` lines always print. | 0 (read only when named) | all: `1` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/floorcap_report` | Print the cumulative floor-gas cap counters every N cycles. | 0 = never | all: `25` | `rad_m1/rad_m1.cpp` (constructor) |
| `rad_m1/clean_nl_print` | Pass-by-pass solver diagnostics (debugging). | false | not named | `rad_m1/rad_m1.cpp` (constructor) |
| `time/ndiag` | Print the cycle line every N cycles. | 1 | all: `50` | `driver/driver.cpp` |
| `time/nan_check_cycles` | Check the state for NaN every N cycles and stop if found (0 = off). | 100 | all: `100` | `driver/driver.cpp` |

## Appendix: advanced and internal keys

These keys are named in the templates to pin the scheme (keep their values) or are solver
internals. Do not change them without a dedicated test (energy and momentum conservation, a
bitwise restart check, and a comparison with the production key set).

Scheme options pinned by the templates:

- `clean_stag_upwind_central` (true), `clean_stag_central` (true): central advective flux within the TVD bound.
- `clean_face_hctau` (true), `clean_face_jse` (true), `clean_face_opac` (`harmonic`): face values of the rays' flux, half ranges and opacity.
- `clean_stag_thinf` (false), `clean_stag_flim` (false), `clean_consistency` (true): forms of the thin-limit correction and flux bound.
- `clean_stag_fwb` (true), `clean_stag_fcone` (true), `clean_stag_fwb_ring` (0): well-balanced face rows, cone limit |F| <= cE, ring projection.
- `clean_stag_xdonor` (`lag`), `clean_stag_xsrc` (true), `clean_stag_xdonor_freeze` (4): donor form of the off-diagonal stress.
- `clean_w_fse`, `clean_w_c0`, `clean_w_c1`: weight of the kinetic share.
- `clean_xt_exact` (`cons`), `clean_wb_fct` (true), `time2_face_fct` (true): exact exchange bookkeeping and the Zalesak limiters.
- `clean_stag_opac_update` (true), `clean_stag_opac_passes` (3): opacities at the iterate's temperature.
- `clean_stag_lat_sums`, `clean_stag_opac_drag`, `clean_stag_tie_gid` (all 0): fresh-run defaults that the templates switch off.
- `implicit_pos_floor`, `implicit_pos_floor_cons`, `implicit_pos_gas`, `implicit_pos_gas_cons` (true): E floor and gas positivity paid from the radiation.
- `vet_seg` (`jump`), `vet_interp` (`auto`), `vet_bc_inner`, `vet_bc_diff`, `vet_gd_tangential`, `vet_gd_wall_interp`, `vet_gd_halo_compact`, `vet_gd_allow_clamp`: ray-sweep details.

Solver internals: `clean_conv` (`budget`), `clean_conv_max` (1e-5), `clean_conv_kappa`,
`clean_tol`, `clean_tol_nl`, `clean_maxit` (200), `clean_maxnl` (200), `clean_maxpic`,
`clean_lin_pass`, `clean_lin_cnorm`, `clean_init_cnorm`, `clean_gas_newton_switch`,
`clean_gas_newton_switch_min`, `clean_stall_fac`, `clean_troot_rtol`, `clean_pc_sweep` (`rb`),
`clean_stag_mg_color`, `clean_stag_fs_line`, `clean_stag_pcr`, `clean_halo_direct`,
`implicit_predictor`, `implicit_predictor_tcap`, `implicit_resid_fatal`, `time/dt_retry_views`,
`time/dt_retry_host`, `time/restart_refill_ghosts` (false: restarts are bitwise).

Not recommended: `clean_pass_damp` (keep 0), `clean_stag_kpat`, `time2_k1_redist`,
`time2_k1_stiff` (keep 0).

Named in the templates but without effect with them (fixed rays, spherical polar): the
rotation keys other than `vet_gd_rotate_dt`, the cubed-sphere keys (`clean_cs_lz_tel`,
`cs_cross_mono`, `cs_vertex_mono`, `mesh/cs_*`), `vet_source`, `vet_scatter_j` and the other
relaxed-source keys (unless `clean_ray_source_relaxed = true`), `vet_col_lat_taucut`,
`implicit_maxit`, `implicit_lin_maxit`, `time2_enth_be`, `time/dt_retry_max`. Keys that are
never read are listed by the `never read` warning after the first cycle.
