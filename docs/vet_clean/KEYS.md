# Input keys of the clean VET radiation scheme

This page lists the keys a user sets for a run: the scheme, the rays, the time step, the
boundaries, opacities and units, the layout mode, and output. The production values are those
of the star templates in `inputs/rad_vet_clean/` (`star_bsg`, `star_m35mms`, `star_t9l62`:
3-D spherical-polar wedges, one MeshBlock along x1, fixed rays, problem generator
`he_star_vet`). Solver internals and keys without effect are listed only by name in the
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
| `<rad_vet>/time2_k1_stiff` | `0` | read only when named; older files stored 1.0 |
| `<rad_vet>/clean_stag_fwb_stage` | `false` | older files stored true |
| `<rad_vet>/clean_stag_fwb_t2` | `false` | older files stored true |
| `<rad_vet>/clean_w_fse_tau` | `true` | older files may store false |
| `<rad_vet>/time2_gtime_mom` | `true` | older files may store false |
| `<rad_vet>/clean_pred_stage` | `false` | production value |
| `<rad_vet>/clean_stag_central` | `true` | production value (default false) |
| `<rad_vet>/clean_stag_kpat` | `false` | removed (only false) |
| `<rad_vet>/clean_stag_eicarry` | `true` | production value (default false) |
| `<rad_vet>/clean_gas_secant` | `max_stage2` | production value (default false) |
| `<rad_vet>/clean_pred_xch` | `true` | production value (default false) |
| `<rad_vet>/clean_stag_fcone_gas` | `false` | production value |
| `<rad_vet>/clean_wb_exch_local` | `true` | read only when named |
| `<rad_vet>/vet_gd_wall_sym`, `vet_gd_theta_seam_sym` | `true` | read only when named; a theta-periodic deck without `vet_gd_theta_seam_sym` stops with a FATAL |
| `<hydro>/floor_pay` | `rad_column` | older star restart files do not store it (restart default: none) |
| `<hydro>/floor_pay_skip`, `dfloor_flux_limit`, `rad_signal_speed_diff` | `true` | read only when named |
| `<mesh>/layout_invariant_mode` | `tile` | see "Layout and reproducibility" |
| `<mesh>/layout_tile` | `128` | the lateral MeshBlock size of the templates |

## Scheme selection

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/transport` | Radiation transport; the clean scheme needs `implicit`. | `implicit` on a fresh run with a `<rad_vet>` block | all: `implicit` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/scheme` | Radiation solver. Only `clean` exists; any other value is FATAL. | `clean` with implicit transport | all: `clean` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/closure` | Closure. `vet_col` with `vet_gd = true`: variable Eddington tensor from the global ray sweep (required on spherical polar and the cubed sphere). | `vet_col` with `vet_gd` on a fresh sp/cs run; `m1` otherwise | all: `vet_col` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_gd` | Use the global ray sweep for the Eddington tensor. | true with `closure = vet_col` on a fresh sp/cs run | all: `true` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_gd_replace` | The global sweep replaces the per-column sweep. | false (read only when named) | all: `true` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/clean_closure` | How the closure follows the solution in time: `extrap1` (one ray sweep per step; recommended), `iterate`, `predict`, `eddington` (D = I/3, no rays; for tests). | `extrap1` | all: `extrap1` | `rad_vet/rad_vet.cpp` (constructor) |

## Rays

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/vet_gd_gl_nmu`, `vet_gd_gl_nphi`, `vet_gd_gl_stagger` | Product direction set: nmu Gauss-Legendre polar nodes times nphi azimuths, alternate rings shifted with `stagger`. 2 pi/nphi must divide the phi extent (FATAL otherwise). Takes precedence over the other sets. | no set named: `vet_gd_geo = 3` (92 directions) on a fresh 3-D run | all: 7, 18, `true` (126 directions) | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_gd_ls` | Level-symmetric set LQ_N, N = 6, 8, 10, 12 (N = 10: 120 directions). | 0 (off) | all: `10` (no effect: the product set is named) | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_gd_geo` | Geodesic set, 10 n^2 + 2 directions. | 3 when no set is named | not named | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_gd_rotate_dt` | Rotate the direction set every dt_rot of physical time (`auto` = 20 x the first dt; 0 = fixed set). | `auto` on a fresh run that names no rotation key | all: `0` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_rotate_avg` | With rotation: average the closure over the last K sets (0 = off). | 4 with a rotating set on a fresh run | all: `0` | `rad_vet/rad_vet.cpp` (`RotAvgParse`) |
| `rad_vet/clean_rot_ring` | Rotating ray ring of the clean scheme: `off`, `live`, `jump`. | read only when named | all: `off` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_scatter`, `vet_scatter_form` | Split the opacity into absorption and electron scattering in the ray source (`ma`: Ma et al. 2026, eq. 7). Changes only the rays' source function. | true on a fresh run with table opacity and a tabulated EOS; `ma` | all: `true`, `ma` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_scatter_kappa_e` | Electron-scattering opacity: `const` or `eos` (from the EOS ionisation). | `eos` with a tabulated EOS; `const` otherwise | all: `const` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/clean_ray_cut_r` | Adaptive ray depth cut on 3-D sp/cs meshes with one MeshBlock along x1: below the cut the closure is the diffusion limit. | 0.01 where the cut applies | all: `0.01` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/vet_gd_wall_sym` | Phi-periodic wedges: centre the set on the wedge, mirror-paired phi-wall maps, the seam face from both sides. Needs `vet_gd_wall_interp = true`. Identical results on wedges that are whole azimuthal periods; required on narrower ones. No effect without a phi wall (cubed sphere, full 2 pi, 2-D). | false (read only when named) | all: `true` | `rad_vet/rad_vet.cpp` (constructor), `rad_vet/rad_vet_vetgd.cpp` |
| `rad_vet/vet_gd_theta_seam_sym` | Theta-periodic wedges: the theta seam face's other side is its exact equatorial mirror (both frames, as the phi seam). Required on a theta-periodic mesh (FATAL otherwise); the ray set must be mirror-symmetric about the equator (FATAL otherwise; e.g. GL with even nmu and stagger is not). No effect without a theta seam. | false (read only when named) | all: `true` | `rad_vet/rad_vet.cpp` (constructor), `rad_vet/rad_vet_vetgd.cpp` |

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
| `time/cfl_number` | Courant number of the hydro step; the radiation does not limit dt. | 0.3 on a fresh clean VET run | all: `0.3` | `mesh/build_tree.cpp` |
| `hydro/rad_signal_speed_diff` | With `rad_signal_speed = true`: the radiation term of the signal speed weighted by 1 - exp(-x), x = tau_d c_eq / c the grid-scale trapping parameter; untrapped radiation (x << 1) does not limit dt. Time step only (accuracy, not stability). | false (read only when named) | all: `true` | `hydro/hydro.cpp` |
| `rad_vet/time_scheme` | `hesdirk2`: two implicit stages per step with the Heun hydro (second order); `be`: one backward-Euler solve per step (first order). | `hesdirk2` on a fresh run that supports it; `be` otherwise | all: `hesdirk2` | `rad_vet/rad_vet_init.cpp` (`ImplicitInit`) |
| `rad_vet/time2_tableau` | Tableau of the two stages: `sdirk2` (Alexander's L-stable SDIRK2, recommended), `hesdirk2`, `trbdf2`. | `hesdirk2` | all: `sdirk2` | `rad_vet/rad_vet_time2.cpp` (`Time2Init`) |
| `rad_vet/time2_auto_be` | After a step whose stages fell back to backward Euler, the next 1, 2, 4, ... N steps take backward Euler at once (0 = off). | 16 | all: `16` | `rad_vet/rad_vet_time2.cpp` (`Time2Init`) |
| `rad_vet/clean_init` | `steady`: before step 1 the radiation is relaxed to the scheme's steady solution in the given gas; `given`: start from the given field. | `given` | all: `steady` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/clean_init_dtfac`, `clean_init_tol`, `clean_init_n` | Pseudo time step (in units of the first dt), tolerance and pass limit of the steady start. | 1e6, 1e-8, 40 or 50 (problem dependent) | all: `1.0e6`, not named, not named | `rad_vet/rad_vet.cpp` (constructor) |
| `time/dt_retry` | A cycle whose implicit solve fails (not converged, non-finite, rejected stage, unpaid floor energy) is undone and taken again with a smaller dt. Bitwise: a retry restores the full cycle. | true | all: `true` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_fac`, `dt_retry_fac_t2` | dt factor of a retry; of a retry caused by a rejected stage. | 0.5; = `dt_retry_fac` | all: `0.5`, `0.9` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_grow`, `dt_retry_ngood` | After a cut, dt grows back by this factor per good cycle, starting after this many good cycles. | 1.25, 10 | all: `1.1`, `20` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_t2_be` | After this many stage-rejection cuts in one cycle, a rejected stage is redone with backward Euler (0 = at once, < 0 = never). | 0 | all: `0` | `driver/dt_retry.cpp` (`DtRetry`) |
| `time/dt_retry_force_cycle`, `dt_retry_force_fac` | Testing: force one retry at this cycle with this factor (the `dt_retry` test). | off | not named | `driver/dt_retry.cpp` (`DtRetry`) |

## Boundaries

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/implicit_bc_x1min`, `implicit_bc_x1max` | Radiation x1 boundaries: `marshak`, `outflow`, `reflect`, `periodic`, and on x1 an imposed `flux`. The ray sweep needs `x1max = marshak` or `reflect`. | `auto` (from the mesh boundary) | all: `flux`, `marshak` | `rad_vet/rad_vet_init.cpp` (`ImplicitInit`) |
| `rad_vet/implicit_flux_x1min` | Flux through the inner face (= L/(4 pi r_in^2) for a star). | 0 | per deck | `rad_vet/rad_vet_init.cpp` (`ImplicitInit`) |
| `rad_vet/implicit_flux_x1max` | Incoming bath at the outer Marshak face (0 = vacuum). | 0 | all: `0.0` | `rad_vet/rad_vet_init.cpp` (`ImplicitInit`) |
| `rad_vet/clean_bc` | Open x1 end: `halfrange` (the end cell's ray half-range flux) or `marshak`. | `halfrange` | all: `halfrange` | `rad_vet/rad_vet.cpp` (constructor) |
| `mesh/ix1_bc` ... `ox3_bc` | Hydro boundaries (the templates: `user` in x1, `periodic` in theta and phi). | required | see the decks | `mesh/mesh.cpp` |

Transverse Marshak faces are not supported by the ray sweep (FATAL).

## Opacity, units and coupling

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/opacity` | `const`, `powerlaw`, `user` or `table` (Rosseland and Planck tables named in `<problem>`). | `const` | all: `table` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/kappa_p`, `kappa_e`, `kappa_f`, `kappa_s` | Constant opacities per unit mass: Planck (emission), energy mean (absorption), flux mean, scattering. | 0, = kappa_p, = kappa_p, 0 | not named | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/c_light`, `arad` | Speed of light and radiation constant in code units. The start-up log compares them with the values `<units>` implies. | required | per deck | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/chat_over_c` | Reduced speed of light; the clean scheme needs 1 (FATAL otherwise). | 1.0 | all: `1.0` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/temp_unit_kelvin`, `rho_unit_cgs`, `kappa_unit` | Unit conversions for the table lookup: K per code temperature, g/cm^3 per code density, code opacity per cm^2/g. | 1.0 | per deck | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/coupling`, `gas_feedback` | Gas-radiation exchange on; exchange acts back on the gas. | true | all: `true` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/force_reference`, `force_reference_work` | `wb_arad`: well-balanced radiative envelope; the problem generator supplies the reference radiative acceleration, the gas receives only the residual force; its work is paid as `split` (required by the clean scheme). | `none`; `auto` | all: `wb_arad`, `split` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/e_floor` | Floor of E (code units). | FLT_MIN | all: `1.0e-4` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/floorcap_rho`, `floorcap_v` | Floor-gas cap: in cells with rho below `floorcap_rho`, |v| is capped at `floorcap_v` and T at max(T_rad, floor); the removed energy goes to E (momentum is removed and reported). | off | all: `1e-14`, `3e7` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/floorcap_lz` | The floor-gas velocity cap scales only the momentum without L_z: spherical polar r / theta, Cartesian the part off phi_hat = (-y, x)/R. The phi part is kept, also where it alone exceeds `floorcap_v`, so L_z is exact. FATAL on the cubed sphere. | false (read only when named) | not named | `rad_vet/rad_vet.cpp` (constructor) |
| `hydro/dfloor`, `pfloor`, `tfloor` | Hydro floors. `hydro/floor_pay = rad_column` pays the internal energy they create from the radiation (needs `hydro/eos = general`). | FLT_MIN | per deck | `eos/eos.cpp`, `hydro/dual_energy.cpp` |
| `hydro/floor_pay_skip` | Pay also the energy that the density floor and the velocity ceiling change. Needs `floor_pay`. | false (read only when named) | all: `true` | `hydro/dual_energy.cpp` |
| `hydro/floor_from_ekin` | General-EOS hydro: the primitive recovery keeps E; the velocity ceiling's kinetic energy becomes heat, and a kinetic-energy-dominated cell below e(rho, tfloor) pays the deficit from its own r and theta kinetic energy where that covers it. The phi momentum (L_z) is never scaled (where the phi part alone exceeds `vceil`, the r / theta part goes to 0 and |v| = |v_phi| may stay above `vceil`; Cartesian phi_hat = (-y, x)/R); cells it cannot pay take the floors as before (paid by `floor_pay`). Removes r / theta momentum there (not booked). FATAL with another EOS or `dual_energy`. | false (read only when named) | not named (Caltech t9l62 arm) | `eos/eos.cpp` (constructor) |
| `hydro/dfloor_flux_limit` | Scale back the outflow of a nearly empty cell so its density stays above the floor (both neighbours see the same flux). | false (read only when named) | all: `true` | `hydro/hydro.cpp` |

## Layout and reproducibility

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `mesh/layout_invariant` | Global and shell sums in exact, order-independent form: results independent of the MPI rank count. | true | all: `true` | `mesh/mesh.cpp` (constructor) |
| `mesh/layout_invariant_mode` | `block`: sums per MeshBlock (independent of the rank count for fixed MeshBlocks); `exact`: every cell term exact (about +8 % cost); `tile`: fixed N x N lateral tiles as the summation and multigrid units, so results do not depend on the MeshBlock layout either (with ray walls in theta: one MeshBlock along x2). | `block` | `tile` (recommended) | `mesh/mesh.cpp` (constructor) |
| `rad_vet/vet_gd_band_ring` | Lend the band intensities per shell from a few slabs: same results, less GPU memory (RUN.md section 6). | on where it applies, off on the raw sp band (read only when named) | A100 sp wedges: `true` | `rad_vet/rad_vet.cpp` (constructor) |
| `mesh/layout_tile` | Tile size N for `tile` mode; required with it. Needs one MeshBlock along x1 and N dividing `<meshblock>` nx2 and nx3. | none | `128` | `mesh/mesh.cpp` (constructor) |
| `mesh/nghost` | Ghost cells; the face line solves across MeshBlocks need >= 3. | 2 | all: `3` | `mesh/mesh.cpp` (constructor) |
| `meshblock/nx1`, `nx2`, `nx3` | MeshBlock size. One MeshBlock along x1 is needed for the ray cut and tile mode. | = mesh | per deck | `mesh/mesh.cpp` (constructor) |

Changing the mode changes the round-off of the sums: switch at a fresh start, or accept a
round-off change at the restart.

## Output and verbosity

| key | meaning | default | production | read in |
|---|---|---|---|---|
| `rad_vet/clean_verbose` | 0: the `CLEAN cycle` lines, `CLEAN WARNING` lines and event lines (`CLEAN WORST`, `GNSW`, `BICGSTAB`, `KPOS` without a face slope); 1: also the per-solve diagnostics (`CLEAN MCAP`, `OLDVEC`, `FWB`, `FCONE`, `LZ`/`LZB`, `WBFCT`, `KPOS`, `KREAL`, `XPAT`, the `XT` write-back closure, the ray cut). `CLEAN WARNING` lines always print. | 0 (read only when named) | all: `1` | `rad_vet/rad_vet.cpp` (constructor) |
| `rad_vet/floorcap_report` | Print the cumulative floor-gas cap counters every N cycles. | 0 = never | all: `25` | `rad_vet/rad_vet.cpp` (constructor) |
| `time/ndiag` | Print the cycle line every N cycles. | 1 | all: `50` | `driver/driver.cpp` |
| `hydro/floor_audit` | Print the floor events and their payment every N cycles (`floor_audit cycle` lines). Needs `eos = general`. | 0 | all: `1` | `hydro/dual_energy.cpp` |
| `hydro/dt_report` | Print every N cycles which cell and term set dt (`DTREP`, `DTREG`; RUN.md section 5). No effect on the state. | 0 (read only when named) | not named | `hydro/hydro.cpp` |
| `time/nan_check_cycles` | Check the state for NaN every N cycles and stop if found (0 = off). | 100 | all: `100` | `driver/driver.cpp` |

## Appendix: advanced and internal keys

These keys are named in the templates to pin the scheme (keep their values) or are solver
internals. Do not change them without a dedicated test (energy and momentum conservation, a
bitwise restart check, and a comparison with the production key set).

Scheme options pinned by the templates:

- `clean_stag_upwind_central` (true), `clean_stag_central` (true): central advective flux within the TVD bound.
- `clean_face_hctau` (true), `clean_face_jse` (true), `clean_face_opac` (`harmonic`): face values of the rays' flux, half ranges and opacity.
- `clean_stag_thinf` (false), `clean_consistency` (true): forms of the thin-limit correction and flux bound.
- `clean_stag_fwb` (true), `clean_stag_fcone` (true), `clean_stag_fwb_ring` (1: 3-D spherical polar with periodic phi only): well-balanced face rows, cone limit |F| <= cE, ring projection.
- `clean_stag_fcone_gas` (false): the cone cut's momentum is not given to the gas.
- `clean_stag_fcone_lz` (read only when named): fcone takes a cell's excess from the face families without an L_z arm first (spherical polar x1, x2; cubed sphere x1; Cartesian x3); that cut is counted and printed signed per family (`CLEAN FCONE_LZ`). The faces with an arm are scaled only where they alone exceed cE, and the gas takes their cut at the face, so L_z is exact.
- `clean_lz_spin` (read only when named; needs `clean_stag_fcone_lz`; a no-op on the cubed sphere, where `clean_cs_lz_tel` covers it): the L_z X of the arm faces' fcone cut, the `clean_stag_fcone_emax` clamp and the `clean_wb_mcap` cap goes to the gas as one rigid rotation dOmega = X/I, I = sum rho V R^2 (no per-cell kick; E pays the kinetic energy; `CLEAN LZSPIN`).
- `clean_stag_fcone_emax` (read only when named; needs `clean_stag_fcone`): every inner face's F within c max(E_L, E_R), an opaque face's too. L_z-exact with `clean_stag_fcone_lz` + `clean_lz_spin`.
- `clean_wb_mcap` (read only when named): the write-back gives the gas the face momentum terms at the Newton passes' cap h chi_f E; the withheld momentum is printed per family (`CLEAN MCAP write-back`). L_z-exact only with `clean_lz_spin`.
- `clean_stag_fwb_cart` (read only when named; uniform Cartesian meshes only, no-op elsewhere): the well-balancing residual minus its line mean along every face direction, so `clean_stag_fwb` gives no net linear momentum.
- `clean_gas_secant` (`max_stage2`): in the stage-2 solves the gas-exchange slope is the chord through the last pass's T where it beats the tangent by more than a derived round-off bound (16 eps of the terms of e(T) - e(T_p)), else the tangent.
- `clean_pred_xch` (true): the stage-2 start blends this step's stage-1 slope and the last step's stage-2 slope by the stiffness of each cell's exchange.
- `clean_wb_exch_local` (true): the write-back's cell exchange is solved exactly from the transported E (E > 0 and gas >= e(rho, tfloor) by construction); needs an EOS `tfloor` > 0.
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
`clean_clos1_budget`, `clean_vel_settled`, `clean_bud_maxnorm` (speed keys, read only when named, default off: the closure pass solved
in the budget norm; the velocity update only in a pass that can end the solve; the budget test in the per-cell max norm), `clean_stag_mg_color`, `clean_stag_fs_line`, `clean_stag_pcr`, `clean_halo_direct` (only true),
`implicit_predictor`, `implicit_predictor_tcap`, `implicit_resid_fatal`, `time/dt_retry_views`,
`time/dt_retry_host`, `time/restart_refill_ghosts` (false: restarts are bitwise).

Removed (only false / 0 accepted): `clean_stag_kpat`, `time2_k1_redist` (old restart files
may store 0.1: name it 0). Not recommended: `time2_k1_stiff` (keep 0).

Removed study and debug keys stop the run with a message when named with another value than
their old default (`CleanPreflight` in `rad_vet/rad_vet_clean.cpp` lists them); keys of the
removed solvers that old decks and restart files carry are accepted without effect.

Named in the templates but without effect with them (fixed rays, spherical polar): the
rotation keys other than `vet_gd_rotate_dt`, the cubed-sphere keys (`clean_cs_lz_tel`,
`cs_cross_mono`, `cs_vertex_mono`, `mesh/cs_*`), `vet_source`, `vet_scatter_j` and the other
relaxed-source keys (`clean_ray_source_relaxed` was removed: only false), `vet_col_lat_taucut`,
`implicit_maxit`, `implicit_lin_maxit`, `time2_enth_be`, `time/dt_retry_max`. Keys that are
never read are listed by the `never read` warning after the first cycle.
