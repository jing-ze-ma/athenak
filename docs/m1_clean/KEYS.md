# Input keys of the clean M1 radiation scheme (production reference)

This page lists every key of the `<rad_m1>`, `<time>`, `<mesh>` and `<meshblock>` blocks, and
the `<hydro>` keys that matter for the radiation coupling and the floors, that the three
production decks use:

| deck | star | mesh (r x theta x phi) | MeshBlock |
|---|---|---|---|
| `inputs/rad_m1_clean/star_bsg.athinput` | blue supergiant (Ma et al. 2026 setup) | 256 x 256 x 256 | 256 x 128 x 128 |
| `inputs/rad_m1_clean/star_m35mms.athinput` | T35L5.2 (Schultz, Bildsten and Jiang 2022) | 336 x 256 x 256 | 336 x 128 x 128 |
| `inputs/rad_m1_clean/star_t9l62.athinput` | T9L6.2 (Jiang et al. 2018; Schultz et al. 2020) | 512 x 256 x 512 | 512 x 128 x 128 |

All three are 3-D spherical-polar wedges with periodic theta and phi boundaries, one MeshBlock
along x1, fixed (non-rotating) rays and the `he_star_m1` problem generator. The code is the
release commit 512cf768a.

The list was made by a script that reads the three decks (skipping `<comment>` blocks and
commented-out lines), takes the union of the keys, and searches `src/` for each key's
`pin->Get*` / `GetOrAdd*` / `DoesParameterExist` call. 249 keys are in scope. For each
key the table gives:

- **meaning**: what the key does, in one or two sentences;
- **type**;
- **code default**: the value the code uses when the key is absent. Many defaults depend on
  whether the run is a **fresh start** or a **restart**, and on the mesh; this is stated;
- **production**: the value in the decks (`all:` when the three agree);
- **read at**: the file and line of the read, and the function.

Markers: **[never read]** = no code reads the key with these decks (see the list at the end);
**[read, no effect here]** = read, but it has no effect with these decks or under
`scheme = clean` (the reason is given).

Words used below:

- *E*, *F*: the radiation energy density and flux (code units = cgs here: erg/cm^3, erg/cm^2/s).
- *staggered*: E lives at cell centres, the flux F on the cell faces (one F per face).
- *rays*: the short-characteristics sweep over the whole mesh (`vet_gd`) that gives the
  Eddington tensor (the closure) and the face half ranges (the outgoing and incoming parts of
  the free-streaming flux).
- *budget tolerance*: the solver measures convergence in each cell relative to the energy that
  cell exchanges in the step (`clean_conv = budget`).
- *the old solver*: the implicit M1 solver that came before the clean scheme. Some keys
  (`implicit_*`) belong to it; a few of them still set options that the clean scheme uses.
- *named*: the key appears in the input (deck, restart file or command line).

## The restart rule: name every production key

A restart is started with `athena -r <file>.rst [-i <deck>] [block/key=value ...]`. The
parameters are put together in this order (`src/main.cpp`, lines 495-557):

1. the input block **stored in the restart file** (the parameters of the run that wrote it);
2. the `-i` deck, if one is given: every key it names replaces the stored value;
3. the command line: every `block/key=value` replaces the value so far. On a restart, a key
   that is not in the stored block is added (the code prints a `### NOTE`). On a **fresh**
   start a command-line key that is not in the deck is FATAL (typo guard), so put every key
   of a fresh run in the deck.

**On a restart, any key that is NOT given in the deck or on the command line takes the value
stored in the restart file, NOT the code default.** If the restart file does not have the key
either, many clean-scheme keys then take their *restart* default, which is usually the old
behaviour (for example "off"), so that old runs continue bitwise. The "code default" column
says "on a restart without the key" where this applies.

Three kinds of defaults exist:

- `GetOrAdd` defaults: the default is written into the run's parameters and therefore into
  every restart file; later restarts keep it even if the code default changes.
- "read only when named": nothing is written; a restart file that does not have the key
  gets the current code default.
- "fresh / restart" defaults: a fresh run takes and records the new value; a restart whose
  file lacks the key takes the old value.

Example of what goes wrong: a production restart file had stored `time2_k1_stiff = 1.0` from an
earlier deck. The code default is 0 (off), but the restart silently ran with 1.0, because the
new deck did not name the key.

Rules that follow:

- When you restart from an old restart file, restart with `-i <production deck>` (or name
  every production key on the command line). The grid keys (`<mesh>` x1min..x3max, nx1..nx3,
  nghost, the stretch keys and `<meshblock>` nx1..nx3) must stay as the file was written:
  a change is FATAL (`CheckRestartGridKeys`, `src/main.cpp`).
- To see the parameters a restart will run with, `athena -r <file>.rst -i <deck> -n` prints
  the merged parameter list and quits (keys the modules add later are not shown).
- After the first cycle the code prints `### WARNING: N input parameter(s) were never read by
  the code`: a key there is a typo or a key this build does not have.

### Keys that must always be named

These are named in all three decks. Earlier production restart files carry other values, so
leaving any of them out of a restart deck changes the scheme silently.

| key | value | why |
|---|---|---|
| `<rad_m1>/time2_k1_stiff` | `0` | read only when named; old files stored 1.0 |
| `<rad_m1>/clean_stag_fwb_stage` | `false` | earlier production files stored true |
| `<rad_m1>/clean_stag_fwb_t2` | `false` | earlier production files stored true |
| `<rad_m1>/clean_w_fse_tau` | `true` | earlier files may store false |
| `<rad_m1>/time2_gtime_mom` | `true` | earlier files may store false |
| `<rad_m1>/clean_pred_stage` | `false` | production key set |
| `<rad_m1>/clean_stag_central` | `true` | production key set (default false) |
| `<rad_m1>/clean_stag_kpat` | `false` | production key set |
| `<rad_m1>/clean_stag_eicarry` | `true` | production key set (default false) |
| `<mesh>/layout_invariant_mode` | `tile` | recommended for production (see below) |
| `<mesh>/layout_tile` | `128` | the lateral MeshBlock size of these decks |

`<mesh>/layout_invariant_mode = tile` with `<mesh>/layout_tile = N` makes the fixed-order sums
use fixed N x N lateral tiles of the mesh instead of MeshBlocks. The history sums and the
multigrid preconditioner are then the same for every MeshBlock layout made of whole tiles,
at the cost of block mode. Conditions (FATAL otherwise): one MeshBlock along x1, N divides
`<meshblock>` nx2 and nx3, and the default multigrid preconditioner (`src/mesh/mesh.cpp:186`,
`src/rad_m1/rad_m1.cpp` near line 790). `layout_tile` has no default: it must be given with
`tile`. The templates in `inputs/rad_m1_clean/` use tile mode; the 2026-11-23 production
runs ran block mode, which differs from them only in these two lines. Changing the mode
changes the round-off of the sums, so switch at a fresh start or accept a round-off change at
the restart.

## Time integration, time step and start

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_init` | Start: steady = before step 1 the radiation is relaxed to a steady solution of the given gas (needed when the initial state is out of radiative equilibrium); given = start from the given field. | string | given (read only when named) | all: `steady` | `src/rad_m1/rad_m1.cpp:965 RadiationM1 constructor` |
| `rad_m1/clean_init_dtfac` | Pseudo time step of the steady start, in units of the first dt. | real | 1.0e6 (or problem/he_ic_equil_dtfac) | all: `1.0e6` | `src/rad_m1/rad_m1.cpp:985 RadiationM1 constructor` |
| `rad_m1/clean_pred_stage` | sdirk2: the stage-2 solve starts from E_old + h K1 of this step's stage 1 instead of the previous step's stage 2. Production: false (on costs t9l62 about 11 percent: stage-2 passes 2.0 -> 3.1). One of the must-name keys (above). | bool | false | all: `false` | `src/rad_m1/rad_m1.cpp:222 RadiationM1 constructor` |
| `rad_m1/time2_auto_be` | After a step whose stages fell back to the backward-Euler redo, the next steps skip the stage attempt and take the redo at once; the stages are retried after 1, 2, 4, ... up to N skipped steps (0 = off). | int | 16 (read only when named) | all: `16` | `src/rad_m1/rad_m1_time2.cpp:648 RadiationM1::Time2Init` |
| `rad_m1/time2_cfl_post` | The next hydro CFL dt is recomputed after the radiation update, on the gas the radiation solve wrote (avoids a step that runs past its CFL after a radiation kick). | bool | true (read only when named) | all: `true` | `src/rad_m1/rad_m1_time2.cpp:644 RadiationM1::Time2Init` |
| `rad_m1/time2_gtime_mom` | sdirk2 only: the hydro stage-2 flux state takes the momentum at the stage time (factor f = 1/g); gives second-order-in-time wall momentum. Must be named true for restarts from files that stored false. | bool | true (read only with time2_tableau = sdirk2) | all: `true` | `src/rad_m1/rad_m1_time2.cpp:147 RadiationM1::Time2Init` |
| `rad_m1/time2_start` | How a step without a stored slope (first step, restart without the slope) starts: be = one Lie-split backward-Euler step (first order for that step); slope = a backward-Euler probe gives K1 and the step runs the stages. | string | be (read only when named) | all: `be` | `src/rad_m1/rad_m1_time2.cpp:389 RadiationM1::Time2Init` |
| `rad_m1/time2_tableau` | Coefficients of the two-stage scheme: hesdirk2, trbdf2 or sdirk2. Production uses sdirk2 (no explicit first slope K1, Heun hydro; both stages implicit with diagonal g = 1 - 1/sqrt(2)). | string | hesdirk2 (read only when named) | all: `sdirk2` | `src/rad_m1/rad_m1_time2.cpp:123 RadiationM1::Time2Init` |
| `rad_m1/time_scheme` | Radiation time scheme. hesdirk2 = two implicit stages per step with the Heun hydro (second order); be = one backward-Euler solve per step (first order). | string | hesdirk2 on a fresh run that supports it (implicit transport, <time>/integrator = rk2, no polar boundary); be otherwise and on a restart without the key | all: `hesdirk2` | `src/rad_m1/rad_m1_init.cpp:134 RadiationM1::ImplicitInit` |
| `hydro/rad_signal_speed` | The hydro CFL uses the radiation-modified signal speed. | bool | on with implicit M1 transport (read only when named) | all: `true` | `src/hydro/hydro.cpp:456 Hydro constructor` |
| `time/cfl_number` | Courant number of the hydro time step. The radiation is implicit and does not limit dt. | real | 0.3 on a fresh clean M1 run if not named (recipe::CleanDefaults); otherwise required | all: `0.3` | `src/mesh/build_tree.cpp:415 Mesh::BuildTreeFromScratch (default added in src/utils/recipe.cpp:68)` |
| `time/dt_retry` | A cycle whose implicit solve fails (not converged, non-finite, rejected stage, unpaid floor energy) is undone and taken again with a smaller dt, instead of being accepted or stopping the run. | bool | true (also on a restart without the key, unless <job>/recipe = none) | all: `true` | `src/driver/dt_retry.cpp:443 DtRetry constructor` |
| `time/dt_retry_fac` | dt factor of a retried cycle (dt -> fac dt). | real | 0.5 | all: `0.5` | `src/driver/dt_retry.cpp:464 DtRetry constructor` |
| `time/dt_retry_fac_t2` | dt factor of a retry caused by a rejected hesdirk2/sdirk2 stage (reason time2_reject). | real | = dt_retry_fac (read only when named) | all: `0.9` | `src/driver/dt_retry.cpp:475 DtRetry constructor` |
| `time/dt_retry_grow` | After a cut, dt is allowed to grow back by this factor per good cycle (up to the CFL dt). | real | 1.25 | all: `1.1` | `src/driver/dt_retry.cpp:465 DtRetry constructor` |
| `time/dt_retry_max` **[never read]** | NOT READ. The maximum number of cuts per cycle is fixed at 8 in the code (src/driver/dt_retry.hpp:174, maxcut_); the deck value 8 matches it. | int | (no key) | all: `8` | (no read in `src/`) |
| `time/dt_retry_ngood` | Number of good cycles after a cut before dt starts to grow again. | int | 10 | all: `20` | `src/driver/dt_retry.cpp:466 DtRetry constructor` |
| `time/dt_retry_t2_be` | After this many stage-rejection cuts in one cycle, a further rejected hesdirk2 stage is redone with backward Euler instead of cutting dt again; 0 = redo with backward Euler at the cycle dt without cutting; < 0 = never. | int | 0 (read only when named) | all: `0` | `src/driver/dt_retry.cpp:483 DtRetry constructor` |
| `time/evolution` | dynamic = evolve hydro and radiation in time. | string | none (required) | all: `dynamic` | `src/driver/driver.cpp:109 Driver constructor` |
| `time/integrator` | Hydro time integrator. Must be rk2: the radiation time scheme hesdirk2/sdirk2 uses the Heun (rk2) hydro as its explicit part. | string | rk2 | all: `rk2` | `src/driver/driver.cpp:125 Driver constructor` |
| `time/ndiag` | Print the cycle line (cycle, time, dt) every ndiag cycles. | int | 1 | all: `50` | `src/driver/driver.cpp:128 Driver constructor` |
| `time/nlim` | Maximum number of cycles (-1 = no limit). | int | -1 | all: `-1` | `src/driver/driver.cpp:127 Driver constructor` |
| `time/tlim` | Stop time of the run [s]. | real | none (required) | bsg: `4.9248e6`; m35mms: `8.64e5`; t9l62: `2.682e6` | `src/driver/driver.cpp:126 Driver constructor` |

## Closure and rays

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_closure` | How the clean solve takes the closure in time: extrap1 (default) builds the rays' source from converged states; alternatives predict, iterate, extrap, extrap2, eddington, sphere. | string | extrap1 | all: `extrap1` | `src/rad_m1/rad_m1.cpp:861 RadiationM1 constructor` |
| `rad_m1/clean_fill_cs_remap` **[read, no effect here]** | Cubed sphere only: the fill takes the half ranges of a seam neighbour in the local face families. Inert on spherical polar. | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:212 RadiationM1 constructor` |
| `rad_m1/clean_fill_passes` | Number of fill layers for the closure in ray-shadowed cells. | int | 8 | all: `8` | `src/rad_m1/rad_m1.cpp:209 RadiationM1 constructor` |
| `rad_m1/clean_predict` | The closure is evaluated on the predicted state. | bool | true | all: `1` | `src/rad_m1/rad_m1.cpp:213 RadiationM1 constructor` |
| `rad_m1/clean_ray_cut_margin` | Optical depth margin of the ray cut below the first non-diffusive cell. | real | 10.0 | all: `10` | `src/rad_m1/rad_m1.cpp:2797 RadiationM1 constructor` |
| `rad_m1/clean_ray_cut_r` | Ray depth cut (on by default for sc rays on 3-D spherical polar with one MeshBlock along x1): below the cut the closure is the diffusion limit. A column's cut sits clean_ray_cut_margin optical depths below its first cell with R = \|grad S\|/(chi S) >= clean_ray_cut_r. | real | 0.01 (read where the ray cut applies) | all: `0.01` | `src/rad_m1/rad_m1.cpp:2796 RadiationM1 constructor` |
| `rad_m1/clean_ray_cut_taumin` | The ray cut is never placed where the optical depth to the top is below this. | real | 0.0 | all: `0` | `src/rad_m1/rad_m1.cpp:2798 RadiationM1 constructor` |
| `rad_m1/clean_ray_cut_wait` | The cut moves down at once but moves up only after this many builds that would raise it (stored in the restart file). | int | 10 | all: `10` | `src/rad_m1/rad_m1.cpp:2799 RadiationM1 constructor` |
| `rad_m1/clean_ray_min` | The rays describe a cell only where their J is at least this fraction of the cell's E; elsewhere the cell is "ray-shadowed" and uses its own half ranges. | real | 1.0e-3 | all: `0.001` | `src/rad_m1/rad_m1.cpp:166 RadiationM1 constructor` |
| `rad_m1/clean_ray_source_relaxed` | Use the vet_source / vet_scatter_j relaxed form as the rays' source of the clean scheme (off = the state's own source). | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:127 RadiationM1 constructor` |
| `rad_m1/clean_rays` | Rays of the clean scheme: sc = the vet_gd short-characteristics sweep (the only form; others are FATAL). | string | sc | all: `sc` | `src/rad_m1/rad_m1.cpp:377 RadiationM1 constructor` |
| `rad_m1/clean_ring_jmin` **[read, no effect here]** | Ring J ratios are measured only where both sets hold J >= clean_ray_min E. Inert with clean_rot_ring = off. | bool | true | all: `1` | `src/rad_m1/rad_m1.cpp:693 RadiationM1 constructor` |
| `rad_m1/clean_rot_fuse` **[read, no effect here]** | Sweep the live ring's two sets in one fused sweep (bitwise). Inert with clean_rot_ring = off. | bool | true (fresh run; a restart keeps the value stored in the file) | all: `true` | `src/rad_m1/rad_m1.cpp:183 RadiationM1 constructor` |
| `rad_m1/clean_rot_ring` | Rotating ray ring of the clean scheme: off, live, jump. Production: off (fixed rays). | string | read only when named (live where the rays rotate with a cross-fade) | all: `off` | `src/rad_m1/rad_m1.cpp:185 RadiationM1 constructor` |
| `rad_m1/clean_rot_ring_fill_rst` **[read, no effect here]** | At a restart from a file without the clean state, the ring's older sets are filled from a sweep of the restart state. Inert with clean_rot_ring = off. | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:698 RadiationM1 constructor` |
| `rad_m1/clean_source_convex` | Keep the scattering ray source between J and B (absorption share kappa_P/chi at most 1), so it cannot go negative in cool envelopes where kappa_P > chi. | bool | true | all: `1` | `src/rad_m1/rad_m1.cpp:709 RadiationM1 constructor` |
| `rad_m1/closure` | Closure of the moment equations. vet_col + vet_gd = true: a variable Eddington tensor from a short-characteristics ray sweep over the whole mesh (the vet_gd rays). | string | vet_col with vet_gd = true on a fresh clean run on spherical polar or the cubed sphere (recipe::CleanDefaults); m1 otherwise | all: `vet_col` | `src/rad_m1/rad_m1.cpp:1841 RadiationM1 constructor` |
| `rad_m1/time2_vet_col_source` | Alias of time2_vet_source: the rays' source during the stages (blend = the state's source where the exchange is resolved, relaxed where it is stiff). Effect under scheme = clean not traced. | string | blend on a fresh vet_col/vet_gd run | all: `blend` | `src/rad_m1/rad_m1_time2.cpp:2653 RadiationM1::Time2VetSourceInit (alias)` |
| `rad_m1/vet_col_every` | Alias of vet_every: rebuild the vet_col closure every N steps. Must be 1 with the time2 vet_col modes; with vet_gd_replace = true no vet_col sweep is run. | int | 1 | all: `1` | `src/rad_m1/rad_m1.cpp:1889 RadiationM1 constructor (via M1KeyName, alias of vet_every)` |
| `rad_m1/vet_col_lat_cut_thin` | Thinness criterion of the ray depth cut (also used by clean_ray_cut): a column is cut no deeper than where it is thinner than this. | real | 1.0 on a fresh run; 0 on a restart without the key | all: `1.0` | `src/rad_m1/rad_m1.cpp:2596 RadiationM1 constructor` |
| `rad_m1/vet_col_lat_offdiag` | Off-diagonal (shear) terms of the lateral closure: operator = put into the operator; none = off. | string | read only when named | all: `operator` | `src/rad_m1/rad_m1.cpp:2605 RadiationM1 constructor` |
| `rad_m1/vet_col_lat_taucut` **[read, no effect here]** | Optical-depth cut below which the lateral sweep is skipped. INERT under scheme = clean: the code sets the cut to +infinity (rays from the first shell, src/rad_m1/rad_m1.cpp:2596 area). | real | read only when named | all: `30.0` | `src/rad_m1/rad_m1.cpp:2589 RadiationM1 constructor` |
| `rad_m1/vet_gd` | Use the global ray sweep (vet_gd) for the Eddington tensor. | bool | true with closure = vet_col on a fresh clean sp/cs run (recipe::CleanDefaults); else read only when named | all: `true` | `src/rad_m1/rad_m1.cpp:367 RadiationM1 constructor` |
| `rad_m1/vet_gd_allow_clamp` | Allow the lateral reach of the ray sweep to be clamped to the band the MeshBlocks hold, instead of stopping the run (FATAL) when the band is too small. Binary compatibility with older production runs, not physics. | bool | false (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:2295 RadiationM1 constructor` |
| `rad_m1/vet_gd_gl_nmu` | Number of polar (Gauss-Legendre) angles of the product ray set. Takes precedence over vet_gd_ls. Production: 7 (x 18 azimuths = 126 directions). | int | no set named: geodesic set vet_gd_geo = 3 on a fresh 3-D run | all: `7` | `src/rad_m1/rad_m1.cpp:2106 RadiationM1 constructor` |
| `rad_m1/vet_gd_gl_nphi` | Number of azimuths of the product ray set; 2 pi/nphi must divide the phi extent of the wedge (FATAL otherwise). | int | see vet_gd_gl_nmu | all: `18` | `src/rad_m1/rad_m1.cpp:2109 RadiationM1 constructor` |
| `rad_m1/vet_gd_gl_stagger` | Stagger the azimuths of alternate polar rings of the product set. | bool | true (with the product set) | all: `true` | `src/rad_m1/rad_m1.cpp:2112 RadiationM1 constructor` |
| `rad_m1/vet_gd_ls` **[read, no effect here]** | Level-symmetric ray set LQ_N (N = 6, 8, 10, 12). INERT here: vet_gd_gl_nmu > 0 is named, and the product set takes precedence (src/rad_m1/rad_m1_vetgd.cpp:790 VgdBuildSet). | int | 0 (read only when named) | all: `10` | `src/rad_m1/rad_m1.cpp:2115 RadiationM1 constructor` |
| `rad_m1/vet_gd_replace` | The vet_gd sweep replaces the per-column vet_col sweep (no vet_col formal solution is built). | bool | false (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:2562 RadiationM1 constructor` |
| `rad_m1/vet_gd_rotate_align` **[never read]** | Align the rotation clock with the time step. NOT READ here: read only when the set rotates on the time clock (vet_gd_rotate_dt > 0 or auto). | bool | true | all: `true` | `src/rad_m1/rad_m1_vetgd.cpp:1446 RadiationM1::VetRotInit` |
| `rad_m1/vet_gd_rotate_dt` | Rotate the ray set on a physical-time clock every dt_rot (0 = fixed set). Production: 0 (rotation off). | real or auto | auto (20 x the first dt) on a fresh vet_gd run that names no rotation key; no rotation on a restart without the key | all: `0` | `src/rad_m1/rad_m1.cpp:2093 RadiationM1 constructor` |
| `rad_m1/vet_gd_rotate_xfade` **[read, no effect here]** | Cross-fade the closure from the current to the next ray set. Inert with a fixed set (the code ANDs it with "set rotates"). | bool | true with a rotating set on a fresh clean run; read only when named | all: `true` | `src/rad_m1/rad_m1.cpp:2168 RadiationM1 constructor` |
| `rad_m1/vet_gd_tangential` | Use the full Eddington tensor including the tangential components (false = D_rr only). | bool | true (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:2570 RadiationM1 constructor` |
| `rad_m1/vet_gd_wall_interp` | Interpolation of the ray sweep at the theta/phi walls of the wedge (see VetGdSweep). | bool | false (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:2430 RadiationM1 constructor` |
| `rad_m1/vet_interp` | Upwind interpolation of the short characteristics: auto = corner segments where the opacity contrast of the stencil is large, log (bilinear in ln chi, ln S) elsewhere. | string | auto on a fresh run; log on a restart without the key | all: `auto` | `src/rad_m1/rad_m1.cpp:2257 RadiationM1 constructor` |
| `rad_m1/vet_rotate_advance` **[read, no effect here]** | The averaged window of rotated sets advances. Acts only with vet_rotate_avg >= 2: inert here. | bool | true on a fresh vet_gd run | all: `true` | `src/rad_m1/rad_m1_vetgd.cpp:1394 RadiationM1::VetRotInit` |
| `rad_m1/vet_rotate_avg` **[read, no effect here]** | Average the closure over the last K rotated sets (0 = off). Inert with a fixed set. | int | 4 with a rotating set on a fresh run; read only when named otherwise | all: `0` | `src/rad_m1/rad_m1.cpp:3312 RadiationM1::RotAvgParse` |
| `rad_m1/vet_rotate_avg_ema` **[read, no effect here]** | Exponential moving average over the rotated sets (0 = off). Inert with a fixed set. | real | 128 with a rotating set on a fresh run; read only when named otherwise | all: `0` | `src/rad_m1/rad_m1.cpp:3351 RadiationM1::RotAvgParse` |
| `rad_m1/vet_rotate_avg_ema_change_form` **[never read]** | Form of the EMA change test. NOT READ here: read only when vet_rotate_avg_ema >= 1. | string | linear | all: `linear` | `src/rad_m1/rad_m1.cpp:3372 RadiationM1::RotAvgParse` |
| `rad_m1/vet_scatter` | Split absorption and scattering opacity in the ray source (kappa_a vs kappa_e). | bool | true on a fresh vet-ray run with table opacity and a tabulated EOS; else read only when named | all: `true` | `src/rad_m1/rad_m1.cpp:1384 RadiationM1 constructor` |
| `rad_m1/vet_scatter_form` | Form of the scattering split: ma (Jiang 2021 eq. 6) or absorption. | string | ma | all: `ma` | `src/rad_m1/rad_m1.cpp:1387 RadiationM1 constructor` |
| `rad_m1/vet_scatter_j` | J of the scattering source: relaxed (end-of-step pair of B). current is FATAL with form ma. | string | relaxed | all: `relaxed` | `src/rad_m1/rad_m1.cpp:1399 RadiationM1 constructor` |
| `rad_m1/vet_scatter_kappa_e` | Electron-scattering opacity: const (kappa_es = 0.34-type constant) or eos (from the EOS ionisation). | string | eos on a fresh run with a tabulated EOS; const otherwise | all: `const` | `src/rad_m1/rad_m1.cpp:1427 RadiationM1 constructor` |
| `rad_m1/vet_seg` | Integration of a characteristic segment through a cell: jump = split at the face with the opaque-to-thin jump term (production); mean, split, emit = older forms. | string | jump on a fresh run; mean on a restart without the key | all: `jump` | `src/rad_m1/rad_m1.cpp:1988 RadiationM1 constructor` |
| `rad_m1/vet_seg_plateau` **[read, no effect here]** | Plateau weights of the jump segment. As far as the code shows, the flag is used only by the vet_sc solver (src/rad_m1/rad_m1_vet.cpp), not by the vet_gd sweep used here. | bool | false (read only with vet_seg = jump) | all: `false` | `src/rad_m1/rad_m1.cpp:1999 RadiationM1 constructor` |
| `rad_m1/vet_source` | Source function of the rays: relaxed_esc = the local exchange with escape, relaxed where the exchange is stiff. | string | relaxed_esc on a fresh run; gas on a restart without the key | all: `relaxed_esc` | `src/rad_m1/rad_m1_vet.cpp:1962 RadiationM1::VetInit` |
| `rad_m1/vet_source_noesrc_h2be` **[read, no effect here]** | Backward-Euler solves read E without the explicit source. INERT here: it acts only together with vet_source_noesrc, which is off with vet_gd (src/rad_m1/rad_m1_time2.cpp:3558 VsrcNoes). | bool | true (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:2427 RadiationM1 constructor` |

## Staggered face scheme

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_consistency` | The face F solve carries the rays' kinetic momentum flux and their face flux as flux factors, so that the discrete radiation is well-balanced against the rays (needed for well-balancing with the vet_gd rays). | bool | true (false only on sp/cs without the vet_gd rays) | all: `1` | `src/rad_m1/rad_m1.cpp:381 RadiationM1 constructor` |
| `rad_m1/clean_face_area` | The sp/cs streaming share in area form also under gflux. | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:120 RadiationM1 constructor` |
| `rad_m1/clean_face_hclat` | Lateral faces: the bracket's flux is the mean of the two cells' H. | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:117 RadiationM1 constructor` |
| `rad_m1/clean_face_hctau` | The bracket's ray face flux interpolated in optical depth between the cell centres (the old solver's implicit_thin_hcface = tau). Production: true. | bool | false | all: `true` | `src/rad_m1/rad_m1.cpp:108 RadiationM1 constructor` |
| `rad_m1/clean_face_jse` | Face half ranges in the shape of the cells' half-range fluxes at the face (the old solver's implicit_thin_recon = jse). | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:111 RadiationM1 constructor` |
| `rad_m1/clean_face_jse_cell` | Apply the jse factor to the donor cell's half ranges instead of the reconstructed face value. | bool | false | all: `false` | `src/rad_m1/rad_m1.cpp:114 RadiationM1 constructor` |
| `rad_m1/clean_face_opac` | x1 face opacity: harmonic (distance-weighted harmonic mean of the two cells) or arith. | string | harmonic (fresh run; a restart keeps the value stored in the file); arith on a restart without the key | all: `harmonic` | `src/rad_m1/rad_m1.cpp:255 RadiationM1 constructor` |
| `rad_m1/clean_layout` | Variable layout of the clean scheme: staggered = cell-centred E with a face flux F on every face (the only layout in this build). | string | staggered on a fresh run; cell on a restart without the key (FATAL: removed) | all: `staggered` | `src/rad_m1/rad_m1.cpp:227 RadiationM1 constructor` |
| `rad_m1/clean_stag_central` | Central form: the rest of an upwind face takes the plain central mean, and the radiation work's drain on the E rows and the gas kick share the same Patankar weight. Needs clean_xt_exact = cons. Production: true. One of the must-name keys (above). | bool | false (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:413 RadiationM1 constructor (stag() lambda: read only when named)` |
| `rad_m1/clean_stag_eicarry` | A stage-start gas internal energy below 0 enters the gas row as it is (paid by the implicit exchange) instead of being taken as 0. Production: true. One of the must-name keys (above). | bool | false (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:576 RadiationM1 constructor (stag())` |
| `rad_m1/clean_stag_fcone` | Every cell's radiation flux vector is kept inside its realisable cone \|F\| <= c E by construction (the faces are scaled toward 0). | bool | true with clean_stag_fwb (fresh run; a restart keeps the value stored in the file) | all: `true` | `src/rad_m1/rad_m1.cpp:498 RadiationM1 constructor` |
| `rad_m1/clean_stag_lat_sums` | Lateral face stress from the cells' stress sums (the cubed-sphere form) instead of 4-cell central differences. Production: 0 (false). | bool | true (fresh run; a restart keeps the value stored in the file) with xdonor and no sbimp; false on a restart without the key | all: `0` | `src/rad_m1/rad_m1.cpp:845 RadiationM1 constructor` |
| `rad_m1/clean_stag_rayratio` | The bracket reads the rays only as ratios (face flux over the donor's ray J) times the state's E. | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:688 RadiationM1 constructor` |
| `rad_m1/clean_stag_thinf` | The staggered bracket's ray correction in the face half ranges in realisable form. Production: false (all three stars). | bool | true (fresh run; a restart keeps the value stored in the file) | all: `false` | `src/rad_m1/rad_m1.cpp:432 RadiationM1 constructor` |
| `rad_m1/clean_stag_upwind_central` | An upwind face carries the central mean of the cells' a E within the TVD bound (the old solver's cell flux). Implied by clean_stag_central = true. | bool | false (read only when named) | all: `true` | `src/rad_m1/rad_m1.cpp:409 RadiationM1 constructor (stag())` |
| `rad_m1/clean_stag_xdonor` | The off-diagonal stress flux as an M-matrix (donor) term: lag (production), p1, off. | string | lag (fresh run; a restart keeps the value stored in the file); off for an old restart without the key | all: `lag` | `src/rad_m1/rad_m1.cpp:537 RadiationM1 constructor` |
| `rad_m1/clean_stag_xdonor_freeze` | From this pass of a solve on, the xdonor cross flux is held fixed (its donor side and cap no longer switch), which stops period-2 pass cycles (0 = off). | int | 4 | all: `4` | `src/rad_m1/rad_m1.cpp:523 RadiationM1 constructor` |
| `rad_m1/clean_stag_xsrc` | With clean_stag_xdonor: the face sources beyond their explicit share are also in donor form. | bool | true with xdonor (fresh run; a restart keeps the value stored in the file) | all: `true` | `src/rad_m1/rad_m1.cpp:546 RadiationM1 constructor` |
| `rad_m1/clean_thin_donor` | The thin share of every face from its donor cell alone (no lagged ray flux). | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:2207 RadiationM1 constructor` |
| `rad_m1/clean_w_c0` | Lower edge of the thin-weight ramp (w = 0 below). | real | 0.3 with clean_w_fse (0.1 without) | all: `0.3` | `src/rad_m1/rad_m1.cpp:103 RadiationM1 constructor` |
| `rad_m1/clean_w_c1` | Upper edge of the thin-weight ramp (w = 1 above). | real | 0.6 with clean_w_fse (0.5 without) | all: `0.6` | `src/rad_m1/rad_m1.cpp:104 RadiationM1 constructor` |
| `rad_m1/clean_w_face_min` | A face's thin weight from the thinner of its two cells (default: the thicker). | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:2219 RadiationM1 constructor` |
| `rad_m1/clean_w_fse` | The thin (free-streaming) weight w of a face from the field's anisotropy g = max(\|H\|/J, eps C1/1.5), eps = min(kappa_P/chi, 1), w = smoothstep between clean_w_c0 and clean_w_c1. | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:99 RadiationM1 constructor` |
| `rad_m1/clean_w_fse_tau` | Scattering isotropises only over the cell's optical depth: eps -> 1 - (1 - eps)(1 - exp(-tau_cell)). Must be named true for restarts from files that stored false. | bool | true (with clean_w_fse) | all: `true` | `src/rad_m1/rad_m1.cpp:102 RadiationM1 constructor` |

## Implicit solver

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_bicg_stall_fac` | BiCGStab stall test factor. | real | 0.99 | all: `0.99` | `src/rad_m1/rad_m1.cpp:285 RadiationM1 constructor` |
| `rad_m1/clean_bicg_stall_its` | BiCGStab stall window (0 = off). | int | 0 with clean_lin_pass (50 without) | all: `0` | `src/rad_m1/rad_m1.cpp:283 RadiationM1 constructor` |
| `rad_m1/clean_conv` | Convergence measure: budget = each row relative to the cell's own energy budget of the step, tolerance clean_conv_kappa x the time-error estimate, between clean_tol_nl and clean_conv_max; strict = E, T relative to themselves at clean_tol_nl. | string | budget | all: `budget` | `src/rad_m1/rad_m1.cpp:926 RadiationM1 constructor` |
| `rad_m1/clean_conv_est_f` | Also stop on the convergence-rate estimate of the face F changes (budget passes only). | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:908 RadiationM1 constructor` |
| `rad_m1/clean_conv_kappa` | Factor on the time-error estimate that sets the budget tolerance. | real | 1.0 | all: `1.0` | `src/rad_m1/rad_m1.cpp:933 RadiationM1 constructor` |
| `rad_m1/clean_conv_max` | Ceiling of the budget tolerance. The MeshBlock layout changes the answer at about this level. | real | 1.0e-5 | all: `1.0e-5` | `src/rad_m1/rad_m1.cpp:153 RadiationM1 constructor` |
| `rad_m1/clean_gas_newton_switch` | Stall detection window W (0 = off; 3..64): a stalled solve takes its remaining passes on the bracketed gas temperature root without the Newton slope. | int | 4 (fresh run; a restart keeps the value stored in the file); 0 on a restart without the key | all: `4` | `src/rad_m1/rad_m1.cpp:761 RadiationM1 constructor` |
| `rad_m1/clean_gas_newton_switch_min` | First pass at which the stall switch may act. | int | 6 (fresh run; a restart keeps the value stored in the file); 20 on a restart without the key | all: `6` | `src/rad_m1/rad_m1.cpp:762 RadiationM1 constructor` |
| `rad_m1/clean_init_cnorm` | Per-cell convergence test of the steady start's solves (0 = off). | real | 1.0e-10 | all: `1e-10` | `src/rad_m1/rad_m1.cpp:129 RadiationM1 constructor` |
| `rad_m1/clean_lin_cnorm` | Per-cell residual test of the steps' strict solves (0 = off). | real | 1.0e-10 | all: `1e-10` | `src/rad_m1/rad_m1.cpp:130 RadiationM1 constructor` |
| `rad_m1/clean_lin_pass` | A linear solve that ends at clean_maxit leaves its iterate to the next pass (no stall test). | bool | true (fresh run; a restart keeps the value stored in the file) | all: `true` | `src/rad_m1/rad_m1.cpp:164 RadiationM1 constructor` |
| `rad_m1/clean_maxit` | Maximum BiCGStab iterations per linear solve. Production: 200. | int | 800 | all: `200` | `src/rad_m1/rad_m1.cpp:152 RadiationM1 constructor` |
| `rad_m1/clean_maxnl` | Maximum Newton/Picard passes per solve. | int | 200 | all: `200` | `src/rad_m1/rad_m1.cpp:151 RadiationM1 constructor` |
| `rad_m1/clean_maxpic` | Passes that update the lagged slopes. | int | 12 | all: `12` | `src/rad_m1/rad_m1.cpp:388 RadiationM1 constructor` |
| `rad_m1/clean_nl_norm` | Newton change norm: global or cell (per-cell relative). | string | global | all: `global` | `src/rad_m1/rad_m1.cpp:888 RadiationM1 constructor` |
| `rad_m1/clean_pc_sweep` | Smoother of the multigrid preconditioner: rb (red-black), rbr, sym, rbrsym, wave. Production: rb. | string | sym (fresh run; a restart keeps the value stored in the file); rb on a restart without the key | all: `rb` | `src/rad_m1/rad_m1.cpp:818 RadiationM1 constructor` |
| `rad_m1/clean_retry_budget` | A step redone by dt_retry keeps the budget tolerance instead of the strict one. | bool | false | all: `0` | `src/rad_m1/rad_m1.cpp:942 RadiationM1 constructor` |
| `rad_m1/clean_stag_fs_line` | Line solver of the face F solve along x1 on a GPU: pcr or thomas. | string | pcr (fresh run; a restart keeps the value stored in the file) | all: `pcr` | `src/rad_m1/rad_m1.cpp:275 RadiationM1 constructor` |
| `rad_m1/clean_stag_mg_color` | Line colouring of the multigrid preconditioner: mirror (mirror-symmetric) or parity. | string | mirror (fresh run; a restart keeps the value stored in the file); parity on a restart without the key | all: `mirror` | `src/rad_m1/rad_m1.cpp:780 RadiationM1 constructor` |
| `rad_m1/clean_stag_opac_drag` | An opacity refresh also rebuilds the implicit drag share of the face momentum. Production: 0 (false). | bool | true (fresh run; a restart keeps the value stored in the file); false on a restart without the key | all: `0` | `src/rad_m1/rad_m1.cpp:254 RadiationM1 constructor` |
| `rad_m1/clean_stag_opac_passes` | Number of passes whose opacities are updated (second order in time for N >= 2). | int | 3 | all: `3` | `src/rad_m1/rad_m1.cpp:249 RadiationM1 constructor` |
| `rad_m1/clean_stag_opac_update` | Opacities re-evaluated at the iterate's temperature during the solve (for clean_stag_opac_passes passes, then frozen). | bool | true (fresh run; a restart keeps the value stored in the file); false on a restart without the key | all: `true` | `src/rad_m1/rad_m1.cpp:246 RadiationM1 constructor` |
| `rad_m1/clean_stag_pcr` | Preconditioner line solves by parallel cyclic reduction on a GPU. | bool | true | all: `1` | `src/rad_m1/rad_m1.cpp:234 RadiationM1 constructor` |
| `rad_m1/clean_stall_fac` | Stall criterion: the window's smallest change is not below this fraction of the previous window's. | real | 0.5 | all: `0.5` | `src/rad_m1/rad_m1.cpp:763 RadiationM1 constructor` |
| `rad_m1/clean_tol` | Relative residual target of the linear (BiCGStab) solve. | real | 1.0e-10 | all: `1e-10` | `src/rad_m1/rad_m1.cpp:128 RadiationM1 constructor` |
| `rad_m1/clean_tol_nl` | Lower bound of the Newton/Picard tolerance (relative changes of E, T, v). | real | 1.0e-11 | all: `1e-11` | `src/rad_m1/rad_m1.cpp:165 RadiationM1 constructor` |
| `rad_m1/clean_troot_rtol` | Residual stop of the gas-radiation exchange temperature root. | real | 1.0e-13 | all: `1e-13` | `src/rad_m1/rad_m1.cpp:896 RadiationM1 constructor` |
| `rad_m1/clean_xt_exact` | Exchange bookkeeping with one linearisation per cell and pass, so that the work between gas and radiation is exact (cons); x2 = also every temperature root from the local 2x2 system. | string | cons (fresh run; a restart keeps the value stored in the file); off on a restart without the key | all: `cons` | `src/rad_m1/rad_m1.cpp:661 RadiationM1 constructor` |
| `rad_m1/implicit_predictor` | The stage solve starts from the predicted E (step) rather than from the old state (none). | string | step on a fresh run (fixed closures); none on a restart without the key | all: `step` | `src/rad_m1/rad_m1.cpp:217 RadiationM1 constructor and src/rad_m1/rad_m1_init.cpp:1051 ImplicitInit` |
| `rad_m1/implicit_predictor_tcap` | Cap on the predicted temperature ratio: the predicted E is capped at tcap^4 times its start value. | real | 1.1 (read only when named; 2.0 with <job>/recipe = none) | all: `1.1` | `src/rad_m1/rad_m1_init.cpp:1820 RadiationM1::ImplicitInit (used in src/rad_m1/rad_m1_clean.cpp:2779)` |
| `rad_m1/implicit_resid_fatal` | Without dt_retry: a backward-Euler solve that ends not converged with a last change above this value stops the run. With dt_retry on (production) the cycle is retried instead. | real | 1e-2 on a fresh run; 0 = off on a restart without the key; 1e2 with <job>/recipe = none | all: `1.0e-2` | `src/rad_m1/rad_m1_init.cpp:904 RadiationM1::ImplicitInit` |

## Positivity, floors and caps

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_wb_fct` | Zalesak limit on the energy flux of the faces written back, so that E > e_floor by construction. | bool | true (fresh run; a restart keeps the value stored in the file) (needs clean_gas_tfloor and clean_xt_exact); false on a restart without the key | all: `true` | `src/rad_m1/rad_m1.cpp:745 RadiationM1 constructor` |
| `rad_m1/cs_cross_mono` **[read, no effect here]** | Cubed sphere only: monotone cross terms at the seams. Inert on spherical polar. | bool | true (read only when named) | all: `true` | `src/rad_m1/rad_m1_init.cpp:1823 RadiationM1::ImplicitInit` |
| `rad_m1/cs_vertex_mono` **[read, no effect here]** | Cubed sphere only: monotone limiter at cube vertices. Inert on spherical polar. | bool | true on the cubed sphere (read only when named) | all: `true` | `src/rad_m1/rad_m1_init.cpp:1805 RadiationM1::ImplicitInit` |
| `rad_m1/e_floor` | Floor of the radiation energy density E [code units]. | real | FLT_MIN | all: `1.0e-4` | `src/rad_m1/rad_m1.cpp:1231 RadiationM1 constructor` |
| `rad_m1/floorcap_rho` | Floor-gas cap after the implicit solve, in active cells with rho < floorcap_rho [g/cm^3]: \|v\| capped at floorcap_v and T_gas at max(T_rad, floor); the removed kinetic and excess internal energy go to E_rad (energy conserved per cell, momentum removed). | real | off (read only when named) | all: `1e-14` | `src/rad_m1/rad_m1.cpp:1617 RadiationM1 constructor` |
| `rad_m1/floorcap_v` | Velocity cap of the floor-gas cap [cm/s]. | real | 0 = no velocity cap (read only with floorcap_rho) | all: `3e7` | `src/rad_m1/rad_m1.cpp:1619 RadiationM1 constructor` |
| `rad_m1/implicit_pos_floor` | Radiation E kept at or above e_floor after the solve (raises counted). | bool | true on a fresh run; old value (off) on a restart without the key | all: `true` | `src/rad_m1/rad_m1_init.cpp:1180 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_pos_floor_cons` | The E floor's raise is taken from the radiation (then gas) of the same column instead of being created. | bool | true (read only when named) | all: `true` | `src/rad_m1/rad_m1_init.cpp:1839 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_pos_floor_gas_keep` | Fraction of its internal energy the gas keeps when it pays an E floor raise. | real | 1.0 on a fresh clean run | all: `1` | `src/rad_m1/rad_m1_init.cpp:1853 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_pos_floor_solve` **[read, no effect here]** | Charge the solve undershoot of the floor raise. Under scheme = clean the code only checks it against implicit_pos_floor (src/rad_m1/rad_m1_init.cpp:1234); no other use was found. | bool | true with implicit_pos_floor on a fresh run | all: `true` | `src/rad_m1/rad_m1_init.cpp:1201 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_pos_gas` | Gas internal energy kept positive in the coupled solve (energy moved from radiation to gas, counted in the m1-positivity report). | bool | true on a fresh run; old value (off) on a restart without the key | all: `true` | `src/rad_m1/rad_m1_init.cpp:1179 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_pos_gas_cons` | The energy implicit_pos_gas needs is taken from the radiation of the same radial column (energy conserved). | bool | true on a fresh clean run; off on a restart without the key | all: `true` | `src/rad_m1/rad_m1_init.cpp:1850 RadiationM1::ImplicitInit` |
| `rad_m1/time2_face_fct` | Zalesak-type flux-corrected limiter on the stage increment's face fluxes, so that E stays positive; the solve stores the energy flux that moved E. | bool | true on a fresh clean run; off on a restart without the key | all: `true` | `src/rad_m1/rad_m1_time2.cpp:585 RadiationM1::Time2Init` |
| `rad_m1/zalesak_active` | The active-set form of the Zalesak limiter rounds. | bool | true on a fresh clean run | all: `true` | `src/rad_m1/rad_m1_init.cpp:2584 RadiationM1::ImplicitInit` |
| `hydro/dfloor` | Density floor [g/cm^3]. Mass created by this floor is NOT in any history flux column (see OUTPUTS.md, energy ledger). | real | FLT_MIN | bsg: `6.15245e-14`; m35mms: `1.0e-17`; t9l62: `3.6e-19` | `src/eos/eos.cpp:32 EquationOfState constructor` |
| `hydro/dfloor_keep_velocity` | When the density floor acts, keep the velocity (the momentum is rescaled; its kinetic energy change is in the event log column dfloor_ke). | bool | false | all: `true` | `src/eos/eos.cpp:34 EquationOfState constructor` |
| `hydro/floor_pay` | Who pays the internal energy the hydro floors create: rad_column = the radiation of the cell, then of its radial column; what is still unpaid is created and reported (CLEAN WARNING). The potential energy of mass created by dfloor is not paid. | string | rad_column on a fresh run with eos = general and <rad_m1>; none otherwise | all: `rad_column` | `src/hydro/dual_energy.cpp:53 DualEnergy::Read` |
| `hydro/pfloor` | Pressure floor [code units]. | real | FLT_MIN | bsg: `8.9014e-3`; m35mms: `1.0e-12`; t9l62: `1.0e-12` | `src/eos/eos.cpp:258 EquationOfState constructor` |
| `hydro/tfloor` | Gas temperature floor in code units (the decks: 1000 K / temp_unit_kelvin). The clean gas row also uses it as the lower bound of T. | real | FLT_MIN | bsg: `134104228155.32465`; m35mms: `133780565496.86447`; t9l62: `123359972487.0939` | `src/eos/eos.cpp:259 EquationOfState constructor` |
| `hydro/vceil` | Velocity ceiling [cm/s]; the kinetic energy removed is reported in the event log (vceil_de). | real | 0 = off | bsg: `1.0e8`; m35mms: `3.0e8`; t9l62: `3.0e8` | `src/eos/eos.cpp:43 EquationOfState constructor` |

## Well-balancing, radiation force and angular momentum

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_cs_lz_tel` **[read, no effect here]** | Cubed sphere only: L_z by telescoping. On spherical polar it must be false (true is FATAL); the decks name false. | bool | true on a fresh staggered cubed-sphere run (read only when named elsewhere) | all: `false` | `src/rad_m1/rad_m1.cpp:601 RadiationM1 constructor` |
| `rad_m1/clean_stag_fwb` | Face F rows in deviation form about the rays' own state, so that the rays' state is an exact steady state of the face equations (radiation well-balancing). | bool | true (fresh run; a restart keeps the value stored in the file) with the F solve and clean_consistency; false on a restart without the key | all: `true` | `src/rad_m1/rad_m1.cpp:445 RadiationM1 constructor` |
| `rad_m1/clean_stag_fwb_ring` | The fwb residual of the phi faces minus its ring mean: no torque on a ring, L_z exact without a ledger. Production: 0 (false). | bool | true (fresh run; a restart keeps the value stored in the file) with fwb on a uniform 3-D spherical-polar mesh with periodic x3 | all: `0` | `src/rad_m1/rad_m1.cpp:481 RadiationM1 constructor` |
| `rad_m1/force_reference` | wb_arad: the gas hydrostatic balance includes a reference radiation acceleration a_ref = kappa F_r/c of the initial column, delivered through the well-balanced effective potential; the radiation module hands the gas only the residual force. | string | none | all: `wb_arad` | `src/rad_m1/rad_m1.cpp:1710 RadiationM1 constructor` |
| `rad_m1/force_reference_work` | How the work of the reference force is paid: split (by the problem generator) or full. | string | auto (= split when the problem generator supplies the reference work, as he_star_m1 does) | all: `split` | `src/rad_m1/rad_m1.cpp:1735 RadiationM1 constructor` |
| `rad_m1/time2_fref_wacc` | With force_reference = wb_arad: the exact payment of the reference force work also on the backward-Euler steps of a hesdirk2 run (false creates energy on those steps). | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:1767 RadiationM1 constructor` |
| `hydro/etotgrav` | The conserved energy includes the gravitational potential energy rho Phi (total energy equation). Phi = GM(1/r_in - 1/r) >= 0 in he_star_m1. | bool | true on a fresh run where the problem generator supports it (he_star_m1 does) | all: `true` | `src/hydro/hydro.cpp:257 Hydro constructor` |
| `hydro/wall_wb_onesided` | One-sided well-balanced reconstruction at closed/reflecting x1 walls. | bool | true with wellbalance_dynamic + wb_x1 (fresh and restart) | all: `true` | `src/hydro/hydro.cpp:338 Hydro constructor` |
| `hydro/wb_cache_every` **[never read]** | NOT READ by this build. | - | (no key) | all: `0` | (no read in `src/`) |
| `hydro/wb_option` | Background assumed by the well-balanced reconstruction: isothermal, isodensity, polytropic, ... | string | polytropic on a fresh sp/cs star start; otherwise required with wellbalance_dynamic | all: `polytropic` | `src/hydro/hydro.cpp:482 Hydro constructor` |
| `hydro/wb_rmax` | Radius [cm] above which the x1 well-balanced reconstruction is switched off. | real | 0 = WB everywhere | bsg: `4.231839e+12`; m35mms: `1.097072e+12`; t9l62: `7.504290e+12` | `src/hydro/hydro.cpp:353 Hydro constructor` |
| `hydro/wb_x1` | Apply the well-balanced reconstruction along x1 (radial). | bool | true on a fresh sp/cs star start | all: `true` | `src/hydro/hydro.cpp:307 Hydro constructor` |
| `hydro/wellbalance_dynamic` | Well-balanced hydro reconstruction about a hydrostatic background. | bool | true on a fresh sp/cs star start (problem early defaults); false otherwise | all: `true` | `src/hydro/hydro.cpp:306 Hydro constructor` |

## Layout, reproducibility and parallel exchange

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_halo_direct` | Direct (list-based) ghost exchange of the clean solve. | bool | true | all: `1` | `src/rad_m1/rad_m1.cpp:948 RadiationM1 constructor` |
| `rad_m1/clean_rst_mshare` | At a restart from a file without the clean state, the memory of the M1 share starts at its steady flux of the first solution. Acts only on such restarts. | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:702 RadiationM1 constructor` |
| `rad_m1/clean_stag_tie_gid` | A face shared by two MeshBlocks takes T, rho and Phi of the lower-gid block (the cubed-sphere rule). Production: 0 (false). | bool | true on a fresh run on a uniform mesh; false on a restart without the key | all: `0` | `src/rad_m1/rad_m1.cpp:299 RadiationM1 constructor` |
| `rad_m1/vet_gd_halo_compact` | Compact ghost exchange of the ray sweep (bitwise the same; faster). | int | 2 on a fresh vet_gd start; 0 on a restart without the key | all: `2` | `src/rad_m1/rad_m1.cpp:2373 RadiationM1 constructor` |
| `rad_m1/vet_gd_halo_ipc` | Exchange the ray-sweep halos between GPUs of one node by CUDA IPC (speed only, bitwise the same). bsg names true (single node, 2 GPUs), m35mms/t9l62 false. | bool | on for clean on CUDA without a theta wall; off on other builds | bsg: `true`; m35mms: `false`; t9l62: `false` | `src/rad_m1/rad_m1.cpp:2530 RadiationM1 constructor` |
| `mesh/layout_invariant` | Sums (reductions, history, shell sums) in a fixed order, so the result does not depend on the MPI rank count. | bool | true (fresh and restart; not for an SMR/AMR restart) | all: `true` | `src/mesh/mesh.cpp:175 Mesh constructor` |
| `mesh/layout_invariant_mode` | block = per-MeshBlock sums in a fixed order (results the same for any rank count with the same MeshBlocks); exact = every cell term exact (+8 percent cost); tile (with layout_tile = N) = block mode on fixed N x N lateral tiles, so the result is also the same for every MeshBlock layout made of whole tiles. Recommended for production: tile with layout_tile = the lateral MeshBlock size (128 here). | string | block | bsg: not named; m35mms: `block`; t9l62: `block` | `src/mesh/mesh.cpp:186 Mesh constructor` |
| `meshblock/nx1` | MeshBlock size along x1. Must equal mesh/nx1 (one MeshBlock along x1) for the ray cut, the history face picks and layout_invariant_mode = tile. | int | = mesh/nx1 | bsg: `256`; m35mms: `336`; t9l62: `512` | `src/mesh/mesh.cpp:870 Mesh constructor` |
| `meshblock/nx2` | MeshBlock size along x2. With layout_invariant_mode = tile, layout_tile must divide it. | int | = mesh/nx2 | all: `128` | `src/mesh/mesh.cpp:872 Mesh constructor` |
| `meshblock/nx3` | MeshBlock size along x3. With tile mode, layout_tile must divide it. | int | = mesh/nx3 | all: `128` | `src/mesh/mesh.cpp:877 Mesh constructor` |
| `time/dt_retry_host` | Where the checkpoint copies used by dt_retry live: false = device memory; true = pinned host memory (saves GPU memory, costs a copy per cycle); m1 = only the radiation copies on the host. Bitwise the same result. | bool or m1 | false (read only when named) | all: `false` | `src/driver/dt_retry.cpp:515 DtRetry constructor` |
| `time/dt_retry_views` | Which arrays dt_retry checkpoints: minimal, all, or clean (minimal without arrays the clean scheme never writes). | string | minimal | all: `minimal` | `src/driver/dt_retry.cpp:491 DtRetry constructor` |
| `time/restart_refill_ghosts` | false: a restart keeps the ghost cells stored in the restart file, so the restart is a bitwise continuation of the straight run. true refills them (not bitwise). | bool | false | all: `false` | `src/driver/driver.cpp:382 Driver::Initialize` |

## Boundaries

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_bc` | Open x1 end of the clean scheme: halfrange = the end cell's ray half-range flux; marshak = the Marshak form. | string | halfrange | all: `halfrange` | `src/rad_m1/rad_m1.cpp:223 RadiationM1 constructor` |
| `rad_m1/clean_bc_area` | The open end's half-range flux in area form (net outward stream h+ + h- A_c/A_f), so an isotropic field gives no net flux. | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:124 RadiationM1 constructor` |
| `rad_m1/implicit_bc_x1max` | Radiation outer x1 boundary: marshak = an open top with an incoming bath of implicit_flux_x1max (0 = vacuum). | string | auto | all: `marshak` | `src/rad_m1/rad_m1_init.cpp:2003 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_bc_x1min` | Radiation inner x1 boundary: flux = a given inward flux implicit_flux_x1min (the stellar luminosity entering the base). | string | auto (from the mesh boundary) | all: `flux` | `src/rad_m1/rad_m1_init.cpp:2000 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_flux_x1max` | Incoming flux/bath at the outer Marshak face (0 = vacuum). | real | 0.0 | all: `0.0` | `src/rad_m1/rad_m1_init.cpp:2092 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_flux_x1min` | Radiation flux [erg/cm^2/s] through the inner face = L/(4 pi r_in^2). he_star_m1 also reads it for the initial column. | real | 0.0 | bsg: `2.2942145641e+13`; m35mms: `1.8587122274e+14`; t9l62: `6.9466003350e+13` | `src/rad_m1/rad_m1_init.cpp:2039 RadiationM1::ImplicitInit` |
| `rad_m1/marshak_q` | Marshak boundary factor q (outgoing flux = q c E at the open face, old form). | real | 0.5 | all: `0.5` | `src/rad_m1/rad_m1_init.cpp:1402 RadiationM1::ImplicitInit` |
| `rad_m1/vet_bc_diff` | Inner diffusion intensity of the formal solution taken at the face (face) or from the first cell (cell). | string | face on a fresh run; cell on a restart without the key | all: `face` | `src/rad_m1/rad_m1.cpp:1966 RadiationM1 constructor` |
| `rad_m1/vet_bc_inner` | Intensity entering the rays at the inner x1 face: diff (diffusion intensity), emit, or auto (emit at a Marshak face, diff at a flux face). Here the inner face is a flux face, so auto = diff. | string | auto | all: `auto` | `src/rad_m1/rad_m1.cpp:3143 RadiationM1::VetBcInit` |
| `rad_m1/vet_bc_reflect` **[read, no effect here]** | Rays at a reflecting x1 end. Inert here (no reflecting x1 end). | string | auto (a restart without the key: open) | all: `auto` | `src/rad_m1/rad_m1.cpp:3163 RadiationM1::VetBcInit` |
| `mesh/ix1_bc` | Inner x1 boundary of the hydro: user = the problem generator's (he_star_m1: closed wall). | string | none (required) | all: `user` | `src/mesh/mesh.cpp:586 Mesh constructor` |
| `mesh/ix2_bc` | Inner x2 (theta) boundary; periodic for these wedges. | string | none (required) | all: `periodic` | `src/mesh/mesh.cpp:625 Mesh constructor` |
| `mesh/ix3_bc` | Inner x3 (phi) boundary; periodic. | string | none (required) | all: `periodic` | `src/mesh/mesh.cpp:660 Mesh constructor` |
| `mesh/ox1_bc` | Outer x1 boundary: user = the problem generator's (he_star_m1: outflow). | string | none (required) | all: `user` | `src/mesh/mesh.cpp:587 Mesh constructor` |
| `mesh/ox2_bc` | Outer x2 boundary; periodic. | string | none (required) | all: `periodic` | `src/mesh/mesh.cpp:626 Mesh constructor` |
| `mesh/ox3_bc` | Outer x3 boundary; periodic. | string | none (required) | all: `periodic` | `src/mesh/mesh.cpp:661 Mesh constructor` |

## Opacity, EOS, units and coupling

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/arad` | Radiation constant a in code units: equilibrium E = arad T^4 with T the EOS code temperature (here a_cgs times T_unit^4). | real | required when the gas and radiation are coupled; 0 otherwise | bsg: `2.3392742567e-47`; m35mms: `2.3619941258e-47`; t9l62: `3.2670356041e-47` | `src/rad_m1/rad_m1.cpp:1606 RadiationM1 constructor` |
| `rad_m1/c_light` | Speed of light in code units (here cgs: 2.99792458e10 cm/s). | real | none (required) | all: `2.99792458e10` | `src/rad_m1/rad_m1.cpp:50 RadiationM1 constructor` |
| `rad_m1/chat_over_c` | Reduced speed of light chat/c. The clean scheme needs 1.0 (FATAL otherwise). | real | 1.0 | all: `1.0` | `src/rad_m1/rad_m1.cpp:56 RadiationM1 constructor` |
| `rad_m1/coupling` | Emission/absorption coupling between gas and radiation. | bool | true when any opacity is non-zero | all: `true` | `src/rad_m1/rad_m1.cpp:1528 RadiationM1 constructor` |
| `rad_m1/gas_feedback` | The radiation energy and momentum exchange acts back on the gas (radiation force and heating). | bool | true | all: `true` | `src/rad_m1/rad_m1.cpp:1529 RadiationM1 constructor` |
| `rad_m1/implicit_enthalpy` | Face value of the advective radiation enthalpy flux (4/3) v E: upwind (donor cell) or plm (limited linear). The matrix always keeps the donor-cell part. | string | plm for the fixed-tensor closures on a fresh run; upwind on a restart without the key | all: `plm` | `src/rad_m1/rad_m1_init.cpp:1478 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_gas_src` | Form of the gas source of the coupled solve. conservative2 = the conservative form (energy exact). | string | conservative2 on a fresh coupled run; matched on a restart without the key | all: `conservative2` | `src/rad_m1/rad_m1_init.cpp:1883 RadiationM1::ImplicitInit` |
| `rad_m1/kappa_unit` | Code opacity per (cm^2/g) of the tables. | real | 1.0 | all: `1.0` | `src/rad_m1/rad_m1.cpp:1295 RadiationM1 constructor` |
| `rad_m1/opacity` | Opacity source: const, powerlaw, user or table. Production: table (Rosseland and Planck tables given in <problem>). | string | const | all: `table` | `src/rad_m1/rad_m1.cpp:1273 RadiationM1 constructor` |
| `rad_m1/rho_unit_cgs` | g/cm^3 per code density unit, used for the table lookup. | real | 1.0 | all: `1.0` | `src/rad_m1/rad_m1.cpp:1294 RadiationM1 constructor` |
| `rad_m1/scheme` | Radiation solver. Only clean exists in this build; any other value is FATAL. | string | clean with transport = implicit on a fresh run | all: `clean` | `src/rad_m1/rad_m1.cpp:89 RadiationM1 constructor` |
| `rad_m1/temp_unit_kelvin` | Kelvin per code temperature unit, used for the table lookup: T[K] = T_code x temp_unit_kelvin. With an ideal gas it must equal mu m_u / k_B of <units>/mu. | real | 1.0 | bsg: `7.4568864364e-9`; m35mms: `7.4749268633e-09`; t9l62: `8.1063567271e-09` | `src/rad_m1/rad_m1.cpp:1293 RadiationM1 constructor` |
| `rad_m1/transport` | Radiation transport: explicit, implicit_x1 or implicit. The clean scheme needs implicit. | string | implicit on a fresh run with a <rad_m1> block (recipe::CleanDefaults); explicit otherwise (a restart file without the key) | all: `implicit` | `src/rad_m1/rad_m1.cpp:67 RadiationM1 constructor` |
| `hydro/eos` | Equation of state: general = the general EOS interface (here with general_eos = gamma, an ideal gamma law through the general code path). Gas only: the radiation is not in the EOS. | string | required | all: `general` | `src/hydro/hydro.cpp:60 Hydro constructor` |
| `hydro/gamma` | Adiabatic index (5/3). | real | required | all: `1.6666666666666667` | `src/eos/general_hyd.cpp:30 GeneralHydro constructor` |
| `hydro/general_eos` | Which general EOS: gamma (ideal gamma law) or table. | string | gamma | all: `gamma` | `src/eos/eos.cpp:364 EquationOfState::BuildGeneralEOS` |

## Mesh and grid

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `mesh/cs_seam_blockend` **[never read]** | Cubed sphere only. NOT READ on spherical polar (the reader returns before it). | string | extend / clamp | all: `extend` | `src/bvals/bvals.cpp:345 CsSeamBlockendExtend` |
| `mesh/cs_seam_rho_guard` **[never read]** | Cubed sphere only. NOT READ on spherical polar. | real | 4.0 | all: `4` | `src/mesh/mesh.cpp:273 Mesh constructor` |
| `mesh/cs_vertex_3cell` **[never read]** | Cubed sphere only. NOT READ on spherical polar. | bool | true on a fresh cs start | all: `true` | `src/mesh/mesh.cpp:249 Mesh constructor` |
| `mesh/f_stretch_r` | Radial stretch parameter of the exponential map. | real | required with use_grid_stretch_r | bsg: not named; m35mms: `-0.455726942889053`; t9l62: `-3.509253962786431` | `src/mesh/mesh.cpp:300 Mesh constructor` |
| `mesh/f_stretch_r_p_amp` | Amplitude of the radial plateau of the polynomial stretch (bsg: a factor-2 refined band). | real | 0 (plateau read only when its amplitude is given) | bsg: `-0.66666666666666667`; m35mms: not named; t9l62: not named | `src/mesh/mesh.cpp:341 Mesh constructor (key built as f_stretch_r_p + _amp)` |
| `mesh/f_stretch_r_p_w` | Edge width of the plateau (1e-12 = sharp edges). | real | required with _amp | bsg: `1.0e-12`; m35mms: not named; t9l62: not named | `src/mesh/mesh.cpp:341 Mesh constructor` |
| `mesh/f_stretch_r_p_xa` | Start of the plateau in the unit coordinate (0..1). | real | required with _amp | bsg: `0.1875`; m35mms: not named; t9l62: not named | `src/mesh/mesh.cpp:341 Mesh constructor` |
| `mesh/f_stretch_r_p_xb` | End of the plateau in the unit coordinate. | real | required with _amp | bsg: `0.6875`; m35mms: not named; t9l62: not named | `src/mesh/mesh.cpp:341 Mesh constructor` |
| `mesh/nghost` | Number of ghost cells (the clean F solve across MeshBlocks needs >= 3). | int | 2 | all: `3` | `src/mesh/mesh.cpp:119 Mesh constructor` |
| `mesh/nx1` | Number of radial cells of the mesh. | int | required | bsg: `256`; m35mms: `336`; t9l62: `512` | `src/mesh/mesh.cpp:120 Mesh constructor` |
| `mesh/nx2` | Number of theta cells. | int | required | all: `256` | `src/mesh/mesh.cpp:121 Mesh constructor` |
| `mesh/nx3` | Number of phi cells. | int | required | bsg: `256`; m35mms: `256`; t9l62: `512` | `src/mesh/mesh.cpp:122 Mesh constructor` |
| `mesh/use_grid_stretch_r` | Exponentially stretched radial grid with factor f_stretch_r (excludes use_grid_stretch_r_poly). | bool | false | bsg: not named; m35mms: `true`; t9l62: `true` | `src/mesh/mesh.cpp:299 Mesh constructor` |
| `mesh/use_grid_stretch_r_poly` | Polynomial radial stretch (can make non-monotonic spacing, e.g. a refined radial band). | bool | false | bsg: `true`; m35mms: not named; t9l62: not named | `src/mesh/mesh.cpp:305 Mesh constructor` |
| `mesh/use_spherical_polar` | Spherical-polar coordinates (r, theta, phi). | bool | false | all: `true` | `src/mesh/mesh.cpp:298 Mesh constructor` |
| `mesh/x1max` | Outer radius [cm]. | real | required | bsg: `5.565600e12`; m35mms: `1.0641150000e+12`; t9l62: `5.6321590000e+13` | `src/mesh/mesh.cpp:113 Mesh constructor` |
| `mesh/x1min` | Inner radius [cm]. | real | required | bsg: `1.391400e12`; m35mms: `6.7463500000e+11`; t9l62: `1.6850984345e+12` | `src/mesh/mesh.cpp:112 Mesh constructor` |
| `mesh/x2max` | Upper theta [rad]. | real | required | bsg: `2.0943951023931953`; m35mms: `1.74532925199432953`; t9l62: `2.3561944901923448` | `src/mesh/mesh.cpp:115 Mesh constructor` |
| `mesh/x2min` | Lower theta [rad]. | real | required | bsg: `1.0471975511965976`; m35mms: `1.39626340159546358`; t9l62: `0.78539816339744831` | `src/mesh/mesh.cpp:114 Mesh constructor` |
| `mesh/x3max` | Upper phi [rad]; with the product ray set 2 pi/vet_gd_gl_nphi must divide the phi extent. | real | required | bsg: `1.0471975511965976`; m35mms: `0.34906585039886590`; t9l62: `3.1415926535897931` | `src/mesh/mesh.cpp:117 Mesh constructor` |
| `mesh/x3min` | Lower phi [rad]. | real | required | all: `0.0` | `src/mesh/mesh.cpp:116 Mesh constructor` |

## Hydro options that matter for the coupling

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `hydro/am_transport` | Angular-momentum module; it also books the hydro x1 boundary fluxes that fill the history columns M_hflx, E_hflx, Lz_hflx. | bool | true on a fresh dynamic sp (multi-D) or cs run | all: `true` | `src/hydro/hydro_am.cpp:74 AmModule::DefaultOn` |
| `hydro/fofc` | First-order flux correction where a cell would hit a floor. Production: false. | bool | true for he_star_m1 on a fresh run (problem default); false on a restart without the key | all: `false` | `src/hydro/hydro.cpp:604 Hydro constructor` |
| `hydro/merge_lat_f` | Merge factor f of each radial range (1..32). | list of ints | required with merge_lat_rc | bsg: `2,1,2`; m35mms: not named; t9l62: not named | `src/hydro/hydro_mergelat.cpp:80 MergeLat constructor` |
| `hydro/merge_lat_mom` **[never read]** | NOT READ: the merged momentum is always the rigid (P and L conserving) form. | - | (no key) | bsg: `rigid`; m35mms: not named; t9l62: not named | (no read in `src/`) |
| `hydro/merge_lat_rc` | Merged lateral cells: radii [cm] below which f x f lateral cell groups are treated as one hydro cell (lifts the lateral CFL limit near r_in). bsg only. | list of reals | off (read only when named) | bsg: `2.43495e12,3.82635e12,5.5656e12`; m35mms: not named; t9l62: not named | `src/hydro/hydro_mergelat.cpp:79 MergeLat constructor` |
| `hydro/reconstruct` | Hydro reconstruction. | string | plm | all: `plm` | `src/hydro/hydro.cpp:608 Hydro constructor` |
| `hydro/rsolver` | Hydro Riemann solver (hllc). | string | required | all: `hllc` | `src/hydro/hydro.cpp:654 Hydro constructor` |
| `hydro/sp_x2_periodic_image` | Theta-periodic wedge: x2 reconstruction with the periodic image (seam mass/energy fix). | bool | true on a fresh theta-periodic sp mesh | all: `true` | `src/coordinates/coordinates.cpp:59 Coordinates constructor` |
| `hydro/split_c2p_once` | One conserved-to-primitive conversion after the operator-split sources. | bool | true (fresh and restart) | all: `true` | `src/hydro/hydro.cpp:282 Hydro constructor` |

## Output and diagnostics

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/floorcap_mesh_vol` | Weight the floorcap counters by the mesh cell volume and report \|dp\| and L_z (diagnostic only; the state is the same). Production: 0. | bool | true (fresh run; a restart keeps the value stored in the file); false on a restart without the key | all: `0` | `src/rad_m1/rad_m1.cpp:1652 RadiationM1 constructor` |
| `rad_m1/floorcap_report` | Print the cumulative floorcap counters every N cycles ("<rad_m1> floorcap cycle ..." line). | int | 0 = never (read only when named) | all: `25` | `src/rad_m1/rad_m1.cpp:1641 RadiationM1 constructor` |
| `hydro/merge_lat_diag` | Print merged-cell diagnostics every N cycles. | int | read only when named | bsg: `1000`; m35mms: not named; t9l62: not named | `src/hydro/hydro_mergelat.cpp:75 MergeLat constructor` |
| `time/nan_check_cycles` | Every N cycles the driver checks the state for NaN and stops the run if it finds one (0 = off). | int | 100 | all: `100` | `src/driver/driver.cpp:129 Driver constructor` |

## Experimental and not recommended

These keys are in the production decks, but not as features to tune. Each is either
deprecated, a study option that production keeps off, a removed path that must stay off, or
an option that production names at its off value only so that a restart file which stored
another value cannot turn it back on. Keep the production values. Do not switch any of them
on without a dedicated test (energy and momentum conservation, a bitwise restart check, and a
comparison against the production key set).

| key | meaning | type | code default | production | read at |
|---|---|---|---|---|---|
| `rad_m1/clean_pass_damp` | From pass N of a solve on, a pass's change of E that turns back is damped. Deprecated: do not use (keep 0). | int | 0 = off | all: `0` | `src/rad_m1/rad_m1.cpp:527 RadiationM1 constructor` |
| `rad_m1/clean_stag_flim` | A face F above its donor's c E is dissipated in the stress and the solved F is bounded. Off by decision (fwb + fcone bound F instead and keep free streaming); named false so that restarts from files that stored true run without it. | bool | false | all: `false` | `src/rad_m1/rad_m1.cpp:438 RadiationM1 constructor` |
| `rad_m1/clean_stag_fwb_stage` | fwb's residual of the stage-old state, once per stage. Off: each solve's own E is used (stage-old E made a flux dipole at an opaque photosphere). Restart files of earlier production stored true: always name false. | bool | false | all: `false` | `src/rad_m1/rad_m1.cpp:457 RadiationM1 constructor` |
| `rad_m1/clean_stag_fwb_t2` | fwb_stage at the stage-time state (needs fwb_stage; no-op without it). Restart files stored true: always name false. | bool | true with fwb_stage on a fresh run, else false | all: `false` | `src/rad_m1/rad_m1.cpp:463 RadiationM1 constructor` |
| `rad_m1/clean_stag_kpat` | Stage increment's gas kick in Patankar form. Production: false (named off: it fails other stars' luminosity tests). One of the must-name keys (above). | bool | false (read only when named) | all: `false` | `src/rad_m1/rad_m1.cpp:573 RadiationM1 constructor (stag())` |
| `rad_m1/clean_stag_thinf_held` **[never read]** | Thin flux factor of a ray-shadowed donor blended with its own half range. Off by decision. NOT READ here: read only when clean_stag_thinf = true. | bool | false | all: `false` | `src/rad_m1/rad_m1.cpp:508 RadiationM1 constructor` |
| `rad_m1/implicit_lin_maxit` **[read, no effect here]** | Old solver's linear iteration limit. No effect on the clean solve (clean_maxit is used); a NOTE is printed if it differs and clean_maxit is not named. | int | 200 (old module) | all: `200` | `src/rad_m1/rad_m1_init.cpp:205 RadiationM1::ImplicitInit` |
| `rad_m1/implicit_maxit` **[read, no effect here]** | Old solver's Picard pass limit. No effect on the clean solve (clean_maxnl is used); a value different from clean_maxnl only prints a NOTE when clean_maxnl is not named. | int | 200 (old module) | all: `200` | `src/rad_m1/rad_m1_init.cpp:203 RadiationM1::ImplicitInit` |
| `rad_m1/time2_enth_be` **[read, no effect here]** | Old-solver option (stage enthalpy rule for backward-Euler solves). Read, but the flag is not used anywhere in this build. | bool | true on a fresh run | all: `true` | `src/rad_m1/rad_m1_time2.cpp:188 RadiationM1::Time2Init` |
| `rad_m1/time2_k1_redist` | Column-local redistribution of K1; the clean scheme does not redistribute. Named 0 to pin restarts. | real | 0 under scheme = clean (read only when named) | all: `0` | `src/rad_m1/rad_m1_time2.cpp:567 RadiationM1::Time2Init` |
| `rad_m1/time2_k1_stiff` | Move the exchange part of K1 in cells with relaxation k dt > s. Off in production. Must always be named 0: a restart file that stored 1.0 otherwise turns it on silently. | real | 0 = off (read only when named) | all: `0` | `src/rad_m1/rad_m1_time2.cpp:542 RadiationM1::Time2Init` |
| `rad_m1/time2_stage_local_be` **[never read]** | Removed path. The clean preflight (src/rad_m1/rad_m1_clean.cpp:434 CleanPreflight) stops the run if it is anything but false/0. Not otherwise read, so it is listed in the "never read" start-up warning. | bool | (removed) | all: `false` | `src/rad_m1/rad_m1_clean.cpp:434 RadiationM1::CleanPreflight (table)` |

## Keys the code does not read, or reads without effect

These keys are in at least one production deck. "Never read" keys are listed by the start-up check `### WARNING: N input parameter(s) were never read by the code` (printed once, after the first cycle). They are harmless; they are kept in the decks so that the decks match the runs made with them.

**Never read with these decks:**

- `hydro/merge_lat_mom` (bsg: `rigid`; m35mms: not named; t9l62: not named): NOT READ: the merged momentum is always the rigid (P and L conserving) form.
- `hydro/wb_cache_every` (all: `0`): NOT READ by this build.
- `mesh/cs_seam_blockend` (all: `extend`): Cubed sphere only. NOT READ on spherical polar (the reader returns before it).
- `mesh/cs_seam_rho_guard` (all: `4`): Cubed sphere only. NOT READ on spherical polar.
- `mesh/cs_vertex_3cell` (all: `true`): Cubed sphere only. NOT READ on spherical polar.
- `rad_m1/clean_stag_thinf_held` (all: `false`): Thin flux factor of a ray-shadowed donor blended with its own half range. Off by decision. NOT READ here: read only when clean_stag_thinf = true.
- `rad_m1/time2_stage_local_be` (all: `false`): Removed path. The clean preflight (src/rad_m1/rad_m1_clean.cpp:434 CleanPreflight) stops the run if it is anything but false/0. Not otherwise read, so it is listed in the "never read" start-up warning.
- `rad_m1/vet_gd_rotate_align` (all: `true`): Align the rotation clock with the time step. NOT READ here: read only when the set rotates on the time clock (vet_gd_rotate_dt > 0 or auto).
- `rad_m1/vet_rotate_avg_ema_change_form` (all: `linear`): Form of the EMA change test. NOT READ here: read only when vet_rotate_avg_ema >= 1.
- `time/dt_retry_max` (all: `8`): NOT READ. The maximum number of cuts per cycle is fixed at 8 in the code (src/driver/dt_retry.hpp:174, maxcut_); the deck value 8 matches it.

**Read, but no effect with these decks (or under `scheme = clean`):**

- `rad_m1/clean_cs_lz_tel` (all: `false`): Cubed sphere only: L_z by telescoping. On spherical polar it must be false (true is FATAL); the decks name false.
- `rad_m1/clean_fill_cs_remap` (all: `0`): Cubed sphere only: the fill takes the half ranges of a seam neighbour in the local face families. Inert on spherical polar.
- `rad_m1/clean_ring_jmin` (all: `1`): Ring J ratios are measured only where both sets hold J >= clean_ray_min E. Inert with clean_rot_ring = off.
- `rad_m1/clean_rot_fuse` (all: `true`): Sweep the live ring's two sets in one fused sweep (bitwise). Inert with clean_rot_ring = off.
- `rad_m1/clean_rot_ring_fill_rst` (all: `true`): At a restart from a file without the clean state, the ring's older sets are filled from a sweep of the restart state. Inert with clean_rot_ring = off.
- `rad_m1/cs_cross_mono` (all: `true`): Cubed sphere only: monotone cross terms at the seams. Inert on spherical polar.
- `rad_m1/cs_vertex_mono` (all: `true`): Cubed sphere only: monotone limiter at cube vertices. Inert on spherical polar.
- `rad_m1/implicit_lin_maxit` (all: `200`): Old solver's linear iteration limit. No effect on the clean solve (clean_maxit is used); a NOTE is printed if it differs and clean_maxit is not named.
- `rad_m1/implicit_maxit` (all: `200`): Old solver's Picard pass limit. No effect on the clean solve (clean_maxnl is used); a value different from clean_maxnl only prints a NOTE when clean_maxnl is not named.
- `rad_m1/implicit_pos_floor_solve` (all: `true`): Charge the solve undershoot of the floor raise. Under scheme = clean the code only checks it against implicit_pos_floor (src/rad_m1/rad_m1_init.cpp:1234); no other use was found.
- `rad_m1/time2_enth_be` (all: `true`): Old-solver option (stage enthalpy rule for backward-Euler solves). Read, but the flag is not used anywhere in this build.
- `rad_m1/vet_bc_reflect` (all: `auto`): Rays at a reflecting x1 end. Inert here (no reflecting x1 end).
- `rad_m1/vet_col_lat_taucut` (all: `30.0`): Optical-depth cut below which the lateral sweep is skipped. INERT under scheme = clean: the code sets the cut to +infinity (rays from the first shell, src/rad_m1/rad_m1.cpp:2596 area).
- `rad_m1/vet_gd_ls` (all: `10`): Level-symmetric ray set LQ_N (N = 6, 8, 10, 12). INERT here: vet_gd_gl_nmu > 0 is named, and the product set takes precedence (src/rad_m1/rad_m1_vetgd.cpp:790 VgdBuildSet).
- `rad_m1/vet_gd_rotate_xfade` (all: `true`): Cross-fade the closure from the current to the next ray set. Inert with a fixed set (the code ANDs it with "set rotates").
- `rad_m1/vet_rotate_advance` (all: `true`): The averaged window of rotated sets advances. Acts only with vet_rotate_avg >= 2: inert here.
- `rad_m1/vet_rotate_avg` (all: `0`): Average the closure over the last K rotated sets (0 = off). Inert with a fixed set.
- `rad_m1/vet_rotate_avg_ema` (all: `0`): Exponential moving average over the rotated sets (0 = off). Inert with a fixed set.
- `rad_m1/vet_seg_plateau` (all: `false`): Plateau weights of the jump segment. As far as the code shows, the flag is used only by the vet_sc solver (src/rad_m1/rad_m1_vet.cpp), not by the vet_gd sweep used here.
- `rad_m1/vet_source_noesrc_h2be` (all: `true`): Backward-Euler solves read E without the explicit source. INERT here: it acts only together with vet_source_noesrc, which is off with vet_gd (src/rad_m1/rad_m1_time2.cpp:3558 VsrcNoes).

## Keys outside these blocks

The decks also have `<job>`, `<units>`, `<problem>` (the `he_star_m1` problem generator:
initial column, tables, gravity, boundaries, seed) and `<outputN>` blocks. They are not listed
here. The history output and the run-time messages are described in `OUTPUTS.md`.

## Caveats

- The "code default" is the value the code adds when the key is absent, read from the
  `GetOrAdd*` call at the listed line. Several defaults are decided in more than one place
  (for example `<time>/cfl_number`, `closure`, `vet_gd` and `transport` are added by
  `recipe::CleanDefaults` in `src/utils/recipe.cpp` before the modules read them). The
  start-up log prints the defaults it applied in a block of lines such as
  `<rad_m1>/vet_rotate_advance = true (default on a fresh vet_gd run)`.
- "No effect here" is based on the code paths of this commit and these decks (spherical
  polar, fixed rays, one MeshBlock along x1, `scheme = clean`). The same key may act on
  another mesh (for example the cubed sphere) or with rotating rays.
