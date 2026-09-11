# Red-giant switches since 40aef0ad, and what they do to other problem generators

Audit of the uncommitted working tree on branch `polar-average-perf` against `HEAD`
(40aef0ad), 2026-09-11.  The rule being checked: every behaviour change must sit behind a
switch that DEFAULTS TO THE OLD BEHAVIOUR, so that `solar_convection` and the
deep-hot-Jupiter problem generators are bit-for-bit unchanged.

## 1. The switches

Conduction parameters live in the `<hydro>` block (the module's block name); the EOS
floor flags likewise.

| block | name | default | what it does | used by |
| --- | --- | --- | --- | --- |
| `<hydro>` | `wb_rmax` | `0.0` (off) | above this x1 the dynamic well-balanced x1 reconstruction and its pressure channel are switched off cell by cell | red giant (3.3e12) |
| `<hydro>` | `dfloor_keep_velocity` | `false` | density floor scales the momentum by `fv = d_old/dfloor` and removes the matching KE instead of leaving `m` and `E` alone | red giant R3/I7 (pending) |
| `<hydro>` | `efloor_from_ekin` | `false` | the energy/pressure floor is paid out of the cell's kinetic energy before any energy is created | red giant (true) |
| `<hydro>` | `eos_floor_consistent` | `false` | tabulated EOS continued linearly in `e` below `T_min` instead of pinning on the bracket | red giant (explicitly false) |
| `<hydro>` | `rad_kappa_rmax` / `rad_kappa_above` | `0.0` / `0.0` | radiatively inert region above a radius (conduction + two-stream together) | earlier red-giant corona runs; I6 sets both to 0 |
| `<hydro>` | `rad_gate_rho` / `rad_gate_dex` | `0.0` / `0.5` | the DENSITY form of the same gate, `kappa_eff = G kappa + (1-G) kappa_above` | red giant I5/I6 (1e-16, 0.5) |
| `<hydro>` | `rad_implicit_x1` | `false` | radial radiative diffusion solved implicitly (tridiagonal, column-local); drops the x1 conduction dt | red giant I2+ (true) |
| `<hydro>` | `rad_cap_ang` | `0.0` (off) | conservative per-face cap on the explicit angular radiative diffusion; drops dt2/dt3 | red giant I2+ (0.5) |
| `<mesh>` | `cs_corner_poison` | `false` | DEBUG: overwrite cube-vertex corner ghosts with 1e30 | debugging only |
| `<problem>` | `nan_report` | `false` | in-kernel capture of the first non-finite flux/cell in HLLC, the WB x1 reconstruction, the hydro fluxes, the conduction flux and the grey RT apply | debugging |
| `<problem>` | `runaway_scan` (+ `runaway_rmin`, `runaway_ratio`, `runaway_ratio_state`) | `false` | ledger of cells far above their neighbours in T, after each operator | debugging |
| `<problem>` | `rt_bface` | `false` | emissivity-weighted far-endpoint Planck source (`BFace`) | red giant only |
| `<problem>` | `rt_use_cons` | `false` | RT reads `e`, `rho` from `u0` instead of the stale `w0` | red giant I6 (true) |
| `<problem>` | `rt_cell_report`, `rt_report_r`, `rt_report_every` | `false`, 3.887e12, 100 | per-cell radiative balance report | red giant |
| `<problem>` | `rt_explicit` | `false` | explicit RT source, bisection control | none |
| `<problem>` | `rt_semi_lin` | `true` = OLD `lambda = 4E/e` linearization | `false` selects the new equilibrium-relaxation step | red giant (false) |
| `<problem>` | `rt_newton` | `false` = OLD | Newton refinement of the backward-Euler balance (general EOS only) | red giant (true) |
| `<problem>` | `rt_rescue_eq` | `false` = OLD 99.9 % drop | non-positive-energy rescue lands on radiative equilibrium instead | red giant (true) |
| `<problem>` | `rt_src_direct` | `false` = OLD flux difference | RT source formed directly (absorption - emission) | red giant (true) |
| `<problem>` | `rt_top_clamp` | `false` = OLD (read the hydro ghost) | the band solver's top slot `i = ie+1` takes its state from the top ACTIVE cell `ie` | red giant (true) |

Red-giant-pgen-only parameters (no other pgen compiles `red_giant.cpp`, so they are safe
by construction): `opac_tmin`, `opac_floor`, `bg_hydrostatic/bg_rho/bg_temp/bg_rtop`,
`bc_use_cons`, `mlt_*`, `vpert_*`.

Two new untracked headers, `src/utils/eint_from_cons.hpp` and `src/utils/runaway_scan.hpp`
(`runaway_scan::on` defaults false).  Neither is in `src/CMakeLists.txt` (headers only).

## 2. Inventory: changes NOT behind an old-behaviour switch

| file | hunk | change | category | pgens affected |
| --- | --- | --- | --- | --- |
| `two_stream_rt.hpp` | ~1041, 1076, 1097, 1162 | the band solver's top slot `i = ie+1` now reads cell `ie` (`ii = (i>ie)?ie:i`) instead of the hydro ghost | unconditional behaviour change | `deep_hot_jupiter_rt` (and red giant) |
| `two_stream_rt.hpp` | 1839-1866 | `rt_src_direct` default `true`: source formed directly, not as `-(Ft-Fb)/dx` | unconditional (switchable, wrong-way default) | `deep_hot_jupiter_rt` |
| `two_stream_rt.hpp` | 1877-2110 | the semi-implicit step rewritten (relax to `e_eq = e (A/Em)^(1/4)`); old form only under `rt_semi_lin = true` | unconditional (wrong-way default) | `deep_hot_jupiter_rt` |
| `two_stream_rt.hpp` | 1877-2110 | `rt_newton`/`rt_rescue_eq` default `true` (general EOS only) | unconditional (wrong-way default) | any general-EOS two-stream run |
| `two_stream_rt.hpp` | 916-934, 735-780 | startup prints (`rt_use_cons` state, `rt_top_re` no-op warning) | diagnostic | all two-stream pgens |
| `hllc_hyd.hpp` | 83-93 | general-EOS branch falls back to `p_floor/(G1-1)` when the reconstructed `w(IEN)` is non-finite or `<= 0` | bug fix (no-op on healthy data; a NEGATIVE reconstructed eint now behaves differently) | all general-EOS HLLC runs |
| `hydro_fluxes.cpp` | 275-300 | non-finite reconstructed `dl/dr(IDPR)`/`wl/wr(IEN)` repaired consistently before the `fmax(.,pfloor)` mask | bug fix (no-op on finite data) | all hydro runs |
| `wb_background.hpp` | 91-110 | `WBGuard` also flattens states more than 1e6 away from the anchor | bug fix (no-op on healthy data) | all well-balanced pgens (dhj included) |
| `wb_background.hpp` | 327-350 | polytropic walk limits the temperature DECREASE over a segment to `0.5 T0` | bug fix (identical when the clamp does not bite) | all `wb_option = polytropic` runs |
| `hydro_tasks.cpp` | 191 | new `ImplicitConduction` task inserted between `srctrms` and `sendu_oa` | unconditional, numerically a no-op (early return unless `rad_implicit_x1`) | all hydro pgens |
| `hydro_newdt.cpp` | 104-200 | extra `MinLoc` reduction + a one-cell diagnostic kernel when dt collapses | diagnostic | all hydro pgens |
| `coordinates.cpp` | 819-941 | `cs_raisev` `par_for` -> `parallel_reduce` (floor-energy accounting); arithmetic unchanged at default flags, `w0` velocity writes moved after the floor block | diagnostic + refactor | cubed-sphere pgens |
| `general_hyd.cpp`, `general_c2p_hyd.hpp`, `coordinates.cpp` | several | `efloor_de` accumulation (a 4th reduction slot) | diagnostic | all general-EOS pgens |
| `eventlog.cpp`, `mesh.hpp` | all hunks | `efloor_de` column in the event log + `MPI_Allreduce` for it | diagnostic | all pgens using `file_type = log` |
| `driver.cpp` | 469-620 | NaN check separates active from ghost cells and prints the state; extra reduction runs ONLY once something is already non-finite | diagnostic | all pgens |
| `mesh.cpp` | 983-1030 | hydro dt-collapse report; conduction dt report says "implicit"/"capped" | diagnostic | all pgens |
| `bvals_cc.cpp`, `prolong_prims.cpp` | 1 hunk each | `cs_corner_poison` hook; C2P signature threading | switch-gated / signature only | none at default |

`src/utils/atm_column.hpp` and `src/utils/correlated_k.hpp` are UNCHANGED since HEAD.

## 3. `rt_bface` audit

`BFace(k_own, k_far, b_own, b_far, on)` returns `b_far` immediately when `on` is false, so
the switch is exact.  All seven call sites pass `bface_on` (grey down/up sweeps and the
grey `Em`; the two ck sweeps and the ck `Em`; the generic-chain down/up sweeps and its
`Em`), and each one's `b_far` argument is the expression the pre-fix code used, so
`rt_bface = false` is the old code bit for bit.  The two 4x "guard" tests in the ck path
only decide whether the extra neighbour-opacity lookup is paid; the result is discarded
when the switch is off.  The parameter is read with `GetOrAddBoolean("problem",
"rt_bface", false)` at all four places the pgen configures the solver.

**FINDING: no red-giant input turned it on.**  Neither the live
`/orion/ptmp/jinma/Athenak/red_giant/I6_fixes/rg.athinput` nor the repo input
`inputs/hydro/red_giant_cs.athinput` set `rt_bface`, so every run since the fix was gated
has been running WITHOUT it.  `rt_bface = true` has now been added to
`inputs/hydro/red_giant_cs.athinput`; the live run's input was deliberately not touched.

**VERDICT, 2026-09-11 (implemented).**  All five number-moving items are now switched
with the OLD behaviour as the default, in the shared header and at both red-giant
configure sites: `rt_src_direct = false`, `rt_newton = false`, `rt_rescue_eq = false`,
`rt_semi_lin = true`, and the band solver's top-slot clamp behind the NEW
`problem/rt_top_clamp` (default `false` = read the hydro ghost, exactly as HEAD).
`inputs/hydro/red_giant_cs.athinput` sets all five to the new behaviour next to
`rt_bface = true`, so the red giant is unchanged and every other pgen reverts to HEAD.
Verified with the harness in `/orion/ptmp/jinma/Athenak/regress/`:

- `deep_hot_jupiter_rt`, correlated-k band path, default input, 50 cycles, 16 x 7:
  `run_ck_new` (build `b_new_dhj`) vs the HEAD run `run_ck_hd` -- **BIT-IDENTICAL**,
  0 differing cells of 524288 in dens, velx, vely, velz, eint, bcc1, bcc2, bcc3.
- `red_giant`, restart from `V9f_rtfix/rst/rg.00002.rst`, 100 cycles, 32 x 7, the
  `I6_fixes` input plus `rt_bface = true`: `run_rg_new` (build `b_new_rg`, with the five
  opt-in lines) vs `run_rg_pin9` (`build_impl/athena_pin9`, whose defaults ARE the new
  behaviour, without them) -- **BIT-IDENTICAL** at t = 503073, cycle 16421, 0 differing
  cells of 2949120 in dens, velx, vely, velz, eint.

## 4. Regression vs HEAD

Worktree `HEAD` = 40aef0ad at `/orion/ptmp/jinma/Athenak/regress/wt_head` (removed after
the test).  Four builds, identical CMake options (MPI + OpenMP, gcc 13 / OpenMPI 4.1, no
Kokkos arch flags, Release), 16 ranks x 7 threads, 100 cycles, final `hydro_w` dump
compared cell by cell.

| case | input | result |
| --- | --- | --- |
| `solar_convection` (general EOS table, 96x64x64, its own two-stream RT in the pgen) | `inputs/hydro/solar_convection.athinput`, `nlim = 100` | **BIT-IDENTICAL** — max abs diff 0, max rel diff 0 in dens, velx, vely, velz, eint (393216 cells each) |
| `deep_hot_jupiter_rt`, PICKET-FENCE path (ideal EOS, spherical polar) | `run/dhj_sph_rt_test/deep_hot_jupiter.athinput`, `nlim = 100`, `bbot = 0` added | **BIT-IDENTICAL** — 0 differing cells in dens, velx, vely, velz, eint (458752 each) |
| `deep_hot_jupiter_rt`, CORRELATED-K band path (general EOS table + MHD + polytropic WB) | `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`, `nlim = 50`, table paths absolutised | **DIFFERS**: dens max rel 1.5e-2, eint 2.0e-2, velx max abs 4.9e3 (310k-440k of 524288 cells per variable) |

Attribution of the band-path difference, by rebuilding the working tree with individual
changes reverted (`/orion/ptmp/jinma/Athenak/regress/wt_mod`):

| build | what was reverted | vs HEAD |
| --- | --- | --- |
| `b_mod_dhj` | `rt_semi_lin`, `rt_newton`, `rt_rescue_eq`, `rt_src_direct` back to the old behaviour | still DIFFERS, but ~7x smaller (dens 1.9e-3, eint 1.2e-3) |
| `b_mod2_dhj` | the above + both `wb_background.hpp` clamps | byte-for-byte the SAME as `b_mod_dhj` -> the WB clamps never bit in this run |
| `b_mod3_dhj` | the above + the band solver's top-slot clamp (`ii = (i>ie)?ie:i`, and the top-ghost `B` read) | **BIT-IDENTICAL** (0 differing cells in all 8 variables) |

`b_mod3_dhj` reproduces HEAD exactly, so the band-path difference is ENTIRELY those five
items: the four RT switch defaults (most of it) and the top-slot change (the remainder).
Nothing else in the audit moves a number -- the NaN repairs and both WB guards are
confirmed inert on healthy data.

The picket-fence path is untouched because `rt_src_direct` only acts when `rt_ck` or
`rt_grey` is on, and the rewritten equilibrium step lives in the split-band `rt_apply`.

## 5. Verdict: what is safe for other problem generators

- Everything in the conduction module, the EOS floors, the well-balanced cut-off and the
  cubed-sphere corner debug is properly gated and defaults to the old behaviour.
- `solar_convection` is unaffected: it does not include `two_stream_rt.hpp` (it carries its
  own `two_stream_RT`), and the regression above is exact.
- The one real exposure is `src/pgen/deep_hot_jupiter_rt.cpp`, which shares
  `two_stream_rt.hpp`.  It inherits four changes it never asked for: the top-slot clamp
  (unconditional), and the three switches whose defaults are the NEW behaviour
  (`rt_src_direct`, the equilibrium step under `rt_semi_lin = false`, and
  `rt_newton`/`rt_rescue_eq`).  RESOLVED 2026-09-11: all four defaults were flipped to the
  old behaviour, the top-slot clamp was put behind `problem/rt_top_clamp` (default false),
  and the red-giant input opts in to all five -- see the verdict at the end of §3.  Measured: on
  `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` they move the state by 1.5 % in density
  and 2 % in internal energy after only 50 cycles.
- The bug fixes (HLLC/reconstruction NaN repair, `WBGuard`, the polytropic walk clamp) are
  no-ops on healthy data by construction and are confirmed no-ops by the regression.
