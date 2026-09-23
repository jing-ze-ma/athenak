# Runtime-switch audit: the cubed-sphere deep hot Jupiter (prod4 and its hydro twin)

2026-09-23. This is a read-only audit: no code was edited and no jobs were run. The question: is
there any runtime switch that the cs dhj run can reach and that should be on or off, but was
forgotten?

* **Inputs.** `inputs/production/deep_hot_jupiter_cs_prod4.athinput` ("prod4"). Its non-comment
  lines are identical to the running copy `bench/cs_mhd_prod4/deep_hot_jupiter.athinput`. Also
  `inputs/production/deep_hot_jupiter_cs_hyd4.athinput` ("hyd4").
* **Code.** HEAD `d837c355`. prod4's binary was built from `395db5bc`
  (`bench/cs_mhd_prod4/BUILD_COMMIT.txt`). `git diff --stat 395db5bc HEAD -- src` touches only
  `src/rad_m1/*`, `src/CMakeLists.txt` and `src/mesh/build_tree.cpp` (the per-rank restart
  fix). So for everything below, HEAD is the prod4 binary. I also re-checked the key facts in
  `bench/cs_mhd_prod4/src_snapshot`.
* **How the keys were found.** I enumerated every `GetOrAdd*` / `Get*` /
  `DoesParameterExist` call in the following files:
  * `src/pgen/deep_hot_jupiter_rt.cpp` and `src/pgen/pgen.cpp`;
  * `src/utils/two_stream_rt.hpp` (namespace `inline` flags) and
    `two_stream_column_ck.hpp`;
  * `src/diffusion/conduction.cpp` and `resistivity.cpp`;
  * `src/mhd/mhd.cpp` and `src/hydro/hydro.cpp`;
  * `src/eos/eos.cpp` and `eos_table.cpp`;
  * `src/coordinates/coordinates.cpp`, `src/mesh/mesh.cpp` and `build_tree.cpp`;
  * `src/driver/driver.cpp`.

  `src/bvals/` and `src/tasklist/` read no parameters.
* **Evidence sources.** Verdicts reuse the existing evidence and cite it:
  * the memory notes;
  * `docs/dev/switch_inventory_2026-09-22.md` (cited as "inv");
  * the `tests_*` and `bench/*` READMEs.

Verdicts: **OK** / **CONSIDER ON** / **CONSIDER OFF** / **UNTESTED** / **READ FROM THE WRONG
BLOCK**. The last one also covers two related cases: keys the dhj path never reads, and keys it
reads but that have no effect, so a value in the input is silently ignored.

## 1. Main finding: the dhj semi-implicit apply is locked to the legacy linearisation

The following are all verified in the code:

* **`rt_semi_lin` is not read by the dhj pgen.** Three things were checked:
  * `rt_semi_lin` has no `GetOrAdd` in `deep_hot_jupiter_rt.cpp`. Only `box_convection.cpp:1826`
    and `red_giant.cpp:1682/1958` read it.
  * Nothing in the dhj path assigns it.
  * It therefore stays at its namespace default, `true` (`two_stream_rt.hpp:484`). That is the
    old step, `lambda = 4 E/e` and `de = (S/lambda)(1 - e^{-lambda dt})`
    (`two_stream_rt.hpp:5498-5502`). Its own note calls it "kept only to reproduce runs made
    before the fix ... leaves heating explicit and unbounded".
* **`rt_newton` is unreachable, and would do nothing even if set.** It is not read by the dhj
  pgen (namespace default `false`). It also acts only in the `else` branch of `if (semilin)`
  (`two_stream_rt.hpp:5579`), so on dhj it could only take effect together with
  `rt_semi_lin = false`.
* **`rt_relax_sub`, `rt_relax_xcrit` and `rt_relax_submax` have no effect on dhj.** They are read
  (`deep_hot_jupiter_rt.cpp:676-678`), but they live in the same `else` branch. They are the
  "knobs to attack the gap" that `tests_ck_sph/README.md:192` names.
* **The physics mismatch.** `tests_ck_implicit/README_T1_stall.md` measured, in day-side top
  cells, `d ln T/d ln e = 5.4` (H2-dissociation energy carried in e) and `d ln E/d ln T = 9.8`
  (Wien-side bands). Together these give `d ln E/d ln e = 53`. The legacy step assumes 4, so its
  relaxation rate is about 13x too low there. `rt_newton` would take T(e) and c_v from the EOS,
  which gives about 21.6 (it keeps E ~ T^4). That is closer, but still not the band exponent.
* **What it costs today.** `tests_ck_implicit/T0_RESULTS.md` (job 11942457, the prod4 binary,
  prod3 rst 00567):
  * the sweep-to-gas gap `rt_desum` is **-4.7 % to -12.9 %** (mean -8.2 %), and steady;
  * `nclip = 0` and `efix = 0`, so this is the relaxation itself, not the limiter.

## 2. Switch table

Values are "(unset)" when the input does not name the key, so the default applies. "=" means
the same as prod4.

### 2a. `<problem>`: radiative transfer and correlated-k (read by the dhj pgen)

| switch | default (HEAD = prod4 binary) | prod4 | hyd4 | evidence | verdict |
|---|---|---|---|---|---|
| rt_ck | false | true | = | core; the ck path | OK |
| rt_split | true | **false** | = | `rt_ck` forces it true (`deep_hot_jupiter_rt.cpp:708-716`); the line is dead | OK (misleading line: remove) |
| rt_layer_legacy | false | (unset) | = | tests_ck_topslab (c2c layers are correct) | OK |
| ck_spherical | true on cs/sp with rt_ck | true | = | tests_ck_sph (gap 1.8e-2 to 2e-10); tests_gpu_ab_0922 | OK |
| ck_beam_sph | same as ck_spherical | true | = | tests_ck_sph/README_beam.md; REAL_RAYS.md (planet total +0.4 %; twilight only 0.5-0.7 of real rays, fix = short characteristics) | OK |
| ck_sweep_form | 1 (tm) under ck_spherical | 1 | = | tests_ck_sweep_form (tm = sd to 3.5e-17); bench/swform_0922 (tm 21 % faster) | OK |
| ck_sweep_cache | 2 under ck_spherical | (unset) | = | tests_ck_sweep_cost; tm requires 2 | OK |
| ck_nquad | 1 | 1 | = | docs/correlated_k_rt.md:93 (1 vs 2 differ 4e-4 after one cycle; 2 doubles the chains) | OK |
| ck_pcut_bar | 10 | 10 | = | inv; tests_ck_sph | OK |
| ck_int_at_cut | true | (unset) | = | forced false under the tau blend (`two_stream_rt.hpp:2227`), and rad_tau_hi = 300 turns the blend on | OK |
| ck_star_teff, ck_table, ck_data_dir, ck_swflux | 6000 / repo paths | set | = | inv | OK |
| ck_implicit (+ ck_impl_tol, dtol, norm_eps, maxit, dtmax, demax, verbose, debug, refresh_kappa, tau_min, arat, colskip, once, frozen_op, frozen_cof, warm, reuse_jac, seed) | false (all knobs inert when off) | (unset) | = | README_T1.md / README_T1_stall.md: does not converge at arat 2; 4.66x the semi-implicit cost at arat 1e30 | OK (keep off); see the default flip of ck_impl_arat |
| rt_de_max | 0.5 | 0.5 | = | T0_RESULTS B: nclip 0 at the production state | OK |
| rt_use_cons | true | true | = | tests_cleanup_0922/DEFAULTS.md | OK |
| rt_semi_implicit | true | (unset) | = | memory cs-vertex-dt-collapse (explicit overcools) | OK |
| rt_ali_diag | true | (unset) | = | inv (B-star leak fix); no dhj A/B | OK (UNTESTED on dhj, low) |
| rt_floor_consistent | false | (unset) | = | tests_ck_topslab 3a: a measured no-op on the dhj reproducer | OK |
| **rt_relax_sub / rt_relax_xcrit / rt_relax_submax** | 1 / 1.0 / 32 | (unset) | = | **no effect on dhj**: only the `rt_semi_lin = false` branch uses them (section 1) | **READ FROM THE WRONG BLOCK** (read, no effect) |
| rt_cell_report, rt_report_every, rt_outer_verbose, rt_apply_debug(_n), ck_dump_*, diag_gid/diag_cycle_*, photosphere_dump | off | off / empty | = | diagnostics only | OK |
| rt_implicit_column | 0 (anything else is fatal) | (unset) | = | grey-only | OK |

### 2b. two_stream_rt switches the dhj pgen NEVER reads (a value in the input is silently ignored)

| switch | effective value on dhj | box / RG read it as | evidence | verdict |
|---|---|---|---|---|
| **rt_semi_lin** | **true (legacy step)** | box default false; RG `!col3` | section 1; T0 gap -4.7..-12.9 %; T1_stall exponent 53 vs 4 | **CONSIDER OFF** (needs a 2-line pgen read first) |
| **rt_newton** | false | box true; RG col3 | section 1; being measured by another agent (tests_rt_newton/). **Has no effect unless rt_semi_lin = false** | **CONSIDER ON** (together with the line above) |
| rt_rescue_eq | false | box true, RG col3 | only acts on a Newton rescue (T0: efix 0) | OK now; goes with rt_newton |
| rt_src_direct | false | box true, RG col3 | an RG fix for round-off noise in transparent cells (the e = 6e-12 corona). The dhj top-slab e is far larger (not measured) | UNTESTED (low) |
| rt_top_clamp | false | box/RG true | a fix for a poisoned top ghost; the dhj ghost comes from the hydrostatic user BC | OK |
| rt_bface | false | RG | the corona/star join; dhj has no join | OK |
| rt_top_vacuum, rt_top_re, rt_tint_override, rt_bottom_flux | false / 0 | box / RG | a planet top is the unresolved column (the default); the rest are box/RG boundary models | OK |
| rt_grey, rt_plane_parallel, rt_rad_force, rt_force_tau_gate, rt_col3_*, rt_impl_* | off | box / RG | grey / taper / mode-3 only; the ck path does not use them | OK (n/a) |

### 2c. `<problem>`: the planet setup (pgen.cpp + dhj pgen)

| switch | default | prod4 | hyd4 | evidence | verdict |
|---|---|---|---|---|---|
| hot_jupiter, user_srcs | false | true | = | required | OK |
| Teq, omega, grav, ap, Rgas, met, bbot | required | set | bbot 0 (still required by pgen.cpp:57) | tests_hyd4 H5 | OK |
| grav_point_mass | false | true | = | memory grav-point-mass-flag (constant g understates H 2.7x) | OK (default flip candidate, section 4) |
| rot_potential | false | true | = | prod4 header (2); with WB on, rotpot_src is true | OK |
| stellar_tide | false | (unset) | = | fatal on cs (`deep_hot_jupiter_rt.cpp:920-925`); 0.5 % at the limb | OK (n/a on cs) |
| bc_outer_maxwell | true | true | removed (MHD-only) | memory dhj-highB-outer-bc | OK |
| seed / seed_amp | 0 / 0 | (unset) | = | ensemble instrument | OK |

### 2d. `<mhd>` / `<hydro>`: hydrodynamics, well-balancing and the cubed-sphere source

| switch | default | prod4 (`<mhd>`) | hyd4 (`<hydro>`) | evidence | verdict |
|---|---|---|---|---|---|
| reconstruct | plm | plm | = | tests_gpu_ab_0922 (ppmx/wenoz null on dhj, 1.44x cost); master-todo LIST B | OK |
| reconstruct_x1 | plm | (unset) | = | fatal with WB + wb_x1 (`mhd.cpp` ~480); ppm4 costs 10x | OK |
| rsolver | required | hlld | hllc | memory dhj-jet-shallow (cs solvers = sp); tests_hyd4 H2 | OK |
| hlld_bx_zero_tol | 1e-4 (MHD only) | (unset) | n/a | a544a761, 0aaf7f39 | OK |
| fofc (+ fofc_rsolver) | false | (unset) | (unset) | memory cs-dhj-fofc-where (fires only at the night top; prod3_fofc died at rot 11.4) | OK |
| bs_emf | false | (unset) | n/a | tests_cs_loop_growth (true = instability) | OK |
| cs_lowbeta_fallback | 0.5 on cs | 0.5 | n/a | memory cs-mhd-lowbeta-fix | OK |
| polar_emf_diss, polar_hlle_rows, sp_* | sp only | (unset) | = | n/a on cs | OK |
| etotgrav | false | true | = | prod4 header (2) | OK |
| wellbalance_dynamic / wb_x1 / wb_option | false / false / required | true / true / polytropic | = | prod4 header (2); memory wb-hydrostatic-scheme-cs | OK |
| wb_cache_every | 0 | 0 | = | tests_tm_prof_growth (wb-cache arms clean to rot 14.6-14.9); the user keeps 0 | OK |
| wb_x2, wb_rho, wb_rmin/rmax, wellbalance_static(_reconst) | off | (unset) | = | inv | OK |
| **cs_wellbalanced_src** | false; **read from `<mhd>` only** (`coordinates.cpp:43`) | true | **cannot be set** | tests_cs_regions §2: with **plm** on `strat` it is 12-14 % WORSE in the interior and 7-9 % worse at the vertex, neutral on `rot`/`loop`; it pays only with ppmx/wenoz, which are null on dhj. Costs 4 % (memory cs-wb-source-cached). The hydro side is being measured (tests_hyd4/README_wbsrc.md) | prod4: **CONSIDER OFF** (UNTESTED on dhj); hyd4: **READ FROM THE WRONG BLOCK** (fix in progress) |
| scratch_level, ausm_* | 0 / ausm-only | n/a | (unset) | n/a | OK |

### 2e. `<mhd>` / `<hydro>`: radiative conduction (conduction.cpp)

| switch | default | prod4 | hyd4 | evidence | verdict |
|---|---|---|---|---|---|
| isotropic_conduction | absent | radiative | = | inv | OK |
| rad_kappa_src | freedman | table | = | inv | OK |
| rad_met / rad_flux_inner / rad_kappa_fac / rad_flux_limit / rad_cs_exact | 0 / 0 / 1 / true / true | 0 / -1 (the pgen sets sigma T_int^4) / (unset) / (unset) / (unset) | = | inv | OK |
| rad_tau_lo / rad_tau_hi | 0 / 0 | 30 / 300 | = | inv (the blend) | OK |
| rad_kappa_above, rad_gate_rho, rad_gate_dex | 0 / 0 / 0.5 | (unset) | = | stellar gates | OK |
| rad_blend_radial / rad_blend_transverse | true / true | (unset) | = | transverse has no effect with rad_angular false | OK |
| rad_angular | true | false | = | tests_gpu_ab_0922 (noang vs sph at most 7 K; drops the dt2/dt3 terms) | OK |
| rad_implicit_x1 | false | true | = | inv | OK |
| rad_cap_ang | 0 | 0.5 | = | silently does nothing with rad_angular false (`conduction.cpp` comment ~169) | OK (cosmetic: remove) |
| rad_implicit_ang, rad_sts_*, rad_ang_solver, rad_adi_*, rad_tr_* | off / sts | (unset) | = | fatal or inert with rad_angular false | OK |
| **rad_x1_uform / rad_x1_kiter** | false / 1 | (unset) | = | adopted by both box productions; test_rad_x1_uform_cpu (halves the frozen-K flux error); never tried on dhj | **UNTESTED** (CONSIDER ON, low) |

### 2f. `<mhd>`: resistivity

| switch | default | prod4 | hyd4 | evidence | verdict |
|---|---|---|---|---|---|
| ohmic_resistivity | absent | eos | removed | memory resistivity-xe-table | OK |
| max_eta | required | 5e12 | n/a | tests_gpu_ab_0922 "Ohmic cap" (sph 17.4 s CFL-limited); physics caveat in the header | OK |
| use_rkg_sts | false | false | n/a | memory dhj-ideal-xe-floor-relaxation (STS off is 1.54x faster at 1e13) | OK |

### 2g. `<mhd>` / `<hydro>`: EOS and floors (eos.cpp, eos_table.cpp)

| switch | default | prod4 | hyd4 | evidence | verdict |
|---|---|---|---|---|---|
| eos / general_eos | required / gamma | general / table | = | inv | OK |
| eos_h2, eos_ionization, eos_metal_ionization, eos_metal_condensation, eos_radiation | true, true, false, false, false | true, true, true, true, false | = | memory general-eos-project | OK |
| eos_xh, eos_yhe, eos_a_metal, eos_metal_mh, eos_metal_tcond | 0.7381, 0.2485, 16, from problem/met, 0 | set / defaults | = | brief: not re-derived | OK |
| eos_logd/logt min/max, eos_dlog(d/t) | -14/2, 1.5/8, 0.05 | -14/0, 1.5/6, 0.05 | = | memory eos-table-dump (NaN below 71 K is unreachable under tfloor 200 K) | OK |
| dfloor / pfloor / tfloor_kelvin | FLT_MIN / FLT_MIN / none | 5e-14 / 1e-5 barye / 200 K | = | tests_ck_topslab 3b; tests_tm_prof_growth floors (no C2P fail, tfloor plateaus). My arithmetic: at rho = dfloor and mu = 2.3, pfloor 1e-5 barye is T ~ 6 K, so the 200 K tfloor binds first everywhere | OK |
| vceil / vceil_thermalise | 0 / false | (unset) | = | eos_vceil = 0 in every counter table | OK |
| **dfloor_keep_velocity / dfloor_keep_temperature** | false; hydro-only (fatal in `<mhd>`) | n/a | (unset) | used by the RG productions; the hyd4 smoke fires tfloor 20-60x earlier cs runs (tests_hyd4/README.md), not yet isolated | hyd4: **UNTESTED** |
| **efloor_from_ekin / efloor_as_tfloor** | false | (unset) | (unset) | implemented only in the hydro c2p (`general_c2p_hyd.hpp:246,297`). In `<mhd>` they are read and **silently ignored** (`general_c2p_mhd.hpp` has neither) | prod4: **READ FROM THE WRONG BLOCK** (harmless while unset; should fatal like dfloor_keep_*); hyd4: **UNTESTED** |
| **eos_floor_consistent** | false | (unset) | (unset) | **no consumer anywhere in src/** (it only toggles `floors_legacy`, `eos.cpp:118,136`) | **READ FROM THE WRONG BLOCK** (a dead switch in both blocks) |
| sfloor | FLT_MIN | (unset) | = | the table refuses it | OK |

### 2h. `<mesh>`, `<meshblock>`, `<time>`

| switch | default | prod4 | hyd4 | evidence | verdict |
|---|---|---|---|---|---|
| nghost | 2 | 2 | = | enough for plm without FOFC | OK |
| use_cubed_sphere, cs_vertex_fill | false, true | true, (unset) | = | memory cs-wire-fill | OK |
| cs_vertex_fill_cc | false | (unset) | = | 212857c1: 11x better statically, no dynamic effect; tests_mpi_gate_0922 on/off | OK |
| use_grid_stretch_r_poly + f_stretch_r_c1..c4 | false | true + fit | = | docs/dev/dhj_resolution_plan.md; memory radial-stretch-refit | OK |
| meshblock 128x16x16 (24 blocks, 12 per GPU) | - | set | = | tests_rank_imbalance (a 10/14 split is +12.6 % at rot 283 but slower from a cold start, bench/prod4_swap_0922). `lb_nmb_eachrank` is an uncommitted patch, not a switch in this binary | UNTESTED (at a developed state) |
| integrator / cfl_number | rk2 / required | rk2 / 0.3 | = | tests_tm_prof_growth (8): CFL 0.3 vs 0.15 growth ratio 1.06 | OK |
| dt_min | 0 | 1e-3 | = | memory dhj-improvement-list item 6 | OK |
| nan_check_cycles, restart_refill_ghosts | 100, false | (unset) | = | driver.cpp:321-334 (false = bitwise restart) | OK |
| `<units>` mu | - | 1.0 | = | general EOS works at mu_ref 1 (eos.cpp tfloor_kelvin note) | OK |

## 3. Ranked items worth acting on, each with the one test that settles it

1. **Semi-implicit apply form (section 1).**
   * **Change.** Expose `problem/rt_semi_lin` and `problem/rt_newton` in `deep_hot_jupiter_rt.cpp`,
     keeping the current values as defaults. That is two lines and bitwise (gate:
     `test_rad_dhj_ck_cpu`).
   * **Test.** A same-state apudev A/B with the T0 recipe: prod4 binary + that patch, restart
     from `bench/cs_mhd_prod3/rst/dhj.00567.rst`, 2000 cycles, `rt_cell_report` every 20.
   * **Arms.** (a) `rt_semi_lin = true` (the production step); (b) `rt_semi_lin = false`;
     (c) `rt_semi_lin = false` plus `rt_newton` and `rt_rescue_eq`.
   * **Read off.**
     * the `rt_desum` gap (now -8.2 % mean);
     * the efix and rescue counts;
     * top-cell T on the day and night sides;
     * max |de/e|;
     * GPU cycles/s.
   * **Adopt (c)** if the gap drops below 1 % at no more than 1.1x the cost.
   * **For the rt_newton agent:** a dhj arm that sets only `rt_newton = true` tests nothing.
     Also, a value set in the input is silently ignored today.
2. **cs_wellbalanced_src with plm (prod4 = on).** The regions gate says it is neutral to worse
   with plm. The hydro benefit is being measured (tests_hyd4/README_wbsrc.md).
   * **Test.** From the first prod4 restart at 0.5 rot, on apudev, same binary, on vs off for
     about 0.3 rot.
   * **Read off.** dt, radial KE, ME decay, dfloor/efloor rates, cycles/s. If there is no
     measurable benefit, drop it from the next production.
3. **hyd4 floors (hydro-only switches, never tried on dhj).** `dfloor_keep_temperature`,
   `efloor_as_tfloor` and `efloor_from_ekin`, with tfloor firing 20-60x earlier runs.
   * **Test.** Restart from hyd4 smoke rst 00002, 2000 cycles, 4 arms: off /
     dfloor_keep_temperature / efloor_as_tfloor / both.
   * **Read off.** tfloor/dfloor/efloor counts, total-E drift, dt.
4. **Silent-no-op traps (no physics change for prod4; cleanup).**
   * Refuse `<mhd>` `efloor_from_ekin` and `efloor_as_tfloor`, as is already done for
     `dfloor_keep_*`.
   * Delete the dead `eos_floor_consistent`.
   * Make the dhj pgen fatal, or warn, on the `rt_*` keys it does not read.
   * Remove the `rt_split = false` and `rad_cap_ang` lines from the input.
   * **Test.** `athena -n` parse checks, plus the bitwise dhj gate.
5. **rad_x1_uform + rad_x1_kiter 2 below the 10-bar handover (low).**
   * **Test.** A 1000-cycle apudev A/B from prod rst.
   * **Read off.** The conduction flux at the cut against `rad_flux_inner`, dt, cost.

## 4. Default-flip list (all users: dhj cs/sp, box_convection He/B-star, red_giant, cs/sp tests)

Where the defaults live: `two_stream_rt` namespace defaults reach **only the dhj pgen**, because
box_convection and red_giant set every one of them from their own `GetOrAdd`. Flipping them
therefore changes only dhj runs. Those are:
* `test_rad_dhj_ck_cpu`, `test_rad_dhj_ck_mpicpu` and `test_rad_dhj_srclim_cpu` (input
  `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`);
* `inputs/tests/dhj_ck_spherical.athinput` and `dhj_ck_implicit.athinput`;
* the `c_dhj_*` rows of `tests_cleanup_0922/gates.sh`.

**Proposed flips**

| # | switch | now | proposed | evidence | bitwise changes | test still needed |
|---|---|---|---|---|---|---|
| F1 | problem/rt_semi_lin (dhj; must first become readable there) | **true (ON)** | **false** | box and RG both run the equilibrium-relaxation step; inv calls the legacy form "pure back-compat"; T1_stall exponent 53 vs 4 | every dhj input and test listed above | item 1 (the dhj A/B); evidence for dhj itself is missing |
| F2 | problem/rt_newton (dhj) | false | true, only together with F1 | box/RG true (general EOS) | same as F1 | being measured (tests_rt_newton/README.md); **no effect without F1** |
| F3 | problem/rt_rescue_eq (dhj) | false | true, with F2 | box/RG true; only acts on a Newton rescue | same as F1, only where a rescue fires | none beyond item 1 |
| F4 | problem/ck_impl_arat | 2.0 | 1e30 (the two-level split off) | README_T1_stall (GPU 11943285): non-converged 400/400 -> 122/400, gap max 1.7e-4 -> 5.4e-6, cost 5.41x -> 4.66x | only runs with ck_implicit = true: tests_ck_implicit gates, dhj_ck_implicit.athinput | none for the flip; the remaining night corner column is a separate problem |
| F5 | problem/grav_point_mass (dhj) | false | true on sp/cs; must stay false or fatal on Cartesian, like ck_spherical | inv DEFAULT ON; memory grav-point-mass-flag; every dhj production sets it | test_rad_dhj_ck* (their input does not set it) | none (the evidence is physics); the user has not decided |
| F6 | problem/rot_potential (dhj) | false | true | inv DEFAULT ON; every dhj production sets it | test_rad_dhj_ck* input already sets it (no change); other dhj inputs that do not set it would change | none; the user has not decided |
| F7 | mhd/cs_wellbalanced_src | false | **keep false** (listed because another agent is testing it) | tests_cs_regions: neutral or worse with plm; helps only with ppmx/wenoz | - | tests_hyd4/README_wbsrc.md (in progress) |

**Defaults that are ON and should be OFF.** Only F1 (`rt_semi_lin`, a legacy form that is on by
default and reaches only dhj). `rad_angular` (default true) should be off for dhj, but box and RG
need it on, so it belongs in the input, not in the default.

**Not flips, but refusals or deletions (section 3, item 4).** `eos_floor_consistent` has no
consumer. `<mhd>/efloor_from_ekin` and `<mhd>/efloor_as_tfloor` are silently ignored and should
fatal.

**Considered and NOT proposed (the evidence says keep the default).**
* `rt_floor_consistent`: a measured no-op.
* `cs_vertex_fill_cc`: no dynamic effect.
* `wb_cache_every`: 0 is exact, and the user decided to keep it.
* `rt_relax_sub`: no effect on dhj, and the box does not use it.
* `rt_src_direct`, `rt_top_clamp`, `rt_bface`: RG fixes with no dhj evidence.
* `mhd/fofc`, `reconstruct`, `rad_x1_uform`: box productions on, red giant off, so no universal
  value.
