# AUDIT 2026-09-30 (viper): which run-time keys can be made default?

User 09-30: "we have made many changes. what can be made default?" Read-only audit, no code changed, no jobs run.

Sources: `git log -p --since=2026-09-24 fork/rt-integration -- src` (tip baed412e; 262 src commits),
docs/dev/default_flips_0927.md (three rounds of flips), TASK-2026-09-28-speed-default-audit.md, NOTE-2026-09-29/30-*,
SESSION-2026-09-30-viper.md, memory notes; unmerged branches he-presn-m1 and lhllc-x1-phi;
inputs /viper/ptmp2/jinma/w121prod_0929/w121prod_{1x,10x,3x,3x_ck10}.athinput,
docs/handover/caltech-2026-09-26/inputs/wasp121_*/ (older, pre-T4; superseded by w121prod_0928),
docs/handover/bench-2026-09-29-hebox/inputs/hebox_bench.athinput, inputs/hydro/he_box_m1_1d.athinput,
he-presn-m1:inputs/radiation/he_presn_m1_wedge.athinput.

Classes: **A** make default now; **B** default after one named check; **C** keep opt-in (physics/problem choice,
known trade-off, or user decision); **D** diagnostic, stays off. "Done" = already the code default.

Restart convention (all recent keys): values are written into the restart via GetOrAdd, so a restart keeps what it
recorded; keys read only when named (DoesParameterExist) are NOT recorded, so a flip of such a key reaches old
restarts too.

## 0. Already default (no action; listed so nobody re-flips them)

| key | default | production | evidence |
|---|---|---|---|
| `<mesh>/cs_seam_flux` | positive | positive (w121 1x/10x 300 rot, 3x) | cs-seam-flux-positivity-0926; the multi-rotation check it needed = the 300-rot w121 runs |
| `<mesh>/cs_seam_resample`, `cs_seam_rho_guard` | linear, 4 | same (named) | w121prod_0927 |
| `problem/wall_closed` (dhj) | true (c6ac0848) | not named = true | w121-bottom-wall-energy-0927: Etot_bot = Lrad_bot, Mdot_bot 0 |
| dhj rot_potential + etotgrav energy | fix always on, no key (8fa76784) | on | rotpot-centrifugal-double-count-0928; 1x rot 300 deep converged (ana_rot300) |
| `problem/flux_hst` / `flux_hst_rkavg` | true / true (semantics change of hst columns only) | not named | ae767d20; hydro hst/bin/rst bitwise |
| `<rad_m1>/implicit_opac_newton` (+ `_guard 0.5`, `_guard_mode 2`) | true where opac_update + gas_newton, fresh runs | He box default; He presn names true | NOTE-2026-09-29-m1-perf (-58..-65 %), NOTE-2026-09-29-m1-opn-guard (He box bitwise, wedge 0 clips) |
| `<rad_m1>/implicit_gas_newton`, `closure_lag step`, `closure_thin_relax 1.5`, `vet_tensor full`, `precond mg` levels 3 | on | named = default in He inputs | default_flips_0927 round 1 |
| `<rad_m1>/force_reference_work` | auto (= split where the pgen registers it: box, wedge) | inputs name split (equivalent) | ke_dt_0926, m1-keydefault-0927 |
| `<rad_m1>/implicit_opac_update` | true with `force_reference = wb_arad` | true | ke_dt_0926 |
| `<rad_m1>/implicit_tol`, `time2_lin_tol_fac` | 1e-8, 1 | 1e-8, not named | round 3 (audit_m1_0928: factor 10 = 1.6-4x time error) |
| `<rad_m1>/time_scheme` | hesdirk2 on Cartesian multi-D rk2 fresh runs | hesdirk2 (named) | m1-hesdirk2-cfl-recommendation |
| `<rad_m1>/vet_col_order2`, `fk_min 0`, `reflect_top`, `surface_face` (sp), `vet_col_surface_q` | on (fresh runs) | not named | runs_5q_sporder2b, runs_5e_vetcol2 |
| dhj `ck_implicit`, T4 set, `ck_impl_arat 1e30`, `nosync`, `floorbound`, `kkt_demax`, `kkt_row`, `frozen_cof`, `ck_nquad 2`, `ck_store_split`, `ck_beam_par` | on | on (named) | default_flips_0927 rounds 1-2; round 3 kept store_split / beam_par |
| `<mhd>/cs_resist_x1_centred` | true | (MHD runs) | 50c0a706; uniform grids bitwise |

## 1. Table of open candidates

| key | code default now | production / test inputs | evidence | scope | risk of flipping | class |
|---|---|---|---|---|---|---|
| `problem/ck_impl_maxit` | 8 | 1x 16; 10x 24; 3x 16; 3x_ck10 24 | /viper/ptmp2/jinma/cknewton_0928/RESULTS.md sect. 2-3: a NOT-CONVERGED call keeps an unchecked iterate (up to 12.5 % dT/T at p > 1e-6 bar); 1x 4/500 -> 0/500 at 16 (max 14 passes); 10x max 17 passes (dtmax 0.25) | dhj ck_implicit | ceiling only: bitwise wherever a call converges within 8 passes; extra passes only in non-converged calls; old restarts keep their recorded 8 | **A** (default 24) |
| `problem/ck_impl_conserve` | 0 (read only when named) | 1 in all w121 inputs; C256 must use it (user 09-30) | solution bitwise (ae767d20); fixes the hst gap +7.2e27 erg/s (ana_rot300); GPU cost +1.2 % (ck_conserve_0929/tm, job 12026982: on 48.85/48.91/49.05 s vs off 48.50/47.93/48.55 s) | dhj ck_implicit | hst/ledger columns change, solution never; not recorded in rst, so old restarts also switch (reporting only) | **A** |
| `problem/bbot` fatal if absent (hydro dhj must name `bbot = 0.0`) | GetReal, no default | hydro inputs carry `bbot = 0.0` | pgen.cpp ReadHotJupiterBbot | hot_jupiter pgens | none if the default is 0 (field off); cosmetic | **A** (minor) |
| `<hydro>/lhllc_x1_phi_min`, `<mhd>/lhlld_x1_phi_min` | not on rt-integration (branch lhllc-x1-phi fbe533f3, off = bitwise) | none yet | lhllc-radial-oddeven-0930; w121prod_0929/x1phi_0930/RESULTS.md: V1 (=1) removes the radial 2-cell mode (radial KE 2.0e32 -> 8.8e30), jet -0.6 %; CPU tst 15/15 rc 0 (x1phi_tst.log). USER 09-30: V1 default after the 10x check | every lhllc / lhlld run | changes results of all lhllc/lhlld runs; **He box and He presn also use lhllc** (x1 = vertical/radial there, where low-Mach dissipation was the reason for lhllc: KE -67..-98 % with hlld, mhdvh_0928) | **B**: 10x check (X10L/X10V1) + lhlld smoke (mhd.sub m0-m3) + combined gate, AND a He box KE check (or scope the default to dhj) |
| `<rad_m1>/vet_col_source` | not on rt-integration (he-presn-m1 5a908941: default relaxed; `gas` = old bitwise) | He presn wedge (default relaxed) | root cause of the He presn vet_col checkerboard (commit message: cfl 0.2 grew, gas-only transient source) | every vet_col run (sp wedge, He presn) | changes all vet_col results; he-presn-m1 also raises NHISTORY/NREDUCTION 20 -> 22 (global) | **B**: merge he-presn-m1 with a combined gate + rerun the vet_col accuracy tests (T-S4, Milne) with relaxed |
| `problem/ck_impl_dtmax` | 0.5 | 1x/3x 0.5; 10x and 3x_ck10 0.25 | cknewton_0928: 10x 0/1500 NOT-CONV with 0.25 + maxit 24; 1x 0.25 NOT recommended (3/500) | dhj | per-metallicity; fixed point unchanged (step cap) | **B -> likely C**: Caltech 3x A/B (3649979/80) + viper w3xk (12034178); no production key change without the user |
| `problem/ck_impl_tol` 1e-7, `ck_impl_rsec` 20 (10x only) | 1e-8, 0 | 10x only | 10x Newton history | dhj 10x | looser tol = tolerance-level change | **C** (per input; part of the same A/B) |
| `problem/sponge_bottom` | true (pgen.cpp, all hot_jupiter) | false in all w121prod inputs (user 09-27); caltech-09-26 inputs still true | user 09-27 (supersedes 09-26 "keep both") | hot_jupiter | changes every input that does not name it; recorded in rst | **C**, flag for the user: the code default contradicts the production choice |
| `problem/sponge_top` | true | true | sparc-sponge-campaign-0925 | hot_jupiter | - | done / C |
| c2 levers `ck_impl_once`, `jreuse 0.2`, `pred`, `ck_impl_every 4` | false, 0, false, 1 (`pred_fac` default already 0.5) | named on in all w121 | next-prod-ck-c2 (user production set) | dhj | tolerance-level / cadence changes of results | **C** (user kept them un-flipped, default_flips_0927; no-accuracy-sacrifice) |
| `problem/ck_impl_xstep` | 0 (round 3) | 8 in all w121 | audit_dhj_0928: +3.2/3.9 %, 1.3-2x kick spread | dhj | accuracy | **C** (relaxation-phase only, user 09-28) |
| `<rad_m1>/implicit_one_pass` | 0 | 0 | NOTE-2026-09-29-viper-onepass-cost: -6/-8.5 % at cfl 0.3, 0 accepted at 0.9 | M1 | 5-22x round-off spread | **C** (user 09-29: opt-in) |
| `<mhd>/hlld_bx_zero_tol` | 1e-4 | original | user 09-30 (SESSION-2026-09-30 sect. 4.7) | MHD | - | **C** (user: keep original) |
| `implicit_op_team_red`, `implicit_vimp_fold`, `implicit_halo_ovl_faces` | false | not named | round 3, audit_m1_0928 (< 3 % gain, round-off change) | M1 | - | **C** |
| `implicit_halo_ipc` | false | - | NOTE-2026-09-29-m1-perf: bitwise, no gain (+-2 %) | M1 CUDA | - | **C** |
| multi-rate `implicit_mr_*`, `coupling_mix_*` | off | off | m1-multirate-rejected-0926 | M1 | accuracy | **C** |
| `time2_one_pass_safety` | 30 | bench names 30 | only acts with one_pass | M1 | - | C |
| `problem/ck_impl_osc` | 0 | 0 | cknewton_0928: not needed with dtmax 0.25 + maxit 24 | dhj | - | **C** |
| `problem/ck_interp_T` | linear | linear | c2d3088a (pchip option, not recommended by cknewton_0928) | dhj | - | **C** |
| `problem/ck_dif_dtau`, `ck_dif_margin` | 0 | 0 | ck-fast2 lever, off | dhj | accuracy | **C** |
| `<hydro>/rad_angular` | true (conduction) | false in w121 | irrelevant: route B has `isotropic_conduction = none` | conduction users | - | C (no effect in production) |
| He box `cfl_number` | input | bench 0.3 (timed setting); user 09-29: 0.6 science, 0.9 relaxation | hebox-cfl-saturated-verdict-0929 | He box | - | C (input) |
| pgen/problem keys: `wg_*` (rad_m1_wedge), `he_*`, `mlt_*` (he-presn-m1), `cav_*`, `b0_*`, `m1_b0_*`, `wg_b0`, `remap_file`, `seed_restart`, `atm_seed*`, `ic_profile`, `allow_potential_change`, `gas_teq/vr/vth`, `lat_*`, `albedo`, `bbot_gauss`, `tfloor_kelvin` | problem setup | per input | - | one pgen | - | **C** |
| `flux_hst_floor`, `flux_hst_wall`, `budget_dt`, `ck_impl_ncref/ncloc/stalldbg/verbose/debug`, `dtloc_every`, `dbg_opac_part`, `dbg_hydro_off`, `implicit_timers`, `implicit_dump_op`, `vet_col_dump_every`, `report_newton_fb` (he-presn-m1), `implicit_eos_cache_check_every`, `implicit_krylov_dev` | off | w121 names `flux_hst_floor = true` | flux_hst_floor = two extra active-cell kernels per ConToPrim; ck_impl_verbose costs 1.28 ms/cycle (Caltech 09-27) | - | cost | **D** |

## 2. Rationale per class

**A (make default now).**
- `ck_impl_maxit 24`: a ceiling, not a scheme change; converged calls are bitwise, non-converged calls stop leaving a
  never-checked iterate (the only way the default changes results is by removing a documented error). All four
  production inputs already exceed 8. 16 is enough for 1x/3x at dtmax 0.5 (max 14 passes); 24 covers 10x (max 17).
- `ck_impl_conserve 1`: solution bitwise, fixes the energy ledger the C256 analysis depends on; +1.2 % GPU is the
  price of correct reporting. Every production and the C256 plan use 1 (user 09-30).
- `bbot` default 0 (cosmetic): hydro inputs must currently carry a dummy `bbot = 0.0`.

**B (after one check).**
- `lhllc_x1_phi_min / lhlld_x1_phi_min = 1`: the user already decided (after the 10x check + lhlld smoke). New
  finding here: the key is generic (`<hydro>`/`<mhd>`), so a global default also changes the He box and He presn,
  whose x1 direction is the convective direction where lhllc was chosen for low dissipation. Either scope the default
  to deep_hot_jupiter_rt (pgen sets it when not named) or add a He box KE A/B before the merge.
- `vet_col_source = relaxed`: a real numerical fix, but it lives on he-presn-m1 (40 commits, global NHISTORY change);
  needs the merge gate plus the vet_col accuracy tests with the new source.
- `ck_impl_dtmax 0.25`: wait for the Caltech 3x A/B; the 1x evidence argues against a global default, so it
  probably stays per-input (C).

**C (keep opt-in).** User decisions (one_pass, hlld_bx_zero_tol, xstep relaxation-only, c2 levers not flipped,
sponges), round-3 speed rule (< 3 % gain with a result change), rejected schemes (multi-rate), problem setup keys.
`sponge_bottom` is flagged: the code default (true) is the opposite of the production choice (false since 09-27);
the user may want the default to follow.

**D (diagnostics).** Off by default; `flux_hst_floor` is named on in production and costs two kernels per C2P.
