# TASK 2026-11-10 viper -> ORION: raycheap15 -- why the cs transport-only deck is 1st order (or worse) (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap15.md`. Goal (user 11-10): implicit M1 truly 2nd order in time
on the CUBED SPHERE. Known: hesdirk2 is p 1.95-1.99 on clean 1-D (raycheap12) and on the 2-D Cartesian pulse with
recon_lag = picard + closure_lag = pass (raycheap13 D3); the cs transport-only deck (raycheap13 Part 3, S0 = raycheap11
T1 recipe) is p 0.16 / 0.46, and S1-S3 (lag keys) and S6 (implicit_lag_restage) are BITWISE = S0.

## viper code reading (fork/raycheap2-1109 dea3e3732)
Why the lag keys are inert on cs, and what is left:
1. **Face coefficients are frozen per stage solve by the thin switch, whatever recon_lag says.**
   rad_m1_implicit.cpp:14504 `doface = aphll && (it == 0 || rst1 || !(rfreeze || impl_thsw != M1_THSW_NONE))`, and the
   lateral faces rad_m1_latface_body.hpp:59 `if (thsw && !first) {return;}`. The cs deck has implicit_thin_switch = fse,
   so S1 (recon picard) cannot reach them. S6 re-evaluates them on pass 1 (rst1) and is still bitwise, i.e. on this deck
   the face coefficients, thin switch, transverse theta and hr-sf block do not depend on the stage iterate at all (they
   are functions of the fixed tensor + opacity). So no "per pass under thin switch" key is needed: it would be bitwise.
2. **The closure is a fixed tensor; closure_lag has nothing to move.** closure = vet_col with vet_gd_replace: the
   chi(f) relaxations are off for fixed-tensor closures (rad_m1_implicit.cpp:1216-1235: implicit_closure_thin_relax is
   forced 0 for vet_col, so S4's `closure_thin_relax 0` was a no-op; S4's change is thin_switch none alone). The tensor
   (D, the vgd half-range data, the hc thin flux DG, the Marshak q) comes from ONE formal solution per step:
   time2_vet_col = predict builds at P = U^n + theta dt K1 (rad_m1_implicit.cpp:13866-13880, rad_m1_time2.cpp:2005-2150),
   theta < 1 where time2_pred_cap (default 0.5) limits (rad_m1_time2.cpp:311-333) -- there the build falls back toward
   the 1st-order lag. Per-pass rebuild exists: `vet_gd_rebuild_every = k` (rad_m1_implicit.cpp:14168, VetGdIterRebuild
   rad_m1_vetlat.cpp:412; rebuilds slot 0 + LAT slots from the Picard iterate E, not the hr data / DG).
3. **With coupling = false the formal-solution source is E itself** (rad_m1_vetcol.cpp:698 `thermal = fl_on && coupling
   && !opac_zero`; src = E of iw EN), so vet_source / relaxed_esc / time2_vet_source do NOT act on this deck. But in a
   hesdirk2 STAGE solve iw EN = U^n + t2inc (the old vector, rad_m1_implicit.cpp:13573), so a lag-mode build sees
   E^n + a21 dt K1, a dt-dependent state. S5's non-convergence (d growing with refinement) says its single frozen build
   differs between levels and/or a per-step term not O(dt^k) accumulates (d ~ number of steps).
4. **Stiffness of the flux drag.** The cs ladder (dt const 0.2956 x cfl/0.15, hydro skipped) has, with ambient chi = 1
   and clump chi = 50 (Ar4: rho 1e4 / 5e5, kappa 1e-4): c dt chi = 0.30 -> 0.037 ambient, **14.8 -> 1.85 in the clump**;
   lateral radiative Courant number at r = 1 c dt/(r dxi) ~ 9.0 -> 1.1, radial 1.9 -> 0.24. The 2-D Cartesian pulse that
   is 2nd order sits at c dt chi = 307 -> 19 (deep diffusion limit). hesdirk2 has stage order 1 (raycheap11): in the
   TRANSITION c dt chi ~ 1 a stage-order-1 SDIRK shows order reduction (Prothero-Robinson), and the cs ladder crosses
   exactly that regime in the clump. Nothing in the Cartesian pulse tests that.
5. Per-pass already (not suspects for a per-step lag): the cs cross terms M1CsTanDer of the iterate (lagged into TRHS,
   recomputed every pass, ImplicitTransverseTerms each pass), the E halo after every pass (ImplicitTransverseHalo(1),
   rad_m1_implicit.cpp:16134), the Krylov operator (implicit_op_stencil). Non-smooth per-step operations that could cap
   p regardless of dt: cs_cross_mono (default true, prod-1108), cs_vertex_mono (true), mesh/cs_vertex_3cell (true, used
   by the M1 off-diagonal VET terms), time2_k1_redist (DEFAULT 0.1 here: one radial MeshBlock, prod-1108), time2_posfix 4,
   implicit_pos_floor(_solve), vet_gd_allow_clamp, the write-back |F| <= cE scale-back (M1ApplyLimits, no key; the inner
   face carries F_in = 1 with E ~ 1, i.e. f ~ 1 there). cs_seam_flux / cs_seam_blockend are hydro (<mesh>) keys: inert
   with dbg_hydro_skip_t2, not run.

## Binary
NO REBUILD: `fork/raycheap3-1110` head `930aad615` = dea3e3732 + tests_raycheap/conv_loc.py only (check
`git diff --stat dea3e3732 930aad615 -- src CMakeLists.txt` is empty; then no keys-absent gate is needed). Use your
/orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732. Take conv_loc.py from that branch (it imports vis/python/bin_convert).

## Base deck
S0 of raycheap13 Part 3 (your Ar4 cs deck = evo_A + gas_rho 1e4, kappa_p/e/f 1e-4, clump_rho 5e5; + coupling = false,
time2_vet_col = predict, dbg_hydro_skip_t2 = true; tlim 20), 24 ranks, cfl 0.15 / 0.075 / 0.0375 / 0.01875. Overrides
below on top of S0 (rad_m1/ prefix unless noted). If a key is refused or a run stops, report and continue.

## Part 0 (no runs; existing raycheap13/14 S0 + S5 logs and bins)
- From each S0 / S5 run.log, quote the end-of-run lines `<rad_m1> hesdirk2: ...` in full (vet clips, k1_redist cells /
  unpaid, time2_pred_cap limited cell-builds of N, fallbacks) and `<rad_m1> m1-positivity: ...` (E floor raises,
  |F| > c E scaled back cell-solves, sum(|F|/(cE) - 1)). Per level, so the scaling with the number of steps shows.
- If the S0 bins still exist: `python3 tests_raycheap/conv_loc.py <c0.15> ... <c0.01875>` (region split: seam band,
  cube-vertex patch, inner/outer 2 shells, the r in [2.6, 3.4] panel-centre patch that holds the clump, rest; per-region
  p; top MeshBlocks; Richardson). Also S5: are its m1_lat / tensor dumps level-independent at t_end (bitwise / max rel)?

## Part 1: is p measuring the time error? (S0, 6 levels)
S0 plus cfl 0.009375 and 0.0046875 (c dt chi clump 0.93 / 0.46, lateral Courant 0.56 / 0.28). conv.py on all 6 and on
the 4 finest; conv_loc.py on all 6 (its Richardson line: p from the 3 finest, error of every level vs E_R).
Expectation if the coarse ladder is just pre-asymptotic: p rising towards 2 at the finest pairs.

## Part 2: one candidate at a time (4 levels each; plus the 2 extra levels for any arm that reaches p >= 1.5)
| arm | overrides on S0 | tests |
|---|---|---|
| A1 | problem/clump_rho=1.0e4 (= gas_rho: no clump) | stiff transition in the clump (c dt chi 15 -> 2) |
| A2 | time2_tableau=trbdf2 | stage order 2 for the same deck |
| A3 | A1 + A2 | |
| B1 | vet_gd_rebuild_every=1 | tensor at the Picard iterate (closure "per pass"; one gd build per pass: if the finest level would exceed ~3 h, run the 3 coarse levels and say so) |
| B2 | time2_vet_col=rebuild | stage-2 build at Y1 |
| B3 | time2_pred_cap=0 | the predictor theta fallback to lag |
| B4 | time2_vet_col=lag vet_col_every=1000000 vet_source_noesrc=true (= S5 + noesrc) | is the frozen build dt-free now? |
| C1 | cs_cross_mono=false cs_vertex_mono=false | cs limiters |
| C2 | mesh/cs_vertex_3cell=false | vertex 3-cell stencil |
| C3 | time2_k1_redist=0 time2_posfix=0 | K1 / stage positivity fixes |
| C4 | implicit_pos_floor=false implicit_pos_floor_solve=false vet_gd_allow_clamp=false | floors / clamps |
| C5 | implicit_flux_x1min=0.3 | inner face f ~ 0.3 instead of ~1 (realizability clip out of play) |
| C6 | implicit_thin_switch=none implicit_recon_lag=picard implicit_closure_lag=pass + A2 | everything per pass + order 2 |

## Part 3: Cartesian 2-D control (your raycheap13 D3 deck, 5 levels, 1 rank)
The pulse is 2nd order in the deep diffusion limit only. Same deck, scattering lowered so the ladder crosses the
transition as the cs clump does: K1 `problem` / rad_m1 kappa_s = 128 (c dt chi 31 -> 1.9), K2 kappa_s = 12.8
(3.1 -> 0.19); each with D3 keys (recon_lag picard, closure_lag pass), and K1T = K1 + time2_tableau=trbdf2. Say where
kappa_s lives in the deck you edit. If K1 drops to p ~ 1 and K1T recovers ~2, the culprit is the hesdirk2 stage order
in the stiff transition, not the cubed sphere.

## Report (NOTE)
Per arm: conv.py d, p (finest pairs), conv_loc region p for S0-6 and every arm with p < 1.5, Picard mean / max /
non-converged per level, the hesdirk2 / m1-positivity end lines, wall per level. Pass criterion: p >= 1.8 at the
finest pairs. CPU only (as before), no code changes; if a key FATALs, report its message and skip the arm.
