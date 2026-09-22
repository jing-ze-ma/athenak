---
name: index-convection-boxes-he-bstar
description: Memory links for the He-star and B-star convection boxes, two-stream RT development, and solar convection.
metadata:
  type: reference
---

## He / B-star box index lines

- **[OLD SEED = whole-column v1 kick at box k 1-4 -> directly excites the (2,2) f-mode; ENTROPY SEED (vpert_var=eint, k 8-32, CZ-confined) f523f4d4; productions relaunched 09-16 evening as prod_w6/box_w6](seed-box-modes.md)**
- **[f-MODE GROWTH = STALE WB CACHE (09-16): wb_cache_every 10 -> 1 kills it (gamma 0.86 -> 0.05, 6 sigma); all production inputs carry 10; dt line reproduced on the merged binary](fmode-wb-cache-culprit.md)**
- **[He w8 REACHES 0.9 v_MLT at 13 turnovers (growth +0.15/turnover, saturation close); box_w9 CANCELLED AND DELETED 09-17 09:30 (user); B star steady at 23 turnovers](he-w8-onset-slow.md)**
- **[B STAR prod_w7: linear phase ends at 10 turnovers (stripe mode), DEVELOPED cellular convection by 15 (0.7-0.9 v_MLT); page v7 v2 has dump 15](bstar-w7-linear-saturation.md)**
- **[He w8 NO PULSATION at 7.7 turnovers: only the 148 s cavity mode, sponge-damped in 1 turnover; rt_surface cadence 75 s cannot see 150 s](he-w8-no-pulsation-7turn.md)**
- **[rt_rad_force TRANSVERSE TERM WRONG ON cs/sp (coordinate spacing, not arc length): fixed 4bcdc855 on he4-presn-global; dhj/rg sp runs with the force on are affected](rt-rad-force-transverse-cs-bug.md)**
- **[SPHERICAL TWO-STREAM FIXED 09-17 evening (8bca3dfa: geometry at the faces only; thick test 1e-4, box bitwise); He4 star now survives 1 turnover; next = the upper-FeCZ drain + floor safety](two-stream-sphere-diffusion-error.md)**
- **[f-MODE GROWTH IS O(dt) NUMERICAL (ratio 0.49); CONFIRMED gamma linear in dt through zero; band, force/heating centering, per-stage weights, integrator (rk3) ALL null; Strang/ImEx implementations break the lid (WIP on rt-imex); driver = the split lag itself (numerical kappa-mechanism)](fmode-dt-taper-test.md)**
- **[He BOTTOM SPONGE built 09-16 (e91b12a1; 48 cells/20 s); the 150 s He mode is a box cavity mode (bottom wall x lid cutoff); inflow still needed](he-bottom-sponge.md)**
- **[RT COST PROFILE 09-16: mode-3 column = 78-82 % of kernel time; Strang 1.03x, ImEx 1.42x; speed-ups must come from the column](rt-cost-profile-0916.md)**
- **[rt-integration MERGE DONE 09-15 (bench/wt_merge 6c0df36b): all bitwise; open: ck sweep lost the mode-3 Jacobian (user choice); GPFS INODE QUOTA AT THE HARD LIMIT](rt-integration-branch.md)**

## Solar / stellar convection

- **[He shallow-ramp instability = NUMERICAL: box organ-pipe mode pumped by the two-stream/diffusion flux mismatch at the cut (10 arms, no convergence); sweep FIXED; with the accurate sweep EVERY He ramp grows (100/300 too) -> convergence matrix being redone; He production ON HOLD](rt-shallow-ramp-missplit.md)**
- **[CONDUCTION u-FORM + K ITERATION DONE (rad_x1_uform, rad_x1_kiter 2): worst-face error 13 % -> 0.1 % (B star), 50 % -> 5 % (He); standard for the next production](conduction-frozen-k-error.md)**
- **[He MODE = dt-LINEAR FORCING by the per-cell semi-implicit factor on the non-local exchange; no cheap switch (sub-cycle/re-sweep/clips all refuted); IMPLEMENTED rt_implicit_column (box branch, bitwise off): dt-linearity GONE with the thick/thin split, single-pass solve works (He amplitude 3x down, L swing 20x down, spread 57 %) but B-star floor rises 2.5x -> OFF for production; next = Jiang-style block-tridiagonal implicit column (intensities as unknowns, exact exponential fluxes), scoped; mode 1 proven a dead end (dt pump removed but a larger dt-independent ring remains); mode 3 BUILT: exact, bitwise off, dt-dependence gone, but the dt-independent ring is LARGER (2.2e4) -> whole-column standard scheme: pump dt-independent + tiny, BUT the 0.707 flux loss PERSISTS in every ramp config = a dead band in the per-cell apply block (hunting); mode 3 conservative; B star floor 0.6 % v_MLT](rt-source-dt-forcing.md)**
- **[TWO-STREAM LAYER STAGGER FIXED + PORTED/PUSHED to polar-average-perf d8f00d49 incl. ck + HLLD tol; red giant radial v changes 25 % rms](two-stream-layer-stagger-fix.md)**
- **[He STRANGE MODES: 3 Msun at the critical mass; the shallow-ramp instability is plausibly PHYSICS; taper must not reach the He driving layer](he-star-strange-modes.md)**
- **[THIN-REGION TAPER + FORCE DONE (box branch, default off; IC must go through mk_taper_ic.py)](rt-thin-region-taper.md)**
- **[MULTI-D SHORT CHARACTERISTICS for the thin layers: NEXT PHYSICS PROJECT (09-15) - replaces the interim flux-limited horizontal diffusion above tau ~3; ~1000 lines, plane-by-plane GPU sweep, ALI coupling](multid-short-characteristics-plan.md)**
- **[B-STAR COLUMN PULSATION = NUMERICAL O(dt) anti-damping of the per-cell source (gamma linear in dt, zero at the production dt); ALI restores physical damping -> FINAL CANDIDATE whole column + rt_bottom_flux + ALI; prod_whole RESTAGED (tapered horizontal operator, 8 GPUs 19 h), NOT submitted (user); MODE 3 QUIET at taper 3/10 (floor 7.9e2 < standard 1.4e3, Ftop exact) at 2.2x prod_whole; acceleration + prod_m3 staging running; USER GO 09-15 ~14:30 for MODE-3 PRODUCTION (prod_m3: mode 3 whole column + rt_bottom_flux + accelerated Newton + horizontal RKL1 tapered 3/10; ~1.8x prod_whole -> 2-link chain on 8 GPUs); waiting on the final-binary gates + the m3-relaxed column; PCR (parallel cyclic reduction) started in parallel on branch rg-box-pcr (bench/wt_m3acc), target mode 3 <= 1.2x standard](bstar-column-pulsation.md)**
- **[STANDING: PORT EVERY RT/CONDUCTION FIX TO rg + dhj; MODE 3 -> rg is step S3, DEFERRED by the user 09-15 ("later"), recipe in the file; do it on rt-integration](port-rt-fixes-to-rg-dhj.md)**
- **[B-STAR TOP-FLUX DEFICIT OPEN: two-stream Ftop 0.71-0.88 F_bot while the gas loses exactly F_bot: the 0.707 loss = LAMBDA-ITERATION LIMIT of the per-cell scheme (applies 1/x of the source in thick cells; thick-region T frozen at the IC in every run incl. rg/dhj); ALI diagonal fix b825409c (default on) recovers half; mode 3 = full cure](bstar-top-flux-deficit.md)**
- **[B-STAR OBSERVABLES: SLFV targets 6-130 umag, nu_char ~1.2/d (= our turnover), Schultz+22 method, extend to 40-60 turnovers](bstar-observables.md)**
- **[B-STAR BOX: convection real, e-fold 3.3 turnovers, PRODUCTION RUNNING prod_taper 11701870 (taper+force, fixed sweep, 20 turnovers, ~7 h); prod_uncut cancelled at 6.9 turnovers](bstar-box-growth-and-production.md)**
- **[rad_cap_ang IS A SPURIOUS ENERGY SINK: 20 % of F, all capped box results compromised; run rad_implicit_ang](rad-cap-ang-spurious-cooling.md)**
- **[BOX RESTART BUG FIXED 929d752d: potential not restart state + BC read w0; bitwise now; mom1 1e-7 residual open](box-restart-bug-fixed.md)**
- **[HLLD Bx-ZERO THRESHOLD BUG (beta > 2e4 drops the rotational waves; unstable at M_A <~ 1): fix 1e-8 on lhlld branch; dhj deep field ran without Alfven structure](hlld-bx-zero-threshold-bug.md)**
- **[LHLLD ALIGNED with MM21 (magnetic-floor chi, unit momentum weight), checkerboard gone; BALSARA VORTEX reproduces Leidi (Mach-independent at M_Alf 10); OPEN low-Mach strong-field BLOW-UP](lhlld-solver.md)**
- **[FeCZ COLUMN RELAXATION: rt2turn 3 % E loss = IC/two-stream imbalance out the top; He column 2 % + relaxed; 1-D He column KAPPA-UNSTABLE; formal F_top diagnostic biased](fecz-column-relaxation.md)**
- **[TO DISCUSS: PCG/multigrid for the horizontal implicit operator (user 09-15). IMPLICIT TRANSVERSE OPERATOR WORKING + rad_blend_transverse: horizontal RT transport was OFF above tau 100 in all runs; on everywhere = 3x cost (thin region stiff); PCG relevant again](implicit-transverse-raddiff-plan.md) — PCG deferred (1.25x cost)**
- **[LTE Prad in the EOS above the handover: Prad/Pgas 8.6 at the B-star box top, stiffness + heat capacity wrong, force ok](rt-thin-region-lte-prad.md) — user decision: cut top at tau 1-3 or separate rad energy**
- **[RT HANDOVER BUG FIXED decad388 (wt_rgbox, unpushed, NOT in main branch): semi-implicit step damped the tau-blend handover, 7-13 % of F deposited; SHARED by red giant + dhj](rt-handover-semi-implicit-bug.md) — port to polar-average-perf**
- **[RT PLANE-PARALLEL PORT 5f2be7d2 (wt_rgbox, unpushed): grey two-stream on Cartesian, spherical bitwise, HIP = CPU to 3e-15](rt-plane-parallel-port.md)**
- **[FeCZ BOX PROJECTS: B star 15 Msun + He stars 3/5/8 Msun staged, real X=0 OPLIB table built](fecz-box-projects.md) — 10 turnovers/24 h on 8 GPUs each; RT port in progress**
- **[FeCZ LOW-LUMINOSITY GAP: no multi-D sim of the B-star (8-20 Msun) iron zone](fecz-low-lum-gap.md) — candidate box project; SLFV test**

- **[run/sun viz + WAVES NOT SHOCKS](sun-convection-viz.md) — gravity waves, 87% solenoidal; data on orion**
- [solar_convection binary provenance](solar-convection-binary-provenance.md) — the pinned ideal192/table192 binaries cannot be rebuilt from git
- [solar_convection with the general EOS](solar-convection-general-eos.md) — OPEN: the atmospheric runaway survives the corrected sponge; sub-photospheric comparison is done and solid. RESUME HERE for this thread.
- [solar_convection restart gravity bug](solar-convection-restart-gravity-bug.md) — FIXED/pushed ef9561e2: restarts ran with zero gravity, dt collapsed to 1e-135
- [solar_convection scaling wall](solar-convection-scaling-wall.md) — MEASURED: saturates at ~36 ranks; run it on ONE node, multinode buys 6%
- [Solar convection test](solar-convection-test.md) — tuning solar_convection pgen (two-stream RT) for solar-like convection + stable atmosphere
- [Sponge inert at 192^3](sponge-inert-at-192.md) — RESOLVED: the sponge read pcoord->x1v, a 1x1 placeholder on Cartesian meshes; the 64^2 sponge results are invalid
