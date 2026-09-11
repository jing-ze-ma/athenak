## RULE ZERO (the user has asked THREE times, 2026-09-09 and 09-10)

- **[SAVE TOKENS](save-tokens.md) — user rule 2026-09-11: never spend tokens on unneeded work; narrow agent prompts, haiku/sonnet for lookups, no agents before the direction is chosen**

- **[DELEGATE simple tasks to OPUS 5 subagents](delegate-simple-tasks-to-opus.md) — FROM THE FIRST TOOL CALL of every session: reading handover notes, git/merge/build/test cycles, creating run dirs, submitting and polling jobs, log reads, analysis scripts. Fable thinks, designs, decides; Opus (`Agent`, `model: "opus"`) executes. Do not wait to be reminded**

## Cubed sphere — start here

- **[RT SEMI-IMPLICIT source](rt-semi-implicit-changes-dhj-answer.md) — 3 days old (048dff30, 09-08); the 09-09 switch to explicit was for ABLATION REPRODUCTION only; the explicit/semi-implicit difference is the EXPLICIT scheme's overcooling error (see cs-vertex note). READ before any A/B claim on cs**
- **[IN FLIGHT, STOPPED 2026-09-11 ~03:30](inflight-2026-09-09-viper.md) — sp_mhd_prod3 CLEAN at rot 23.8 on apu (check it first); cs collapse CLOSED; instrument + cap + handover bundle PUSHED (c26b01ac); explicit ens arms CANCELLED. Next = dt_min 1e-3, density floor vs WB, longer si ensemble. START HERE**
- **[GPU BUILD RACE: flag never compiled](gpu-build-race-flag-not-compiled.md) — verify flags with `strings` + a startup print before trusting any arm**
- **[cs RADIAL = sp RADIAL, DONE 979edada](cs-radial-unification.md) — centroid everywhere, x1 Grid-PLM always; every cs baseline before it is stale**
- **[cs STRETCHED-grid SOURCE TERM off 0.57-2.15x: FIXED 13a97399](cs-stretched-source-term-bug.md) — index-space dr in the curvature source; affects every stretched cs dhj run**
- **[cs STRETCHED-grid RESISTIVITY was ANTI-DIFFUSIVE: FIXED c5c85e3b](cs-stretched-resistive-rcm-bug.md) — r_cm missed the stretch; pre-fix cs resistive MHD suspect. START HERE for the deep sheet**
- **[cs seam CO-LOCATION: FOUND and FIXED](cs-seam-colocation-fixed.md) — the clamped-window cubic; 5-11x on the evolved field. START HERE for the seam halo**
- **[cs seam order: 2nd order CONFIRMED](cs-seam-order-limiter.md) — global L1 2.01; never read order off region bins. The seam is DONE**
- **[cs seam: de-staggering fix REFUTED](cs-seam-destag-refuted.md) — superseded by [[cs-seam-colocation-fixed]]**
- **[resistive seam 1st order: MECHANISM CLOSED](cs-resistive-seam-order.md) — operator innocent (1.98), the halo INPUTS were guilty**
- **[cs cross-level SEAM halo: CLOSED](cs-crosslevel-seam-halo-first-order.md) — a01ace75 + 180a9b3e. START HERE for cs+SMR halos**
- **[cs CUBE-VERTEX corner halo x RADIAL GHOST: FIXED](cs-cube-vertex-corner-radial-ghost.md) — 397b4ad3. START HERE for cs halo bugs**
- **[cs seam conservation CLOSED](cubed-sphere-seam-conservation.md) — 985faa22; exact to round-off, and history.cpp used the CARTESIAN volume on every grid**
- **[cs seam EMF CLOSED](cubed-sphere-seam-emf.md) — 42323a66; exact at 1/4/16 blocks per panel. div B is the WRONG gate. START HERE for cs+MHD**
- **[cs RESISTIVITY](cubed-sphere-resistivity.md) — 4bfacdd8; gnomonic two-pass curl, 2nd order; refinement since c4735f18. START HERE**
- **[cs SMR: refined MHD CONVERGES, MPI-clean](cubed-sphere-smr.md) — 74cbc8df; edge FLUX buffers lacked the seam transform. "Rank dependence" RETRACTED. START HERE for cs+SMR**
- **[cs MHD convergence](cubed-sphere-mhd-convergence.md) — e63b571a; mhd_corner_e had no cs form. Residual is the radial-BC phase lag**
- **[SHOCKS through cs seams](cs-shocks-through-seams.md) — structure ON a seam NaN'd SILENTLY; the resample is now clamped. FOFC DECIDED not needed. START HERE**
- **[cs blast vs a CARTESIAN grid](cs-blast-vs-cartesian.md) — the seam costs NOTHING (4.1% vs Cartesian 4.5%); only the vertex is ~2x worse. START HERE for shocks**

## Cubed sphere — validation and limits

- **[cs_test FF decay + rot_axis (4990eb41)](cs-test-ffdecay-rotaxis.md) — cs matrix vs sp, no 1.5-order region on cs. START HERE for cs-vs-sp simple tests**
- **[cs PURE HYDRO validation](cs-hydro-validation.md) — mass and energy to MACHINE PRECISION under a shock; L1 2.3-2.4. Records the FALSE first-order reading. START HERE**
- **[cs MHD + RESISTIVE validation by region](cs-mhd-validation.md) — RE-MEASURED: ideal MHD 2nd order everywhere; the resistive energy drift is PHYSICAL Ohmic heating, NO BUG**
- **[cs ANGULAR MOMENTUM: sp exact, cs cannot be](cs-angular-momentum.md) — rigid rotation 6.9e-4/rot at nx2=32; in dhj runs Lz tracks sp to 1e-4 over 100 rot**
- **[cs: radiation + srcterms now REFUSED](cs-unsupported-physics-guards.md) — b283ad3a; they used to run and return a wrong answer**
- **[cs NARROW-BLOCK resample degeneracy](cs-narrow-block-resample-degeneracy.md) — the along-seam stencil INVERTS below 3 cells (MeshBlock under 6). Test at the minimum block size**
- **[cs GENERAL/TABULATED EOS stale cache: FIXED](cs-general-eos-stale-cache.md) — p/Gamma_1/T cached BEFORE the gnomonic correction, 6-15% wrong**
- **[cs MHD blast from a VECTOR POTENTIAL](cs-mhd-blast.md) — b283ad3a; div B zero by construction. Gate with TWO numbers — div B alone passes a ZERO field**
- **[cube-vertex corner fill: PROMOTED, default ON](cs-wire-fill-wip.md) — `<mesh>/cs_vertex_fill`; ownership is the exact chart DIAGONAL. Merged**
- **[cs cube-vertex REAL fill: prototyped, NOT built](cs-cube-vertex-real-exchange.md) — sampling beats extrapolation 3.4-5.5x, needs no new exchange**

## Cubed sphere — history

- **[Cubed sphere ON GPU](gpu-this-capture-device-lambda.md) — af539941; hydro and MHD match CPU. The constant-memory CLOSURE LIMIT**
- **[cs GPU multi-block fault: FIXED](cs-gpu-multiblock-fault.md) — 4302a008; one kernel captured both send and recv buffers**
- **[Cubed sphere + MPI: FIXED](cubed-sphere-mpi-hang.md) — 1c2e29d6; unmatched receive at a cube vertex. Run the RECIPROCITY AUDIT first**
- **[cs MHD seam FIXED](cubed-sphere-mhd-seam.md) — 4f19a244; three plumbing bugs, the transform was always exact**
- **[Cubed sphere: MHD](cubed-sphere-mhd.md) — e0c74357; dxedge was zero, mhd_fluxes had no gnomonic rotation**
- **[Cubed sphere: COMMITTED as 9492a946](cubed-sphere-committed.md) — five sessions of hydro work is in git**
- **[cs along-seam resample, 2nd ORDER](cubed-sphere-seam-interp.md) — the "flat interior residual" is RETRACTED**
- **[cs seam basis transform, FIXED](cubed-sphere-seam-basis.md) — a tangent-BASIS transform, not a signed permutation**
- [Cubed sphere panel frames](cubed-sphere-panel-frames.md) — forced by panel_neighbors; a duplicate copy had 3/4 swapped
- [Cubed-sphere hydro fix](cubed-sphere-hydro-fix.md) — halo axis bug FIXED; numbers predate the panel 3/4 swap
- [Cubed sphere: x1 as radial](cubed-sphere-x1-radial.md) — DONE, bit-identical; axis roles cyclically permuted
- [Cubed sphere: hydro state](cubed-sphere-hydro-state.md) — 4 early bugs FIXED; its rigid-rotation section is RETRACTED
- [Cubed sphere for the hot Jupiter](cubed-sphere-for-hot-jupiter.md) — the dt argument is DEAD (radial binds)

## Hot Jupiter: physics, runs, campaigns

- **[dhj JET is SHALLOW, deep equator WESTWARD to rot 164 on sp](dhj-jet-shallow-westward-deep.md) — all cs solvers = sp; the deep super-rotation never spins up. hllc upper jet HALF of ausm+up/sp. START HERE for any jet question**

- **[cs VERTEX dt COLLAPSE: CLOSED 09-11](cs-vertex-dt-collapse-0907-defaults.md) — explicit RT source overcools stiff cells by x/(1-e^-x) (verified 0.04%), a NIGHTSIDE cold patch sinks, the vertex column drains, hydro reheats it, RT cascade, conduction reports; semi-implicit default is RIGHT and self-heals; rad_tmax_kappa is no rescue; si/s05 = marginal dt_min case. START HERE**
- **[RESTART dropped rot_potential: centrifugal force DOUBLED after every chain restart, FIXED](restart-rot-potential-bug.md) — 0.4 % of g; all first-round restart ablation arms INVALID; hllc/ausmpup/prod2 chains ran on it. READ before trusting any restarted 09-07-default run**
- **[WB RESTART cache bug, FIXED 5c0b98e4 (orion)](wb-restart-cache-bug.md) — zero background for 9 cycles after every restart; not the collapse trigger**
- **[sp POLAR BLOW-UP BISECTED to 863e8337 and FIXED](sp-pole-bisect-culprit-863e8337.md) — `mesh/polar_x3_shift=rotate` default (9a9396f7); A/B shift 6.4e34 vs rotate 7.2e31; sp_mhd_prod3 clean past rot 23 (old prod died 10.3). START HERE for sp MHD**
- **[cs_mhd_prod DIED at rot 41.6](cs-mhd-prod-nan-rot41.md) — NaN everywhere in one interval, no precursor; the code does not stop on NaN. START HERE for the old cs production**
- **[cs_hyd_rs: hydro cs twin x {hllc, lhllc, ausmpup, ppmx}](cs-hyd-rs-run.md) — the jet comparison. START HERE for the solver comparison**
- **[cs_mhd_prod2: cs MHD production FROM SCRATCH on the 09-07 defaults](cs-mhd-prod2-run.md) — WB+rot_potential+deep RT; A/B vs cs_mhd_prod. START HERE for current production**

- **[Deep RT: radiative conduction + tau_R BLEND fad2f5db](radiative-conduction-deep-interior.md) — `rad_tau_lo/hi`; `rad_kappa_src=table` makes the handover consistent to 10%. START HERE for deep RT**
- **[cs MHD dhj VALIDATED at rot 20](cs-dhj-diagnostics-rot20.md) — T structure = sp, jet WEAK (~1 km/s) on every grid, B does not slow it. START HERE for any cs dhj science plot**
- **[sp_mhd_prod DIES at rot 10.3: polar-row RADIAL field at the BOTTOM](sp-pole-bottom-radial-blowup.md) — BISECTED: 2a64e2c7 clean, 18f5dd21 dirty, resistivity innocent. START HERE for sp MHD**
- **[cs MHD dhj BOTTOM-BOUNDARY DRIFT: FIXED 248f1b77](cs-mhd-bottom-inflow.md) — cs ghost bcc as plain face averages, 8% energy drift by rot 16, NaN at rot 32; redo pre-fix runs. START HERE**
- **[sp MHD ENERGY EXCESS vs hydro, ABSENT on cs](sp-mhd-energy-excess.md) — +1%/2 rot; floors on the upper nightside, cause OPEN. READ before comparing sp MHD with anything**
- **[Correlated-k design + build](correlated-k-design.md) — THE MAIN THREAD. COMPLETE: physics, docs, tests, MPI. START HERE**
- **[dhj RT on the CUBED SPHERE: WORKING](dhj-cubed-sphere-port.md) — tracks spherical polar; four bugs found, three PRE-EXISTING**
- **[cs NON-ORTHOGONAL audit: BOTH FIXED](cs-nonorthogonal-audit.md) — e415b91a; history KE/ME ~1% wrong on cs, newdt understated the fast speed**
- **[cs MHD bug FIXED: the C2P floor corrupts u.e](cs-mhd-c2p-floor-corrupts-ue.md) — magnetic energy from the NON-ORTHOGONAL triple; vertex order -0.04 -> +1.97, but the dhj blow-up stays OPEN**
- **[cs MHD low beta: per-step operator 2nd ORDER incl. the vertex](cs-mhd-low-beta-divergent.md) — the instability is NONLINEAR; HLLE fallback is legitimate. START HERE**
- **[cs MHD low-beta: fixed ON THE TEST PROBLEMS ONLY](cs-mhd-lowbeta-fix.md) — cb1afd28, `<mhd>/cs_lowbeta_llf`, default ON for cs. INSUFFICIENT: dhj still dies at 0.21 rot**
- **[cs SEAM GHOST e_int: RETRACTED as a mechanism](cs-seam-ghost-eint.md) — the seam degrades EVERY quantity ~15-19x; E-vs-B inconsistency is ~14%. READ THE RETRACTION FIRST**
- **[GENERAL-EOS AUDIT](general-eos-audit.md) — the (rho,e)->T inversion is bracketed and monotone, no live bug; off-table states silently extrapolated**
- **[3D RESISTIVITY BUG: curl B missing theta terms on sp + Cartesian, FIXED](resistivity-3d-curl-missing-terms.md) — every sp resistive dhj run since June is wrong; cs unaffected. START HERE**
- **[sp POLE-EDGE CURRENT: far-side segment WRONG SIGN, FIXED](sp-pole-edge-current-sign.md) — J_r grew ~1/dtheta; the m=2 polar pressure pattern**
- **[sp POLE: resistive J_r divided by a ZERO dual area, FIXED](sp-pole-edge-area-zero.md) — NaN in cycle 1 on any uniform-theta grid**
- **[sp POLE FIXED: polar_emf_diss + x3 face-state theta shift](sp-pole-fixes.md) — de667f32/863e8337; polar-row L1(B) rate 0.76 -> 2.1. START HERE for the pole**
- **[sp_test: RIGID ROTATION + RESISTIVE DECAY rates](sp-test-rigidrot-resist.md) — hydro 2.6, resistive B 2.1; POLAR ROWS do not converge for a tangential field**
- **[RESISTIVITY AUDIT](resistivity-audit.md) — eta_b ghosts, resistive dt, pole-face curl all CLEAN; one gap fixed behind `use_polar_average_eresist`**
- **[sp POLAR FIELD BLOW-UP: FIXED 1cabe85c, default ON 3147de4f](sp-polar-field-blowup.md) — GS05 corner-EMF grows a face checkerboard. READ BEFORE ANY cs-vs-sp MHD CLAIM**
- **[cs ROTATION SOURCE BUG: cs ran the CARTESIAN beta-plane](cs-rotation-source-bug.md) — f75ad783, removed ae831665. READ BEFORE TRUSTING ANY cs-vs-sp COMPARISON**
- **[cs 1-ULP amplification: old nx=64 claim STALE](cs-ulp-amplification.md) — round-off grew to O(1) and killed nx=64 at 0.309 rot on pre-fix binaries; retest ALIVE. START HERE**
- **[WHAT RECONSTRUCTION ACTUALLY RUNS](reconstruction-what-actually-runs.md) — sp ignores `reconstruct` entirely; cs honours it but x1 ignores the radial stretch. READ before any reconstruction claim**
- **[cs VERTEX = LIMITER CLIPPING](cs-vertex-limiter-clipping.md) — ppmx and WENO-Z cut the worst-cell vertex force 6.5x; ppmx is the production candidate (nghost 4)**
- **[cs VERTEX ORACLE: halo innocent](cs-vertex-oracle-halo-innocent.md) — exact ghosts change the residual <1%; it is the CELL BALANCE. Seam continuation DECIDED against**
- **[sp POLAR-ROW RESIDUAL FIXED 7.7x: `mesh/polar_quadratic_recon` 6410be99](sp-polar-row-reconstruction.md) — default OFF; limiters see the pole as an extremum. READ before any sp accuracy claim**
- **[sp GEOMETRIC-SOURCE RESIDUAL: polar rows 14x interior](sp-geometric-source-residual.md) — WB source does NOT touch it; Cartesian update 34882bd7 fixes only the transverse field. "WB will fix it" RETRACTED**
- **[cs WELL-BALANCED SOURCE CACHED: 46% -> 4% of a run](cs-wb-source-cached.md) — 88064f67; bitwise identical, restart-clean; cs_test cannot restart**
- **[GS07 corner EMF UNCONDITIONAL on cs since b71ef392](cs-gs07-emf-gate.md) — same order as the plain average, 1.8x its L1(B); `bs_emf` is the only way back**
- **[cs DEEP TOROIDAL SHEET: CLOSED](cs-deep-toroidal-sheet.md) — it was the stretched resistive bug c5c85e3b; rcmfix has NO sheet. START HERE**
- **[cs dhj production RETRY: ARM MATRIX](cs-dhj-production-retry.md) — nx=32 clean to 5.5 rot, nx=64 dies at 0.3 rot. Pre-fix binaries. START HERE**
- **[cs MHD instability CHARACTERIZED](cs-mhd-instability-characterized.md) — NOT gnomonic-specific; the WB source does NOT fix it. TWO RETRACTIONS (vertex order deficit, rate vs resolution). START HERE**
- **[cs MHD MINIMAL REPRODUCER: blows up from NOISE](cs-mhd-minimal-reproducer.md) — f3d35a96, cs_test iprob=13; rest atmosphere + force-free field dies in 2.5 min on one core. START HERE**
- **[cs MHD dhj blow-up: OPEN, at a CUBE VERTEX](cs-mhd-dhj-blowup.md) — 10 mechanisms ELIMINATED; refining makes it WORSE; cs dies 0.15 rot, sp clean. START HERE**
- **[cs dhj NaN: CAUSE FOUND](cs-raisevel-missing-floor.md) — GnomonicEquiangleRaiseVel never re-applied the floors. START HERE for cs+dhj**
- **[cs_dhj_long went NaN in ONE rotation](cs-dhj-long-run.md) — RESOLVED; a job that exits 0 is not a job that ran**
- **[ck_grav_prod: THE PRODUCTION CAMPAIGN](ck-grav-prod-run.md) — point-mass gravity, nx1 234, 2 GPUs, ~70-110 h. HELD, NOT LAUNCHED. START HERE**
- **[ck_mhd_b3: MHD + EOS resistivity, bbot 3 G](ck-mhd-b3-run.md) — ck_grav_prod's input, 6 lines changed. SMOKE PASSED, NOT LAUNCHED**
- [ck_limb: the CONSTANT-g run](ck-limb-run.md) — superseded by ck_grav_prod; still the calibration reference (1.76 H buffer)
- [ck_grav_size: the sizing run](ck-grav-size-run.md) — ANSWERED: r99 = 1.794e10, H = 22 cells
- [ck_hydro_long: RETIRED](ck-hydro-long-run.md) — sim "days" are Earth days; its mass drift was a transient
- [Ideal-gas + EOS x_e resistive runs](xe-resistivity-long-runs.md) — carry the stellar-heating bug; do not restart on a fixed binary
- **[dhj blow-up: CLOSED, it was the race](dhj-ck-eos-blowup.md) — no thermodynamic blow-up exists; every pre-fix arm comparison is retracted**
- [Exo-FMS cross-validation](exofms-cross-validation.md) — LW agrees 1.7%; fixed a ~25% stellar-heating loss (b4e0953c changes answers)

## Hot Jupiter: grid, EOS, atmosphere

- **[HYDROSTATIC WB scheme SUPPORTS cs + point-mass gravity; wb_option=polytropic 140dbf9e](wb-hydrostatic-scheme-cs.md) — cs == sp to 4 digits, 40-100x residual gain; recipe polytropic + wb_x1 + wb_cache_every=10 + rot_potential. Not in production**
- **[Radial grid stretch](radial-grid-stretch.md) — 1e19d4e7, `mesh/use_grid_stretch_r_poly`. Post-processing MUST apply the map to x1v**
- **[Radial stretch REFIT](radial-stretch-refit.md) — production coefficients leave a 3.7x spread in cells/H; refit gives 1.34x and 1.46x dt on cs**
- **[dt binding direction](dt-binding-direction.md) — dt is a straight MIN over x1/x2/x3; the radial stretch flips binding to radial, costing 1.76x**
- **[Point-mass gravity flag](grav-point-mass-flag.md) — `problem/grav_point_mass` (c37ebe75); constant g understates H by 2.7x at 1e-6 bar**
- **[Stellar tide flag](stellar-tide-flag.md) — c4aa730f, `problem/stellar_tide`; validated, deliberately OFF (mirror-symmetric)**
- [Stellar tide at the domain top](stellar-tide-at-domain-top.md) — missing terms are 0.5% at the limb; r_L1 = 4.28 R_p
- **[Isobar, not shell](dhj-isobar-vs-shell.md) — the 1e-6 bar level must be measured on the ISOBAR. Supersedes the floor conclusions**
- [Sizing x1max from the isobar](isobar-buffer-calibration.md) — buffer in SCALE HEIGHTS against the accepted run (1.76 H)
- [Floors for the 1e-6 bar level](dhj-floors-for-1e-6-bar.md) — lower pfloor to 1e-3; do NOT extend x1max or raise dfloor; tfloor_kelvin=50 crashes
- [Grid resolution design](dhj-grid-resolution-design.md) — measure cells/H AT the feature; the cost is the timestep, not the cells
- **[Upper-atmosphere mottling](upper-atm-mottling.md) — the speckled T maps are the H2 dissociation front plus sinking plumes, NOT noise**
- **[H2 chemistry quenches above 1e-3.6 bar](h2-chemistry-quench.md) — the EOS assumes instant dissociation; that front is not physical**
- **[sp HYDRO vs MHD at 100 rot](sp-hydro-vs-mhd-comparison.md) — the field takes 31% of the zonal KE, 89% below 20 bar; ~2% in radius. START HERE for "what does B change"**
- **[COMPOSITION + FIELD maps in the atlas](dhj-composition-maps.md) — dayside atomic/ionized, nightside molecular/neutral; plasma beta min 0.96 at tau=2/3**
- **[tau=2/3 PHOTOSPHERE diagnostic](dhj-photosphere-diagnostic.md) — a4d30d43, `problem/photosphere_dump`; surface spans 7040 km. Records the x2v CENTROID trap**
- **[EOS table dump](eos-table-dump.md) — f882159f, `<hydro>/eos_table_dump`; recovers p and T from dumps. Table is NaN below 71 K**
- **[EOS inversion NaN trap](eos-inversion-nan-trap.md) — a naive root find on the dumped table is SILENTLY 0.35 dex low in T**
- [General EOS optimization](general-eos-optimization.md) — 1.41x on CPU but only 1.062x on GPU; further work capped at ~8%
- [FastChem vs the general EOS](fastchem-vs-general-eos.md) — buys nothing, but validates the composition. DECIDED against
- [UHJ band structure in the literature](uhj-band-structure-literature.md) — 11 Kataria bins x 8 g; our grid choice is right

## Solar / stellar convection

- **[run/sun viz + WAVES NOT SHOCKS](sun-convection-viz.md) — four measures say gravity waves (max compression 0.048, 87% solenoidal); the 512^3 data is on ORION, not viper. START HERE**

## Bugs and performance

- **[ORION MERGE round-off: GATED caad9247, default path bitwise old](orion-merge-roundoff-8da093f5.md) — `floors_legacy`; the general-EOS c2p was the hunk; same-TU kernels still differ (hipcc inlining trap)**

- **[All-Mach solvers lhllc/ausmpup: algebra matches the papers](lowmach-solver-audit.md) — rk2 == rk3, NO CFL~M restriction; wb_column tgrad=-0.3 is SUPERADIABATIC (trap)**

- **[REFLECTING WALLS leaked mass under a blast: FIXED](reflect-wall-mass-leak.md) — mirror the wall state; invisible in every v_r = 0 test; sp always, cs since 979edada**
- **[GPU nondeterminism, FIXED 6e600f12](dhj-run-to-run-nondeterminism.md) — a missing team_barrier before the x1 interface-pressure floor**
- [Polar MPI host-mirror bug, FIXED](polar-mpi-host-mirror-bug.md) — 3882e37f; every multi-rank polar run died on cycle 1
- [Restart bug with tfloor_kelvin, FIXED](athenak-restart-tfloor-bug.md) — ba2f0943; general-EOS runs could never be restarted
- [dhj ideal input out-of-bounds, FIXED](dhj-ideal-input-oob.md) — 0d1f6f6a; use a Debug build for this class
- [MHD flux stale x1 limits](mhd-fluxes-stale-x1-limits.md) — 3429f59f; LATENT, changes no answer. Records the full race audit
- [RT chain-parallel split kernel](rt-chain-parallel-split.md) — 38311a8a, the live architecture; cost numbers predate the coalescing fix
- **[MAGNETIC ENERGY is not in the dumps](mhd-energy-not-from-dumps.md) — bcc is the mean of the FIELDS, not of the ENERGIES; 8.8x off in B_r. Quote the history**
- [Overall GPU profile of a dhj run](dhj-overall-gpu-profile.md) — RT 35.7%, fluxes 27.8%, ConsToPrim 12.2%, boundaries 9.6%
- [Meshblock decomposition on one GPU](meshblock-decomposition-gpu.md) — 32 blocks beats 2 by 1.35-1.45x; the second APU buys ~1%
- [nx1 ceiling from MHD LDS](nx1-ceiling-lds.md) — capped at nx1=264; x1 cannot be split because RT is a column solve
- [par_for_outer team size](par-for-outer-team-size.md) — ~3% of wall. DECIDED not worth it, do not re-propose
- RT kernel history (all HISTORICAL, superseded by [[rt-chain-parallel-split]]): [interleaving](inflight-rt-kernel-optimization.md), [band scaling](inflight-rt-band-scaling.md), [occupancy](rt-kernel-occupancy-limit.md)
- [IN FLIGHT: nx2=64 + GPU timings](inflight-nx64-timings.md) — paused 2026-08-16; sections 5/6 still to re-measure on GPU

## Working practice

- **[DELEGATE simple tasks to OPUS 5 subagents](delegate-simple-tasks-to-opus.md) — FROM THE FIRST TOOL CALL of a session; the user has had to repeat this TWICE. Fable thinks/designs; reading, builds, spec'd edits, jobs, monitoring go to `model: opus`**
- **[VALIDATE THE INSTRUMENT](validate-the-instrument.md) — a gate reporting "nothing" usually means "I did not look". READ BEFORE TRUSTING ANY CLEAN RESULT**
- **[Measure impact before claiming it](measure-impact-before-claiming.md) — report the defect, not its consequences, until an A/B has shown them**
- [Localise by dilution](localise-by-dilution.md) — a FLAT error profile does not exonerate a boundary; vary the DOMAIN SIZE at fixed dx
- [Don't propose confounded tests](confounded-tests-rejected.md) — the user vetoes experiments that could not discriminate
- **[apudev: ONE-TIME TESTS ONLY](apudev-one-time-tests-only.md) — chained/production jobs go to `apu`; user rule 2026-09-10**
- [Use fork, not origin](use-fork-not-origin.md) — all git fetch/push goes to jing-ze-ma/athenak
- **[Restart output dt trap](restart-output-dt-trap.md) — last_time comes from the restart (a dump silently never fires); output dt is SIMULATED time**
- [Never write in run/](never-write-in-run-dir.md) — 165 GB of output; read-only, and never `git add -A`
- [Viper HIP build recipe](viper-hip-build-recipe.md) — module loads and cmake flags for the MI300A APU nodes
- [Reusable CPU baseline binary](bench-baseline-worktree.md) — a worktree at bench/base_wt reproduces every baseline digit-for-digit

---

**CURRENT STATE, 2026-09-11 (~05:30), HELD on VIPER: waiting for bug-fixes from ORION; start nothing.** (was: HEAD = fork = c26b01ac on `polar-average-perf`.
Read docs/handover/HANDOVER-2026-09-11.md and [[inflight-2026-09-09-viper]] (bottom section) first. Running: sp_mhd_prod3
(bench/sp_mhd_prod3, apu chain 11569324-27), clean at rot 23.8. Nothing else. The cs collapse is CLOSED; semi-implicit RT is
the production default. Next steps are in the inflight note (dt_min 1e-3, density floor vs WB, longer si ensemble, rotate the
origin token). Rules: delegate to Opus from the first tool call; apudev one-time tests only; never write in run/.
