## Cubed sphere — start here

- **[cs RADIAL = sp RADIAL, DONE 979edada: centroid everywhere, x1 Grid-PLM always, full audit](cs-radial-unification.md) — every cs baseline before it is slightly stale; rcmfix resubmitted as 11442863**
- **[cs STRETCHED-grid SOURCE TERM off by 0.57-2.15x: FIXED 13a97399 + the full sp-vs-cs source AUDIT](cs-stretched-source-term-bug.md) — index-space dr and r in the angular-momentum curvature source; cs_test was unstretched too. Affects every stretched cs dhj run, hydro included**
- **[cs STRETCHED-grid RESISTIVITY was ANTI-DIFFUSIVE: FIXED c5c85e3b](cs-stretched-resistive-rcm-bug.md) — r_cm never got the stretch, sign flip near the top; the production dhj config. Every pre-fix cs resistive MHD result is suspect. START HERE for the deep sheet**
- **[cs seam CO-LOCATION: the binding term, FOUND and FIXED](cs-seam-colocation-fixed.md) — the clamped-window cubic; 5-11x on the evolved field. START HERE for the seam halo**
- **[cs seam order: NOTHING limits it, 2nd order CONFIRMED](cs-seam-order-limiter.md) — global L1 2.01; never read an order off the fixed-cell region bins. The seam is DONE**
- **[cs seam: de-staggering fix REFUTED](cs-seam-destag-refuted.md) — superseded by [[cs-seam-colocation-fixed]]; kept for the trail of what was refuted**
- **[resistive seam 1st order: MECHANISM CLOSED](cs-resistive-seam-order.md) — the operator is innocent (1.98), the halo INPUTS were guilty; six committed instruments**
- **[cs cross-level SEAM halo: CLOSED by two fixes](cs-crosslevel-seam-halo-first-order.md) — a01ace75 + 180a9b3e; a level boundary no longer degrades the seam. START HERE for cs+SMR halos**
- **[cs CUBE-VERTEX corner halo x RADIAL GHOST: FIXED](cs-cube-vertex-corner-radial-ghost.md) — 397b4ad3; the whole "resistivity across a radial block interface" defect. START HERE for cs halo bugs**
- **[Cubed sphere: seam conservation CLOSED](cubed-sphere-seam-conservation.md) — 985faa22 +; mass and energy exact to round-off, and history.cpp used the CARTESIAN volume on every grid**
- **[Cubed sphere: seam EMF CLOSED](cubed-sphere-seam-emf.md) — 42323a66 +; exact at 1/4/16 blocks per panel. div B is the WRONG gate. START HERE for cs+MHD**
- **[Cubed sphere: RESISTIVITY](cubed-sphere-resistivity.md) — 4bfacdd8; gnomonic two-pass curl, cs_test iprob=11, 2nd order; refinement supported since c4735f18. START HERE**
- **[Cubed sphere SMR: refined MHD CONVERGES and is MPI-clean](cubed-sphere-smr.md) — 74cbc8df; the edge FLUX buffers had no seam transform. "Rank dependence" RETRACTED. START HERE for cs+SMR**
- **[Cubed sphere MHD convergence](cubed-sphere-mhd-convergence.md) — e63b571a; mhd_corner_e had no cs form. Residual is the radial-BC phase lag, one law in dt**
- **[SHOCKS through cs seams](cs-shocks-through-seams.md) — structure ON a seam NaN'd everything SILENTLY; the unlimited along-seam resample, now clamped. FOFC DECIDED not needed. START HERE**
- **[cs blast vs a CARTESIAN grid](cs-blast-vs-cartesian.md) — the seam costs NOTHING (4.1% vs Cartesian's own 4.5%); only the cube VERTEX is ~2x worse. START HERE for "is cs good enough for shocks"**

## Cubed sphere — validation and limits

- **[cs_test FF decay + rot_axis (4990eb41): cs matrix vs sp, no 1.5-order region on cs](cs-test-ffdecay-rotaxis.md) — START HERE for cs-vs-sp simple tests**
- **[cs PURE HYDRO validation](cs-hydro-validation.md) — mass and energy to MACHINE PRECISION even under a shock; space 2.3-2.4 L1. Records the FALSE first-order reading. START HERE**
- **[cs MHD + RESISTIVE validation, BY REGION](cs-mhd-validation.md) — RE-MEASURED: ideal MHD 2nd order EVERYWHERE (the old 1.84/1.76 was a radial floor); the resistive energy drift is PHYSICAL Ohmic heating, NO BUG**
- **[cs ANGULAR MOMENTUM: why sp is exact and cs cannot be; RE-MEASURED 09-06](cs-angular-momentum.md) — rigid rotation 6.9e-4/rot at nx2=32 (~2.7th order); in the dhj runs the rotating-frame Lz tracks sp to 1e-4 of the frame L over 100 rot**
- **[cs: radiation + srcterms are now REFUSED](cs-unsupported-physics-guards.md) — b283ad3a; they used to run and return a wrong answer. Two traps that make a guard which can never fire**
- **[cs NARROW-BLOCK resample degeneracy](cs-narrow-block-resample-degeneracy.md) — the along-seam stencil INVERTS below 3 cells, i.e. any MeshBlock under 6. Test at the MINIMUM legal block size**
- **[cs GENERAL/TABULATED EOS stale cache: FIXED](cs-general-eos-stale-cache.md) — p/Gamma_1/T cached BEFORE the gnomonic correction, 6-15% wrong. The gate is VACUOUS without a flow**
- **[cs MHD blast from a VECTOR POTENTIAL](cs-mhd-blast.md) — b283ad3a; div B zero by construction. Gate the construction with TWO numbers — div B alone is satisfied by a ZERO field**
- **[cube-vertex corner fill: PROMOTED, default ON](cs-wire-fill-wip.md) — branch cs-wire-wip, `<mesh>/cs_vertex_fill`; ownership is the exact chart DIAGONAL. Merged; widening coarse/fine MEASURED not worth it**
- **[cs cube-vertex REAL fill: prototyped, NOT built](cs-cube-vertex-real-exchange.md) — sampling beats extrapolation 3.4-5.5x and needs NO new exchange, just a widened destination range**

## Cubed sphere — history

- **[Cubed sphere ON GPU](gpu-this-capture-device-lambda.md) — af539941; hydro and MHD match CPU. The constant-memory CLOSURE LIMIT; incremental builds hide warnings**
- **[cs GPU multi-block fault: FIXED](cs-gpu-multiblock-fault.md) — 4302a008; one kernel captured both send and recv buffers. Second over-large-functor bug that session**
- **[Cubed sphere + MPI: FIXED](cubed-sphere-mpi-hang.md) — 1c2e29d6; unmatched receive, cube-vertex cause. Run the RECIPROCITY AUDIT first for any boundary bug**
- **[Cubed sphere: MHD seam FIXED](cubed-sphere-mhd-seam.md) — 4f19a244; three plumbing bugs, the transform was always exact**
- **[Cubed sphere: MHD](cubed-sphere-mhd.md) — e0c74357; dxedge was zero, mhd_fluxes had no gnomonic rotation, bcc was a non-orthogonal triple**
- **[Cubed sphere: COMMITTED as 9492a946](cubed-sphere-committed.md) — five sessions of hydro work is in git**
- **[Cubed sphere: along-seam resample, 2nd ORDER](cubed-sphere-seam-interp.md) — the "flat interior residual" is RETRACTED**
- **[Cubed sphere: seam basis transform, FIXED](cubed-sphere-seam-basis.md) — a tangent-BASIS transform, not a signed permutation**
- [Cubed sphere panel frames](cubed-sphere-panel-frames.md) — the six frames are forced by panel_neighbors; a duplicate copy had 3/4 swapped
- [Cubed-sphere hydro fix](cubed-sphere-hydro-fix.md) — halo axis bug FIXED; its numbers predate the panel 3/4 swap
- [Cubed sphere: x1 as radial](cubed-sphere-x1-radial.md) — DONE, bit-identical; axis roles cyclically permuted
- [Cubed sphere: hydro state](cubed-sphere-hydro-state.md) — 4 early bugs FIXED; its rigid-rotation section is RETRACTED
- [Cubed sphere for the hot Jupiter](cubed-sphere-for-hot-jupiter.md) — the dt argument is DEAD (radial binds)

## Hot Jupiter: physics, runs, campaigns

- **[cs_mhd_prod DIED at rot 41.6: NaN everywhere in one interval, no precursor, ran on NaN to rot 54](cs-mhd-prod-nan-rot41.md) — NaN-hunt restart 11529359 (0.05-rot dumps from rot 40). The code does not stop on NaN. START HERE for the old cs production**
- **[cs_hyd_rs: pure-hydro cs twin x {hllc, lhllc, ausmpup, ppmx} for the jet comparison](cs-hyd-rs-run.md) — launched 09-07 08:25, jobs 11529253/56/59 + chains. START HERE for the solver comparison**
- **[cs_mhd_prod2: cs MHD production FROM SCRATCH on the 09-07 defaults (WB+rot_potential+deep RT)](cs-mhd-prod2-run.md) — lead 11526708 (05:19, 6.5 h), chain 11526681/2; A/B against cs_mhd_prod. START HERE for the current production**

- **[Deep RT: radiative conduction + the tau_R BLEND fad2f5db](radiative-conduction-deep-interior.md) — `rad_tau_lo/hi`; sp column: dt+cost equal, RT chain SHORTER, 1-cell cut spread; the sp dt collapse (angular widths) FIXED; `rad_kappa_src = table` (kappa_R tabulated from the ck table at start-up) makes the handover CONSISTENT to 10 %; Freedman blend was off 1.7x. START HERE for deep RT**
- **[cs MHD dhj VALIDATED at rot 20: T structure = sp, jet WEAK (~1 km/s) on every grid, B does NOT slow it (beta>=100)](cs-dhj-diagnostics-rot20.md) — report URL + the EOS-inversion / cs lat-lon / wind-projection pipeline and its traps. START HERE for any cs dhj science plot**
- **[sp_mhd_prod (18f5dd21) DIES at rot 10.3: polar-row RADIAL field at the BOTTOM; BISECTED: 2a64e2c7 CLEAN, ideal 18f5dd21 DIRTY, resistivity innocent](sp-pole-bottom-radial-blowup.md) — round 2 arms wall_r (4990eb41) + x3shift_r (863e8337), jobs 11526902/3, gate ME1 at rot 4.7. START HERE for sp MHD**
- **[cs MHD dhj BOTTOM-BOUNDARY DRIFT: FOUND + FIXED 248f1b77](cs-mhd-bottom-inflow.md) — the radial BC built cs ghost bcc as plain face averages, the ghost energy drifted every call (8 % by rot 16), bottom cells drained, NaN at rot 32. EVERY cs MHD dhj run before it is suspect; redo from scratch. START HERE**
- **[sp MHD ENERGY EXCESS vs hydro (+1 %/2 rot, 3-4x horizontal KE), ABSENT on cs](sp-mhd-energy-excess.md) — floors on the upper NIGHTSIDE where sp carries 7x the field of cs; deep field identical; every sp MHD variant has it; cause OPEN. READ before comparing sp MHD with anything**
- **[Correlated-k design + build](correlated-k-design.md) — THE MAIN THREAD. COMPLETE: physics, docs, tests, MPI. START HERE**
- **[dhj RT on the CUBED SPHERE: WORKING](dhj-cubed-sphere-port.md) — the full production physics tracks spherical polar; four bugs found, three PRE-EXISTING**
- **[cs NON-ORTHOGONAL audit: BOTH FIXED](cs-nonorthogonal-audit.md) — e415b91a; history.cpp KE/ME were ~1% wrong on cs, newdt understated the fast speed by 1/sin_cell. Pre-fix cs history files carry that 1%**
- **[A real cs MHD bug, FIXED: the C2P floor corrupts u.e](cs-mhd-c2p-floor-corrupts-ue.md) — ConsToPrim's magnetic energy uses the NON-ORTHOGONAL triple; a floor then OVERWRITES u.e. Vertex order -0.04 -> +1.97. But it does NOT save the dhj run (0.2023 -> 0.2129) — that blow-up stays OPEN**
- **[cs MHD low beta: per-step operator 2nd ORDER everywhere incl. the vertex](cs-mhd-low-beta-divergent.md) — the instability is NONLINEAR on a consistent scheme; HLLE fallback is a legitimate stabiliser. START HERE**
- **[cs MHD low-beta: fixed ON THE TEST PROBLEMS ONLY](cs-mhd-lowbeta-fix.md) — cb1afd28, `<mhd>/cs_lowbeta_llf`, a per-face HLLD->HLLE fallback below beta 0.5, default ON for cs. Real but INSUFFICIENT: the production dhj run still dies at 0.21 rot**
- **[cs SEAM GHOST e_int: RETRACTED as a mechanism](cs-seam-ghost-eint.md) — the seam degrades EVERY quantity by the same ~15-19x; eint is not special, the 1/beta signature was in the null row too. Real E-vs-B inconsistency is ~14 %. Fix b9efbbfa single-rank only. READ THE RETRACTION FIRST**
- **[GENERAL-EOS AUDIT](general-eos-audit.md) — the (rho,e)->T inversion is BRACKETED, monotone, warm-start-independent; no new live bug. Off-table states are silently log-linear-extrapolated**
- **[3D RESISTIVITY BUG since 2026-06-20: curl B missing its theta terms on sp + Cartesian, FIXED (uncommitted)](resistivity-3d-curl-missing-terms.md) — every sp resistive dhj run since June is wrong; cs unaffected. START HERE for any sp resistive claim**
- **[sp POLE-EDGE CURRENT: far-side loop segment had the WRONG SIGN, J_r grew ~1/dtheta, FIXED](sp-pole-edge-current-sign.md) — the m=2 polar pressure pattern; new edge-by-edge EMF check in sp_test iprob=11**
- **[sp POLE: resistive J_r divided by a ZERO dual area, FIXED](sp-pole-edge-area-zero.md) — NaN in cycle 1 on any uniform-theta grid; production hid it with the theta stretch**
- **[sp POLE FIXED: third-difference polar_emf_diss + x3 face-state theta shift](sp-pole-fixes.md) — de667f32/863e8337; polar-row L1(B) rate 0.76 -> 2.1, v_r 2.3; tangential stays ~1.5 (cart REPAIRED but no better, 4th-order faces REJECTED). dhj checkerboard gate PASSED to rot 1.2. START HERE for the pole**
- **[sp_test: RIGID ROTATION + RESISTIVE DECAY tests + rates](sp-test-rigidrot-resist.md) — hydro 2.6, resistive B 2.1; the POLAR ROWS do not converge for a tangential field (O(1))**
- **[RESISTIVITY AUDIT](resistivity-audit.md) — eta_b ghosts, resistive dt, and the pole-face curl are CLEAN (with file:line); the one live gap is the polar EMF average excluding the resistive EMF, fixed behind `use_polar_average_eresist`**
- **[sp POLAR FIELD BLOW-UP: FOUND + FIXED 1cabe85c `<mhd>/polar_emf_diss`, DEFAULT ON 3147de4f](sp-polar-field-blowup.md) — GS05 corner-EMF upwind terms grow a (-1)^(i+k) face checkerboard in the polar row. READ BEFORE ANY cs-vs-sp MHD CLAIM**
- **[cs ROTATION SOURCE BUG: cs ran the CARTESIAN beta-plane, sp ran full Coriolis](cs-rotation-source-bug.md) — f75ad783; the switch was REMOVED in ae831665, full rotation is now unconditional. "Same physics, cs dies and sp does not" was FALSE. READ THIS BEFORE TRUSTING ANY cs-vs-sp COMPARISON**
- **[cs 1-ULP amplification: nx=64 was DYNAMICALLY UNSTABLE on PRE-FIX binaries; RETEST cs_n64_fixed 11517116 ALIVE past 0.5 rot, dt flat 2.3 s -- old claim STALE](cs-ulp-amplification.md) — round-off grows to O(1) and KILLS nx=64 at 0.309 rot while it saturates ~1% at nx=32; the perturbation script must skip the TREE metadata. START HERE**
- **[WHAT RECONSTRUCTION ACTUALLY RUNS: sp ignores `reconstruct` in ALL directions (stretched PLM hard-wired); cs honours it but x1 ignores the radial stretch](reconstruction-what-actually-runs.md) — READ before any reconstruction claim**
- **[cs VERTEX = LIMITER CLIPPING: ppmx (the existing extremum-preserving PPM) and WENO-Z cut the worst-cell vertex force 6.5x; the Cartesian update (d166238b) does nothing there](cs-vertex-limiter-clipping.md) — same root cause as the sp pole; ppmx is the production candidate (nghost 4), gates in memory**
- **[cs VERTEX ORACLE: exact ghosts change the vertex residual <1 %; it is the CELL BALANCE, no ghost stencil can help](cs-vertex-oracle-halo-innocent.md) — the through-vertex seam-continuation idea DECIDED against; Cartesian-momentum update is the only lever left**
- **[sp POLAR-ROW RESIDUAL FIXED 7.7x: `mesh/polar_quadratic_recon` 6410be99, default OFF](sp-polar-row-reconstruction.md) — limiters see the pole as an extremum; `reconstruct` is IGNORED on sp. READ before any sp accuracy claim**
- **[sp GEOMETRIC-SOURCE RESIDUAL: polar rows 14x interior; face-sum WB source (default OFF) does NOT touch it; Cartesian polar-row update 34882bd7 (default OFF) fixes the transverse-field part only](sp-geometric-source-residual.md) — my "WB will fix it" RETRACTED**
- **[cs WELL-BALANCED SOURCE CACHED: 46 % -> 4 % of a run, bitwise identical, restart-clean](cs-wb-source-cached.md) — 88064f67; records the tlim-clipping restart-test TRAP and that cs_test cannot restart**
- **[GS07 corner EMF is UNCONDITIONAL on cs since b71ef392 (switch REMOVED)](cs-gs07-emf-gate.md) — gate: same order as the plain average, 1.8x its smooth-field L1(B); `bs_emf` is the only way back, for diagnosis**
- **[cs DEEP TOROIDAL SHEET: CLOSED, it was the stretched resistive bug c5c85e3b](cs-deep-toroidal-sheet.md) — rcmfix on 979edada has NO sheet; the old cs_prod_mhd_rot (sheet, eruption, death) is archived. START HERE**
- **[cs dhj production RETRY: ARM MATRIX + production reruns](cs-dhj-production-retry.md) — nx=32 clean to 5.5 rot, nx=64 dies at 0.3 rot (1-ULP seed amplifies). Pre-fix binaries. START HERE**
- **[cs MHD instability CHARACTERIZED](cs-mhd-instability-characterized.md) — a6401406, df2edc8e, c3992145. NOT gnomonic-specific (spherical polar does it too). The WELL-BALANCED source is BUILT and does NOT fix it, refuting the quadrature-mismatch hypothesis. TWO RETRACTIONS: the vertex order deficit and "the rate rises with resolution". START HERE**
- **[cs MHD MINIMAL REPRODUCER: it blows up from NOISE](cs-mhd-minimal-reproducer.md) — f3d35a96, cs_test iprob=13. An atmosphere AT REST + a force-free field dies in 2.5 min on ONE CORE, no RT/rotation/EOS/resistivity/shear. Growth rate RISES with resolution. Halo and stratification EXONERATED by A/B. START HERE**
- **[cs MHD dhj blow-up: OPEN, at a CUBE VERTEX](cs-mhd-dhj-blowup.md) — 10 mechanisms ELIMINATED incl. timestep, CT/monopoles, seam transform+resample (exact), block decomposition, general EOS, radial stretch. Refining makes it WORSE. Minimal reproducer: ideal gas + ideal MHD + uniform grid, cs dies 0.15 rot, sp clean. START HERE**
- **[cs dhj NaN: CAUSE FOUND](cs-raisevel-missing-floor.md) — GnomonicEquiangleRaiseVel never re-applied the floors after the metric correction. START HERE for cs+dhj**
- **[cs_dhj_long went NaN in ONE rotation](cs-dhj-long-run.md) — RESOLVED; kept for the 2x2 that pointed at the grid. A job that COMPLETES with exit 0 is not a job that ran**
- **[ck_grav_prod: THE PRODUCTION CAMPAIGN](ck-grav-prod-run.md) — point-mass gravity, nx1 234 stretched, 2 GPUs, ~70-110 h. HELD, NOT LAUNCHED. START HERE**
- **[ck_mhd_b3: MHD + EOS resistivity, bbot 3 G](ck-mhd-b3-run.md) — ck_grav_prod's input, 6 lines changed. SMOKE PASSED, NOT LAUNCHED**
- [ck_limb: the CONSTANT-g run](ck-limb-run.md) — superseded as science by ck_grav_prod; still the calibration reference (1.76 H buffer)
- [ck_grav_size: the sizing run](ck-grav-size-run.md) — ANSWERED: r99 = 1.794e10, H = 22 cells
- [ck_hydro_long: RETIRED](ck-hydro-long-run.md) — sim "days" are Earth days; its mass drift was a transient
- [Ideal-gas + EOS x_e resistive runs](xe-resistivity-long-runs.md) — COMPLETED but carry the stellar-heating bug; do not restart on a fixed binary
- **[dhj blow-up: CLOSED, it was the race](dhj-ck-eos-blowup.md) — no thermodynamic blow-up exists; every pre-fix arm comparison is retracted**
- [Exo-FMS cross-validation](exofms-cross-validation.md) — LW agrees 1.7%; found and fixed a ~25% stellar-heating loss (b4e0953c changes production answers)

## Hot Jupiter: grid, EOS, atmosphere

- **[HYDROSTATIC well-balanced scheme SUPPORTS cs + point-mass gravity; NEW wb_option=polytropic 140dbf9e](wb-hydrostatic-scheme-cs.md) — cs == sp to 4 digits; 40-100x residual gain, floor = the IC's own discrete inconsistency; avg-anchor correction REFUTED; table cost 1.03x of plain (polytropic + wb_cache_every=10); VERIFIED hydro/MHD/resistive on cs, sp, Cartesian x1/x3, RT on, rotation on (WB = the CONVERGED answer), GPU == CPU to 1e-11. Recipe: polytropic + wb_x1 + wb_cache_every=10 + problem/rot_potential=true (barotropic rotating IC). Not yet in production**
- **[Radial grid stretch](radial-grid-stretch.md) — 1e19d4e7, `mesh/use_grid_stretch_r_poly`. Post-processing MUST apply the map to x1v**
- **[Radial stretch REFIT](radial-stretch-refit.md) — the production coefficients leave a 3.7x spread in cells/H; refit vs H_rho gives 1.34x and 1.46x dt on cs**
- **[dt binding direction](dt-binding-direction.md) — dt is a straight MIN over x1/x2/x3; the radial stretch flips binding from azimuthal to radial, costing 1.76x**
- **[Point-mass gravity flag](grav-point-mass-flag.md) — `problem/grav_point_mass` (c37ebe75); constant g understates H by 2.7x at 1e-6 bar**
- **[Stellar tide flag](stellar-tide-flag.md) — c4aa730f, `problem/stellar_tide`; validated, deliberately OFF (mirror-symmetric between limbs)**
- [Stellar tide at the domain top](stellar-tide-at-domain-top.md) — the missing terms are 0.5% at the limb, +6.8 cells substellar; r_L1 = 4.28 R_p
- **[Isobar, not shell](dhj-isobar-vs-shell.md) — the 1e-6 bar level must be measured on the ISOBAR. Supersedes the floor conclusions**
- [Sizing x1max from the isobar](isobar-buffer-calibration.md) — buffer in SCALE HEIGHTS against the accepted run (1.76 H), plus two mu validations
- [Floors for the 1e-6 bar level](dhj-floors-for-1e-6-bar.md) — lower pfloor to 1e-3 (free); do NOT extend x1max or raise dfloor; tfloor_kelvin=50 crashes
- [Grid resolution design](dhj-grid-resolution-design.md) — measure cells/H AT the feature; the cost is the timestep, not the cells
- **[Upper-atmosphere mottling](upper-atm-mottling.md) — the speckled T maps are the H2 dissociation front plus sinking plumes, NOT noise or floors**
- **[H2 chemistry quenches above 1e-3.6 bar](h2-chemistry-quench.md) — the EOS assumes instant dissociation; recombination is 8 decades too slow, so that front is not physical**
- **[sp HYDRO vs MHD, controlled at 100 rot](sp-hydro-vs-mhd-comparison.md) — the field takes 31% of the zonal KE, 89% of it below 20 bar, and shows at the photosphere as ~2% in radius. START HERE for "what does B change"**
- **[COMPOSITION + FIELD maps in the atlas](dhj-composition-maps.md) — dayside atomic + iron ionized, nightside molecular + neutral; the atlas now carries BOTH runs, with log|B|, B_r and plasma beta (min 0.96 at tau=2/3)**
- **[tau=2/3 PHOTOSPHERE diagnostic](dhj-photosphere-diagnostic.md) — a4d30d43, `problem/photosphere_dump`; the emitting surface spans 7040 km between near- and far-IR. Records the x2v CENTROID trap**
- **[EOS table dump](eos-table-dump.md) — f882159f, `<hydro>/eos_table_dump`; recovers p and T from binary dumps. Table is NaN below 71 K**
- **[EOS inversion NaN trap](eos-inversion-nan-trap.md) — a naive root find on the dumped table is SILENTLY 0.35 dex low in T. Validate against ck_dump_file**
- [General EOS optimization](general-eos-optimization.md) — 1.41x on CPU but only 1.062x on GPU; further EOS work capped at ~8%
- [FastChem vs the general EOS](fastchem-vs-general-eos.md) — MEASURED: buys nothing, but validates the composition. DECIDED against
- [UHJ band structure in the literature](uhj-band-structure-literature.md) — Parmentier+2018 and Tan+2024 both use 11 Kataria bins x 8 g; our grid choice is right

## Solar / stellar convection

- **[run/sun viz + WAVES NOT SHOCKS](sun-convection-viz.md) — the "Buoyancy Box" artifact; tau=2/3 granulation from the run's OWN opacity; four measures say gravity waves (max compression 0.048, 87% solenoidal). The 512^3 convection data is GONE from viper (it is on ORION). START HERE**

## Bugs and performance

- **[All-Mach solvers lhllc/ausmpup: algebra matches the papers; rk2 == rk3, NO CFL~M restriction on the stratified rest test](lowmach-solver-audit.md) — two design choices (normal-velocity chi; AUSM Mref=1 hard-wired). wb_column tgrad=-0.3 is SUPERADIABATIC (trap)**

- **[REFLECTING WALLS leaked mass under a blast on sp (always) + cs (since 979edada): FIXED, mirror the wall state](reflect-wall-mass-leak.md) — invisible in every v_r = 0 test; closed blasts back to round-off on both grids**
- **[GPU nondeterminism, FIXED 6e600f12](dhj-run-to-run-nondeterminism.md) — a missing team_barrier before the x1 interface-pressure floor; new _gpu regression test**
- [Polar MPI host-mirror bug, FIXED](polar-mpi-host-mirror-bug.md) — 3882e37f; every multi-rank polar run died on cycle 1
- [Restart bug with tfloor_kelvin, FIXED](athenak-restart-tfloor-bug.md) — ba2f0943; general-EOS runs could never be restarted
- [dhj ideal input out-of-bounds, FIXED](dhj-ideal-input-oob.md) — 0d1f6f6a; use a Debug build (Kokkos bounds checking) for this class
- [MHD flux stale x1 limits](mhd-fluxes-stale-x1-limits.md) — 3429f59f; LATENT, changes no answer. Also records the full race audit
- [RT chain-parallel split kernel](rt-chain-parallel-split.md) — 38311a8a, the live architecture; its cost numbers predate the coalescing fix
- **[MAGNETIC ENERGY is not in the dumps](mhd-energy-not-from-dumps.md) — bcc is the mean of the FIELDS, the code integrates the mean of the ENERGIES; 8.8x off in B_r. Quote the history**
- [Overall GPU profile of a dhj run](dhj-overall-gpu-profile.md) — RT 35.7%, fluxes 27.8%, ConsToPrim 12.2%, boundaries 9.6%
- [Meshblock decomposition on one GPU](meshblock-decomposition-gpu.md) — 32 blocks beats 2 by 1.35-1.45x; the second APU buys ~1%
- [nx1 ceiling from MHD LDS](nx1-ceiling-lds.md) — radial resolution capped at nx1=264; x1 cannot be split because RT is a column solve
- [par_for_outer team size](par-for-outer-team-size.md) — ~3% of wall. DECIDED not worth it, do not re-propose
- RT kernel history (all HISTORICAL, superseded by [[rt-chain-parallel-split]]): [interleaving](inflight-rt-kernel-optimization.md), [band scaling](inflight-rt-band-scaling.md), [occupancy](rt-kernel-occupancy-limit.md)
- [IN FLIGHT: nx2=64 + GPU timings](inflight-nx64-timings.md) — paused 2026-08-16; sections 5/6 still to re-measure on GPU

## Working practice

- **[VALIDATE THE INSTRUMENT](validate-the-instrument.md) — a gate reporting "nothing" usually means "I did not look". READ THIS BEFORE TRUSTING ANY CLEAN RESULT**
- **[Measure impact before claiming it](measure-impact-before-claiming.md) — report the defect, not its consequences, until an A/B has shown them; two retractions**
- [Localise by dilution](localise-by-dilution.md) — a FLAT error profile does not exonerate a boundary; vary the DOMAIN SIZE at fixed dx and look for 1/L
- [Don't propose confounded tests](confounded-tests-rejected.md) — the user vetoes experiments whose result could not discriminate
- [Use fork, not origin](use-fork-not-origin.md) — all git fetch/push goes to jing-ze-ma/athenak
- **[Restart output dt trap](restart-output-dt-trap.md) — last_time comes from the restart (a dump silently never fires), and output dt is SIMULATED time (38 GB in two minutes)**
- [Never write in run/](never-write-in-run-dir.md) — 165 GB of output; read-only, and never `git add -A`
- [Viper HIP build recipe](viper-hip-build-recipe.md) — module loads and cmake flags for the MI300A APU nodes
- [Reusable CPU baseline binary](bench-baseline-worktree.md) — a worktree at bench/base_wt reproduces every baseline digit-for-digit; the other bench binaries are GPU builds

---

**CURRENT STATE, 2026-09-07 (~09:45), HANDOVER TO ORION.** HEAD 7c652768 on `polar-average-perf`, PUSHED
to the fork (jing-ze-ma/athenak). Viper is in MAINTENANCE 09-07 12:00 -> 09-12 12:00; the next
session runs on ORION, which cannot see viper's bench/ or scratch. Everything needed is in git:
docs/handover/HANDOVER-2026-09-07.md (READ FIRST), docs/handover/scripts/ (the analysis scripts),
docs/handover/claude-memory-2026-09-07/ (a copy of this memory directory as of the handover).

### Running on VIPER through/after the maintenance (check when back on viper, NOT from orion)
1. cs_mhd_prod2 (11526708 -> chain 11526681/2): cs MHD production on the new defaults, rot ~11 at 09:00.
2. cs_hyd_rs/{hllc,lhllc,ausmpup,ppmx} (11529253/56/59/11529330 + chains): pure-hydro solver and
   reconstruction comparison; the jet at rot 15-20 is the question. LHLLC keeps 70x the deep radial KE.
3. sp_pole_bisect wall_r (11526902) and x3shift_r (11526903): gate ME1 at rot 4.7 / polerow at rot 4.
4. cs_mhd_prod_nan (11529359): passed rot 41.6 without NaN -> the original NaN is not reproducible.
5. Login-node serial: lowmach/isoL_N64_* (2000 t_BV Edelmann slow-mode test), ~50 % at 09:00.

### NEXT STEPS on ORION (in order)
1. Clone/pull the fork branch; `git submodule update --init`; build CPU (PROBLEM=deep_hot_jupiter_rt) and GPU
   per viper-hip-build-recipe's orion equivalent; the ck tables are at /orion/u/jinma/ATHENAK/athenak/data/exo_fms_ck.
2. Regression: `cd tst && python run_test_suite.py --test test_suite/rad/test_rad_dhj_ck_cpu.py` and
   test_rad_cs_raddiff_cpu.py (both PASS on viper at ebd57244).
3. GPU binary needs the NaN guard (d3d74f2b) -- rebuild before any GPU run.
4. Determinism test (cs MHD restart twice from one rst, 200 cycles, bitwise) -- the rot-41.6 NaN hunt showed
   restarts diverge at 5e-6 in 0.1 rot ([[cs-mhd-prod-nan-rot41]]).
5. Open science: the jet comparison (results on viper after 09-12), the deep interior with the radiative
   outlet, and the sp polar blow-up bisection result.
