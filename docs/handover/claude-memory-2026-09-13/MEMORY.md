## RULE ZERO

- **[SAVE TOKENS](save-tokens.md) — no unneeded work; narrow prompts, cheap models for lookups**

- **[DELEGATE simple tasks to OPUS 5 subagents](delegate-simple-tasks-to-opus.md) — from the first tool call; Fable decides, Opus executes**
- **[STANDING RULE: delegate to Opus 5 agents](delegate-heavy-work-to-opus.md) — EVERY session from the first tool call: Opus agents do reading/editing/building/launching/MONITORING/analysis; I do diagnosis, briefs, VERIFY every diff and number; ONE deliverable per agent, name the hunks to read, "STOP, no timers" (token blow-up 09-11)**
- **[SAVE TOKENS (standing rule 09-11): especially on Fable; no unnecessary tool calls, no surveys, agents one deliverable each, no timers left armed](save-tokens-everywhere.md)**

## Cubed sphere — start here

- **[RT SEMI-IMPLICIT source](rt-semi-implicit-changes-dhj-answer.md) — 048dff30; explicit overcools, semi-implicit is right**
- **[IN FLIGHT, STOPPED 09-11](inflight-2026-09-09-viper.md) — cs collapse closed; c26b01ac pushed**
- **[GPU BUILD RACE: flag never compiled](gpu-build-race-flag-not-compiled.md) — verify flags with `strings` plus a startup print**
- **[cs RADIAL = sp RADIAL](cs-radial-unification.md) — 979edada; earlier cs baselines stale**
- **[cs STRETCHED-grid SOURCE TERM off 0.57-2.15x](cs-stretched-source-term-bug.md) — fixed 13a97399, index-space dr**
- **[cs STRETCHED-grid RESISTIVITY was ANTI-DIFFUSIVE](cs-stretched-resistive-rcm-bug.md) — fixed c5c85e3b, r_cm missed the stretch**
- **[cs seam CO-LOCATION found and fixed](cs-seam-colocation-fixed.md) — the clamped-window cubic, 5-11x**
- **[cs seam order: 2nd order confirmed](cs-seam-order-limiter.md) — global L1 2.01, done**
- **[cs seam de-staggering fix REFUTED](cs-seam-destag-refuted.md) — superseded**
- **[resistive seam 1st order: closed](cs-resistive-seam-order.md) — operator innocent, halo inputs guilty**
- **[cs cross-level SEAM halo: closed](cs-crosslevel-seam-halo-first-order.md) — a01ace75 + 180a9b3e**
- **[cs CUBE-VERTEX corner halo x radial ghost: fixed](cs-cube-vertex-corner-radial-ghost.md) — 397b4ad3**
- **[cs seam conservation closed](cubed-sphere-seam-conservation.md) — 985faa22, exact**
- **[cs seam EMF closed](cubed-sphere-seam-emf.md) — 42323a66; div B is the wrong gate**
- **[cs RESISTIVITY](cubed-sphere-resistivity.md) — 4bfacdd8; two-pass curl, 2nd order**
- **[cs SMR: refined MHD converges, MPI-clean](cubed-sphere-smr.md) — 74cbc8df; rank dependence retracted**
- **[cs MHD convergence](cubed-sphere-mhd-convergence.md) — e63b571a; residual is radial-BC phase lag**
- **[SHOCKS through cs seams](cs-shocks-through-seams.md) — resample clamped; FOFC not needed**
- **[cs blast vs a CARTESIAN grid](cs-blast-vs-cartesian.md) — seam costs nothing, vertex ~2x worse**

## Cubed sphere — validation and limits

- **[cs_test FF decay + rot_axis](cs-test-ffdecay-rotaxis.md) — 4990eb41; no 1.5-order region**
- **[cs PURE HYDRO validation](cs-hydro-validation.md) — machine precision under a shock; L1 2.3-2.4**
- **[cs MHD + RESISTIVE validation by region](cs-mhd-validation.md) — 2nd order; energy drift is physical Ohmic heating**
- **[cs ANGULAR MOMENTUM: sp exact, cs cannot be](cs-angular-momentum.md) — 6.9e-4/rot at nx2=32**
- **[cs: radiation + srcterms now REFUSED](cs-unsupported-physics-guards.md) — b283ad3a**
- **[cs NARROW-BLOCK resample degeneracy](cs-narrow-block-resample-degeneracy.md) — stencil inverts below 3 cells**
- **[cs GENERAL/TABULATED EOS stale cache: fixed](cs-general-eos-stale-cache.md) — cached pre-correction, 6-15% wrong**
- **[cs MHD blast from a VECTOR POTENTIAL](cs-mhd-blast.md) — b283ad3a; gate with two numbers**
- **[cube-vertex corner fill promoted, default ON](cs-wire-fill-wip.md) — `<mesh>/cs_vertex_fill`**
- **[cs cube-vertex REAL fill prototyped, not built](cs-cube-vertex-real-exchange.md) — sampling beats extrapolation 3.4-5.5x**
- [cs determinism test on orion](cs-determinism-test-orion.md) — DONE 193377: CPU restarts are BIT-IDENTICAL, even 16x7 vs 28x4; the viper 5e-6 divergence is GPU-specific or a binary difference
- [cs orthogonal-KE AUDIT: red_giant ghost fills + IC + inner passes FIXED in the working tree 09-09 (T14: correct, not the killer); latent sites listed in docs/handover/NOTE-2026-09-09-cs-covariant-basis.md (pushed 40aef0ad)](cs-orthogonal-ke-audit-2026-09-09.md)

## Cubed sphere — history

- **[Cubed sphere ON GPU](gpu-this-capture-device-lambda.md) — af539941; constant-memory closure limit**
- **[cs GPU multi-block fault: fixed](cs-gpu-multiblock-fault.md) — 4302a008**
- **[Cubed sphere + MPI: fixed](cubed-sphere-mpi-hang.md) — 1c2e29d6, unmatched receive at a vertex**
- **[cs MHD seam fixed](cubed-sphere-mhd-seam.md) — 4f19a244; transform always exact**
- **[Cubed sphere: MHD](cubed-sphere-mhd.md) — e0c74357; dxedge was zero**
- **[Cubed sphere committed](cubed-sphere-committed.md) — 9492a946**
- **[cs along-seam resample, 2nd order](cubed-sphere-seam-interp.md) — flat-residual claim retracted**
- **[cs seam basis transform fixed](cubed-sphere-seam-basis.md) — tangent-basis transform**
- [Cubed sphere panel frames](cubed-sphere-panel-frames.md) — duplicate copy had 3/4 swapped
- [Cubed-sphere hydro fix](cubed-sphere-hydro-fix.md) — halo axis bug fixed
- [Cubed sphere: x1 as radial](cubed-sphere-x1-radial.md) — done
- [Cubed sphere: hydro state](cubed-sphere-hydro-state.md) — 4 bugs fixed; rigid-rotation part retracted
- [Cubed sphere for the hot Jupiter](cubed-sphere-for-hot-jupiter.md) — dt argument dead

## Hot Jupiter: physics, runs, campaigns

- **[cs WB arm died rot 12.3; KE gap is sp MHD, not cs](cs-wb-arm-and-ke-gap.md) — every cs arm = sp hydro = 2-3e34**
- **[dhj JET is SHALLOW, deep equator WESTWARD to rot 164 on sp](dhj-jet-shallow-westward-deep.md) — all cs solvers = sp; hllc jet half of ausm+up**
- **[cs VERTEX dt COLLAPSE closed 09-11](cs-vertex-dt-collapse-0907-defaults.md) — explicit RT source overcools; semi-implicit is right**
- **[RESTART dropped rot_potential: centrifugal force doubled](restart-rot-potential-bug.md) — fixed; first-round ablation arms invalid**
- **[WB RESTART cache bug fixed 5c0b98e4](wb-restart-cache-bug.md) — zero background for 9 cycles**
- **[sp POLAR BLOW-UP bisected to 863e8337 and fixed](sp-pole-bisect-culprit-863e8337.md) — `mesh/polar_x3_shift=rotate`**
- [sp pole bisect round 2: null](sp-pole-bisect-round2-null.md)
- **[cs_mhd_prod died at rot 41.6](cs-mhd-prod-nan-rot41.md) — NaN with no precursor**
- **[cs_hyd_rs: hydro cs twin x 4 solvers](cs-hyd-rs-run.md) — jet comparison**
- **[cs_mhd_prod2: cs MHD production from scratch](cs-mhd-prod2-run.md) — WB, rot_potential, deep RT**
- **[Deep RT: radiative conduction + tau_R blend](radiative-conduction-deep-interior.md) — fad2f5db; `rad_kappa_src=table`, 10%**
- **[cs MHD dhj validated at rot 20](cs-dhj-diagnostics-rot20.md) — T = sp, jet weak on every grid**
- **[sp_mhd_prod dies at rot 10.3](sp-pole-bottom-radial-blowup.md) — polar-row radial field at the bottom**
- **[cs MHD dhj bottom-boundary drift fixed 248f1b77](cs-mhd-bottom-inflow.md) — cs ghost bcc as face averages**
- **[sp MHD EXCESS SOLVED: polar-row HLLD->HLLE swap](sp-polar-hlle-swap-is-the-excess.md) — `mhd/polar_hlle_rows=false` 7d87f3c5; sp_mhd_nopole testing**
- [sp MHD excess needs no field](sp-mhd-excess-is-hydro-path.md) — the bbot=0 step
- [sp MHD energy excess, earlier symptoms](sp-mhd-energy-excess.md) — superseded by the above
- **[Correlated-k design + build](correlated-k-design.md) — complete**
- **[dhj RT on the CUBED SPHERE working](dhj-cubed-sphere-port.md) — tracks sp; three bugs pre-existing**
- **[cs NON-ORTHOGONAL audit: both fixed](cs-nonorthogonal-audit.md) — e415b91a; history KE/ME ~1% wrong**
- **[cs MHD: the C2P floor corrupts u.e](cs-mhd-c2p-floor-corrupts-ue.md) — vertex order -0.04 to +1.97**
- **[cs MHD low beta: per-step operator 2nd order](cs-mhd-low-beta-divergent.md) — instability is nonlinear**
- **[cs MHD low-beta fixed on the test problems only](cs-mhd-lowbeta-fix.md) — cb1afd28; dhj dies at 0.21 rot**
- **[cs SEAM GHOST e_int retracted as a mechanism](cs-seam-ghost-eint.md) — seam degrades everything 15-19x**
- **[GENERAL-EOS AUDIT](general-eos-audit.md) — bracketed and monotone, no live bug**
- **[3D RESISTIVITY: curl B missing theta terms on sp + Cartesian](resistivity-3d-curl-missing-terms.md) — fixed; cs clean**
- **[sp POLE-EDGE CURRENT: far-side segment wrong sign, fixed](sp-pole-edge-current-sign.md) — J_r grew ~1/dtheta**
- **[sp POLE: resistive J_r divided by a zero dual area, fixed](sp-pole-edge-area-zero.md) — NaN in cycle 1**
- **[sp POLE fixed: polar_emf_diss + x3 face-state theta shift](sp-pole-fixes.md) — de667f32; L1(B) rate 0.76 to 2.1**
- **[sp_test: rigid rotation + resistive decay rates](sp-test-rigidrot-resist.md) — hydro 2.6, B 2.1**
- **[RESISTIVITY AUDIT](resistivity-audit.md) — clean; one gap behind `use_polar_average_eresist`**
- **[sp POLAR FIELD BLOW-UP fixed 1cabe85c, default ON 3147de4f](sp-polar-field-blowup.md) — GS05 corner-EMF checkerboard**
- **[cs ROTATION SOURCE BUG: cs ran the Cartesian beta-plane](cs-rotation-source-bug.md) — f75ad783**
- **[cs 1-ULP amplification: old nx=64 claim stale](cs-ulp-amplification.md) — retest alive**
- **[WHAT RECONSTRUCTION ACTUALLY RUNS](reconstruction-what-actually-runs.md) — sp ignores `reconstruct`; cs x1 stretch-aware since 4a16f07b**
- **[cs VERTEX = LIMITER CLIPPING](cs-vertex-limiter-clipping.md) — ppmx/WENO-Z cut the worst-cell force 6.5x**
- **[cs VERTEX ORACLE: halo innocent](cs-vertex-oracle-halo-innocent.md) — it is the cell balance**
- **[sp POLAR-ROW RESIDUAL fixed 7.7x](sp-polar-row-reconstruction.md) — `mesh/polar_quadratic_recon` 6410be99, off**
- **[sp GEOMETRIC-SOURCE RESIDUAL: polar rows 14x interior](sp-geometric-source-residual.md) — WB claim retracted**
- **[cs WELL-BALANCED SOURCE CACHED: 46% to 4%](cs-wb-source-cached.md) — 88064f67**
- **[GS07 corner EMF unconditional on cs since b71ef392](cs-gs07-emf-gate.md) — `bs_emf` is the way back**
- **[cs DEEP TOROIDAL SHEET closed](cs-deep-toroidal-sheet.md) — was the stretched resistive bug c5c85e3b**
- **[cs dhj production RETRY arm matrix](cs-dhj-production-retry.md) — nx=32 clean to 5.5 rot, nx=64 dies**
- **[cs MHD instability characterized](cs-mhd-instability-characterized.md) — not gnomonic-specific; two retractions**
- **[cs MHD MINIMAL REPRODUCER](cs-mhd-minimal-reproducer.md) — f3d35a96, cs_test iprob=13**
- **[cs MHD dhj blow-up OPEN, at a cube vertex](cs-mhd-dhj-blowup.md) — 10 mechanisms out; refining makes it worse**
- **[cs dhj NaN: cause found](cs-raisevel-missing-floor.md) — RaiseVel never re-applied the floors**
- **[cs_dhj_long went NaN in one rotation](cs-dhj-long-run.md) — resolved; exit 0 is not a run**
- **[ck_grav_prod: the production campaign](ck-grav-prod-run.md) — point-mass gravity, nx1 234; held**
- **[ck_mhd_b3: MHD + EOS resistivity, bbot 3 G](ck-mhd-b3-run.md) — smoke passed, not launched**
- [ck_limb: the constant-g run](ck-limb-run.md) — superseded; calibration reference
- [ck_grav_size: the sizing run](ck-grav-size-run.md) — r99 = 1.794e10, H = 22 cells
- [ck_hydro_long retired](ck-hydro-long-run.md) — sim "days" are Earth days
- [Ideal-gas + EOS x_e resistive runs](xe-resistivity-long-runs.md) — carry the stellar-heating bug
- **[dhj blow-up closed: it was the race](dhj-ck-eos-blowup.md) — no thermodynamic blow-up**
- [Exo-FMS cross-validation](exofms-cross-validation.md) — LW agrees 1.7%
- [Correlated-k shared module](correlated-k-shared-module.md) — DONE 2026-09-07 in TWO stages: correlated_k.hpp (opacity+Rosseland/conduction coupling) and two_stream_rt.hpp + atm_column.hpp (the solver); pgen 6082 -> 3648 lines; gated BITWISE identical
- [dhj conservation check](dhj-conservation-check.md) — mass +0.21%, internal E +0.56% over t=4.32e5; the hst cadence and weighting both make hst useless here
- [dhj dt limited by Alfven floor](dhj-dt-limited-by-alfven-floor.md) — MEASURED: dt is the radial Alfven CFL at r/Rp~1.28 in a cell sitting ON dfloor; raising dfloor x10 buys 2.7x
- [dhj general EOS + correlated-k blow-up](dhj-general-eos-ck-blowup.md) — OPEN: read docs/HANDOFF_dhj_ck_eos.md IN THE REPO first (9 refuted hypotheses); TWO blow-ups, only the general-EOS one is live
- [dhj high-B crash is the outer BC](dhj-highB-outer-bc.md) — FIXED bc6b5774: RKG super-time-stepping already fixes the diffusive dt; the high-B crash was the outer-x1 Maxwell term, not the Alfven speed
- [dhj ideal vs general cost](dhj-ideal-vs-general-cost.md) — general is 3.09x ideal per simulated second (10 G, done); max_eta=1e13 beats 1e14 at both fields. Its 3 G runs were CANCELLED and superseded by the viper-base set
- [dhj ideal + EOS x_e: floors, max_eta, STS](dhj-ideal-xe-floor-relaxation.md) — CLOSED, shipped as 6003c1ce: max_eta=1e13 + STS off is 1.54x faster at the same dt; STS is NOT a standing win at 1e14
- [dhj restart bugs](dhj-restart-loses-gravity-potential.md) — both FIXED (30d21859, 4cfdc329): restarts lost the gravitational potential and the cell-centered field bcc0
- [dhj viper dt discrepancy](dhj-viper-dt-discrepancy.md) — LARGELY DISSOLVED: viper bottoms at ~2.6 not ">3", and V_old (viper commit 2af153a3) reproduces it; only the t~1e6 dip is still unchecked
- [Exo-FMS correlated-k tables](exo-fms-ck-tables.md) — installed and md5-verified in the repo at data/exo_fms_ck/; ROMIO cannot write AthenaK output on /tmp
- [Perna resistivity for UHJs](resistivity-perna-uhj.md) — good to ~3x below 5000 K; the min_xe→max_eta rename that aborts the dhj runs; a tabulated x_e would be 6x FASTER
- [x_e table for ideal-gas runs](resistivity-xe-table-for-ideal.md) — DONE 5d1c3436: ohmic_resistivity=eos under eos=ideal; 7% FASTER than perna, plus two latent bugs fixed
- [sp polar bisection on orion](sp-pole-bisect-orion.md) — jobs 193408-11; the production grid RECONSTRUCTED (nx2=64 f_stretch_theta=3 verified against the 8.1 deg polar cell); dt still 7x off, so the CONTROLS decide

## Hot Jupiter: grid, EOS, atmosphere

- **[HYDROSTATIC WB scheme supports cs + point-mass gravity](wb-hydrostatic-scheme-cs.md) — 140dbf9e; cs == sp, 40-100x**
- **[Radial grid stretch](radial-grid-stretch.md) — 1e19d4e7, `use_grid_stretch_r_poly`**
- **[Radial stretch REFIT](radial-stretch-refit.md) — 3.7x to 1.34x, 1.46x dt**
- **[dt binding direction](dt-binding-direction.md) — radial stretch costs 1.76x**
- **[Point-mass gravity flag](grav-point-mass-flag.md) — c37ebe75; constant g understates H 2.7x**
- **[Stellar tide flag](stellar-tide-flag.md) — c4aa730f, validated, off**
- [Stellar tide at the domain top](stellar-tide-at-domain-top.md) — 0.5% at the limb
- **[Isobar, not shell](dhj-isobar-vs-shell.md) — measure the 1e-6 bar level on the isobar**
- [Sizing x1max from the isobar](isobar-buffer-calibration.md) — 1.76 H buffer
- [Floors for the 1e-6 bar level](dhj-floors-for-1e-6-bar.md) — pfloor 1e-3; tfloor_kelvin=50 crashes
- [Grid resolution design](dhj-grid-resolution-design.md) — cost is the timestep
- **[Upper-atmosphere mottling](upper-atm-mottling.md) — H2 dissociation front plus plumes, not noise**
- **[H2 chemistry quenches above 1e-3.6 bar](h2-chemistry-quench.md) — front not physical**
- **[sp HYDRO vs MHD at 100 rot](sp-hydro-vs-mhd-comparison.md) — field takes 31% of the zonal KE**
- **[COMPOSITION + FIELD maps in the atlas](dhj-composition-maps.md) — beta min 0.96 at tau=2/3**
- **[tau=2/3 PHOTOSPHERE diagnostic](dhj-photosphere-diagnostic.md) — a4d30d43; spans 7040 km**
- **[EOS table dump](eos-table-dump.md) — f882159f; table NaN below 71 K**
- **[EOS inversion NaN trap](eos-inversion-nan-trap.md) — naive root find 0.35 dex low in T**
- [General EOS optimization](general-eos-optimization.md) — 1.41x CPU, 1.062x GPU
- [FastChem vs the general EOS](fastchem-vs-general-eos.md) — decided against
- [UHJ band structure in the literature](uhj-band-structure-literature.md) — grid choice right

## Red giant FOFC (viper GPU)

- **[RESTARTS NOT BITWISE: user BCs read w0 before the first c2p](restart-bc-reads-w0-bug.md) — red giant fixed 3fd3a836; dhj same defect, fix in progress**
- **[RED-GIANT DEATH MECHANISM: seam-row energy floor + EOS inversion garbage T=2.6e18](red-giant-seam-floor-eos-garbage.md) — not FOFC, not chaos**
- **[RG_fofc SURVIVED on GPU AND on orion CPU](rg-fofc-gpu-survived.md) — FOFC replaces the vceil; floors are vertex corner ghosts**

## Solar / stellar convection

- **[run/sun viz + WAVES NOT SHOCKS](sun-convection-viz.md) — gravity waves, 87% solenoidal; data on orion**
- [solar_convection binary provenance](solar-convection-binary-provenance.md) — the pinned ideal192/table192 binaries cannot be rebuilt from git
- [solar_convection with the general EOS](solar-convection-general-eos.md) — OPEN: the atmospheric runaway survives the corrected sponge; sub-photospheric comparison is done and solid. RESUME HERE for this thread.
- [solar_convection restart gravity bug](solar-convection-restart-gravity-bug.md) — FIXED/pushed ef9561e2: restarts ran with zero gravity, dt collapsed to 1e-135
- [solar_convection scaling wall](solar-convection-scaling-wall.md) — MEASURED: saturates at ~36 ranks; run it on ONE node, multinode buys 6%
- [Solar convection test](solar-convection-test.md) — tuning solar_convection pgen (two-stream RT) for solar-like convection + stable atmosphere
- [Sponge inert at 192^3](sponge-inert-at-192.md) — RESOLVED: the sponge read pcoord->x1v, a 1x1 placeholder on Cartesian meshes; the 64^2 sponge results are invalid

## Bugs and performance

- **[FOFC on GPU VERIFIED](fofc-gpu-verified.md) — 48/48 HIP Release; Debug abort was a pre-existing dt-diagnostic over-read, fixed b18f3205**
- **[ORION MERGE round-off gated caad9247](orion-merge-roundoff-8da093f5.md) — `floors_legacy`; general-EOS c2p was the hunk**
- **[All-Mach solvers lhllc/ausmpup: algebra matches the papers](lowmach-solver-audit.md) — no CFL~M restriction**
- **[REFLECTING WALLS leaked mass under a blast: fixed](reflect-wall-mass-leak.md) — mirror the wall state**
- **[GPU nondeterminism fixed 6e600f12](dhj-run-to-run-nondeterminism.md) — missing team_barrier**
- [Polar MPI host-mirror bug fixed](polar-mpi-host-mirror-bug.md) — 3882e37f
- [Restart bug with tfloor_kelvin fixed](athenak-restart-tfloor-bug.md) — ba2f0943
- [dhj ideal input out-of-bounds fixed](dhj-ideal-input-oob.md) — 0d1f6f6a; use a Debug build
- [MHD flux stale x1 limits](mhd-fluxes-stale-x1-limits.md) — 3429f59f; latent
- [RT chain-parallel split kernel](rt-chain-parallel-split.md) — 38311a8a, live
- **[MAGNETIC ENERGY is not in the dumps](mhd-energy-not-from-dumps.md) — bcc is the mean of the fields**
- [rt_chain_ck SWEEP SPLIT measured 1.13x, shelved](rt-ck-sweep-split-measured.md) — patch in bench/prof/rt_new
- [Overall GPU profile of a dhj run](dhj-overall-gpu-profile.md) — RT 35.7%, fluxes 27.8%
- [Meshblock decomposition on one GPU](meshblock-decomposition-gpu.md) — 32 blocks beats 2 by 1.4x
- [nx1 ceiling from LDS](nx1-ceiling-lds.md) — MHD capped at nx1=264; hydro escape `hydro/scratch_level=1` 955c39df
- [par_for_outer team size](par-for-outer-team-size.md) — ~3% of wall, not worth it
- RT kernel history, superseded by rt-chain-parallel-split: [interleaving](inflight-rt-kernel-optimization.md), [band scaling](inflight-rt-band-scaling.md), [occupancy](rt-kernel-occupancy-limit.md)
- [IN FLIGHT: nx2=64 + GPU timings](inflight-nx64-timings.md) — paused 2026-08-16
- [AMR null-tree bug](amr-broken-on-branch.md) — fixed; the cubed-sphere refactor left `Mesh::ptree` unconstructed, and the same half-migration may bite elsewhere
- [Pre-existing branch breakage](branch-preexisting-breakage.md) — non-MPI build fails in bfield_bcs.cpp; hydro/mhd_linwave test scripts pass a `vflow` param that no longer exists
- [Event log MPI deadlock](eventlog-mpi-deadlock.md) — FIXED 27da6380: file_type=log hangs ANY clean MPI run (early return skipped the schedule advance -> Allreduce every cycle); likely upstream too
- [FOFC 1D segfault](fofc-1d-segfault.md) — pre-existing crash on the branch: FOFC + 1D hydro + PLM, not caused by the EOS work
- [FOFC COMPATIBILITY COMPLETE 09-11 incl. MHD on cs + static/dynamic WB tests (new pgen wb_atm; dynamic WB: fofc raises the Mach ceiling 39 -> >80) + MPI/SMR tests; COMMITTED b3345f55 + merged 7cee2eb7 on 09-12, not pushed; 49 fofc tests pass, bit-identical off; left: red-giant fofc run, GPU check, wiki, commit split from red-giant hunks](fofc-compatibility-plan.md)
- [hst uses Cartesian cell volume](hst-cartesian-volume-on-spherical.md) — .hst mass/tot-E are NOT physical integrals on the spherical dhj mesh
- [pcoord arrays unallocated in Cartesian](pgen-cartesian-coordinate-arrays.md) — which geometry arrays are placeholders on a Cartesian mesh, and the solar/cooling_convection audit for unguarded uses
- [RT source is semi-implicit](rt-source-semi-implicit.md) — FIXED 048dff30: exponential relaxation at rate 4E/e; cured the red giant dt collapse; the wall control is what ruled out the boundary
- [RT clip warning is rank-local](rt-srclim-warn-rank-local.md) — BUG, unfixed: the clip count is rank 0's share, and the warning is silent when rank 0 never clips
- [RT source in a transparent cell: RESOLVED](rt-transparent-cell-cancellation.md) — direct source in all 3 kernels + the REAL cause: the opacity lookup CLAMPED rho to the table edge (1e-14) and returned molecular opacity; table now extended to 1e-24 at the floor; background holds 411 K
- [two_stream_rt has THREE sweep kernels (grey/ck/generic)](two-stream-rt-three-kernels-trap.md) — a sweep edit must hit all three; red_giant uses rt_chain_grey; verify with rt_apply_debug. Also: <hydro>/wb_rmax added+verified
- [VIPER floors_legacy GATE 09-11 (caad9247): orion merge changed GPU answers 1 ulp via the general-EOS c2p; legacy kernels in separate .cpp under floors_legacy (all floor switches off); MERGED locally 7cee2eb7 on 09-12](viper-floors-legacy-gate.md)
- [VIPER MERGE LANDED 09-11 08:05: a10e367d on top of 9 topical commits; gates bit-identical both sides; not pushed](viper-merge-trial-2026-09-11.md)

## Working practice

- **[VALIDATE THE INSTRUMENT](validate-the-instrument.md) — "nothing" usually means "I did not look"**
- **[Measure impact before claiming it](measure-impact-before-claiming.md) — report the defect, not consequences**
- [Localise by dilution](localise-by-dilution.md) — vary domain size at fixed dx
- [Don't propose confounded tests](confounded-tests-rejected.md) — user vetoes non-discriminating tests
- **[apudev: ONE-TIME TESTS ONLY](apudev-one-time-tests-only.md) — chained/production jobs go to `apu`**
- [HIP DualView sync idiom + hipcc contraction](hip-dualview-sync-idiom.md) — modify_host/sync_device only; kernel splits lose GPU bitwise
- [Job submission permitted](job-submission-permitted.md) — from the main session
- [Use fork, not origin](use-fork-not-origin.md) — push to jing-ze-ma/athenak
- **[Restart output dt trap](restart-output-dt-trap.md) — last_time from the restart; dt is sim time**
- [Never write in run/](never-write-in-run-dir.md) — read-only
- [Viper HIP build recipe](viper-hip-build-recipe.md) — cmake flags for MI300A
- [Reusable CPU baseline binary](bench-baseline-worktree.md) — bench/base_wt
- [Fork git setup](athenak-fork-git-setup.md) — origin = jing-ze-ma/athenak fork over SSH (all git ops); upstream = IAS-Astrophysics; the user's terminal cannot copy OUT
- [Freya build procedure](freya-build-procedure.md) — how to build AthenaK on Freya (MPI+OpenMP, SPR); clean only build/ contents
- [Freya/Orion job submission](freya-job-submission.md) — p.shared vs p.exclusive, keep --cpus-per-task with OMP, and an orion node has 112 PHYSICAL cores (the 224 SLURM reports are hyperthreads)
- [GPFS quota wall](gpfs-quota-wall.md) — df lies; a per-user quota SILENTLY truncates AthenaK dumps without any error
- [Orion build traps: OpenMP-off binary + module load](orion-build-openmp-and-module-traps.md) — a Kokkos_ENABLE_OPENMP=OFF build is SILENT and 3.7x slow at 96x7; `module load` does not survive between tool calls, so make falls back to gcc 7.5
- [run/ is untouchable](run-directory-untouchable.md) — ~40GB of simulation data; never add, clean, or modify it
- [Scratchpad invisible to compute nodes](scratchpad-not-visible-to-compute-nodes.md) — /tmp is node-local; stage SLURM jobs on /orion
- [Test output location](test-output-location.md) — all run output/analysis goes in /orion/ptmp/jinma/Athenak/ (separate GPFS, no quota); never /orion/u or /tmp
- **[VIPER HIP CONVENTIONS (standing, 09-11): DualView modify_device/sync_host only; kernel splits change GPU round-off; CPU bitwise != GPU bitwise; build on viper before calling GPU code done](viper-hip-code-conventions.md)**
- [Work pace preference](work-pace-preference.md) — move faster / less deliberation on routine tasks

## General EOS

- [EOS electron regression test](eos-electron-regression-test.md) — DONE: mhd_eos_electrons.py, plus how to dry-run a test without letting run_tests.py delete build/
- [EOS picked by run name](eos-selection-by-run-name-bug.md) — analysis scripts guessed the EOS from the directory name and silently used the wrong one; fixed, plus the audit that found the scope
- [EOS x_e inert in the resistivity](eos-xe-resistivity-capped.md) — RESOLVED: not the EOS at all; two out-of-bounds reads of eta_b, one of which made ohmic_resistivity=constant a silent no-op
- [eostest input inventory](eostest-input-inventory.md) — what the ideal/general comparison inputs cover and how to run them; .bin output is single precision
- [General EOS: MHD+polar difference](general-eos-mhd-polar-bug.md) — resolved, not a bug: the polar boundary forces HLLD→HLLE at the pole, and HLLE is the solver whose Roe average has no general-EOS analogue
- [General EOS project](general-eos-project.md) — non-ideal EOS for Newtonian hydro/MHD; branch `general-eos`; ALL stages including the analytic EOS are done and verified. START HERE.
- [General EOS: Stage 3 loose ends](general-eos-stage3-loose-ends.md) — the last three TODO(stage3)s closed; why the WB background hands back a TEMPERATURE and not a pressure
- [General EOS: Stage 3 perf cleanup](general-eos-stage3-perf-cleanup.md) — one T solve per cell instead of ~11-30; fully verified and committed as 7d1fed3f on 2026-08-12
- [General EOS: Stage 3 tabulated EOS](general-eos-stage3-table.md) — the analytic EOS itself: H2 + Saha + radiation via a (log rho, log T) table; done and verified, commit 03febafa
- [General EOS Stage 4: (rho,e) table](general-eos-stage4-rho-e-table.md) — SCOPED, NOT STARTED; re-baselined 2026-08-17 to ~1.5x (3.07x -> ~2.0x), proxy ceiling measured twice; task 1 is the (rho,e) domain shape
- [General EOS table cost](general-eos-table-cost.md) — RE-MEASURED 2026-08-17: 3.07x per cycle / 2.32x per simulated second, NOT 4.4x; profile (49% libm) still stands. Restrict any such comparison to the hydro-limited window
- [General EOS: table linear-wave tests](general-eos-table-linear-wave.md) — how general_eos=table is regression tested; the tabulated EOS cannot be called from host code

## Red giant envelope

- [1.5 R grid STARVED the deep: 0.92 cells/H_p at the wall -> dies at 3.2e5](red-giant-15R-grid-starved-deep.md) — ambient medium was fine; refit with >=5 cells/H_p in the star and more nx1; measure on the actual grid, not the fitter
- [Red giant: the ANALYSIS used the wrong opacity](red-giant-analysis-opacity-trap.md) — eoslib.get_kapr is solar_convection's analytic formula, 30-44x too high in the interior; the RUN is clean; retracts the "wall cannot shed the flux" story
- [OPAC CLAMP UNNECESSARY OPEN-TOP TOO (I3, 09-10 23:30): I2 vs I3 indistinguishable to 6.2e5, top 50-90 K warmer without it, fewer cold cells; open: 31.6 K floor cells in both after 6e5, L_out/L ~2.7 after the cap burst](red-giant-clamp-off-open-top-ok.md)
- [COMMON DEATH 5.71e5 (V9f/V10/V9_dense): EXPLICIT radiative conduction runaway in an underdense supersonic downdraft at 1.058 R on a cube EDGE, inside the star; corona/top innocent; next test = sub-cycled/implicit radial conduction from rst 5e5](red-giant-common-dt-collapse-571e5.md) — START HERE 09-10 18:00
- **[red giant dt COLLAPSE is RADIATIVE CONDUCTION, still OPEN](red-giant-conduction-dt-collapse.md) — RESUME HERE: the RT fix moved it 2.9e4 -> 1.86e5 but did not remove it; MEASURED: the limiting cell is IN THE BLEND WINDOW at the photosphere (tau 1-3), both times; tau blend 10/100 NaNs; STS for conduction is the promising cure**
- [CORONA V8 DIED 2.96e5 (09-10 07:50): the LAST ACTIVE CELL below the kappa cutoff cools to the 31.6 K floor in 2e4 s, the corona accretes onto it and evaporates; the RT "Newton rescue = 99.9% drop" message precedes every open-top/corona failure (V7, V8) -> suspect; V9 (2e4 K) follows R7's atmosphere-cooling path](red-giant-corona-join-cell-cooling.md)
- [Red giant DEEP motion = RCB pile-up transient](red-giant-deep-onset-is-rcb-pileup.md) — L switched on with no convection piles heat above the RCB; sup grows linearly from t=0, ignites at 7e5 at 200x MLT driving, overshoots 10x; surface convects first and correctly; remedy under test: mlt_alpha=3
- [DFLOOR-ENERGY DIAGNOSIS WRONG 09-11 11:10: the 1e16 K dt-dip cells are VACUUM cells (rho_old<=0, v=-0) whose energy GROWS at v=0; dfloor_keep_temperature (e*=fv) vacuous and UNCOMMITTED; proposal: vacuum reset + instrument the energy source; user decides](red-giant-dfloor-energy-diagnosis-wrong.md)
- [Red giant: dt collapse SOLVED](red-giant-dt-collapse-solved.md) — the t=1.9e5 collapse was the EXPLICIT radiative diffusion at a cooling photospheric cell; a DEEPER tau blend (10/100) hands it to the semi-implicit two-stream and dt stays 30.7 s
- [DUST OPACITY kills lid-free runs: PROVEN, opac_tmin=2500 survives; explicit RT is a negative control](red-giant-dust-opacity-kills-lidfree.md) — START HERE 09-09 14:40: T6 extended to 3e6 (194594); prod12 candidate = lid-free+open+grains off; prod11 lidded control at 4.1e6
- [Red giant envelope project](red-giant-envelope-project.md) — 2026-09-07: pgen builds and RUNS on the cubed sphere with the CORRELATED-K two-stream (problem/rt_ck); stellar opacities merged (OPLIB+AESOPUS); grey 1-D validation passed; COMMITTED+PUSHED (35ea57ba, 1d99e22e, 361be207)
- [PHOTOSPHERIC RUNAWAY: photosphere IS inside (1.077 R, 20 cells under the top; the 'outside' claim was a misread); 1-cell ionization front; every blend/two-stream-semi variant dies ~1.84e6, EXPLICIT two-stream sails through and dies at the top cell instead; D1_catch budget pending](red-giant-explicit-conduction-explosion.md) — START HERE 09-10 evening
- [FLOOR DEATH 1.43114e6 CLOSED 09-10 02:40: efloor_from_ekin alone cures it (FL2) AND wb_rmax=3.3e12 alone cures it (FL4, floor energy still created but harmless); ship BOTH; ported to main tree + random seed added; smoke tests V5/V6](red-giant-floor-energy-creation-fix.md)
- [Red giant: the 0.34 L deficit is a SPIN-UP transient, not a bug](red-giant-flux-deficit-is-spinup.md) — CHECKED 09-08 23:20 at t=5.5e6: the convection front IS advancing outward ~25 cells/1e6 s, deep runaway saturated; chain 194259+194262-5 healthy, just needs wall time
- [GREY OPAQUE-LID BUG: grav=0 on the grey path -> dtau_top = inf in every grey run incl. prod11](red-giant-grey-opaque-lid-bug.md) — fix doubles L_out at t=0 but the thin atmosphere cools 3364->2200 K; being isolated on the plain star (R8); corona itself was inert and fine
- [Red giant: the GREY two-stream, rebuilt](red-giant-grey-two-stream.md) — bd156d77 `problem/rt_grey`: one band, the conduction module's opacity; removes the thermal runaway and the runs live past every ck death point. The drain/outer-BC follow-ups are CLOSED by rt_top_re
- [V8f 6e5 K CORONA DROPPED (09-10 18:20): it is a Parker wind (T 3x above the critical T at the join, H = 1.9 R); the 14-scale-height atmosphere drains, over-pressured corona shocks down, dfloor cells with no velocity ceiling spike -> hydro dt 0.95 s at 3.2115e5; enable eos_vceil in any open-top run](red-giant-hot-corona-is-a-parker-wind.md)
- [I4 ACTIVE 400 K MEDIUM: passed V10 death, top never floored, then DIED 6.31e5 of a hydro collapse at an ejected-shell front (all open-top runs launch a shell at 5.7e5; the pressure-supported corona lets I2/I3 ride it out); next = DENSITY-GATED opacity (corona inert, no radius cutoff)](red-giant-i4-active-medium.md)
- [I5 DENSITY GATE DONE (09-11 01:40): ran to 9e5, cutoff artefact CONFIRMED (L_out 1.2-3.3 vs 2.6-5.3 with the cutoff, L_cut 0.06 vs 1.36); 0.35 L at restart = clamp off, not the gate; corona inert; NEW: photospheric cells reach 1e16-1e17 K under flat dt -> runaway source (hydro or RT) still open(red-giant-i5-density-gate.md)
- [IMPLICIT RADIAL DIFFUSION + ANGULAR CAP DONE 09-10 23:00 (hydro/rad_implicit_x1, rad_cap_ang=0.5): I2 ran 5e5 -> 9e5 with dt flat, tot-E flat, capping only 5.69-5.98e5 then silent; +21% cost; binary build_impl/athena_pin2](red-giant-implicit-radial-diffusion.md)
- [MLT closure lessons: per-column MLT flickers deep and cannot hand over at the top](red-giant-mlt-closure-lessons.md) — R1/R1b/R1c all failed at the photosphere; fix = shell-averaged 1D MLT with f = min(F_MLT, L - F_rad,diff) (problem/mlt_mean)
- [MOLECULAR-OPACITY KNEE RUNAWAY (09-10 00:45): open-top cells cool radiatively through the kappa_R knee at 3200-2500 K (emission grows as T falls), park at 200-500 K, become the Riemann precursor; lid keeps the atmosphere 1000 K hotter; TEST opac_tmin~3200](red-giant-molecular-opacity-knee-runaway.md) — START HERE
- [T7 (grains off + WB fix) STILL died at 1.04e6 with no precursor: the 100-cycle NaN check hides the origin (dt ignores NaN); T8 194637 reruns with the check every cycle — OPEN](red-giant-nan-check-hides-origin.md) — START HERE 09-09 17:00
- [Red giant: open outer BC sealed the star, FIXED bd2b974d](red-giant-open-outer-seals-the-star.md) — `problem/rt_top_re` makes the column above the domain RADIATE and the 6650 K seal is gone; its "the defect is the handover" conclusion is SUPERSEDED by the spin-up note above
- [3 L EMERGENT FLUX = OPACITY-CUTOFF ARTEFACT (09-11 01:00): the inflating surface crosses rad_kappa_rmax where kappa=0 and radiates 9000 K gas as a blackbody; lidded top is a uniform 9200 K/1e-8 at 1.3-1.5 L; "+34 L dE/dt" is float32 dump noise; I4 (no cutoff) is the test](red-giant-open-top-kappa-cutoff-artefact.md)
- [OPEN-TOP prod12 RECIPE FAILS FROM SCRATCH (V7, 09-10 06:50): cold collapse of the thin atmosphere at 5.8-6.0e5, EARLIER than V4; opac_tmin=3200 thins the open atmosphere 19x; corona V8 is the remaining open-top option; lidded B12 is healthy](red-giant-open-top-prod12-recipe-fails.md)
- [prod11: ROOT CAUSE FIXED and the LIDDED REVIVAL GATE PASSED (B12, 09-10 06:55: fixed binary + opac_tmin=3200 from rg.00124 to 3.45e7, dt flat, 0 vertex cold cells; cold population grows slowly, recycles, no dt response); revive from B12_fix_tmin3200/rst](red-giant-prod11-died-3e7.md)
- [PRODUCTION prod11 LAUNCHED 06:40: 1.1 R, all fixes, chain 194515-7](red-giant-prod11-launched.md) — START HERE next session: check dt, L_out, deep quiet, surface onset; old chain cancelled; option A (hot corona) proposed, not approved
- [Quadrature 0.904x, RCB dipole = seed?, pin binaries](red-giant-quadrature-seed-binary.md) — ck_nquad=1 is 10% below diffusion (use 2); the RCB dipole is closure-blind (R0==MLT to 4 digits) -> testing vpert_rmin; never rebuild under a running job
- [R3 RESULT 09-11 07:20: dfloor_keep_velocity correct but NOT the loop (R3_bface died like R2 at 5.8667e5); cells at 1e5 x escape speed, KE-cancellation energy creation + NaN T in RT from negative ei -> NEXT: Newtonian hydro/vceil + RT NaN guard, R4](red-giant-r3-dfloor-keep-velocity-not-the-loop.md)
- [RESTART x730 radial-KE kick: FOUND+FIXED (WB cache zero for ncycle%N cycles after restart; affects ALL well-balanced restarts, dhj too), COMMITTED 5c0b98e4, pushed; viper must pull](red-giant-restart-radial-ke-injection.md) — START HERE 09-09 16:20: prod11 chain 194624/5 with wb_cache_every=1; T7_nodust_fix 194627 is the prod12 gate
- [RG_fofc LAUNCHED 09-11 15:30: red giant from scratch nghost=3, fofc on, NO vceil, floor fixes on, chain 196662-4; deaths to beat 5.7e5/5.87e5/6.13e5](red-giant-rg-fofc-run.md)
- [RT SWITCHES DONE 09-11 06:30: rt_bface + 4 RT defaults now OLD by default, rt_top_clamp new; ck dhj bit-identical to HEAD; red giant inputs need SIX opt-ins from pin10 on (repo input has them, ptmp run inputs do not)](red-giant-rt-bface-switch-and-head-regression.md)
- [RT NEIGHBOUR-PLANCK BUG FIXED+GATED 09-10 09:00: BFace kappa-weighting at 16 sites; smooth star bit-identical; B13 (lidded, NO clamp, NO wb_rmax) holds dt flat past the old collapse with cold cells 0-27 vs 156 -> the cold-collapse family WAS this bug; V8f join cell holds 3300 K](red-giant-rt-neighbour-planck-bug.md)
- [RUNAWAY SOURCE (09-11 03:30): RT read STALE w0 -> fixed (rt_use_cons, floor energy creation gone, R2); loop reopens via the DENSITY FLOOR keeping momentum (v 4e11 in a floored cell) -> dfloor_keep_velocity in progress (R3/I7); I6 DIED the same way at 5.877e5 (v 5e16 in a dfloor cell)(red-giant-runaway-source-rt-stale-w0.md)
- [PROVEN: the vpert seed across the RCB causes the deep transient](red-giant-seed-across-rcb-proven.md) — vpert_rmin=1.6e12 flattens the dipole with/without MLT; nquad=2 killed both no-MLT runs at 1.44e6 (bisecting); MLT gap = chi limiter
- [Red giant session 2026-09-08](red-giant-session-2026-09-08.md) — the 10 commits, the working production config, and the orion SLURM core-vs-hyperthread trap; its jobs are all finished/cancelled
- [SPONGE KE BUG on the cubed sphere: real, fixed (887c241e), but NOT the killer — T11 died the same way at 8.72e5](red-giant-sponge-cs-kinetic-energy-bug.md) — START HERE 09-09 17:50: the lid-free deaths at ~1e6 were the sponge draining the cube-vertex cells; not RT, not the ghost, not dust (that was the EARLIER deaths)
- [STALE-PRIMITIVE AUDIT (09-11 02:40): 6 more sites of the RT bug class; LIVE: arm rt_use_cons, RedGiantBC ghosts seeded from stale w0 (at the dying faces), rg_wallflux unguarded ei; LATENT: rg_relax, open-inner BC (incl. orthogonal KE on cs), MLT caps, perturb-reconstruction ordering](red-giant-stale-primitive-audit.md)
- [TOP-CELL NEGATIVE DENSITY 09-11 12:40 = root of the RT -nan (rescue was an accidental energy sink) and likely of the 1e16 K vacuum cells; RT guard in tree UNCOMMITTED; with it the implicit conduction fails at 6.1265e5 -> needs FOFC step 3 (cs) or a vacuum/top-cell treatment; user decides](red-giant-top-cell-negative-density.md)
- [Red giant: the TOP cools away and kills the run](red-giant-top-cooling-runaway.md) — historical: all three jobs died of a HYDRO dt collapse; the optically thin top cell falls 3364 -> 1100 K and the cooling GROWS as it cools; the emergent-flux decay is this, not the interior
- [two-stream CANNOT own the interior (09-10 19:00): deposition throttled to src*t_rad deep; core L is a conduction wall term; killer cell is tau_cell~1 in a tau>100 column -> diffusion operator invalid there; cure = implicit radial diffusion (one block radially, no MPI) or local-tau blend](red-giant-two-stream-cannot-own-interior.md)
- [VCEIL WORKS 09-11 07:00: hydro/vceil=5e7 + RT NaN guard; R4 passed the 5.87e5 death; the "hang" was runaway_scan's rank-local line budget skipping an Allreduce (FIXED, pin11, ghost-cell counter fix too); I8_vceil 196221 clean at 6.26e5 -> 9e5; next physics = denser corona](red-giant-vceil-r4-passes-then-hangs.md)
- [VERTEX CHIMNEY dissected 09-09 22:30: sponge+open top lethal, sponge-off survives, wall caps; T15 = prod12 gate sponge-off to 3e6](red-giant-vertex-chimney.md) — START HERE for prod12
- [HYDRO TRIGGER CAUGHT (09-09 23:30) + CENSUS 00:15: a COLD-COLLAPSED cell in the thin layer above the 1-cell front (any angle; vertex only first) + vacuum-Riemann energy creation; radiation only amplifies; open-top atmosphere is grid-scale noisy (x3 adjacent-cell rho) and grid-locked](red-giant-vertex-floor-hydro-trigger.md) — START HERE
- [CHART-FREE SEED: vertex asymmetry GONE (V4 copies bit-identical) but V4 still DIED at 9.45e5 by an INTERIOR cold-collapse -> conduction-dt (no floor event); seed buys 2x time only; random component added (V6)](red-giant-vpert-seed-chart-imprint.md)
- [Red giant: the ambient medium is killed by the WELL-BALANCED scheme](red-giant-wb-kills-ambient-medium.md) — WB on: NaN in 2 cycles, ~400 g spurious force; WB off: background free-falls at exactly g and survives. RT/BC/opacity fixes were all red herrings for THIS failure. Fix: <hydro>/wb_rmax
- [WB POLYTROPIC WALK OVERFLOW = ROOT CAUSE of the prod11 NaN AND the "WB kick at a cold cell": FIXED in the working tree 09-10 (4 hunks), bit-identical on healthy data, B11 passes the death with WB ON and the sawtooth is gone; binary build_guard](red-giant-wb-polytropic-walk-overflow.md)

## Session states (orion)

- [Session state 2026-08-16](session-state-2026-08-16.md) — superseded: four unpushed commits on general-eos; the open thread was the solar atmospheric runaway
- [Session state 2026-08-17](session-state-2026-08-17.md) — previous entry point (origin/general-eos = d1c289dd)
- [Session state 2026-08-18](session-state-2026-08-18.md) — superseded entry point: no code changed today; all cost questions answered; the ONE open thread is the viper dt discrepancy
- [Session state 2026-08-24](session-state-2026-08-24.md) — superseded: the ckrepro reproducers ran clean; hst is unusable here (cadence + Cartesian volume)
- [Session state 2026-08-25](session-state-2026-08-25.md) — superseded: job 190939 was left running that day; the viper doc's blow-up does NOT reproduce on orion CPU; input pushed as cd2e8815
- [Session state 2026-09-04](session-state-2026-09-04.md) — START HERE: artifact PUBLISHED with the time player + ray-marched convection plumes; nothing pending
- [Session state 2026-09-07 (orion)](session-state-2026-09-07-orion.md) — START HERE: pulled viper's handover (5e9b1273); read docs/handover/HANDOVER-2026-09-07.md; steps 1-2, the determinism test and the EOS table all DONE; left: GPU build with d3d74f2b (needs viper)
- [Session state 2026-09-09 evening](session-state-2026-09-09-evening.md) — superseded by the NIGHT note; keeps the T11/T10 details and the earlier agent transcripts
- [Session state 2026-09-09 night](session-state-2026-09-09-night.md) — superseded by the 09-10 EARLY note: prod12 gate FAILED (see the explosion note), nothing running except prod11 chain 194625; the user must choose: bigger domain, stay lidded, or implicit diffusion
- [Session state 2026-09-10 05:30](session-state-2026-09-10-0300.md) — superseded by the 09:10 note; keeps the FL/V4/V5/V6/B-series detail
- [Session state 2026-09-10 09:10](session-state-2026-09-10-0910.md) — superseded by the 22:00 note: WB-walk overflow + RT neighbour-Planck bugs FIXED and gated (build_guard); lidded run needs no clamp (B13); V9f/V10 open-top tests running to ~12:40; commit plan; user decides prod11 revival / open-top config
- [Session state 2026-09-10 22:00](session-state-2026-09-10-2200.md) — superseded by the 09-11 03:40 note: implicit radial diffusion + angular cap gated; I2 alive past both deaths; I3 (clamp off) and I4 (400 K active medium, from scratch) launched with watcher agents; tree uncommitted (~29 files)
- [Session state 2026-09-10 early](session-state-2026-09-10-early.md) — superseded by the 03:00 note; keeps the FL/V4/catch job list and the agent transcript paths
- [Session state 2026-09-11 03:00](session-state-2026-09-11-0340.md) — superseded by the 12:50 note: implicit radial diffusion + angular cap gated, V8f dropped (Parker wind), I4 active medium
- [SESSION STATE 2026-09-11 12:50](session-state-2026-09-11-0800.md) — START HERE: a10e367d pushed; FOFC 0+1 committed c0163568 (not pushed); I8 done 9e5; CHECK C1 first; uncommitted RT guard (keep) + vacuous dfloor_keep_temperature (drop); decision: top-cell negative density -> FOFC step 3
- [SESSION STATE 2026-09-12 01:00](session-state-2026-09-12-0100.md) — START HERE: FOFC committed locally + viper merged 7cee2eb7 (NOT pushed); gates pass except dhj ck 1-ulp from the RT guard e1db81d8; RG_fofc running, no events at 1.5e5

---

**CURRENT STATE, 2026-09-12 (~00:30), STOPPED on VIPER.** HEAD = fork = fb69c554 on `polar-average-perf`
(bc972cc8 ME fix + MHD implicit conduction; 1bd31298/27ca5b13 cs seam packs bitwise; 77de5618 radimpx1 split, 1 ulp on HIP;
fb69c554 HIP conventions note for orion). cs prod restart 62.0 -> 40.5 ms/cycle; cs grid gap vs sp CLOSED (37.2 vs 38.5).
RUNNING on apu: cs_mhd_prod3 (11611223 + links 11611229/30/31, rot 6 clean, binary swapped to 27ca5b13 for links 2+) and
sp_mhd_prod3 (11569325, rot ~120). s01 AND s05 switch arms all reached 8 rot incl. controls: switches harmless, rescue
UNPROVEN (controls no longer die). NEXT: read cs_mhd_prod3 (drained columns, rt_eiclamp, 1-ME vs sp) at rot ~20-42;
rt_chain_ck is 34% of cs GPU time = next speed lever; the handover list (dt_min 1e-3, density floor vs WB). Details in
[[inflight-2026-09-09-viper]] bottom. Rules: save tokens, delegate to Opus, apudev one-time only, never write in run/,
sbatch from the main session is permitted.

**ORION STATE 2026-09-12 01:00:** the orion snapshot's own bottom paragraph is still dated 09-07 and reads:
"**CURRENT STATE, 2026-09-07 (~09:45), HANDOVER TO ORION.** HEAD 7c652768 on `polar-average-perf`, PUSHED
to the fork (jing-ze-ma/athenak). Viper is in MAINTENANCE 09-07 12:00 -> 09-12 12:00; the next session runs
on ORION, which cannot see viper's bench/ or scratch. Everything needed is in git:
docs/handover/HANDOVER-2026-09-07.md (READ FIRST), docs/handover/scripts/ (the analysis scripts),
docs/handover/claude-memory-2026-09-07/ (a copy of this memory directory as of the handover)."
The live orion entry point is [[session-state-2026-09-12-0100]]: FOFC committed locally + viper merged
7cee2eb7 (NOT pushed); gates pass except dhj ck 1-ulp from the RT guard e1db81d8; RG_fofc running, no
events at 1.5e5. Standing orion rules: [[delegate-heavy-work-to-opus]], [[save-tokens-everywhere]],
[[viper-hip-code-conventions]], [[run-directory-untouchable]], [[test-output-location]].
