# Viper session 2026-09-29 03:15 to 09-30 06:40 CEST (supersedes SESSION-2026-09-28-viper.md) — START HERE on viper

Read this, then `MEMORY.md`. Memory snapshot of this handover: **`docs/handover/claude-memory-2026-09-30/`**
(copy it into the local memory dir on a machine whose memory is older). `rt-integration` tip at writing:
**4067ce18** (pushed to fork). Always `git fetch fork` before pushing; never force.

## 0. Read this first (newest first)

- **09-30 05:14** viper runs **w3xk** = WASP-121b 3x with ONLY `ck_impl_dtmax 0.25` + `ck_impl_maxit 24` (tol 1e-8,
  no rsec), fresh start, jobs **12034178 (+12034179 afterany)**, pending on apu (est. start 09-30 08:01).
  Compare with the Caltech 3x (3x keys) on dt decline, night-side downflows, top-cell ck cycle.
- **09-30 ~04:30** Caltech 3x dt fell 12 -> 7.5 s (rot 18-30): radial CFL from supersonic night-side downflows.
  Caltech A/B from its rot-36 rst (3649979 3x keys / 3649980 10x keys). User: viper repeats it only if the Caltech
  result is ambiguous (TASK-2026-09-30-w121-3x-ck-solver-ab.md; snapshot staged at Caltech `to_viper_0930/`).
- **09-30 03:15** DeltaAI merges gated on viper (bitwise, 15/15 CPU tests) + FOFC inline-lambda fix merged b6a6eff3.
- **09-30 ~01:25** viper MHD max_eta arms **held** (see 1). The 3x runs on Caltech only (viper copy cancelled 03:06).
- **09-30 ~00:50** opacity-Newton guard merged (23ec1fb8): fixes negative implicit-row diagonals from the Newton
  default (the He presn solver breakdown). He presn 128x128 wedge resubmitted with the guard (12030056-58).

## 1. Running / queued (checked `squeue -u jinma` 09-30 06:38 CEST)

| where | job(s) | what | dir | state / ETA |
|---|---|---|---|---|
| viper apu | 12034178 (+12034179 afterany) | w3xk: 3x with dtmax 0.25 + maxit 24, TROT 300, binary athena.gpu.w3xk = cba4e797 (src = b6a6eff3) | /viper/ptmp2/jinma/w121prod_0929/w3xk | PD, est. start 09-30 08:01 |
| viper apu | 12026299 (32 nodes), 12028087 (16 nodes), 30 min each | C256 x nx1 256 benchmarks from the remapped 1x rot-300 rst | /viper/ptmp2/jinma/w121_c256_0929 | PD, est. start 07:34 / 07:51 |
| viper apu | 12030056 -> 12030057 -> 12030058 (afterany), 6 h each | He presn wedge 128x128 (ad3d_128), guarded Newton, binary athena_he_gpu72_0b8b6c0d | /viper/ptmp2/jinma/hepresn3d_0929/M1/ad3d_128 | PD, est. start 07:56 |
| viper apu | 12028523 b10_e14_f13, 12028524 b10_e13_f16 (24 h), 12028525 b3_e13_f16 (6 h) | w121 MHD max_eta long arms, binary athena_dhj_gpu72_ae767d20 | /viper/ptmp2/jinma/w121_mhd_0929/viper_arms/<arm> | **HELD** (user 09-30 01:25) pending the DeltaAI short arms; release with `scontrol release 12028523 12028524 12028525`. 12028522 (b3_e14_f13) cancelled by the user. |
| Caltech | 3640302 | WASP-121b 3x production (793e03c3, 3x ck keys), started 09-29 17:46 PDT, dt ~12 -> 7.5 s | /resnick/groups/carnegie_poc/jingze/w121prod_0930/ | ETA rot 300 was 09-30 06:00 PDT (at dt 12 s; later with the decline). Not verifiable from viper. |
| Caltech | 3649979 / 3649980 | 3x ck solver A/B from rot-36 rst (A 3x keys, B 10x keys), nlim 417000, verbose | /resnick/scratch/jingze/w121_3x_ckab_0929 | results -> NOTE on Caltech. Not verifiable from viper. |
| Caltech | (task) | C256 remap + benchmark 8/12 H200 nodes (TASK-2026-09-29-caltech-c256.md); no production | | pending on Caltech |
| DeltaAI | smoke 3271618 then 12 short arms + noise twin (ghx4-interactive, 2 h cap) | MHD max_eta scan (TASK-2026-09-29-deltaai-w121-mhd-etamax.md, NOTE-2026-09-29-deltaai-w121mhd-queue.md) | | state not visible from viper |

Stopped/finished since 09-29 03:15: WASP-121b 1x and 10x reached rot 300 (w121prod_0929, follow-on links
cancelled 09-29 18:30); viper 3x copies 12030240/41 and 12030356/57 cancelled; w3xk first submission
12034146/47 cancelled (key change); He presn R9 (12021350) blew off; C64 x nx1 256 smoke 12027341 done
(~33 KB/cell GPU memory -> C256 needs >= 8 H200 nodes).

## 2. Merges into rt-integration since b3a3f53f (all pushed)

| sha | what changes |
|---|---|
| **d26b7364** m1-perf-0928 | `<rad_m1>/implicit_opac_newton` DEFAULT ON where implicit_opac_update + implicit_gas_newton (He box, sp wedge, MHD); restarts lacking the key keep it off. He box -58/-60 % ms/cycle (1/2 GPUs), Picard 2.03 vs 6.76 (viper check 3aaaefc6). Also implicit_halo_ipc (opt-in), CUDA M1Evt real events. `implicit_one_pass` stays opt-in (user; NOTE-2026-09-29-viper-onepass-cost.md). |
| **ae767d20** dhj-ck-conserve + **fe7301e5** inputs | `problem/ck_impl_conserve` (default 0, solution bitwise): reported ck face fluxes consistent with the applied energy in bound (floorbound/KKT) cells; `problem/budget_dt` per-shell energy/mass ledger (default off); `problem/flux_hst_rkavg` default TRUE = RK-weighted fluid flux hst columns (semantics change). fe7301e5 sets ck_impl_conserve = 1 in the WASP-121b inputs; C256 must use >= ae767d20 with the key. |
| **9e0b204a** units audit | `problem/bbot_gauss` (Gauss) beside `problem/bbot` (Heaviside-Lorentz code units: B_G = bbot sqrt(4 pi); 3 G = 0.846, 10 G = 2.821; prod4 "3 G" was 10.6 G); startup prints both units; dhj cgs guard; rad_m1 unit-constant print. |
| **354a88b7** | He box / He star inputs `tfloor = 5.0e3` (code T, inert) -> `tfloor_kelvin = 5000`; He box smoke bitwise. |
| **23ec1fb8** m1-opn-guard | `implicit_opac_newton_guard` 0.5 + `_guard_mode` 2 (default): drops a face's opacity-Newton term that would take the row diagonal below 0.5x its no-Newton value or flip a neighbour entry sign (M-matrix). guard 0 = old rows bitwise. He box 0 drops, bitwise, same speed; He presn wedge 0 floor clips / 0 NC, Picard 8.1 (NOTE-2026-09-29-m1-opn-guard.md). |
| **ad992f7b** (DeltaAI) | c2p_track Before/After kernels defined once in utils/c2p_track.cpp (CUDA dhj MHD startup segfault from inline header lambdas). |
| **793e03c3** (DeltaAI) | Hydro/MHD ClearSend/ClearRecv also wait on the cubed-sphere seam flux requests (uniform meshes leaked Isend requests -> Cray MPICH abort after ~8100 cycles). Waits only, bitwise. C256 production needs >= 793e03c3. |
| **b6a6eff3** fix-fofc-inline-lambda | FofcTrialRemoveGrav (fofc_etotgrav.cpp) and pgen_eos::HostGamma1FromP (pgen_eos_utils.cpp) defined once in a .cpp; viper gate of the DeltaAI merges + this fix: dhj hydro, He box, dhj MHD fofc 0/1 bitwise vs pre-merge on HIP, 15/15 CPU tests (NOTE-2026-09-30-viper-gate-deltaai-merges.md). Left: multigrid.hpp templates (nvcc risk unclear); CUDA fix untested on GPU. |

Also on rt-integration: 2fd94098 (build_inc.sh mkdir lock), cba4e797 (header order in pgen_eos_utils.cpp).
Not merged: he-presn-m1 (branch, many commits; NHISTORY/NREDUCTION 20 -> 22 is a global change to flag
before merge), dhj-remap-h (script only, 9896d452), dhj-olr-dump (diag, 2efb2851).

## 3. Science results

**WASP-121b 1x rot 300** (w121prod_0929/ana_rot300/RESULTS.md): deep CONVERGED (rot 200-300: 100 bar +0.10 K/10 rot,
bottom -0.04; rotpot fix works). Sub-adiabatic 3-250 bar, only the 2 wall cells superadiabatic: **no deep CZ**.
Jet 14.2 km/s at 4e-4 bar, 3.65 at 1 bar; 0.1-bar day-night 558 K. hst gap +7.2e27 erg/s = ck clamp reporting in the
top 2 cells (p < 1e-7 bar) -> fixed in the ledger by ck_impl_conserve (solution unchanged).

**10x rot 300** (ana_rot300_10x, budget_rot300_10x): the flat -3.8e28 erg/s total-E loss is **REAL deep cooling**
(10-430 bar, ~-0.9 K/100 rot): upward fluid enthalpy flux 1.15e28 (wall) -> 4.4e28 (10 bar) vs Lrad_bot 9.2e27
drains the initial deep adiabat (1x warms slightly instead). Open: resolved circulation vs numerical mixing (-> C256).

**Synthetic observables vs data** (synth_rot300*/RESULTS.md; literature obs_lit/LITERATURE.md):

| quantity | 1x | 10x | observed |
|---|---|---|---|
| phase-curve peak offset | 17-21 deg E | 6-6.5 deg E | ~3 deg E (NIRSpec) |
| NRS1 / NRS2 night (ppm) | 845 / 1291 | 278 / 539 | 136 / 630 |
| dayside NRS2 deficit | -667 ppm | -641 ppm | — (both too low: coarse 3.5-4.4 um ck band / missing opacity) |
| eps (redistribution) | 0.35 | 0.13 | 0.25 |
| T_day / T_night (K) | — | 2789 / 1337 | 2717 / 1562 |
| limb Fe RV morning / evening (km/s) | +2.2 / -13.2 | evening -12.2 | -4.1 / -6.9 |

Bond albedo 0.277 is imposed by the input. 1x too much transport, 10x too little -> 3x test (below). Missing
candidates: drag, nightside clouds.

**3x ck table** (memory w121-3x-ck-table-0930.md): no upstream 3x premixed table; interpolate ln(k/X) per H nucleus
in log Z between 1x and 10x. Validation: Sonora 0 + 1.0 -> 0.5 grey means median 1-3 %, 90 % < 8 %; Exo-FMS
1x + 100x -> 10x median 8-12 %. Worst: optical and 3.5-4.4 um. Files ckdata3/ck/Premixed_3x_g8_11_hiT2.txt +
CE_tables/FastChem_ck_3x_int_hiT2.txt; eos_xh 0.7188, eos_yhe 0.2420, met 0.4771; RCE IC, nx1 76.

**He box cfl verdict** (hebox_cfl2_0927/RESULTS.md): cfl 0.6 = cfl 0.3 within noise (v1'/v_MLT 1.265 vs 1.257);
0.9 marginal deficit (KE1 -3.6..-4.7 %). **User: M1 uses cfl 0.6 in the science window, 0.9 in relaxation.**

**He presn wedge (4 Msun, he_star_m1, branch he-presn-m1)**: vet_col source fix (relaxed) cured the inner-wall
instability; frozen MLT R9 blew off (FeCZ heats before convection grows) -> adaptive shell-mean deficit closure.
3-D: convection grows (e-fold 0.54 turnover) but is inefficient; from 3.5-4 turnovers a **porous runaway** (rho
rms/mean 1.4, <F>/F_diff 2.4, L_top/L 1.16, mass loss 4.5 %/turnover) — physics, robust to cfl, ranks, BE and
theta walls. The solver breakdown at ~5.1 turnovers: **root cause = unguarded opacity-Newton face term** (negative
diagonals in fast thin plume cells) -> guard 23ec1fb8. **Retracted**: the theta-periodic-seam artefact claim
(reflecting-wall run wr3d gives the same box-scale mode within 1-5 %).

## 4. User decisions 09-29 / 09-30

1. m1-perf merged with Newton default (d26b7364); `implicit_one_pass` stays opt-in.
2. He box / M1: cfl 0.6 science, 0.9 relaxation (hesdirk2).
3. He presn: 4 Msun on a sp wedge with a lean pgen (he_star_m1), 2nd order + fully coupled; frozen MLT then ramp
   (option 1) -> replaced by the adaptive shell-mean deficit closure; input back to GUARDED Newton (ae9e5c77).
4. Opacity-Newton guard: go + merged; vimp positivity fallback diagnosis LATER.
5. ck_impl_conserve merged and = 1 in WASP-121b inputs; C256 must use it.
6. C256 x nx1 256 continuation of the 1x (viper bench 16 + 32 nodes; Caltech remap + bench, no production).
7. 1x rot 300 -> resistive MHD max_eta scan; dfloor 1e-13 yes; hlld_bx_zero_tol original (1e-4); split DeltaAI
   (12 short + twin) / viper (long arms); b3_e14_f13 cancelled; viper arms held.
8. Inputs in physical units (bbot_gauss, tfloor_kelvin); merged.
9. Deep-mixing question goes into the C256 analysis (deep enthalpy flux mean/eddy split, deep drift vs C32).
10. 3x: nx1 76; runs on Caltech only; viper w3xk with only dtmax 0.25 + maxit 24; wait for the Caltech A/B.
11. All agents Opus 5.5 medium (worker.md); apudev max 3 chained jobs.

## 5. Open items / next steps

- **C256**: collect 12026299 / 12028087 (viper) and the Caltech benchmarks -> user picks machine and node count.
  The C256 analysis includes the deep-mixing split.
- **3x**: Caltech 3640302 to rot 300 -> synthetic observables vs the table above. Solver-key question: w3xk (viper)
  vs Caltech A/B (3649979/80); viper A/B repeat only if ambiguous. No production key change without the user.
- **MHD scan**: DeltaAI short arms first; then decide on releasing viper arms 12028523-25 (`scontrol release`).
  Analysis `ana_eta.py <armdirs>` in /viper/ptmp2/jinma/w121_mhd_0929.
- **He presn**: 128x128 run (12030056-58) — does the porous runaway persist at resolution, does the guarded solver
  pass ~5 turnovers; He-star / porosity physics question (steady porous state vs runaway) open for the user.
- **Parked**: vimp positivity fallback (milder trigger, 16 stages); ck_impl_conserve extension to the linearised
  cadence steps (optional; 10x residual 9.1e26); dayside NIRSpec deficit (coarse ck band); drag / clouds.
- Merge of he-presn-m1 needs the NHISTORY/NREDUCTION flag and a combined gate.

## 6. Standing rules added

- **ALWAYS SMOKE-TEST before any production/chained launch or relaunch (user 09-30, all machines, agents too):** a
  few-cycle run on the dev partition (viper apudev) with the SAME binary, input and command-line KEYS as the job,
  in a separate dir; submit/release the chain only when it cycles without FATAL. (He presn 12037023/44 died in 9 s
  after hours in the queue: KEYS named rad_m1/implicit_opac_newton, absent from the input -> parser fatal.)

- apudev: at most 3 chained jobs (user 09-29); longer runs go to apu.
- Every input key in physical units named in the key/comment, startup prints code and physical values;
  old "bbot 3 G" statements are code units (10.6 G).
- All delegated agents: subagent_type worker = Opus 5.5 medium (supersedes the Sonnet split).
- M1 cfl 0.6 in the science window, 0.9 in relaxation.
- Standing earlier: 2 GPUs per viper node; HSA_XNACK=1 + HSA_NO_SCRATCH_RECLAIM=1; GPU timings only; ROCm 7.2
  (`*_gpu72` builds); incremental builds (builds/build_inc_viper.sh); rsolver lhllc (hydro) / lhlld (MHD) for dhj;
  noise gates at equal physical time.

## 7. Update 09-30 ~10:20 CEST (read this first; supersedes the sections above where they differ)
- **Merged since 06:40:** defaults flip 9bf6db23 (ck_impl_maxit 24, ck_impl_conserve 1, bbot default 0, sponge_bottom
  false; NOTE-2026-09-30-defaults.md, AUDIT-2026-09-30-default-candidates.md); lhllc/lhlld radial-face fix aa9d51b7
  (NOTE-2026-09-30-lhllc-x1-phi.md): lhllc left an undamped 2-cell RADIAL odd-even velocity mode in the deep (low-Mach
  phi ~3e-3 radially); `lhllc_x1_phi_min`/`lhlld_x1_phi_min` = 1 restores the full dissipation on x1 faces; default 1 for
  deep_hot_jupiter_rt ONLY (global 0; the He box would lose ~4 % v1'/v_MLT). Deep v_r, per-shell Mdot and the deep
  enthalpy flux of every earlier lhllc dhj run are dominated by that mode (deepmix_c32 baseline, 10x deep cooling).
- **Queued on viper apu (est. ~10:45-13:10):** C256 benchmarks 12028087 (16 N), 12026299 (32 N); WASP-121b FRESH starts
  WITH the radial fix and original solver keys: w1xf 12037717(+18), w10xf 12037719(+20), w3xk 12034178(+79) (3x),
  binary athena.gpu.x1f / .w3xk = 44be5b99 (same src as aa9d51b7), dirs /viper/ptmp2/jinma/w121prod_0929/w{1xf,10xf,3xk};
  He presn with **hllc** (user: near-sonic convection, Mach ~0.4): 64x64 12037044(+45) and 128x128 12037023-25
  (M1/ad3d_hllc, M1/ad3d_128). HELD: viper MHD arms 12028523-25.
- **Caltech:** 3x (original keys, pre-fix) 3640302 running, ETA rot 300 ~15:00 CEST; ck A/B: no limit cycle
  (494540fd) -> the 3x dt decline is not a solver effect; compare with w3xk (radial fix) once it runs.
- **Next:** C256 benchmark -> machine choice; He presn hllc 64 vs 128; fixed arms vs pre-fix (dt, odd-even share, deep
  flux via deepmix.py); DeltaAI MHD max_eta scan; then the 3x/10x MHD drag set (PLAN-2026-09-30-w121-mhd-drag-metallicity.md).

## 8. Update 09-30 ~23:30 CEST (read this first; the previous session ended here)
**The old session stayed open (screen) to let its agents finish; if an agent dir has NO RESULTS.md, that agent was cut off. Their detached drivers keep computing (synth2_trans/driver.sh -> "ALLDONE" in driver.log; synth2_emis/run_prt.sh -> 3 pRT runs, logs/), Slurm jobs keep running, but NO agent will write RESULTS.md: the new session runs the analysis itself (synth2_emis/products.py, synth2_trans/ana.py, hepresn_wind_0930/analyze_wind.py; read each script header first) or relaunches an agent with the original brief goals.**

| agent | results land in | ETA |
|---|---|---|
| WASP-121b emission post-processor (pRT, angle-resolved, instrument bands) | /viper/ptmp2/jinma/w121prod_0929/synth2_emis/ (RESULTS.md, out/, logs/) | ~00:10-00:20 |
| WASP-121b transmission / Fe-Na limb RV (hi-res chords, CCF) | /viper/ptmp2/jinma/w121prod_0929/synth2_trans/ (RESULTS.md, out/, driver.log ALLDONE) | ~01:15-01:30 |
| He-star steady wind, step 1 (branch he-wind-bc 697be6eb, not pushed) | /viper/ptmp2/jinma/hepresn_wind_0930/ (RESULTS.md, W*/G* runs; GPU job hewind 12046760) | open |
| M1 positivity, step 2 (branch m1-positivity 38500817; fork has 33150622) | FINISHED: tests_m1/runs_6b_m1pos/RESULTS.md on the branch; runs /viper/ptmp2/jinma/m1pos_0930 (`python3 analyze.py`) | arms running |

**He presn (4 Msun wedge):** 128x128 NaN at t 26,825 s (5.7 tt), 64x64 hllc NaN at 31,626 s (6.7 tt): implicit M1 E<0 in the
porous Fe-bump outflow. dt collapse there = hydro radial CFL of the outflow channels (v_r ~970 km/s > v_esc, rho 0.5-4 %
of shell mean). Fix on m1-positivity (keys default off, OFF bitwise on wedge + He box, He box GPU cost unchanged):
implicit_g0_exchange, implicit_g0_limit=1, implicit_opac_newton_guard_mode=6, implicit_pos_gas, implicit_pos_floor.
First trigger = Picard pass-0 solve E<0 from the lagged -c dt v g0 face term. Restart gates from the last good rst:
controls (off, BE-off) NaN at 26,825 s again; ON arms pass. Speed/stability at 128: **BE + options (beon128) 12,600
sim-s/wall-h, 0 stage rejects; hesdirk2 + options 6,000-7,300 sim-s/h with 74-76 % of stages rejected -> BE redo.**
Recommendation (not decided): BE + options in the transient/relaxation, hesdirk2 only at steady state.
Steadiness: mass loss peaked ~-7 %/tt at ~4.9 tt and fell to ~-0.01 %/tt by 6.4-6.9 tt at BOTH resolutions (burst,
not a steady wind); KE and L_top/L_in still drifting after the closure ramp ended at 28,200 s (L_top/L 0.96 at 6.9 tt).
Caveat: the mass-loss decline coincides with the closure ramp start (23,500 s).
**Queued (user 09-30):** beon128 extension to t = 47,000 s (10 tt): job 12046848 (afterany 12046528 + smoke 12046798;
the script exits unless the smoke log shows `rc=0 fatal=0`). Other arms: onv128 12046525, onn128 12046526, onnh64
12046530, offh64 12046529 (control), He box accuracy 12046701 (`hebox/ana.sh 60000 74000`).
**Open for the user:** implicit_vimp on/off in relaxation; making the five keys the He-star default (after the 10 tt
run + He box accuracy); hydro-side eint<=0 in near-void cells (rho/<rho> 2e-5); step 3 (3-D with wind BCs + keys).
**WASP-121b fixed arms (radial fix) at 23:00:** w1xf rot ~74 dt 14 s; w3xk rot ~46 dt 7-7.7 s (declined from 12 s at rot
18-20 like Caltech, radial fix does not prevent it); w10xf rot ~43 dt 9.7 s; all clean. C256 32-node bench 12026299
still pending (16-node: 0.118 s/cycle, dt 3.7 s -> ~1 h / ~31 GPU-h per rotation). MHD arms 12028523-25 HELD.
**New standing rule (user 09-30):** ALWAYS smoke-test (same binary + input + KEYS) before any launch/relaunch; in
worker.md and memory. User granted `scontrol hold/release` (allow rule in .claude/settings.local.json; the auto-mode
classifier still blocked it twice in this session - may need /hooks or a restart to load).
Orion task TASK-2026-09-30-orion-he-m1pos DROPPED (no transfer path). WASP-121b analysis scripts pushed:
docs/handover/scripts/w121/.
**Secondary hot spot (user question 09-30 23:40):** WASP-121b 1x rot-300 (dumps 146-150) has a compact equatorial hot
knot at lon +110..+120 deg (past the evening terminator) at 1e-5..1e-4 bar, ~3300 K vs ~1500 K night; gone by 1e-3 bar;
absent at 10x. Maps /viper/ptmp2/jinma/w121prod_0929/hotspot2/ (Tmaps_rot300.png, maps_w1x/w10x.npz, maps.py; script copy
docs/handover/scripts/w121/hotspot2/). Masked-knot phase-curve test queued with the emission agent (RESULTS.md section in
synth2_emis/); my estimate before it: at most tens of ppm, only in strong band cores (CO 4.5-4.8 um).
**Knot heating check (09-30 ~23:55, hotspot2/knot.py, equator |lat|<8, 1e-4 bar):** along the equator the knot is only a
~100-200 K bump on the jet's cooling track (T 2950 K at +102 -> 3060 K at +110..+118 -> 2690 K at +126); its ~1500 K
contrast is vs the cold OFF-equator night at the same lon (a hot equatorial tongue carried by the jet). Cause: strong
DOWNWELLING (v_r -2 km/s vs |v_r| < 0.5 km/s around) = adiabatic compression; NOT a shock (u stays 13-14 km/s, zonal Mach
2.4-2.6, no deceleration at 4-deg binning) and NOT H recombination (mu 1.26 = atomic through the knot; recombination, mu
rising to 2.3, starts downstream at +126..+174). At 1e-5 bar downwelling -6.5 km/s: near the top sponge (1e-7..1e-6 bar),
possible top influence. The jet is supersonic everywhere on the night side (Mach 2-6).
**Emission post-processor DONE (09-30 ~23:50):** /viper/ptmp2/jinma/w121prod_0929/synth2_emis/RESULTS.md (scripts in
docs/handover/scripts/w121/synth2_emis/, uncommitted there). pRT 3.4 R1000 ck + FastChem (ions, dissociation), PHOENIX
star, real throughputs, 72 phases, angle-dependent per-column rays (no scattering). **Main finding: GCM closure bias**
(memory ck-spherical-closure-overemits-0930): ray RT on the GCM's own table gives 0.80 (1x) / 0.87 (10x) of the GCM OLR;
the ck_spherical face closure reproduces it (0.96/1.02). New vs observed (1x | 10x | obs): NRS1 day/night 3219/805 |
3479/326 | 3924/136; NRS2 3830/1257 | 4261/632 | 4924/630; offsets NRS1/2 17.8/18.6 | 6.4/7.0 | 3.4/2.7 deg E; T_day/
T_night 2488/1622 | 2667/1300 | 2717/1562 K; eps 0.370 | 0.138 | 0.246. Offsets unchanged by the post-processing (< 1 deg).
Masked knot: max -47 ppm (CO 4.5-4.8 um), -36 NRS1, -28 NRS2 at phase 0.194; offsets move 0.3-0.9 deg. Open: fix the
closure (user), scattering/reflected light (user), 416 MB Chrome cache /viper/u2/jinma/.cache/selenium (deletion was
blocked for the agent; user to decide).
**10-01 ~00:20 ck_spherical closure over-emission CONFIRMED in the real kernel** (cksph_test_0930/RESULTS.md; memory
ck-spherical-closure-overemits-0930): 2x/(1+x) exactly; W121 1x isothermal 1.22-1.25. Prototype problem/ck_sph_dilute on branch
ck-sph-closure-fix (not pushed; tm + ck_implicit only, refuses ck_impl_lin/jac_lin which production uses). User decisions:
adopt + port + GPU A/B; thick-cell blend; variable-Eddington factor; top-ghost boundary value.
**10-01 ~00:20 He presn beon128 is NOT steady: the envelope COLLAPSED.** Mass above 0.7 R 71 % (3.0 tt) -> 60 % (4.5 tt) ->
1 % (7.5 tt); median rho at 0.9 R 3.6e-9 -> 1.6e-13; 24 % floor cells above 0.8 R; only 7.7 % of the mass left the domain
(fell onto a dense shell at 0.5-0.65 R above the closed inner wall). Large inward-flux patches in the empty outer region
(suspect). Slices page https://claude.ai/artifact/HcgYDVQ3CozhxXavZppqcU ; script /viper/ptmp2/jinma/m1pos_0930/slices/slices.py.
Open: physical (porosity -> radiation no longer supports the inflated envelope) vs closure ramp-off / closed wall / floor.
**10-01 ~00:30 Transmission post-processor DONE** (/viper/ptmp2/jinma/w121prod_0929/synth2_trans/RESULTS.md; scripts copied to
docs/handover/scripts/w121/synth2_trans/): Kurucz Fe I (300-line mask) + Na D Voigt, H- + Rayleigh continuum, R 600k chords,
15 epochs, CCF at R 70k. Fe RV morning / evening / whole: 1x -0.66 / -11.86 / -6.75, 10x +0.48 / -11.65 / -6.48 vs
Seidel+25 -4.12 / -6.90 (whole ~-5.2 to -6): evening ~5 km/s too blue, morning not blueshifted; only the whole-transit value
agrees. Decomposition: winds alone -4.9 / -6.0, rotation +5.0 / -5.8. Lines form at 10^-4.5 (1x) / 10^-5.2 (10x) bar median;
domain top high enough (< 3 % above it). Robust to Fe x0.1/x10, extension, chemistry, LSF. Lever = circulation (drag / less
super-rotation), not post-processing. Line depths (Fe CCF ~1500 ppm, Na ~3000 ppm) not checked against measured amplitudes.
**Session end 10-01 00:20.** All agents reported except the He steady-wind agent (step 1, branch he-wind-bc; last activity
arms_apu.sh + smoke_FH at ~00:15): check /viper/ptmp2/jinma/hepresn_wind_0930/RESULTS.md; if missing, the agent was cut off ->
run analyze_wind.py on its W*/G*/arm dirs. Next-session checks: beon128 ext to 47,000 s (does the collapsed state persist?),
He box accuracy (m1pos_0930/hebox/ana.sh 60000 74000), onn/onv128 + offh64 ends (m1pos_0930/analyze.py), C256 bench 12026299.

## 9. Update 10-01 ~00:40 CEST
- He steady-wind agent was cut off; the session ran analyze_wind.py: /viper/ptmp2/jinma/hepresn_wind_0930/RESULTS.md.
  **No steady wind in any 1-D arm** (A wall, B frozen, D wall+frozen, E relax300): the beta-law IC wind drains in ~7,000 s,
  r > 1.2 R becomes near-vacuum infall; Gamma above R 0.2-0.5, below R <1 on average -> continuum driving cannot launch
  a wind; envelope bursts (M(r<R) -25 %) fall back. Smoke F/H (hesdirk2 0.6; nx1 660) rc=1 in 1 s, never launched.
- offh64 12046529 was held ("user env retrieval failed"); released 00:35.
- beon128 ext at 37,713 s, clean (0 NaN); ETA 47,000 s ~01:30. onn128 32,845 / onv128 31,713 / onnh64 32,981 s, clean.
- offh64 12046529 cancelled 00:45 (user: control redundant; off128 controls already NaN at 26,825 s).
- 10-01 02:20: beon128 ext reached 47,000 s (10 tt), clean; envelope collapse PERSISTS (M(>0.7R) 0.6 %, rho(0.9R) still falling).
  He box gate: positivity keys accuracy-neutral (N-C < P-C round-off). onv128 vs onn128: vimp on ~33 % slower, no gain.
  New threads 10-01: ck closure variants (cksph_test_0930/VARIANTS.md; recommend face 5 VEF + top 1, GPU timing running),
  He HSE check (he_hse_1001), BSG reproduction of Ma+2026 (bsg_1001; step 3 running on branch bsg-setup).
