# Viper session 2026-09-27 night to 09-29 (supersedes SESSION-2026-09-27b-viper.md) — START HERE on viper

Read this, then `MEMORY.md`. **Caltech / DeltaAI: read section 1 before running ANY dhj / WASP-121b case.**
`rt-integration` at the end of this session: **6d690e09** (pushed). Always `git fetch fork` before pushing; never force.

## 1. Bugs fixed (all merged on rt-integration)

| fix | effect | commit |
|---|---|---|
| **Centrifugal energy double count** (dhj, `rot_potential` + `etotgrav`): the horizontal centrifugal work was added in SourceFunc although the etotgrav flux of Phi_tot already does it | ~1e31 erg/s taken from 0.1-10 bar and put below 10 bar: the 100-bar warming (+28 K/10 rot) and 1-10 bar cooling of every rot_potential+etotgrav run (w121 1x/10x, likely prod3/prod4) were this bug. Fixed: 100-bar drift -1.7 K/10 rot, energy budget closes. **Runs before the fix are unusable below 1 bar after rot ~20; the bug also froze a sub-adiabatic layer below ~60 bar (thermal time ~1e6 rot), so do NOT restart old runs: start fresh.** | 8fa76784 (+ tst tolerance 87d5a244) |
| He box `wall_noflux` subtracted rho*Phi with the EFFECTIVE WB potential from an energy holding the TRUE one | wall-cell e 2.9x too large at the top wall; globally ~1e-9 of L; + guards refusing wb_phi_eff combos that are not energy-consistent | ab556bf6 |
| FOFC under etotgrav tested e + rho*Phi as internal energy (hydro + MHD) | energy-floor hits invisible to FOFC in red giant / He4; now the trial energy is what the real update produces (He4: FOFC flags 3e3 -> 1.25e7, floor events down) | 96f8a7f0 |
| energy guards: ghost rho*Phi in built-in copy BCs, SMR/AMR + etotgrav FATAL, restart potential-key check (`problem/allow_potential_change`), rank-summed RT limiter warning | | 96f8a7f0 |
| cs resistive curl on stretched radial grids (x1-face derivative, `<mhd>/cs_resist_x1_centred`, default on) | WASP-121b-type stretch: field error 1.2e-2 -> 2.3e-4; uniform grids bitwise | 50c0a706 |
| ck Jacobian fp atomics -> ordered per-block sums (`jac_lin/lin/frozen_op` off paths; T4 unaffected) | GPU repeat runs bitwise; +2.6-3.6 % on those old paths only | 2a7c083d |
| M1 on MHD: box two-stream / wb_arad_force and the simple M1 tests ported to <mhd>; closure echo | | 83bd1623 |

## 2. Defaults changed (merged; docs/dev/default_flips_0927.md rounds 2 and 3)

- dhj: `ck_implicit` true (fresh rt_ck runs; restarts keep their recorded value), `ck_impl_floorbound` / `ck_impl_kkt_demax` true (fused path), `ck_impl_xstep` 0 (relaxation-phase option; the WASP inputs name 8).
- M1: `implicit_precond` rbgs_fwd on spherical polar; `implicit_op_team_red`, `implicit_vimp_fold`, `implicit_halo_ovl_faces` off, `implicit_one_pass` 0 everywhere; `time2_lin_tol_fac` 1 (was 10: stiff radwave time order 0.12-0.48).
- Kept on after the MI300A audits: `ck_store_split`, `ck_beam_par`, `implicit_halo_overlap`, `implicit_mg_levels 3`, `implicit_gas_newton` (user).
- Diagnostics: `problem/flux_hst_floor` (Efloor, Efloor_rt, Mfloor hst columns), `ck_impl_ncref`, `ck_impl_osc`.

## 3. Verdicts

- **He box cfl (M1 coupling fix):** KE dt-converged cfl 0.15-0.9 (growth 1.878e-4/s every arm; old cfl-dependent KE was numerical). Strict round-off: cfl 0.9 relaxation-only (KE profile 60x noise), 0.6 ~7x, 0.3 converged. Box not convecting yet: continuation to saturation RUNNING (hebox_cfl2_0927, R3/H9/P9).
- **WASP-121b bottom BC:** fresh 1x on the fixed code, 25 rot: deep stays on its starting adiabat (entropy within 0.005); fixed-flux closed wall is enough, no fixed-entropy BC. Thermal time below 10 bar ~8e5 rot => the IC adiabat IS the interior entropy; **user: keep the current IC adiabat**.
- **10x ck Newton:** input keys `ck_impl_dtmax 0.25` + `ck_impl_maxit 24` -> 0 non-converged (was 144/500); 1x `maxit 16` -> 0. A fresh start with the old keys: 1x 19/75, 10x 69/75 non-converged.
- **200 K cells:** at ~1e-9 bar at the night-side top edge (not in the region that counts); ck cools optically thin unirradiated cells below the floor; no bug; floor energy ~1e-4 of the absorbed starlight.
- **Top-shell "leak":** ConToPrim floor energy (~1e25 erg/s), not a leak.
- **Caltech ck refactors (ck-next, ck-lin2-fma):** exact reformulations; HIP noise gate PASSES at equal physical time; keep. MI300A: -2.7 % / +1.5 %. ck-lin2-hip (V1) dropped.
- **Noise gates must stop arms at equal PHYSICAL time (time/tlim), never nlim** (an nlim gate faked a 4.6x failure).
- **MHD He box v_h order loss:** an additive ~0.1 cm/s error, negligible; hlld is 3-60x more dissipative at low Mach (keep lhlld); lhlld's low-Mach correction has a magnetic floor (MHD damps slow flows more than hydro) — for the MHD phase.
- **Rm_dx (prod4):** explicit eos eta unresolved at the grid scale; user: too pessimistic; follow-ups (drag time, effective numerical eta) deferred to the MHD phase.
- **Validation of the tip:** no regressions; pre-existing: MPI+AMR segfault/deadlock (dyngrmhd, z4c, jeans, rad lwave AMR), multigrid 900 GiB allocation on GPU (upstream cb39c07c), flaky in_cshock2d; GR-family tests skipped (user: GR EOS broken).

## 4. WASP-121b fresh start — ready, awaiting the user's go

Tip 6d690e09; inputs = w121prod_0927 1x/10x plus: 1x `ck_impl_maxit 16`; 10x `ck_impl_maxit 24`, `ck_impl_dtmax 0.25`;
`ck_impl_xstep 8` named (spin-up); optional `flux_hst_floor = true`. Current IC; closed wall (default). 1 node x 2 GPUs per arm;
~24 rot/h at 1x (bench). Plan: coarse spin-up to rot 300, then remap to nx1 256 (remap port RUNNING: /viper/ptmp2/jinma/remap_0929).
Recheck the deep entropy step at ~rot 100; switch xstep 8 -> 0 with a drift check at the end of the spin-up; decide sponge_bottom.

## 5. Infrastructure

- **Incremental builds:** `/viper/ptmp2/jinma/builds/build_inc_viper.sh <target> <commit>` (README there; dhj/box/none/rg GPU, dhj/box/none CPU); ~3 min for a pgen change vs ~8-25 min full. Note `docs/handover/NOTE-2026-09-28-incremental-builds.md` for Caltech/DeltaAI.
- **Cross-cluster timing:** `docs/handover/TASK-2026-09-28-cross-cluster-timing.md` + `bench-2026-09-28/` (pinned 11c9a5be); viper results in `RESULTS_viper.md` (1x 19.2 / 12.7 ms/cycle at 2 / 4 GPUs).
- **GPU tuning (RUNNING, /viper/ptmp2/jinma/gputune_0928):** ROCm 7.2 (+ HIP malloc async off) looks ~9 % faster per cycle; equal-time noise test in progress; build-script change needs the user's OK.
- viper: no CPU partitions for us; CPU checks on the login node only if <= 30 min (nice, <= 16 cores); apudev 2 nodes; tst `--gpu` accepts HIP flags; style script no longer clobbers cpplint.py.
- DeltaAI: `TASK-2026-09-28-deltaai-bringup.md`; token cleanup note `NOTE-2026-09-28-token-cleanup.md` (the old GitHub token is revoked).

## 6. Open, later

Massive stars: global He-star model geometry (sp wedge with vet_col, or cubed sphere after M1 stage S5). MHD phase: max_eta cap / drag time, effective numerical eta, lhlld magnetic floor, GS05 corner-emf time order, small wall Mdot in MHD. Pre-existing bugs listed in section 3.

## 7. Update 09-29 ~03:15 (read this first; supersedes section 4 where they differ)

- **WASP-121b fresh start RUNNING on viper** since 09-29 02:08 CEST: `/viper/ptmp2/jinma/w121prod_0929` (1x job 12018387 +
  12018388 afterany, 10x 12018389 + 12018390), ROCm 7.2 binary of cda4da33 (md5 7bf3f753), inputs = w121prod_0927 +
  `rsolver = lhllc`, 1x maxit 16, 10x maxit 24 + dtmax 0.25, xstep 8, flux_hst_floor; ~25 rot/h (1x), ~14 rot/h (10x).
  The same run is also on Caltech (restarted with lhllc, smoke clean) and queued on DeltaAI; the user decides which
  copies to keep (NOTE-2026-09-29-viper-w121prod.md, NOTE-2026-09-28-deltaai-w121prod.md).
- **dhj uses lhllc (hydro) / lhlld (MHD) from now on** (user 09-29; NOTE-2026-09-29-dhj-rsolver.md).
- **ROCm 7.2 adopted on viper** (gcc/16 rocm/7.2 openmpi_gpu/5.0 + `Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF`; ~10 %
  faster; build targets `*_gpu72` in /viper/ptmp2/jinma/builds; jobs must load the same modules).
- **Radial remap merged** (e804daf1): `problem/remap_file` + docs/handover/scripts/dhj_remap.py; HOWTO in
  /viper/ptmp2/jinma/remap_0929/HOWTO.md (76 -> nx1 256 at rot 300; 17.6 min/rot on 2 GPUs at 256).
- **He box:** the M1 box DOES convect (v1'/v_MLT 1.2-1.5 saturated, same as two-stream box_w8; the earlier "no
  convection" was a measurement error: plane-mean flux + cm/s read as km/s; heconv_0929). Resolution sufficient (no
  doubling). Saturated cfl comparison R3/H9/P9 pending in /viper/ptmp2/jinma/hebox_cfl2_0927 (R3 job 12016075 near
  t = 74000; results in RESULTS.md, script scripts/analyze2.sh, metrics res/conv.py).
- **Benchmarks:** WASP-121b (bench-2026-09-28; viper, DeltaAI done: GH200 1.4-1.7x MI300A) and He box
  (bench-2026-09-29-hebox; viper done, hold lifted for Caltech/DeltaAI).
- **Agents:** simple mechanical tasks on Sonnet ("sonnet" = Claude Sonnet 5, claude-sonnet-5, checked 09-29), judgement on
  Opus 5.5. Noise gates stop at equal physical time (tlim). CPU checks on the login node (<= 30 min).
