# TASK 2026-11-17 Orion -> viper (USER): staggered clean scheme, GPU stage 1 (A100 + MI300A)

The USER decided to hand over the staggered layout of the clean scheme ("stag2") once the
positivity and drag/preconditioner rounds passed. Both have passed, and so have the follow-up
fixes (beam gap, F-solve stop rule, M7 opacity key, cubed sphere package). This TASK asks for
the GPU build, a short correctness check against Orion CPU references, and **timing against the
cell layout, which is the point of the staggered layout** (on the CPU it is 1.2x to 12x cheaper
per cycle than the cell layout on the same decks).

## Code
Fork `jing-ze-ma/athenak`, branch **`clean-cc-stag2-gpu1` = 9689f4216**. It is your
push-16 candidate **0678599eb** plus the staggered work, replayed on top (77 commits):
- rounds 1-17 of the staggered scheme (d6b6ed4e2..788e3391c);
- the audit / production-review fixes (CleanGasEint use, area-form transverse stress,
  Doppler with the gas's face F, eint floor check, ...), fix2 (dt_retry checkpoint of the clean
  state, d memory on acceptance, D_ab ghost parity at walls/poles), force_reference = clean
  on acceptance (one rule for both layouts, = your 4ac0c0e49);
- positivity: flux factors (`clean_stag_gflux`), the write-back fix, the beam-gap fix (CLN_JR);
- drag/preconditioner: `clean_stag_vimp`, `clean_stag_upwind`, `clean_stag_pc = mg`;
- perf: F-solve factors once (S2), face weights once per pass (S3), fused ghost write-back
  (S6a), F-solve stop rule (S1, new `CLEAN FSOLVE` log line);
- `clean_stag_opac` (M7, opt-in, experimental);
- the cubed-sphere staggered package (included, not gated here; its rhse_cs WB is still 3x the
  cell layout's and being fixed);
- cleanup: study keys removed, the recommended set made the staggered default.
Integration log with every conflict resolution: `task_1117_stag2_gpu1/INTEGRATION.md`.

**Keys off, the branch is bitwise 0678599eb** (Orion reg harness, 41/41 run dirs, x1line
variants 5/5, all 7 RT tests, wb levels, 0 CLEAN WARNING; and a cell-layout restart file
written by 0678599eb continues bitwise on the branch). Your cell-layout defaults are unchanged.

## The staggered key set (now the defaults of `clean_layout = staggered`)
Every deck in `decks/` has the keys written out (the explicit decks reproduce the defaults
bitwise on Orion):
```
<rad_m1>
scheme = clean
clean_layout = staggered
clean_stag_split = true        # energy flux apart from the face momentum state
clean_consistency = true       # default where the rays are not opt-in (not cs; sp only with vet_gd)
clean_stag_fsolve = true       # the face momentum of its own F solve
clean_stag_fkin = true         # the rays' kinetic momentum flux in the F solve ...
clean_stag_fwork = true        # ... with the F solve in the Newton loop (fkin <-> fwork)
clean_stag_gflux = true        # the rays' face flux as flux factors (positivity)
clean_stag_vimp = true         # gas velocity implicit in the face momentum (stiff drag; off on cs)
clean_stag_upwind = true       # donor-cell advective flux + lagged correction, open-end a E
clean_stag_pc = mg             # the old production mg on the E solve (line: x1-line solve)
clean_pos_limit = false        # USER: off for staggered (cell default stays true)
clean_rays = fv1               # spherical polar with vet_gd (unchanged default)
# clean_stag_opac = 0          # experimental, opt-in (M7: opacities of the iterate's T)
```
A named key that its coupling drops prints one notice (`<rad_m1> clean_stag_fkin and
clean_stag_fwork need each other, ...: ignored`). Log lines: `CLEAN cycle ... solves N` per
solve, `CLEAN FSOLVE` per CleanSolve with transverse F-solve passes.

## Build
- A100 (Raven): `./build_cuda_raven.sh <worktree at clean-cc-stag2-gpu1>` (your production
  recipe, gcc/13 cuda/12.6 openmpi_gpu/5.0, AMPERE80, `-ffp-contract=off` -> `--fmad=false`).
- MI300A: `./build_hip_mi300a.sh <worktree>` (the stage-1 MI300A recipe: gcc/16 rocm/7.2
  openmpi_gpu/5.0, hipcc, GFX942_APU, `-ffp-contract=off`).
- NMAX 432 is enough for every deck here.
- Orion CUDA compile check (cuda/12.1, sm_80, nvcc_wrapper, --fmad=false): see "GPU readiness".
  Orion has no A100 and no ROCm; the NVVM 12.6/12.8 miscompile class cannot be run here.

## S1-S3 correctness (one GPU unless noted; minutes per run)
```
LAUNCH='srun -n %N --gpus-per-task=1 <binding>' KEEPGOING=1 ./run_all.sh <build>/src/athena S1 S2 S3
python3 compare.py > gates.md
```
| group | runs | what they test |
|---|---|---|
| S1 | cb, cbg, ws66, shd, shdg, hcyl (x0.5), nltek1 (x1) | Cartesian beams (flux factors, beam gap closed: cb 0.46 = cell), shadows, thick diffusion, non-LTE |
| S2 | rhse_sp2, sp3, sp3_mb2 (2 GPUs), flz q, pulse0 | sp with gas: hydrostatic WB with force_reference = clean, 1 vs 2 radial blocks, rotating box (L_z), closed box energy |
| S3 | rw +-amp (P 3.3, tau 10, N64), drag (beta 0.1, R 3, tau 10, icfl 4), shock M7 N256, fs | radiation wave, stiff drag (vimp), Lowrie-Edwards radiative shock (upwind), radial free streaming |

Gates (compare.py, CPU numbers in ref/ref.json, final dumps in ref/dumps/, md5 in ref/md5.txt):
- rc 0, no FATAL/NaN; **0 CLEAN WARNING, 0 DT_RETRY** (the CPU has 0 on every S run);
- CLEAN cycle lines within 10 % of the CPU's (same solve sequence);
- final-dump field maxrel (vector components share one scale) <= max(1e-6, 10 x the CPU's own
  1- vs 4-thread maxrel). The CPU spread is <= 2e-10 everywhere except nltek1 (2e-10), drag
  (7e-14) and fs (9e-6 -> gate 9e-5).
- These are tolerance gates, not bitwise: GPU reductions change the Krylov iterates.

Orion CPU reference (4 OpenMP threads, branch binary md5 bd94f61ea366c9bf3538a9be28cd9714):

| run | CLEAN lines / solves | CPU 1-vs-4-thread maxrel | CPU wall s |
|---|---|---|---|
| s1_cb_x0.5 | 19 / 19 | 4.2e-11 | 4.3 |
| s1_cbg_x0.5 | 19 / 19 | 5.7e-15 | 5.4 |
| s1_ws66_x0.5 | 19 / 19 | 2.5e-14 | 6.4 |
| s1_shd_x0.5 | 11 / 11 | 1.0e-13 | 5.5 |
| s1_shdg_x0.5 | 11 / 11 | 6.0e-16 | 10.1 |
| s1_hcyl_x0.5 | 65 / 65 | 4.4e-16 | 22.3 |
| s1_nltek1_x1 | 119 / 119 | 2.1e-10 | 19.9 |
| s2_rhse_sp2 | 1316 / 4600 | 6.5e-14 | 10.0 |
| s2_sp3 | 381 / 1730 | 2.4e-14 | 21.3 |
| s2_sp3_mb2 (2 ranks) | 381 / 1743 | 1.9e-14 | 15.7 |
| s2_flzq | 42 / 188 | 7.9e-11 | 23.1 |
| s2_pulse0 | 37 / 73 | 5.2e-11 | 12.3 |
| s3_rwp / s3_rwm | 257 / 772, 776 | 0 | 4.9, 4.5 |
| s3_drag | 131 / 153 | 7.0e-14 | 4.8 |
| s3_shock_m7 | 398 / 3585 | 1.3e-11 | 8 (main loop) |
| s3_fs | 57 / 530 | 9.1e-6 | 12.7 |

`compare.py --runs <Orion 1-thread runs>` gives "all gates passed".

## T timing (the main deliverable; report only, no gates)
Per-cycle cost of the staggered layout against the cell layout on the same deck, same GPU:
```
LAUNCH=... KEEPGOING=1 TIMEOUT=30m ./run_all.sh <build>/src/athena T     # one run at a time
python3 timing.py > timing.txt          # s/cycle over cycles 6-N, staggered / cell ratio
```
| pair | deck | size | nlim |
|---|---|---|---|
| t_shd_x2_{stag,cell} | Cartesian shadow (shd), 1 block | 1024 x 512 | 10 |
| t_sp1x1h_{stag,cell} | rhse envelope, 3-D sp wedge, GL 10x12 vet_gd rays, gas, gravity (the stage-1 G5 sp1x1h deck) | 548 x 32 x 64 | 25 |

Please run both pairs on **one A100 and one MI300A**, and send s/cycle, the staggered/cell
ratio, the initialisation time, the kernel table (kokkos-tools kernel timer or nsys) for each
staggered run (kernels of interest: cln_rhs, cln_apply, cln_fwc, cln_pcsetup, cln_thomas,
cls_*, cln_fs0..3 (F solve), cln_velocity, cln_real, cln_hcopy/cln_hback, cln_fv*), and the
device memory high-water mark. Orion CPU smoke of the T decks (nlim 2, 32 threads, all three at once on one node): t_shd_x2
staggered 4.1 s against cell 67 s of main-loop time, t_sp1x1h_stag 15.0 s, 0 warnings (the cell
sp1x1h was not run at full size on the CPU: its steady start alone takes minutes). Orion CPU
(16 threads, 548 x 8 x 16 wedge, 3 cycles): staggered
2.5 s against cell 5.0 s of main-loop time, and its initialisation is shorter (13.6 s against
488 s wall, dominated by the cell layout's steady start). Bench s/cycle on the CPU (4 threads,
same node, all arms together; int/bench/tab_*.md): cb 0.18 vs 1.26, shd 0.61 vs 7.5, shdg 2.1 vs
8.9, hcyl 0.72 vs 1.68, nltek1 0.30 vs 0.71, sp3 0.12 vs 0.16, flz q 0.85 vs 1.33 (27 vs 40
cycles), rhse_sp2 0.015 vs 0.017, pulse0 equal.

## Known limits (staggered; none block the GPU stage)
1. **Cubed sphere**: included (vimp switched off there with a notice; consistency opt-in), not
   gated in this TASK; rhse_cs WB is 3x the cell layout's (9.8e-9 vs 3.4e-9) and cs6 vs cs24
   disagree: being fixed on Orion.
2. **clean_stag_opac (M7) is experimental, opt-in.** With a T-dependent opacity (kappa ~ T^-3
   relaxation at implicit_cfl 4) the default set needs dt_retry at cycle 0 (9 CLEAN WARNING,
   recovers, rc 0) where the cell layout needs none; with clean_stag_opac = 0 it fails there,
   and it is worse on the T-dependent Marshak cases c/d. On constant opacity the key is a
   bitwise no-op. All S decks have constant opacity.
3. **pulse0 closed-box energy** is closed: 4.4e-16 on this branch (cell 8.9e-16). For the record:
   2.0e-15 with pc = mg and 2.4e-14 with pc = line on the earlier 3cff6ccbb; the 6.6e-14 of
   round 18 belonged to the line preconditioner.
4. **rhse_sp2 well-balancing**: 1.2e-8 with pc = mg (2.0e-9 with pc = line), cell 1.5e-9:
   solver-tolerance level, not a regression class (that is 1e-6..1e-3).
5. **flz q wall v_r** points outward (+0.003/+0.010/+0.009) where both parents of the merge
   had it inward (fix branch -0.0022, round 18 -0.0048; cell -0.0032); |v_r| is at the cell
   layout's level and dE/E, dLz/Lz stay at round-off. Under investigation on Orion.
6. **clean_rot_ring restarts are not bitwise** (cell DIFF 3e-3 identical to 0678599eb; staggered
   1e-1). Pre-existing; the RT tests here run without rotation.
7. **Restart files with clean extrap/ring history**: the cell layout keeps push 16's format
   (49 closure components; 0678599eb files continue bitwise); staggered restart files carry the
   50th (CLN_JR) and are new.
8. F-solve stop rule (S1): with gas coupling the transverse passes go up to 256 (stall exit);
   watch the cost of cln_fs* on gas-coupled multi-D decks (t_sp1x1h).

## GPU readiness (Orion)
- Review of every new staggered device kernel against the A100 NVVM classes (no runtime-indexed
  local arrays, no pointer/reference selects between locals or captured structs, no std::array
  in device lambdas), the HIP level-1 team-scratch race class, and implicit this-capture:
  **no team scratch, no TeamPolicy, no atomics, no member use in any device lambda.** The line
  solves (E-solve Thomas/mg smoother cln_thomas, F solve cln_fs2f/cln_fs2) are one thread per
  line. Fixed during integration: CleanGhostFill's captured `bool rf[6]` (indexed by side) is
  now an int bitmask.
- Known risks to watch (no change made, low risk):
  - cls_galerkin (rad_m1_clean_solver.cpp, CleanStagCoarseSetup): local `am[3]`, `ap[3]` indexed
    by the trip-3 loop variable d (the cell layout's clc_galerkin form accumulates in the View);
  - cls_restrict: `pos[3]` indexed by d (same form as the accepted cell residual);
  - CleanMomentumChange: `Real *mag`, and in the cs package `cf[side]`/`wf[side]` arrays written
    through pointers with a runtime side;
  - viper's inherited CleanRayFlux pointer select (cs code path).
  If the A100 12.6/12.8 build shows wrong numbers on S1/S2, these are the first suspects.
- CUDA compile on Orion (cuda/12.1, sm_80, nvcc_wrapper, `--fmad=false`, NMAX 4096), state at
  hand-off: rad_m1_clean_solver.cpp (mg, fused ghost kernel), rad_m1_clean_cell.cpp,
  rad_m1_clean_closure.cpp, rad_m1_time2.cpp and the rest of rad_m1 compiled with 0 errors
  (tree of 841d51a76; the last three commits touch only rad_m1_clean.hpp/.cpp comments, the
  Marshak boundary share and the dt_retry visitor). **rad_m1_clean.cpp (most of the new
  staggered kernels) was still in cicc after 45 min**: please report its compile time on Raven
  (cuda 12.6) and MI300A; a compile that long is itself worth knowing. Orion will post the
  result in a NOTE.

## Rules
- Fresh starts only. Stop on FATAL, NaN or rc != 0 in S1-S3 and send the last 40 log lines;
  do not adapt decks or keys. KEEPGOING=1 for T.

## What to send back
1. Build: modules, Kokkos configure summary, `--fmad=false` in flags.make (CUDA), errors and
   warnings in the clean files (`*.clean_diag.txt`), binary md5s, for A100 and MI300A.
2. gates.md for A100 and MI300A.
3. timing.txt (+ kernel tables, memory high-water marks) for both GPUs.
4. Run-dir paths, `md5sum runs/*/bin/*.bin runs/*/tab/*.tab | sort`, run.log of any failure.
Reply: NOTE-2026-11-17-viper-orion-stag2-gpu1.md on fork/bsg-files-1009.
