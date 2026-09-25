# runs_5p_coarse2: a global coarse correction for the implicit M1 preconditioner

- **Date:** 2026-09-25, viper. **Branch:** `m1-coarse2` (from rt-integration 038148c6),
  worktree `/viper/ptmp2/jinma/wt_coarse2`.  Run tree `/viper/ptmp2/jinma/coarse2_0925`
  (binaries `bin/`: `base_*` = 038148c6, `g4_*` / `c4_*` = final GPU / CPU; final timing
  `runs/t4` job 11977352, GPU gates `runs/bw` + `runs/ts4` job 11977353, profiles
  `runs/prof4`; operator dumps `runs/d`; CPU gates `gates/`, `gates_eval.txt`,
  `refgates_eval.txt`, `rstgc.log`, `tst_rad_m1.log`).
- The brief's Part 1 (scaling baseline at 1-8 GPUs) was dropped by the user (apu queue);
  1 and 2 GPUs (apudev) only.  4-8 GPU runs: `job_scale.sh` (commands below), not run.

## 1. What was built: `<rad_m1>/implicit_precond = mg_gc`

Two-level preconditioner, multiplicative, coarse FIRST:

    x_g = (P^T A P)^{-1} P^T r;   r' = r - A P x_g;   z = P x_g + M_mg r'

- P: piecewise constant over `implicit_gc_bands2 x implicit_gc_bands3` bands of the
  WHOLE (x2,x3) mesh (default 1 x 1 = the sideways mean of each x1 layer), one unknown
  per band and global x1 index.  M_mg = the merged `mg` with `implicit_mg_levels`
  (1 = rbgs_fwd alone, then no halo or coarse levels).
- A is the Krylov operator itself (the stored stencil, edges, folded vimp).  A first
  version built P^T A P from the 7-point preconditioner rows (TA..CKP) and was WORSE than
  no coarse space with velocities (CPU gate box: mg_gc 3 39.8 it vs mg 3 12.3): the
  vimp terms are missing from those rows, and a mean-mode error in them is amplified.
- Coarse first: P^T r needs no operator product; A P x_g needs no halo (x_g is global).
  Coarse-first and coarse-last have the same error-propagation spectrum.
- Per application: 1 band-reduction kernel (it also makes the BiCGStab p / s update),
  1 host round trip with ONE `MPI_Allreduce` of N1 x nbands numbers, the banded LU solve
  (half bandwidth 3 nbands - 1) on the host of every rank, 1 fine kernel for r' (5 row
  sums per cell for one band), and P x_g added in the existing prolongation kernel.
  The coarse matrix: 1 per-cell kernel + 1 band reduction + 1 allreduce + LU, once per
  Picard pass.  Band sums are deterministic (fixed-order team sums, host sum in a fixed
  order); 1 vs 2 ranks differ by round-off only.
- Needs a single-level mesh and `implicit_op_stencil` (fatal otherwise).  The rest of
  the code: 5 lines of dispatch in rad_m1_implicit.cpp; everything else in
  rad_m1_precond.cpp.
- Also added, opt-in (read only when named, default off = bitwise):
  `implicit_bcg_rho_direct` sums (rhat, r) directly in the update kernel of the
  krylov_fuse = 3 BiCGStab.  The (rhat,r) recurrence stalls the 1-2 hardest solves of
  the CPU gate box (pass 0 of cycles 1-2 run to lin_maxit = 200, with rbgs_fwd AND
  mg_gc; mg avoids it; 0 NON-CONVERGED, Picard absorbs it); rho_direct cures those
  (rbgs_fwd 26.9 vs 43.6 it/solve there).  It costs one blocking reduction per
  iteration (GPU: box rbgs_fwd +3.5 %, wedge +10 %, job 11972909), and the GPU
  production cases never stall, so it stays off.

## 2. Offline (study2.py on cycle-20 dumps, BiCGStab to 1e-10 from x0, runs/d)

| preconditioner | box vet_sc | wedge vet_col |
|---|---|---|
| rbgs_fwd | 22 | 15 |
| rbgs_fwd + global mean (mg_gc levels 1) | 21 | 9 |
| mg 3 (block-local) | 9 | 7 |
| global mean + mg 3 (mg_gc 3) | 7 | 4 |
| 4 x 4 bands + mg 3 | 9 | 4 |
| full global V-cycle | 3 | 2 |

The box's leftover after the mean is k = 1-2 sideways error, which piecewise-constant
bands do not represent (bands made it worse on the box).  The wedge's slow error IS
the sideways mean.

## 3. GPU (apudev, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, binary g4, interleaved, 2 reps,
runs/t4, job 11977352; box 2 GPUs = weak, nx2 208; wedge 2 GPUs = strong, same grid)

| case | rbgs_fwd | mg 3 | mg_gc 1 | mg_gc 2 | mg_gc 3 |
|---|---|---|---|---|---|
| box 1 GPU ms/cycle | 50.99, 51.03 | **47.83, 47.87** | 54.35, 54.39 | 50.45, 50.39 | 48.15, 48.21 |
| box 1 GPU it/solve (max) | 21.2 (49) | 9.8 (22) | 17.0 (43) | 10.7 (26) | 7.9 (19) |
| box ms/iteration (slope 1e-10 vs 1e-8) | 0.54 | 0.97 | 0.84 | 1.11 | 1.28 |
| box 2 GPUs ms/cycle | 57.26, 57.44 | **54.66, 54.74** | 62.20, 61.96 | 58.55, 58.41 | 55.68, 55.91 |
| box weak eff. 1 -> 2 GPUs | 0.89 | 0.88 | 0.88 | 0.86 | 0.86 |
| wedge 1 GPU ms/cycle | 37.62, 37.77 | 38.28, 38.23 | **30.79, 30.88** | 31.03, 31.02 | 31.46, 31.31 |
| wedge it/solve (max) | 13.1 (21) | 7.3 (10) | 1.1 (2) | 1.0 (1) | 1.0 (1) |
| wedge ms/iteration | 0.73 | 1.39 | ~1.8 (*) | (*) | (*) |
| wedge 2 GPUs ms/cycle | 24.20, 24.20 | 24.75, 24.98 | **18.96, 18.80** | 19.12, 18.92 | 19.37, 19.25 |
| wedge strong eff. 1 -> 2 GPUs | 0.78 | 0.77 | 0.82 | 0.81 | 0.81 |

(*) 1.10 vs 0.81 it/step: the slope is not measurable; the kernel trace (runs/prof4)
gives ~0.1 ms of gc work per application (GCPre 28 us/launch + 25 us GPU idle for the
host round trip).  0 NON-CONVERGED in all 50 runs.

- **Wedge: mg_gc 1 = -18 % per cycle on 1 GPU (30.8 vs 37.7 ms) and -22 % on 2 GPUs
  (18.9 vs 24.2)**, 1.1 it/solve vs 13.1: the solve almost disappears; the cycle is
  now VetColBuild (6.7 ms), the Picard/closure kernels and the halo.
- **Box: mg_gc 3 = parity with mg 3** (+0.7 % on 1 GPU, +2 % on 2 GPUs) at 7.9 vs 9.8
  it/solve: the coarse space saves 1.9 it/solve and costs ~0.3 ms per iteration
  (box profile: GCPre 19 us/launch, 2 per application, + 1.4 ms/cycle GPU idle for the
  host round trips).  mg stays the box default.
- Cost history (box 1 GPU mg_gc 3): 52.4 ms (per-thread 45-slot coarse build,
  1.4 ms/launch; serial band sums, 52 us) -> 51.3 (per-cell one-band build) -> 48.2
  (row-chunked team band reductions, 19 us).

## 4. Gates (all on the final binaries c4 / g4)

- **Default bitwise:** CPU gates.py arms box_a1/m2, slab_a1/m2, loose and tight, base
  038148c6 vs new: 8/8 rst bitwise, hst identical (gates_eval.txt).  GPU: slab 1 GPU,
  box and box vet_sc 2 GPUs 9/9 files bitwise each (runs/bw); the wedge case has no
  output files (0/0).  tst/test_suite/rad_m1: 3 passed (tst_rad_m1.log).
- **Same fixed point** (gates_gc.py, gates.py criteria, lin_tol 1e-11 tight arm):
  mg_gc 1 / 3 / 2x2 bands vs rbgs_fwd: box hst_cons 1.5-1.6e-10, slab 3.1e-9,
  rst_cons 2-6e-10; 1 vs 2 ranks PASS in every variant.  The strict gate's F1top miss is
  the level every other valid preconditioner shows with the same binary
  (refgates_eval.txt): line vs rbgs_fwd 1.8e-10 (box) / 1.8e-9 (slab), rbgs 1.8e-9
  (slab), mg 8.3e-11 / 2.6e-9.
- **Restart** (rst_gc.py: CI restart test + mg_gc 3, 2 bands, rho_direct, 2 ranks,
  Eddington and vet_sc): 16 cycles vs 8 + restart + 8 BITWISE.
- **0 NON-CONVERGED** in every CPU gate arm and GPU run.
- **Wedge T-S4** (sph_atm_vc 800 cycles, runs/ts4): mg_gc 3 vs default dE/E 5.2e-10,
  dF/|F| 1.5e-9, transverse |F2,F3|/|F1| 4.9e-10 (1 GPU, 4x4); 2x2 bands 5.9e-10 /
  1.5e-9; 2 GPUs (8x8) 6.6e-10 / 3.0e-9: unchanged to solver tolerance.

## 5. Verdict

- sp wedge (vet_col): `implicit_precond = mg_gc`, `implicit_mg_levels = 1` (-18 % /
  -22 % per cycle at 1 / 2 GPUs, better strong scaling 0.82 vs 0.78).
- Cartesian box (vet_sc): keep `implicit_precond = mg`, `implicit_mg_levels = 3`;
  mg_gc 3 only matches it.  The box needs the k = 1-2 sideways modes (a global V-cycle
  reaches 3 it/solve offline), i.e. coarse levels across blocks, not band means.
- The default (rbgs_fwd) is unchanged and bitwise.

## 6. 4-8 GPU scaling (not run; scripts ready)

    cd /viper/ptmp2/jinma/coarse2_0925
    sbatch -p apudev -N 1 --ntasks-per-node=1 --gres=gpu:1 scripts/job_scale.sh 1
    sbatch -p apudev -N 1 scripts/job_scale.sh 2
    sbatch -p apu -N 2 scripts/job_scale.sh 4
    sbatch -p apu -N 4 scripts/job_scale.sh 8
    python3 scripts/tsum.py runs/s<NP> 10 30

(strong: box 84x208x208 and wedge 96x128x128, 16 blocks; weak: box 4 and wedge 16
blocks per GPU; rbgs_fwd / mg 3 / mg_gc 1 / mg_gc 3 interleaved, 2 reps.)

## Files

`study2.py` (offline, needs runs_5m_precond/study.py), `gates_gc.py` (CPU gates),
`refgates.py` (reference pairs, same binary), `rst_gc.py` (restart gate),
`cpu_gates.sh`, `job_t4.sh` (timing + profiles), `job_g.sh` (GPU gates),
`job_scale.sh`, `build.sh`, `tsum.py`, `ksum.py`.
