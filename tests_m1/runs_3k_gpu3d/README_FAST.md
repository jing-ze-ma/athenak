# Speed-up #5: the GPU cost of the implicit-M1 Krylov iteration

Date 2026-09-23, viper apudev (MI300A). NOT COMMITTED.

- **Patch:** `bench/m1_fast_0923/fast.patch` (md5 89af182b). It applies cleanly (`git apply --check`) on
  HEAD 66d5c59a and on 54060691. It touches `src/rad_m1/rad_m1.{hpp,cpp}` and
  `src/rad_m1/rad_m1_implicit.cpp` only.
- **Development base:** a6202226. The patch was then carried onto 66d5c59a, and every final gate below
  uses 66d5c59a.
- **Snapshots:** the `/viper/u2` inode quota was full at the start (590 inodes left), so the
  snapshots and build trees went to local scratch on viper13. The final source tree is kept at
  `/viper/ptmp2/jinma/m1_fast_0923/`. The binaries are in `bench/m1_fast_0923/bin/`:
  - `base2` = HEAD 66d5c59a (GPU md5 07e6a5c7, CPU 1a724602);
  - `fin` = 66d5c59a + fast.patch (GPU 508529bc, CPU 16a849d9);
  - `base`, `new1` … `new5` = the development builds on a6202226.
- **Scripts:** jobs `time.sh`, `fin.sh`, `prof.sh`; CPU gates `cpu_gate*.sh`; tools in `tools/`
  (`split.py` = profile split, `cmp.py` / `stats.py` = gates, `summarize.py` = timing).
- **Environment:** every GPU job exports `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.

## Baseline

Baseline = **HEAD 66d5c59a defaults + `implicit_predictor = step`**. That is pcr, bcg_sync 1,
lres_test false, conv_est true, lin_ew_max 1e-2, which is the "recommended" configuration of
README_PICARD. `inp/box3d_gpu.athinput` spells out the pre-0923 values, so every arm restates
these defaults on the command line (`REC` in the scripts).

- Box: 84x104x104, 4 MeshBlocks.
- Timings: 120 cycles, ms/cycle measured from cycle 20 to cycle 120.
- Every timing comparison uses the same binary, with the arms interleaved and repeat b run in
  reverse order.

## 1. Profile of the baseline (Step 1)

rocprofv3 kernel trace, 40 cycles, window = cycles 5-40. Tool: `tools/split.py`. Each kernel's
category also carries the idle gap in front of it.

Runs:
- Eddington: `runs/P_base`, job 11945096, a6202226 + REC;
- vet_sc full: `runs/P_base_V`, job 11945535.

Profiled cycle times:
- Eddington: 131.2 ms/cycle, 1660 launches, idle 12.3 ms;
- vet_sc: 129.1 ms/cycle, 1741 launches, idle 10.2 ms.

Unprofiled, the Eddington baseline runs at 113-115 ms/cycle.

| ms/cycle (kernel + idle before) | Eddington | vet_sc full |
|---|---|---|
| Krylov halo: `ImplicitHaloCopy` 109-118 us x 280-300 + bvals pack/unpack | **42.4** | **42.2** |
| Krylov operator: `ImplicitOffDiagOp` 189 us (E), 91 us (V) + `ApplyOp` | **30.6** | 15.4 |
| Krylov vectors + reductions (BiCGStab kernels 33-59 us) | 20.9 | 16.9 |
| Preconditioner (PCR 47.9 us x ~130 + copies) | 7.6 | 7.0 |
| SC sweep (`VetSweepMB` 172 us x 85) | - | 18.5 |
| Picard level (assembly, EOS/gas Newton, opacity, closure) | 7.1 | 6.6 |
| hydro + conduction weights + BCs + sources + C2P | 22.7 | 22.6 |

Findings:
- `ImplicitHaloCopy` is 77 % of the halo cost. It loops over (m,n,k,j) rows with i serial
  inside, so the 3-D box has ~13k threads walking strided rows.
- Under the **Eddington closure** `ImplicitOffDiagOp` spends 27 ms/cycle computing terms that are
  identically zero, since D_ab = (3 chi - 1)/2 n_a n_b = 0 at chi = 1/3.

## 2. Levers (all behind `<rad_m1>` switches, all default off)

| key | what it does | exact? |
|---|---|---|
| `implicit_halo_direct = true` | The implicit exchanges become ONE on-rank ghost-copy kernel (neighbour table per block and direction). It is used only if every rank has all neighbours on-rank at the same level, with no seam or pole. Otherwise the ordinary exchange runs, with the shell copy done one thread per cell. | bitwise |
| `implicit_od_cache = true` | Computes `M1OffDiv` once per cell (3 values) instead of 12 times per cell, and fuses the 7-point row with the off-diagonal faces into one kernel. | bitwise |
| `implicit_krylov_fuse = 1` | The PCR preconditioner reads r and writes z itself, and builds the p / s update in its load phase (1 launch instead of 3). | bitwise |
| `implicit_krylov_fuse = 2 / 3` | 2: the dot products ride in the operator kernel. 3: `ImplicitBiCGStabTwo`, with **2 blocking reductions per iteration**; the convergence test is taken half an iteration late, and rho comes from a recurrence. | round-off |
| `implicit_precond = rbgs / rbgs_fwd` | Symmetric or forward red-black **transverse line Gauss-Seidel**: PCR line solves per colour, block-local, no communication. | round-off (different iterates, same converged answer) |
| `implicit_op_stencil = true` | Writes the frozen operator of each Picard pass once, as a 19-point stencil (7-point row + Eddington off-diagonal terms). It uses the 7 face/centre slots only when all edge coefficients are 0. Reductions use 256-thread blocks: the default takes 33 kB of LDS, which allows one block per CU. | round-off |
| `implicit_precond_float = true` | Runs the preconditioner's line solves in float (MIXED PRECISION; the operator, the vectors and the reductions stay double). | round-off |

**Gains:** 1 GPU, Eddington, ms/cycle for repeats a and b, cycles 20-120. NON-CONVERGED 0 in
every arm.

| arm (job, binary) | ms/cycle | Picard/step | inner its/step |
|---|---|---|---|
| baseline (11945172, new2 = HEAD path) | 113.4, 114.6 | 2.858 | 62.88 |
| + halo_direct | **80.4, 79.7** (-34) | 2.858 | 62.88 |
| + od_cache alone | 104.7, 103.5 (-10) | 2.858 | 62.88 |
| + krylov_fuse 1 alone | 111.9, 111.9 (-2) | 2.858 | 62.88 |
| halo + od_cache + fuse 2 | 69.6, 71.1 | 2.858 | 63.13 |
| halo + od_cache + fuse 3 (2 syncs/its) | 65.1, 65.4 (-5 vs fuse 2) | 2.858 | 63.33 |
| ... + rbgs (symmetric) | 55.7, 57.5 | 2.158 | 39.97 |
| ... + rbgs_fwd | 54.8 (a) | 2.233 | 41.50 |
| ... + rbgs_fwd + op_stencil (11945309, new5) | **48.7, 48.5** (RBF alone: 53.3, 54.3) | 2.233 | 41.49 |
| ... + rbgs_fwd + op_stencil + precond_float (11945266, new4) | 47.0, 49.8 (vs 47.6, 48.6 without) | 2.233 | 41.78 |

- **Preconditioner:** rbgs cuts inner iterations by 34 % and Picard passes by 22 %, because
  Eisenstat-Walker needs fewer passes.
- **Mixed precision** is neutral within noise on 1 GPU (above) and gives -1 ms on 2 GPUs (C2, job
  11945267: 44.9, 45.4 vs 46.1, 46.2). The line solve is latency- and barrier-bound, not
  bandwidth-bound. **Not recommended.**
- **Launch bounds and the edge-skip** (new4 → new5, both binaries in job 11945309) are neutral:
  48.7, 48.8 vs 48.7, 48.5.

**RECOMMENDED (`FAST`):** `implicit_halo_direct = true`, `implicit_od_cache = true`,
`implicit_krylov_fuse = 3`, `implicit_op_stencil = true`, `implicit_precond = rbgs_fwd`.

## 3. Gates, on HEAD 66d5c59a (`base2`) vs 66d5c59a + fast.patch (`fin`)

**Default path, byte-identical to HEAD.**
- GPU: `tests_m1/runs_3k_gpu3d/he_slab_m1_3d.athinput` untouched (pure HEAD defaults), 120
  cycles. hst identical, BITWISE:
  - 1 GPU: `F1_dB` vs `F1_dF`, job 11945578;
  - 2 GPUs: `F2_*`, job 11945580.
- CPU: seeded 2-D slab, 200 s (`cpu/f_def_*`, `f_vdef_*`):
  - Eddington and vet_sc: hst identical;
  - all 6 bin payloads identical (`tools/pmd5.py`).
- CPU with REC and the new keys off (`f_ref` vs `f_ref_n`): hst identical.
- Development: halo_direct, od_cache and krylov_fuse 1 are bitwise equal to the baseline on the
  GPU (job 11945172) and on the CPU. The CPU runs covered 1 rank, 1 rank with 2 blocks (the direct
  kernel) and 2 ranks (the flat-copy fallback) (`cpu/` of gate 2).

**Converged answer, CPU 2-D slab, 200 s** (`cpu_gate3.sh`, relative to base2 + REC).
- `ctl` = base2 with `rad_flux_inner` changed by ±1 ulp; `tol10` = implicit_tol 1e-10.
- NON-CONVERGED is 0 in all 19 runs.

| run | F1top/Fin mean t>100 | KE1 end | KE2 end | dt mean | totE end |
|---|---|---|---|---|---|
| Eddington ref (1.000002, 1.0852e26, 1.728e23, 0.16120, 8.0476e32) | | | | | |
| ctl +1 / -1 ulp | 5.5e-13 / 6.3e-13 | 2.0e-9 / 5.8e-9 | 2.6e-8 / 5.3e-9 | 4e-12 | 5.9e-13 / 8.5e-13 |
| tol10 | 3.1e-13 | 7.7e-9 | 1.7e-8 | 3.8e-12 | 1.4e-12 |
| **FAST** | 5.4e-13 | 8.4e-9 | 8.1e-9 | 1.2e-12 | 1.9e-12 |
| FAST, 2 blocks on 1 rank | 4.8e-13 | 4.0e-9 | 3.1e-9 | 3.2e-12 | 7.6e-14 |
| FAST, 2 ranks | 1.1e-12 | 5.6e-9 | 5.8e-9 | 3.0e-12 | 1.2e-12 |
| vet_sc full ref (1.000004, 1.7783e26, 1.735e23, 0.16053) | | | | | |
| ctl +1 ulp | 5.5e-12 | 6.4e-9 | 1.5e-8 | 1.7e-11 | 2.3e-12 |
| **FAST** | 3.7e-12 | 1.6e-8 | 5.5e-8 | 6.6e-12 | 2.0e-12 |
| FAST, 2 ranks | 5.9e-12 | 9.5e-9 | 4.3e-8 | 6.8e-11 | 3.0e-12 |

- Inner iterations, CPU slab: Eddington 42 216 → 25 611, vet 39 717 → 25 621. CPU wall:
  146 → 110 s and 151 → 125 s.
- **Verdict: PASS.** Every change is at the ±1-ulp / tol10 level, on 1 and 2 ranks.

**GPU 3-D box, 120 cycles (t = 19.36), FAST vs base2.**
- Eddington, 1 GPU: totE 5.3e-12, KE1 1.2e-9, F1top 2.4e-12; the ±1-ulp control gives 5.8e-12,
  3.2e-10 and 5.3e-12.
- vet_sc, 1 GPU: totE 4.9e-13, KE1 4.6e-10; the control gives 8.3e-13 and 2.3e-10.
- 2 GPUs: the same level.
- The largest normalised hst column differences (up to 0.55) are all in the 2-mom and 3-mom
  columns. Those are net transverse momenta of about 4.7e9, against 1.4e21 for 1-mom, i.e.
  round-off noise at 3e-12 of the momentum scale.

**Final GPU timing** (jobs 11945578-81, fin and base2 interleaved, ms/cycle for repeats a, b):

| | baseline (HEAD + predictor) | FAST | speed-up |
|---|---|---|---|
| Eddington, 1 GPU | 112.9, 114.3 | **48.3, 48.4** | 2.35x |
| Eddington, 2 GPUs | 94.3, 92.8 | **47.9, 47.0** | 1.97x |
| vet_sc full, 1 GPU | 119.8, 120.0 | **67.1, 67.4** | 1.78x |
| vet_sc full, 2 GPUs | 106.5, 106.3 | **69.7, 69.8** | 1.53x |
| (HEAD pure defaults, no predictor, E, 1 GPU) | 129.3, 127.5 | - | - |

## 4. Where FAST's cycle goes now, and the gap to hydro

Profiles:
- Eddington: `runs/P_new5_STF`, job 11945490. The source is identical to fast.patch before the
  rebase.
- vet_sc: `runs/P_new5_VSTF`, job 11945534.

| ms/cycle (profiled) | Eddington: 50.9 total, 588 launches | vet_sc: 67.2 total |
|---|---|---|
| operator (`ImplicitStencilOp` 87 us x 92 + build 1.3) | 9.8 | 9.3 |
| preconditioner (PCR half-sweeps 35.6 us x 183) | 8.2 | 7.5 |
| halo (`ImplicitHaloDirect` 11.5 us x 102) | 1.3 | 1.3 |
| vectors + reductions (incl. sync gaps) | 2.7 | 2.3 |
| SC sweep | - | 18.3 |
| Picard level | 6.1 | 6.1 |
| hydro + rest | 22.8 | 22.4 |

- **Eddington:** radiation is ~28 ms/cycle against hydro's 22.8, i.e. 1.2x hydro. At baseline it
  was 108 ms, 4.8x.
- **vet_sc full:** radiation is ~45 ms/cycle, 2.0x hydro. 18 ms of that is the SC formal
  solution.
- **2 GPUs** are no faster than 1 (47.5 vs 48.3 ms). 454k cells per GPU are latency-bound: 588
  launches per cycle and 2 host syncs per Krylov iteration.

## 5. Ranked TO-DO for the next session (numbers from the profiles above)

1. **SC sweep (vet_sc):** `VetSweepMB` costs 171 us x 85 per cycle = 18.3 ms, 27 % of the FAST vet
   cycle. This is the other agent's file (rad_m1_vet.cpp); it is the largest remaining item for
   vet_sc.
2. **Inert conduction in the M1 input (the "hydro" 22.8 ms):** about 6.3 ms/cycle of work whose
   weight is 0 on every face (`rad_tau_lo = 1e5`):
   - `Conduction::BuildRadWeights` 2.57 ms x 2;
   - `AddIsotropicHeatFluxRadiative` 3 x 0.26 ms;
   - `Conduction::NewTimeStep` 0.35 ms.

   Skipping it when M1 owns rad_flux_inner would cut the M1 run cycle by ~13 %.
3. **Reduction kernels:** the fused operator + reduction runs at 85-87 us, against 46 us for the
   same stencil as a par_for. A plain dot product takes 42 us, against 17 us for a 5-vector
   update. The operator reductions alone cost (87 - 46) us x 92 = 3.7 ms/cycle. Try a hand-rolled two-stage (team-partial) reduction.
4. **Iterations** (41.5 inner its and 2.23 passes per step; the preconditioner costs 8.2 ms):
   - Mode analysis: the slow modes are horizontally smooth. Line Jacobi gives λ ≈ (s + a2 θ²)/(s + 4 a2),
     with a2 ≈ 15 in the thick layer.
   - A semicoarsening (x2-x3) multigrid V-cycle with rbgs smoothing is the robust cure. Its
     estimated net gain is only ≤ 4 ms at the current ~0.3 ms per iteration, because it adds
     ~15 launches per application.
5. **Picard level (6.1 ms):**
   - `Opacity` 0.90 ms plus `ImplicitSolve#8` 0.74 ms (inferred to be the opacity pass
     `m1_impl_opac`), once per step;
   - stencil build 1.3 ms;
   - EOS / assembly kernels at 0.2-0.5 ms each.
6. **Host syncs:** 2 per iteration, ~3.5 ms/cycle of idle before the PCR and vector kernels.
   - A pipelined BiCGStab (1 sync) would save ~1.7 ms.
   - On 1 rank, device-side scalars plus a convergence check every N iterations would remove
     nearly all of it.
7. **Multi-GPU:** there is no scaling from 1 to 2 GPUs on this box. Production should load
   ≥ 1 M cells per GPU, and the timings should be re-measured there.

Not done:
- **Restricting the implicit solve by optical depth** (candidate d). It is the thick layer where
  horizontal diffusion over one step spans several cells (a2 ≈ 15), and that coupling is exactly
  what limits the Krylov convergence. So excluding cells does not remove the hard part.
- **Float operator / vectors.** They would move the converged answer by ~cond x 1e-7 unless the
  true residual is refined in double.
- **HIP graphs.** They need device-resident scalars first (item 6).
