# runs_5n_mgfuse: fusing the mg coarse-level kernels (phase 1) -- negative result

- **Date:** 2026-09-25, viper. **Branch:** `m1-mgfuse` (rt-integration 9de6060b + the
  fusion 2415449f + its revert).  Run tree `/viper/ptmp2/jinma/mgf_0925` (profiles
  `runs/prof_ref` = merged mg 091d1422 code, `runs/prof_v1` = 2415449f; timing `runs/t`,
  job 11970346; CPU gates `gates/`, `gates_eval.txt`, `rstmg.log`, `tst_rad_m1.log`).
- **Verdict: no code change.** The fusion (the restrictions of level l formed in the
  load phase of level l's own line sweeps, `implicit_mg_fuse`, 2415449f) is exact
  (bitwise the unfused map) but SLOWER, and adding it slowed the unfused kernel too, so
  the source is reverted to rt-integration.  `implicit_precond = mg` stays as merged.

## 1. Profile of the merged mg (rocprofv3 kernel trace, 1 GPU, cycles 10-40, 3 levels)

| | box vet_sc | wedge vet_col |
|---|---|---|
| mg applications per cycle | 37.9 | 25.9 |
| M1PCRMG (4 per application) | 1.48 ms, 9.8 us/launch | 1.32 ms, 12.7 us |
| res0 (fine residual + restriction) | 1.00 ms, 26 us | 1.17 ms, 45 us |
| res1 (level 1 -> 2) | 0.67 ms, 18 us | 0.46 ms, 18 us |
| prolongation (one kernel) | 0.48 ms, 13 us | 0.45 ms, 17 us |
| **mg kernels total / GPU idle before them** | **3.71 / 0.27 ms/cycle** | **3.50 / 0.16** |
| + one extra ImplicitHaloDirect per application | ~15 us | ~24 us |

**The launch gaps are not the cost**: the GPU idles only 0.16-0.27 ms/cycle in front
of all 186-270 mg launches.  The cost is the execution of small latency-bound kernels
(a 7-round PCR per coarse column, a residual that reads 8 fine components), about
100 us of kernel time per application, i.e. ~0.2 ms per BiCGStab iteration of GPU
work plus the extra halo.

## 2. The fused variant (2415449f), same binary, interleaved, 2 repeats (job 11970346)

| case | rbgs_fwd | mg unfused | mg fused |
|---|---|---|---|
| box 1 GPU ms/cycle | 52.69, 52.67 | **48.83, 48.86** | 49.94, 49.96 |
| box 2 GPUs | 60.42, 60.28 | **57.34, 57.30** | 57.75, 57.92 |
| wedge 1 GPU | 43.18, 43.40 | **42.45, 42.47** | 43.46, 43.14 |
| wedge 2 GPUs | **24.99, 24.89** | 25.98, 25.91 | 26.10, 26.12 |
| it/solve box 1 / 2 GPUs | 21.5 / 23.0 | 9.8 / 10.1 | 9.8 / 10.1 |
| it/solve wedge 1 / 2 GPUs | 13.0 / 13.1 | 7.3 / 7.3 | 7.3 / 7.3 |
| ms per iteration, box 1 GPU (slope lin_tol 1e-10 vs 1e-8) | 0.69 | 1.12 | 1.24 |
| ms per iteration, wedge 1 GPU | 1.27 | 2.04 | 2.17 |

0 NON-CONVERGED in all runs.  Kernel view (`ksum.py`): fused M1PCRMG 29.6 us/launch on
the box (vs 9.8 merged) and 41 us on the wedge; the mg group 5.03 ms/cycle fused vs
4.48 unfused-in-the-same-binary vs 3.71 merged.  The restriction inside the team's
load phase runs on the nx1 threads of one column team with scattered reads, far
slower than the flat par_for it replaced; and carrying the restriction views and
branches in the kernel slowed even the unfused sweep (15 us/launch vs 9.8).

## 3. Gates (all pass; for the record, on 2415449f)

- default (rbgs_fwd) bitwise: gates.py arms box_a1/m2, slab_a1/m2, base 9de6060b vs new;
- mg fused vs unfused: bitwise (rst and hst), 1 and 2 ranks, loose and tight;
- mg vs rbgs_fwd: as runs_5m_precond (F1top 8e-11 box / 2.6e-9 slab, the level of the
  line and rbgs reference pairs); restart (rst_mg.py) bitwise; tst rad_m1 3 passed.

## 4. Recommended next step (phase 2, NOT implemented): a global coarse space

The runs_5m_precond diagnosis: after rbgs_fwd, 98-99.7 % of the box error is the
sideways MEAN of each x1 layer, and a global V-cycle needs 2-3 iterations offline vs 9/6
for the block-local mg.  Cheapest global form:
- coarse space = piecewise-constant fields over gb2 x gb3 GLOBAL (x2,x3) bands, one
  unknown per band and global x1 index (1 x 1 band = the per-layer sideways mean; for
  the wedge, where k = 1-2 dominate at the base, a few theta/phi bands or 2-4 phi
  Fourier modes);
- restriction P^T r: per-row partial sums reduced in a fixed order on the device, one
  MPI_Allreduce of N1 x nbands numbers per application;
- coarse operator P^T A7 P (built once per Picard pass, same allreduce): block-
  tridiagonal in x1 (nb x nb blocks, diagonal x1 couplings), factorised and solved
  redundantly on every rank on the host (N1 nb^2 flops);
- multiplicative, coarse FIRST: x_g = A_g^{-1} P^T r, then the existing mg on
  r - A7 P x_g (computable without a halo, x_g is global), z += P x_g in the existing
  prolongation kernel.  Cost ~2 extra kernels + one host round trip per application.
An untested draft of exactly this is kept in `phase2_global_draft_UNTESTED.patch`
(it also contains the phase-1 fusion; keys implicit_mg_global, implicit_mg_gbands2/3).
Measure 1, 2 and 4-8 GPUs, both geometries, before adopting.

## Files

`ksum.py` (kernel time per cycle from a rocprofv3 trace), `job_t.sh` (timing),
`job_prof.sh` (profiles), `build.sh`, `cpu_gates.sh`, `gates_fuse.py`.
