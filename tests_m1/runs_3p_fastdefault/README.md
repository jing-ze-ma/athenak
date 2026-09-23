# Implicit-M1 fast path and vet_mb_halo = 3 become the defaults

Date 2026-09-23, viper. Base HEAD 5cd0a86d (predictor default + restart state, 66f1f889).
Build and run tree: `/viper/ptmp2/jinma/m1_fastdef_0923/` (binaries deleted afterwards;
inputs, scripts, logs and hst files kept). Patch = `/viper/ptmp2/jinma/fastdef.patch`.

## What changed

`src/rad_m1/rad_m1_implicit.cpp`:
- The `FAST` set of `tests_m1/runs_3k_gpu3d/README_FAST.md` is now the default:
  `implicit_halo_direct = true`, `implicit_od_cache = true`, `implicit_krylov_fuse = 3`,
  `implicit_op_stencil = true`, `implicit_precond = rbgs_fwd`.
- These defaults apply only when all of the following hold:
  - closure = eddington, vet_sc or tau. Closure m1 (Levermore), minerbo and kershaw keep the old
    path, as 66d5c59a did; the fast path was never gated for them;
  - transport = implicit with implicit_solver = bicgstab;
  - implicit_bcg_sync = 1 and implicit_line_solver = pcr;
  - one MeshBlock along x1;
  - no implicit_lin_cnorm > 0.
- `implicit_op_stencil` is also off by default under a periodic x1 wrap, where it would be fatal.
- Any key named in the input keeps its value. So would a restart written before this commit,
  because GetOrAdd echoed the old `false` values into it.

`src/rad_m1/rad_m1_vet.cpp` (the default value only):
- `vet_mb_halo` defaults to 3 with the ray kernel and no vet_mb_lag.
- If the band does not fit a MeshBlock, it is narrowed to the largest K that fits.
- It is still read only when given, so the restart echo does not change.

## Gates

**CPU** (`cpu_gate.sh`, gcc/openmpi build of the same source).
- `*_old` = the same binary with the old values written explicitly into the input
  (`inp/*_off.athinput`). `*_new` = the input without those keys (`inp/*_def.athinput`).
- NON-CONVERGED = 0 in all 27 runs.

**Explicit old values give bitwise HEAD.** 2-D seeded slab, 200 s: HEAD 5cd0a86d against the
new binary run with the old values (`e_head`/`e_old`, `v_head`/`v_old`). The hst files are
identical.

**2-D seeded slab, 200 s:** relative differences against `*_old`. `ctl` = old values with
rad_flux_inner changed by +1 ulp.

| arm | F1top/Fin mean | KE1 end | KE2 end | dt mean | totE end |
|---|---|---|---|---|---|
| Eddington ctl | 5.5e-13 | 2.0e-9 | 2.6e-8 | 4.0e-12 | 6.0e-13 |
| Eddington new | 5.4e-13 | 8.4e-9 | 8.1e-9 | 1.2e-12 | 1.9e-12 |
| Eddington new, 2 blocks, 1 rank | 4.8e-13 | 4.0e-9 | 3.1e-9 | 3.2e-12 | 7.6e-14 |
| Eddington new, 2 ranks | 1.1e-12 | 5.6e-9 | 5.8e-9 | 3.0e-12 | 1.2e-12 |
| vet_sc full ctl | 5.5e-12 | 6.4e-9 | 1.5e-8 | 1.7e-11 | 2.3e-12 |
| vet_sc full new | 3.7e-12 | 1.6e-8 | 5.5e-8 | 6.6e-12 | 2.0e-12 |
| vet_sc new, 2 blocks, 1 rank | 5.7e-12 | 1.2e-8 | 5.7e-8 | 7.7e-11 | 2.7e-12 |
| vet_sc new, 2 ranks | 5.9e-12 | 9.5e-9 | 4.3e-8 | 6.8e-11 | 3.0e-12 |
| tau new | 8.6e-13 | 3.0e-9 | 1.3e-8 | 6.6e-13 | 5.2e-12 |

Inner iterations, CPU slab:
- Eddington 42 216 → 25 611;
- vet_sc 39 717 → 25 621;
- tau 40 711 → 25 803.

**3-D box on CPU:** 84x32x32 in 4 blocks of 84x16x16, 60 cycles, t = 9.68. Relative
differences at the end:

| arm | totE | KE1 | F1top |
|---|---|---|---|
| Eddington ctl | 1.9e-13 | 1.2e-9 | 7.7e-12 |
| Eddington new | 5.8e-12 | 5.6e-10 | 4.8e-12 |
| Eddington new, 2 ranks vs new, 1 rank | 4.7e-13 | 2.2e-10 | 2.5e-12 |
| vet_sc full ctl | 4.9e-12 | 4.7e-10 | 1.5e-11 |
| vet_sc full new | 2.8e-12 | 4.8e-10 | 1.0e-11 |
| vet_sc new, 2 ranks vs new, 1 rank | 2.0e-12 | 1.2e-9 | 8.5e-12 |

**vet_mb_halo 3 vs 1 is bitwise.** The hst files are identical in three cases:
- 2-D, 2 ranks (`v_m2`/`v_m2h1`);
- 3-D, 2 ranks (`v3_m2`/`v3_m2h1`);
- GPU, 1 and 2 GPUs (`T*_Vnew_a`/`T*_Vh1_a`).

**Restart bitwise, with the new defaults (predictor + fast path + halo 3):**
- Method: run to t = 40 with a restart file at t = 20, then restart from t = 20.
- CPU, 2-D slab:
  - Eddington, 1 block (`re_*`);
  - vet_sc full, 2 blocks on 1 rank (`rv_*`).
- GPU job 11945786, 3-D box: t = 8 with the restart file at t = 4, Eddington and vet_sc, on
  1 and on 2 GPUs.
- In every case these are identical: the final rst file, both final bin files, and every hst
  row after the restart (`tools/hstrst.py`).
- NON-CONVERGED = 0.

**Verdict: PASS.** Every change is at the ±1-ulp control level, on 1 and 2 ranks and in 2-D and
3-D. Halo 3 and restarts are bitwise.

## GPU timing

Jobs 11945784 (1 GPU) and 11945785 (2 GPUs), apudev:
- one binary (md5 0df0b801);
- arms interleaved, with repeat b in reverse order;
- `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.

Setup: 3-D He box 84x104x104 in 4 blocks, `inp/box3d_def.athinput` (= `he_slab_m1_3d`), 120
cycles, ms/cycle from cycle 20 to 120 (`tools/summarize.py`). Old = the old values written in the
input (`box3d_off`). NON-CONVERGED = 0 in every arm.

| ms/cycle (a, b) | old defaults | new defaults | speed-up | new, vet_mb_halo = 1 |
|---|---|---|---|---|
| Eddington, 1 GPU | 113.9, 114.7 | **48.9, 49.4** | 2.35x | - |
| Eddington, 2 GPUs | 93.6, 93.5 | **47.2, 47.2** | 1.98x | - |
| vet_sc full, 1 GPU | 103.0, 103.0 | **53.3, 55.1** | 1.90x | 53.2, 53.0 |
| vet_sc full, 2 GPUs | 90.6, 90.3 | **51.9, 51.4** | 1.75x | 54.1, 53.8 |

- Picard passes per step:
  - Eddington 2.858 → 2.233;
  - vet_sc 2.983 → 2.733.
- Inner iterations per step:
  - Eddington 62.9 → 41.5;
  - vet_sc 59.4 → 40.0.
- **vet_sc with the fast path and dd4df74c together** runs at 51-55 ms/cycle. The hydro-only
  cycle is 22.8 ms.
- **Halo 3** saves 2.3 ms/cycle on 2 GPUs. On 1 GPU it is neutral within noise.
