# TASK for viper: speed defaults that change results for little or unproven gain (Caltech audit 09-27)

The user asked us to check whether any earlier speed-up changes results without gaining much speed. The Caltech
session ran a **read-only audit** at 1a19b375. It used memory notes, the handovers, docs/dev (default_flips_0927,
switch_inventory_2026-09-24, m1_wedge_0926, ke_dt_0926), the tests_m1 READMEs, merge messages and the GetOrAdd
defaults. **Nothing was changed.** The user wants viper to decide on each item. Most items need MI300A numbers,
and MI300A is the production machine.

Terms used below:
- "round-off" means GPU results differ at round-off level; the noise ratio is against the 1e-14 kick twins.
- "tolerance level" means results differ above round-off but stay below the solver tolerance.

## 1. Changes results, gain ~0 or negative: candidates to make opt-in

| switch | result effect | measured gain | suggested check |
|---|---|---|---|
| M1 `implicit_one_pass 8` (+ auto), where it is still the default (every setup except `force_reference = wb_arad`, which is 0 since m1-keydefault-0927) | 5-22x the round-off spread | ~0 (`sph_atm`); -1..0 % (wedge); auto-disable timed on H200 only | default 0 everywhere? |
| M1 `implicit_op_team_red` (default on since dbfbe0d6) | round-off | -12 % on the static T-S4 wedge (1 GPU); **-0.7 %** on the H200 box and on the moving MI300A sph_wedge (2 GPU); box -2 % / -4 % (1 / 2 GPU) | re-measure moving-gas box and wedge on MI300A; opt-in, or on only where it pays |
| M1 `implicit_halo_overlap` / `ovl_faces` (on with >1 rank on Cartesian) | round-off | -1.4..-3.4 ms at 4-8 GPUs across nodes; **off is faster by up to 0.9 ms at 2 GPUs on one node** | default by rank layout (off within one node) |

## 2. Changes results, real gain, thin evidence

| switch | result effect | gain | what is missing |
|---|---|---|---|
| M1 `time2_lin_tol_fac 10` (hesdirk2 stage Krylov tolerance 10x looser) | tolerance level, looser than lin_tol | -8 % (Eddington), -10 % (vet_sc), MI300A box only | stiff-radwave error with factor 1 vs 10 (fast3 found lin_tol 1e-9 gives 20x the error); timing on the wedge and on H200 |
| dhj `ck_impl_xstep 8` (504d8e30) | 1-10 bar within noise; **1e-3..1 bar at 1.2-1.5x a single noise member** | -11 % (MI300A WASP 1x), -19 % vs xstep 2 (H200) | a second noise member (2-3 member bar) |
| dhj `ck_store_split` (55762ada) | round-off; noise **1.21x** on H200 (1.11-1.25), 0.91x on MI300A | only **-4.3 %** at production nx1 76 on MI300A (whole merge) | re-gate the noise with 2-3 members |
| dhj `ck_beam_par` (WASP inputs; default under T4 since defaults-0927) | GPU not bitwise; noise 1.03x (nx1 76), 1.15x (nx1 256) | -6.2 % on 2 H200 at nx1 256, **measured before ck-next reworked ck_beam_tau** | MI300A nx1 76 timing with the current binary; keep only if >= ~3 % |
| dhj ck-next 6a9ffb0d, ck-lin2-fma 8be0ad67 | ck-next: GPU round-off (noise 0.78x); lin2-fma: bitwise on CPU and CUDA | -9..-11 % each on H200 | HIP gate + MI300A nx1 76 timing (see `TASK-2026-09-27-cklin2-hip-gate.md`) |

## 3. Never measured individually

`ck_impl_colskip` (default true since 09-22; can change results at tolerance level), `gas_newton` (default since
09-27; no timing A/B), `mg_levels 3` vs 2, and `vimp_fold`. One A/B each if a clean inventory is wanted.

## 4. Checked and fine (no action)

Bitwise, or a large measured gain (>= ~5 %) with round-off or tolerance-level differences:
- ck-jlin, nosync, the T4 set (fuse, jac_lin, cvsec), sweep_cache/tm, Kokkos 4.6.02;
- M1 fast path (halo_direct, od_cache, halo_mpi);
- op_stencil, krylov_fuse, bcg_sync, pcr;
- predictor step, lin_ew_max, precond mg (Cartesian), rbgs_fwd (sp);
- fast3/fast4/fast5 kernels, closure_lag step, vet_col predict.

Also fine, though the gain is small: `implicit_fast_kernels` (+1.5 %) and `predictor_order 2` (+3 %). Keep them.

## Also from Caltech 09-27 (for context)

- The Caltech timing harness had `problem/ck_impl_verbose=true`. Its CkCadStore diagnostics add **1.28 ms/cycle**
  through two contended global atomics per cell. If any viper deck or benchmark sets it, the same cost applies.
  Production default is false.
- Per-cycle ck profile on 1 H200 (nx1 256): lin1p + lin_sum are 45 % of ck; the storing kernels are 32 %.
- Hydro findings the user declined as not worth the effort (recorded only):
  - `splitpc2p` / `splitsc2p` ConToPrim tasks are redundant for dhj (2 of 5 c2p + raisev per cycle, ~3.2 ms/cycle);
  - PackAndSendCC runs with vector length 1 (~1.8 ms/cycle).
