# runs_5m_precond: a coarse correction for the implicit M1 BiCGStab preconditioner

- **Date:** 2026-09-25, viper. **Branch:** `m1-precond` (from rt-integration d2572b4f),
  worktree `/viper/ptmp2/jinma/wt_precond`.
- **Run tree:** `/viper/ptmp2/jinma/pc_0924`: `bin/` = WIP 3deb6b00, `bin2/` = final
  (fused prolongation), `base/` = d2572b4f; GPU runs in `runs/t1` (wedge),
  `runs/t2` (box levels), `runs/t3` (final timing), `runs/ts4` (T-S4 gate), `runs/bw`
  (GPU default bitwise); dumps + offline study in `runs/d_box1`, `runs/d_wed1`;
  CPU gates in `gates/`, `gates2/` (final binary), `rstmg*/`.
- **New switch:** `<rad_m1>/implicit_precond = mg`, with `implicit_mg_levels` (default 2,
  levels including the fine one) and `implicit_mg_halo` (default true).  Code in
  `src/rad_m1/rad_m1_precond.cpp`; the other values of `implicit_precond` are untouched
  (bitwise).  Debug: `implicit_dump_op = N` writes the frozen system of cycle N
  (`m1op.c<N>.r<rank>.bin`) for `study.py` (`m1op.py` reads it).
- Cases: box = the 3-D seeded He slab, vet_sc full tensor, 84x104x104, 4 blocks
  (`inp/box_vsc_m.athinput`; 2 GPUs: nx2 = 208, 8 blocks); wedge = the sp He wedge
  vet_col, 96x128x128, 16 blocks (`inp/hewedge_d.athinput`; 2 GPUs same grid).
  40 cycles, lin_tol 1e-10, iterations from `implicit_picard_log`, steps 10-40
  (`tsum.py`).

## 1. Diagnosis (offline, `study.py <dump> diag`, the system of cycle 20)

- The row excess (diagonal minus the sum of |off-diagonals|, over the diagonal) is tiny:
  box 5e-3 at depth falling to 4e-4 and slightly negative (-1.4e-3) at i = 70-77; wedge
  1.6e-2 at the base to 1.8e-3 at depth.  The system is nearly singular for modes that
  are smooth sideways: c dt rho kappa is large, the x1 coupling is 0.8-0.9 of the diagonal
  (wedge base: sideways 0.8, x1 0.17, the angular cells are small there).
- The x1 line solve is exact; red-black transverse Gauss-Seidel (block-local) damps high
  sideways wavenumbers only.  The error left after n rbgs_fwd iterations is:
  - box: 99.7 % in the sideways mean (|k| = 0) after 5 iterations, 98 % after 10, then
    |k| = 1-2 after 20 (runs/d_box1/diag.txt);
  - wedge: |k| = 1-2 (52-85 %), residual largest at the base (i ~ 18) where the sideways
    coupling dominates (runs/d_wed1/diag.txt).
  - No block-edge concentration (rms error at block-edge columns 1.1-1.3x the interior):
    the slow modes are global smooth sideways modes, not block-boundary artefacts.
- Offline BiCGStab iterations to 1e-10 on the dumped systems (`study1.txt`, `study2.txt`):

  | preconditioner | box | wedge |
  |---|---|---|
  | line | 29 | 20 |
  | rbgs_fwd (default) | 22 | 16 |
  | rbgs (symmetric) | 21 | 13 |
  | 2 forward sweeps | 18 | 13 |
  | rbgs_fwd, global colouring (no block cut) | 18 | 11 |
  | rbgs_fwd + block-local semicoarsened levels, 2 / 3 / all | 10 / 9 / 9 | 9 / 6 / 7 |
  | same, coarse levels across blocks (`mgb`) | 8 | 7 |
  | full global V-cycle, x1 never coarsened | 3 | 2 |

  So the lever is a coarse correction for the smooth sideways error; more or better
  sweeps give little.  Line solves in x2/x3 (ADI) were not built: the error is not
  high-frequency in x2/x3, and one more line direction is again a smoother.

## 2. The implemented preconditioner (`implicit_precond = mg`)

z0 = rbgs_fwd(r); r1 = R (r - A7 z0); z = z0 + P V(r1), with P piecewise constant over
2 x 2 (x2,x3) aggregates at the same x1 index inside a MeshBlock (x1 never coarsened),
R = P^T, Galerkin coarse rows R A7 P (7-point again, block-local), V = the same thing
recursively (one forward red-black line sweep per level, sawtooth, no post-smoothing).
The only communication is the halo of z0 for the fine residual (`implicit_mg_halo`).
The prolongation of all levels is one kernel (z(k,j,i) += sum_l z_l(k>>l, j>>l, i)).
Per application: 2 fine PCR launches (as rbgs_fwd), 1 halo, 1 residual/restriction,
2 PCR + 1 restriction per coarse level, 1 prolongation; coarse rows rebuilt once per
Picard pass.

## 3. GPU results (apudev, HSA_XNACK=1, HSA_NO_SCRATCH_RECLAIM=1, same binary
`bin2`, 2 interleaved repeats, runs/t3, job 11968168)

| case | preconditioner | it/solve (max) | it/step | ms/cycle (2 reps) | change |
|---|---|---|---|---|---|
| box 1 GPU | rbgs_fwd | 21.5 (50) | 29.8 | 56.86, 56.85 | |
| | mg, 2 levels | 12.9 (31) | 18.4 | 54.06, 54.08 | -4.9 % |
| | **mg, 3 levels** | **9.8 (22)** | 13.9 | **52.80, 52.75** | **-7.1 %** |
| | mg 3, halo off | 12.0 (26) | 17.1 | 55.86, 55.71 | -1.9 % |
| | mg 3 + implicit_precond_float | 9.8 (21) | 13.9 | 52.43, 52.23 | -7.9 % |
| | mg 4 / 5 levels (runs/t2, `bin`) | 8.6 / 8.3 | | 53.88 / 54.35 | -5.2 / -4.4 % |
| box 2 GPUs | rbgs_fwd | 23.0 (56) | 30.4 | 64.55, 64.54 | |
| | **mg 3** | **10.1 (24)** | 14.4 | **60.33, 60.53** | **-6.4 %** |
| | mg 3, halo off | 12.1 (29) | 17.2 | 62.11, 61.62 | -4.1 % |
| wedge 1 GPU | rbgs_fwd | 13.0 (21) | 13.0 | 43.39, 43.43 | |
| | mg 2 | 8.5 (12) | 8.5 | 42.41, 42.36 | -2.4 % |
| | **mg 3** | **7.3 (10)** | 7.3 | **41.97, 41.84** | **-3.5 %** |
| | mg 3, halo off | 8.8 (12) | 8.8 | 44.36, 44.27 | +2.1 % |
| wedge 2 GPUs | rbgs_fwd | 13.1 (21) | 13.1 | 24.82, 24.92 | |
| | mg 3 | 7.3 (10) | 7.3 | 25.37, 25.34 | **+2.0 %** |
| | mg 3, halo off | 8.8 (12) | 8.8 | 25.66, 25.66 | +3.2 % |

0 NON-CONVERGED in every run.

**Cost per iteration** (slope between lin_tol 1e-10 and 1e-8, same preconditioner):
box 1 GPU rbgs_fwd (56.86-47.39)/(29.8-15.9) = 0.68 ms, mg 3 (52.80-45.69)/(13.9-7.8) =
1.17 ms; wedge 1 GPU rbgs_fwd (43.39-33.36)/(13.0-5.1) = 1.27 ms, mg 3
(41.97-33.17)/(7.3-2.9) = 2.0 ms.  So mg costs about +0.5 ms (box) / +0.7 ms (wedge)
per BiCGStab iteration (two preconditioner applications), i.e. above the ~0.2 ms the
brief set as the break-even for a 2x iteration cut.  The iterations are cut 2.1x (box)
and 1.8x (wedge), and the net is -7 % (box) and -3.5 % / +2 % (wedge 1 / 2 GPUs).
The extra cost is kernel launches, not the halo: halo off saves ~0.04 ms/iteration
but costs 2-3 iterations per solve.  At lin_tol 1e-8 the gain remains (box 47.39 ->
45.69-45.81, -3.4 %; wedge 33.3 -> 33.2, flat).

The iteration count is now 7-10 per solve; the ideal (global V-cycle, 2-3 offline)
needs cross-block coarse levels (communication per level), not attempted.

## 4. Gates

- **Default bitwise:** CPU `gates.py` arms box_a1, box_m2, slab_a1, slab_m2, loose and
  tight, base d2572b4f vs new: rst bitwise, hst identical (gates/eval.txt).
  tst/test_suite/rad_m1 (ATHENAK_M1_DATA=/viper/ptmp2/jinma/faces_0924/m1data): 3 passed
  (slab, opcheck incl. implicit_op_check = 2, restart; tst_rad_m1.log).  GPU, 3
  Cartesian cases base vs new (slab 1 GPU, box and box vet_sc 2 GPUs): 9/9 files
  bitwise each (runs/bw, job 11969009).
- **Same fixed point** (`gates_mg.py`, gates.py criteria, final binary: gates2/):
  mg vs rbgs_fwd hst_cons 8e-11 (box) / 2.6e-9 (slab, F1top), rst_cons 2e-10 / 5e-10,
  hst_dyn <= 2e-8; mg 1 vs 2 ranks PASS both.  The F1top misses of the strict gate are
  what any other valid preconditioner shows: line vs rbgs_fwd 1.8e-10 (box) / 1.8e-9
  (slab), rbgs vs rbgs_fwd 1.8e-9 (slab), same binary (gates/, reference pairs).  So mg
  reaches the fixed point to solver tolerance as well as the existing preconditioners.
- **He slab / box 0 non-converged:** all CPU gate arms and all GPU runs.
- **Restart** (`rst_mg.py`: the CI restart test + mg levels 4, 2 ranks, Eddington and
  vet_sc): 16 cycles vs 8 + restart + 8 BITWISE (rstmg2.log).
- **Wedge T-S4** (`sph_atm_vc`, 800 cycles, runs/ts4, GPU): mg vs default dE/E 1.0e-9,
  dF/|F| 1.7e-9 (1 GPU, 4x4); 1.2e-9 / 3.2e-9 (2 GPUs, 8x8); transverse flux
  |F2,F3|/|F1| 7.6e-11 / 2.7e-11 (default 5.5e-11): unchanged to solver tolerance.

## 5. Verdict

- Use `implicit_precond = mg`, `implicit_mg_levels = 3` for Cartesian box runs:
  -7 % (1 GPU) / -6 % (2 GPUs) per cycle, 2.1x fewer iterations, max per solve 22 vs 50.
  With `implicit_precond_float = true` -8 %.
- Keep rbgs_fwd on the sp wedge for now: -3.5 % on 1 GPU but +2 % on 2 GPUs; the wedge
  cycle is dominated by VetColBuild and the bvals halo, and its solve is only ~10 ms.
- Remaining lever: the mg application costs ~0.25 ms (launch-bound coarse PCR kernels);
  fusing the coarse restriction into the PCR kernel or doing all coarse levels in one
  team kernel per block would bring the per-iteration cost toward the 0.2 ms target.

## Files

- `study.py`, `m1op.py`: offline study on dumped systems; `tsum.py`: ms/cycle and
  iterations per run and arm; `gates_mg.py`: CPU gates; `rst_mg.py`: restart gate.
