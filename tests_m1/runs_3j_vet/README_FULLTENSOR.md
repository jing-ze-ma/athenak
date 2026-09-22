# FULL Eddington tensor from the SC VET closure (vet_tensor = full), 2026-09-22

The patch is not committed. It is in bench/m1_fullt_0922/fulltensor.patch and contains VET plus the full
tensor, relative to rt-integration 1715190b (`git apply --check` passes). The same content
relative to 92d36db3 is in fulltensor_vs_92d36db3.patch; it lacks only the later eig_min limit edit.
Gate builds: bench/m1_fullt_0922/{ref,build} (gcc/14, Release, serial, PROBLEM=box_convection;
ref = HEAD 92d36db3 + vet_sc.patch). A 1715190b + patch compile is in new1715/. The HIP build is athena_gpu.
Runs are in bench/m1_fullt_0922/runs and gpu_runs. The analysis scripts are runs/gate2.py and runs/plane.py.

## What changed

* `<rad_m1>/vet_tensor = uniaxial | full` (default uniaxial) and `vet_eig_min` (default 0). Both are
  read only when present, so inputs without them write byte-identical restarts.
* `VetFullTensor()` (rad_m1_vet.cpp) runs after the SC sweep. It stores the guarded D = K/J in
  `vet_cell(M1_VET_D11..+5)` and fills the ghosts: periodic in x2/x3, edge copy in x1.
* Solver: new helpers `M1DDiag(iw, vd, full, m, d, k,j,i)` and `M1POff(..., vd, full)`, and
  `M1OffDiv(..., vd, full)` (rad_m1_implicit.hpp). Every D_ab of the multi-D solve now goes
  through them: the x1 line row (D11 at i, i+-1), the face-flux kernel, the x2/x3 face fluxes, the
  7-point cell terms and Krylov coefficients (D22, D33), the off-diagonal divergence (D12, D13,
  D23) in the RHS and in the `operator` Krylov term, the transverse limiter theta, and the
  enthalpy coefficients a_d plus de0 in step (b). With full = false the arithmetic is identical
  term for term. ImplicitTridiagSolve, Thomas and PCR are untouched. No `this` capture: all views
  are copied to locals. The one hipcc -Wgpu-maybe-wrong-side warning, at `coupling` in the gas
  update, is pre-existing.
* Realizability guard: J <= 0 or trace <= 0 gives delta/3. D is normalised to trace 1. The
  eigenvalues come from Jacobi; if any is below vet_eig_min, the map
  lam -> emin + (max(lam,emin) - emin)(1-3emin)/(sum - 3emin) is applied. It keeps the
  eigenvectors, the trace 1 and every eigenvalue in [emin,1]. The raw SC tensor is PSD with
  trace 1 by construction because I >= 0 and the weights are positive, so at emin = 0 the guard
  should never fire. D11 >= 0 keeps the x1 row a column-dominant M-matrix with diagonal >= 1.
* Diagnostics: `<rad_m1>/dbg_opac_patch` with `dbg_opac_x1lo/x1hi/x2lo/x2hi` multiplies the
  opacities inside a box (default off, read only when given). `VetDump` also writes a
  `.plane` file with K/J, the (chi,n) projection and the guarded D.

## Gate 1: byte-identity (CPU, ref vs new binary)

| arm | result |
|---|---|
| default, tlim 3, bin/rst every 1 s | 16/16 files cmp-identical, log identical except timings |
| uniaxial vet_sc, tlim 3 | 16/16 files cmp-identical, log identical except timings |
| uniaxial vet_sc, 200 s (V2_ref vs V2_new) | every output file cmp-identical |

## Gate 2: seeded 2-D He slab, 200 s, lag = step, vpert 1e-3 (runs/gate2.py)

| arm | Picard mean/max | nonconv | posfb | min E | F1top/Fin (min/max) | dt | KE_1 / KE_2 |
|---|---|---|---|---|---|---|---|
| V2 uniaxial (ref) | 4.14 / 6 | 0 | 0 | 1.651e5 | 1.0000 (1.0000/1.0003) | 0.1612 | 1.29e26 / 1.81e21 |
| T3_tau (ref) | 4.13 / 5 | 0 | 0 | 1.651e5 | 1.0000 (1.0000/1.0003) | 0.1612 | 1.23e26 / 1.87e21 |
| **F_full** | 4.14 / 6 | 0 | 0 | 1.651e5 | 1.0000 (1.0000/1.0003) | 0.1611 | 1.38e26 / 1.90e21 |

Guard: 0 triggers in 3.34e6 cell-calls. max |D - K/J| = 4.4e-16 (trace round-off). Min D11 =
0.33314. The same run built on 1715190b (runs/N1715_full) has a user.hst identical to F_full.
In this slab the tensor has almost no structure the uniaxial form misses: max |K/J - D_uni| =
1.8e-4 (bottom cell), K12/J <= 7.5e-5. F_full differs from V2 by 7 % in KE_1 and 5 % in KE_2.

## Gate 3: a shadowing absorber (opacity x10 in x1 [-1.05e7, 5.5e6], x2 [1.273e8, 2.122e8], tau ~ 1 to 0.3; 100 s)

| arm | Picard | nonconv/posfb | min E | F1top/Fin (min/max) | dt | KE_1 / KE_2 |
|---|---|---|---|---|---|---|
| P_uni | 6.97 / 8 | 0 / 0 | 9.05e4 | 0.9998 (0.9981/1.0009) | 0.0845 | 1.19e30 / 1.26e30 |
| **P_full** | 6.97 / 7 | 0 / 0 | 8.53e4 | 1.0000 (0.9979/1.0009) | 0.0997 | 1.20e30 / 1.33e30 |
| P_emin (full, floor 0.30) | 6.97 / 7 | 0 / 0 | 9.33e4 | 1.0000 (0.9980/1.0008) | 0.0981 | 1.22e30 / 1.32e30 |

The patch drives violent physical motion (KE ~1e30) in both arms. Both stay stable.

Tensor structure (runs/plane.py on the vet `.plane` dumps):

| dump | max \|K12/J\| | max \|K22-K33\|/J | max \|K/J - D_uniaxial\| |
|---|---|---|---|
| step 0 | 3.0e-2 | 1.8e-2 | 3.0e-2 |
| call 700 | 6.4e-2 | 5.0e-2 | 3.9e-2 (about 12 % of 1/3) |

The solve reads D exactly (max |D - K/J| = 0 in the dumps). Min eigenvalue is 0.266.

Guard at floor 0: 0 triggers. At floor 0.30 it acted in 7.26e4 of 2.02e6 cell-calls (3.6 %)
and moved D by at most 0.043. The run stays stable.

## Gate 4: cost on the GPU (apudev job 11943125, 1 MI300A, 1 MeshBlock, same binary, interleaved x2)

| case | E eddington s/step | U uniaxial s/step | F full s/step | SC ms/call (U / F) | SC share of U / F step |
|---|---|---|---|---|---|
| 2-D 84x32, 400 steps | 0.0756, 0.0788 | 0.0886, 0.0858 | 0.0884, 0.0884 | 9.10 / 9.13 | 10.3-10.6 % / 10.3 % |
| 3-D 84x32x32, 100 steps | 0.225, 0.225 | 0.243, 0.244 | 0.239, 0.239 | 14.3 / 14.4 | 5.9 % / 6.0 % |

Timings are the log's "cpu time used"/nlim and the fenced VetReport SC seconds (gpu_runs/*/run.log,
log.out.11943125).

* Full vs uniaxial costs +0.03 ms per SC call; whole-run differences are within run-to-run noise.
* VET vs Eddington costs +13 % per step in 2-D and +6-8 % in 3-D. In 2-D the SC is launch-bound
  (84 per-layer launches of 32 threads, 9 ms), slower than on the CPU.
* The GPU binary predates two edits that affect no timed arithmetic: the guard-counter fix and
  the vet_eig_min limit.

CPU side note (login node, noisy): SC is 10.2 ms/call full vs 9.3 ms uniaxial (F_full, V2_new logs).

## Open

The multi-block/MPI sweep is still missing. The top Marshak factor q is still 0.5 against the SC
H/J = 0.57. Plane dumps are for k = ks only.
