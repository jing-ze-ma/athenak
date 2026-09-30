---
name: rad-m1-implicit-gpu-survey
description: Literature/code survey 2026-09-21 on implicit radiation transport on GPUs for stage 3 of the M1 module; recommended architecture and the key references
metadata: 
  node_type: memory
  type: project
  originSessionId: e67ee046-aa0b-46f4-881b-84099bfda99a
  modified: 2026-09-21T01:41:22.778Z
---

Survey by an Opus agent (not independently verified by me; items it marked unverified: HERACLES'
solver, VETTAM's KSP/PC, Enzo/ZEUS-MP, Castro GPU radiation).

- No published code does implicit two-moment/VET transport on GPUs at scale. Closest analogue:
  ARK-RT (Bloch et al. 2021, arXiv:2011.13926; Kokkos, same AP-HLL M1 as ours): Newton outer +
  BiCGStab on an ASSEMBLED matrix (Trilinos); AMG ~10 it, ILU/Schwarz ~20, damped Jacobi ~250;
  iteration count doubles at >= 4 MPI ranks; implicit/explicit gain 160x on CPU but only 11x on
  GPU; AMG is SLOWER on GPU than CPU, damped relaxation is the best GPU preconditioner.
- AREPO-IDORT (user's paper): Jacobi-like iteration, CPU only, GPU listed as future work.
  AthenaK's own GPU radiation (arXiv:2608.22020 per the agent) is explicit transport.
- Upstream origin/feature/multigrid: constant-coefficient Poisson only, red-black GS/SOR on
  device, AMR yes, but coarse tail below one MeshBlock is host-serial (round trip per V-cycle);
  variable-coefficient scaffolding (coeff arrays, CalculateMatrixPack) declared, unimplemented.
  No Krylov infrastructure anywhere in AthenaK.
- Published clean route: Schur-eliminate F so the solve is a scalar diffusion-like operator in E
  (Olivier et al. arXiv:2404.17473); DSA theory -> tau-independent iteration counts if the
  diffusion preconditioner matches the discrete transport operator. Nonlinear GMG for M1 with a
  realizability guard: Bloch et al. 2022 arXiv:2208.14703 (prototype only).
- A tau-gated explicit-M1/implicit-diffusion hybrid is NOT in the literature (novel; hazard at
  the handover surface).
- AMD GPUs: PETSc HIP matrices "in development" (Kokkos backend only), hypre BoomerAMG GPU
  partial, AmgX CUDA-only, Ginkgo has HIP. External libraries = validation oracle at best.

Recommended stage-3 architecture: matrix-free BiCGStab on the Schur-reduced scalar E operator
(face coefficients = our ap_hll c/(3 tau_face)), Picard outer loop with lagged Eddington tensor
and opacities, preconditioned by the existing radial block-tridiagonal column solve + ADI
transverse lines (line relaxation; already on GPU and partitioned across MeshBlocks/ranks).
Weak point: no coarse space -> iterations grow with MeshBlock count; eventual fix = geometric
multigrid (extend feature/multigrid with a line smoother) used as the preconditioner. Global
reductions are the GPU cost: Chebyshev smoothing is reduction-free. See [[rad-m1-design]].

ADDENDUM (same day): strongest precedent for the recommended route = Enzo FLD (arXiv:0901.1110):
inexact Newton + Schur complement to a scalar radiation equation + hypre-MG-preconditioned CG,
near-constant iterations to 4096 cores (uniform grid, CPU). Published iteration counts: RAMSES FLD
(arXiv:1102.1216) CG+diagonal ~11/step, and for pure diffusion STS ~10x cheaper than CG; Jiang
2021 Jacobi-like, no preconditioner, 20-1000 iterations; IDORT Nitermax 10, tol 1e-8, DO costs
~10x hydro; VETTAM Picard (chosen over Newton) + PETSc GMRES/ASM, ~10x hydro per step, and an
upstream wavespeed correction at AMR level boundaries "crucial" in the diffusion regime; CASTRO
GMRES preconditioned by hypre SMG/PFMG (operator nonsymmetric from v/c terms); MPI-AMRVAC FLD is
CARTESIAN ONLY because point relaxation fails for spherical operators -> a line/column smoother
is not optional for us. AthenaK GPU radiation (arXiv:2608.22020): explicit, radiation lowers
throughput 3.4-3.7x vs GRMHD alone = the in-code number an implicit scheme must beat. HERACLES'
solver never verified (A&A 403): do not cite it.

UPSTREAM STATUS checked 2026-09-21: multigrid is ON upstream main since 2026-07-31 (PR #776,
8a6a8efa, Velasco-Romero): src/multigrid/* + src/gravity/mg_gravity.*; only user = self-gravity
(Poisson); CalculateMatrixPack still declared, never defined; feature/multigrid is 43 header lines
ahead of main, last commit 2026-07-17. Also new on main: RKL2 super-time-stepping for diffusion
(#779, 2026-08-10, Fielding). rt-integration branched from main at a1b2afab (2025-12-19) and is 72
commits behind: using the multigrid means merging upstream main first (src/multigrid absent here).
