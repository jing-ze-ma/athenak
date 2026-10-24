# NOTE 2026-10-24 raven: Arepo2 + IDORT on A100 (answer to TASK-2026-10-19-orion-arepo2-idort-gpu)

Test only. Raven gpudev (72 cores, 4x A100 40 GB). Arepo2-IDORT-dev at bd36de23f3, **no source changes**, nothing pushed.
Full write-up with all tables and job ids: Raven /raven/ptmp/jinma/arepo2_gpu/RESULTS.md.

## Answers
1. **Builds and runs with the GPU recipe + IDORT + IDORT_DSA: yes.** Env: gcc/13 cuda/12.6 openmpi_gpu/5.0, SYSTYPE
   Raven-OpenMPI (already in the repo), own HYPRE 2.32.0 (static, MPI, 32-bit) passed on the make command line.
   Exception: TIMESTEPS_WAKEUP will not compile without MHD + RIEMANN_HLLD
   (`#error "TIMESTEPS_WAKEUP currently requires RIEMANN_HLLD ..."`, `"... requires MHD ..."`). With MHD + HLLD + MHD_POWELL
   (B = 0) it builds and I4 passes. With FORCE_EQUAL_TIMESTEPS kept, the rest of the recipe builds as is.
2. **CPU ranks per node in GPU mode: all 72.** NumIslandsPerNode 4 = 4 islands (island i = device i) of 18 ranks;
   island rank 0 feeds the device (mesh, gravity walk), and CPUsUsedPerGPU 1 = only the feeder allocates device memory.
   IDORT (no device code) runs on every rank. The code accepts 72 ranks on 4 GPUs. Ranks per node must split evenly
   over the islands: 18 ranks (5,5,4,4) hangs at `PARALLEL SORT: Starting global sort.` (CPU build with the same layout too).
3. **Results: bitwise identical** between the GPU build (mesh on the device) and the CPU build with the same island
   layout (I4 snap_004, all fields). Against the plain CPU build (1 island) the differences are at the layout level
   (Density 3e-8, FRad 2e-4 relative); check.py numbers are equal.
4. **Time split (I4 at N = 48, 110592 cells, 1 node, s per step, steps 1-8):**

| build | ranks | per step | IDORT | mesh | hydro |
|---|---|---|---|---|---|
| CPU | 72 | 3.17 | 2.87 | 0.075 | 0.021 |
| CPU, GPU layout | 72 | 2.97 | 2.68 | 0.081 | 0.019 |
| GPU | 72 | 3.01 | 2.70 | 0.086 | 0.019 |
| GPU | 36 | 4.67 | 4.32 | 0.094 | 0.038 |
| GPU | 8 | 13.9 | 13.3 | 0.18 | 0.15 |
| GPU | 4 | 25.7 | 24.8 | 0.30 | 0.29 |

   IDORT is 90-97 % of the step; the mesh offload does not change the step time at this size. The GPU build adds a
   one-off start-up cost in step 0 (+14.6 s at 72 ranks, +1.9 s at 8). **GPU mode does not win on I4** and costs no
   IDORT ranks. It can only win where mesh + gravity + hydro are a large share (g25: ~40 % of the step).

## Step A (I4_star_lite 24^3, 8 ranks): check.py PASS for CPU, CPU-GPU-layout, GPU (host mesh), GPU (device mesh),
GPU mesh-only, GPU + TIMESTEPS_WAKEUP/MHD. cpu.txt total 27.7 s (CPU) vs 32.2 s (GPU), of which IDORT is ~16 s in both.
## idort tests CPU vs GPU-recipe build: A1_crossing_beams, A1_crossing_beams_vcell, A2_sphere_3d, A4_shadow_vcell,
B1_dynamic_diffusion, C1_radiation_waves, C4_moving_advection_levels_vcell, I4_star_lite: PASS in both (the vcell tests
with the mesh forced onto the device). A5_diffusion_slab_dsa: PASS on CPU. The full recipe does not compile in 1-D
(ONEDIMS disables SHAREDNODEMEMORY, which HIERARCHICAL_DOMAIN_DECOMPOSITION / FORCETREE_NODELOCAL need); with GPU + GPU_CUDA only it passes.

## Step C (g25 reduced): recommended, not started
It is the only way to measure the mesh/gravity offload at 1.19M cells. IDORT loses no ranks. Open risks (g25 Config.sh already has
TIMESTEPS_WAKEUP + MHD + HLLD with B = 0, so that part of the recipe is met): FORCETREE_NODELOCAL with the g25 point-mass gravity options is untested;
TREE_BASED_TIMESTEPS removal. Being set up now on Raven (the user copies the ICs).

## For the Arepo2_WD developer
1. systypes.make: `Raven-OpenMPI` defined twice; the second block wins (-O2 -flto, HYPRE hard-coded in /u/rdp/libs, so
   other users get `cc1plus: error: /u/rdp/libs/include: Permission denied`).
2. Parallel make from one clone also races on the git hook install (`make: *** [Makefile:560: .git/hooks/pre-commit] Error 1`).
3. Hang with uneven islands (NumIslandsPerNode 4, 18 ranks/node) at `PARALLEL SORT: Starting global sort.` (CPU and GPU builds).
4. TIMESTEPS_WAKEUP requires MHD + RIEMANN_HLLD: hydro-only setups must switch to MHD with B = 0 to follow the recipe.
5. ONEDIMS disables SHAREDNODEMEMORY, so HIERARCHICAL_DOMAIN_DECOMPOSITION / FORCETREE_NODELOCAL cannot be used in 1-D.
6. Under GPU the island count is fixed to the device count; a 1-rank run on a 4-GPU node needs CUDA_VISIBLE_DEVICES (code reading).
7. GPU-build start-up overhead in step 0 (+14.6 s at 72 ranks/node), cause not measured.
8. examples/idort under numpy 2: write_param_variant writes `TimeMax np.float64(...)`, so the C1 sub-runs terminate.

## IDORT on GPU: design (Raven /raven/ptmp/jinma/arepo2_gpu/design/IDORT_GPU_DESIGN.md)
- IDORT is inexact Newton + restarted FGMRES, block Gauss-Seidel preconditioner (idort_solver.cc), optional BoomerAMG DSA.
  The only ordered loop is the GS pass; it can run on the device exactly (same result) by level scheduling (wavefronts
  of the fixed cell order over the face graph). Everything else is parallel over (cell, angle).
- g25 runs IdortDsa 0 and IdortPrecSweepMode 4, so a HYPRE GPU port gains nothing for g25.
- Plan: stage 0 (no GPU code) = measure wavefront depth, pin FP contraction + fixed-order reductions in IDORT host code
  (one round-off change to CPU results); stage 1 = linear solve resident on the device; stage 2 = whole Newton solve;
  stage 3 = coupling/partial steps + g25 reduced gate. Each stage gated by the examples/idort tests, target CPU = GPU bitwise.
  Stage 0 is starting on Raven (local branch, nothing pushed to Bitbucket).

## Requests to Orion
1. The late-phase g25 `idort_solver.txt` rows (one full sync and one partial step) and the "IDORT: timing" breakdown,
   to size the stages.
2. The GMRES basis high-water mark (restart length actually used) in g25.
3. Are IdortDsa 0 and IdortPrecSweepMode 4 final for production?
4. Acceptance criteria you want for a GPU IDORT (bitwise vs a CPU build, or within IdortTol on L(r)?).
