# NOTE 2026-11-16 Orion -> viper (USER): build the GPU tests from the integration branch

The user asks you to use the clean scheme's integration branch for the GPU tests:

  fork jing-ze-ma/athenak, branch **clean-cc-int**, head **4680ce500**

It contains everything in origin/clean-cc dae5c4133 (push 8 + CUDA fix) plus: ray-set rotation for the
clean FV rays, the vet_sc band fix (cb runs multi-block), deterministic axisymmetric FV sweep,
non-blocking halo in the operator apply, opacities at the Newton temperature (kappa(T)), ic_rad off under
clean, wb_arad reference work, the asymptotic-preserving lld ray sweep (clean_rays = lld), FV/lld rays on
Cartesian meshes (opt-in; SC stays the Cartesian default), cubed-sphere seam slopes, the FV bracket limiter,
and the temporary test keys clean_line_solver and clean_dissip_power.

How this changes TASK-2026-11-16-viper-cc_gpu1-a100 (stage 1):
- Build clean-cc-int 4680ce500 instead of dae5c4133 (same build_cuda_raven.sh; the CUDA fix is included).
- G5 timing (the main deliverable) needs no references: start it now on 4680ce500.
- G1-G4: the shipped references are for dae5c4133. Orion regenerates them on 4680ce500 and pushes an
  updated ref/ within about an hour; run G1-G4 after that (or run now and compare later).

The decisions package (second task: line solver, dissipation power, rotation, direction sets, fv1 vs lld)
follows on the same commit, with its own TASK file, shortly.
