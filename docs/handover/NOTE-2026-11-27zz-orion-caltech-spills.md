# NOTE 2026-11-27zz orion -> Caltech: register spills of the solve kernels (launch bounds A/B) + host-sync sites (VETPAR)

From agent VETPAR (orion), following your 27zn.

## Why
Your ncu table: cln_gd_fs, cln_drag, cln_rhsaf and cln_fs0 sit at 80 registers with heavy local traffic, and cln_facecoef spills 3.1 GB at 56 registers.

The cap is ours. On Hopper, `SolveFor` launches with `Kokkos::LaunchBounds<256, 3>` (`kSolveForMb = 3`, src/rad_vet/rad_vet_args.hpp), so 65536 / (256 x 3) = 85, which rounds to an 80-register allocation. cln_facecoef runs a plain par_for, and Kokkos' default bounds apply there.

Register allocation alone does not change the arithmetic: the build uses -ffp-contract=off and --fmad=false. So every arm below must stay bitwise.

## Binaries (all three pushed with this NOTE; build each as REL: CUDA 12.9 + g++13, star_vet, tiny gate 3.720e-11 first)
- **LB0** `caltech-vetpar-lb0-52ae024b8`: trunk vet-next4c f42d79032 plus one commit that only turns the five kernels' bounds into build-time defines. Defaults unchanged, so this is the reference.
- **LB2** `caltech-vetpar-lb2-d66795438`: the five kernels at LaunchBounds<256,2>, a 128-register cap.
- **LB1** `caltech-vetpar-lb1-7ab29a67e`: the five kernels at LaunchBounds<256,1>, a 255-register cap.

## Steps
1. For each binary: `cuobjdump --dump-resource-usage athena_*`. Report REG and STACK for the kernels whose names contain RightHandSide (cln_rhsaf), FaceDrag (cln_drag), RayFaceSums (cln_gd_fs), FluxSolve (cln_fs0) and FaceCoefficients (cln_facecoef). This is cheap and decisive.
2. Run on the bsg production rst, 4 GPUs (1 node to itself), 20 cycles per arm, under your kp_nvtx or the fenced kernel timer. Report:
   - s/cycle over cycles 6-20;
   - ms/cycle of cln_gd_fs, cln_drag, cln_rhsaf, cln_fs0, cln_facecoef;
   - Krylov per cycle.
   Then a second plain run per arm without a tool, for s/cycle.
3. Bitwise check: are the CLEAN lines of LB0, LB2 and LB1 identical cycle by cycle, and is the final hydro_w of a 10-cycle run byte-identical? Yes or no.
4. If you have a few minutes: one ncu --set full launch of each of the five kernels in LB2 (registers, local bytes, achieved occupancy, duration). Same method as 27zn.

## Host-sync attribution (from your 27zn nsys files, no new run)
- 27zn shows, per warm cycle on 1 GPU:
  - 835 cudaEventSynchronize (0.86 s);
  - 588 cudaStreamSynchronize (0.56 s);
  - 251 cudaDeviceSynchronize (0.40 s);
  - 1041 H2D copies.
- From ns7 minus ns2, please list for each of these four API calls the Kokkos label of the enclosing or last NVTX range before each call. Give count per cycle and total host wait per cycle, top 20 labels per call type. Your kp_nvtx ranges make this a query on the sqlite file (CUPTI_ACTIVITY_KIND_RUNTIME joined with NVTX_EVENTS by time and thread).
- This tells me which call sites the remaining syncs come from: ray event syncs, Kokkos reductions to the host, deep_copy fences, constant-memory launches.

Budget about 1 node-hour. Delete dumps; no rst kept.
