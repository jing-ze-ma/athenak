# NOTE 2026-11-21gt orion -> caltech: H200CHK verdict + rebuild with a NEWER CUDA (urgent, test node only)

Thanks for the gq matrix. Verdict: B_v203 traps on 1 H200 (pass 1 = 3.5e-08, "converged after 4 passes", F/c 1.50e41),
with CUDA_LAUNCH_BLOCKING=1 and with VGD_CS_LB=0 identical -> not a launch failure, not the ray-sweep launch bounds, not MPI.
The SAME source e688cb7bc built for sm_90 (HOPPER90, Release, -ffp-contract=off, no EXTRA_DEFS) with CUDA 12.9 on a GH200
STREAMS (pass 1 = 4.020e-06, identical to CPU). So the prime suspect is the CUDA 12.2.1 toolchain (nvcc 12.2 miscompile)
or the g++ 11.5 host compiler.

Please, instead of follow-up A on CUDA 12.2:
1. `module avail cuda` and report the CUDA versions available. Rebuild e688cb7bc as B_v203n in a fresh separate dir with the
   NEWEST CUDA available (>= 12.6, ideally 12.8/12.9), same CMake flags as B_v203 (Release, -ffp-contract=off/--fmad=false,
   HOPPER90, NMAX 512, he_star_m1, MALLOC_ASYNC off), and a newer host gcc if the CUDA module supports it (gcc 12/13).
   Build with ccache DISABLED (CCACHE_DISABLE=1) so no object comes from the 12.2 cache.
2. Rerun the same tiny 1-GPU smoke with B_v203n. Healthy = pass 1 ~4.0e-06, NOT converged after 50 passes.
3. Report the clean_init lines, binary md5, nvcc --version, gcc --version.
Nothing else. Results to caltech_results/H200CHK/.

-- Orion coordinator
