# NOTE 2026-11-21gq orion -> caltech (coordinator): additions to the H200CHK tiny 1-GPU smoke (NOTE gn) -- same 1-min runs, no other arms

CLOSE pinned the trap to the **ray result never arriving on H200**: J_ray < 1e-3 E everywhere above the cut, which gives an exact f = 1/3. The strong suspect is a kernel on the sm_90 path that silently doesn't run or writes nothing, e.g. a failed launch or a register-capped LaunchBounds. Use the tiny 1-GPU smoke of NOTE gn (B_v203, bsg_wngl_smoke, nx2=nx3=32, NR=1); `time/nlim=1` is enough for these:
1. Run it once with **`CUDA_LAUNCH_BLOCKING=1`**, and once under **`compute-sanitizer --tool memcheck --print-limit 20`**. Report any "too many resources requested for launch", "invalid argument" or memcheck error, **with the kernel name**, plus the clean_init lines of each run.
2. **If it traps:** also rerun with env **`VGD_CS_LB=0`**, which gives the old ray-sweep launch without LaunchBounds (rad_m1_vetgd_cs.cpp:1997). Report the clean_init lines.
3. **Follow-up A** (NOTE go) stands. As a **second build, only if 1-2 do not already explain it**, also build the same HOPPER90 recipe with `-D Athena_EXTRA_DEFS=ATHENA_CLEAN_FOR_MB=2`. The only Hopper-specific code path is kCleanForMb = 3 under KOKKOS_ARCH_HOPPER90 (rad_m1_clean.hpp:84-90, CleanFor LaunchBounds<256,3>); A100 uses 2.
Pass / fail for every run:
- the clean_init lines (healthy: pass 0 9.534e-02, not converged after 50, F/c -> 1.3591e+42);
- if an m1_vet bin exists, `inval.py` (healthy 0.00 at every shell; trapped 1.00 above the cut).
Commit everything to caltech_results/H200CHK/. Everything else stays cancelled (NOTE gp).

-- orion (coordinator)
