# NOTE 2026-11-21ed orion -> caltech (USER): question for an H200 profile (H200PROF); arms to follow

The USER asked for an H200 profile. Two questions before we post the arms (about 2 arms of ~40 min total, 1 node; production has priority):
1. Is **Nsight Compute (ncu)** available on the H200 nodes (path or module)? Do the nodes allow **GPU performance counters** for a normal user, i.e. no ERR_NVGPUCTRPERM? A quick `ncu --query-metrics` / one-kernel test would tell. nsys works already (PF_* arms).
2. Which node could take ~40 min of profiling around/after 17:00 PST without touching fresh production or R_LF (02-16)? 02-10 was suggested.
The arm file (HG + RSG on B_v203, NR=4: nsys all ranks, 4 cycles, plus an ncu pass on rank 0 over the top-10 kernels) follows in a later NOTE.

-- orion (for the USER)
