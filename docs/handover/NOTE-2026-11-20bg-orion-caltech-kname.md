# NOTE 2026-11-20bg orion -> caltech (USER): one CPU-only post-processing job on the 3 PF sqlites (kname.py)

Please run `python3 docs/handover/caltech_runs_1121/batch4/kname.py <dir>/a/nsys_r0.sqlite 1 3 80 > kname.txt` (pure python3 + sqlite3) for
PF_rsg_v20g, PF_rsg_h200 and PF_hg_h200. Put each output in `caltech_results/<arm>/kname.txt`. It keys kernels by function name, because the
nvtx labels came out as "?". Low priority: CPU only.
For any future nsys arm, if a Kokkos nvtx connector (libkp_nvtx_connector.so) can be built from kokkos-tools, set KOKKOS_TOOLS_LIBS to it on rank 0.

-- orion (for the USER)
