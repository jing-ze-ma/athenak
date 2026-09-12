---
name: fofc-gpu-verified
description: FOFC (orion b3345f55) verified on the viper MI300A GPU 2026-09-12 -- HIP Release 48/48 test cases pass, drift vs CPU is round-off; the HIP Debug aborts were a PRE-EXISTING Cartesian x1v over-read in the hydro dt-collapse diagnostic, fixed b18f3205
metadata:
  type: project
---

Harness: bench/fofc_gpu (run_cases.py imports the 9 tst FOFC test modules and calls
their own test functions against any binary; compare.py diffs two result files;
submit.sh = the apu job; builds build_hip_{dbg,rel}, build_cpu{,_mpi,_dbg} at e3410692).
apudev MaxTime is 15 min, so the ~40 min job ran on apu (11620418).

Results: CPU 48/48 (reproduces orion's numbers exactly), HIP Release 48/48; 13 cases
differ from CPU only at round-off (div B ~1e-14, MPI-vs-serial ~1e-16, vceil/FOFC event
counts by 0.1%). HIP Debug: 17 Cartesian hydro cases aborted on
`x1v out of bounds [0,3] extents [1,1]` -- NOT FOFC: hydro_newdt.cpp's dt-collapse
diagnostic (runs on the first NewTimeStep of every run) read pcoord->x1v on Cartesian
grids where it is a (1,1) dummy. Reproduces with fofc=false and with linear_wave_hydro.
Value feeds only the printed "r=" field; Release physics never affected. Fixed b18f3205
(CellCenterX from the block size on Cartesian). MHD has no such diagnostic.

**How to apply:** treat orion's "build Debug once before trusting a GPU FOFC run" item as
DONE; a bounds-checked Debug run of ANY Cartesian hydro input is the regression check.
The style gate (`run_test_suite.py --style`) is unusable on this branch: thousands of
pre-existing cpplint errors in pgen/*, hydro.hpp, mhd.hpp, and flake8 needs `python`.
See [[validate-the-instrument]].
