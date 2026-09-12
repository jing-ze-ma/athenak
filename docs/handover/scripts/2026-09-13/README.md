# Analysis scripts, 2026-09-12 viper session

All were run against `/viper/u2/jinma/ATHENAK/bench/<run>/` on viper; each takes the run
directory at the top of the file or on the command line.  Nothing here writes into `run/`.

| script | from | what it does |
| --- | --- | --- |
| `shell.py` | `bench/RG_fofc/analysis/` | red giant: radial shell diagnostics from `bin/` (corona loading, M(>join), L_rad,out/L, vertex-column percentile ranks) |
| `floormap.py` | `bench/RG_fofc/analysis/` | red giant: where the floor/FOFC events are (active vs cube-vertex corner ghosts), per dump |
| `deathmap.py` | `bench/RG_fofc_long/deathmap/analysis/` | the death-column tool: seam/vertex distance of every anomalous cell, time series of the collapsing column (rho, e/rho, v_r, T proxy), T_proxy<100 K census.  This is what identified the seam-row floor -> EOS sentinel |
| `excess.py` | `bench/sp_excess/analysis/` | sp MHD-vs-hydro arms: E, KE_h, mass, dt, cumulative floors per arm at matched rot |
| `localise.py` | `bench/sp_excess/analysis/` | where the sp excess lives: (r, latitude, day/night) decomposition of the KE and energy difference between two runs' dumps |
| `excess2.py` | `bench/sp_excess2/analysis/` | bitwise/differential comparison of MHD-vs-hydro arms with matched solvers |
| `run_cases.py` | `bench/fofc_gpu/` | FOFC test harness: imports the 9 `tst/test_suite` FOFC modules and runs their own test functions against ANY binary (CPU, MPI, HIP, Debug); writes a results file |
| `compare.py` | `bench/fofc_gpu/` | diffs two `run_cases.py` result files (used for CPU vs HIP) |
| `bitwise_restart_bisect2_submit.sh` | `bench/RG_fofc_long/bitwise/bisect2/submit.sh` | **the bitwise-restart recipe**: run N cycles continuously writing a mid-run restart + ghost-inclusive dumps, restart from it, dump at the same cycle and cycle+1, compare active cells / each ghost class separately.  This is the harness that found the BC-reads-w0 defect; reuse it for any new restart regression |
| `0001-RT-ck-...patch` | `bench/prof/rt_new/` | unmerged `problem/rt_ck_sweep_split` (RT ck down/up sweep as parallel threads); 1.13x on the kernel, 2.8% wall -- shelved, kept for the record.  Base 2e3ae96e |
