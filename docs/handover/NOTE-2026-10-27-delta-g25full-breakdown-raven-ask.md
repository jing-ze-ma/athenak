# NOTE 2026-10-27 Delta -> viper (Raven) + deltaai (USER): g25_full breakdown. CORRECTION to b0301516: Delta's g2m arm ran IDORT on the CPU (IdortGpu not set, default 0), so its 55.8 s is not a GPU-IDORT number. Ask Raven for the 31511100 files

Follow-up to b0301516 (hsmlfix PASS + 3-arm timing, job 23248921).

## 1. Correction: the g2m arm ran no GPU radiation solves
- Every `IDORT: timing` line in g2m/run.log ends `gpu io 0.0000 (0 solves)`, including all 9 full (all-cell) solves.
- Cause: Delta's job script (built from data/g25 README_C.md / run_g25.sh) adds the GPU block
  `NumIslandsPerNode 4, CPUsUsedPerGPU 1, MeshCandidateSource 1, MeshCellRetirement 1, GpuMaxMemSize 32000, MaxMemSizeNode 30000`,
  but **not `IdortGpu`**, whose default is 0. That recipe does not mention IdortGpu at all.
- The deck meets every GPU-path condition in idort_island.cc: no DSA, IdortPrecSweepMode 4, Krylov restart 100 (< MAX_VEC - 1 = 127), no ONEDIMS_SPHERICAL.
  Full solves have 6.09 M cells / 16 islands = 380 k cells per island, above IdortGpuMinCells 1e5.
  So with `IdortGpu 1` the whole Newton solve should run on the device.
- **So b0301516's g2m 55.8 s/top step = GPU mesh + GPU gravity + CPU IDORT.** It does not measure the GPU recipe.
  c1m and c0 are unaffected. The hsmlfix PASS is unaffected (it is in the init density, before any IDORT solve).
- Raven 31511100 shows "IDORT full sync" GPU 6.15 vs CPU 15.23 s, so Raven's g2m clearly ran device solves. Please confirm Raven set `IdortGpu 1` (and IdortGpuMinCells); if so, README_C.md should list it.

## 2. Breakdown of job 23248921 (s per Delta top step; steps 1-256 / 10)
Delta "top step" = 7275.58 s sim = 25.6 substeps (TimeMax 72755.8 / 10). The full IDORT solves fall every 32 substeps (steps 32, 64, ..., 256): 8 of them in 10 top steps.
IDORT numbers are from idort_solver.txt; cpu.txt has no IDORT line, so IDORT sits inside cpu.txt "misc".

| section | g2m (IDORT on CPU) | c1m | c0 |
|---|---|---|---|
| **total** | **55.8** | 62.9 | 67.2 |
| IDORT step (t_step) | 34.5 | 40.2 | 39.4 |
| of which the 8 full solves | 30.7 | | |
| IDORT t_comm (halo) | 7.9 | 10.7 | 9.9 |
| IDORT t_orth (GMRES) | 8.0 | 9.9 | 10.0 |
| IDORT t_reduce (part of orth + rest) | 6.5 | 9.0 | 9.2 |
| IDORT t_sweep | 7.0 | 6.2 | 6.2 |
| IDORT t_transport | 2.8 | 2.7 | 2.6 |
| voronoi (mesh) | 7.8 | 11.6 | 11.1 |
| treegrav | 6.3 | 2.8 | 8.2 |
| of which nltreeimbal | 2.7 | 0.8 | - |
| hydro | 0.9 | 1.2 | 1.2 |
| rest (domain, kicks, timesteps, i/o, logs, misc outside IDORT) | ~6 | ~7 | ~7 |

- Per full solve: 30-55 s (wall max), 6-9 Newton outers, 67-115 matvecs. ~95 GMRES iterations and ~110 matvecs per top step.
- The CPU IDORT solve on Delta is communication-bound: halo + orth + reductions are ~16-21 s of the ~31-35 s solve; the cell work (block + lin + transport) is ~4 s.
- GPU gravity is 2.3x slower than c1m gravity, mostly imbalance wait. Raven sees the same direction (2.91 vs 1.47 s).
- GPU mesh: 1.5x faster than CPU (7.8 vs 11.6 s). Raven's mesh is 5.1 s in both arms.

## 3. Ask to Raven (viper): the files of job 31511100 (4 nodes, 288 ranks, g2m + c1m)
Please push (or copy next to data/g25) for **both arms**:
1. `param.txt-usedvalues` (to diff against Delta's: IdortGpu, IdortGpuMinCells, Krylov / sweep settings, island layout).
2. `output/cpu.txt` and `output/idort_solver.txt` (Delta compares IDORT t_comm / t_orth / t_reduce / t_sweep per full solve and per matvec).
3. `output/timebins.txt` or the first ~300 lines of `run.log`, and your **definition of "top step"** (sim-time length and how many full IDORT solves fall in one).
   Delta's top step holds 0.8 full solves; if Raven's holds 1, the 20.9 vs 55.8 comparison is off by that factor too.
4. The Config.sh of the g2m build (to check that Delta's d9c0 / hsmlfix builds match e93504ccef's options).

## 4. Next on Delta (pending the USER's go, nothing submitted)
- Rerun g2m only with `IdortGpu 1`: same builds_hf binary, same layout, 10 top steps. Estimate ~12 min, 16 GPUs, interactive x2 = ~6.5 GPU-h.
- The job script will add IdortGpu to the GPU block, and check that the full solves report `> 0 solves` in the gpu io field.

**SU:** this NOTE used no jobs. Delta's bsg-files-1009 total stays 42.08 GPU-h; 952 of 997 left.

All changes stay local. No Arepo source or diffs in git.
