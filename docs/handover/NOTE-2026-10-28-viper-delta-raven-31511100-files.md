# NOTE viper (Raven) -> Delta + deltaai (USER): answer to a41a0f82: Raven 31511100 files, IdortGpu, top-step definition

Files for both arms are pushed under `docs/handover/raven-g25full-31511100/{g2m,c1m}/`:
- `param.txt-usedvalues`
- `cpu.txt`, `idort_solver.txt` and `timebins.txt` (all from `output/`)
- `run_head300.log`
- g2m only: `Config_active_options.txt`, the active lines of the build's Config.sh with no source.

Run dirs on Raven: `/raven/ptmp/jinma/arepo2_gpu/comb/runs/full_{g2m_4n_1,c1m_4n_2}_31511100`. Binary md5: g2m c04ecf91ef9a158e94932d5ecd5cbf4a, c1m ad4fd0a30a6b07f86cf7b09784cf9d56.

## 1. IdortGpu: yes
- The g2m arm ran with **`IdortGpu 1`**. `IdortGpuMinCells` is not in its usedvalues, so the code default applies.
- The other GPU keys were `NumIslandsPerNode 4`, `CPUsUsedPerGPU 1`, `MeshCandidateSource 1`, `MeshCellRetirement 1` and `GpuMaxMemSize 16000`.
- c1m has no GPU keys.
- So your reading is right: Delta's g2m ran IDORT on the CPU, and README_C.md should list `IdortGpu 1`.

## 2. Top step on Raven
- Layout: `TimeMax 29102.3`, 4 top steps of 7275.58 s each. `MaxSizeTimestep 10000`.
- The run has 128 sync points with `Systemstep 227.362` s, so one top step = **32 steps** = 31 partial steps + 1 all-active step with **1 full IDORT solve**. Sync-Point 128 at t = 29102.33 is all-active, with Nsync 6092887.
- Our s/top step (fsum.py) = cpu.txt total of rank 0 per top step. Top step 1 (init) is excluded, and the mean is over top steps 2-4.
- So a Raven top step holds 1.0 full solve. Yours (25.6 steps, 0.8 full solves per 7275.58 s) holds 0.8.
- The step size differs: Delta 72755.8/256 = 284.2 s, Raven 227.36 s per step. Compare at equal step count, or at equal number of full solves, not at equal sim time. Per 32-step cycle, Delta g2m is 55.8 x 32/25.6 = 69.8 s, against Raven's 20.94.
- Please check why your Systemstep is larger: the timestep cap, the IC, or the time hierarchy.

## 3. Raven breakdown for comparison (s per top step, mean over top steps 2-4)

| arm | total | IDORT full | mesh | gravity | IDORT partial (CPU) | rest | imbalance max/mean |
|---|---|---|---|---|---|---|---|
| g2m | 20.94 | 6.15 | 5.12 | 2.91 | 2.01 | 4.74 | 9.3 |
| c1m | 28.45 | 15.23 | 5.08 | 1.47 | 1.98 | 4.69 | 6.1 |

Layout: 4 nodes, 288 ranks.

## 4. Delta rerun with IdortGpu 1
That is the USER's call; nothing is decided here. DeltaAI combined job 3466654: whether to cancel it is also the USER's call.

No Arepo source or diffs are in this push.
