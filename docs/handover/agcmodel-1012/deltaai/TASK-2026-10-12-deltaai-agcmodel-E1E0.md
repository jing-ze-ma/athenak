# TASK raven -> DeltaAI: AG Car A 3D long-ramp runs E1 (deep inner wall) then E0 (old-wall control), <= 1 h link chains, 1 node x 4 GH200

From: Raven (worker AC, agcmodel-1012).

## Race
- Same runs as TASK-2026-10-12-caltech-agcmodel-E1E0.md (Caltech H200) and the Raven gpu queue (31106649 E1, 31105722 E0).
- **First result wins.** Post STARTED with the job id as soon as link 1 runs. If you cannot start within ~2 h, SKIP and say so.
- The physics, files, decks, analysis and NOTE content are the same as in the Caltech TASK. Read its sections "Why", 2 and 4.
  The files are in `docs/handover/agcmodel-1012/deltaai/` (identical copies).

## Environment
- `export MPICH_GPU_SUPPORT_ENABLED=1`, plain srun, 1 rank per GPU, 4 ranks (4 MeshBlocks; E0 480x64x64, E1 536x64x64).

## 1. Build (report md5)
- Fork branch **`prodfinal-1011` @ `ed32f0fa`**, built with `build_inc_deltaai.sh hes_gpu` (PROBLEM=he_star_m1).
- Use the newer hash instead if the coordinator's push message names one.

## 2. Runs: chained links, each <= 1 h (`-t 00:55:00`)
- **E1 first, then E0.** ghx4-interactive allows one job at a time, so submit each link when the previous one ends.
- Link 1: `srun -n 4 <bin> -i agcA_E1.athinput -d <run> -t 00:55:00`.
- Link n: `srun -n 4 <bin> -r <run>/rst/agcar3d.<newest>.rst -d <run> -t 00:55:00`. The code writes a final rst at the -t limit.
- **STOP rules:** rc != 0, FATAL, NaN or `Picard solve DIVERGED` (record cycle, t and the last 40 lines). DONE at tlim 2.4192e7 s.
- **Estimated length:**
  - E1: 2.6-3e4 cycles, about 4-10 links.
  - E0: 1.9e4 cycles, about 3-8 links.
  - The s/cycle rises in the ramp.
- **Cap:** 16 links per run. When the cap is hit, stop and report how far it got.
- Keep the rst at t = 4e6, 1.2e7 s and the last one.

## 3. Report: NOTE-2026-10-12-deltaai-agcmodel-E1E0.md (bsg-files-1009, docs/handover/agcmodel-1012/)
- Per link:
  - t, cycle, dt, s/cycle, wall;
  - the Picard line;
  - FATAL / NaN / NOT ADMISSIBLE counts.
- Interim hstA.py rows at t = 4e6, 8e6, 1.2e7, 1.6e7 s.
- Final:
  - the hstA.py, fluxA.py (`RAMP_S=4.032e6 RAMP_T=2.016e7`; E1 `ICFILE=<F>/ic_agcar_A_ge_deep70.txt`, E0 `ICFILE=<F>/ic_agcar_A_ge.txt`)
    and plume.py outputs, verbatim, plus GPU-h;
  - push the hst files, logs, prof.npz and plume.npz.
