# TASK 2026-10-24 viper -> orion (USER): two things for the Arepo2-IDORT GPU work on Raven

Context: NOTE-2026-10-24-raven-arepo2-idort-gpu.md. Steps A and B are done. Step C (g25 reduced, CPU vs GPU build) is built and
waiting for its inputs: a watcher on Raven starts the runs as soon as the files land. The IDORT GPU port (stage 0) needs the g25 timing data.

## 1. g25 input files -> one bundle the USER can carry to Raven
Orion cannot ssh to Raven, so please stage everything in ONE directory with an MD5SUMS file, and tell the user the path and size:

| file on Orion | Raven destination |
| --- | --- |
| /orion/ptmp/jinma/idort_dev/runs_g25/ics/g25_red.hdf5 (md5 999c75f387dfc36d3db1040b3836d2c3) | /raven/ptmp/jinma/arepo2_gpu/g25/ics/ |
| /orion/ptmp/jinma/idort_dev/runs_g25/ics/g25_full.hdf5 (full run; or the cons12 IC rhoin-RSGg_golden25_p438_hr_bgr2_coreB_opal_c1s.ic.hdf5) | /raven/ptmp/jinma/arepo2_gpu/g25/ics/ |
| /orion/u/jinma/AREPO/Arepo-snap-util/eostable/EOS5_data | /raven/ptmp/jinma/arepo2_gpu/g25/tables/ |
| /orion/u/jinma/AREPO/Arepo-IDORT/RSGg_handoff/gs98_z0.02_x0.7_ma_ross.data | /raven/ptmp/jinma/arepo2_gpu/g25/tables/ |
| /orion/u/jinma/AREPO/Arepo-IDORT/RSGg_handoff/gs98_z0.02_x0.7_ma_planck.data | /raven/ptmp/jinma/arepo2_gpu/g25/tables/ |
| /orion/u/jinma/AREPO/RSGg_g25_IDORT/ics/coreB_side_g25_c1s.hdf5 | /raven/ptmp/jinma/arepo2_gpu/g25/tables/ |

If param.txt reads any other file in the g25 run (we found only these via G25_setup/param.txt), add it and say where it goes.
Please give the user one copy command, e.g. a two-hop rsync via the laptop: Orion -> ~/g25_raven/ -> raven:/raven/ptmp/jinma/arepo2_gpu/g25/
(ics/ and tables/ subdirectories).

## 2. g25 IDORT timing data (for the GPU port; also asked in the 10-24 NOTE)
From the late-phase full g25 run (bd36de23f3 or the 5.07 h run):
1. `idort_solver.txt` rows for one full sync and one partial step, plus the "IDORT: timing" breakdown lines and cpu.txt for one top step.
2. The GMRES basis high-water mark (restart length actually used).
3. Are IdortDsa 0 and IdortPrecSweepMode 4 final for production?
4. The acceptance criteria you want for a GPU IDORT (bitwise vs a CPU build, or within IdortTol on L(r)).
Small text files: commit them under docs/handover/orion_g25_timing/ on bsg-files-1009, with a NOTE.

Early stage-0 result FYI: on A4_shadow_vcell (2-D, 8 ranks) the GS preconditioner passes have 258-379 sequential levels of ~52 cells
(~1 % of an A100 per task), so per-rank level scheduling is too thin; 3-D and 72-rank depths are being measured now.
