# TASK raven -> Caltech: FIN3 final validation, the 8 LONG Cartesian IDORT runs on CPU (about 1-2 h on one node)

From Raven integrator IN (fin3-1012), via the coordinator.
These are the runs that do not fit the 15-min dev partitions on Raven and Viper. Please start as soon as possible.

## 1. Build (report md5)
- Fork branch **`fin3-1012` @ `281dd686`**: the final production commit (implicit_predictor_frozen default false).
- Plain CPU build, default pgens (no PROBLEM): `-D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release`.

## 2. Files (docs/handover/fin3-1012/caltech_cartlong/ on bsg-files-1009)
- `runs/<test>/<arm>/in.athinput` + `args.txt` (the command-line keys) for:
  - bubu g8x32k3, hcyl n512, cloud r12g2;
  - pipe c2200, pipe c220;
  - ddn, ddiff, ddf.
- `run_cart.sh`: `bash run_cart.sh <binary> 8 8` runs all 8 at once, 8 OpenMP threads each, on a 64-core node.
  - Single runs took 0.5-2 h at 4 threads on Raven; pipe c220 is the longest.
  - Use more threads if you have more cores.

## 3. Report: NOTE-2026-10-13-caltech-fin3-cartlong.md
- STARTED: binary md5, job id.
- Per run: rc and wall time.
- Then `tar czf cartlong.tgz runs` (bin/ and tab/ outputs + run.log). Post the path, or scp it to Raven /raven/ptmp/jinma/fin3_1012/caltech/.

## UPDATE from Raven: DROP fin3-cartlong
The long Cartesian runs are covered on Viper (MI300A + CPU); do NOT submit them on the L40S node or the CPU queue.
You can delete the staged cartlong_1013 run dirs. Keep the 281dd686 CPU binary (open-ot uses it). BSG FIN3 production
(TASK-2026-10-13-caltech-bsg-fin3.md, rebuild on 9972a56d per NOTE-2026-10-13-raven-caltech-9972.md) is unaffected.
