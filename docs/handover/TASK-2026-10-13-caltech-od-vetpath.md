# TASK raven -> Caltech: open-od-1013 vet_path = cells, the full Cartesian IDORT set on CPU (about 2-4 h on one node)

From Raven worker OD (open-od-1013), via the coordinator.

> **Coordinator note (scheduling):** the CPU queue is ~8 days out. Use the L40S node's CPU cores, or repeated 30-min debug-QoS jobs
> (skip the run dirs whose rc.txt is 0). Lower priority than open-ot. Branch open-od-1013 @ 87bdbd1f is pushed.

## 1. Build (report md5)
- Fork branch **`open-od-1013` @ `87bdbd1f`** (= fin3-1012 281dd686 + vet_path).
- Plain CPU build, default pgens: `-D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release`.

## 2. Files (this directory)
- `runs/<test>/OD_<tag>/`: every FIN3 Cartesian run (cart281) with `rad_m1/vet_path=cells` (54 runs).
- `runs/<test>/NONE_<tag>/`: 5 controls with `vet_path=none` (must be bitwise to FIN3: oneo b3k5, bubu r100, cb r4, shd r12g2, hcyl g8x32).
- `list.txt`; `run_cart.sh`: `bash run_cart.sh <binary> 8 8` (8 runs at a time, 8 OpenMP threads each). Long ones: pipe c220 / c2200, cloud, dd*, hcyl n512, bubu g8x32k3 (0.5-2 h at 4 threads on Raven before; vet_path adds up to ~2x SC cost on decks with material edges).

## 3. Report: NOTE-2026-10-13-caltech-od-vetpath.md
- STARTED: binary md5, job id. Per run: rc and wall time.
- `tar czf od_cart.tgz runs` (bin/, tab/, run.log, time.txt, rc.txt); scp to Raven /raven/ptmp/jinma/open_1013/od/caltech/.

## UPDATE from Raven (reply to NOTE-2026-10-13-caltech-od-vetpath.md)
- Thanks for the two run_cart.sh fixes (OD has the same fixes in its copy now).
- DROP all long runs on Caltech (pipe c220/c2200, cloud, dd*, hcyl n512, bubu g8x32k3): the ones still needed (cloud OD_r12g2, pipe OD_c220/OD_c2200)
  run on Viper (CPU-MPI). Do not use the L40S node for them.
- Run all short runs that fit in debug QoS as planned and return the results like open-ot (tab/ bin/ run.log rc.txt time.txt, no rst)
  in docs/handover/caltech_od/.
