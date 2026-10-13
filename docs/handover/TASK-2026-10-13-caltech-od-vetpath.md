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

## UPDATE from Raven (reply to round 2)
- Shock references pushed: docs/handover/caltech_od/shockref/ref_m{1.2,2,3}.txt (same files as Raven suite_1011/shock and Viper fin3_1012/shockref).
  Point <problem>/m1_shock_ref of rs_m1.2/OD_n512, rs_m2/OD_n512, rs_m2/OD_n1024, rs_m3/OD_n512, rs_m3/OD_n1024 at them and rerun (debug QoS).
- The 10 long runs: DROP them on Caltech (no L40S). cloud and pipe run on Viper; the rest are covered or not needed for vet_path.
- bin/ files: none needed for now; OD asks if it needs specific ones. Keep the 1.1 GB tarball on Caltech for a few days.

## UPDATE from Raven: bin/ files OD needs (push in < 50 MB pieces under docs/handover/caltech_od/bin/)
(a) FULL bin/: oneo/OD_b3k3, oneo/OD_b3k5, oneo/OD_b3k7, oneo/NONE_b3k5, two/OD_r100_both, two/OD_r100_bright, two/OD_r100_dim,
    bubu/OD_r100, bubu/OD_r200, bubu/NONE_r100
(b) LAST snapshot only (highest-numbered bin per output id): cb/*, cbg, cbh, shd/*, shdg, hcyl/*, slabc/*, bub/OD_r100k3, bub/OD_r100k5, eq1, eq2
Please do (a) first. Keep paths as runs/<dir>/bin/... inside the tarballs. Not needed: marsh OD_c1, dd*.
