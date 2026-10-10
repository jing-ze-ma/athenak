# TASK viper -> Caltech: H200 timing of the vet_gd halo speed-up (vgdspeed-1009), short job (user GO 10-09)

## What changed

Branch `vgdspeed-1009` (fork) changes only src/rad_m1/{rad_m1.hpp, rad_m1.cpp, rad_m1_vetgd.cpp, rad_m1_vetcol.cpp}.
It is bitwise identical to rt-integration b2f2e897 on AG Car A hr, BSG hr reduced and He giant fresh: CPU 4 ranks, and
4 A100 on Raven, including a run-to-run repeat.
- The per-shell compact halo uses cached gather/scatter lists per (pass, shell). They are rebuilt only when the
  direction set rotates (`vet_gd_halo_list_mb`, default 4096 MB per rank, capped at half the free GPU memory).
- The shell kernel runs only over the directions of the pass's branch (`vet_gd_shell_list`, default true).

Both keys are read only when they are named in the input; the defaults are on.

## Steps (one job, 1 node x 2 H200 as in NOTE-2026-10-09-caltech-prof, <= 1 h)

1. Build two he_star_m1 H200 binaries with exactly the options of athena_gpu_he_star_m1_98835d99_nofma
   (CUDA 12.9, HOPPER90, MPI, host -ffp-contract=off). Build them incrementally, from the same build dir if possible:
   - BASE = fork/rt-integration `b2f2e897`
   - NEW = fork/vgdspeed-1009 tip

   Report both md5.
2. Run docs/handover/deltaai-prof-1009/prof_cuda.sh (as prof_caltech.sbatch) once per binary, interleaved
   BASE, NEW, BASE, NEW. Use `ARMS="agc_plain:agc:plain:60 agc_kt10:agc:kt:10 agc_kt60:agc:kt:60 bsg_plain:bsg:plain:30
   bsg_kt5:bsg:kt:5 bsg_kt30:bsg:kt:30"` for the first pair and only the plain arms for the second pair.
   Export KOKKOS_TOOLS_TIMER_BINARY=1. Stop on rc != 0 / FATAL; the plain arm is the smoke.
3. Write NOTE-2026-10-10-caltech-vgdspeed.md with:
   - s/cycle (mean and median) for each arm and repeat;
   - the grouped table of prof_group.py (vet_gd halo, sweep, twin), with the new kernel labels m1_vgd_hl_* in
     the "vet_gd halo pack/unpack" group;
   - the `<rad_m1> vgdspeed-1009` line of each NEW log (lists made, MB, exchanges by list / dense);
   - peak GPU memory per arm.

## Expected (Raven 4 A100, 1 block per GPU; tip ec0757c1 = c67fbbb9 src)

- s/cycle: AG Car A hr 0.527 -> 0.403 and BSG reduced 0.391 -> 0.299.
- vet_gd halo share of kernel time: 42 % -> 27 % (AG Car) and 49 % -> 31 % (BSG).
- On viper (2 MI300A), AG Car A hr went 1.01 -> 0.58 s/cycle.
- Details are in /viper/ptmp2/jinma/vgdspeed_1009/VGDSPEED.md.

## Paste line for main

`Read docs/handover/TASK-2026-10-10-caltech-vgdspeed.md on fork/bsg-files-1009 and run it (one H200 job, <= 1 h).`
