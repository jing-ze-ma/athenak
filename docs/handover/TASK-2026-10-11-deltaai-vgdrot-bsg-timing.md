# TASK raven -> DeltaAI: vet_gd rotation-spike speed-up, BSG timing + bitwise (1 short ghx4-interactive job)

The Raven gpudev slot for this is ~35 min out (the gpu partition is ~1 day). **First result wins.** Same environment rules as
your AG Car A tasks (MPICH_GPU_SUPPORT_ENABLED=1, plain srun, 1 rank per GPU).

## Code (fork), PROBLEM=he_star_m1, build_inc_deltaai.sh hes_gpu
- BASE: `rt-integration` **d28dacd8** (production code).
- NEW: `vgdrot-1011` **a749711a**. vet_gd: the per-shell exact halo-list builds are batched (CUB scans, fewer collectives), the next rotation's
  tables are prefetched, the per-cycle allocations are removed, and the per-shell sweep kernels take kernel-argument launches. Opt-in
  vet_gd_halo_ipc modes stay OFF. Intended to be bitwise.

## Case: BSG true-repro hr, reduced to 4 blocks (1 node x 4 GH200), as your BSG gates
Production BSG input (your bsg_hrdet_1009 / Raven bsg.athinput set) plus:
`mesh/nx2=128 mesh/nx3=128 mesh/x2min=1.3089969389957472 mesh/x2max=1.832595714594046 mesh/x3max=0.5235987755982988 meshblock/nx2=64 meshblock/nx3=64`
32 cycles, all dumps off except hst every cycle (`time/ndiag=1`), plus one bin dump at the last cycle for the bitwise check.
Run BASE, NEW, BASE, NEW interleaved in one job (2 reps each).

## Report: NOTE-2026-10-11-deltaai-vgdrot-bsg-timing.md
- s/cycle per run: steady (exclude the rotation cycles c % vet_gd_rotate_every == 0 and the last), rotation cycles, average.
- NEW vs BASE: hst bitwise yes/no, last bin dump bitwise yes/no. Rep-to-rep bitwise for each binary.
- GPU-h.

## WITHDRAWN (raven)
Raven gpudev job 31067404 finished first (BSG base -> a749711a: 0.2072 -> 0.1653 s/cycle average, bitwise). No need to run this; if you already have a result, post it anyway as a GH200 cross-check.
