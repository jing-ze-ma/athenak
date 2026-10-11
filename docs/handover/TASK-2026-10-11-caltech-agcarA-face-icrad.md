# TASK raven -> Caltech: AG Car A GPU A/B, face boundary + IC radiation start (H200, short)

This comes from the Raven session that took over viper's work. Raven's A100 queues (gpudev and gpu) are blocked on
priority tonight. The same tests are queued on Raven as one combined job. **Whichever site finishes first wins**: post
your NOTE as soon as you have results; Raven will post a NOTE when its own run lands, and then this task is closed.
If your H200 queue is still about a day out, reply "SKIP" in the NOTE rather than queueing.

## Code (fork github.com/jing-ze-ma/athenak)
- A: branch `deepatm-1011` **e84d2513**. `rad_m1/vet_bc_diff = face` is now the default for vet_sc, vet_gd, vet_col_lat
  and the vet_col core (the diffusion intensity at the inner/bottom boundary is taken at the face, not the cell);
  `cell` reproduces the old behaviour bitwise.
- B: branch `icrad-1011` **7045bc14**. New opt-in `rad_m1/ic_rad = fs | relax | fs_relax` (fresh starts only):
  starts the thin envelope from a radiation field consistent with the VET transport instead of the IC's Eddington-like E.
  Default off, bitwise.
- Build the he_star_m1 GPU binary for each, as for your AG Car B binary (HOPPER90, MPI, CUDA 12.x).

## Input
AG Car A, production RT keys (be, vet_col + vet_gd), with files (md5):
- IC: ic_agcar_A_ge.txt, md5 4ab9e3f3dd2539fe12eb029ba9188fca
- Rosseland: rosseland_ext2_gs98_x0.36_z0.02.txt, md5 b28e97c5b5cac3afe0ea3d3fc60b7cfd
- Planck: planck_ext2orion_gs98_x0.36_z0.02.txt, md5 493f9d014bbc7d7ee0c6abb813330244
- Input deck: Raven's agcar_A_va.athinput (= agcar_A_xtf + named keys + hst every cycle). If you do not have the AG Car A
  files with these md5s, say "NEEDFILES" in your NOTE; the user will transfer them. Do not substitute other tables.

## Runs (1 node, 4 GPUs, 1 block per GPU, as Raven: 4 blocks of 480x64x64)
Fixed dt: `rad_m1/implicit_cfl=1251.3` (dt = 1239.996 s; dt/2 runs use 625.65 and twice the cycles), hst every cycle.
1. Binary A, 50 cycles: `rad_m1/vet_bc_diff=cell` vs `rad_m1/vet_bc_diff=face`, plus a repeat of `face` (run-to-run check).
2. Binary B, 60 cycles: default IC vs `rad_m1/ic_rad=fs`, `ic_rad=relax`, `ic_rad=fs_relax`.
   Then default IC and the best ic_rad also at dt/2 (120 cycles).

## Report: NOTE-2026-10-11-caltech-agcarA-face-icrad.md
- Every run: rc, FATAL, NON-CONVERGED, Picard mean/max, s/cycle (exclude rotation cycles 10, 20, ... and the last).
- Runs 1: max relative difference face vs cell of the shell means of E, F1, dens, velx vs radius (bands < 0.8,
  0.95-1.05, > 1.05 R_ph, R_ph = 388.3 Rsun) at cycle 50, and of L_top in the hst. Repeat: bitwise or not.
- Runs 2: L_top per cycle for cycles 0-60 (table or plot) for each IC option; the dt vs dt/2 L_top difference over
  cycles 1-30 for default IC vs the best ic_rad; the cost of the IC step.
- Ledger line as usual.
