# NOTE 2026-10-02 (viper): BSG code READY -- build rt-integration 30bf6c03

Follows NOTE-2026-10-02-bsg-queue-now.md (READY-marker recipe).

## Commit
Build he_star_m1 from rt-integration **30bf6c03** (or any later rt-integration commit).
It contains the stall fix (41e8c2c7: `<rad_m1>/implicit_src_stable`, default true for fresh runs),
fixes-1002c (restart grid guard always on, `<rad_m1>/implicit_face_weight` = distance on sp) and
hse-face-1002 (2b0e740d: `<problem>/he_bc_hse_flux` = face where he_bc_outer = hse).
viper binary: ROCm 7.2, -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF, md5 90ce905a317123da73fa102103dfee2c.

## Input
Bundle bsg3d_arm2.athinput with ONE change: `<output1>` dt = 1000.0 (hst cadence <= 1158 s). time_scheme stays be.
The new defaults come from the code; do not name them.

## Check in the fresh-start restart file header (parameter dump, "Default value added at run time")
- `implicit_src_stable = 1`
- `implicit_face_weight = distance`
- `he_bc_hse_flux = face` (arm 2, he_bc_outer = hse)
- also `implicit_det_reduce = 1`, `<hydro>/rad_signal_speed = true`
A restart from an rst written by older code prints
`rad_m1: implicit_src_stable = false (restart file lacks the key; ...)` and keeps the old form: start fresh.

## Port gate: bundle bsg_col_arm2.athinput (tlim 2e4 s), viper CPU at 30bf6c03, 1 rank, gate_compare.py
`dt_end 1.243014e+02  Picard mean 2.994  L_top/L_in end 0.999466 range [0.999374, 0.999467]  FATAL 0  nan 0`
(21 hst rows, first dt 124.3031 s).  Your GPU column should match to round-off.

## 3-D smoke (viper MI300A, apudev 1 node x 2 GPUs = 2 ranks x 2 MeshBlocks, job 12059365)
Fresh, 12 cycles: dt 88.1770 s (cycle 0) -> 88.1768 s (cycle 12), Picard mean 3.0 (max 3, 0 non-converged),
1.63 s/cycle (cycles 2-12); restart leg 6 cycles: Picard 3.0, 2.16 s/cycle; no FATAL/nan.
Production on viper: 2 nodes x 2 GPUs (1 MeshBlock per GPU), chain jobs 12059375-77 (arm 2).

Report back with NOTE-2026-10-0x-<site>-bsg.md as before.
