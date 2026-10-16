# TASK 2026-10-15 caltech -> raven: accretor bigger domain (stream nozzle + donor sink), first test

**Test only, report back by NOTE** (docs/handover/NOTE-2026-10-1x-raven-accretor-bigdom.md on bsg-files-1009).
User decision 10-15: this test runs on Raven, not Caltech. Nothing of it has run anywhere yet.

## What
User request 10-15: make the Plaskett accretor domain larger (r_out 0.85 -> 1.5 d_L1) and inject the L1 stream
inside the domain instead of through the r_out boundary window. Design: docs/dev/accretor_rhd_design.md section 13
on the branch (read 13.1-13.4 first).

- Branch `accretor-bigdom-1015` (github jing-ze-ma/athenak), head **239a7374** (code c22b696e + docs/handover data).
  Base: accretor-rhd-1015 af45fa46 (the orbit1 production code, fdb18517 + docs).
- Build: `-D PROBLEM=ry_per_accretor`, GPU, Release, `-D CMAKE_CXX_FLAGS=-ffp-contract=off` (the "nofma" builds the
  accretor runs use). It compiles on Caltech (CUDA 12.9, H200); it has NOT been executed.
- grid_stretch.hpp changed (new log-tail slots): expect a full rebuild.

## New keys (all default off; absent = old behaviour)
| key | value in the deck | meaning |
|---|---|---|
| mesh/f_stretch_r_tail_nin, _x1max, _ij | 640, 13.7926, 520 | old 640-cell grid kept face for face to face 520 (r 10.4157), then constant-ratio (log) tail to x1max |
| mesh/x1max, nx1 (and meshblock/nx1) | 24.33996 (1.5 d_L1), 840 | cost 840/640 = 1.31x per cycle |
| problem/stream_inject | volume | stream nozzle inside the domain (no r_out window) |
| problem/stream_r | 13.7926 | nozzle at the old r_out: same ballistic state (phi 3.410 deg, v_r -126.46, v_phi 60.68) |
| problem/stream_ncell | 2 (default) | nozzle = 2 cells above the face nearest stream_r (r 13.7974-13.8708) |
| problem/outer_rho | copy | plain outflow at the new r_out (now honoured with thermo = general) |
| problem/env_r_top_max | 13.7926 | WB top capped at the old r_out (else fake envelope inside the donor lobe) |
| problem/donor_mask | roche | donor lobe (Phi < Phi(L1), x > d_L1) reset to floor ambient each stage; hst columns Mmsk/Emsk replace Mwal/Jwal |
| problem/donor_irr | (off in the base test) | lobe also holds E = a T_d^4, F = 0; problem/donor_teff default 26200 K |

## Files
docs/handover/accretor-bigdom-1015/ on the branch: `bigdom.athinput` (the deck), `column_z0.02_p381.txt`
(md5 78f4a632), `rosseland_ext2_gs98_x0.70_z0.02.txt` (ad516c7e), `planck_ext2_gs98_x0.70_z0.02.txt` (99d668f6).
Replace `@HDIR@` in the deck by the directory holding the three tables:
`sed "s#@HDIR@#$PWD#" bigdom.athinput > run.athinput`.
No restart file: the mesh differs from orbit1's, so the run starts fresh (static phase, stream off).

## Test plan (1 GPU)
1. Static phase: `athena -i run.athinput -d <dir>` (deck has stream = false, tlim 0.00946, rst at 0.00946).
   Old 640-cell static phase: 1345 cycles, dt 7.03e-6, ~700 s on 1 H200.
2. Stream phase: restart from that rst with
   `problem/stream=true time/tlim=0.025 output2/dt=0.0011829 output4/dt=0.0047317 output3/dt=0.025`.
   ROL2U keys and time/dt_retry = true are already in the deck. Old: 0.00946 -> 0.0237 took ~950 s on 1 H200.
3. Optional, if the base test passes: same stream phase from the same static rst with `problem/donor_irr=true`
   (donor_teff 26200 default), to t 0.015 is enough: check the donor light enters (Erad, Ltop rise; no FATAL).

## Pass checks
- No FATAL, no NON-CONVERGED; dt_retry events listed if any.
- Min rate (hst col 7, rate over t 0.012-0.024) = **19.795 code/time** (orbit1/start0; analytic 19.80) within 1 %.
- Stream path: same as start0 (hydro_w bin at t ~0.0189 and 0.0237: density ridge from phi ~3.4 deg at r 13.8 to
  the star; ballistic arrives at R_acc phi 16.5 deg). Old impact: hottest dense gas r 9.43, phi 15 deg near t 0.026.
- Mask stable: Mmsk (hst col 9) small and smooth, no gas piling at the lobe edge, nothing hot leaking out of it.
- dt not limited by the mask or the nozzle: dt should stay ~7.0e-6 before impact as in start0 (dt 7.03e-6); if dt
  drops, locate the limiting cell.
- Static phase: envelope unchanged vs the 640-cell run for r <= 10.4 (same faces); region r > 13.8 floor ambient
  at radiative equilibrium, no psi > 0 artefacts near the donor.
- Cost per cycle vs orbit1 approach (0.51-0.56 s/cycle on H200), expect ~1.3x on equal hardware.

## Reference numbers (orbit1 / start0, 640 cells, boundary window, fdb18517)
Min rate 19.795; dt 7.03e-6 (static and approach); Mout 4.8e-4 cumulative by t 0.1235 (orbit1 link sums restart at
0 each link); MR 0.574 by t 0.1235; impact passes with ROL2U (21 dt retries, t 0.0243-0.0244, smoke_a);
stream peak rho 1.231e-8 g/cc, T 64767 K, P_rad/P_gas 0.42.

## Your rules
Use Raven's own GPU/node rules (GPUs per node, partition/QoS, node exclusions, GPU environment exports, walltime
and scratch rules) - not Caltech's. Short test jobs only; no production chain.
