# TASK 2026-10-15 caltech -> raven: accretor bigger domain, L1 nozzle + donor sink/irradiation, first test (v2)

**v2 supersedes v1 (a5b01dc7).** **Test only, report back by NOTE**
(docs/handover/NOTE-2026-10-1x-raven-accretor-bigdom.md on bsg-files-1009). User decision 10-15: the long test runs
on Raven; Caltech ran only the regression gate and short smokes (results below).

## What
User request 10-15: Plaskett accretor with r_out 0.85 -> 1.5 d_L1, the L1 stream injected inside the domain at L1
(Lubow-Shu nozzle), the donor's Roche lobe as a sink that radiates the donor's luminosity. Design and all numbers:
docs/dev/accretor_rhd_design.md section 13 (13.5 and 13.6 are the current state) on the branch.

- Branch `accretor-bigdom-1015` (github jing-ze-ma/athenak), head **ddb14cf8** (code 82197717). Base:
  accretor-rhd-1015 af45fa46 (orbit1 code fdb18517 + docs).
- Build: `-D PROBLEM=ry_per_accretor`, GPU, Release, `-D CMAKE_CXX_FLAGS=-ffp-contract=off` ("nofma", as the
  accretor runs). grid_stretch.hpp changed: full rebuild.

## Regression gate (done on Caltech, 1 H200): BITWISE
84b55d97 vs the orbit1 binary fdb18517, new keys absent, 20 cycles: env13, S3a-final, S5 uniform, S5 profile and
the orbit1-lineage restart (smoke_a rst 00009 + orbit1 keys) give identical hst files, identical bin files and
identical cycle/dt lines. Later commits only touch code behind the new keys.

## Deck and keys
docs/handover/accretor-bigdom-1015/ on the branch: `bigdom.athinput` + `column_z0.02_p381.txt` (md5 78f4a632),
`rosseland_ext2_gs98_x0.70_z0.02.txt` (ad516c7e), `planck_ext2_gs98_x0.70_z0.02.txt` (99d668f6).
`sed "s#@HDIR@#$PWD#" bigdom.athinput > run.athinput` in that directory.
No restart file: the mesh differs from orbit1's, the run starts fresh.

| key | deck value | meaning |
|---|---|---|
| mesh/f_stretch_r_tail_nin, _x1max, _ij | 640, 13.7926, 520 | orbit1's grid face for face to r 10.4157, then log tail |
| mesh/x1max, nx1, meshblock/nx1 | 24.33996 (1.5 d_L1), 840, 840 | 1.31x cells |
| problem/stream_inject | l1 | Lubow-Shu nozzle: 2 cells just inside L1 (r 16.1348-16.2206), |phi| < 0.2393 rad |
| problem/stream_mdot | 19.795 | band mass rate (= orbit1 Min rate); peak rho_L1 0.6859 code = 3.505e-8 g/cc |
| problem/stream_l1_t | 155300 K | uniform T_L1; flow speed c_T 45.45 km/s along -x; sigma_y 1.2942 Rsun (A 7.994, B 6.994) |
| problem/outer_rho | copy | plain outflow at r_out, dark Marshak radiation top |
| problem/env_r_top_max | 13.7926 | WB top capped at orbit1's r_out |
| problem/donor_mask | roche | lobe (Phi < Phi_L1, x > d_L1) + the sector behind the nozzle = sink; hst Mmsk/Emsk (cols 9/10) |
| problem/donor_irr, donor_lum | true, 69183.1 Lsun (10^4.84) | lobe holds E = a T_d^4, F = 0; T_d 26019.8 K over 4 pi R_L^2, R_L 12.9432 |
| hst cols 25/26 | Mnzf, Eirr | mass through the nozzle's inner face; radiation energy added by donor_irr |
(stream_inject = volume + stream_r = 13.7926 is the alternative nozzle at orbit1's r_out, not for this test.)

## Caltech smoke (g10: this deck, stream on from t = 0, no static phase, 24 min)
1563 cycles to t 0.0078: no FATAL, 1 NON-CONVERGED, no dt retry. Mnzf rate 21.5 -> 20.84 (still falling; +5 %
over 19.795), Min - Mnzf ~ 0, Mmsk 1.4e-4 total, dt 7.05e-6 -> 4.0e-6 (location not measured), 0.84 s/cycle on
H200, Ltop 6.0e5 -> 1.12e6 code (donor light), Eirr 2.3e4 by t 0.0074.

## Test plan (1 GPU)
1. Static phase: `athena -i run.athinput -d <dir>` (deck: stream = false, tlim 0.00946, rst at 0.00946).
2. Stream phase: restart from that rst with
   `problem/stream=true time/tlim=0.045 output2/dt=0.0011829 output4/dt=0.0047317 output3/dt=0.005`.
   The stream now starts at L1 (16.2 Rsun, not 13.8): the first impact is later than orbit1's t ~0.025, at about
   0.025 + 0.012 = ~0.037 (stop at 0.045 or at first impact + 0.005; rst links fine if the wall is short).
3. Report the cost per cycle and the time-to-impact.

## Pass checks
- No FATAL / NON-CONVERGED storms; list dt_retry events.
- Stream mass rate: Mnzf rate (col 25) over the steady part = 19.795 within 5 % (report the value; the Caltech smoke
  saw +5 %, falling); Min (col 7) - Mnzf small.
- L1 nozzle stable: no growing junk around r 16.1-16.3, phi +-15 deg (hydro_w bins), Mmsk smooth and small.
- Same stream path: density ridge from L1 bending to phi > 0 and hitting the star near phi 15-16.5 deg (orbit1:
  ballistic arrival phi 16.5 deg, hottest impact gas r 9.43 phi 15 deg).
- Mask stable: lobe at floor, no hot gas leaking out of it; Eirr and the donor light (Ltop) reasonable: the band
  should receive sigma T_d^4 = L_d/(4 pi R_L^2) per area at the lobe surface.
- dt not limited by the nozzle or the mask: locate the dt-setting cell after the switch-on (Caltech: 7.05e-6 ->
  4.0e-6 by t 0.006); orbit1 reference before impact 7.03e-6.
- Envelope r <= 10.4 unchanged vs orbit1 in the static phase.

## Reference numbers (orbit1 / start0, 640 cells, window at r 13.79, fdb18517)
Min rate 19.795 (every link); dt 7.03e-6 (static and approach); stream at 13.79: phi 3.41 deg, v_r -126.46, v_phi
60.68, rho 1.23e-8 g/cc, T 64.8 kK; impact passes with ROL2U (21 dt retries at t 0.0243-0.0244, smoke_a); MR 0.574
and Mout 4.8e-4 by t 0.1235 (link sums restart at 0 each link); approach 0.51-0.56 s/cycle on H200.

## Your rules
Use Raven's own GPU/node rules (GPUs per node, partition/QoS, node exclusions, GPU environment exports, walltime
and scratch rules) - not Caltech's. Short test jobs only; no production chain.
