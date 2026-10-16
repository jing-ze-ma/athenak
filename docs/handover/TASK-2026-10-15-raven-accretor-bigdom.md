# TASK 2026-10-15 caltech -> raven: accretor bigger domain, pinned-rate stream nozzle + donor sink/irradiation (v3)

**v3 supersedes v2 (3b8315fc, 3ebb49a2) and v1.** **Test only, report back by NOTE**
(docs/handover/NOTE-2026-10-1x-raven-accretor-bigdom.md on bsg-files-1009). User decision 10-15: the long test runs
on Raven; Caltech ran the regression gate and short smokes (below).

## What
Plaskett accretor with r_out 0.85 -> 1.5 d_L1; the stream injected inside the domain by a 2-cell nozzle at r 14.5
(inside L1, where it is supersonic) with its mass rate PINNED at stream_mdot = 19.795 (orbit1's Min rate); the
donor's Roche lobe is a sink that radiates the donor's luminosity. Design and numbers: docs/dev/accretor_rhd_design.md
section 13 (13.7 is the current nozzle; 13.5 L1 nozzle is NOT used) on the branch.

- Branch `accretor-bigdom-1015` (github jing-ze-ma/athenak), head **efe33e6a** (code 7fb237a8).
- Build: `-D PROBLEM=ry_per_accretor`, GPU, Release, `-D CMAKE_CXX_FLAGS=-ffp-contract=off` ("nofma"). grid_stretch.hpp
  changed: full rebuild.

## Regression gate (Caltech): BITWISE
84b55d97 vs the orbit1 binary fdb18517, new keys absent, 20 cycles: env13, S3a-final, S5 uniform, S5 profile, orbit1
lineage restart: identical hst, bin (whole files) and cycle/dt lines. Later commits touch only code behind new keys.

## Deck
docs/handover/accretor-bigdom-1015/ on the branch: `bigdom.athinput` + `column_z0.02_p381.txt` (md5 78f4a632),
`rosseland_ext2_gs98_x0.70_z0.02.txt` (ad516c7e), `planck_ext2_gs98_x0.70_z0.02.txt` (99d668f6).
`sed "s#@HDIR@#$PWD#" bigdom.athinput > run.athinput`; add `dt_report = true` under <problem> (prints the hydro
dt-setting cell at every hst output). Fresh start (mesh differs from orbit1; no restart file).

| key | value | meaning |
|---|---|---|
| mesh/f_stretch_r_tail_nin, _x1max, _ij; x1max; nx1 | 640, 13.7926, 520; 24.33996; 840 | orbit1 grid to r 10.4157, log tail to 1.5 d_L1 |
| problem/stream_inject, stream_r | volume, 14.5 | nozzle cells r 14.5106-14.5878, inner face at 14.5 (Mach |v|/c_tot 2.20, |v_r|/c_gas 2.36) |
| problem/stream_mdot | 19.795 | pins the face rate (peak rho 0.2726 code from the discrete face sum) |
| problem/env_t_stream, stream_width | 72966 K, 1.720 | one-zone state at 14.5 along the true path (stream_path_1015.py) |
| problem/stream_back_sink, floor_keep_momentum | default true, true | sector behind the nozzle is sink; floor keeps momentum, |v| capped |
| problem/outer_rho, env_r_top_max | copy, 13.7926 | outflow r_out (dark Marshak top); WB top capped |
| problem/donor_mask, donor_irr, donor_lum | roche, true, 69183.1 Lsun | lobe sink (hst Mmsk/Emsk cols 9/10), E = a T_d^4, T_d 26019.8 K |
| hst cols 25/26 | Mnzf, Eirr | mass through the nozzle's inner face; donor radiation energy added |

## Caltech smoke (g12: this deck + dt_report, stream on from t = 0, 1 H200)
To t 0.0021: Mnzf rate 19.802 -> 19.798 per 1e-4 (pinned to 0.04 %), Min - Mnzf constant (nozzle fill), Mmsk 4.6e-6,
no DT_RETRY, no NON-CONVERGED, dt 7.047e-6 flat, set by an envelope cell (r 8.787, at rest), not the nozzle or mask.

## Test plan (1 GPU)
1. Static phase: `athena -i run.athinput -d <dir>` (stream = false, tlim 0.00946, rst at 0.00946).
2. Stream phase: restart from that rst with
   `problem/stream=true time/tlim=0.036 output2/dt=0.0011829 output4/dt=0.0047317 output3/dt=0.005`.
   First impact: flight 14.5 -> 13.79 is 0.0065 code, then as orbit1 (0.0151 from 13.79 to impact): t ~0.031.
   Run to first impact + ~0.005 (tlim 0.036); chain short links from the rst if your wall is short.
3. Report cost per cycle and the time of first impact.

## Pass checks
- No FATAL / NON-CONVERGED storms; list dt_retry events (orbit1 needed ROL2U + 21 retries at its impact).
- **Pinned rate: Mnzf rate (col 25) = 19.795 within 0.5 % from the first output after the switch-on**; Min (col 7) -
  Mnzf constant after the fill.
- Nozzle stable: no growing junk around r 14.4-14.7, |phi - 2 deg| < 25 deg (hydro_w bins); Mmsk small, smooth.
- Same stream path as orbit1: arrives at the star near phi 15-16.5 deg (orbit1 ballistic 16.5 deg, hottest impact gas
  r 9.43 phi 15 deg).
- Mask stable, lobe at floor; donor light (Ltop, Eirr) sensible: sigma T_d^4 per area at the lobe surface.
- dt not set by the nozzle or the mask (dt_report lines); reference before impact 7.03-7.05e-6.
- Envelope r <= 10.4 unchanged vs orbit1 in the static phase.

## Reference numbers (orbit1 / start0, 640 cells, boundary window at 13.79, fdb18517)
Min rate 19.795; dt 7.03e-6 before impact; stream at 13.79 phi 3.41 deg, v_r -126.46, v_phi 60.68; impact at
t ~0.024-0.026 with 21 dt retries (t 0.0243-0.0244); MR 0.574, Mout 4.8e-4 by t 0.1235 (link sums restart per link);
0.51-0.56 s/cycle on H200 during the approach.

## Your rules
Use Raven's own GPU/node rules (GPUs per node, partition/QoS, node exclusions, GPU environment exports, walltime
and scratch rules) - not Caltech's. Short test jobs only; no production chain.
