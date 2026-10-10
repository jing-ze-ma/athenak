# TASK viper -> Caltech: xthinfix-1009 ARM MATRIX, CPU nodes + 1-2 H200 (user 10-10)

This replaces TASK-2026-10-10-caltech-steep-battery.md, whose arms are superseded. Short jobs only; do not touch
production. The viper login node is saturated, so this matrix runs on Caltech.

## Code

Fork branch `xthinfix-1009` commit `413b38af`. It is rt-integration 5b304cf9 (vgdspeed) plus the half-range face-flux
work.

Build:
- the CPU binary (no PROBLEM, MPI, Release): `ATHENA_CPU`;
- H200 he_star_m1 as your AG Car binaries: `ATHENA_GPU`;
- rt-integration `5b304cf9` he_star_m1 H200: `ATHENA_BASE` (reuse the vgdspeed build if it is that commit).

Report the md5s.

Keys used by the arms. A key absent from a FRESH input is FATAL, so the bundled inputs carry them all.
- `implicit_blend_xthin_mode = all | steep | kn`:
  - `steep` is the gate X^8/(X^8+X0^8) on every face.
  - `kn` is the override weighted by the face radiation Knudsen number Kn = |grad E|/(chi E), with Kn0 =
    `implicit_blend_xthin_r0`.
- `implicit_hr_recon = plm`: van Leer half-range face states.
- `implicit_hr_recon_qs = eps`: plm only where |E^k - E^n| <= eps E^n.
- `implicit_hr_recon_qs_rel`: no plm on the vacuum side of a front.
- `implicit_hr_pos = kill | floor`.
- `implicit_accel = anderson`.

## Arms (`docs/handover/xthinfix-1009/arms/ana/od.py`, `battery_arms.sh`)

| arm | meaning |
| --- | --- |
| cen | central |
| hr | production (mode all) |
| hrx0 | xthin 0 |
| sq | steep, X0 60, + plm + qs 0.1 |
| sqf | sq with pos floor |
| sqv | sq with qs_rel 1e-3 |
| sqa | sq with Anderson |
| aq | mode all + plm + qs 0.1 |
| kn01 / kn03 / kn1 | Knudsen switch, Kn0 = 0.1 / 0.3 / 1, dc |
| knp | kn03 + plm + qs 0.1 |

## Steps

1. **CPU battery.**

       B=<checkout of this branch>/docs/handover/xthinfix-1009/arms
       REPO=<xthinfix-1009 413b38af checkout> ATHENA_CPU=... XTF_RUN=<dir> P=<cores> bash $B/battery_arms.sh all

   You can split it as `order` | `gates` | `beams`, then `ana`. About 1700 one-core runs; the 512-cell levels
   dominate.
   - **order:** grey atmosphere (steady, cfl 1e4) 32-512 for all arms, vs Hopf and as per-region
     self-convergence by tau and by Kn (`ana_atm_self.py`). Pulses kappa 0.128 / 12.8 / 1280 and rw_t10 / rw_t1000,
     be and hesdirk2, in X (fixed dt), T and C (fixed c dt/dx). The relaxed criterion asks for convergence at fixed CFL.
   - **gates:** G1 / G5 / G3 for hr, sq, aq, kn03, knp.
   - **beams:** xb20, ba0, ba20, cyl, shd3b, 100 cycles np 1. Record wall, Picard, inner iterations, NON-CONVERGED,
     plm positivity fallbacks and floor clips.
   - Output: `$XTF_RUN/RESULTS/*.txt`.
2. **GPU point.** `arms/gpu_arms.sh` on 2 H200, one job (~25 min). Arms base / all / sq / aq / kn03 / knp / base on AG
   Car A (`agcar_rcxA_ge_accel_st_mg_hr.athinput` from your agcarb_1009/files), 30 cycles each. It reports s/cycle,
   Picard, inner iterations, base-vs-all bin-data identity, and all vs each arm by radius. Smoke one arm first
   (nlim 3).

## Viper results so far (the reference)

Grey atmosphere, per-region self-convergence:
- Arms with first-order override in the thin top (hr, beam, beam_kn, steep60 dc) do NOT converge in Kn > 0.3 and
  diverge in 0.03 < Kn < 0.3 at 512.
- sq converges everywhere: thin 0.95, transition 1.5, thick 1.6. Hopf L1 3.9e-4 at 512.

Thin pulse and rw_t10:
- sq: thin pulse k0.128 = central; rw_t10 converges.
- hr / aq: the thin pulse converges at fixed CFL (C mode 1.72) but not at fixed dt (X mode). rw_t10 fails in both X
  and C.

Beam tests (Picard passes / positivity events, against hr 4.0 / 5.0):
- sq: xb20 65 passes, cyl 9.8 passes. The E <= 0 events happen every pass and are pre-existing (cold absorbers under
  dc as well).
- aq: xb20 7.8 passes, cyl 34 passes.
- sqf: DIVERGES.
- sqv: identical to sq.

## Report

Push `docs/handover/NOTE-2026-10-10-caltech-arm-matrix.md`: md5s, the RAW `RESULTS/*.txt` (incl.
`atm_self_regions.txt`) and `RESULTS_gpu_arms.txt`, any FATAL / NON-CONVERGED, and a verdict of 5 lines at most.
