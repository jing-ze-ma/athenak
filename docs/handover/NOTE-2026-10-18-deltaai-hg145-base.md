# NOTE DeltaAI -> Raven (2026-10-18 05:00 UTC): hg145 baseline twin B0 + background analysis of B and B0 (TASK 23dd32be)

## Verdict
**The KE growth is not a floor effect.** The baseline B0 (dfloor 1e-18, pfloor 1e-11, sponge_dmax 2e-14, mk IC) shows the SAME
KE rise in the SAME radial bands, ~0.5 d EARLIER, and it saturates (~8.4e38 erg from 2.5 d). In fl15 B the same rise starts ~2.0 d
and is still climbing at 3 d (6.3e38). In both runs >= 99 % of the KE growth after 2 d sits at **0.90-1.01 R_ph**, below and at the
photosphere, not in the floor/interface zone (1.01-1.10 R_ph carries ~1 %). This looks like the near-surface convective layer
spinning up (onset delayed ~0.5 d in fl15); floors only move the onset.

## B0 run (job 3403137, gh095, 1 GH200, ONE link 3.1 min)
Deck hg145_mk_nobh.athinput.in (md5 67fc1f7e, @FILES@ only), binary 4c494982 (md5 385860369a7d3d3801df2b4649c0ced7), keys exactly the
Test B keys. rc 0, FATAL 0, NaN 0, DIVERGED 0, summary NON-CONVERGED 0, Picard mean 3.07 max 5, 1815 cycles to t 2.592e5, DONE.
(First submit 3403136 died at launch on gh011 "srun: Error configuring interconnect", 6 s, no run; resubmitted excluding gh011.)
Dir: /work/nvme/bivj/jma20/hg145_1018/B0 (bins every 4320 s).

## KE vs time (hst; user KE_int / hydro 1-KE+2-KE+3-KE, erg)
| t (d) | B (fl15) KE_int | B KE_tot | B0 (base) KE_int | B0 KE_tot |
|---|---|---|---|---|
| 0.00 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| 0.25 | 2.80e+37 | 2.86e+37 | 5.03e+37 | 5.09e+37 |
| 0.50 | 4.72e+37 | 5.23e+37 | 8.34e+37 | 9.03e+37 |
| 0.75 | 5.18e+37 | 6.37e+37 | 9.13e+37 | 1.07e+38 |
| 1.00 | 5.89e+37 | 7.46e+37 | 1.04e+38 | 1.23e+38 |
| 1.25 | 6.75e+37 | 8.52e+37 | 1.21e+38 | 1.41e+38 |
| 1.50 | 7.48e+37 | 9.39e+37 | 1.36e+38 | 1.59e+38 |
| 1.75 | 8.07e+37 | 1.00e+38 | 1.76e+38 | 2.02e+38 |
| 2.00 | 8.83e+37 | 1.08e+38 | 3.33e+38 | 3.64e+38 |
| 2.25 | 1.11e+38 | 1.32e+38 | 5.79e+38 | 6.15e+38 |
| 2.50 | 2.22e+38 | 2.44e+38 | 7.93e+38 | 8.37e+38 |
| 2.75 | 4.19e+38 | 4.45e+38 | 8.28e+38 | 8.83e+38 |
| 3.00 | 6.25e+38 | 6.57e+38 | 8.42e+38 | 9.13e+38 |

## Which band carries the growth (bands.py: dump nearest 2.0 d and the last 3 dumps; volume-weighted <v_r>, mass-weighted rms|v|, KE)
Band sums reproduce the hst KE_tot (B 1.079e38 -> 6.566e38; B0 3.638e38 -> 9.134e38), so the spherical volume weights are right.
| band r/R_ph | B: KE 2.0 -> 3.0 d (share of dKE), rms|v| km/s | B0: KE 2.0 -> 3.0 d (share), rms|v| km/s |
|---|---|---|
| 0.00-0.90 | ~1e33, 0 % | ~1e33, 0 % |
| 0.90-0.99 | 6.3e37 -> 3.9e38 (60 %), 0.19 -> 0.48 | 1.3e38 -> 4.8e38 (64 %), 0.27 -> 0.53 |
| 0.99-1.01 | 4.5e37 -> 2.6e38 (39 %), 0.29 -> 0.70 | 2.4e38 -> 4.3e38 (35 %), 0.66 -> 0.90 |
| 1.01-1.03 | 5.0e35 -> 4.1e36 (1 %), 0.39 -> 0.96 | 3.2e35 -> 6.0e36 (1 %), 0.32 -> 1.13 |
| 1.03-1.10 | 5.7e33 flat, rms 1.36 -> 1.70 | 1.2e32 -> 2.9e32, rms 0.18 -> 0.27 |
| 1.10-1.50 | 2e30 flat | 2e31 flat |
Shell-mean <v_r> at 3 d: B 0.90-0.99 -0.14, 0.99-1.01 -0.19, 1.01-1.03 -0.33, 1.03-1.10 -0.42 km/s;
B0 -0.21, -0.34, +0.57, -0.01 km/s.
Difference worth noting: in fl15 B the thin 1.03-1.10 R_ph layer moves faster (rms 1.7 vs 0.27 km/s) and bg145 shows
max|rho/rho0-1| up to 1.58 there (B0: <= 0.29). That is the floor-zone signature, but it holds negligible KE (~6e33 erg).
B also has 10x less sponge power (int P_sponge/int L_top 1.7e-9 vs 1.7e-8 in B0) and 10x lower Mdot_top (8.1e14 vs 8.1e15).

## bg145.py full stdout (`BASE=hegiant LOGF="run.*.log" VISPY=<wt>/vis/python python bg145.py 7.95704e12 2.22652e13 <run> 4`)
Files: docs/handover/hg145-prod-1018/deltaai-base/{bg_B.txt,bg_B0.txt,bands_B.txt,bands_B0.txt,bands.py}
### B
```
#### B  rows 865  t_end 2.592e+05  Picard 3.04/5 NONCONV 0  fofc at-floor max 0 other max 0
  hst t 4.286e+02 L_top/L0 0.92497 Mdot_top 9.09e+14 Min_top 0 P_sponge/L 2.80e-10 (L_bot-L_in)/L_in -1.19e-02 v1sq_wall 5.48e+02 Mdot_bot 0
  hst t 3.243e+04 L_top/L0 1.14454 Mdot_top 8.68e+14 Min_top 0 P_sponge/L 1.75e-09 (L_bot-L_in)/L_in -1.16e-02 v1sq_wall 3.93e+03 Mdot_bot 0
  hst t 6.486e+04 L_top/L0 1.14581 Mdot_top 8.6e+14 Min_top 0 P_sponge/L 1.74e-09 (L_bot-L_in)/L_in -1.13e-02 v1sq_wall 4.37e+03 Mdot_bot 0
  hst t 1.297e+05 L_top/L0 1.14704 Mdot_top 8.43e+14 Min_top 0 P_sponge/L 1.74e-09 (L_bot-L_in)/L_in -1.08e-02 v1sq_wall 3.87e+03 Mdot_bot 0
  hst t 1.944e+05 L_top/L0 1.15125 Mdot_top 8.27e+14 Min_top 0 P_sponge/L 1.73e-09 (L_bot-L_in)/L_in -1.03e-02 v1sq_wall 4.62e+03 Mdot_bot 0
  hst t 2.592e+05 L_top/L0 1.14408 Mdot_top 8.08e+14 Min_top 0 P_sponge/L 1.76e-09 (L_bot-L_in)/L_in -9.71e-03 v1sq_wall 6.01e+03 Mdot_bot 0
  int P_sponge dt / int L_top dt = 1.74e-09
  bands (r/R_ph): 0.00-0.20 | 0.20-0.50 | 0.50-0.90 | 0.90-0.99 | 0.99-1.01 | 1.01-1.03 | 1.03-1.10 | 1.10-1.50 | 1.50-2.00 | 2.00-2.80 | 2.80-3.11
  c0     t 0.000e+00    +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000  | Mdot(R_ph)/sr +0.00e+00  sponge-edge dv +0.00 km/s d(e ratio) +0.0e+00  min rho 1.0e-15
  c152   t 2.171e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   +0.01  0.002 |   -0.54  0.162 |   -0.20  0.009 |   +0.45  0.623 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +3.82e+20  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.8e-16
  c212   t 3.029e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.08  0.011 |   -0.45  0.148 |   -0.27  0.016 |   +0.49  0.807 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +1.53e+21  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.7e-16
  c424   t 6.057e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.43  0.056 |   -0.42  0.171 |   -0.52  0.047 |   +0.70  1.255 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +1.48e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.6e-16
  c484   t 6.914e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.42  0.058 |   -0.39  0.178 |   -0.58  0.057 |   +0.76  1.345 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +2.04e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.7e-16
  c636   t 9.086e+04    -0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   -0.45  0.066 |   +0.45  0.178 |   -0.73  0.089 |   -0.91  1.505 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +2.69e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.5e-16
  c696   t 9.943e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.47  0.070 |   +0.50  0.178 |   -0.79  0.105 |   -0.99  1.546 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +2.61e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.5e-16
  c756   t 1.080e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.48  0.074 |   +0.52  0.179 |   -0.84  0.122 |   -1.06  1.575 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +2.29e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.4e-16
  c908   t 1.297e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.49  0.080 |   +0.62  0.181 |   -0.97  0.171 |   -1.25  1.580 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +5.49e+21  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.4e-16
  c968   t 1.383e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.49  0.082 |   +0.63  0.203 |   -1.02  0.192 |   -1.32  1.556 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -7.89e+20  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.4e-16
  c1180  t 1.686e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.48  0.093 |   +0.74  0.328 |   -1.19  0.276 |   -1.55  1.344 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.29e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.1e-04  min rho 9.5e-16
  c1240  t 1.771e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.47  0.099 |   +0.76  0.369 |   -1.24  0.301 |   -1.61  1.259 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.44e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.5e-16
  c1392  t 1.989e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.47  0.114 |   +0.83  0.481 |   -1.35  0.367 |   -1.73  1.024 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.63e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.6e-16
  c1452  t 2.074e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.46  0.120 |   +0.84  0.527 |   -1.39  0.404 |   -1.78  0.958 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.70e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.2e-04  min rho 9.6e-16
  c1512  t 2.160e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.46  0.125 |   +0.89  0.576 |   -1.44  0.465 |   -1.82  0.913 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.76e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.6e-16
  c1664  t 2.377e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.45  0.146 |   +0.96  0.713 |   -1.53  0.672 |   -1.95  0.780 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.86e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.6e-16
  c1724  t 2.463e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.45  0.158 |   +0.96  0.767 |   -1.56  0.750 |   -2.02  0.740 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.91e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.0e-04  min rho 9.6e-16
  c1815  t 2.592e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.44  0.194 |   +0.91  0.807 |   -1.61  0.858 |   -2.16  0.724 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.96e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.2e-04  min rho 9.6e-16
```
### B0
```
#### B0  rows 865  t_end 2.592e+05  Picard 3.07/5 NONCONV 0  fofc at-floor max 0 other max 0
  hst t 4.286e+02 L_top/L0 0.98384 Mdot_top 8.81e+15 Min_top 0 P_sponge/L 2.43e-09 (L_bot-L_in)/L_in -1.19e-02 v1sq_wall 5.48e+02 Mdot_bot 0
  hst t 3.243e+04 L_top/L0 1.12619 Mdot_top 8.69e+15 Min_top 0 P_sponge/L 1.70e-08 (L_bot-L_in)/L_in -1.16e-02 v1sq_wall 3.93e+03 Mdot_bot 0
  hst t 6.486e+04 L_top/L0 1.12855 Mdot_top 8.6e+15 Min_top 0 P_sponge/L 1.69e-08 (L_bot-L_in)/L_in -1.13e-02 v1sq_wall 4.37e+03 Mdot_bot 0
  hst t 1.297e+05 L_top/L0 1.13212 Mdot_top 8.44e+15 Min_top 0 P_sponge/L 1.68e-08 (L_bot-L_in)/L_in -1.08e-02 v1sq_wall 3.87e+03 Mdot_bot 0
  hst t 1.944e+05 L_top/L0 1.12540 Mdot_top 8.28e+15 Min_top 0 P_sponge/L 1.70e-08 (L_bot-L_in)/L_in -1.03e-02 v1sq_wall 4.62e+03 Mdot_bot 0
  hst t 2.592e+05 L_top/L0 1.13154 Mdot_top 8.07e+15 Min_top 0 P_sponge/L 1.70e-08 (L_bot-L_in)/L_in -9.71e-03 v1sq_wall 6.01e+03 Mdot_bot 0
  int P_sponge dt / int L_top dt = 1.70e-08
  bands (r/R_ph): 0.00-0.20 | 0.20-0.50 | 0.50-0.90 | 0.90-0.99 | 0.99-1.01 | 1.01-1.03 | 1.03-1.10 | 1.10-1.50 | 1.50-2.00 | 2.00-2.80 | 2.80-3.11
  c0     t 0.000e+00    +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   +0.00  0.000  | Mdot(R_ph)/sr +0.00e+00  sponge-edge dv +0.00 km/s d(e ratio) +0.0e+00  min rho 1.0e-14
  c152   t 2.171e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   +0.01  0.002 |   -0.67  0.120 |   +0.02  0.000 |   -0.05  0.093 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +1.62e+21  sponge-edge dv +0.00 km/s d(e ratio) -5.9e-05  min rho 1.0e-14
  c212   t 3.029e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.09  0.013 |   -0.62  0.094 |   +0.03  0.001 |   +0.07  0.112 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +2.91e+21  sponge-edge dv +0.00 km/s d(e ratio) -6.0e-05  min rho 1.0e-14
  c424   t 6.057e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.58  0.077 |   -0.57  0.177 |   +0.08  0.009 |   +0.10  0.152 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +1.80e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.1e-05  min rho 1.0e-14
  c484   t 6.914e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.58  0.079 |   -0.53  0.186 |   +0.10  0.013 |   +0.11  0.159 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.000  | Mdot(R_ph)/sr +2.49e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.1e-05  min rho 1.0e-14
  c636   t 9.086e+04    -0.00  0.000 |   +0.00  0.000 |   +0.00  0.000 |   -0.60  0.088 |   +0.50  0.184 |   +0.14  0.025 |   +0.14  0.173 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +3.22e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.1e-05  min rho 1.0e-14
  c696   t 9.943e+04    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.63  0.091 |   +0.55  0.184 |   +0.15  0.031 |   +0.15  0.178 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +3.08e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.1e-05  min rho 1.0e-14
  c756   t 1.080e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.65  0.093 |   +0.56  0.185 |   +0.17  0.038 |   +0.17  0.182 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +2.61e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.2e-05  min rho 9.9e-15
  c908   t 1.297e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.66  0.101 |   +0.66  0.213 |   +0.22  0.061 |   +0.20  0.193 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr +2.95e+21  sponge-edge dv +0.00 km/s d(e ratio) -6.2e-05  min rho 9.9e-15
  c968   t 1.383e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.66  0.106 |   +0.67  0.240 |   +0.24  0.073 |   +0.22  0.198 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -4.91e+21  sponge-edge dv +0.00 km/s d(e ratio) -6.8e-05  min rho 9.9e-15
  c1180  t 1.686e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.64  0.134 |   +0.81  0.340 |   +0.33  0.123 |   +0.27  0.214 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -1.96e+22  sponge-edge dv +0.00 km/s d(e ratio) -2.5e-05  min rho 9.9e-15
  c1240  t 1.771e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.64  0.143 |   +0.85  0.371 |   +0.35  0.140 |   +0.28  0.219 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.16e+22  sponge-edge dv +0.00 km/s d(e ratio) -9.5e-05  min rho 9.9e-15
  c1392  t 1.989e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.63  0.166 |   +0.96  0.452 |   +0.47  0.207 |   +0.31  0.233 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.52e+22  sponge-edge dv +0.00 km/s d(e ratio) -1.0e-04  min rho 9.9e-15
  c1452  t 2.074e+05    -0.00  0.000 |   -0.00  0.000 |   +0.00  0.000 |   -0.62  0.190 |   +0.98  0.502 |   +0.57  0.261 |   +0.32  0.239 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.63e+22  sponge-edge dv +0.00 km/s d(e ratio) -5.1e-05  min rho 9.9e-15
  c1512  t 2.160e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.62  0.222 |   +1.05  0.534 |   +0.73  0.340 |   +0.33  0.245 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.72e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.6e-05  min rho 9.9e-15
  c1664  t 2.377e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.61  0.288 |   +1.15  0.661 |   +1.13  0.596 |   +0.36  0.263 |   -0.05  0.001 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.91e+22  sponge-edge dv +0.00 km/s d(e ratio) -1.1e-04  min rho 9.9e-15
  c1724  t 2.463e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.65  0.303 |   +1.14  0.716 |   +1.23  0.689 |   +0.37  0.271 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.95e+22  sponge-edge dv +0.00 km/s d(e ratio) -6.5e-05  min rho 9.9e-15
  c1815  t 2.592e+05    -0.00  0.000 |   -0.00  0.000 |   -0.00  0.000 |   -0.67  0.339 |   +1.06  0.747 |   +1.30  0.809 |   +0.38  0.287 |   -0.05  0.002 |   -0.02  0.000 |   -0.02  0.000 |   -0.02  0.001  | Mdot(R_ph)/sr -2.89e+22  sponge-edge dv +0.00 km/s d(e ratio) +1.9e-04  min rho 9.9e-15
```

Nothing running on DeltaAI. Cost: B0 ~0.05 GPU-h.
