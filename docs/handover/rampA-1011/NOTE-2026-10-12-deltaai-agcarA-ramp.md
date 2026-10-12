# NOTE DeltaAI: AG Car A full-size 3D through the frozen-MLT ramp (TASK-2026-10-12-deltaai-agcarA-ramp) -- DONE: FATAL at w = 0.11, collapses like B

**STARTED**: job **3363376** RUNNING on gh138 since 20:52 CDT (fresh from t = 0), link 1 of at most 2 (each <= 1 h, -t 00:55:00).
- Binary prodfinal-1011 **0197f020** exactly: `athena_hes_gpu_0197f0202e10`, md5 d36ed4ff9aed3bc092ab03e42d425c93.
- Files md5 4ab9e3f3 / b28e97c5 / 493f9d01 OK; deck = agcar_A_ramp.athinput.in with only @FILES@ filled (no other key changed).
- 1 node x 4 GH200, 4 ranks, 1 block/GPU. ghx4-interactive allows 1 submitted job per user, so link 2 (if needed) and the
  consfix-c2drag job are submitted one after another as each job ends.

## FATAL during the ramp (posted 21:45 CDT; analysis follows)
**Job 3363376 stopped with rc 134 at cycle 6695, t = 7.1827e6 s (ramp weight w = 0.114): Picard solve DIVERGED on all 4 ranks**
```
<rad_m1> Picard NON-CONVERGED after 200 passes: resid=2.001984e-02 (tol 1.000000e-08) lin_resid=3.331631e-08 worst cell of rank 0 (m,k,j,i)=(0,66,25,456)
### FATAL ERROR in .../src/rad_m1/rad_m1_implicit.cpp at line 13526
<rad_m1> Picard solve DIVERGED: resid=2.001984e-02 lin_resid=3.331631e-08 > implicit_resid_fatal=1.000000e-02 at time 7.182736e+06 cycle 6695 (rank 0)
```
(worst cell i = 456 of 480, i.e. near the top of the envelope; NC-local on rank 2 at (0,18,26,344).) STOP rule applied: no link 2, no rerun.
- Before that: **hesdirk2 stage 1 NOT ADMISSIBLE** first at cycle 5827 (t ~ 6.28e6 s, w ~ 0.41), then more often: 104 lines in total
  (per 100 cycles: 6 / 5 / 4 / 4 / 15 / 17 / 26 / 27 for cycles 5800-6699).
- **Cost rose 5x through the ramp**: s/cycle per 100 cycles 0.20 (c500), 0.26 (c2000), 0.34 (c3000), 0.40 (c4500), 0.32 (c5500),
  0.49 (c6000), 1.01 (c6500), ~1.5 at the end. dt 1252 s to t ~ 2e6 s, 940 s at 3.6e6, 990-1080 s through the ramp.
- fofc_report lines: 1995, with at-floor > 0 in none.
- Reached: cycle 6650 output at t = 7.134e6 s; newest rst agcar3d.00014.rst (t = 7.0e6 s); 15 bin dumps (every 5e5 s, t = 0 ... 7.0e6).
- The task's minimum (t >= 6.5e6 s) is reached; tlim 8.064e6 s is not.
- Wall 43.8 min x 4 GH200 = 2.9 GPU-h.
Next: hstA.py / reduceA.py / fluxA.py on the 15 dumps (running now on DeltaAI).

## FINAL (22:05 CDT)
**Verdict: AG Car A collapses through the ramp in full-size 3D, as B does.** Resolved convection never takes over the MLT flux
(L_conv <= 0.02 L0 at every radius and time); as w goes 1 -> 0.16 the total L at 0.55-0.90 R_ph falls to 0.30-0.54 L0 and
L_top/L0 0.96 -> 0.55; the envelope falls in (band <v_r> -23 / -30 / -34 km/s at 0.70-0.95 R_ph at t = 7e6 s) and the
photosphere contracts (1.14 -> 1.09 R_ph). The implicit solve gets harder all the way (Picard ~5 -> 42 in the hst,
stage-1 fallbacks from w ~ 0.41) and diverges at w = 0.114 (FATAL above).

### Job
- 3363376, gh138 (1 node x 4 GH200, 1 block/GPU), rc 134 (FATAL x4, Picard DIVERGED), NaN 0, wall 43.8 min, 6695 cycles.
- s/cycle: 0.20-0.40 to w = 1 (t <= 4.0e6), then rising to ~1.5; mean 2514 s / 6650 cycles = 0.378 s/cycle,
  **2.08e7 zone-cycles/s** (7.86e6 cells).
- No end-of-run Picard line (crash). Picard from user.hst: mean / max **11.72 / 69 (hst rows every 1e4 s)**.
- `NOT ADMISSIBLE` (stage fallback) lines: 104 (first at cycle 5827); `time2_pred_cap ... limited` lines: 0 (no end summary);
  fofc_report lines with at-floor > 0: 0 of 1995.
- GPU-h: **2.9**.

### hstA.py (verbatim; the run ends at 7.18e6 s, so the rows stop at 7.0e6)
```
/work/nvme/bivj/jma20/rampA_1011/run t_end 7.1806e+06 rows 719
      t     w_mlt  Ltop/L0  Lmid/Lin  LMLT/Lin  Lbot/Lin  KE_int    dt     Picard  Mdot_top
  5.010e+05  1.000  1.0007  0.5567   0.8069   1.0005  7.34e+42    1252  3.19  0.00e+00
  1.001e+06  1.000  1.0016  0.5578   0.8069   1.0012  2.57e+43    1252  3.12  0.00e+00
  1.500e+06  1.000  1.0017  0.5577   0.8069   1.0020  5.52e+43    1252  4.81  0.00e+00
  2.001e+06  1.000  0.9984  0.5544   0.8069   1.0023  9.59e+43    1252  5.56  0.00e+00
  2.501e+06  1.000  0.9927  0.5501   0.8069   1.0027  1.47e+44    1252  5.00  0.00e+00
  3.001e+06  1.000  0.9888  0.5500   0.8069   1.0038  2.07e+44    1095  8.28  0.00e+00
  3.500e+06  1.000  0.9779  0.5438   0.8069   1.0237  2.69e+44     942  9.40  0.00e+00
  4.001e+06  1.000  0.9609  0.5388   0.8069   1.1013  3.05e+44     942  11.73  0.00e+00
  4.251e+06  0.993  0.9498  0.5377   0.8011   1.1219  3.12e+44     938  19.05  0.00e+00
  4.500e+06  0.967  0.9328  0.5355   0.7804   1.1252  3.16e+44     939  13.15  0.00e+00
  4.751e+06  0.924  0.9113  0.5311   0.7454   1.1315  3.14e+44     945  15.05  0.00e+00
  5.000e+06  0.865  0.8852  0.5230   0.6976   1.1421  3.06e+44     943  12.85  0.00e+00
  5.251e+06  0.791  0.8519  0.5093   0.6385   1.1376  2.93e+44     956  15.77  0.00e+00
  5.500e+06  0.707  0.8110  0.4916   0.5707   1.1253  2.79e+44     961  10.05  0.00e+00
  5.751e+06  0.615  0.7680  0.4777   0.4963   1.1121  2.65e+44     985  9.50  0.00e+00
  6.001e+06  0.519  0.7294  0.4713   0.4187   1.1017  2.52e+44    1005  10.60  0.00e+00
  6.250e+06  0.422  0.6938  0.4676   0.3405   1.0902  2.40e+44    1000  17.70  0.00e+00
  6.501e+06  0.328  0.6556  0.4611   0.2643   1.0771  2.29e+44    1007  22.10  0.00e+00
  6.751e+06  0.240  0.6090  0.4499   0.1937   1.0652  2.18e+44    1036  33.70  5.81e+11
  7.000e+06  0.162  0.5515  0.4345   0.1311   1.0540  2.09e+44    1071  41.94  1.49e+13
```

### fluxA.py (verbatim; all 15 dumps, t = 0 ... 7.0e6 s, every 5e5 s)
Run with a copy patched for NumPy 2.4 (only `float(z['time'][k])` -> `float(np.ravel(z['time'])[k])`: time is stored (15, 1)).
reduceA.py needed the grid dump agcar3d.x1grid.txt in the run dir (the code writes it to the job's cwd; copied in).
```
r/R_ph = 0.552 0.629 0.700 0.799 0.901 0.950 1.000 | L_MLT,IC/L0 there: 0.75 0.81 0.71 0.52 0.71 0.48 0.00
t 0.000e+00 w 1.000 | KE(<R) 0.00e+00 | rph 1.000/1.000 | vr bands   +0.0   +0.0   +0.0   +0.0   +0.0   +0.0
   rad   0.249  0.195  0.294  0.480  0.286  0.517  1.000
   conv  0.000  0.000  0.000  0.000  0.000  0.000  0.000
   kin   0.000  0.000  0.000  0.000  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   1.000  1.001  1.000  1.000  1.000  1.000  1.000
t 5.010e+05 w 1.000 | KE(<R) 7.34e+42 | rph 1.000/1.001 | vr bands   +0.0   +0.1   +0.1   +1.4   -0.1   -0.1
   rad   0.250  0.195  0.294  0.481  0.287  0.517  1.000
   conv  0.053  0.030  0.029 -0.001  0.000  0.000  0.000
   kin   0.000  0.000  0.000  0.000  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   1.053  1.030  1.029  1.000  1.001  1.000  1.000
t 1.001e+06 w 1.000 | KE(<R) 2.57e+43 | rph 1.001/1.001 | vr bands   +0.0   +0.4   +0.2   -4.3   -0.8   -0.1
   rad   0.249  0.197  0.293  0.480  0.289  0.518  1.001
   conv  0.025 -0.003  0.042  0.001  0.000  0.000  0.000
   kin  -0.000  0.000  0.000  0.000  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   1.024  0.999  1.041  1.002  1.003  1.001  1.001
t 1.500e+06 w 1.000 | KE(<R) 5.52e+43 | rph 1.002/1.006 | vr bands   +0.1   +0.9   +1.0   -2.5   +0.5   -0.1
   rad   0.248  0.195  0.290  0.469  0.293  0.518  1.001
   conv  0.025  0.014  0.008  0.015  0.000  0.000  0.000
   kin  -0.002  0.002  0.001  0.000  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   1.022  1.015  1.006  1.005  1.007  1.001  1.001
t 2.001e+06 w 1.000 | KE(<R) 9.59e+43 | rph 1.004/1.012 | vr bands   +0.1   +1.4   +2.8   +4.5   +5.7   -0.1
   rad   0.247  0.194  0.290  0.481  0.293  0.514  0.998
   conv  0.029  0.017  0.009  0.003  0.000  0.000  0.000
   kin  -0.010  0.004  0.004  0.001  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   1.018  1.020  1.010  1.005  1.007  0.997  0.998
t 2.501e+06 w 1.000 | KE(<R) 1.47e+44 | rph 1.013/1.025 | vr bands   +0.2   +1.8   +6.2   +7.5   +3.4   -0.1
   rad   0.247  0.194  0.288  0.479  0.302  0.510  0.992
   conv  0.030  0.016  0.010  0.005  0.001  0.000  0.000
   kin  -0.025  0.004  0.009  0.002  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   1.003  1.019  1.013  1.006  1.016  0.994  0.992
t 3.001e+06 w 1.000 | KE(<R) 2.07e+44 | rph 1.032/1.042 | vr bands   +0.3   +2.5   +8.3  +12.6   +9.9   -0.1
   rad   0.246  0.193  0.286  0.475  0.304  0.508  0.989
   conv  0.027  0.014  0.009  0.008  0.002  0.001  0.000
   kin  -0.046  0.000  0.011  0.004  0.000  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   0.979  1.013  1.012  1.007  1.020  0.992  0.989
t 3.500e+06 w 1.000 | KE(<R) 2.69e+44 | rph 1.057/1.072 | vr bands   +0.3   +2.5   +9.0  +12.8  +13.2   -0.1
   rad   0.246  0.193  0.285  0.472  0.304  0.501  0.978
   conv  0.023  0.011  0.007  0.008  0.003  0.001  0.000
   kin  -0.067 -0.009  0.009  0.006  0.001  0.000  0.000
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   0.953  1.001  1.007  1.006  1.021  0.985  0.978
t 4.001e+06 w 1.000 | KE(<R) 3.05e+44 | rph 1.080/1.092 | vr bands   +0.1   +2.0   +9.7  +12.1  +12.0   -0.1
   rad   0.246  0.193  0.285  0.469  0.303  0.499  0.959
   conv  0.018  0.009  0.004  0.006  0.004  0.001  0.000
   kin  -0.084 -0.022  0.002  0.006  0.001  0.000  0.001
   mlt   0.751  0.805  0.706  0.520  0.714  0.483  0.000
   sum   0.931  0.985  0.997  1.001  1.021  0.983  0.960
t 4.500e+06 w 0.967 | KE(<R) 3.16e+44 | rph 1.098/1.114 | vr bands   +0.0   +1.2   +9.9  +12.2  +11.9   -0.1
   rad   0.246  0.194  0.284  0.467  0.302  0.499  0.929
   conv  0.015  0.007  0.000  0.002  0.005  0.001  0.000
   kin  -0.092 -0.038 -0.010  0.003  0.001  0.000  0.001
   mlt   0.726  0.779  0.683  0.503  0.690  0.467  0.000
   sum   0.896  0.942  0.957  0.975  0.999  0.968  0.930
t 5.000e+06 w 0.864 | KE(<R) 3.06e+44 | rph 1.114/1.127 | vr bands   -0.1   +0.6   +8.6  +10.0  +10.8   -0.4
   rad   0.247  0.195  0.288  0.465  0.306  0.491  0.877
   conv  0.014  0.007  0.000 -0.003  0.005  0.002 -0.000
   kin  -0.090 -0.047 -0.023 -0.004  0.001  0.000  0.001
   mlt   0.649  0.696  0.610  0.450  0.617  0.417  0.000
   sum   0.820  0.851  0.875  0.909  0.928  0.910  0.877
t 5.500e+06 w 0.707 | KE(<R) 2.79e+44 | rph 1.131/1.140 | vr bands   -0.2   -0.9   +4.9   +5.2   +8.6   -0.6
   rad   0.247  0.195  0.289  0.471  0.311  0.463  0.798
   conv  0.015  0.007 -0.002 -0.003  0.003  0.001 -0.000
   kin  -0.082 -0.052 -0.032 -0.009  0.000  0.000  0.000
   mlt   0.531  0.569  0.499  0.368  0.505  0.341  0.000
   sum   0.711  0.720  0.754  0.826  0.818  0.806  0.798
t 6.001e+06 w 0.518 | KE(<R) 2.52e+44 | rph 1.139/1.145 | vr bands   -0.1   -0.8   -1.4   -4.2   -4.3   -1.0
   rad   0.249  0.196  0.291  0.471  0.335  0.448  0.717
   conv  0.014  0.009  0.003 -0.006  0.000 -0.000 -0.000
   kin  -0.069 -0.052 -0.038 -0.012  0.000  0.000 -0.000
   mlt   0.389  0.418  0.366  0.270  0.370  0.250  0.000
   sum   0.584  0.571  0.622  0.722  0.705  0.698  0.717
t 6.501e+06 w 0.327 | KE(<R) 2.29e+44 | rph 1.128/1.143 | vr bands   +0.0   -0.6   -9.9  -14.2  -14.8   -1.7
   rad   0.250  0.197  0.292  0.473  0.355  0.432  0.645
   conv  0.017  0.013  0.005 -0.006 -0.003 -0.001 -0.000
   kin  -0.053 -0.048 -0.038 -0.011  0.000 -0.000 -0.001
   mlt   0.246  0.264  0.231  0.170  0.234  0.158  0.000
   sum   0.460  0.426  0.490  0.626  0.585  0.589  0.644
t 7.000e+06 w 0.162 | KE(<R) 2.09e+44 | rph 1.087/1.111 | vr bands   +0.2   -0.4  -23.4  -29.6  -33.7   -2.6
   rad   0.250  0.197  0.290  0.471  0.379  0.416  0.546
   conv  0.019  0.016  0.009 -0.009 -0.005 -0.001  0.000
   kin  -0.039 -0.043 -0.034 -0.008 -0.000 -0.000 -0.005
   mlt   0.122  0.131  0.115  0.084  0.116  0.078  0.000
   sum   0.352  0.300  0.380  0.538  0.490  0.494  0.540
```

### Files (pushed here, docs/handover/rampA-1011/deltaai-run/, 9.5 MB)
agcar3d.user.hst, agcar3d.hydro.hst, run.3363376.log, prof.npz, slice_*.npz (dumps 0, 2, 6, 8, 10, 12, 14).
On DeltaAI: /work/nvme/bivj/jma20/rampA_1011/run (bins every 5e5 s to 7.0e6, rst 00000-00014), out/.
