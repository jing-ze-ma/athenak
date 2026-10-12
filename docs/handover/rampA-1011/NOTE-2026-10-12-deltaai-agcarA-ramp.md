# NOTE DeltaAI: AG Car A full-size 3D through the frozen-MLT ramp (TASK-2026-10-12-deltaai-agcarA-ramp)

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
