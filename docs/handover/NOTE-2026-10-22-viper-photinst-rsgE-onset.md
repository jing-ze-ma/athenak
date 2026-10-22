# NOTE 2026-10-22 viper -> caltech: RSG E N=32 reproduces the cliff-cell runaway; same signature as HG (re b85d1750)

Run: RSG E deck at cs N=32 (fresh, f55be371, viper 1 node x 2 MI300A), dir /viper/ptmp2/jinma/photinst_1022/runs/base, scr/ana.py.
Metric: lateral sd/mean of e_int/rho on the cliff shell r = 0.9950 R_ph (= cell 319, the one-cell MLT deposit, 68 % of L).

| t (d) | 21 | 24 | 27 | 30 | 33 | 36 | 39 | 42 | 45 | 48 |
|---|---|---|---|---|---|---|---|---|---|---|
| e/rho sd | 0.031 | 0.036 | 0.034 | 0.041 | 0.049 | 0.061 | 0.077 | 0.098 | 0.125 | 0.169 |

- Onset ~27-30 d at N=32; e-fold 11-12 d (36->48 d: 11.8 d), consistent with your RSG E prod growth (8.8 d at N=128).
- At 45-48 d on that shell: rho sd 2.7-3.5 % (was 0.2-0.4 %); corr(rho, T) = -0.99 (isobaric two-phase, as your HG);
  corr(v_r, T) = +0.97 (hot cells move out = expansion, no overturning signature); corr(T_cliff, F_r 3 cells above) falls from
  +0.4..+0.6 to ~0 by 45 d. Trad sd on the shell 1.5-2.4 % (the e/rho sd is dominated by the H ionisation energy).
- Surface still small at 48 d: T_eff sd 0.06-0.08 %, panel contrast |fc/ed - 1| <= 1.6e-4, R_ph sd 8e-5.
- dt = 1000 s throughout, 0 dt retries through 40 d: the retries you saw in RSG E prod come AFTER the onset, not before.
=> RSG E has the same column-local runaway in the converged discrete system as HG (one-cell deposit at the under-resolved H cliff).

Next on viper (restart-from-growing-state, per your suggestion): restart at t = 40.0 d (e/rho sd 0.08, growing); arms
mlt_fsub_smooth = 4 and 8 (your 2f50377d), time2_vet_col = rebuild; base continuation to 60 d as reference. We will add your
Kirchhoff / T-integral face-diffusivity fix on the same restart as soon as you push it. Our fresh-start arms (lag passes,
gas Newton off, kP=kR, MLT off, table x4) are cancelled after your restart-arm correction.
