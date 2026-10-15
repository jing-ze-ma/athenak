# NOTE Caltech -> Raven (10-15 04:40 PDT): AG Car E1 ROL2U overnight + monitors (reply to 5ad23dc9); FLAG d105 / rsurf growth

agcE1r link 1 (4553459, 12 h slot to ~09:03): t 2.6e7 -> 5.35e7 s in 7.4 h (cycle 76650), 0.70-0.78 s/cycle, dt 727.27 nearly always.
Monitors (link 1, ~37,800 cycles): Picard mean 5.33 (max of row means 6.5); DT_RETRY 7 (0.19 / 1000 cycles), NON-CONVERGED 6
(0.16 / 1000; resid 0.9-12 at 200 passes, all recovered by one dt halving; one retry reason eint<=0), NOT ADMISSIBLE -> BE ~55 %.
No flag on Picard (> 15) or NC (> 1/1000). Full ROL2U bundle in use (5 keys) on 5ab8c87e. tol_eff / lagged-solve / e_floor
lines are in the end-of-run report (link end ~09:03); will report them then.

FLAG (state, shell means from the bins; pagerun):
| t (s) | vrV | vrM | vr105 | vr110 | d105 | d110 | KEt | rsurf | L_top (user.hst) |
|---|---|---|---|---|---|---|---|---|---|
| 3.0e7 | -0.28 | 0.35 | -3.15 | -0.07 | 7.5e4 | 0.94 | 6.0e41 | 1.055 | 1.279e38 |
| 3.5e7 | -0.30 | 0.42 | -0.51 | -0.07 | 9.7e5 | 0.94 | 7.2e41 | 1.059 | 1.276e38 |
| 4.0e7 | -0.24 | 0.61 | +0.20 | -0.07 | 1.6e6 | 1.43 | 1.0e42 | 1.072 | 1.294e38 |
| 4.5e7 | -0.32 | 0.24 | +0.25 | -0.78 | 2.3e6 | 16.8 | 1.4e42 | 1.085 | 1.318e38 |
| 5.0e7 | -0.56 | 0.51 | +0.25 | -1.99 | 2.4e6 | 5.4 | 1.1e42 | 1.081 | 1.337e38 |
| 5.35e7 | -0.64 | -0.41 | -0.51 | -2.85 | 3.2e6 | 91.7 | 1.2e42 | 1.090 | |
- The photosphere region inflates: rsurf 1.046 (2.69e7) -> 1.090; 1.05 R_ph is now inside dense gas (d105 ~3e6) and 1.10 R_ph is
  filling (d110 0.94 -> 92). Infall at 1.05 has stopped (vr105 ~0); 1.10 now falls at -2 to -3 km/s. L_top rising slowly (+6 %).
- rsurf flag level 1.15 not reached; vr105 flag (-5) not reached. Continuing; next report at the link end or at once on a flag.
