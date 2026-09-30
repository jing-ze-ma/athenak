---
name: hebox-cfl-saturated-verdict-0929
description: He box M1 saturated cfl verdict 09-29 (R3 vs H9 vs P9, t 60000-74000): marginal, cfl 0.9 ~5 % low in v1'/v_MLT and KE1
metadata:
  type: project
---
/viper/ptmp2/jinma/hebox_cfl2_0927, window t 60000-74000, FeCZ tau 4.4-64 (RESULTS.md last section).
v1'/v_MLT: R3 1.257 (snapshot std 0.03), H9 1.192 (0.115), P9 1.208 (0.084); |H9-R3| 0.065 vs noise |H9-P9| 0.016.
KE1: H9 -4.7 %, P9 -3.6 % vs R3 (noise 1.2 %); KEh and Fc/F within noise; mean T/rho profiles differ ~1e-5.
Verdict MARGINAL: fails the strict round-off-spread rule ([[no-accuracy-sacrifice]]) but only one noise pair, and the
difference is below the snapshot scatter. Both cfl-0.9 members sit below R3, so it looks systematic.
Pending: H6 (cfl 0.6), which was due ~05:50 09-29. The user has not decided yet; this bears on [[m1-hesdirk2-cfl-recommendation]].
H6 IN (06:35): cfl 0.6 = R3 within noise (v1'/v_MLT 1.265 vs 1.257, KE1 +1.6 %); deficit only at 0.9. Recommendation to the
user: cfl 0.6 for the science window (1.9x cheaper than 0.3), 0.9 for relaxation. User decision pending.
**USER DECISION 09-29 ~06:40: M1 runs use cfl 0.6 in the science window, cfl 0.9 in relaxation** (hesdirk2). Supersedes the
09-26 "cfl 0.9" decision for the science phase. For the He presn wedge re-check with its own gate 1b dt-order data.
