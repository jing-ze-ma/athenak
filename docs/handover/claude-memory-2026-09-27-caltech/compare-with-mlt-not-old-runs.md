---
name: compare-with-mlt-not-old-runs
description: User rule 09-17 - judge convection onset/amplitude against v_MLT (and F_conv/F), never against an earlier run like box_w4 (which was 19x v_MLT, over-driven by the old kick)
metadata:
  type: feedback
---

Judge a convection run against the PHYSICAL target: v_MLT at the driving layer (He star
1.86e4 cm/s at z=-1.02e8, i~80; B star 2.44e5 cm/s) and the convective flux F_conv/F. Not
against an older run: box_w4 at 1 turnover was 3.6e5 cm/s = 19 v_MLT (old whole-column kick,
f-mode contaminated), so "200x below w4" meant nothing.

**Why:** user 09-17: "why are you comparing with w4? you should compare with mlt".

**How to apply:** in onset checks print rms v_z / v_MLT per dump and band, and F_conv/F.
He w9 at 3 turnovers: 0.62 v_MLT at i=84 (w8 0.19). B star dump 10: 0.4-0.6 v_MLT at
tau 300-1000. See [[he-w8-onset-slow]], [[bstar-w7-linear-saturation]].
