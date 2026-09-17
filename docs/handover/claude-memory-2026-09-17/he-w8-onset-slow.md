---
name: he-w8-onset-slow
description: He box_w8 (velocity seed 1e-3 c_s, kmin/kmax 4..6 per component) grows but slowly; seed landed in k7-12 because kmin/kmax apply per component; box_w9 relaunched 09-17 with vpert 1e-2, k 3..4
metadata:
  type: project
---

He box_w8 onset (analysis_0916/onset_he_w8/RESULTS.txt, dumps 1-3 = 3 turnovers, 2026-09-17):
the seed is NOT erased (unlike the entropy seeds w6/w7). k4-6 v_z power at i=84 grows
1.02e6 -> 1.55e6 -> 2.37e6 (cm/s)^2, +0.42/turnover in power (amplitude e-fold 4.7
turnovers); deep planes i=70/49 grow +0.5-0.6/turnover, like box_w4. But amplitude is ~200x
below box_w4 at the same turnover (w4's old whole-column kick started at 3.6e5 cm/s), and
the dominant power sits in k7-12 and is flat: `vpert_kmin/kmax` are drawn PER COMPONENT
(box_convection.cpp ~L2500), so 4..6 gives |k| = 5.7-8.5, not 4-6.

**Why:** at +0.42/turnover the 40-turnover production would spend ~25 turnovers in linear growth.

**How to apply:** box_w9 (launched 09-17 ~02:30, same everything) uses vpert 1e-2 c_s and
kmin/kmax 3..4 (|k| 4.2-5.7). w8 left running for comparison; user decides which to keep.
For |k| in a band [a,b], set per-component kmin ~ a/sqrt2, kmax ~ b/sqrt2. See [[seed-box-modes]].
