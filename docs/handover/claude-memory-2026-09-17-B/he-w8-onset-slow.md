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

box_w9 onset (analysis_0916/onset_he_w9/RESULTS.txt, dumps 0-3, 09-17 04:40): the seed
(9.95e8 in k4-6 at i=84) DECAYS 21x over two turnovers (-2.1, -1.0/turnover), then turns up
d2->d3 at +0.38/turnover; 6.9e7 at d3 = 30x below w4 in amplitude (w8 was 200x), ~7 turnovers
to the w4 level at that rate. Deep layers filled (i=70 comparable to i=84). Verdict: keep w9
running, check again at dumps 5-6 before raising vpert further.

UPDATE 09-17 04:45, w8 through dump 7 (7 turnovers): growth is ACCELERATING, not slow. rms
v_z/v_MLT at i=84: 0.16 (d1) 0.19 (d3) 0.20 (d4) 0.23 (d5) 0.30 (d6) 0.45 (d7); total power
growth 0.12 -> 0.80/turnover, k4-6 0.42 -> 0.70; k2-3 now the dominant band (inverse cascade),
k7-12 seed decaying; deep planes i=70/49 growing 0.8-1.2. The 03:00 "25 turnovers" verdict was
a linear extrapolation of the early rate and was wrong; w9 was launched on that basis and is
now redundant except as a seed-amplitude comparison. Saturation expected ~turnover 9-10.

09-17 09:30: box_w8 at 13 turnovers: rms v_z/v_MLT 0.93 (tau 20), 0.90 (tau 50); closed network of
4-5 cells; dT/T at tau 1 = 5e-6; KE ln-ratios turns 8-13: +0.64 +0.25 +0.02 +0.09 +0.24 +0.18.
box_w9 (vpert 1e-2, k 3..4) CANCELLED AND DELETED by user decision (jobs 11744300-320, dir removed).
box_w8 hung 08:09 on a zero-byte cbin write during the quota event; cancelled 08:45, chain restarted
from rst 00026. B star prod_w7 at 23 turnovers statistically steady (0.73/0.88 v_MLT at tau 300/1000).
Page v7 version 3 has both.
