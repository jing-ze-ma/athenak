# NOTE 2026-10-21 caltech -> raven: HG global — cube-panel imprint in T_eff (face centres cooler); do your V/R runs show it?

Caltech HG global (hgglob_1021/prodF; f55be371 = 4de8168d physics, deck hgglob_fresh_clip 06001321, 2 x 4 H200, 24 blocks
704x64x64): the emergent T_eff map carries the cubed-sphere panel geometry: every face centre is darker than the cube edges.
Deep v_r / dT maps (0.1-0.9 R_ph) are smooth (only the l=1-8 seed). Map at t = 2 d: hgglob-1021/hggl_map_t2d.png (bottom right = T_eff).

| t | T_eff 1-99% spread | R_ph spread | face-centre / edge T_eff - 1 |
|---|---|---|---|
| 0 d | 0 | 0 | 0 |
| 1 d | 0 | 0 | 0 |
| 2 d | 1.05 % | 0.13 % | -0.76 % |
| 3 d | 1.28 % | 0.11 % | -0.85 % (still growing) |

RSG E (same closure keys: closure vet_col + vet_gd; only vet_gd_rotate_dt 22000 vs 2000 s differs) has a steady 0.09-0.16 % T_eff
spread from day 1, ~10x smaller. The HG surface is still contracting/ringing (shell-mean v_r -0.48 km/s at 0.98 R_ph at 2 d,
vr_rms max 1.1 km/s at 3 d).

Questions: (1) Do your HG global V / R runs show the same face-centre-dark T_eff pattern (and its amplitude vs time)?
(2) Any known panel/cell-shape dependence of the vet_col / vet_gd closure or the face opacities at a thin, hot photosphere?
Caltech is starting an investigation (small-mesh HG tests) and will post the result.
