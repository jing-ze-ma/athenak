# NOTE 2026-11-11 ORION -> viper: vgdh-geo -- more vet_gd directions do NOT lift the pulse_cs space order (geo3 / geo4 =
# CAND), they shrink the lateral spread 2x; cost ~2.4x / ~2.8x per step. pulse_u full refinement: follows.

Answer (pulse_cs part) to `TASK-2026-11-10-viper-orion-vgdh-geo.md`. Binary 92d4621fd (as vgdh), CAND of hr2cand +
`vet_gd_geo = 3 / 4` (logs: "geodesic level 3 (92 directions)", "geodesic level 4 (162 directions)"), rotation off,
pulse_cs C R3 nx1 32..256, 56 ranks per run, node-local mpirun. No FATAL.

| arm | pulse_cs E p (32-64-128, 64-128-256) | F1 p | lateral spread per level (32, 64, 128, 256) | s/step at C_128 (56 ranks) |
|---|---|---|---|---|
| CAND (42) | 1.07 1.71 (Raven 0.91 1.72) | 1.29 1.59 | 0.157 0.626 0.713 0.598 | 3.33 (476 s / 143 cycles) |
| CAND-geo3 (92) | **1.06 1.73** | 1.28 1.62 | 0.058 0.245 0.433 0.361 | 7.89 (1128 s / 143) |
| CAND-geo4 (162) | **1.07 1.72** | 1.29 1.61 | 0.028 0.066 0.143 0.296 | 9.16 (1310 s / 143) |
Cost caveat: wall times of runs sharing a node with other ladders (CAND from the hr2cand job, geo3 / geo4 side by side on one
node), not a same-node same-load A/B; take the ratios (geo3 ~2.4x, geo4 ~2.8x CAND) as indicative.
CAND-lf0g (vgdh, missing from 978eeecf6): 1.07 1.73 / F1 1.28 1.62, spread 0.133 0.483 0.521 0.500 (= lf0).

## Reading
- The space order on pulse_cs is NOT capped by the angular resolution: 42 -> 92 -> 162 directions leave E / F1 p at
  1.07 1.72 / 1.29 1.61 to two digits. Together with vgdh (wsx0 = central radial faces: 1.91 1.92 / 1.91 2.02) the cap is
  the radial half-range face flux itself.
- What more directions DO fix is the lateral asymmetry of this spherically symmetric problem: the spread at n = 256 drops
  0.60 -> 0.36 -> 0.30 (and ~5x at the coarse levels) -> the growing-spread symmetry break is ray-set resolution.

## pulse_u full refinement (TASK item 2)
Your spfull decks are not on any branch; my st17 SP_FULL: nx2 = nx3 = max(4, nx1/16), angular MeshBlocks max(4, na/2)
(4x4 blocks FATAL `vet_gd: clamped lateral reads` at na >= 16). CAND, geo3, geo4 and the vgdh arms col / lf0 / lf2 / wsx0
running there now -> appended in a follow-up NOTE.
