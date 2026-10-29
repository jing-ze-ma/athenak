# NOTE caltech -> viper (2026-10-29): defaults-1029 on fin51-r-fix (USER 10-29)

fin51-r-fix gains 5 commits on top of 2ff6a8ee:

| commit | key | default now |
|---|---|---|
| vtxfix-1029 | `<rad_m1>/cs_vertex_mono` (monotone two-point lateral M1 at the 8 cube vertices) | **true on the cubed sphere** (no-op elsewhere) |
| vtxfix-1029 | `<rad_m1>/implicit_pos_floor_cons` (the write-back E floor takes its raise from the same column; no energy created) | **true** |
| autobe-1029 | `<rad_m1>/time2_auto_be = N` (after a hesdirk2 stage fallback skip the doomed stage attempt on the next 1, 2, 4, ... N steps, BE redo at once; a successful retry resets) | **16** |
| stiffk1-1029 | `<rad_m1>/time2_k1_stiff = s` | 0 = off (opt-in; helps HG, hurts RSG E) |

Naming a key with the old value (false / 0) restores the old behaviour bitwise.
**Runs that restart on this tip without naming the keys change behaviour** (not bitwise).

Why (Caltech RSG E / HG, 10-29):
- RSG E 1681-1704 d plume: M1 E-floor cells at the (-x,+y,+z) cube vertex; the floor created
  5.4e42 erg, heated the vertex cells to ~45 kK, lateral blast -> evacuated column -> plume.
  With both vtxfix keys: zero floor energy, no event (re-run from 1675 d).
- hesdirk2 stage 1 is not positivity preserving for k dt >> 1; HG fell back on ~81 % of steps,
  RSG E on 42-100 %. time2_auto_be: HG 0.40 -> ~1.8-2 d/h (with k1_stiff 1.0), RSG E 12 -> 27-38 d/h;
  HG accuracy vs the prod dump 391: front T rms 2.4e-4, L rms 1.4e-3, photosphere L 0.9 %.
- The auto_be skip state is not in restart files (a restart costs one stage attempt).

Caltech production switched: RSG E prodH and HG (keys named explicitly), S1 full sphere and
accretor run12 (defaults) at their next restart.
