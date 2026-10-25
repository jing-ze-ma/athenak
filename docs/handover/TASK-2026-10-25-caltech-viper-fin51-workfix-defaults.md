# TASK 10-25 Caltech -> viper (USER): fin51-r-fix <- work fix, then defaults-1025

Two steps, both on the fork (jing-ze-ma/athenak). Please gate, merge and report back here.

## 1. Fast-forward fin51-r-fix 958cac78 -> c60714bd (branch rsgfresh-1025)

c60714bd = fin51-r-fix 958cac78 + two commits (opt-in, key off = bitwise):

- ba4bb54d kick581-1024 `<rad_m1>/implicit_work_fix`: the implicit write-back subtracted
  the full gas work after the solve where the E row already carried it (E - work << 0,
  floored). This was the RSG E 580-581 d energy kick.
- c60714bd the same fix after the implicit_vimp positivity fallback: dv^k is still built
  (ImplicitVimpBuild, Jacobian unused), so the write-back subtracts only work - w^k.

Caltech gates (H200):
- With the key off, the results are bitwise identical to 958cac78.
- RSG E kick replay from 580 d, 12 ranks: with the fix, clean through 584.5 d (Ltop/Lin
  max 0.762, NON-CONV 0) and bitwise identical to the W12 reference.
- vimp-fallback gates, 4/4 pass: R4 Ltop/Lin 0.759 vs 2.57 without the fix (W4), and
  W12b 0.760. The created energy is ~5-7e39 erg (solve undershoot).
- Production: RSG E (fresh start 10-25, ~150 d) and HG global (~47 d) both run this
  binary (md5 140d0e79) with `rad_m1/implicit_work_fix=true`. Both are clean.

Recommendation: Raven run12 (on 958cac78) adds `rad_m1/implicit_work_fix=true` at its
next link.

Note: rt-integration was stale (d28dacd8, 10-10). Caltech fast-forwarded it to c60714bd
today. fin51-r-fix stays the production branch.

## 2. Then gate and merge defaults-1025 (83d97f8c, on c60714bd)

USER 10-25: the RSG E / HG production keys become code defaults. A named key still wins.

| key | old default | new default |
|---|---|---|
| `rad_m1/implicit_work_fix` | false | true |
| `rad_m1/time2_posfix` | 0 | 4 for hydro; stays 0 for MHD, which is not wired (no FATAL) |
| `rad_m1/time2_vet_col` (sph, fresh start) | predict | extrap2 (a restart that lacks the key keeps lag) |

Not yet built or gated anywhere. Suggested gates:
1. New binary, no keys == c60714bd with the production XKEYS (bitwise; sph RSG/HG-like deck).
2. New binary with the old values named (posfix=0, work_fix=false, vet_col=predict)
   == c60714bd with no keys (bitwise).
3. An MHD deck runs (posfix stays 0) and a Cartesian deck is unchanged.
4. A restart from an old rst that lacks the keys behaves as intended.

If they pass, fast-forward fin51-r-fix to defaults-1025.
