---
name: vet-col-enough-no-sc-on-sp
description: 09-24 measurement: vet_col's column-only Eddington tensor vs a 3-D long-characteristics solve on the most porous He4 sp wedge snapshot -- <1 % below tau 1, <=5-6 % at the corrugated photosphere/thin top; SC on sp not justified
metadata:
  type: project
---

/viper/ptmp2/jinma/vetcmp_0924/README.md (wedge9d dump 00020, 3.5 turnovers, rms rho/<rho> 0.8-1.8 at the photosphere;
192 rays; S = B, kappa_R incl. electron scattering; no S-from-E variant, scattering in the source untested).
|df_rr/f| 99 %: tau>10 <=0.15 %, tau 1-10 0.56 %, tau 0.1-1 5.1 % (max 5.8 %), thin top 3.5-4.4 %; no cell > 10 %.
Lateral anisotropy / off-diagonals (zero in vet_col): 99 % <= 0.4 % below tau 1, 1-3.5 % above.
Flux at each column's tau=2/3: ratio median 0.997 (min 0.917), direction off by median 0.95 deg; slice sums 5-6 % low at the photosphere.
Error comes from the large-scale corrugation of the tau=1 surface, not from cell-scale porosity.

**Why:** decides whether to build the full short-characteristics VET on sp (option A, 10-15 agent-days).

**How to apply:** default answer = vet_col is enough for the He sp work; revisit SC on sp only if a problem needs <5 %
photospheric flux accuracy or strong lateral anisotropy (irradiated dhj-like cases). Related: [[m1-hesdirk2-cfl-recommendation]].
