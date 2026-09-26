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

**Update 09-25 (dhj, /viper/ptmp2/jinma/scsp_0925):** p-psi SC prototype on prod4 rot 80 (full sphere 128x64x128,
grey kappa_R, S=B, no beam). vet_col 99 % |df_rr/f| at tau 0.1-1: day 0.4 %, terminator 1.6 %, NIGHT 8 % (max 10 %);
flux at tau 2/3 within 1.5 % median. p-psi SC needs a tangent ray at EVERY shell (sparse stride 8 -> 30-40 % errors);
then <0.1 %. Lateral anisotropy is smeared spatially by any SC (misses LC D_r,lat 0.13 by ~0.06). Cost ~60x vet_col
(1176 directions/cell). sp operator uses only f_rr (impl_offdiag forced none; rad_m1_sph.cpp:22-24). VERDICT: NO-GO on
C++ SC for sp; revisit only if the off-diagonal operator becomes stable or a full LC run (lc.py, needs a CPU node;
sbatch -p general refused for the account) shows night-side f_rr errors that matter.

**DO study 09-25 (/viper/ptmp2/jinma/do_study_0925, Jiang 2014/2014b/2021 fetched there):** Jiang uses GLOBAL Cartesian
angles on every grid (n.grad I = div(nI), no angular-flux terms; curvature only via face areas/volumes) and block-Jacobi
implicit iterations (no sweep -> no cycle on poles/phi/cs; slow in thin gas, ~100-1000 iterations). AthenaK src/radiation
is GR-only, explicit, no face areas -> not reusable as an sp/cs solver (geodesic angle set + moment sums are).
Ranking: (1) keep vet_col, extend it to cs (sp-only today, rad_m1_vetcol.cpp:220; 4.2 ns/cell vs 59 ns M1 step);
(2) only if lateral tensor / night 8 % matters: DO formal solver with global geodesic angles (92-162), exact radial
sweeps (acyclic per global direction) + lateral Jacobi, ~3-25x vet_col, ~2k lines, prototype in Python first;
(3) angular fluxes for option A: no; (4) replace M1 by DO: no.

**Lagged geometric-source SC ("lagsrc", user's Jiang-2021 idea) 09-25: NO** (/viper/ptmp2/jinma/scsp_0925/README_lagsrc.md).
Pure lagging does not converge (period-2 oscillation, 3-13 % of I clipped, thin top 30-77 % wrong); with the own-cell
term implicit + ascending-mu ordering it converges (first-order upwind = direct solve in 1 sweep). Night photosphere
99 % f_rr error stays 1.0-1.2 % for 128-512 directions (vs vet_col 8.1 %, naive A 0.23-0.35 %, p-psi 0.08 %); the error is
spatial (spreading S_geo over dtau~1 segments with linear SC weights), not angular. Conservative form helps the thin
top (1.5-2 %) but needs >= 8 sweeps (more reads than p-psi). Only p-psi with a tangent ray per shell reaches 0.1 %.
Untried: second-order SC segments.
