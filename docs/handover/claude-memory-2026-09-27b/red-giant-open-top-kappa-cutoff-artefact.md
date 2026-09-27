---
name: red-giant-open-top-kappa-cutoff-artefact
description: "THE 3 L EMERGENT FLUX after the cap burst (I2/I3, 09-11 01:00) is an OPACITY-CUTOFF ARTEFACT: the surface layer inflates (mass above i=340 x4.5, E x10 between 5.6e5 and 7.0e5), the tau=2/3 surface and the tau=100 cut are pushed out to rad_kappa_rmax=3.6528e12 where kappa=0 by fiat, and the RT radiates 9000-10000 K ionized gas as a bare blackbody (L_cut columns have icut 392-396 at 6600-7000 K). Not shock heating (only one Mach-1.5 column at 6.2e5), not the vertex chimney (T-ranking on a flat 9400 K plateau). I4 (no cutoff, active medium) is the direct test"
metadata:
  type: project
---
Evidence (I2): fraction of columns with icut >= 390: 0 -> 0.1% -> 9.8% -> 13.5% -> 4.1% at
5.6/6.2/7/8/9e5; per-column sigma T_ph^4 sum 0.66 -> 3.16 -> 2.41 L tracks L_out 0.62 -> 5.1
-> 3.7 to ~30%; T_ph median 3050 -> 3890 K; i(tau=2/3) p95 356 -> 396. Physical diffusion
flux at tau=100 reconstructed from the dumps stays 0.05-0.12 L while the run's L_cut grows to
2.1 (int_at_cut=false, I_up = B(T_cut) at two_stream_rt.hpp:1287). Lidded B12/B13 at 3.4e7:
top layer a UNIFORM 9200 K at 1e-8 with L_out 1.3-1.5 L -- the same state I2 reaches by 8-9e5,
held by the lid. So the model's surface wants to be ~9000 K/1e-8 (over-luminous 1.3-1.5 L),
and the thin 2200 K IC atmosphere is not its equilibrium.
The agent's "standing +10-50 L dE_tot/dt" (dumps, incl. PE; lidded too; no step at the burst)
is AT THE SINGLE-PRECISION NOISE FLOOR of the .bin dumps (E_tot 1e50 in float32 -> 1e43
resolution -> 50 L over 1e5 s) and the hst has 6 digits: NOT established. A real budget
needs an in-code 16-digit reduction. Do not chase it from dumps.
Fix candidates: rad_kappa_rmax = x1max (or rad_kappa_above = table) so the inflated surface
keeps its opacity -- exactly I4's configuration; with the 2e4 K inert corona the join must
move up or go. Scripts: scratchpad q12,q34,q45,q6-q9.py (session f9284861) run from
_analysis_0910/. See [[red-giant-clamp-off-open-top-ok]], [[red-giant-implicit-radial-diffusion]].
