---
name: cs-vs-sp-dhj-rot50
description: 2026-09-13 cs_mhd_prod3 vs sp_mhd_nopole (polar swap off) at rot 20/50 (bench/sp_mhd_nopole/analysis/cs_vs_sp): FLOW agrees (total KE within 15%, zonal jet identical: eq +540/+515 m/s, peak ~+800 at 60 deg) but sp loses DEEP field 1.6x faster (i<20, mid-latitude bulk, 79% of the loss at |lat| 20-53; polar rows 12x decay but only 5% of the loss) and has 1.6x the radial KE (deep, equatorial + high-lat, growing on both). sp fires 26x the ENERGY floors: dayside photosphere i 58-110 (the substellar dissociation front), NOT the night-side top (cs pattern). Likely setup, not grid: cs has rad_implicit_x1 + rad_cap_ang + rt_use_cons, sp does not; sp has 1.55x more steps/rot. cs seams/vertices innocent. Dump ME is 17-26% below the face ME at rot 50 (grid-scale field)
metadata:
  type: project
---
Discriminator to run: an sp arm with the cs RT switches (rad_implicit_x1, rad_cap_ang 0.5, rt_use_cons) at fixed
everything else. Remaining sp polar-row artefact (j 0,1,62,63: field decays 12x, KE_r density 2x, 3.6% of volume):
candidates polar_emf_diss and the polar-row reconstruction, not the Riemann solver (mask 0 in this run).
