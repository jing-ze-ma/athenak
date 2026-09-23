---
name: cs-mhd-prod2-run
description: cs_mhd_prod2 (bench/cs_mhd_prod2), the cubed-sphere MHD production run FROM SCRATCH on ebd57244 with the 2026-09-07 defaults (WB polytropic + rot_potential + radiative diffusion tau blend with ck-table kappa_R); launched 09-07 05:20, lead 11526708 (6.5 h, cut by the 12:00 maintenance) then chain 11526681 -> 11526682
metadata:
  type: project
---

**What:** bench/cs_mhd_prod2. Input = cs_mhd_prod's (128x32x32 cs, stretch, bbot 3 G,
point-mass gravity, floors) + `<mhd>` wellbalance_dynamic/wb_x1/wb_option=polytropic/
wb_cache_every=10 + isotropic_conduction=radiative, rad_tau_lo/hi 30/300,
rad_kappa_src=table, rad_flux_inner=-1 (pgen sets 4.735e6 = sigma 537.6^4), rad_met 0
+ problem/rot_potential=true. rst every 2 rot (6.10e5 s) because the lead job only has
6.5 h before the 09-07 12:00 -> 09-12 12:00 maintenance. Binary build_dhj_gpu ebd57244
(05:08). Start-up prints the Rosseland table vs Freedman and the wall flux.

**Jobs:** 11526708 lead (06:30 wall, running 05:19 on vipa1287); 11526681 and 11526682
(23:55 each) chained afterany, pick up the newest rst. 11526680 was a false start
without an input file (the input writer asserted on a comment mentioning rot_potential);
it failed in seconds and the chain was re-pointed with scontrol update Dependency.

**First minutes:** 43 ms/cycle on 2 GPUs (cs_mhd_prod was ~40), dt 20-29 s, no floors
reported. Compare against cs_mhd_prod (11516277, same physics minus the three new
pieces, healthy past rot 19) -- that is the A/B for WB + rot_potential + deep RT:
deep v_r/v_theta drift, deep T profile evolution, bottom-boundary behaviour, jet.

**Trap:** `sbatch` succeeded on a directory whose input file did not exist; check the
input exists before submitting a chain, and check `squeue` for the false start.
