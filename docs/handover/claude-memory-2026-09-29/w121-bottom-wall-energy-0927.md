---
name: w121-bottom-wall-energy-0927
description: flux_hst on w1x rot 60-63 shows the inner wall injecting ~1.5e30 erg/s fluid energy flux (160x L_int) with zero net mass flux
metadata:
  type: project
---
dhj-fluxhst merged 09-27 (392f1e50; columns in dhj.user.hst). GPU check from w1x rst 00120 (rot 60->63,
/viper/ptmp2/jinma/fluxhst_0927/gates/gpu_new): Lrad_bot 9.39e27 = sigma T_int^4 4pi r_in^2 (OK);
Lir_top/(Lsw_abs+Lrad_bot) = 0.835 (atmosphere still gaining); Etot_bot 1.51e30 erg/s, almost all fluid flux,
Mdot_bot oscillates +-1e18 g/s with mean ~0; hydro dE/dt ~2e30 matches. Likely source of the deep drift
(100 bar +13-16 K / 10 rot, [[next-prod-ck-c2]]).
**Why:** the 300-rot w1x/w10x relaxation may be relaxing toward a wall-heated state, not the T_int one.
**How to apply:** diagnose the inner BC (p*v work of the oscillating wall) before trusting the rot-300 deep profile.

**ROOT CAUSE (09-27 night, wallE_0927, branch dhj-wall-0927 c84c8e57, not merged):** ix1_bc=user
(HydrostaticEquilibrium) never writes the inner fluid ghosts -> frozen t=0 reservoir; the is-face Riemann flux
exchanges +-3.5e20 g/s locally, inflow 0.5 % higher specific energy -> net 1.28e30 erg/s. Fix
`problem/wall_closed = true` (mirror state at the face, hydro+MHD, cs+sp): Etot_bot = Lrad_bot exactly,
Mdot_bot 0, box dE/dt 1.71e30 -> 3.9e29. BUT not the main cause of the 100-bar warming (+5.4 K vs +4.5 K / 3 rot).
ck r^2 NOT missing (ck_spherical, ck_beam_sph); 40 % starlight gap = transmitted limb light through the thin shell.
