---
name: mhd-relax-hydro-spinup-ok
description: 09-25 prod4 vs hyd4: at bbot=3 MHD is indistinguishable from hydro after rot 6; field establishes in ~10 rot (20 for >100 bar) -> hydro spin-up then C256 MHD for 10-20 rot is enough
metadata:
  type: project
---
/viper/ptmp2/jinma/mhdrelax_0924/README.md (prod4 to rot 80.8, hyd4 to 99.1; same IC, field on from t=0, bbot=3,
eta from EOS, max_eta 5e12). ME: x100 in 0.5 rot, peak 3.7e33 at rot 10, then wanders x4 (no steady level); 55-75 % at
10-100 bar; >100 bar band fills in ~20 rot; bands above 1 bar decay steadily. MHD-hydro: no difference > 1.1 noise
units in u_eq, T_day/night, zonal T at 1e-6..100 bar after rot 6 (only a rot 2-6 wind-up transient, +-38 K night side
aloft). The KE jump at rot 58-70 is in both runs (hydrodynamic). Effective drag time (H/v_A)(u/v_A) ~10-20 rot at
0.1-10 bar; v_A/u <= 0.06.

**Why:** decides the production plan "long low-res hydro spin-up, remap to C256, switch on the analytic field".
**How to apply:** final MHD phase N = 10-20 rot suffices at bbot=3; caveat: not tested as a switch-on into a spun-up
flow, and C256 changes Rm. Related: [[next-prod-ck-c2]].
