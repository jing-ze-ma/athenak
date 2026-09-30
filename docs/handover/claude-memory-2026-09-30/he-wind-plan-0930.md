---
name: he-wind-plan-0930
description: He presn porous runaway = Fe-bump wind; user GO 09-30 for steady-wind setup (step 1, branch he-wind-bc) + M1 positivity (step 2, branch m1-positivity) in parallel agents
metadata:
  type: project
---
09-30: 128x128 He presn (ad3d_128) = 64x64 to 5.5 tt (runaway onset 3.5 tt identical); NaN at t 26825 (5.7 tt);
last good rst 00011 (t 25850). dt drop = hydro radial CFL in channel cells at r/R 0.98, v_r 970 km/s outflow
(> v_esc 660), rho 0.5-4 % of shell mean. Interpretation: Fe-bump radiation-driven wind (~3e-5 Msun/yr full sphere,
my scaling). Setup cannot be steady: inner x1 closed wall, outer outflow at tau~1e-2 -> domain drains.
User GO 09-30 ~22:00: step 1 = mass-supplying inner BC + outer boundary far out + EOS/opacity tables below
rho 1e-14 + 1-D column steady-wind test (branch he-wind-bc, runs /viper/ptmp2/jinma/hepresn_wind_0930);
step 2 = positivity-preserving implicit transport on all faces (ap_hll/blend exist for x1 only), gas-eint positivity
limiter, energy-conserving BE floor (branch m1-positivity, runs /viper/ptmp2/jinma/m1pos_0930); decisive gate =
restart ad3d_128 rst 00011 past 26825 to >= 7 tt. Step 3 (3-D 64x64 from t=0 with both) after the gates, needs user.
22:25: 64x64 hllc (ad3d_hllc, 12037045) also failed: L_top/L unphysical from t 29,200 (6.2 tt, Picard 150+), NaN
t 31,626 (6.7 tt); last clean rst 00012 (t 28,200); passed to the positivity agent as 2nd gate case. Mass loss
had slowed (-6.6 % at 5.5 tt, -6.9 % at 6.0 tt) before the solver failure.
23:10 m1-positivity REPORT (branch 38500817, unpushed beyond 33150622; binary m1pos_0930/bin/athena_gpu_he_star_m1_7192a176):
first trigger = Picard pass-0 solve E<0 from lagged -c dt v g0 face term (g0 = rho(kE E - kP aT^4), cancellation in
stiff rad-dominated plume cells) + opac-Newton rhs; NOT transverse advection (ap_hll x2/x3 not built). Keys (off=bitwise,
verified IDENT wedge+box): implicit_g0_exchange, implicit_g0_limit=1 (pass 0 only), guard_mode=6 (row dominance),
implicit_pos_gas, implicit_pos_floor. He box GPU cost: no change (agent's numbers 25.0-26.0 vs 25.0-25.3 s).
Controls off128/be128 NaN at 26,825 s; ON arms past it (onnh64 past 31,626 s) but 320-680 hesdirk2 stages
NOT ADMISSIBLE -> BE redo (vimp on AND off), Picard 100-175. Open for user: vimp off in relaxation; defaults;
hydro-side eint<=0 in near-void cells (rho/<rho> 2e-5). He box accuracy job 12046701 queued (ana.sh 60000 74000).
10-01 00:20: beon128 envelope COLLAPSED (M(r>0.7R) 71%->1%), not steady; slices https://claude.ai/artifact/HcgYDVQ3CozhxXavZppqcU
