---
name: fecz-column-relaxation
description: 2026-09-13 FeCZ box initial columns vs the two-stream - rt2turn (B star, handover 100/300) lost 3.1 % of E in 3500 s (15-20 % of F*A) through the TOP as excess emergent flux, drained from the deep interior, global subsidence -3e4 cm/s; not a wall or ramp leak (rt2turn/analysis/eloss.txt). He 3 Msun column only 1.3-2 % off and self-heals in 0.03 turnover; relaxed IC ic/ic_fecz_he3_rt_relaxed.txt (one-line ic_profile change in he3_2turn). 1-D He column is ACOUSTICALLY UNSTABLE (growing ~90 s mode, Fe-bump kappa mechanism) - no stationary 1-D state; relax by v-annihilation segments (hestar_fecz/relax1d_3msun/relax_loop.sh). The formal-solution F_top diagnostic (eloss_profiles.py) is BIASED a few % when dtau/cell ~1; the hst energy budget 1-(dE/dt)/(A F) is authoritative.
metadata:
  type: project
---
B-STAR DONE (bstar_fecz/relax1d/, 1e5 s 1-D, 32 min GPU): energy-budget F_top/F_bot 1.0156 at t=0 (formal solution said 1.089 - biased), 1.00002 at 1e5 s; column rings NOT (breathing mode decays on its own, single run). Relaxed IC ic/ic_fecz_rt_relaxed.txt: fix is photosphere and above (tau<=1 was 15 % too hot, 2-3x too dense; tau 1e-2 T x0.92 rho x0.32; tau 1 T x0.85 rho x0.56; FeCZ/base <1.5 %). Verified 0.9985 over 100 cycles. rt2turn FINISHED 2 turnovers in 81 min (1 link): E recovered to -0.9 %, KEh 5.4e33 still rising. BOX RESTART BUG: restarting a box RT run from rst corrupts the state at cycle 0 (rt_de_max clipping, |v1| 5e6) - diagnosis agent launched 09-13 evening; rt2turn's chain relies on restarts. Recipe: nx2=nx3=1, vpert 0,
box_convection.cpp has NO velocity damping (cool_layer would defeat the purpose); 16 x 25 s segments restarted from the dumped
column with v=0 (analysis/make_ic.py). He relaxed column: tau 1e-2 dT -3.8 % drho -15 %, tau 1-10 <1 %, tau 100 dT -2 %
drho +6 %, bottom drho +18 %; start-up KE 20x lower. smoke_rt2 (handover 10/100) was flat because the diffusion operator,
which the IC was integrated with, owned everything below tau 100. See [[fecz-box-projects]], [[rt-thin-region-lte-prad]].
