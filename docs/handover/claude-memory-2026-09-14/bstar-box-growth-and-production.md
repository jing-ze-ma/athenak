---
name: bstar-box-growth-and-production
description: 2026-09-13 B-STAR FeCZ BOX, clean (operator, relaxed column) - convection is REAL but radiatively damped: KEh e-fold 3.3 turnovers (18x slower than adiabatic BV), at 2 turnovers KEh and F_conv are 1e-3 of MLT; saturation ~11-13 turnovers -> PRODUCTION NEEDS 15-20 TURNOVERS (1.5-2e5 s) = ~11 h on 2 GPUs (6 x 2 h links, restart verified) or <3 h on 8 GPUs; 10x seed amplitude would save ~7 turnovers. Uncut vs top-cut (tau 2.16): interior identical to 0.3 % (FeCZ velocities/flux 1-7 %), the thin top only adds motion above and in the 1.5e9 cm below the cut -> UNCUT = production configuration (keeps the photosphere = the observable), but its top layer has the LTE-Prad EOS defect (heat capacity 17x) -> surface signals untrustworthy until fixed. PRODUCTION CONFIG FINAL 09-14: prod_taper LAUNCHED 09-14 job 11701870 (safety 11701871, 3 h) (taper+force, operator, FIXED sweep + cut BC + limiter + solver, column ic_fecz_rt_relaxed_fix_taper.txt (F_top/F_bot 1 to 3e-8), vpert 1e-2, build_hip_fix, 7.5 h) dt-forcing floor 0.6 % of v_MLT. The untapered prod_uncut was CANCELLED 09-14 at 6.9 turnovers (old sweep, phantom top; its 7 turnovers of data remain in bench/bstar_fecz/prod_uncut). Old note: bench/bstar_fecz/prod_uncut, job 11700420 (ONE 12.5 h link, safety 11700421, 2 GPUs, plus one 3 h safety link), vpert 1e-2 (3x seed = 0.3 v_MLT; 3e-2 would have been ~v_MLT, cancelled), tlim 2e5 s, bin every 1e4 s (0.63 GB), rst every 1e4 s, binary build_hip_prod; ETA ~11 h.
metadata:
  type: project
---
Runs: rt2turn_relaxed (cap reference, 11699793), rt2turn_tr (operator, uncut, 11699959, 2 turnovers in ~65 min),
rt2turn_topcut (operator, top tau 2.16, 11700053, 2 turnovers in 36 min, dt 1.3-1.5 s). Analyses: rt2turn_topcut/analysis/
{growth,compare_cut}.*, rt2turn/analysis/turnover2.*, rt2turn_relaxed/analysis/early.*. Production binary
bench/wt_rgbox/build_hip_prod/src/athena (HEAD 84b19667, restart bitwise in dens/ener; mom 1e-7 residual). Handover 100/300
is below the B-star FeCZ (tau 580-3300): fine. Thin-region fix design: above the handover drop aT^4 from the EOS and apply
kappa F/c from the two-stream as a momentum source, tapered across the blend. See [[rad-cap-ang-spurious-cooling]],
[[rt-thin-region-lte-prad]], [[fecz-column-relaxation]], [[implicit-transverse-raddiff-plan]].
