# cs-tan-1106: full Eddington tensor on the cubed sphere (running report)

Branch `cs-tan-seam-1106` (f5d002cad = cs-tan-1106 + seamE 55db30f95; RSG arms ran 5b03746b2 = + seamE cc81d3159). Not merged, no defaults changed, no production touched. Key: `rad_m1/vet_gd_tangential = true` (read only when named; false on cs if absent).

## Why it was off
Energy budget, He 128^2 sp wedge, deliberate, from Viper's session log (docs/handover/NOTE-2026-11-07-viper-caltech-vet-gd-tangential.md on bsg-files-1009): no-cap full tensor R = +1.34, L_top +127 %, source not diagnosed. Not reproduced here so far (below), but not on Viper's state.

## Design
Q (trace-free tangential part) in LAT4/LAT5 about a per-cell Cartesian frame exchanged once with the slots; no new tensor slots; memory +6 reals per column. Face term n.div(QE) = d_n(n.Q.n E) + d_t(t.Q.n E) with n, t frozen at the face midpoint. Operator (19-point) form for implicit part; fixed point is the lagged second-order form.

## RSG real state, pass 1: rst 1295, global fallback, cs_cross_mono + guard_mode 7, 3000 cycles = 2.97e6 s, on vs off, 4 GPUs each (MEASURED)
- FATAL 0/0, NON-CONVERGED 0/0, NOT POSITIVE 709 on / 695 off, NOT ADMISSIBLE 357/353, BE steps 736 on.
- Fallback: 2537 on / 2538 off in 3000 cycles = 0.85 per cycle, independent of the key. Q active in 41561 of 54827 assemblies = 75.8 %.
- Energy gate: R/int L_in +8.29e-5 on, +9.18e-5 off. Slope/L_in all +9.2e-5 on, +1.06e-4 off; last half +8.95e-5 on, +9.10e-5 off. R_on-R_off slope -1.4e-6 of L_in. Background drift present in both.
- Cost: s/cycle 2.94 on vs 2.74 off (+7 %); swap pass same sign on swapped nodes (elapsed at cycle 50: 249 s on, 230 s off). Measured on shared stretches of the same pair, nodes differ.
- Q rms: tau_top 0.01-0.1 0.085, 0.1-1 0.076, 1-3 0.053, 3-10 0.034, 10-30 0.020.
- Flux at cycle 592500, F1/(cE)>0.25 band (tau ~ 1 PROXY): see real/flux_592500.txt; seams and depth: real/diag_592500.txt. No reference solution: these are changes, not improvements. Gas divergence between arms grows with time; cycle 590000 is the cleanest direct measure (median dE/E 1.2 %, dF1/|F| 1.8 %, drho/rho 0.015 %).

## RSG: why its fallback rate is high (MEASURED, fbdiag of rsg2 arms, 200 triggers each)
E<=0 cells sit deep: median tau_top 425 (p10 8, p90 1400), r/R_ph 0.94-0.99; only 19-22 % above tau_top 30 (the vet_gd cut). 16-20 % still solve to E<=0 with D_r,lat and vimp masked off locally. HG l1: median tau_top 3.7, r/R_ph 0.996. So RSG's fallbacks come mostly from below the tangential-tensor region and are not caused by the key (2537 on vs 2538 off in pass 1). Separate problem; not addressed on this branch.
Local fallback (implicit_fb_local) on RSG: m1_vlat_fall at cycle 590020 111 on / 107 off vs 199 / 186 with the global path: ~45 % fewer global drops.

## HG real state vs seam's l1 (key off reference, same recipe, rst 139, same 6 GPUs) (MEASURED)
- 1550 cycles to tlim 5.0285e7 (l1 1200 cycles to 5.0245e7). FATAL 0/0, NON-CONVERGED 0/0, NOT POSITIVE 8 vs 12, NOT ADMISSIBLE 12 vs 7.
- Q active in 23031 of 23102 transverse assemblies (99.7 %); vet_col_lat fallbacks 10 (l1 19 in 1200 cycles).
- Energy gate: R/int L_in -3.9e-6 on, -5.2e-6 off; slopes/L_in -3.7e-6 / -5.3e-6 (all), -7.3e-6 / -6.7e-6 (last half). R_on-R_off slope +3e-6 of L_in.
- Cost 2.44 vs 2.31 s/cycle (+6 %), same GPUs, different times.
- Q rms: tau_top 0.01-0.1 0.067, 0.1-1 0.064, 1-3 0.051, 3-10 0.037, 10-30 0.019.
- Flux at the same time 5.0227e7 (115000 s after split), band = F1/(cE)>0.25 (tau ~ 1 proxy): median |dE/E| 2.7 %, |dF1|/|F| 2.5-2.7 %, |dFt|/|F| 3.8 %, p99 dF1 0.6-0.9; drho/rho median 0.17 % (gas divergence ~15x smaller). Top shell sum F1 -0.08 %, T_eff map rms 0.47 %, p99 1.5 %. Seam cells 2.9 % vs interior 2.7 % vs vertex 2.3 %. Shell mean dE/E +0.1..+0.7 % 36-94 cells outside the band; deep shells <= 3e-6. Files real/hg_flux_01744.txt, real/hg_diag_626053.txt.

## sp checks (MEASURED, short)
- R136a1 rst 5 copy, 2 ranks on 1 GPU: key on 600 cycles; key off 300 cycles (stopped for priority). Common 1504 s: FATAL 0/0, L_top within -0.44..+0.09 %, R/int L_in -3.6e-5 on, -1.8e-4 off; slopes of opposite sign (noise). NOT POSITIVE 2 on / 0 off, NOT ADMISSIBLE 1 / 0. Viper's sp source not reproduced in this short window.
- BSG column 256x8x8 seeded 0.05, 6000 cycles: R_on-R_off +1e-6 of int L_in dt, L_top <= 9e-6. Q small there. The first unseeded column was laterally uniform (vacuous).

## Gates
- keys-absent bitwise IDENTICAL: CPU HG hgb +4, accretor +2, RSG rsb +4, dhj W121 +4 (f59b131fa); GPU HG hgb +4 (f59b131fa and the merged 5b03746b2). Not run: GPU accretor/RSG/dhj MHD/R136a1 bitwise, CPU R136a1 (run started, not checked), layout invariance of the key (cases written: he_hg_tan, rsg_tan, he_col_gd_tan, not run), beam/thin-sphere accuracy test.
- MMS (cs, CPU, C16-C128): lagged Q face term 2nd order on interior, seam and vertex faces (1.7e-2 -> 3.1e-4); operator form 2nd order interior, 1st order at seams (by design, the fixed point is the lagged form).
- implicit_op_check with the key: operator variants agree (<= 7e-16); cs_label FAIL pre-existing on the base binary (fixed on seamE).

## Accuracy vs the formal solution, HG real state (MEASURED, 02:20-02:45)
l1 recipe, rst 139, 6 GPUs, 160 cycles key off and key on. Reference: the vet_gd short-characteristics J of the same build (output6 -> m1_thsw, m1_tj). Band = F1/(cE)>0.25 proxy.
| offset from band | off med\|E-J\|/J | on med\|E-J\|/J | on/off |
|---|---|---|---|
| -40 | 1.7e-4 | 1.7e-4 | 0.996 |
| -10 | 1.94e-2 | 2.00e-2 | 1.03 |
| -3 | 0.103 | 0.106 | 1.035 |
| 0 | 0.132 | 0.136 | 1.036 |
| +3 | 0.140 | 0.143 | 1.02 |
| +20 | 0.154 | 0.155 | 1.004 |
Same at an earlier output. The full tensor does NOT reduce E's deviation from the SC J; slightly larger near tau ~ 1. Small HG state (Q rms 1e-4): ratio 1.000.
Caveats: reference approximate (42 directions, ray noise; the hc thin correction already pulls E toward it); E only, no reference flux H output; Q enters the rotation average from zero after the restart (8 sets x ~18 cycles), so it may not be fully in at cycle 160. File real/hg_acc_n121.txt. No 2-resolution study of one state and no beam test: not done.

## RSG pass 2 (rsg2): local fallback + cs_cross_mono + guard 7 + fb_diag, rst 1295, on (02-09) vs off (02-16), stopped 04:00 by exact PID (MEASURED)
- 3700 cycles on / 4000 off (3.7-4.0e6 s). FATAL 0/0, NON-CONVERGED 0/0, NOT POSITIVE 851 / 871, NOT ADMISSIBLE 422 / 484.
- Fallback: m1_vlat_fall 1703 on (0.47/cycle) / 1732 off (0.44/cycle) vs 0.85 with the global path. FBDIAG_SUM stage2 passes with vlat active: on 15935/16593 (96.0 %), off 16021/16726 (95.8 %). (Finalize line not printed: no graceful stop; FBDIAG_SUM used instead.)
- Floors: m1_pos_flr 3044 on / 3031 off, ~2.6x pass 1 (1173): the local mask costs floors in both arms.
- Energy gate: R/int L_in +7.6e-5 on, +6.6e-5 off; slope/L_in all +8.5e-5 / +7.7e-5, last half +6.8e-5 / +6.2e-5; R_on-R_off slope -1.0e-6 of L_in. Same background drift as pass 1.
- Cost: 2.97 vs 2.75 s/cycle (+8 %) with the nodes swapped relative to pass 1: the +7-8 % is the key, not the node.
- Flux at cycle 593500 (band proxy): median |dE/E| 2.3 % (3 inside), 1.9 % (band), 1.3 % (3 outside); |dF1|/|F| 1.9 / 2.0 / 1.3 %; drho/rho 0.6-0.7 % (gas diverged by now). T_eff map rms 0.23 %, p99 0.73 %, top shell F1 +0.28 %. File real/rsg2_flux_06333.txt.

## Summary for the user
- Why off: Viper's He sp energy budget (source not diagnosed). Not reproduced on cs HG/RSG (budget difference <= 1e-5 of L_in over 0.17-4e6 s) nor in short sp R136a1 / BSG-column runs; Viper's He state not tested here.
- Stable on real HG and RSG states; Q at tau_top 1-3 rms ~0.05; active 99.7 % on HG, 76 % (global fallback) or ~96 % of passes (local fallback) on RSG.
- Cost +6-8 %; memory +6 reals per column.
- Effect: ~2-3 % median local change in E and F near tau ~ 1, p99 15-60 %; T_eff map rms 0.2-0.5 %; no systematic L shift beyond ~0.4 %.
- Accuracy: against the SC J of the same build, E is 2-4 % (relative) further from J with the key on near tau ~ 1. No evidence the full tensor is more accurate; flux reference not available. A beam/thin-sphere test with an analytic reference is still owed.
- Decisions: (1) adopt or not, given no shown accuracy gain; (2) if pursued, a flux-reference test (SC H output) and a Viper He-state rerun; (3) RSG fallbacks come from deep cells (median tau_top ~425), unrelated to the key: separate issue.
- Seam agent's report (a2ebd3e33b4efd30b) asked to be saved as seamE_1106/REPORT.md; not written by this agent.
