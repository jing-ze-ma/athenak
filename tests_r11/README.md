# tests_r11 (2026-09-18, Fable, no Opus): the top boundary is NOT the cause; a stronger
# closure is worse; the 3-D run to 5 turnovers is the next measurement

Binary athena_v7 = a6d66c3c + `problem/rt_top_vacuum` exposed in red_giant.cpp (it was never read
there: the box sets it true, the sphere always used the unresolved-column model). Arms are the
tests_r9 1-D configuration (r11_arm.sh = r9_arm.sh with the new binary, apudev).

| arm | change | L_out/L at 0.1 / 0.5 / 0.9 turnover | end |
| --- | --- | --- | --- |
| base11 | none | 1.278 (t=0) -> 0.981 / 0.978 / 0.970 | 1 turnover, no collapse |
| vac11 | rt_top_vacuum = true | 1.278 -> 0.983 / 0.979 / 0.970 | identical to base to 1e-4 |
| a15_11 | none, 2 turnovers | 0.980 / 0.977 / 0.960 (1.2) / 0.985 (1.6) | fluctuates 0.96-0.99 |
| a3_11 | mlt_alpha = 3.0, 2 turnovers | 0.979 / 0.984 / 0.941 (1.2) / **0.359 (1.6)** | outer shells fall inward at 4-8e6 cm/s, T at 0.97-1.01 R rises to 2e5 K, eint of the 0.97 R cell x25 |

Reading: on this star the column above the domain has dtau ~ kappa p_gas/g ~ 6e-3, so the top
boundary model is irrelevant (vacuum is a null). The 2-3 % shortfall of L_out is the power going
into the envelope (expansion work + stored radiation energy) of a layer at Gamma_rad ~ 1.0-1.05:
a 2 % flux error there is a 2 % net acceleration, i.e. ~1e6 cm/s per turnover, which is the
observed drift. A stronger shell-mean closure does not cure it: it removes radiative flux from
the top layers, which then contract and heat. Only resolved convection carrying its share can hold
this envelope; whether it starts within a few turnovers is a 3-D question. Submitted: r11_3d.sh
(job 11793509, apu, 4 h, one piece, rst every half turnover, tlim 5 turnovers, arm tests_r11/g3d).
Gates when it returns: alive, L_out/L, lnKE1/lnKEh, F_res/F_req and v_r rms/v_MLT from the last
mltfaces/rt_profile, shell drift per turnover.
The tests_r10 thick-top matrix numbers (0.67 at t=0 for thick tops, the thin-top mode-3 runaway)
concern the unresolved-column model with p = gas pressure and kappa_const and are a separate item.
