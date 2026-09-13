---
name: rt-handover-semi-implicit-bug
description: 2026-09-13 FOUND + FIXED (decad388 on rg-box-convection in bench/wt_rgbox, NOT pushed, NOT in polar-average-perf) - the two-stream RT tau-blend HANDOVER term ((w_{i+1}F_{i+1}-w_iF_i)/dx, two_stream_rt.hpp ~2165) went through the semi-implicit relaxation (damped by (1-e^-x)/x per cell) while conduction took its half in full (conduction.cpp 548/1043) -> the two halves no longer cancel and ~7-13 % of F is deposited for ever inside the blend ramp. SHARED by the spherical path - every red_giant and deep_hot_jupiter_rt run to date had it (red_giant.cpp 983 makes the blend mandatory; rt_semi_implicit defaults true); at nx1=480 the ramp is 2-3 cells so the whole handover lands in one damped cell near tau 10-100 (photosphere region). Invisible in tot-E (7e-9 of the red giant's E over 3e6 s). Fix - split src into the (1-w)-weighted local balance (relaxed) and the handover (explicit, after Newton/rescue); inert where w=0
metadata:
  type: project
---
Conservation test inputs/hydro/box_convection_grey_cons.athinput (tau_bottom 1000, blend 10/100, starts from the analytic Eddington
column via tools/box_convection/mk_eddington_column.py): dE/dt 0.098 F -> 0.0064 F at 4e3 s, 0.020 F -> 0.0005 F at 4e4 s
(explicit-source reference 0.0057 F). Grey atmosphere now hits the EXACT grey surface T/Teff 0.8102 (exact 0.8112; Eddington
0.8409) and slope 0.7506 (exact 0.75): the residual 3 % vs the Eddington formula is the FORMULA's closure error. Spherical path
NOT bitwise after the fix (red giant 200 cycles: KE_r 0.6 %, momenta 4 %). TODO for the user - port decad388 to
polar-average-perf; re-check the red-giant photosphere/corona structure (the deposit sat at its tau 10-100) and the dhj RT.
Rad tests 4/4. See [[rt-plane-parallel-port]], [[fecz-box-projects]].
