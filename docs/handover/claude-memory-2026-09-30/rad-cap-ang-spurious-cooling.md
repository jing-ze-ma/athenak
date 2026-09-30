---
name: rad-cap-ang-spurious-cooling
description: 2026-09-13 THE TRANSVERSE CAP rad_cap_ang IS A SPURIOUS ENERGY SINK in the FeCZ boxes - narrow 336x8x8 1-rank arms (bench/bstar_fecz/deepcool/): cap 0.5 + vpert -> F_top/F_bot 1.21-1.24 plateau within 200 s, E -0.58 % in 600 s, box-wide <v1> -3e4 (reproduces rt2turn/rt2turn_relaxed's 3 % loss + subsidence); vpert 0 -> 1.000; cap 0 + rad_implicit_ang (vpert on) -> 1.000; 1-D + vpert -> 1.000. Not the wall, not MPI, no floors. Mechanism: unrelaxed horizontal T structure + nonlinear kappa(T,rho), F ~ T^3/kappa dT/dz raises the horizontally averaged flux; energy leaves through the top. EVERY capped box result (rt2turn convective fluxes, He smoke amplitudes) is compromised; cap inflated He KEh 3-4x. RULE: FeCZ boxes run rad_cap_ang = 0 + rad_implicit_ang = true (build_hip_tr).
metadata:
  type: project
---
vpert (box_convection.cpp 666-724): v1 only, compressive, 16 modes, zero at walls. bc_mode 3 wall (803-845): cancels mass
flux exactly; energy via an enthalpy estimate at the wall cell (deliberate, keeps rad_flux_inner); residual dm(h_ghost-h_cell)
for nonzero wall v1; silently skips when rho-dm<=0. Arms 09-13 night: rt2turn_tr 11699959 (operator, uncut), rt2turn_topcut
11700053 (operator, top at tau 2.16, rt_top_vacuum=false), rt2turn_relaxed 11699793 (cap reference). See
[[implicit-transverse-raddiff-plan]], [[fecz-column-relaxation]], [[rt-thin-region-lte-prad]].
