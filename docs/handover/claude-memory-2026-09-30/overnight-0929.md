---
name: overnight-0929
description: User asleep from 09-29 ~09:30: scope of autonomous work (He presn 3-D scout only), what needs the user, morning summary
metadata:
  type: project
---
User went to bed 09-29 ~09:30 ("can I trust you with this").
ALLOWED: 3-D scout agent (grid refine, column gate, GPU build, smoke, submit approved R9 relaxation cfl 0.9 ~10 turnovers,
1 apu node) in /viper/ptmp2/jinma/hepresn3d_0929; verify its numbers; stop a clearly broken run and diagnose; commits on
he-presn-m1 only; push handover notes only; 10-min watcher on R9.
NOT without the user: WASP-121b jobs (12018387/89 + afterany 12018388/90 untouched), merges to rt-integration, default
changes, new threads, order-gate arms, production. Queue decisions for the morning.
Morning: one short summary (WASP-121b, scout state + first convection diagnostics, decision list).

**~10:05 R9 FAILED (verified):** job 12021350 (apudev) blew off: M_tot 3.40e25 -> 1.66e25 g by t 9000 s (through the top),
NaN at t 16677 s, killed. No convection (v_r'/v_MLT 1-3e-2). Cause (agent): FeCZ thermal time ~500 s << convective growth
~5000 s; without the frozen MLT flux the Gamma~1 FeCZ heats and unbinds the envelope. Grid B184 (nx1 184, tanh bumps in
StretchRPoly, default off/bitwise), seed on eint+E (he_seed_nk/he_seed_rad), commits c3fb3795 2e817ac7, binaries
athena_he_{cpu,gpu72}_c3fb3795; smoke clean, dt 97 s at cfl 0.9, 0.35 s/cycle on 2 GPUs. NOTHING RESUBMITTED (overnight
rule). DECISION FOR THE USER: (1, agent-recommended) keep frozen MLT (esrc) on while convection grows, ramp off over a
few turnovers (new pgen key); (2) velocity seed 0.1-0.3 v_MLT; (3) accept blow-off with top sponge. Caveat: the
two-stream era's MLT sub-grid flux + ramp arm D (09-17) also drained, but that was with the two-stream defects.
Also open: wall-cell push, photosphere 8e-3, thin top 0.12, interior IC residual 1.4-1.6e-3.

**USER 09-29 ~10:15: option 1 chosen** (frozen MLT on while convection grows, then time-based cosine ramp off: keys mlt_ramp_start/mlt_ramp_time). Opus agent running it in hepresn3d_0929/M1; it stops and reports if it blows off again (no options 2/3 without a go).
