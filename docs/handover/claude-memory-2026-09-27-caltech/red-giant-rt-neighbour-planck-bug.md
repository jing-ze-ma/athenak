---
name: red-giant-rt-neighbour-planck-bug
description: 2026-09-10 08:15 — TWO-STREAM RT BUG (all three sweep kernels, two_stream_rt.hpp ~1231/1463/1583 + the alp/bet/gm layer source): a cell's emission term Em and the layer source function use the NEIGHBOUR's Planck function B(i+1)/B(i-1) regardless of the neighbour's opacity; next to a hot transparent corona the last active cell is charged the CORONA's emission (1e9x its own) and is driven to the T floor in ONE step (the "Newton rescue 99.9% drop"); next to any hotter cell (a 1-cell ionization front, 8983/2874 K -> 48x) a cold cell over-cools. Candidate root of the whole photospheric cold-cell-collapse family (V4, V7, B4/B5, T6, R7)
metadata:
  type: project
---
Captured at t=0 in V8r_rtdiag (job 195456, <problem>/rt_cell_report): cell i=439 (last active
below the join) T 3364 K, kap*rho 2.4e-27, I_up 2.85e9, B(i) 2.31e9, B(i+1) 2.34e18 (6e5 K,
kap 0): abs_up 3.6e-17 correct, but emi 1.49e-8 = 4 pi kap rho 0.5 B(i+1) -> src -3.6e-8,
e_eq 4.7e-12 vs e 7.5e-10, Newton crosses -e, rescue de = -0.999 e; T -> 8-34 K by cycle 50
and can NEVER recover (emission fixed by the corona's T). With B(i) instead: equilibrium
T = 3364/2^0.25 = 2828 K, radiative time 3.3e7 s = the physical answer.
Rescue itself is innocent (solves the poisoned equation). Also: rt_top_re's back-radiation is
multiplied by (1-exp(-kappa_above p/g)) = 0 whenever rad_kappa_above=0 -> silently defeated.
FIX (implementing 08:20, agent): weight the layer source and Em by kappa*rho of the two cells
(reduces to B_own for an inert neighbour, identity for equal opacities), in ALL THREE kernels
([[two-stream-rt-three-kernels-trap]]); optionally clamp the implicit step to the local
radiative-equilibrium energy from the actual incident intensity instead of the fixed 99.9%.
GATES: the grey 1-D validation must still match; then B13 (rg.00124, fixed RT, NO opac_tmin:
does the cold-collapse death at 3.2256e7 vanish?), V8 rerun (join cell ~2800 K, corona holds?),
V7 rerun (open top from scratch).
Related: [[red-giant-corona-join-cell-cooling]], [[red-giant-molecular-opacity-knee-runaway]],
[[red-giant-top-cooling-runaway]], [[red-giant-open-top-prod12-recipe-fails]].

## FIXED AND GATED 2026-09-10 09:00 (agent a369a053; BFace hunk read by me; 16 call sites)
two_stream_rt.hpp: BFace(k_own,k_far,b_own,b_far): returns b_far (old expression, bit for bit)
unless k_far < 0.1 k_own, then blends to the kappa*rho-weighted mean, b_own for an inert far
cell. Applied at every two-point-B site in all three kernels (grey/ck/generic sweeps, Em,
rt_top_re probe). rt_rescue_eq (default true): a would-be e<=0 Newton step lands on the local
equilibrium (Em(T_eq)=A) instead of the fixed 99.9% drop; counters split eq/floor. rt_top_re
now warns when rad_kappa_above=0 (its back-radiation is then exactly zero). rt_cell_report
diagnostic kept (off). Binary pinned: red_giant/build_guard_pin/athena.
GATE 1: smooth star t=0 column BIT-IDENTICAL (321 rows, F_out 1.22780715e10). At the V8 join
the net face flux changes 1e-8 (the fake emission telescoped out of every FLUX diagnostic),
the cell source changes 7 orders (de/e -0.999 -> -6.6e-7). NOTE: an UNCONDITIONAL kappa-
weighted variant moved the smooth column's F_out by 3.47% (1.2278e10 -> 1.2704e10) -> the
two-point average is baked into the validated grey solution; which is closer to the analytic
Eddington atmosphere is UNCHECKED (requested 09:05).
GATE 2: B13_rtfix (rg.00124, RT fix, NO opac_tmin, NO wb_rmax, WB-walk fix): 3.205e7 -> 3.24e7
dt flat 30.141; cold cells 2/2/0/12/27 at 3.210/3.218/3.224/3.230/3.240e7 vs B4 4/85/156(dead)
and B12 (clamp) 113-129; nT<1000 = 0, floor marks 0. => the cold-cell collapse family in the
lidded run WAS this RT bug; the opacity clamp is no longer needed there.
GATE 3: V8f (6e5 K corona): join cell 3274-3384 K through 1.5e5 (V8: 31 K in 1.5e3 s), 0
rescues in 6600 cycles, atmosphere expands INTO the corona (<v_r> +3e4..+1.2e5), corona T
6.0e5 -> 4.2e5 as it loads stellar gas, L_out/L 0.90 -> 0.80 (V8 0.65). Still running to 6e5.
ANALYTIC CHECK 09:15: smooth star t=0 (L 1.914e36, T_eff 4000 K, F_expected 1.229268e10):
old/surgical F_out/F_exp = 0.99881 (bit-identical to each other); the UNCONDITIONAL kappa-
weighted form = 1.03348 (+3.3%, 28x worse; it pulls every layer's source toward its hotter
deeper endpoint since kappa falls outward) and worsens the tau 1-1000 residuals 2-8x. The
thin-top (tau<1) residual of ~0.15 is the hemispheric quadrature, common to all. => SHIP THE
SURGICAL FORM (default RT_BFACE_R = 0.1).
V8f CLOSE-OUT (agent, after the session exit): clean to 3.06e5 (past V8's 2.96e5), join cell
3270-3390 K throughout, 0 rescues. Then a DIFFERENT failure at ~3.1-3.2e5: ONE near-vacuum cell
INSIDE the star (i=411, 27 cells below the join, rho 6.4e-26, T 2.06e9 K between 4499 K and
1043 K neighbours) -> its OWN Planck function 3e32, A < 0, so rt_rescue_eq is inert (needs
A > 0) and the floor rescue fires; dt 30.6 -> 0.87. = the vacuum/hydro trigger family, not
the neighbour-Planck bug. Cold cells T<1500 creep 1 -> 30 by 3.2e5. V8f hit its 2 h wall at
3.307e5. OPEN for the open-top/corona configs: what evacuates a cell to 6e-26 / 2e9 K in the
star's thin atmosphere (rho 1e-18 floor? no: dfloor was 1e-26 in V8f; V9f/V10 have 1e-18 —
check whether THEY show it).
