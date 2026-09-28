---
name: red-giant-runaway-source-rt-stale-w0
description: "RUNAWAY SOURCE FOUND 09-11 02:30 (R1_scan, per-operator max-T ledger): a 3-operator loop -- the two-stream apply removes MORE energy than the cell has (de -1619 vs e 1517) because it evaluates ei, the positivity guard, the rescue and LimitRTSource on the STALE w0(IEN) (previous stage's ConToPrim) while applying de to u0 that RKUpdate just changed; ConToPrim's floor then CREATES energy (0.1 K cell next to 1e4 K at equal rho); the next RKUpdate dumps +1e3 into that contact -> 1e8-1e17 K. Conduction (implicit+cap) contributed EXACTLY 0; gravity/WB 1e-5. Fix in progress: rt_use_cons (RT reads e, rho, T from u0 at apply time)"
metadata:
  type: project
---
Instrument: src/utils/runaway_scan.hpp (new, header-only, <problem>/runaway_scan, runaway_rmin,
runaway_ratio, runaway_ratio_state), scan calls in red_giant.cpp (1886, 1949, 2782, 2796, 2839)
and hydro_tasks.cpp (after ImplicitConduction 534, after ConToPrim 690). Bit-identical on/off.
R1_scan (job 196042, gate 1e-16, clamp off, implicit+cap, rst 5.0e5, stopped at 5.787e5):
117 cells crossed T/T_nb > 10, all in 5.690-5.787e5. Sum de by operator over 2577 lines:
hydro RKUpdate +1032, RT +631 (huge oscillation), ConToPrim floors +180, gravity/WB -0.013,
implicit conduction +0.000. Traced cell gid 90 (9,3,391) r 3.62e12 (the I5 dt-collapse cell):
cycle 18698: hydro +960.9 (T 9.3e7 K) -> RT -1619.2 (e = -102) -> floor +102.3 (e 6.7e-7,
0.1 K) -> next cycle hydro +480.9 ... Neighbour rho all 1e-14..8e-14 (NOT a vacuum cell);
i-1 neighbour had non-positive e. efloor_de 1.5e8 -> 9.5e8 over the burst.
Code facts: two_stream_rt.hpp ei = w0(IEN) at 1865/1954/1997, guard `ei + de > 0` 1923,
rescue defl = -(1-1e-3)*ei, LimitRTSource(de, w0(IEN)) 1945/2404, T_g precompute from w0 at
995/1139/2213; the user source function (red_giant.cpp ~1880-2840) runs after RKUpdate/
SrcTerms/ImplicitConduction with w0 STALE ("THE SPONGE RUNS BEFORE ConsToPrim" note 2731).
The agent proposed rt_de_max=0.5; REJECTED: it clamps against the same stale w0.
Fix brief: EintFromCons helper shared with ImplicitRadialUpdate; <problem>/rt_use_cons makes
the T_g precompute, ei, deq, Newton, guard, rescue floor, limiter and rho read u0. Gates:
off -> bit-identical to pin5; on -> smooth-star agreement; R2_scan (crossings should vanish,
efloor_de flat) then I6 = I5 + rt_use_cons to 9e5.
Pins: athena_pin4 (scan, code-unit T), pin5 (scan + K units + runaway_ratio). Runs:
R1_scan/. See [[red-giant-i5-density-gate]], [[red-giant-implicit-radial-diffusion]].

## R2 RESULT 09-11 03:30: rt_use_cons CLOSES THE FLOOR-ENERGY STEP, loop reopens via MOMENTUM
R2_scan (pin6, rt_use_cons on): ConToPrim floor energy creation in the watched cells +182 -> -1e-8
(gone); crossings at 5.79e5 117 -> 59; first COLLAPSE 5.728e5 -> 5.867e5; but the new collapse
cell is a DENSITY-FLOOR cell (rho = dfloor 1e-18, r 3.641e12) with v = (4.2e11, -2.9e11, 5e11)
cm/s: general_c2p_hyd.hpp:44-48 floors rho "without changing momentum or energy" -> v inflated by
rho_old/dfloor. Fix (4) in progress: <hydro>/dfloor_keep_velocity (scale m by fv = rho_old/dfloor,
u.e -= (1-fv^3) KE_true; cs: metric-correct KE in GnomonicEquiangleRaiseVel), R3_scan, then I7.
Also done: src/utils/eint_from_cons.hpp (shared), 52 two-stream sites now read eiN/rhoN;
bc_use_cons (RedGiantBC ghosts from u0: smooth star dT 1.3e-6 in the atmosphere, L_out identical,
inner-face budget accumulator moves 5e-3 -- diagnostic only); rg_wallflux EOS guard (never fires
on the smooth star). Pins: pin6 (rt_use_cons), pin7 (+bc_use_cons, wallflux guard).
I6_fixes job 196082 (pin7, rt_use_cons + bc_use_cons, rst 5e5, tlim 9e5) RUNNING.
Queued after R3/I7: rt_bface switch (default off, red giant on) + regression comparison vs HEAD
for solar_convection / dhj + docs/handover/NOTE-2026-09-11-red-giant-switches.md.

## I6 RESULT 09-11 05:10: DIED 5.877e5 the R2 way (the flag was not yet in pin7)
I6_fixes (pin7 = rt_use_cons + bc_use_cons, NO dfloor_keep_velocity): dt COLLAPSE cycle 20122
t=587763 dt 1.56e-7, cell (2,9,9,342) gid 53 r 3.452e12 rho=1e-18 (=dfloor) T 8.3e6 K,
v=(-5e15, 5e16, -1.5e16) cm/s -> the momentum-kept floor cell exactly; then I_dn = nan in the
two-stream column and L_rad,out/L = nan; crawls at dt 0.03 afterwards (scancel is denied to the
assistant by the auto-mode classifier: the USER must `! scancel`). Confirms fix (4) is the
next gate: dfloor_keep_velocity is implemented in the tree (eos.hpp:62-71, ideal/general
c2p, hydro.cpp:136 deferred-floor path, pin8 built 02:59 probably has it) but UNVERIFIED;
agent audit + gates + R3_scan + I7 launched 05:15.
