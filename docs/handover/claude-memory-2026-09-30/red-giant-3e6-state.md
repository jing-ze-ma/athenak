---
name: red-giant-3e6-state
description: RG_fofc_long2 at 3e6 (2026-09-13, bench/RG_fofc_long2/analysis): numerically sound (tot-E 1e-6, mass exact, dt flat, eos_fail=0, efloor_de 5e-38 of E, seams/vertices at rank 50, no hot cells) BUT (1) an EXPONENTIALLY GROWING small-scale radial mode at r < 1e12 (i 5-40, above the inner WALL), e-fold 2.7e5 s, growing since t~5e5, 44x in KE over the link, now 68% of the star's KE, not a pulsation (<v_r>/rms 0.01-0.08); (2) the open top + no sponge DRAINED the corona inward 150x (M>4e12 down 3.6e4x, 500x more went down than out the top), 0.45-3% of cells at the density floor at r>3.83e12. prod12 gate NOT passed on physics
metadata:
  type: project
---

KE_r by shell: r<1e12 2.7e39 (9e5) -> 5.9e43 (3e6) = 100% of the rise; all other shells flat; corona KE fell.
rms v_r at i=20 (r 2.9e11): 3.7e3 -> 1.4e5 (2e6 -> 3e6), i=1 (wall cell) pinned. face-budget "gain/L inner" flips
sign wildly = the wall being hammered; hst mass exactly constant so nothing crosses it. Corona: M(>join) 3.5e28 ->
2.3e26, top-shell median rho 2.7e-16 -> 6e-18 (~6x dfloor), mean v_r INWARD at the top; only 6.3e25 g left the domain
(saturated at 2.16e6). FOFC/tclamp episodes: 1.71-1.73e6, 1.93-2.79e6 (peak 7.6e7 fofc, 9.4e6 tclamp at 2.6e6),
2.91-3.0e6 (rising at tlim); floor cells only at r > 3.83e12, never at seams/vertices. No hydro_fofc dumps in this
run (add <output5> hydro_fofc for the next). **Next:** the inner-wall mode (inner_bc = wall at 0.03 R*, the RCB is at
0.04 R*: is it the MLT/conduction/wall interaction? the T structure is unchanged, only v grew) and the corona drain
(open top without sponge: [[red-giant-vertex-chimney]] decided sponge off; a floor-supported top is the price).

**Refinement:** on the prod12 gate's OWN terms (V7_prod12gate died 5.8-6.0e5 of a COLD collapse of the thin atmosphere:
31.6 K floor cells at i 298-316 blown at 1.4e7 cm/s, dt 30 -> 0.2 s, plus the vertex chimney) RG_fofc_long2 PASSES:
zero floor cells below i 427, eos_tfloor 0, max |v| 5.7e6, dt flat, vertex rho ranks 52-60 at i 300-340. FOFC + the
EOS clamp bought that. The two remaining defects (inner-wall radial mode, corona drained) are OUTSIDE what the gate
tests and are the next physics questions, not numerics.

**INNER MODE DIAGNOSED (bench/RG_fofc_long2/mode/parta.txt, 2026-09-13):** NUMERICAL OVERSTABILITY, not convection.
Background static (rho, e at i 0-80 unchanged <0.5% from 9e5 to 3e6 and = the MLT column); real-EOS nabla 0.39998 vs
grad_ad 0.399983 (superadiabatic by <=1e-5) => fastest convective e-fold 2.2e6 s, 8x SLOWER than measured 2.7e5 s;
one global rate at i 5/20/40 (3.35/2.76/2.66e5 s) while the dynamical time varies 5x; amplitude 5x the MLT velocity,
enthalpy flux 260-411x L/4pi r^2 (energetically impossible); smooth and large scale (<1% variance above 0.75 k_max,
radial correlation +0.91..1.0); buoyancy-driven (C(v_r,T') 0.35-0.94) => entropy created in phase with v_r.
Wall flux = exactly L, fixed, no feedback; deep radiative diffusivity gives no thermal restoring force. Suspects: the WB
scheme (background cached every 10 cycles = delayed feedback on a buoyancy oscillation) and the implicit radial
tridiagonal. Experiment mode/ (11658246): restart 1.0e6 -> 1.3e6, arms base / implicit_off / wb_off / wallflux_off /
cond_off (rad_kappa_rmax gate; rad_kappa_fac=0 is NaN and the pgen overwrites it), growth.py fits the e-fold.
**DRIVER FOUND (mode/ arms, 2026-09-13 ~01:30): the WELL-BALANCED SCHEME.** e-fold at i=20 (3e5-s window from 1.0e6):
base 3.13e5 s (growth x2.20), wbcache0 (background rebuilt every stage) 3.13e5 (x2.20, identical), wallflux_off 2.97e5
(x1.96, still grows), **wb_off (wellbalance_dynamic=false wb_x1=false) 7.17e6 s (x1.05, GONE)**; cond_off pending.
So the WB reconstruction/source (polytropic closure, wb_option) in the deep RADIATIVE zone above the wall does work on
buoyancy oscillations; the cache delay is innocent. Fix candidates: restrict WB radially (wb_rmax / a window) to
exclude r < ~1e12, or a closure that matches the radiative profile; the atmosphere still needs WB (33-58x residual).
**Complete table (mode/, e-fold at i=20 over 1.0e6-1.3e6):** base 3.130e5 s; wbcache0 3.132e5; wallflux_off 3.130e5;
cond_off (ALL radiative conduction gated off, grey two-stream inert too) 3.130e5 -- IDENTICAL; wb_off 7.17e6 (gone).
=> a purely hydrodynamic overstability of the dynamic WB scheme in a stable/marginal stratification; radiation,
conduction, the wall flux and the cache cadence play no role (so the implicit-solve hypothesis is moot). Next arms
(mode2/, agent): wb_rmin (new lower bound) at 1e12 / 6e11, wb_option isentropic / adaptive, plainsrc (problem/
wb_grav_source=plain: plain rho*g source with the WB reconstruction kept).
**FIX WORKS (mode2/, apudev arms to 1.19e6, 2026-09-13 ~04:00): `hydro/wb_rmin` (aebccba3).** WB off below 1e12: growth
x0.95 (e-fold -1.8e6 s); WB off below 6e11 (= exactly the mode region i 5-40): x0.60, DECAYS with e-fold -2.9e5 s.
Control x2.20, WB off everywhere x1.05. => set wb_rmin ~ 6e11 (just above the RCB at 1.3e11 .. the mode top at 5.7e11)
in the red-giant input; the atmosphere keeps WB. Closure (isentropic/adaptive) and plainsrc arms pending.
**Closure / source arms (mode2, 1.0e6 -> 1.19e6, e-fold at i=20):** isentropic 3.26e5, adaptive 3.43e5, plainsrc
(problem/wb_grav_source=plain: plain rho*g source, WB deviation reconstruction kept) 3.98e5 -- all still GROW (x1.4-1.5
vs control x1.83 over the window); only rmin removes it (x0.60). => the overstability is in the WB DEVIATION
RECONSTRUCTION/FLUXES (background rebuilt from the perturbed state), not the source, not the closure. RULE: the
dynamic WB scheme only where the column is meant to be hydrostatic to the closure family (the atmosphere); bound it
below with wb_rmin above any stably stratified interior. CHECK the sp hot-Jupiter runs (same scheme, deep interior).
**sp hot-Jupiter CHECKED (bench/sp_mhd_nopole/analysis, 2026-09-13): NO WB-driven mode there.** Deep KE_r grows 3.0x
(sp, WB on) vs 2.57x (cs, WB off) over rot 4-36 with identical fit rates (3.98e-7 vs 3.94e-7 s^-1 at i 1-10) = the
spin-up of the near-adiabatic interior, saturating after rot ~40 in cs; the sp/cs ratio is a flat 1.2-2.5x offset;
the stably stratified band (i >= 42, the irradiated inversion) DECAYS in both; coherence/correlations equal in both.
Note: the MHD module has NO wb_rmin/wb_rmax at all (only hydro reads them) -- add if ever needed.
**CORRECTION 2026-09-13 ~07:30: the wb_rmin "fix" is NOT confirmed.** RG_fofc_long3 (wb_rmin=6e11, from 9e5): over
1.4e6-1.94e6 rms v_r grows with e-fold 2.96e5 s at i=5 and 3.5e5 s at i=20 -- the SAME rate as RG_fofc_long
(5.4e5/3.9e5 over 1.4-1.76e6) and RG_fofc_long2 (3.1e5/3.9e5 over 1.72-2.28e6) -- and a LARGER amplitude at i=5
(2.3e4 at 1.94e6 vs 4.9e3 at 2.28e6 in long2). Character differs: C(v_r,T') at i=20 ~0 (was 0.9 with WB), at i=5 0.6-0.7.
The 1.9e5-s arms measured the transient decay of the WB-shaped pattern after the switch, then growth resumed.
=> WB shapes the mode but is not its root; the driver sits at/near the inner WALL under the plain scheme too
(i=5 fastest, buoyant; long2's wb_off arm x1.05 over 3e5 = decay + regrowth cancelling). OPEN AGAIN. Suspects: the
reflecting inner wall + gravity (wall-adjacent overstability), marginal physical convection (i 5-40 is inside the
nominal convection zone; the proxy-EOS N^2 estimate has 1e-2 uncertainty in nabla), the rmin edge residual.
