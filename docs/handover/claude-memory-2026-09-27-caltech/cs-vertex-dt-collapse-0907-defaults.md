---
name: cs-vertex-dt-collapse-0907-defaults
description: cs dhj dt COLLAPSE at a cube vertex on the 09-07 defaults, CLOSED 2026-09-11 - the EXPLICIT RT source overcools stiff cells by x/(1-e^-x), a nightside cold patch sinks, the vertex column drains and reheats, conduction reports it; semi-implicit default is right; per-cycle diag instrument built (uncommitted)
metadata:
  type: project
---

Found 2026-09-09 (viper). cs_mhd_prod2 (rot 14.4) and all four cs_hyd_rs arms (hllc 15.2,
lhllc 4.4, ausmpup 15.3, ppmx 3.7 rot) stalled with dt 1e-4..1e-9 s, no NaN, floor counters
overflowing int32, and burned GPU nodes for two days: the NaN guard cannot see a dt collapse.
`time/dt_min` (added 09-09, driver.cpp; 1e-2 in the dhj inputs) now aborts it.

The event: ONE cube vertex (its three panel copies) grows a coherent RADIAL velocity column
over ~95 of 128 radial cells, v_r 1e6 -> 1-2e7 cm/s, eint and B rising with it. The vertex
column was already the dt-binding cell in the healthy state (v+1.6cs ~4e6 at i~68). No
precursor in the floor counters (flat to within 150 cycles). Three of five collapses
happened mid-job with NO restart, so orion's WB restart bug [[wb-restart-cache-bug]] is
NOT the trigger (hllc even restarted at ncycle%10==0 and still died).

Input diff old cs_mhd_prod (ran to rot 41.6) -> the dead runs, all ADDED: wellbalance_dynamic
+ wb_x1 + wb_option=polytropic + wb_cache_every=10, isotropic_conduction=radiative (+rad_*),
problem/rot_potential=true. Conduction operator GATED on cs_test iprob 15 (09-09): energy-conservative to dump
precision, no growing vertex feature; the vertex cell is UNDER-relaxed by ~1 % at nx=16
(5x interior error, 1st order there, 2nd elsewhere), a damping deficit, NOT a heat source.
rad_cs_exact=false is 35-50x worse, so the cross term is right. cs_vertex_fill is FC-only
(a null arm for hydro). So the linear operator cannot make the column; suspects now the
WB radial path, rot_potential, or conduction's tau-blend/limiter on the real profile.
Instrument trap: kappa(T) makes an O(amp) even term; isolate with +amp/-amp parity. Orion's red giant sees the same "vertex collapse" on cs with WB + conduction.

**How to apply:** the 09-07 defaults are NOT production-safe on cs until the ablation
(nocond / nowb / norot / cache1 arms, lhllc config dies fastest at 4.4 rot) says which
piece kills it. Analysis: bench/rst_vmax.py (hydro-restart fix in the session scratchpad);
hot cells = v+1.6cs above the healthy max. See [[cs-mhd-prod2-run]], [[cs-hyd-rs-run]],
[[cs-vertex-limiter-clipping]].

**UPDATE 2026-09-09 late (viper): the collapse is NOT reproducible from a restart.** Twenty
restart arms (lhllc rot 4.0 dump, 0.42 rot before its collapse; hllc rot 15.11 dump, 0.08 rot =
1400 cycles before its collapse; new binary with rt_semi_implicit=false and the restart
rot_potential fix 42c0a6a2, and the ORIGINAL binary f6f0d6f4) ALL survived past the original
death time with dt >= 8 s. Restart faithfulness after the fix: mass/tot-E within 1e-4 of the
continued run, KE components decorrelate at the 1-25 % level within 0.1 rot (chaotic). So the
event is a fast (<0.02 rot), chaos-sensitive trigger with no committed precursor; 2 of the 5
originals (hllc, ausmpup) died within 0.1 rot of a chain restart that carried the doubled
centrifugal force [[restart-rot-potential-bug]], the other 3 (lhllc, ppmx, cs_mhd_prod2) died
inside lead jobs with no restart. Restart-based ablation is therefore USELESS here; only
from-scratch arms with long integration can discriminate, and one realization per arm has
weak power (a survival to 6 rot proves little; a death is informative). From-scratch arms on the
right physics: bench/cs_ablate/{ctl2,nocond,nowb,norot,cache1} 11539561-65; cs_ablate/ctl
11539431 is the semi-implicit-RT control. Restart trap: <outputN>/last_time must be overridden
when restarting off the output grid.

**FROM-SCRATCH ABLATION, 6-rot gate (2026-09-10): NO arm collapsed** -- ctl (semi-implicit
RT), ctl2, nocond, nowb, cache1 all reached rot 6 with min dt 9-20 s; norot at rot 5.4 clean.
ctl2 has the same physics as the original lhllc run up to ROUND-OFF and diverges from it by
14 % in KE at 0.1 rot (cs 1-ulp amplification), so the original death at 4.42 rot is a
chaos-selected event, not a deterministic property of the setup. Single-realization
survival to 6 rot has NO POWER. To discriminate the 09-07 defaults one needs either LONG
integration per arm (~20 rot; P(ctl survives) ~ e^-2.5 if the rate is ~1/8 rot) or an
ensemble of realizations per arm (GPU is deterministic, so seed via a tiny IC perturbation
or a different meshblock decomposition). Restarts are now faithful, so the 6-rot arms can be
CHAINED (bench/cs_ablate/*/rst). Decision pending with the user.

## ENSEMBLE RESULT 2026-09-10 (bench/cs_ens, problem/seed, seed_amp 1e-10, 8 rot)
The collapse is NOT a rare chaotic event: with the EXPLICIT two-stream RT source (rt_semi_implicit=false)
10/10 seeds die in the narrow window rot 4.3-5.7. Semi-implicit source: 1/8 dead (s05 at 5.3), the rest
past rot 7.3-8.0. Conduction off (explicit source): 3/5 dead. The 09-09 "all arms survived" ablation
had run semi-implicit ON by the build race, which is why it discriminated nothing. Details and the
held/cancelled job bookkeeping in [[inflight-2026-09-09-viper]]. The "conduction top suspect" line
above is superseded: the RT source is the dominant trigger, conduction a secondary amplifier.

## MECHANISM 2026-09-10 (bench/cs_ens/analysis, ctl2/s01 vs si/s01): a THERMAL RUNAWAY under the EXPLICIT RT source
- The dt owner at the abort is CONDUCTION reporting a cell at T = 1.33e5 K (50x normal), p 4.06 bar, rho/3.4: kappa_rad ~ T^3 -> dt1 0.047 s;
  hydro dt fell 7x on the same hot cell. Conduction is the MESSENGER (conduction-off still dies 3/5), not the cause.
- Cell (m,k,j,i)=(3,14,15,41), lat +32.0, lon -136.4 (nightside), 3.5 deg from the panel-0 cube vertex, one cell in from the corner.
  All six vertex cells (three corners + images) agree to a few % in T,p in BOTH arms: NO seam/panel asymmetry at the vertex.
- Timeline: dt 8.68 -> 1.13 -> 0.047 -> 0.0097 s in 86 cycles (1.4e-3 rot); three decades inside one ndiag=100 interval. Recovered
  precursor dips at rot 3.99 (2.1 s), 4.86 (1.52 s), 5.69. si/s01 never below 8.59 s in 8 rot; same median dt and cost as ctl2.
- THE PRECURSOR IS THE COLUMN ABOVE, NOT THE CELL: cold front (T<1500 K) eats down the vertex column 3e-7 bar (rot 2) -> 1e-3 bar
  (rot 5, 75/128 cells); |v_r| grows with height, all DOWNFLOW, co-located with the front. Domain-wide: 770 columns >30% cold in ctl2
  at rot 5, ZERO in si. frac(T<1500 K) 0.005 -> 0.076 in ctl2, 0.004 -> 0.001 in si. ctl2 runs 4.8x dfloor, 8.9x efloor, 5.8x tfloor/cycle.
- Code: explicit de = src*dt clipped to +-0.5 e per application (LimitRTSource, two_stream_rt.hpp ~1500) -> floor/overshoot oscillation
  in the low-density upper atmosphere, compounding; semi-implicit de = (src/lam)(1-exp(-lam dt)) cannot overshoot equilibrium. The
  clip warning fired in 14 cells (ctl2) vs 1 (si) at startup. The vertex is WHERE it lands (narrow-block/ vertex cells are the worst
  clipped), not evidently WHY. Global mass/energy/KE and the event counters show NO precursor (eos_fail, c2p_it, fofc all zero).
- NOT MEASURED: the onset (dumps 1/rot, last at rot 5.000, abort at 5.714); RT heating/tau are not in the dumps; nclip over time.
- NEXT: fine-cadence restart of ctl2/s01 from rst rot 5 (bin every 200 s, ndiag=1), ~17k cycles; restarts may not reproduce (chaos).
- Production implication: rt_semi_implicit=true (the default) is the safe choice; it changes the answer ~3% KE ([[rt-semi-implicit-changes-dhj-answer]]).
- RESTART REPRODUCES IT (2026-09-10, bench/cs_ens/ctl2_s01_fine, job 11569445, apudev 15-min one-time test): from rst rot 5.000 with
  bin every 200 s + ndiag=1, the SAME binary, it died at cycle 138854 / rot 5.7162 vs the original 138886 / 5.7144 (round-off shift).
  The earlier "20 restart arms all survived" is EXPLAINED: those ran semi-implicit by the build race. 1092 dumps (17 GB) cover the onset.
  Abort cells: (m,k,j,i)=(6,17,2,44)->(6,17,2,45), the panel-1 VERTEX CORNER cell column; T 7.8e4 -> 3.8e5 K in one cycle, F_free 2e15 -> 1e18.
- ONSET (bench/cs_ens/analysis/fine/REPORT.md, 2026-09-10 evening): NOT PHYSICS. In the last 16 cycles (33 s) the panel-1 vertex corner
  cell i=43 goes T 2367 -> 376251 K (x159) at rho x1.006; compression explains 2e-6 of it, the cells above are COLDER and moving DOWN,
  the downdraft KE is 408x too small, and the required 3.8e5 erg/cm^3/s is 1.45x the planet's luminosity into one cell. A source-term
  integration instability; radiative conduction (kappa_rad ~ T^3/rho, x77 in one cycle) is the amplifier that reports the dt.
  dt fall is 22 cycles wide: hydro owns the first two decades, conduction the last two. INVISIBLE at 200-s dump cadence (T never >2700 K
  in any dump of the abort cells). Precursor = a DRAIN: from ~400 cycles before, an accelerating downdraft evacuates i=42-47 by 30-40x
  (all-time min rho, all-time max |v_r| 1.9e6). An identical drain at rot 5.367 (dt 4.6 s) RECOVERED. No seam asymmetry between the
  three vertex images at any time; the vertex NEIGHBOURHOOD is special (min rho 50-200x lower than mid-panel, 40-60% cold burden).
- OPEN: explicit RT source vs radiative conduction as the integrator that blows up cannot be separated from hydro_w dumps. NEXT
  INSTRUMENT: restart from ~rot 5.70 (no rst exists yet; write one from the rot-5 rst) with PER-CYCLE dump of the RT source term and
  kappa_rad for meshblock 6, plus an rt_semi_implicit=true twin from the same restart.

## TWINS from rst rot 5.700 with the PER-CYCLE DIAGNOSTIC (2026-09-10 night, bench/cs_ens/ctl2_s01_twin)
- exp (explicit) reproduces the abort at cycle 138867 (fine run: 138854), but in block m=3 (k,j,i)=(16,17,42-43),
  ANOTHER image of the same cube vertex; T 8.8e4-1.3e5 K, conduction-limited. Chaos moves the image, not the vertex.
- si (semi-implicit, restarted into the ALREADY-DRAINED column) does NOT hit dt_min by tlim but its dt falls 11 -> 0.16 s
  and gid 6 column (17,2) RUNS AWAY: rt_src at i=43-45 goes -6 -> +20 -> 9e2 -> 2e4 -> 3e5 -> 2.5e6 -> 7e9 erg/cm3/s in
  6 cycles while T stays ~2.4e3 K until the last two; delta_u per cycle == rt_de (within the rk2 bracket), conduction
  divergence IDENTICALLY ZERO at onset; rt_clip pins de/e at 0.375 for 4 cycles. So: the ENERGY COMES FROM THE RT
  SOURCE, not conduction, in both integrators; the semi-implicit form only delays it once the drain exists (from
  scratch it prevents the drain: ensemble 7/8). A src growing 50x/cycle at fixed local T means the flux comes from
  another cell of the same column (a 1-D solve) or the two-stream solve is unphysical there: flux components being
  dumped next (rt_Ft/Fb/Qs/Em), bench/cs_ens/ctl2_s01_twin2/{exp_g3,exp_g6,si_g6}.
- Instrument trap: pin diag_gid to the block that actually aborts in THAT binary (the image moves between runs).

## ROOT CELL FOUND (2026-09-10 night, bench/cs_ens/ctl2_s01_twin2/exp_g3, block 3 column (16,17), flux components dumped)
Stellar term rt_Qs == 0 in the column (nightside): src is purely -(Ft-Fb)/dx. The sequence at the root cell i=52:
1. DRAIN (cycles 138640-138780, 140 cycles, dt ~10 s): accelerating downdraft v_r -5.8e5 -> -8.6e5 cm/s, rho 1e-8 -> 1.6e-10
   (60x), T 1700 -> 200 K = tfloor_kelvin; cells i=50-56 sit AT the T floor at ~1e-10 g/cm3. RT there is a weak HEATER
   (rt_de +0.03..0.5 per stage) -- the RT source does NOT drive the drain's cooling; expansion/advection does (delta_u < 0).
2. REHEAT BY HYDRO (138781-138792): T 200 -> 2433 K with rt_de ~ 0 or negative, delta_u > 0, v_r -8.6e5 -> -1.3e6.
3. EOS JUMP (138793): T 2433 -> 6622 K for e x1.43 -- the low-density dissociation plateau of the table (physical, on-table;
   eos grid logd -14..0, logT 1.5..6 brackets everything).
4. RUNAWAY (138793-138805): T -> 1.3e5 K while rt_de is a huge SINK (-600..-1730 per stage, clipped at 0.5 e) and delta_u
   stays positive: HYDRO supplies the energy at the root, RT removes it. rt_T is the EOS T of (rho,e) at the START of the
   last stage, so rt_T can exceed T(w_eint at the end of the cycle) by a lot when RT removes half the energy per stage.
5. CASCADE: the hot root radiates Fb 1e11 -> 1e14 into the thin column above (i=56-128, at dfloor 5e-14, T 1500 K);
   those cells have e/src ~ 0.4 s ~ dt, so the WHOLE floor column trips in ONE cycle (138802), clipped at de/e 0.375
   (x1.625 per cycle). The heated thin gas is carried back down by the downdraft into the root. Conduction only appears
   once T ~ 1e5 K (kappa_rad ~ T^3) and then owns the dt: dt COLLAPSE lines at (3,16,17,42-43), 45 cycles later.
ANSWER to "explicit RT or conduction in the drained cells": NEITHER at the root. The DRAIN is hydro; the reheat is hydro;
RT is a sink at the root and the clip-limited amplifier in the evacuated column above; conduction is the messenger. The
explicit source's role is UPSTREAM: it produces the cold front / cold columns (770 in ctl2 vs 0 in si at rot 5) that
sink and drain the vertex column. The si twin restarted INTO the drained state runs away too (same column type), so the
semi-implicit default prevents the drain, not the runaway. Next lever: why the vertex column drains under the explicit
source (cold front formation), or make the atmosphere robust to a drained column (floor/EOS behaviour at 1e-10 g/cm3).
CONFIRMED in the si twin (si_g6 column (17,2), root i=68): rho drops 575x (1.7e-9 -> 2.9e-12), v_r reverses +2.2e5 -> -6.9e5,
T 1690 -> 283 K, then a hydro reheat 283 -> 5950 K over 165 cycles with rt_de NEGATIVE throughout; ignition at 138855 when the
cell below reaches 2.5e4 K (rt_de flips positive, clipped); T peaks 2.4e7 K (OFF the EOS table top 1e6 K), then the column REFILLS
and T decays to 1e4 K: si survives because dt bottoms at 0.1 s > dt_min, not because the physics differs. The heating front in
exp_g3 descends from above into the T-floor slab (i=58 -> 55 over 20 cycles) with v_r < 0 everywhere and accelerating (i=53 to
-3.3e7 cm/s). The tfloor cell count FALLS during the reheat (52 -> 2); the RT clip count jumps 2 -> 72 one cycle after the root
passes 1e4 K. The run's event log never fired (log dt 3.05e4 > the run length): set output2/dt small in any short diagnostic run.
Analysis: bench/cs_ens/ctl2_s01_twin2/analysis/ (twin2_analysis.py, root2_*.py, run_out.txt, root2_task*.txt).

## CLOSED 2026-09-11 ~01:00: the EXPLICIT SOURCE OVERCOOLS STIFF CELLS BY x/(1-e^-x); the cold patch is NIGHTSIDE, not the vertex
- Ratio law VERIFIED on identical states (twins, first cycle after the shared rst, 25547 cells): de_explicit/de_semi-implicit =
  x/(1-exp(-x)) with x = lam*beta*dt, lam = 4 Em/e, beta = 0.5 (rk2 stage 2), to 0.04% median. In the healthy state 1292/32768
  cells have x > 1 (forward Euler oscillatory), 22 have x in 3-10 (upper column, T 1400-1800 K, rho 4e-13..3.5e-11, i 69-102);
  the rest at the RT/conduction handover i 31-41, T ~3000 K. Stiffness is physical (t_cool = c_v T/(4 kappa_P sigma T^4) ~ 10 s ~ dt
  at 1500 K, independent of rho); the explicit+clip integrator removes 1.6-3x too much energy per stage there.
- Ensemble dumps (bench/cs_ens/analysis/coldmap): ctl2 is colder than si in the upper layer A (i 69-104) from rot 1 (p10 -103 K,
  before any cold column exists), growing to -846 K by rot 5; cold-cell fraction 0.094 vs 0.00035; deep layer B (i 31-41) -30..-75 K,
  steady. Cold columns (>30% cells < 1500 K): ctl2 418/322/102/770 at rot 2-5, si 59/0/0/0. They form as ONE contiguous patch in the
  ANTISTELLAR panel INTERIOR (V=0, E=0 at rot 2-3; enrichment I 1.51), i.e. at a physical location, and reach the seams/vertices only
  by rot 5. The DRAINED set (rho(i=52) < 1e-9) is edge/vertex-weighted from rot 3 (E 1.94, V 1.63): the vertex is where the cold
  patch's downdraft drains deepest (the vertex neighbourhood already has 50-200x lower min rho), hence WHERE the collapse lands.
- CAUSAL CHAIN (complete): explicit source overcools the stiff nightside upper layer -> a cold patch on the antistellar panel -> cold
  gas sinks; the downdraft is strongest at the cube-vertex neighbourhood -> a column drains 60-600x to the dfloor/tfloor -> hydro
  reheats the floored gas from above (front descending, v_r < 0) -> low-density EOS dissociation plateau: T jumps 2400 -> 6600 K ->
  the hot cell radiates into the evacuated column (t_heat ~ dt), whole column trips in one cycle, clipped at 0.5 e per stage ->
  T ~ 1e5 K -> kappa_rad ~ T^3 -> conduction dt collapse. RT is a SINK at the root; conduction is the messenger.
- PRODUCTION: rt_semi_implicit=true (default) removes the cold bias (si warms, A median 3130 -> 3872 K) and is the right default.
  Residual risk: 4% of cells still have x > 1, and the semi-implicit form only linearises EMISSION; absorption-dominated heating in an
  evacuated column still clips (si twin from the drained state reached 2e7 K and recovered). A per-cell implicit (Newton on T with
  kappa(T) and B(T)) or sub-cycled source for x > 1 cells would remove the last stiffness; not built.
- COMMITTED 2026-09-11: instrument 1d988879, rad_tmax_kappa 5c58426e, pushed to fork.

## RESCUE 1 BUILT 2026-09-11: `<hydro>/rad_tmax_kappa` (K, default 0 = off) caps the T used in kappa_rad (fluxes, newdt, dt_diag;
gradient and flux limiter untouched). Uncommitted with the instrument. Inert at 0/absent/1e9 (bitwise); at 1000 K the cap can RAISE
kappa where kappa_R(T) falls faster than T^3 (T 5-6.5e3 K, table region) -- cap only at >= 3e4 K where the table is clamped.
TEST (bench/cs_ens/ctl2_s01_twin2/exp_cap, explicit, cap 3e4 K, from rst570): passed the original abort by 2870 cycles / 300 s of
simulated time, the drained column REFILLED (rho at i=52 x757, back to 73% of pre-drain), but the hot patch kept GROWING (block-3
cells > 1e4 K: 86 -> 1500, Tmax 1e6 K), dt sat at 0.01-0.1 s, and it died again at dt_min (conduction cells at 3-5e4 K, now DENSER
gas rho 4-9e-8). The cap alone does not rescue the EXPLICIT run: the two-stream clip-limited cascade continues. Not yet tested on the
semi-implicit twin (which healed on its own with dt bottoming at 0.1 s): si_cap would measure whether the cap cheapens the recovery.
si_cap TESTED (semi-implicit + cap 3e4 K, same rst): reached tlim 1.75e6, NO collapse, and is statistically identical to the
uncapped si twin over the common window (1790 vs 1783 cycles, min dt 0.18 vs 0.10 s, none < 0.1 s): the cap is a NO-OP on the
semi-implicit branch; its dt trough is hydro-owned. The si excursion is LARGE but SELF-HEALING: 18.8k of 32.8k block cells > 1e4 K
at the peak (cycle 140100, Tmax up to 6.8e6 K), back to ZERO hot cells by cycle 141500 (t 1.7458e6), block max T 6470 K again,
the vertex column refilled 2e-12 -> 1.7e-8; mass/energy drift < 1e-3. VERDICT on rescue 1: useful only as a guard for explicit
runs and insufficient there; not needed in production (semi-implicit). Keep the flag (inert), default 0.

## si/s05 POST-MORTEM 2026-09-11 (bench/cs_ens/analysis/si_s05, si_s05_pm): the ONE semi-implicit death (rot 5.34)
- Its binary predates the dt-COLLAPSE report: abort cell unknown. NO cold patch, NO cold bias (thermally == si/s01 at rot 1-5).
  It DID drain: 65 drained columns at rot 5 (si/s01 has 35!), antistellar panel INTERIOR (I enr 1.47, V=0), formed in the last
  rotation; deepest column gid 0 (J,K)=(14,11): rho x40 down at i=52, plateau 2-5e-10 over i 52-68 at 2340 K, downdraft -8.8e5.
  Global precursor in the hst (unlike ctl2): radial KE x3.9 and efloor/tfloor counters x3 in the last 0.1 rot.
- NOT REPRODUCIBLE: a faithful replay from rst rot 5.000 (new binary) tracks KE to <3% through rot 5.2, then has no excursion
  (min dt 7.2 s to rot 5.34). Chaos-selected near-miss; the explicit ctl2 death reproduced within 32 cycles.
- READING: drains form under semi-implicit too (35-65 columns/run at rot 5) but heal (twins: dt bottoms 0.1-0.18 s). s05 died at
  dt 8.9e-3 vs dt_min 1e-2: a MARGINAL case of the same self-healing excursion. PRACTICAL: production dt_min=1e-2 is too tight
  for this setting; 1e-3 lets excursions heal at the cost of ~1e3-1e4 slow cycles. Prevention lever remains option 2 (density floor
  vs the WB background). To instrument a si death: from-scratch seed with cyclediag armed over a cycle window, not a restart.
