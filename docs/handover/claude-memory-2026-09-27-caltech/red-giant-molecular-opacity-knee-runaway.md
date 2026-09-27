---
name: red-giant-molecular-opacity-knee-runaway
description: 2026-09-10 00:45 — WHAT COOLS the cells above the front in open-top runs (task_cold.py, 5 cells tracked through 30 dumps): a RADIATIVE thermal runaway on the molecular-opacity knee: kappa_R has a MINIMUM near 3200 K and rises 5-20x by 2500 K (dln(kappa T^4)/dlnT ~ -12), so emission GROWS as T falls; cells cool at constant/RISING density (adiabatic and floor artefacts rejected), park at 200-500 K where opac_tmin=2500 freezes kappa, then become the vacuum-Riemann precursor. The lid keeps the atmosphere ~1000 K hotter so no cell reaches the knee. Fix to test: opac_tmin ~3200 (flatten the knee)
metadata:
  type: project
---

Files: /orion/ptmp/jinma/Athenak/red_giant/_analysis_0910/task_cold*.py/.txt, task_noise*.
Tracked cells (T15, gid 23 (3,3) i=294 etc.): 1.32e6 T 2901 -> 1.34e6 2599 (kappa_R 4.65e-3
-> 1.56e-2) -> 1.36e6 2219 -> 1.38e6 1121 -> 1.40e6 356 K while rho RISES 7.5e-9 -> 1.3e-8 and
div v < 0 (infall, -1e-5..-6e-5/s); dlnT/dlnrho = -0.7..-3.6 (wrong sign for adiabatic,
nabla_ad 0.06-0.3); T19 30-s dumps show the pre-collapse drift is smooth (0.06%/30 s), no
floors touched (rho 1e-9..1e-8, p 10-200). Net cooling / (4 sigma kappa T^4): 0.5% at 3100 K,
15% at 2600, 32% at 2200, >100% at 1100 K. The clamp at opac_tmin=2500 does NOT stop the
collapse, it stops it being 30-90x worse (grain regime), and it is what parks the cells.
End state: cold dense low-p cell wedged between neighbours at 7x its pressure = the
vacuum-Riemann precursor ([[red-giant-vertex-floor-hydro-trigger]]). 8 cells end below the
EOS table floor (eos_logt_min = 2 -> 100 K): p(rho,e) unreliable exactly there.
Why open vs lid: at t=1.4e6 prod11's atmosphere is 3600-3800 K (0.3% of cells < 2500 K),
T15's 2600-3050 K (10.3% < 2500 K, 34 < 1000 K); T6 open is already 2524-2591 K at 1e5
with everything infalling 3-4 km/s. Noise (adjacent-cell rho ratio 2-4 vs prod11 1.02-1.04
late, 1.77 at the same age) is broad-band grid-scale scatter, NOT checkerboard/stripes
(coherent modes 2% of variance), rises monotonically OUTWARD (1.35 at i=286, 2.33 at 296,
3.98 at 306), weak x2 preference; precedes the cold cells (ratio > 1.5 at 2e5, first
collapsed cell 4.4e5 at i=296, persistent from 1.28e6). Panel degeneracy p0=p2 etc. is the
symmetric vpert seed (bit-identical in prod11 too), not an open-top symptom.
CAVEAT (mine): the cooling was inferred from 4 sigma kappa T^4, the RT source itself is not
in the dumps; a grey scheme with the same kappa for emission and absorption has a kappa-
independent equilibrium, so a net cooling of 15-100% of emission needs J << B(T) in those
cells -- why the two-stream gives a cell at tau ~0.3 above an 8000 K front so little J is
unexplained (1-D columns + corrugated front? the direct source?). The front is 1 cell; the
10<tau<100 blend layer holds 0.8-1.3 cells and 29-60% of columns have NONE (lidded too).

## TEST 2026-09-10 01:40: opac_tmin=3200 (agent a8555e50). SPLIT RESULT
opac_tmin is a LOAD-TIME clamp of the table (rows below it replaced by the first row at/above,
so 3200 -> 3350 K node); ONE ktab shared by the two-stream, the conduction (rad_kappa_src=
table_rho) and the initial column -> consistent. Clamping at the knee sets kappa_R over
2500-3350 K to the MINIMUM (less emission AND absorption). Traps: command-line overrides
cannot ADD parameters (sponge had to go into rg.athinput for a t=0 run); output/last_time
is the time of the LAST dump, not an end time.
- T21_tmin3200_r14 (restart T15 1.4e6, jobs 195213/195226): SURVIVED to tlim 2.0e6, dt 30.65
  throughout, 0 collapses; the 16 collapsed cells at 1.42e6 HEAL to 0 by 1.44e6 and stay 0
  through 28 dumps; median T(i 300-318) 2733 -> 3050 K and holds; adjacent-cell ratio 3.3
  -> 2.4. Passed 1.43113e6 and 1.51539e6 and 4.8e5 s beyond. Strongest single fix so far.
- T21_tmin3200_t0 (from t=0, 195211/195212): DIED 7.80183e5 (chain restart reproduces
  7.78037e5), EARLIER than T6 (1.04e6)/T11 (8.7e5)/T15; atmosphere NOT hotter (tracks T6,
  slightly cooler as the lower kappa predicts: 2005 K at 4e5); first collapsed cells at 6e5
  (later than T6's 4.4e5) but death is a no-precursor vertex-chimney event (top cells i
  308-318, |v2,v3| 6-8e5, identical values across panels 1/3/4/5, dt 30 -> 1e-20 in one cycle).
=> the knee runaway is what SUSTAINS the late (>1.2e6) failure; the early open-top
atmosphere cooling to ~2000-2500 K in the first 1e5 s (all open runs) and the vertex-chimney
event at 7.8e5-1.04e6 are separate and NOT cured by the clamp.
Next (launched 01:45): extend r14 to 3.5e6 (the prod12 gate) and a t0 run at opac_tmin=2800
to separate "flatten the knee" from "lower kappa below 3350"; the early atmosphere cooling
and the vertex chimney remain the from-scratch blockers.

## GATE PASSED 2026-09-10 01:30: T22_tmin3200_r14x (195239) 2.0e6 -> 3.5e6 clean
dt flat 30.65 the whole way, zero collapses, NO cold cube-vertex cell at any dump; only mild
interior cells (rel 0.17-0.30 in symmetric groups of 8) that vanish by 3.0e6; median
T(i 300-318) 2988 -> 3060-3140 K. Branch = T15 (2500 K) to 1.4e6, then opac_tmin=3200:
2.1e6 s of continuous clean evolution -> the prod12 gate condition is MET on this branch.
T22_tmin2800_t0 (195241) from scratch DIED 7.7344e5: NaN in the top 3-4 cells on panels
0/1/2/4 at once, the 7.6e5 dump already has 8 fully collapsed SEAM-edge cells (rel 9e-6,
rho 4.9e-13) at i_arr 317 on panels 1/4; atmosphere tracks the 2500 K baseline (2001 K at
2e5); noisier than 3200. => opac_tmin does not rescue from-scratch lid-free runs at any
value; the restart branch survives because of its evolved head start, not the clamp alone.
Candidate prod12 recipe: run with the 2500 K table (or lidded) to ~1.4e6, then continue
open/sponge-off with opac_tmin=3200 (both fixes from the 01:00 tasks should go in first).
