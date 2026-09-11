---
name: red-giant-explicit-conduction-explosion
description: 2026-09-10 — the open-top red giant's lethal event is a SINGLE-CELL EXPLOSION (9e3 -> 3.6e8 K in ~2 cycles, no precursor) from the EXPLICIT radiative-conduction branch in the tau blend at the photospheric front; bisected (T18), root-cause candidate = the conduction dt uses the CELL kappa while the flux uses the FACE kappa (conduction.cpp face_flux vs NewTimeStep); fix rad_dt_face under test, blend 30/300 is the workaround (T19)
metadata:
  type: project
---

T15 (open top, sponge off) died at 1.51539e6: gid 23 interior cell (4,3), i_arr 292
(1.074 R, tau_R 32, blend w 0.51), T flat at 9078 K for 1e5 s then 3.57e8 K within ~2
cycles at constant rho. Never visible in a dump. Lidded runs never show it because their
photosphere sits ON the wall (prod11: tau=1 only in the outermost cell) -- the tau~30 front
does not exist there; the lid hides, not cures.
Bisection T18 (from T15's rst 1.4e6, restarts bit-identical across rank counts; the
CONTROL survives -> the T15 death is marginal/chaotic, dying is informative, surviving is
not): rad_tau 3/10 (fully explicit diffusion) dies 1.43113e6; rt_explicit dies 1.46406e6
(same cell i=292, hydro dt first, 3.2e10 K); rt_newton=false stalls 1.51497e6 (T15's time);
rad_tau 30/300 clean to 1.56e6. => the explicit radiative energy exchange in the blend.
Root cause (my reading of conduction.cpp): face_flux forms kappa from the FACE-averaged
T, rho, p and limits with sigma T_face^4; NewTimeStep forms kappa_ and ffree from the
CELL's own state -> across the ionization front (kappa_R varies by orders over one cell,
kappa_rad ~ T^3/kappa_R) the applied face kappa can exceed the dt's kappa by 10-100x ->
explicit step beyond its limit -> 2-cycle blow-up. Also explains the 1.0477e6 86 kK event
and the older [[red-giant-dt-collapse-solved]] history (each blend shift only moved the front).
Fix under test: hydro/rad_dt_face=true (face-consistent keff, min over faces), binary
bin_dtface, tests F0 (control must die 1.43113e6), F1 (3/10 + fix), F2 (rt_explicit + fix),
F3 (10/100 + fix: cost). Workaround running: T19_blend30 / T19_blend100 gate runs to 3e6.
UPDATE 09-09 21:30 CEST (system clock; the note above says 09-10): the previous session ended
before the fix agent wrote any code (src_vfix conduction.* still identical to the repo); the
fix + F0-F3 tests were RE-DELEGATED in session 293c860c (subagents/agent-ae59c594d8cfc1084).
T19_blend30 = job 195144 (+195145 chained), T19_blend100 = 195146 (+195147); both launched
21:20, at 1.49e6 by 21:30 with dt 30.65, ~170 s wall per 1e5 s sim -> 3e6 by ~22:15.
prod11 = 194625 at 2.87e7 (restarted from its last rst; the 3.05e7 was 194624's tail). Related: [[red-giant-vertex-chimney]].
UPDATE 21:35: T19_blend30 (30/300) EXPLODED at 1.83619e6, cell (0,7,9,300) then 299, r=1.074 R
(3.46e12), rho 1.8-2.3e-9, T 9e3 -> 40 kK -> 115 kK in 13 cycles, tau_R 117-172, w=0.64-0.86
(explicit branch), dt 30.65 -> 1e-3, crawling back (0.7 at 1.8369e6). => a deeper blend only
MOVES the front (predicted). T19_blend100 (100/1000) clean at 1.84e6, dt 30.65 (its 1.41e6
dt=5.6 dip is restart start-up only). The rad_dt_face fix is the real gate (F0-F3).
UPDATE 21:45: T19_blend100 (100/1000) ALSO exploded: 1.86455e6 at (0,4,9,300) (same r=1.074 R
as blend30's cell), T 9k -> 28 kK -> 149 kK in 16 cycles at w=0.14 -> 0.002 (i.e. the
SEMI-IMPLICIT two-stream carries >85% of the weight there, tau_R 107-174); then 1.86692e6
DEAD: (0,7,2,297) T=5e9 K, hydro dt 4.7e-18, with the two-stream "Newton step left e<=0
rescued" line at the same cycle. Cancelled 195146/195147 (dt 1e-18, unrecoverable); blend30
(195144) left crawling at dt 0.14. LESSON: the explosion also fires where w is small -> the
face/cell kappa mismatch cannot be the whole story; kappa_rad ~ T^3/kappa_R jumps by orders
when a cell crosses the H-ionization front within ONE 30 s step, and dt is evaluated on the
post-step state (dt lags). F1-F3 decide whether rad_dt_face alone is enough; if not, the
next candidate is a kappa-limited/implicit treatment of the front or the two-stream source.
UPDATE 21:50 — rad_dt_face REFUTED as the cure (bin_dtface md5 9f24268c, F0-F3, 2 nodes):
F0 control dies 1.43114e6 (= T18_blend_shallow, default path unchanged). F1 (3/10 + fix)
dies 1.43105e6 at (2,2,2,297): T 63 kK, tau 4.4, w=0.22, kface/kcell = 0.55 (the face kappa
is SMALLER than the cell's -> the dt mismatch was not operating), next cycle (2,2,9,301) at
171 kK with w=0 (pure two-stream, tau 0.004). F2 (rt_explicit + fix) dying at 1.4573e6
(dt 0.26), same as T18_rtexpl. F3 (10/100 + fix): dt ratio vs T18_ref exactly 1 -> the
fix costs nothing but also does nothing here. CONCLUSION: the runaway is NOT an explicit
conduction stability overshoot; it fires in the two-stream-dominated (w~0) and even
optically-thin cells, with Newton on or off and rt_explicit. Every radiative energy
exchange variant heats the same photospheric-front cells to 1e5-1e9 K in a few cycles.
Next suspects: a positive-feedback absorption (kappa rising with T at 1e4 K in the
opacity table, cell absorbs the interior flux faster than it re-emits) in the two-stream
source; check J vs B and the per-band heating rate in the exploding cell from a 1-cycle
dump before 1.43105e6 (F1 has rst?), and try capping kappa_P/kappa_R in the source at the
front (rt_srclim) or a sub-cycled/implicit source. The 30/300 and 100/1000 blends and the
lid only move/hide it.
FINAL (agent report 22:05): jobs F0 195151 / F1 195152 / F2 195153 / F3 195154, all ended.
F2 stalled at 1.457327e6 (dt slid to 1.8e-2, no collapse event). F3 (10/100 + fix) DIED at
1.540578e6 (hydro dt 1e-26 at (0,6,6,294), rho 1e-18, T 3e10: an exploded floor cell) where
T18_ref survived to 1.56e6 -> chaotic single-cell event, also unprotected by the fix.
F0 is NOT bit-identical to T18_blend_shallow (7th-digit dt drift from cycle 60100; the
face_flux -> RadFaceKappaCgs refactor changed inlining/FMA); behaviourally identical.
Fix lives ONLY in src_vfix (conduction.cpp/.hpp + mesh.cpp collapse printout kface/kcell);
keep as a diagnostic, not a cure. The agent's closing suggestion to use T19 blend30 as the
live path is WRONG (blend30 exploded at 1.836e6, see above).

## RESULT 2026-09-10 ~13:30 — rad_dt_face built and tested: diagnosis CONFIRMED, fix INSUFFICIENT
bin_dtface (src_vfix; conduction.cpp: shared RadFaceKappaCgs helper, face-consistent keff
min over 6 faces, diag slot 16 kface/kcell; default off = bit-identical: F0 control dies at
1.43114e6 like the baseline). At the F1 exploding cell kface/kcell = 143 and w=0,0 with
dt1..3 = -0.3: the old dt loop RETURNED EARLY (wmax==0 && blend_r) while the ANGULAR faces
(4-point-averaged weight) still carried heat -> no constraint at all there. With the honest
limit the runs still die, now as dt collapses (F1 1.43105e6, F2 1.4573e6): the flux-limited
explicit diffusion is GENUINELY STIFF at the ionization front; F3 (10/100 + fix) even died
at 1.5406e6 where the reference survived. Zero cost where not stiff (dt identical to 7
digits for 3600 cycles). Verdict: keep rad_dt_face as a diagnostic (turns silent 1e8 K
blow-ups into attributable collapses), NOT for production; the cure is to hand the front
to the semi-implicit two-stream = deeper blend. T19_blend30 (195144->195145) and
T19_blend100 (195146->195147) both clean past 1.75e6 with dt 30.65 (blend100 had one
post-restart dip to 5.6 s at 1.4137e6, recovered). Gate = 3e6.
Anomaly kept on file: T15 (continuous, 96 ranks) died at 1.51539e6, T18_ref (restart at
1.4e6, 32 and 96 ranks bit-identical) survived: a restart is not bit-identical to the
continuous run (WB cache rebuild / budget accumulators) and the event is marginal.
NEXT (launched 22:15, session 293c860c, subagent ae9c52e5b4c0c5944): catch-it-in-the-act.
D0_stageA = T18_blend_shallow config 1.4e6 -> 1.4295e6 with an rst at the end; D1_catch =
restart with hydro_w every 30 s + rt_apply_debug on the column of (0,9,9,295)@rank5 to
1.432e6; D2_catch = same with the F1 config. Analysis: per-cycle energy budget of the
exploding cell (two-stream source vs conduction divergence vs hydro), opacity slope
dkappa/dT across 8-15 kK at its rho. Scripts go to red_giant/_analysis_0910/.

## GATE FAILED 2026-09-10 evening — the blend ladder SATURATES; prod12 NOT launched
T19_blend30 died 1.83619e6 (deterministic: the chained restart from 1.8009e6 reproduced the
same cycle 73052, same rank), T19_blend100 died 1.86455e6. Ladder: 10/100 -> 1.515e6,
30/300 -> 1.836e6, 100/1000 -> 1.865e6 (+1.6%): asymptote ~1.85e6. Killer both times: a
conduction-branch thermal runaway in the TOP ~6 cells (r 3.451-3.461e12 = 0.980-0.983 x1max),
T 3e4 -> 6e5 -> 5e9 K in 10-20 cycles, different angular cells per run = a shell-wide
instability of the open top, not one bad cell. Those top cells carry tau_R 110-174: the
domain top at 1.1 R sits INSIDE the opaque envelope (the initial column's tau=2/3 is at
r=3.86e12 = 1.2 R, outside x1max=3.52e12 -- the startup line says so), and with the open
top the envelope keeps expanding/piling into the boundary (top-20-cell mass x10 by
1.24e6), so the stiff explicit front always ends up in whatever blend window is set.
=> The open outer boundary at 1.1 R is ILL-POSED for this star, not a bug to patch.
Options: (1) domain to ~1.25-1.3 R so the photosphere is inside, radial grid refit keeping
>= 5 cells/H_p in the star ([[red-giant-15R-grid-starved-deep]] lesson), open top,
sponge off; (2) stay lidded = prod11 (healthy at 3.0e7, zero collapses); (3) a
semi-implicit radiative-diffusion operator (code project). Recommend (2) now, (1) if
lid-free physics is required. Reproducer for any future attempt: T19_blend30/rst/rg.00010.rst
(1.8009e6) -> collapse at 1.83619e6 in ~14 min on 6 nodes.
UPDATE 22:40 (session 293c860c): user asked about two-stream EVERYWHERE (no conduction blend),
also explicit. Fact: ck_pcut_bar=1e30 already -> the sweep covers the whole column; the
blend only sets which operator DEPOSITS (w=1 conduction, w=0 two-stream; w rises from lo
to hi). So two-stream-everywhere = rad_tau_lo=1e29, rad_tau_hi=1e30 (input only, no cost
change). Launched on the deterministic reproducer (T19_blend30/rst/rg.00010.rst = 1.8009e6,
30/300 dies 1.83619e6 in ~14 min on 6 nodes): T20_ts_semi (195177, semi-implicit source)
and T20_ts_expl (195178, + rt_explicit=true), tlim 1.95e6. Expectation: expl dies (local
radiative relaxation time ~20 s at the front AND deep vs dt 30; T18_rtexpl/F2 already died);
semi is the real question -- blend100 exploded at w=0.14 and F1's 171 kK cell had w=0, so
the two-stream itself is suspect. Implicit-conduction assessment given to the user: radial
implicit per-column tridiagonal = 2-3 days (columns are rank-local: mb nx1 = 320 = full
radius; angular dt is 300-2000 s so stays explicit); full 3D implicit = weeks; RKG STS not
worth it. WAIT for the D1_catch budget before building anything.
T20 RESULT 23:05: T20_ts_semi (two-stream everywhere, semi-implicit) died 1.83707e6 (hydro dt
1e-9 -> frozen 1e-71 at 1.83853e6; 3 Newton e<=0 rescues; cond dt 9e306 = blend provably off);
cancelled 195177. T20_ts_expl (rt_explicit) SAILED THROUGH 1.836e6 at dt 30.65, front layer
i 292-301 NOT running away (Tmax 4-9.5 kK, falling), then died 1.8636e6 by a DIFFERENT route:
NaN at the OUTERMOST active cell (0,7,7,321) gid 35 (i=320 rho 1.6e7, v1 -6e20; open-ghost
guard 15992 fallbacks in the last report; L_rad,out/L 2e27) = an open-top boundary blow-up,
not the in-layer thermal runaway. Cost: semi +5.5%, expl -1.4% vs blend30 per sim second.
So: the blend/semi-implicit deposit IS implicated in the 1.836e6 death; the explicit
two-stream removes it and exposes the open-boundary cell as the next failure.
CORRECTION to the "photosphere outside the domain" note: the collapse printout's two taus
are the cell's LOWER and UPPER face tau_R (diag slots 9/10 = rad_tauf(i), (i+1)); at the
death cells tau = 174 below / 0.16 above -> the whole front sits in ONE cell at i~300, 20
cells below the top, r = 1.077 R INSIDE the domain; the startup "tau=2/3 at 3.86e12, T=nan"
line is the failed IC search. Agent ad2879c5 measuring the tau profile from dumps.
TAU PROFILE MEASURED 23:15 (_analysis_0910/taskF_photosphere.py, kappa_R replica verified to
6 digits vs the printouts): photosphere tau_R=2/3 at r = 3.43-3.46e12 = 1.07-1.08 R, i_arr
289 (IC) -> 298-300 (open runs at death) -> 315 (prod11 lidded, 6 cells under the wall);
20-30 cells BELOW x1max, top cell tau 1e-5..1e-4. The "outside the domain" claim is WRONG.
The H-ionization front is ONE CELL in every column (tau 2 -> 250 across a single face,
kappa_R 0.01 -> 15, T 3300-5000 -> 9200 K, rho DROPS inward 2-3x = density inversion);
H_p/dr 13-16 below, 2.5-7 above (the grid is fine in H_p; the jump is narrower than H_p).
tau=2/3 sits 3-6 cells ABOVE the kappa jump, in the neutral gas. Mass above tau=2/3
2.6e28 g (open) / 1.4e28 g (prod11). prod11 has the same 1-cell front and is healthy, so
a 1-cell front alone is not fatal; the difference is the open top's thin atmosphere above.
