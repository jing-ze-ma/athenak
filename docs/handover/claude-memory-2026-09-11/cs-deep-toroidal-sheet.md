---
name: cs-deep-toroidal-sheet
description: CLOSED 2026-09-06 -- the kG zonal sheet at radial cells 5-6 in cs_prod_mhd_rot (and its rot-10 eruption) was the anti-diffusive stretched resistive EMF (c5c85e3b); rcmfix on 979edada has NO sheet (3-ME flat at 6-8e32, B3 at i=5/6 -6/-8 G vs +574/-614). The old run is archived; the production MHD rerun continues from rcmfix
metadata:
  type: project
---

**Observed 2026-09-05 in bench/cs_prod_mhd_rot** (cs_rot32's input + cs_full_rotation, aa6ddb9b
HIP binary): 3-ME grew exponentially from rot 2.5 (7.9e32 at 2.1 -> 8.6e34 at 6.1 rot,
e-fold ~0.6 rot), then the per-0.5-rot ratio fell 2.66 -> 1.12 by rot 6.6: SATURATING. 96 %
of it sits in the innermost 8 radial cells as a coherent zonal sheet that reverses sign
between i=5 and i=6 (mean B3 +574 / -614 G at rot 6), uniform over all six panels, with a
locked grid-scale mean-v_r pattern (+616, -797, -909 cm/s at i=5,6,7; the hydro rerun has
30-56 there). Beta there is ~1e4, so it is not the low-beta regime. The same two-cell dipole
sits at the same radius in the beta-plane run cs_both_chain at 1/30 the amplitude; the
initial field is smooth there; sp_dhj_ctl at rot 6 has a smooth 15-70 G profile. So: a
numerical deep-interior mode whose LOCATION is set by the setup and whose growth the full
rotation feeds. Cause NOT found (candidates: the eos-resistivity max_eta step, the
well-balanced source, GS07 EMF; an ablation needs restarts -- cadence raised to 1/rot in
that directory from the rot-10 slot on).

**RETRACTION (the user pushed back, rightly):** I called the run "not science-usable". The
measurement says otherwise at rot 6: KE3 tracks the hydro rerun within 10 % at every 0.5
rot; at i >= 16 rms |B| is ordinary (47 G at i=16, 27 G at i=32, tracking cs_both_chain and
sp); rms horizontal and radial velocities and the eint/rho proxy agree with the hydro rerun
to 5-10 % and 1 % from i=16 to the top. The mode is confined and saturating. Keep the run;
gate its science on (a) the saturation holding, (b) tot-E showing no extra Ohmic heating
from the kG sheet (smooth so far), (c) the atmosphere continuing to track hydro.

**Why:** [[measure-impact-before-claiming]] -- I reported a consequence ("unusable") before
measuring it. **How to apply:** for a deep-interior artefact, compare the ATMOSPHERE against
the hydro control before judging the run. See [[cs-dhj-production-retry]],
[[cs-mhd-instability-characterized]].

**STAGED, NOT SUBMITTED (2026-09-05): bench/cs_deep_ablate/{nowb,nogs07,norot,noeta}** -- from
scratch, 6 rot, the production binary (md5 97bf06ab), ONE line changed each (wb source off;
GS07 EMF off; omega = 0; constant eta = 0); apu1, 2 GPUs, 4 h slots; README.md there has the
gate. Control = cs_prod_mhd_rot itself. `cd <arm> && sbatch submit.sh`.

## ERUPTION at rot 10.3-10.5 (2026-09-05 ~23:30) -- "confined and harmless" RETRACTED AGAIN
cs_prod_mhd_rot (NO restart happened; single job since 17:51): mass -0.6 %/2 rot, tot-E -1.6 %,
KE1 x40, KE3 x5, 1-ME x13, 3-ME halved (1.3e35 -> 5.4e34), deep rms v_r 1e3 -> 1e4 cm/s at
cells 4-8, floors x2, dt still 15-18, no NaN. The saturated sheet let go. rst/dhj.00001.rst
is at rot 10.0, just before. **SUBMITTED (user): bench/cs_deep_ablate/{nowb,nogs07,norot,
noeta}_r10, jobs 11424548-51**, each `-r` that file `-i` one-change input, to rot 12, 2 h apu1;
control = the production run's own rot 10-12 (mass 3.4522e26 -> 3.4311e26, KE1 5.5e32 -> 2.0e34).
Gate: which arm does NOT erupt (mass flat, KE1 flat) by rot 11.

**RESTART ARMS READ (2026-09-06 00:20):** NONE of the four ingredients drives the eruption from
the saturated rot-10 state. nowb_r10 tracked the control DIGIT FOR DIGIT to rot 11.7, then went
NaN at rot 11.8 (dt frozen 10.4466, the grind-on signature) -- the WB source is a STABILISER
of the erupting state, not its cause; cancelled. nogs07_r10: sheet dissolves smoothly but
erupts anyway, ~0.3 rot delayed (mass 3.4395e26 at 11.69 vs control 3.4381e26; KE1 9.3e33 vs
1.24e34). noeta_r10: sheet decays FASTER without resistivity, erupts on the control's
schedule. norot_r10 confounded (frame-force removal transient). Conclusion: by rot 10 the
eruption is baked into the stored 1.3e35 erg; the question "what BUILDS the sheet" needs the
FROM-SCRATCH arms (bench/cs_deep_ablate/{nowb,nogs07,norot,noeta}, staged, not submitted).
**FROM-SCRATCH ARMS SUBMITTED (user, 2026-09-06 00:30): jobs 11425897 nowb, 11425898 nogs07,
11425899 norot, 11425900 noeta** (bench/cs_deep_ablate/<arm>, 6 rot, apu1 4 h). Gate: 3-ME vs
rot against the control's 7.9e32 (2.1) -> 1.4e34 (4.1) -> 8.6e34 (6.1), and mean B3 at radial
cells 4-8 in the rot-4 and rot-6 dumps (the +/- pair at i=5/6).
**2026-09-06 01:05: cs_prod_mhd_rot is DEAD -- dt = 0 from rot 12.19 (21:31), cycle counter
grinding for 3.5 h without advancing time; cancelled with its chained slots (11409654-56).
Final state: mass 3.4286e26 (-0.7 %), KE1 2.5e34, 1-ME 3.2e34. The from-scratch arms
(11426490/91/93: nowb, nogs07, noeta; norot NaN'd by rot 0.8 and was cancelled) are the
live question; at rot 0.6 all three have 3-ME 3.4-4.2e32 vs the control's ~3e32 -- too
early. cs_prod_hyd_rot continues (rot ~70+).**

## FROM-SCRATCH VERDICT (2026-09-06 ~03:00): the EOS RESISTIVITY builds the sheet
    rot   control 3-ME   nowb      noeta
    2.5   ~1.0e33        9.7e32
    2.8   ~1.4e33                  6.5e32
    3.3   3.1e33         3.0e33
    3.7   5.6e33 (3.6)             6.9e32   <- NO growth without the eos eta
nowb = control (WB source cleared as builder too); nogs07 lagging in time, unread; norot
NaN'd before rot 0.8 (a NON-rotating cs MHD atmosphere dies -- unexplained, item of its own).
Resistivity cannot create B; it changes the deep-layer diffusion at the eta-cap step and lets
the winding pile up the sheet. sp with the SAME eos eta has no sheet -> a cs-specific defect in
the VARIABLE-eta resistive path is the leading hypothesis (constant-eta cs_test iprob=11 is 2nd
order; variable eta was never gated; the edge-eta averaging is the same 4-cell mean on both
grids). STAGED, NOT SUBMITTED: bench/cs_deep_ablate/{etacap (constant eta = 1e13 uniform),
etalow (eos eta, max_eta 1e12)} -- README there. SUBMITTED: bench/cs_deep_ablate/ppmx (job
11430489, production input + reconstruct=ppmx, nghost=4, 6 rot) -- compare 3-ME growth and the
vertex behaviour against the control.
    4.1   1.4e34         1.45e34 (nowb == control to 3 digits)
    4.7                            7.3e32   <- noeta FLAT at 7e32 through rot 4.7: CONFIRMED

## 2026-09-06 session: nogs07 read, the CAUSE candidate found, rcmfix running
nogs07 (from scratch) tracked the control's 3-ME to rot 1.3 (6.6e32 vs 5.6e32) then went
NaN at rot 1.39 (floors exploding) -- GS07 is not the builder, and the plain-average EMF is
not viable on cs. The resistivity verdict pointed at the cs resistive path, and the audit
found it: the cs resistive dual mesh was ANTI-DIFFUSIVE on a stretched radial grid
([[cs-stretched-resistive-rcm-bug]], c5c85e3b), plus a 0.57-2.15x error in the
angular-momentum curvature source ([[cs-stretched-source-term-bug]]) and the
midpoint-vs-centroid inconsistency ([[cs-radial-unification]]). bench/cs_deep_ablate/rcmfix
(job 11442863, production input on 979edada, 6 rot) is the test; etacap/etalow are
probably moot. Every number in this file was measured with those bugs in.

## CLOSED (2026-09-06 session 2): rcmfix has NO sheet -- cause = c5c85e3b
rcmfix (job 11442863, production input on 979edada, 6 rot in 2h10): 3-ME 6.1e32 (rot 2.0),
5.8e32 (3.2), 6.1e32 (4.0), 7.0e32 (5.2), 7.8e32 (6.0) -- FLAT, exactly the noeta curve, where
the control went 7.4e32 -> 8.4e34. Mean B3 per radial cell at rot 6: control +574/-614 G at
i=5/6 (rms 850/960, max 1.8 kG); rcmfix -5.6/-8.4 G (rms 9/14, max 22/31 G). Profiles at
i>=8 agree with the control to ~20 %. So the sheet, its saturation, the rot-10 eruption and
the rot-12 death were all the stretched-grid resistive bug [[cs-stretched-resistive-rcm-bug]]
(anti-diffusive EMF near the top, wrong dual mesh everywhere), fed by winding: noeta was flat
because it removed the broken operator, not because resistivity "builds" field. etacap/etalow
were never submitted (moot). ppmx (11442987, angular ppmx + Grid-PLM radial) died AGAIN at
0.35 rot (dt 0, dfloor 3.4x rcmfix's at the same cycle) -> ppmx is not viable on cs production.
**The production MHD run is now bench/cs_prod_mhd_rot on 979edada, continued from rcmfix's
rot-6 restart (jobs 11470115 -> 11470116, 24 h each; bins every 2 rot, rst every 10 rot);
the pre-fix run is archived in old_97bf06ab/ there.** Every cs MHD memory measured before
979edada (cs-mhd-* blow-ups, eruption, 3x KE claims) ran with this bug in and needs the
re-measure caveat.
