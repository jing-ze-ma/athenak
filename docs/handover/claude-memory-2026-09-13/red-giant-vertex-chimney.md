---
name: red-giant-vertex-chimney
description: 2026-09-09 22:30 — the open-top vertex collapse DISSECTED: the 24 cube-vertex columns are evacuated CHIMNEYS (rho 70x low at 8.2e5, 1e4-1e5x by 8.54e5, inflowing while the bulk outflows, 3x shear below the sponge); the sponge (exact) + open top turns it lethal; wall caps it; sponge-off survives. Ghost-fill basis, corner ghosts, conduction cross term, vertex fill flag, explicit RT all CLEARED
metadata:
  type: project
---

Runs (all restarts of T11_spongefix/rst/rg.00004.rst = 8e5, binary 887c241e, under
/orion/ptmp/jinma/Athenak/red_giant/): T12_spongeoff SURVIVES 9e5 (excursion peaks 8.4e5,
relaxes); T13_vfilloff null (FillPanelCornersCC unconditional for CC); poison test:
corner ghosts read only via conduction cross term at 2e-6; T14_ghostfix (covariant ghost
fills, ported to the main tree) dies identically 8.721147e5; E2 rad_cs_exact=false
identical; E1b rt_explicit dies 8.718e5; E4 sponge only i>=316 dies 8.766e5 with the cold
cells inside the shallow sponge; E3 outer_bc=wall SURVIVES, vertex |v_h| capped 4.8e5.
E1 (rt_grey=rt_ck=false) is a NULL: the Eddington relaxation fallback stalls dt at once.

Anatomy (E0_fine, 2e3-s dumps, m=0 vertex cell = (k,j)=(7,0) in the bin dump — the
first angular axis is REVERSED vs mb_geometry; (0,0) is a healthy cell!): at 8.2e5 the
vertex column at i=310 has rho 1.3e-12 vs bulk 8.9e-11, v1 = -4.4e5 (INFLOW, bulk +2.5e5).
Order: |v_h| doubles every 2e3 s from 8.24e5 (2e4 -> 1.8e6) -> e/rho sags from 8.28e5 ->
rho min + deepest inflow 8.40e5 -> e/rho crashes to the floor 8.54e5 at i=310-311 (12
cells above zs, 11 below the top) -> flow reverses to outflow, cold front marches up ->
hydro dt collapse 8.7211e5. Profile at 8.54e5: vertex column 1e4-1e5x under-dense above
i~295, |v_h| 4e5 through the whole upper envelope incl. i=280 (bulk 1.6e5) -- the shear
and the evacuation are NOT made by the sponge. Sponge guard fires on that exact cell
1.5e3 s before the floor (symptom). No NaN anywhere: pure dt collapse.

Decision: the sponge exists for REFLECTING walls (its own comment); under an open top it is
pointless and lethal -> prod12 gate = open top + sponge OFF to 3e6 (T15). Open question
being measured: does the chimney form without the sponge / under the wall, and when did it
start (rho_vertex/rho_bulk vs t through the T6..T11 chain)? If it forms anyway, the open
top has a vertex defect (mass leaves the vertex column) that needs its own fix.
Related: [[red-giant-sponge-cs-kinetic-energy-bug]], [[cs-orthogonal-ke-audit-2026-09-09]].

## UPDATE 2026-09-10 morning — T15 (sponge off) STALLED at 1.0455e6, not cured
T15_prod12gate (195093): dt 30.65 -> 0.03 at t=1.0455e6, crawling back (1 s at 1.052e6);
open-ghost guard fired at ~9.2e5 in 64 columns with d_i = 1e-18 (dfloor) and T at the
floor in the TOP active cell. Same 1.04e6 event as T7 (sponge on) and T10: the sponge was
never the primary. NEW PRIME SUSPECT: the open-ghost guard's FALLBACK (red_giant.cpp
column_state ~2642, called at ~2684/2722) puts the dense INITIAL hydrostatic column into
the ghost above a column that has evacuated to the floor = a pressure wall that drives
INFLOW (the vertex columns inflow at -4.4e5 at i=310 while the bulk outflows). Fix under
test: problem/open_guard_fallback = zerograd (d_g=d_i, e_g=e_i when ghost_ok rejects),
binary bin_zg, T16_zg_nospg (from T12's 9e5, sponge off, tlim 1.3e6: does the 1.0455e6
collapse vanish?) and T17_zg_spg (from T11's 8e5, sponge on, tlim 1e6: does the 8.72e5
vertex death vanish?). Also being measured: H/dr and rho at the top vs t (is the open top
a vacuum boundary for the thin atmosphere?), the mass flux through the top, the chimney
history along T6..T11 and the first guard firing in each (which comes first).
Agents' transcripts: .../884370bd-7441-4eee-9239-e5fc61707ea3/subagents/agent-a233e12b6617aac84.jsonl
(analysis) and agent-ae6cf1324b4b472e9.jsonl (zerograd build + T16/T17).

## UPDATE 2026-09-10 ~10:30 — T15 RECOVERED; chimney anatomy measured; pressure-wall REFUTED
T15 (sponge off): dt back to 30.65 by ~1.06e6, at t=1.445e6 healthy -- past T7/T10's 1.04e6
and R2a's 1.44e6 deaths. The 1.0477e6 collapse was RADIATIVE CONDUCTION's explicit branch
(dt1) at an 86 kK, rho 3e-9 cell at i_arr 295 in the tau 26-76 blend, on a panel-seam
block corner that is NOT the cube vertex; hydro dt untouched (31.9). Same class as
[[red-giant-dt-collapse-solved]] (hot cell in the explicit branch).
Guard/pressure-wall hypothesis DEAD: every guard line names i=322 with ANGULAR ghost
indices (0,1,10,11) = the FillPanelCornersCC corner-halo columns, whose "top active cell"
is itself an extrapolated ghost at dfloor. No active cell ever has rho<1e-16 (6144 top
columns, every dump). H/dr >= 2.5 everywhere at the top: the open top is NOT a vacuum
boundary. Mass loss through the top 9e5->1.04e6 = 7.8% of the top-20-cell mass, refilled
from below faster (top-20 mass GREW 1.5e27 -> 1.3e28 g by 1.24e6).
Chimney anatomy (scripts /orion/ptmp/jinma/Athenak/red_giant/_analysis_0910/): born DEEP
at i_arr 290 (1.074 R ~ the photosphere) from the first dump (R_rho 0.96 at 2e4, 0.60 at
1.2e5, 0.36 at 2e5), climbs to i_arr 310 by 2.9e5; CYCLIC, period ~3e5 (m=0 R_rho at 310:
0.044 -> 0.0055 (8.4e5) -> 0.13 (9e5) -> 0.0017 (1.06e6) -> 0.098 (1.2e6) -> 0.003
(1.36e6)); at 1.04e6 the vertex column is a coherent -7.5 km/s DOWNFLOW funnel, 200x
under-dense, T 2000 K vs bulk +4..7 km/s outflow, 3000-9000 K. OPEN-TOP ONLY: E3_wall
heals it (0.0145 -> 0.065 in 6e4 s), prod11 has NO chimney (R_rho 0.92-1.01 at 1e7-2e7).
Sponge off does not remove it (deepest values are sponge-off). Cleared: corner ghosts
(poison 2e-6), seam-end resample, conduction cross term, ghost-fill basis, guard fallback.
Remaining suspects for WHY the vertex column drains under an open top: the cube vertex's
2x truncation error ([[cs-blast-vs-cartesian]]) in a marginally resolved (H/dr 2.5-4)
flowing atmosphere; flux_seam_cc at the only cell with two seam faces. Impact if left:
24 of 6144 photospheric columns (0.4%) are funnels; the conduction dt collapses cost ~2x
wall during each chimney minimum. T16/T17 (zerograd fallback) now only confirm the refutation.

## UPDATE 2026-09-10 ~11:15 — T15 DEAD at 1.5154e6: PHOTOSPHERIC THERMAL RUNAWAY (new, lethal)
T16_zg_nospg / T17_zg_spg (zerograd guard fallback, bin_zg): T17 dies at T11's exact
8.721147e5; T16 dies HARD at 1.0413e6 (T15 recovered from the same event by luck).
Guard fallback CONFIRMED irrelevant (tracks baseline to 0.6%). The edit stays in
src_vfix only (option open_guard_fallback, default column).
T15 (195093, cancelled 11:10, rst every 2e5 up to rg.00008 = 1.4e6?) 2nd collapse at
1.51539e6: conduction dt 1.9e-13 at code (0,5,5,292) = the BULK (3,3) column, r=3.4355e12
(1.074 R, the photosphere), rho 7.1e-9, T = 3.57e8 K, tau_R 32; 70 cycles later (0,5,5,296)
at 1.36e8 K. The 1.0477e6 event (86 kK at (0,2,9,295)) was the same thing, milder. So the
open-top star has TWO independent defects, both absent under the wall (E3, prod11 to 3e7):
the vertex chimney (non-lethal without the sponge) and a photospheric thermal runaway in
bulk columns (lethal). Agent tracing the runaway's growth + a hot-cell census vs t with
wall controls: transcript .../884370bd-7441-4eee-9239-e5fc61707ea3/subagents/agent-a79b55bdf9846202e.jsonl
RECOMMENDATION forming: prod11 (lidded) IS the production; the open top needs the runaway
understood (RT two-stream or conduction blend at tau~30 with an outflow?) before prod12.

## UPDATE 11:50 — the runaway is a SINGLE-CELL EXPLOSION with no history (agent taskC/D/E)
gid 23 (panel 1, interior cell (4,3), 4 cells from any block face), i_arr 292: T 9025 -> 9078 K
flat over 1.40-1.50e6 (86th percentile of its shell), then 9078 -> 3.57e8 K within ~2
cycles at t~1.5147e6 at CONSTANT rho (8.3e-9 -> 7.1e-9): tau_e <~ 65 s. Whole layer
i_arr 280-321 has ZERO cells > 2e4 K in every dump of T15, T11, T12, E3, prod11 (to 3e7),
so it never shows in dumps. Sits at the tau_R~30 base of the photospheric front (open-top
photosphere at 1.077 R, <tau>=1 at i_arr 299; prod11's photosphere is ON the wall, tau=1
only in the outermost cell -> lidded runs have no such front: the confound). The 1.0477e6
86 kK event was at gid 37's cube-vertex cell and HEALED within one dump. Diagnosis: only an
implicit/semi-implicit update can deposit 4e4x the thermal energy in 2 steps -> the
two-stream's Newton/closed-form source solve at blend weight w=0.51 is the suspect
(the "Newton step left e<=0 rescued" messages are the same solver). Bisection T18_{ref,
rtexpl,nonewton,blend(30/300),blend_shallow(3/10)} from T15's 1.4e6 rst, agent transcript
.../884370bd-7441-4eee-9239-e5fc61707ea3/subagents/agent-a32ed6085be03b591.jsonl
