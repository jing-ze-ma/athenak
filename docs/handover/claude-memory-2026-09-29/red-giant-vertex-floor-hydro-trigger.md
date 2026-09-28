---
name: red-giant-vertex-floor-hydro-trigger
description: 2026-09-09 23:30 — CAUGHT IN THE ACT (D1/D2_catch): the open-top "photospheric thermal runaway" is TRIGGERED BY HYDRO at a permanently FLOORED cell (T=0.1 K, i_arr 301, r=1.082 R) that exists in the 8 cube-vertex columns and ONLY there; the vacuum-Riemann update creates +1.7e39 erg in one cycle (995 km/s jet, 29x v_esc), the opacity's dln k/dln T=+10 across 8-11 kK flips the tau blend w 0->1 in one stage, and the explicit conduction then delivers x116. Radiation is the amplifier, not the cause. Fix the VERTEX COLUMN, not the blend/source
metadata:
  type: project
---

Runs (2 nodes, from T15_prod12gate/rst/rg.00008.rst = 1.4e6, no intermediate restart):
D1_catch 195176 (T18_blend_shallow config 3/10, bin_catch md5 6c9a8467, = F0 to every digit:
collapse cycle 60813 t=1431135.185, cell (0,9,9,295) rank 5 = gid 15 panel 0 CUBE-VERTEX
column, r=3.4449e12=1.0765 R) and D2_catch 195184 (+rad_dt_face, = F1: (2,2,2,297), the
OTHER vertex column mb 32). Per-RK-stage budget from the rt_apply column stream (new
problem/rt_dump_rank + a host header line, outside the kernel; bit-reproduction verified):
- Before: cell i_arr 301 of the vertex column sits at the ENERGY FLOOR (T 0.1 K, e/rho 2.8e8
  vs 1.7e12 normal) for the whole 3.1e4 s; census: exactly 8 such cells in 6144 columns =
  the 8 cube vertices (mb 0/5/10/15, 32/37/42/47 at |x2v|=|x3v|=0.9688). Standing radial
  downflow 3-5 km/s into it. T at the victim i=295 is 2802.2 K flat to 0.01 K for 2000 cycles.
- cycle 60810 s1: the floored cell's last energy wiped; s2: cell 300 gets +3925 erg/cc (x24,
  1988 -> 9918 K) with dE_RT = +0.113 (0.003%) and conduction EXACTLY 0 (w=0 both faces) ->
  HYDRO. Column total energy (i 290-310) 1.55e39 -> 3.27e39 erg in ONE cycle (+112%, 92%
  kinetic), the three angular neighbour columns unchanged to <0.01%; v_r at i=301 -3e5 ->
  -1e8 cm/s (29x v_esc, 200x c_s): energy CREATED at a vacuum interface (floors).
- Amplifier: kappa_R(rho 4.4e-10) rises 4150x from 2 kK to 10 kK (dln k/dln T +9.6..+11
  over 8-11 kK, turns over above 11 kK), so that one cell's dtau goes 3e-3 -> 14 and the
  blend weight of the column below flips w=0 -> 1 in ONE stage (two-stream off, explicit
  conduction on at a dt computed with w=0 = unconstrained); next stage tau drops, w -> 0,
  the two-stream flux comes back sign-reversed 6400x (+9.6e9 -> -6.2e13) and gives the
  victim +59% then +1035% (D2: 99.95% two-stream at w<=0.235 -- the amplification does not
  need the blend to flip). The blend is BISTABLE across the ionisation front.
- Kill: victim at 7.2e4 K, tau 9, w=0.984: explicit conduction x116 in one 7.4 s stage ->
  1.77e7 K (recomputed with condcheck.py to 0.3%). Newton "e<=0" fires AFTER the death.
Why every radiative variant died the same way: all downstream of a hydro event none controls.
NEXT: (1) confirm the other deaths (T15 1.515e6, T19 1.836/1.865e6, T20_expl top cell) are
vertex-column floored-cell events (census of floored cells vs t); (2) fix the vertex column
under the open top (why it drains to the floor: [[red-giant-vertex-chimney]] suspects =
cube-vertex 2x truncation, flux_seam_cc two-seam-face cell, corner fill); interim: a
guard/damping on the 8-24 vertex columns above the photosphere, or an energy-conserving
floor. Scripts: _analysis_0910/catch.py budget.py condcheck.py opac.py vertexbudget.py
task_catch.py -> task_catch_out.txt; runs D1_catch/D2_catch (5.5 GB each, 15 s dumps).
Related: [[red-giant-explicit-conduction-explosion]] (the refuted dt story), [[cs-blast-vs-cartesian]].

## CENSUS 2026-09-10 00:15 (task_floor.py): the vertex is NOT the general trigger
Deaths classified geometrically: D1 (1.431e6) and T15's 1.0477e6 stall = VERTEX columns;
T15's death (1.515e6, gid 23) = INTERIOR, 11 cells from a vertex, its column carried a
collapsed cell (e/rho 0.03x median, i=294) for 1.35e5 s; T19_blend30/100 = x2-SEAM edge
cells 3-5 from a vertex, precursor 3-5e4 s; T20_ts_expl = top cell i=321 of a HEALTHY
interior column (no precursor, pure open-boundary hydro NaN). Vertex cells are clean (0/24)
in every dump after 1.44e6; 8-48 cold-collapsed cells persist at i=294-301 (r 1.077-1.082 R,
the layer just above the 1-cell front) in all open runs; lidded prod11/E3 have none.
=> the general precursor = COLD-COLLAPSED CELLS in the thin layer right above the front
under an open top, at any angular position; the vertex is only the first place.
Chart asymmetry: the 3 panel copies of each vertex (adjacent cells 2.3 deg apart) differ
x3 in rho, grouped by chart axis, identical across all 8 vertices to <1%; panel means pair
(p0=p2, p1=p4, p3=p5 to 5 digits) = the top-of-atmosphere flow is GRID-LOCKED, not
turbulent; IC exactly symmetric. Control: the median ADJACENT-cell rho ratio within a block
at that radius is already 2.9-3.6 in open runs (seams 1.1-1.3; prod11 1.02) -> the open-top
atmosphere above the front is grid-scale noisy EVERYWHERE; that is the disease.
Files: _analysis_0910/task_floor*.txt. NEXT: thermal history of a collapsing interior cell
(adiabatic expansion vs radiative cooling; is its kappa at the opac_tmin=2500 K cliff?),
and why the layer is grid-scale noisy (WB/wb_rmax above the front? ghost fill? RT?).
