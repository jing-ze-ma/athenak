---
name: red-giant-prod11-died-3e7
description: 2026-09-10 00:40 — the LIDDED production prod11 DIED DETERMINISTICALLY at t=3.210573e7 (cycle 1061300): jobs 194624 AND 194625 both restarted from rst/rg.00006.rst (3.0e7) and hit the same FATAL non-finite in 26574 cells at i_arr 291-310 on many ranks at once (shell-wide, no dt precursor, nan_check=100); chain is dead, "prod11 never dies" is FALSE; prod11/out.txt is 42 GB of UCX mpool.c:55 WARN spam. Reproducer prod11_catch launched (agent abcfdfcf); do not resubmit blindly
metadata:
  type: project
---

Facts: sacct: 194624 FAILED after 3:03:33 (19:58), 194625 FAILED after 3:02:36 (23:01), both
12 h limit -> not wall; both segments' out.txt end identically: last cycle line
elapsed=1.09e4 cycle=1061200 time=3.210271e7 dt=30.14, then "### FATAL ERROR driver.cpp:533
non-finite conserved density or energy in 26574 cell(s) at cycle 1061300, time=3.210573e7",
first bad cells (0,9,11,298) rank 35, (0,11,11,305) rank 30, (0,5,4,298) rank 74, (0,1,11,298)
ranks 41/42 (spans i 293-312 = the photospheric layer, prod11's front at i_arr 313-315).
194625 REDID 194624's segment because no rst had been written after 3.0e7 (rst dt 5e6).
Last dump #160 = 3.20e7, 1e5 s before the death. MaxRSS 591 MB -> no OOM.
Log: prod11/out.txt is 42 GB; the middle is millions of "[orionNNN:pid:0] mpool.c:55 UCX
WARN object ... tag_send from host memory length 102400 rendezvous zero-copy read from
remote" lines (an earlier segment; the last 20 MB are clean). NEVER grep it whole; use
tail -c. Consider UCX_LOG_LEVEL=error or a different transport in future sub.sh.
Reproducer: prod11_catch/ (from rg.00006 with nan_check 1, rst 2.5e5, bin 2e4) ~3 h to the
event; README in the directory. Analysis of dumps #150-160 for precursors: agent abcfdfcf
(_analysis_0910/). Related: [[red-giant-vertex-floor-hydro-trigger]] (same layer in the open
runs), [[red-giant-molecular-opacity-knee-runaway]].

## ANALYSIS 2026-09-10 01:45 (agent abcfdfcf; _analysis_0910/task_G*.py -> task_G*_out.txt)
CORRECTION: rst/rg.00006.rst is at t = 2.585e7 (its output3/last_time=3e7 is the exit bump;
output2 last_time 2.58e7 / file_number 130 is the truth), so each dead segment was 6.26e6 s.
Top shell i_arr 290-321 at dumps #100-#160: ZERO cold/hot/over-dense/fast cells at every dump,
Tmax 10.7 kK flat, vmax 6-9e5, no sponge-guard/clip line, hst mass/E constant to 7 digits
573 s before the death; front tau=1 at i_arr 313.1 stationary over 2e6 s (0.12 cells/1e6 s
long-term; NOT running into the wall). => NO precursor visible; the NaN is born inside the
last 3020 s (100-cycle check) after dump #160. Circumstantial: the NaN shell (freshest spans
i_arr 305-307) is INSIDE the sponge (first damped cell 298) at Mach 0.85-0.97, T 2650-3150 K
= right at the opac_tmin=2500 clamp; the sponge guard covers only ei<=0.
Reproducer prod11_catch job 195261 (6 nodes, 6 h wall, NOT chained): from rg.00006 with
nan_check 1, rst every 2.5e5 (rst/ from file 100), hydro_w every 2e4 from 3.19e7 (files
from 01000); death expected ~03:35 CEST 09-10. README.md in the directory.
prod11 usable data = 2.0e7-3.2e7 (dump #160); DO NOT resume blind (bit-identical restarts,
any rank count). Next: read 195261's FATAL block (true first cycle + seed cells), bisect
from the nearest 2.5e5 rst with sponge=0 and opac_tmin as the first switches; future
production: nan_check <= 10 and rst cadence << 5e6.

## REPRODUCER RESULT 2026-09-10 03:40 (prod11_catch 195261, nan_check_cycles=1; verified)
FATAL at cycle 1061295, t=3.210557e7: 6 cells on 2 ranks, ONE seed column (rank 59, panel 3,
J=14, K=29, code i=305-307, r 3.478e12, 14 cells under the lid) + its ghost copy on rank 62.
prod11's "26574 cells shell-wide" was the 100-cycle blind spot: 6 -> 26574 cells in 5 cycles.
Not a vertex, not a seam, NOTHING floored (min e/rho 9e11, floor mark 1.42e7, all 11 dumps),
no dt precursor (dt flat 30.14). The column: 1-cell ionization front 8983 K -> 2884 K at
305/306, cold side 4-6x overdense and sinking, cooling monotonically through the 3200-2500 K
molecular-opacity knee over the last 2e5 s = [[red-giant-molecular-opacity-knee-runaway]] in
the LIDDED run. ~25-30% of the shell at i_arr 308 sits on a cold plateau (T~3300 K, top of
the knee); this column is the only one that fell through it. Last dump 01010 at 3.210e7,
nearest rst rg.00124.rst at 3.205002e7 (184 cycles before). Analysis: _analysis_0910/
task_H_out.txt, task_H2_out.txt, task_H_catch.py, task_H2_col.py, task_G4_out.txt (its
"i_arr" column is the dump index; true i_arr = code i).
Bisection launched ~03:45: B1_catch_ctrl (determinism) and B2_catch_tmin3200 from rg.00124.

## BISECTION 2026-09-10 04:15 (B1/B1b/B2 from rg.00124.rst, agent-run, numbers as reported)
B1 ctrl (pinned prod11 binary): dies at the EXACT cycle 1061295 / t 3.210557e7 / cell
(0,7,8,305) rank 59 -> deterministic. The pinned prod11 binary has NO opac_tmin (predates
it), so B2 used bin_nandiag/athena; B1b = that binary, param off, dies 7 cycles later at the
same span (binary swap innocent). B2 opac_tmin=3200: NO NaN to tlim 3.26e7 (4.9e5 s past the
death), dt flat 30.14, 0 floor events. BUT the cold column becomes a coherent 20-cell
DOWNDRAFT (v_r -2.2e5, e/rho 0.45x median, rho drained 2x), cold cells 2 -> 102 (0 vertex),
min T top-45 2222 K (< 2527 K in prod11_catch), L_out/L 0.76 -> 1.30. Shell median unchanged
(plume-local). OPEN: does the clamp remove the runaway or just let the cell convect instead
of NaN-ing; the exact numerical failure at cycle 1061295 (no floor, no dt drop) is NOT
understood -> B3 (per-cycle dumps over the last 30 cycles) launched 04:20.
Dirs: B1_catch_ctrl, B1b_ctrl_nandiagbin, B2_catch_tmin3200; analysis task_H2_B2_*,
task_H_census_*.

## PER-CYCLE CATCH 2026-09-10 04:40 (B3_catch_percycle 195367; bit-reproduces B1)
Trick: `-i extra.athinput` may ADD an <output5> block (main.cpp loads -i after the restart's
embedded input; any block named output* is accepted); last_time = start time; dt=30 < dt_run
=> one 39 MB dump per cycle. 85 dumps, cycles 1061211-1061295 (the death cycle included).
FINDINGS: single-cycle event, NO precursor (column stationary to 4-5 digits, M 0.09-0.12,
subsonic everywhere). What goes non-finite is the VELOCITY of 3 cells (i_arr 305-307 of the
seed column); rho stays finite (changes 0.01%); eint of those 3 cells drops 5 orders in the
same step to EXACTLY eps = 1.42264e7 erg/g = T 31.62 K = 10^1.5 = the EOS table logt_min edge
(the same fingerprint as the "floored" cells of the open-top deaths -> that fingerprint IS
the table clamp). Nothing crosses dfloor/pfloor. The bad column is a standing 1-cell spike
(rho 5.4x, T 2879 vs 9500 K of its 8 neighbours), neighbours untouched at the death cycle.
PERIODIC PRECURSOR: v_r(305) sawtooths: drifts -120 cm/s per cycle, then jumps +9.16e3 in
ONE cycle at cycles 1061218, 1061231, 1061245, 1061258, 1061267, 1061281 (period 13-14
cycles ~ 420 s); the death at 1061295 is exactly 14 cycles after the last jump = on the cycle
the next jump was due. => a periodically applied operator kills the cell. Identify it (agent
launched 04:45): which red_giant operator runs every ~14 cycles / 420 s (rt cadence? MLT
mean? WB cache? sponge?) and what it does to v and e at a cold overdense cell next to a hot one.
No rg.log in these runs (no log block in the embedded input). Files: _analysis_0910/
task_B3_percycle.py, task_B3_out.txt, task_B3_track.py, task_B3_track_out.txt.

## MECHANISM 2026-09-10 05:00 (agent source+dump analysis; established vs inferred marked)
NO periodic operator exists (every cadence parameter tabulated; dt constant to 8 digits); the
+9.2e3 jumps are a LOCAL limit cycle of the standing 1-cell spike (neighbours and the shell
median smooth to 6 digits at the jump cycle), amplitude slowly growing. The spike is a cold
blob (T 2879, rho 5.4x) sitting at v_r -1.8e4 inside a +1.5e5 OUTFLOW.
The kick: the only writer of u0(IM1) is the WB gravity source (red_giant.cpp ~1836-1843):
src = bdt*(A_{i+1}(pr-p)+A_i(p-pl))/V with pr/pl the cached background and p = the cell's
OWN eos.Pressure(d, w0(IEN)); wb_rmax unset -> active at 1.082 R; a cell whose p departs
from the background gets -2p/r: 304 cm/s^2 = 18x g. => the SAME WB-at-a-cold-cell mechanism
as the open-top 1.43114e6 death (FL4 proved wb_rmax alone cures that one).
The NaN: eos.Pressure(d,e) takes log10(e); a NON-POSITIVE e -> NaN -> src NaN -> u0(IM1) NaN,
rho untouched. The sponge guards this case (red_giant.cpp ~2607-2613); the gravity source
does NOT. Then GnomonicEquiangleRaiseVel: ekin NaN -> eint NaN -> e_positive false -> temp=-1
-> EnergyFromPressure(d,pfloor,temp) pins at the bracket 10^(ymin-3) = 0.1 K (eos_logt_min=2
here; the "31.6 K = 10^1.5" was the ANALYSIS table's logt_min) -> eint = rho*1.4226419e7,
rho-independent, temp=0.1 > tfloor so not overwritten, w(IVX)=NaN. Matches the dump exactly.
INFERRED (not proven; dumps are post-c2p): the non-positive e comes from the UNCLIPPED grey
RT energy increment (two_stream_rt.hpp ~2107-2114: rt_de_max=-1 -> no clip, no positivity
guard) in the preceding stage. efloor_from_ekin is absent from this binary.
CONSEQUENCE: prod12's wb_rmax=3.3e12 + efloor_from_ekin address this death too. Still to add:
a guard on eos.Pressure in the gravity source (mirror the sponge's) and decide on an RT
positivity clip. Tests launched 05:05: B4 (bin_nandiag + wb_rmax=3.3e12 from rg.00124), B5
(build_port binary + wb_rmax + efloor_from_ekin = the prod12 binary on the prod11 death).

## GUARD + B6 2026-09-10 05:40 (verified: hunk at red_giant.cpp:1845-1846, B6 0 FATAL)
Main tree: rg_grav curv branch now uses p = (e_>0) ? eos.Pressure(d,e_) : 0.5*(pl+pr)
(build_guard/src/athena). B6_guard (195370, rg.00124, WB on, no floor fix, no opac_tmin):
passes the death, runs to tlim 3.215e7 (cycle 1062769, 1474 cycles past), dt 30.14. The cold
clump persists (rel 0.20, rho 5x, v_r -1.1e5 at i 306-308).
KEY LOG LINE: "sponge guard FIRED (rank 59, cycle 1061446, t=3.21101e7): (0,7,8,305) d=1.4e-8
ei=-nan ke=36.06". The sponge computes ei = u0(IEN) - ke with ke FINITE => u0(IEN) ITSELF is
NaN before the sponge (which runs before ConsToPrim). So an energy-updating operator (RT most
likely; conduction possible) produces a NON-FINITE energy in that cold dense spike; the floor
pass then pins it to the 0.1 K bracket (the 1.42264e7 fingerprint) and the cell carries on.
ROOT CAUSE STILL OPEN: which kernel writes the NaN energy and why (opacity lookup at T~2900 K
rho 6e-8? a division in the grey sweep?). Instrument: a device counter + first-cell report of
non-finite du in the three grey/ck/generic sweep kernels (see [[two-stream-rt-three-kernels-trap]]).
Unguarded eos.Pressure(w0(IEN)) sites listed by the agent (MLT sweep 1896, shell mean 1971/2392,
MLT faces 2215, inner wall 2356/2426/2449/2594, rg_relax 2718, fill_open 2846/2884): a NaN in
the shell-mean reductions would poison a GLOBAL mean.

## ROOT OPERATOR FOUND 2026-09-10 05:10 (B7/B8 nan_report instrumentation, agent; raw lines quoted)
<problem>/nan_report (default off; scans after every operator; in main tree: red_giant.cpp,
two_stream_rt.hpp, hydro_tasks.cpp, hydro.cpp/.hpp, conduction.cpp/.hpp) at cycle 1061446:
FIRST offender = [hydro_flux_x1] face (0,7,8,306): flx(IDN)=3.2e-3 and flx(IM1)=5705 FINITE,
flx(IEN)=NaN; 2 bad faces (306,307) bracketing the cold cell 306 -> 3 bad cells 305-307.
Face state: L d=1.41271e-08 e=44840.6 T=8941 p=8657 v1=-7407 | R d=5.80803e-08 e=62289.4
T=2874 p=8413 v1=+6050 (diverging, 4x rho jump, 3x T jump). Conduction is INNOCENT (F_cond
2.7e8 finite added to an already-NaN flux; w_tau 0.0106). RT grey apply never produced a bad
cell; MLT off; psrc null. Operator order in a stage: Fluxes(+conduction) -> RKUpdate ->
GnomonicSrc -> user src (rg_grav, wall, sponge, RT, relax) -> BCs -> ConToPrim (floor).
=> the GENERAL-EOS RIEMANN SOLVER's energy flux goes NaN at a hot/cold diverging face.
B9_rsolver (agent, 05:15) instruments the solver's energy-flux path (star state / EOS call).

## B4/B5 VERDICT 2026-09-10 05:20 (checked by me in out.txt): the NaN cures do NOT hold the run
B4 (bin_nandiag + wb_rmax=3.3e12 only): passes 3.2106e7 but dt decays 30 -> 1.3 s and it dies
FATAL at cycle 1067883, t=3.225763e7, 1 cell (0,2,2,310) rank 80, gid 80, panel 5, lloc (0,0,0)
= a CUBE-VERTEX/corner cell, i 310 (NOT the old seed column).
B5 (build_port + wb_rmax + efloor_from_ekin + log): no FATAL, but dt collapses gradually
30 -> 0.24 s from ~3.2262e7 (no "COLLAPSE" line: the ratio trigger misses a gradual fall);
rg.log floor counts jump x10 (dfloor 2.3e6, efloor 3.4e6 per interval, efloor_de 2.2e10 at
cycle 1070097); cancelled at 3.2269e7 (05:20) — would never reach tlim.
B2 (opac_tmin=3200) stayed at dt 30.14 to 3.26e7. => in the lidded run the WB/floor/guard
fixes only remove the NaN death path; the underlying cold collapse (photospheric cells
falling through the molecular-opacity knee) proceeds and chokes dt ~1.5e5 s later.
The OPACITY CLAMP is a necessary part of prod12, not optional. Analysis of the B4/B5 dt
collapse (which limiter, which cells, census growth vs B2) launched 05:25.

## EXACT OPERATION 2026-09-10 05:35 (B9_rsolver capture, verbatim in the agent report)
rsolver=hllc, plm, general_eos=table, wb_x1=true, wellbalance_dynamic. hllc_hyd.hpp general
branch makes NO EOS call in the energy path. At face (0,7,8,306): the RECONSTRUCTED right state
(cell 306 side) arrives with e = -nan and p = 1e-12 (= pfloor: hydro_fluxes.cpp:253
`dr(IDPR,i) = fmax(dr(IDPR,i), pfloor)` masks a NaN pressure because fmax(NaN,x)=x; IEN has
no mask). Cell-centred w0(IEN) on both sides FINITE (44840.6 / 62289.4). S_M=+3.0e5 -> qd=0,
but flx(IEN) = qc*fl.e + 0*(-nan) + ... = NaN (IEEE 0*NaN). Momentum flux uses only d,v,p ->
finite. Two defects: (a) HLLC: guard wl/wr(IEN) before el/er (mirror the qa/qb/cp clamps);
(b) UPSTREAM: the x1 reconstruction returns NaN p AND e from finite cell values -> the WB x1
reconstruction is the suspect (the fmax mask should be an explicit non-finite check).
Instrumented: hllc_hyd.hpp, hydro_fluxes.cpp, hydro.hpp, hydro_tasks.cpp (nan_report, off by
default). Next (05:40): instrument the WB x1 reconstruction at that face.

## B4/B5 dt COLLAPSE DIAGNOSED 2026-09-10 06:10 (agent; task_B45_*): CONDUCTION dt in cold cells
Correction: dt is FLAT 30.14 (= 0.2557 x the deep hydro CFL 117.9 s at i_arr=2, unchanged in
all runs) until cycle 1066300 (t 3.2256e7), then drops x73 in <=100 cycles in BOTH B4 and B5
(bit-near twins), and stays 0.24-1.4 s. Hydro CFL in the cold shell is 4000-9700 s, so the
limiter is conduction (by elimination; no COLLAPSE print because the ratio trigger missed it).
Cold population without the clamp: 4 -> 156 cells in 1.4e5 s, a distinct 1350-1400 K
population (then B5: 48 cells < 1000 K, min 95 K); with opac_tmin=3200 (B2): 2 -> 102 in 5e5 s
and NEVER below 2190 K (a floor at ~2200-2250 K). B4's death cell is EXACTLY a cube vertex
(gid 80 panel 5 J=0 K=0) but B5's collapse cell and all min-T cells are INTERIOR (d 9-10);
vertices are one site among many. Cold cells are always ONE cell, v_r -1e5..-3e5, e/rho
8e9-1e11 (100-10000x above the 1.42264e7 mark), between two normal-T neighbours.
=> two separate problems: (1) the WB-walk NaN (fixed, see [[red-giant-wb-polytropic-walk-overflow]]);
(2) the single-cell cold collapse through the opacity knee -> conduction dt; opac_tmin=3200 is
the only tested remedy (B2), and it leaves a 2200 K floor population + a growing downdraft.

## B12 = THE LIDDED REVIVAL GATE, PASSED 2026-09-10 06:55 (job 195378, agent a564796a)
build_guard binary (WB fix + guards) + opac_tmin=3200 from rg.00124 (3.205e7) -> tlim 3.45e7
cleanly: dt 30.143-30.146 s throughout (0.01% spread), 0 COLLAPSE, 0 nan_report, no FATAL.
Tracks B2 cell-for-cell to 3.26e7 (90 vs 102 cold cells, minT identical). Beyond: cold cells
90 -> 118 -> 141 -> 178 -> 200 at 3.45e7 (growth 16 -> 5.8 per 1e5 s, non-monotonic =
recycling, NOT saturating yet); VERTEX cold cells exactly 0 always; two 1-2-cell episodes to
~180 K (i=312) each recovering to ~2000 K within 2e4 s with NO dt response (the clamp removes
the runaway feedback of B4/B5). Seed-column downdraft sits on a 0.45-0.51 plateau, does not
deepen. max|v| 8.2e5 -> 9.8e5 -> 9.0e5; L_out/L 1.29-1.36; L_cut/L 0.77 -> 0.56.
Restarts in B12_fix_tmin3200/rst/ -> prod11 can be REVIVED from there (user decides; wall
per 2.5e6 s ~ 1.2 h on 6 nodes).
