# RSG E<=0 fallback diagnosis (rsgdiag-1107), 11-07 night (requested by the coordinator)

J=/resnick/groups/carnegie_poc/jingze, D=$J/rsgdiag_1107. Branch rsgdiag-1107 in $J/wt_rsgdiag (off seamE-1106 55db30f95, local commits only, no push, no merge).
State: RSG rst 01295 (cycle 589788, t=5.4648e8), the same state as rsg_on1/off1/on3/off3.

## 1. Rates (measured, counters in rsg.log / run.log)
| run | vlat fallback/cycle | NOT POSITIVE/cycle | NOT ADMISSIBLE/cycle | posfix cells p50 |
|---|---|---|---|---|
| prodH (no keys), recent hst windows | 1.6-1.9 | 0.36 (whole log) | 0.18 | 171 |
| rsg_on1 (cross_mono + gm7, tan=true) | 0.85 | 0.25 | 0.13 | 38 |
| rsg_off1 (same keys, tan=false) | 0.87 | 0.25 | 0.12 | 38 |
| HG (seam agent, for comparison) | off 1.5 -> keys 0.13 | | | |
On RSG the keys halve the fallback rate (on HG: 11x fewer). They cut the posfix cells per event 4x but not the event rate. tan on/off makes no difference.

## 2. Where the bad cells are (measured, rsg_on3/off3 fbdiag, fb_local on, ~35 cycles, 3600 records; ana_on3.txt ana_off3.txt)
- Seam distance: 78-93 % more than 5 cells from a seam (82 % of the area is). Seam-adjacent cells (distance 0): 0-0.1 %. Vertex 0-0.3 %. Block edges 5-9 % (6 % of the area). The seam and vertex chain found on HG does not explain RSG.
- Radius: r/R_ph 0.90-1.00 (99 %), i 400-455; r/R_ph 0.95-0.98 holds about half.
- tau_top 100-1e4 for 70 % of cells. dtau_cell is bimodal: >10 for 75 %, 0.003-0.03 for the rest, within one column: one cell has dtau 136-1000 next to cells with 0.003-0.01.
- Solve kind: hesdirk2 stage1 and stage2 only. No BE solve was bad in this window (fb_local arms). Pass 0 is the first bad pass in about 80 % of the solves.
- Persistence: clusters recur over several cycles: (panel,j,k) = (1,35,86-87), (5,46,93), (5,60-63,94-95), (0,116-121,2-4), (4,39-40,111-112), (1,92,75). 18 columns are bad in 6-20 of about 35 cycles.
- Terms-off passes (vlat and vimp dropped): 85-90 % are cells that were already bad earlier in the same solve, with E_prev at the floor (1e-10). Dropping the terms does not cure them (as on HG).

## 3. What the cells are (measured, prodH float32 dumps 06325 = the rst cycle; hotstate.py colprof.py p0cmp.py)
- 100 % of the pass-0 bad cells are cold, H-recombined gas (e/rho 3-8e11, against 1.5e13 in the ionized shell; cold cells are 19 % of the i 420-455 cells).
- They are cold downdraft fingers inside the hot shell: v_r -1 to -3.4e6 cm/s (3rd percentile), rho at the 3rd percentile of the shell (2.9e-10 against <rho> 1.1e-9), E 1-10 against 1e3-4e3 in the hot gas, F/cE 0.05-0.5. Example column (1,35,87): E drops 1740 -> 18 from i 428 to 436 (2 cells) and e/rho 1.4e13 -> 7.6e11.
- They are not front cells: max neighbour E / own E is 1.1 at the median, >10 for 1 %.

## 4. Mechanism (partly measured, partly inferred; GPU rows pending)
Measured (cpu1 = CPU sd10c restart; only the valid fields: T2 increments, gas state, TRIG):
- (a) The stage old vector is negative. Cell (p0,j59,k40,i430): E^n 14.8, a21 dt K1_E = -15.5. Cells (p1,j35-36,k84-85): the gas old vector e + a21 dt K1_EN is below 0 (TEGN = -1; K1_EN -0.43/s against e_gas ~100-140, (1-gamma) dt = 677 s).
- (b) The predictor start is far. Pass-0 start E is 2-12x E^n (4.9 -> 60.6; 14.8 -> 33.5): the implicit_predictor=step extrapolation of the last increment. The cell with dtau 136-1000 sits next to cells with dtau 0.003 in the same cold column, so its T iterate is hot: the H-recombination opacity rises ~1e4x between ~4000 and ~8000 K.
Inferred: in cold cells E is 10-100x smaller than e_gas. An explicit exchange slope or a predicted T that is too hot then takes E or e_gas below 0 in one step, and the global fallback fires. HG cells are radiation-rich and seam-bound; RSG cells are gas-dominated, cold and interior. This also explains why time2_k1_stiff hurt RSG on 10-29: it moves the stiff exchange slope into K1_E, which is tiny in cold cells, so E + a21 dt K1_E < 0 (inferred).

## 5. Tools and fixes
- Bug in the CPU diagnostics (meshfb/seamE `fetch`): create_mirror_view on a HostSpace build returns the device array itself, so every fetched field aliased the last one. CPU fbdiag records were empty (min E 1e300) and ROW/CELL values were wrong. Fixed in rsgdiag-1107 89d801876 (create_mirror). GPU dumps were never affected.
- Row dump extended to every E<=0 cell of the rank (cap 40): binary rd2.
- New key rad_m1/implicit_predictor_ecap = f (default off, bitwise when absent): the predicted E is accepted only within [E0/f, E0 f], and T is predicted only where E is. Commit 246e5e400, binary rd3.

## 6. GPU arms 04:02-04:38 (RSG's 8 GPUs, rst 01295, 4 GPUs each; all ended rc=0 by -t, nothing left running)
All arms: XKEYS + cs_cross_mono + guard_mode 7 + fb_diag 300, no fb_local; binary rd3/rd4 (seamE 55db30f95 + rsgdiag; rd4 = rd3 + tstart key; absent keys = default).
| arm | cycles | vlat fb/cyc | NOT POS/cyc | NOT ADM/cyc | stage1 E<=0 solves | stage2 E<=0 solves | bad cell-passes | Picard passes/solve s1/s2 | posfix cells p50 |
|---|---|---|---|---|---|---|---|---|---|
| base (rd3, + row dump) | 150 | 0.78 | 0.29 | 0.13 | 40 % | 39 % | 3909 | 6.65/5.38 | 46 |
| ecap2 (implicit_predictor_ecap=2) | 200 | **0.41** | 0.22 | 0.085 | 27 % | 11 % | 464 | 6.64/5.62 | 14 |
| prednone (implicit_predictor=none) | 250 | 0.70 | 0.24 | 0.052 | 42 % | 19 % | 743 | **9.96/10.71** | 21 |
| tsold (implicit_tstart_old=true) | 300 | 0.88 | 0.32 | 0.13 | 51 % | 33 % | 11120 | 6.82/5.28 | 79 |
- **ecap2 halves the mesh-wide fallback (0.78 -> 0.41/cycle), cuts bad cells 8x and posfix cells 3x, and leaves the Picard passes unchanged.** The high-dtau (hot-T) bad cells disappear (dtau>100: 52 % -> 6 %). What remains is tau_top 3-30, dtau 1-10, and 24 % within 6 cells of a vertex (panel 0 j115-122 k3-21; panel 4/5 interior).
- prednone removes the same cells but costs +50 % Picard passes (s/cycle 3.8 vs ~3.2): not a fix.
- tsold makes it worse: the gas old-vector T start (hypothesis (a) for gas) is refuted as the main trigger.
- Row split (base, $D/rows_base.txt, 63 pass-0 stage rows): the most negative RHS term is SRCR (the exchange or opacity-Newton RHS) in 64 % of rows, SRCR/|E_P| p50 -15; lateral lag and seam faces are negligible (worst face never a seam face; max single-face lag 1-4 % of E_P). Tk/T(EGN) p50 1.55, p90 1.9-2.0 (just inside the default tcap 2). This matches HG finding (3) and dominates on RSG.
- Energy budget: NOT measured. My arm script disabled the hst output (output1 dcycle 1e6), my mistake. Counters: efloor_de of the same order in all arms; pos_flr/cycle base 0.21, ecap2 0.00. ecap moves only the Picard starting point, so the fixed point and the conserved write-back are unchanged (inferred). A 1-2 h ecap2 vs base A/B with the hst on is the needed gate.
- Noise: one member per arm, 150-300 cycles; the vlat counts are 116 vs 82 events (Poisson ~10 %). ecap2's -47 % is well outside that; prednone's -10 % is not.

## 7. Recommendation
1. implicit_predictor_ecap=2 (default off; branch rsgdiag-1107 246e5e400; it also caps T where E is not predicted) is the RSG fix candidate. Next: a GPU A/B with the hst on (energy budget, Etot drift, Picard, s/cycle), 2 members, then try ecap 1.5. HG check too (its trigger (3) is the same term).
2. Keep cs_cross_mono + guard 7 (they cut posfix cells 4x on RSG: prodH -> rsg_on1).
3. Do not use time2_k1_stiff on RSG (inferred cause in sec. 4) or implicit_tstart_old (measured worse).
4. The CPU-diagnostics aliasing fix (89d801876) should go to seamE-1106/meshfb-1106 before anyone trusts CPU fbdiag/seamrow output.
Open: cpu1 (my sd10c CPU restart on 02-09, 8 x ~50 GB) still runs. My stop was denied by the permission classifier; the coordinator asked the user. It ends at nlim (60 cycles, never reached) or with the allocation.
