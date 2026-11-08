seamE-1106 overnight report. I couldn't write REPORT.md myself (harness refused a report file); please save this as /resnick/groups/carnegie_poc/jingze/seamE_1106/REPORT.md.

Branch seamE-1106, local only, tip 7de51b303, worktree /resnick/groups/carnegie_poc/jingze/wt_seamE (based on meshfb-1106). No merge, no default change, no production touched. Every key is read only when named. Keys-absent runs are bitwise with prod-1102 on GPU (HG small rst, 4 cycles: IDENTICAL at sd4, sd18 and sd19).
All arms restart from HG rst 139 on 6 H200 and are compared per simulated time against the control meshfb_1106/hg/off. Run dirs are in seamE_1106/runs; analysis scripts (simt.py, budget.py, diag.py, fp/cmp.py) are in seamE_1106.

**What drives E<=0 (row dumps unless marked)**
1. Seam cross term (the "seam" problem). A seam-adjacent face's skew cross term uses a one-sided 3-point derivative (M1CsTanDer, rad_m1_implicit.hpp).
   - Across a sharp front (panel 0 j=100 k=0: E = 20, 174, 3390) it has the wrong sign and 5x the size, so F = 6.5 cE on that face.
   - 55 of 60 pass-0 minimum cells are seam-adjacent, and 53 of those have the derivative out of bounds.
   - The fix `rad_m1/cs_cross_mono` clamps it: bad rows 1253 -> 156, seam share 91 % -> 2 %.
2. Opacity-Newton term on the RHS: about -1e4 on a row with E ~150 (confirmed by switching the term off). HG's guard mode 6 has no RHS test. New guard bit 8 (mode 15) helps; `implicit_opac_newton=false` removes the term.
3. Hot T predictor: in 12 of 18 stage rows pass 0 starts at T up to 1.7x the gas temperature, which makes the source term SRCR negative. `rad_m1/implicit_predictor_tcap=1.1` removes this class. In the 5-cycle real-state test it changed nothing.
4. Vertex holes (hydro). With cs_cross_mono on, 19 of 19 NOT POSITIVE cycles have a negative-eint hole in panel 3 vertex cells (555-778 km/s).
   - In production dumps such holes first appear 1.5 h after cs_seam_rho_guard=0 (11-03 22:49) and keep reappearing at new columns.
   - RSG had the same switch and shows none, so the cause is open.
   - `mesh/cs_seam_rho_guard=4` gives 0 holes in every arm.
   - The NOT POSITIVE rate per sim time shows no lasting step at the switch, so the holes are not what sets that rate.
5. Explicit stage-1 increment: EN = E + a21 dt K1_E < 0. K1's face fluxes push about 60 E per step through a cell.
   - `time2_k1_reset=pos` (existing key) creates energy, +2 % of L: rejected.
   - Face-flux FCT and the neighbour-relay version both failed (cascade / too slow).
   - SDIRK2 tableau: stage 2 goes negative, energy budget off by 0.5 % of L, first order in a time ladder (E error 1.6e-3 vs 5.7e-5 for hesdirk2 at dt 10). Dead.
   - `rad_m1/time2_k1_redist=0.1`: a conservative radial-column transfer that works (budget closed). It is a test vehicle: FATAL on x1-split decks, and the face-flux form you asked for doesn't work at this dt.
   - It fires on ~1e5 cells per stage (0.7 % of the mesh) and moves 2.0e37 erg/s, 15 % of L gross, about 1.6 cells. This is the largest open accuracy question.
6. Global fallback: control has D_r,lat active in only 28 % of BE passes. `implicit_fb_local` (meshfb-1106): vlat_fall per cycle 1.5 -> 0.002.

**Key table**

| key | fixes | gates |
|---|---|---|
| cs_cross_mono (k=0) | seam E<=0 | keys-absent bitwise; layout no worse than off; CPU restart identical (audit); Linf first order at smooth extrema (audit). Threshold k>0 (0.1, 0.02, 0.005) brings NON-CONVERGED retries back, so k=0 is used. |
| implicit_opac_newton=false | Newton RHS | same fixed point as Newton on with opac_lag=none (1e-12); lag error 20x smaller at equal dt |
| cs_seam_rho_guard=4 | vertex holes | layout = off |
| implicit_predictor_tcap=1.1 | hot T start | inert in the fixed-point test |
| time2_k1_redist=0.1 | stage EN<0 | budget closed; layout = off; fires widely; column-only |
| implicit_fb_local | global fallback | no layout gate; GPU restart unsettled (below) |
| guard_mode bit 8 | Newton RHS hidden by TRHS | used only with Newton on |

**Candidate r7 (all of the above, k=0) to t=5.018e7 (68000 s sim) vs off**
- BE redos 1 vs 153; NOT POSITIVE 1 vs 184; dt 107 vs 81.
- Wall 1481 vs 2180 s (-32 %); Picard 5.8 vs 10.0 (max 15 vs 45); 0 retries.
- Budget closed: 4.93e36 vs L_in-L_top 4.93e36.
- Noise member r7p (E x 1e-10): differences ~0.01 %, so every change vs off is real.
- At 1 d vs off: E(i280-320) +2.45 %, deep KE -0.41 %, L_top +0.9 %, holes 0, vmax 138 km/s.
- meshfb's local-only arm shows the same signs (E +1.9 %, KE -0.76 %), so most of this comes from keeping D_r,lat on. Which result is right is not shown.
- r8 (r7 + opac_lag_passes=3): wall +3 %, 2 redos; E(i280-320) a further +0.15 % at 1 d, above noise.

**Opacity lag (new finding)**
- HG's deck has `implicit_opac_lag=pred` with 2 passes: kappa is frozen after pass 2, so the converged answer depends on the path.
- On the real state with lag none, Newton on and off agree to 1e-12. With the lag they differ by 4e-3 at dt 60.
- Lag error vs each config's own lag-none truth (E median / p99, i240-360):
  - candidate: 8.3e-5 at dt 78, 1.7e-4 at dt 105, 1.9e-5 at dt 105 with passes=3;
  - production config: 3.6e-5 at dt 30 (truth converged). At dt 78/105 its lag-none truth never converged (300 passes), so those rows can't be trusted.
- The candidate's larger dt raises its lag error about 2x; passes=3 more than recovers it.
- The same lag causes the layout excess: Newton off + lag gives dens 3.9e-5 at b2x2 vs 1.1e-7 with lag none; key-off with lag none gives 9.7e-8.

**Other gates**
- cs_label FATAL was a false alarm (halo depth 2); fixed in opcheck.
- GPU restart test: DIFF even for the keys-absent base, so it is my harness or the hesdirk2 restart path, not the keys. Unsettled.
- Smooth-order test of the limiter: the CPU run timed out, not done.

**Decisions for the user**
- (a) Candidate for a longer test: cs_cross_mono, implicit_opac_newton=false, cs_seam_rho_guard=4, tcap 1.1, time2_k1_redist 0.1, fb_local, probably with opac_lag_passes=3.
- (b) Accept or replace the column redistribution: it is energy-conserving, but it moves 15 % of L gross and only works with one radial block.
- (c) The 0.9 % L_top / 2.5 % near-photosphere E change is mostly D_r,lat kept on. Is that the intended physics?
- (d) Settle the GPU restart method before any merge.

**GPUs**
- 02-17:0-3 + 02-11:2,3 handed to a0bf at 03:57.
- r7p ends by itself at ~04:17 on 01-17:0-3 + 02-11:0,1, then those are free.

---
## Addendum 04:30 (production GPU restart test)

Production restarts are not bitwise on GPU under hesdirk2. They are bitwise under backward Euler, writing a restart doesn't perturb the run, and neither the predictor nor the opacity lag causes the difference.

Setup: prod-1102 binary athena_gpu_he_star_m1_wallwb_0f1bbb708_fmadoff, production HG keys, HG small rst (prod_1103b/he/runs/hgb), 1 GPU per run. Each test compares 4 straight cycles with 2 cycles + restart + 2 cycles, on the bytes after `<par_end>` of the final restart files.

| test | result |
|---|---|
| hesdirk2 (production scheme) | DIFF: about 1.08e9 of 6.67e9 bytes differ, so the state itself, not metadata |
| time_scheme=be | IDENTICAL |
| writing a restart at cycle 2 vs not writing one (4-cycle runs) | IDENTICAL |
| hesdirk2 with implicit_predictor=none | DIFF |
| hesdirk2 with implicit_opac_lag=none | DIFF |
| hesdirk2 with both | DIFF |

The cause sits in the hesdirk2 restart path itself. One candidate is the logged rule that restarts without the slope take a backward-Euler first step; others are stage or closure history that the restart file doesn't hold. Not localized further before the allocation ended, so not proven.

This also means the earlier restart DIFFs with tonight's keys (base, without fb_local, fb_local alone) say nothing about the keys. A restart gate for any of them needs time_scheme=be, or a fix of the hesdirk2 restart first. Results: /resnick/groups/carnegie_poc/jingze/seamE_1106/RESTART_PROD.txt and RESTART_PROD2.txt (run dirs rsp/ and rsp2/).

GPUs: none of my processes is running. 02-17:0-3 and 02-11:2,3 were empty when these tests started at 04:25 and finished by 04:29. r7p on 01-17:0-3 + 02-11:0,1 ended by itself as planned. The three remaining srun processes belong to rsgdiag_1107, not to me.

My earlier handback (the overnight report, to be saved as seamE_1106/REPORT.md) stands as written, except its GPU restart line: the DIFF is now settled as a production hesdirk2 restart property, not caused by the seamE keys.
