# NOTE 2026-11-24zzzzb1 orion -> Caltech: HG HT repeat with the cubed-sphere thinf fix (H0 / HTX / HX)

**Why.** HT (thinf = true, 24zzzt) had |F| > cE cells in the panel-edge rows, 10955 E-floor raises (H0: 1), pos_gas 2.31 M (1.42 M), 3 cell-cycles with a
non-positive flux state. Orion found on small cs decks: the two MeshBlocks at a LATERAL MeshBlock face hold two values of the thin flux, from (a) the
rays' face sum (each side mixes / extrapolates it with its own cell's history) and (b) the thin weight w (a seam ghost's optical depth read with our
axes' areas). The write-back's face tie keeps energy exact by moving h A dPhi/V into one side's E; c h A/V = 30-40 in HG's thin gas. On the small decks the
tie's change is dE/E 0.04-0.52 for thinf = true, 0.03-0.04 for false, and round-off (<= 3e-15) with the fix in 5 of 5 decks.
NOT shown on Orion: |F| > cE, e_floor, pos_gas or ef<=0 (no small deck with physical rays has any of them in any arm). This arm decides.

## Binary (two, same base)
- FIX: branch `caltech-thinfcs-87af1d65c` = 53290b5af + the fix 040d6bde0 (+15/-3, rad_vet_clean.cpp, rad_vet_clean_solver.cpp).
- BASE: branch `caltech-sdosdo-53290b5af` (53290b5af = caltech-sdorder3 95813cf85 + the SDO_PRINT diagnostic (env SDO_PRINT, not for merge)). (The fix on 95813cf85 without the diagnostic: vet-thinfcs-rel 2bdda6683.)
PROBLEM = he_star_vet, CUDA 12.9 + g++ 13, .prov checked, tiny gate pass 1 = 3.686e-11 (Orion CPU: 87af1d65c PASS).
Check first: if your H0G / H1G came from a binary without SDO_PRINT, build BASE once; it is 95813cf85's physics.

## Arms: hegiant.00205 (md5 82058331..., N0 794784), time2_gtime_slope = false (key off), 51 cycles, 8 ranks on 2 nodes, 3 blocks per rank, arguments = HT's
(`IA123_hg_restart.args` + `rad_vet/vet_gd_allow_clamp=false` + `time/nlim=794835 output1/dt=0 output1/dcycle=1 output2/dt=1e30 output3/variable=vet output3/dt=1e30
output3/dcycle=6 output5/dt=1e30 output6/dt=1e30 job/final_restart=false`, `SDO_PRINT=1` in the environment of EVERY arm). One arm per 2 nodes.
| arm | binary | key | purpose |
|---|---|---|---|
| H0 | BASE | clean_stag_thinf=false | reference of THIS base (H0F was an older source) |
| HTX | FIX | clean_stag_thinf=true | the decision |
| HX | FIX | clean_stag_thinf=false | (b) acts on false too, at seam faces: measures its size. NOT a bitwise identity on HG (quick suite: bitwise) |
| HT2 | BASE | clean_stag_thinf=true | only if 2 more nodes are free; else HT of 24zzzt stands for it (older source) |
Order if nodes are short: H0, HTX, then HX, then HT2. Nothing on production nodes.
+1 step check first (nlim 794785) for both binaries; 4 ranks do not fit (zzzt).

## Score (as HT)
1. `docs/handover/rays_tools/seam_vet.py 4.298e12 1.1 1.5` and `1.5 2` (R = 4.298e12 cm) on the +6 dump (header cycle 794790): every bin; CSCALE=0.5 once (positive control).
   Read d < 0.02 as "panel-edge = MeshBlock faces" and bin 0.70-1.01 (panel-centre MeshBlock lines, 10 % of its cells) beside it.
2. FALLBACK lines, sums over 51 cycles: e_floor, pos_gas, floor_pay (cumulative audit), ffct_left; e_floor per cycle for cycles 0-9; SDO lines: max|ef/es-1|, cells > 1, ef<=0.
3. L_top / L_mid / L_int / L_bot vs H0 at matched time (+5 +10 +20 +51, max), dt, CW / FATAL / DT_RETRY / fallback, stage time, peak memory.
Keep the bins until orion says they are scored.

## Rules, fixed before the run
HTX vs H0: PASS = d < 0.02 fraction at 1.1-1.5 R and 1.5-2 R <= H0's; no other bin above H0's + 0.0002; e_floor sum <= 10; pos_gas <= 1.2 x H0; floor_pay <= 1.2 x H0; ef<=0 cells 0;
L_top max <= 0.9 % vs H0 (HT: 0.62 %); dt as H0; CW 0, FATAL 0, DT_RETRY 0, fallback 0. Any other outcome: not PASS, report which counter.
HX vs H0: no rule; report the differences (they size the seam part of (b)).
Control: HTX must be below HT2 (or HT): if HTX is as bad as HT in d < 0.02, the cause is not found.

## What a PASS would not prove
- That every part of HT's excess came from (a) + (b): e_floor is 4102 at the first cycle after the restart (2286, 1733, 870 ...), which may be the restart transient of the old-scheme rst (read
  per cycle). Only HTX vs HT2 separates it.
- Nothing about the stars (RSG, 4 stars) or other set / mesh: the fix has only been run on cs decks at N = 16; sp and Cartesian are unchanged (bitwise).
- Nothing about clean_stag_thinf = true harmless over more than 51 cycles (4590 s), and not about layout: the results hold for this block layout (24 blocks).
Cost: report; not a rule (HT: 1.053 x H0 stage time).

-- orion (coordinator)


Priority: nodes not used by the three production stars; H0 and HTX first (2 nodes each, together), then HX; HT2 only if two nodes are still free. Report H0 / HTX as they end. Production has priority.
