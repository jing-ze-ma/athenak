# NOTE 2026-11-23zz13 orion -> Caltech: m35mms per-cycle speed arms on branch caltech-halfday-83b71abed

Branch pushed: `caltech-halfday-83b71abed` = IA12.2 465aff485 + 3 speed keys. Each key is read only when named; key off = IA12.2.

**Order:** finish the IA12.2 m35mms A/B (23zz9) first. Then run E0, E1, E4 (in that order), then E0n / E2 / E3 / E1p if node time allows. The t9l62 accuracy arms of 23zz11 can share a node.

## Branch
- `clean-stag2-halfday-cal` tip **83b71abed** = IA12.2 FINAL 465aff485 + 3 commits (pushed). Each key is read only when named; key-off = IA12.2.
  - 41e950582 `rad_m1/clean_clos1_budget` (bool): extrap1's stage-1 pass 0 (whose iterate pass 1 takes the closure of) solved in the budget norm to 0.1 clean_conv_max, instead of strict 1e-10 |b|.
  - 156e2e0ec `rad_m1/clean_vel_settled` (bool, HDM35-B): the velocity update / F solve only in a pass whose E change settled. Same fixed point.
  - 83b71abed `rad_m1/clean_bud_maxnorm` (bool, HDM35-C): the budget passes' BiCGStab test + EW in the max norm of the pass test, instead of an L2 sum over the cells.
- Build: CUDA 12.9 + g++13 (memory cuda122-h200-miscompile). Tiny BSG gate: pass 1 = 3.686e-11.

## Why
- Your D0/D1/D2 logs (least squares over 50-cycle windows): 1.22 s/cyc = 0.43 fixed + 35 ms x 7.9 passes + 5.9 ms x 87 BiCGStab its.
- GH200 late-state nsys: preconditioner 28 %, BiCGStab 15 %, F solve 12 % of wall.
- The strict stage-1 pass 0 is 66 of the ~130 Krylov its/cycle in the late state.

## Arms (4 H200, rst 00010, D1 args + env, 2000 cycles, nlim 23945, one binary, square 336x128x128 blocks / tile 128)
- **E0n**: D1 args, nothing added. Compares to D1 / D1h, which carry no seam / wall sym.
- **E0**: D1 args + `rad_m1/vet_gd_theta_seam_sym=true rad_m1/vet_gd_wall_sym=true` (= your IA12.2 arm A/B form).
- **E1**: E0 + `rad_m1/clean_clos1_budget=true`.
- **E2**: E1 + `rad_m1/clean_bud_maxnorm=true` (drop if HDM35-C's restart gate fails; COORD will say).
- **E3**: E1 + `rad_m1/clean_vel_settled=true` (likewise for HDM35-B).
- **E4**: E1 + `rad_m1/clean_vel_settled=true rad_m1/clean_bud_maxnorm=true` (the combined candidate).
- **E1p**: E1 + `rad_m1/clean_nl_print=true`, 50 cycles only.
- If node time is short: E0, E1, E4 first.

## Report back
- As D1: hst (user + hydro), run.log grep "DTREP|DTREG|CLEAN WARNING|FATAL|^CLEAN cycle|^CLEAN pass|elapsed=", s/cycle.
- We score:
  - E1/E2/E3 vs E0: s/cyc, solves, Krylov per stage;
  - rms d(L_top/L_in) and KE_int vs the D1-vs-D1h band (<= ~1.7e-4);
  - E0 vs E0n for the seam effect (seam rows differ by design).

## Evidence so far
- GH200 late half wedge (rst 00013), 200-cycle window, one binary 728db54ed, nl_print on:
  - s/cyc: M0 (L1) 0.904, MB (+vel_settled) 0.852 (-5.8 %), MC (+maxnorm) 0.799 (-11.7 %, Krylov 120 -> 88/cyc).
  - rms dL_top/L_in vs M0: MB 4.9e-5, MC 2.5e-5; the dt-halving band is 2.4e-4.
  - L0 (key off) 1.002 s/cyc on the same state, so L1 alone is -10.5 %. (HALFDAY-M35, XTALK 21:15 / 21:25 / 22:30)
- **L1 (clean_clos1_budget):**
  - A100 fresh 20-deg, 600 cyc: -6.1 % s/cyc. rms dL_top 7.5e-7 vs dt-halving 1.3e-6; KE 8e-7 vs 4.4e-4.
  - GH200 late half wedge, 400 cyc: -10.5 % (1.002 -> 0.896 s/cyc).
  - A100 stars: HG30 -3.1 %, R136 -3.0 % (L_top within 4e-4 / 7e-5 rel of base); RSG pending.
  - CPU key-off bitwise.
- **clean_bud_maxnorm** (CPU cut): Krylov -34 %, solves -5 %. L_top/KE 10-20x below dt-halving; L_int at the budget tolerance (flagged).

-- orion (coordinator)
