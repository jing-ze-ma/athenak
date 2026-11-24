# NOTE 2026-11-24s orion -> Caltech: adopt clean_stag_mg_levels = 4 (PL4 -9 % on late m35mms)

Your PC arms (late m35mms, rc3 + speed keys, 4 H200 each):

| arm | s/cyc | Krylov it / stage |
|---|---|---|
| P0 (rb, mg 2) | 1.060 | 7.15 |
| PR (rbr) | 1.069 | 5.99 |
| PRL4 (rbr + mg 4) | 1.005 | 4.86 |
| **PL4 (mg 4)** | **0.963 (-9.2 %)** | 4.93 |

The preconditioner only changes how the solve converges (same stop tolerance), so the answer stays within the solver tolerance (your same-time L_top row: 1e-5).

**Do:**
1. At each star's next rst (bsg, t9l62, m35mms, M35ZAMS), add `rad_vet/clean_stag_mg_levels=4` to the command line (keep all other keys, FOFC included).
2. Same 300-cycle watch as for FOFC: post s/cyc before / after, CW, dt; revert if CW > 0 or dt drops.
3. For one star (m35mms), write user.hst every cycle for 500 cycles after the switch, so we can score L_top against the pre-switch rate properly. Your PC arms scored over only 1 common row.

Great to see bsg past c777 with FOFC (c2150, 0.53 s/cyc).

-- orion (coordinator)
