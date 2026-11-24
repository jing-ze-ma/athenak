# NOTE 2026-11-24h orion -> Caltech: IA12.3-rc3 = caltech-ia123rc3-af836f5a3 (rc2 + the speed keys' accuracy fix); decks with clos1_budget + bud_maxnorm ON

- rc3 = rc2 5b493acfa + one 6-line commit: the budget tolerance floor uses clean_tol_nl when a budget key is named. That makes clean_clos1_budget / clean_bud_maxnorm accurate on the radiation-wave test (4.4e-3 -> 2.4e-5). It is bitwise to rc2 on the star (the floor never binds there).
- Gates: CPU 20/20 key off (56 dirs bitwise vs rc2) and 20/20 with the keys; MI300A 10/10.
- Decks: `docs/handover/ia123rc3_decks/` (seam keys + floor_from_ekin + clos1_budget + bud_maxnorm ON, vel_settled OFF).

**What to do with it:**
1. If you haven't built rc2 yet: build rc3 instead (incrementally, 1 file changes) and run everything from 24e / 24g on it. The 24g PC arms then use the rc3 decks as K.
2. If rc2 is already built and the 24e arms are running: let them finish. Build rc3 for the PC arms (24g) and for the production starts.
3. Production starts (bsg / t9l62, USER-approved): on rc3 with the rc3 decks, once the 200-cycle arms are clean.
4. m35mms production: at its next rst, move from IA12.2 B to rc3 with the rc3 m35mms deck (tile / notile as now). That gives the speed keys (your E4 arm: about -15..-20 %) + the PCR kernel (-4.7 %). Run 200 cycles first and check CW 0 + L_top within the D1h band vs IA12.2 B.

-- orion (coordinator)
