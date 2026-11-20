# NOTE 2026-11-20dp orion -> caltech (USER): production EMA 128 named; v20.3 gate at the FIRST free node

1. **USER: production HG/RSG = `rad_m1/vet_rotate_avg_ema=128` with `vet_rotate_avg_ema_change_form=step`, NAMED explicitly.** The gate restarts (rsg.01295 / hegiant.00139) predate this, which is why their dumps show only the plain ring.
   - **From v20.3 on, every gate or long arm names both keys** (the v20.3 keys files will carry them).
   - Updated now: `arms_v203.txt` (V203L_rsg, V203L_hg and V203LP_rsg gain both keys) and `arms_emawb.txt` (all four EW arms gain change_form=step; EMA 128 was already named). EW_on vs EW_off, both with EMA 128, is the production baseline.
   - V203_rsg / V203_hg stay without EMA, so they remain the bitwise check against V202F.
2. **v20.3 gate: squeeze it in at the first free node**, ahead of the batch. It is 5 short arms and the only thing open for the v20.3 cut. Order: V203L_rsg, V203L_hg, V203LP_rsg (the key-on rows) first, then V203_rsg / V203_hg.

-- orion (for the USER)
