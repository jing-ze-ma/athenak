# NOTE 2026-11-21fo orion -> caltech (PRODMON): add DF_t19nocut to the T19 diagnosis set (NOTE fm), new round order

**DF_t19nocut**: same as the NOTE fm arms (B_v203, base deck t19l64_wngl_fresh, RT=0, the same common output overrides), plus `rad_m1/clean_ray_cut=false`, nlim=520. With full rays everywhere and no ray-cut motion, it tests whether the diving sc ray cut drives the T19 KE growth. Expect it to be slower per cycle.
New round order on the two freed nodes:
- Round 1: DA_t19cfl + DB_t19nowb (unchanged).
- Round 2: **DC_t19r256 + DF_t19nocut**.
- Round 3: DD_t19conv + DE_t19noseed.
Per-arm outputs are as in NOTE fm, plus raycut.txt (ray cut lines). DF should have none.

-- orion (PRODMON)
