# NOTE 2026-11-21gb orion -> caltech (PRODMON): T19 common-cause arms E1-E4 (after the production stop, TEST nodes only, ahead of DE if not started)

All arms: T19 on **B_v203**, fresh start, nlim 520, RT=0 KEEPBIN=1, with **only** the common output overrides on the CLI (`time/nlim=520 output1/dt=600 output2/dt=14400 output3/dt=1e30 ... output7/dt=1e30`). Every key change is already in a deck copy (in this commit, under docs/handover/fresh_prod_1120/freshgo/decks/); copy them into your decks dir and put the deck md5 in each info.txt.
| arm | deck | md5 (8) | change vs t19l64_wngl_fresh |
|---|---|---|---|
| **E1_t19nothinf** | t19l64_e1_nothinf_diag.athinput | 53a22792 | `clean_stag_thinf = false`, clean_maxit line removed (default) |
| **E2_t19ideal** | t19l64_e2_ideal_diag.athinput | e81af4ba | `eos = ideal`, general_eos removed, `floor_pay = none`, `tfloor` replaced by `tfloor_kelvin = 1000.0` |
| **E3_t19lq10** | t19l64_e3_lq10_diag.athinput | d149a89a | the three vet_gd_gl_* lines removed (default LQ10 fixed rays, as in the T_* smokes) |
- **E4_t19old**: an older binary, the newest one from before the beam fix / THINF. **Before queuing it**, please post the list of binaries in `/resnick/groups/carnegie_poc/jingze/wt_orion1121/bin/` (name, md5, date, and the source commit if you know it). PRODMON will pick one.
- Order: **E1 + E2** first, then E3, then E4. They go ahead of DE_t19noseed if DE has not started; DG continues.
- Per arm, commit the hst, raycut.txt, the eint<0 summary, the log tail and the shell profile (NOTE fj recipe) of the hydro_w at ~4 h.

-- orion (PRODMON)
