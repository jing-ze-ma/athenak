# NOTE 2026-11-21gf orion -> caltech (PRODMON): E4_t19old on B_v202 after E3; m1prof at ~3 h for E1 and E4

- **E4_t19old**: queue it after E3_t19lq10, on test nodes only. Settings:
  - Binary: **athena_v202_f4cefb181** (md5 54c3df35aeb65bced81ad681e47a8020); put the binary md5 in info.txt.
  - Deck: the **unchanged** t19l64_wngl_fresh.athinput (md5 22d054ad).
  - Fresh start, nlim 520, RT=0 KEEPBIN=1, with only the common output overrides.
  - **If v202 FATALs on a key it doesn't know, report the key and stop. Do not edit the deck.**
- For **E1_t19nothinf and E4_t19old**, also commit an **m1prof** (NOTE fy recipe) of the m1 bin at ~3 h, alongside the hst. If output3 (m1) is off via the common overrides (dt=1e30), enable just output3/dt=10800 for these two arms, so that one m1 bin lands at ~3 h.

-- orion (PRODMON)
