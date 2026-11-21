# NOTE 2026-11-21gk orion -> caltech (USER via coordinator): test on BSG instead of T19 -- re-prioritise now

The USER wants the tests run on BSG rather than the very luminous T19.
1. **Now, all in parallel, on free nodes: BB_bsgcutr, BD_bsgthinf0, BC_bsgtops** (NOTE ge decks).
2. **Let E1 / E2 / E3 finish**, but **do NOT start E4_t19old, DGL_t19cutr or any new T19 arm**. DGL only runs if a node would otherwise idle.
3. Then queue these **BSG equivalents**, each on its own node as nodes free. All are B_v203 unless stated otherwise, fresh start, tlim 1.12e6 s (13 d), outputs as in the deck (m1 every 0.5 d). Each deck is committed here, with one change from bsg_wngl_fresh plus tlim:
| arm | deck | md5 (8) | change |
|---|---|---|---|
| **BE_bsgideal** | bsg_diag_ideal.athinput | 45622ec5 | eos = ideal, general_eos removed, floor_pay = none, tfloor -> tfloor_kelvin = 1000.0 (as E2) |
| **BF_bsglq10** | bsg_diag_lq10.athinput | 17d0dc00 | vet_gd_gl_* removed (default LQ10 fixed rays) |
| **BG_bsgold** | bsg_diag_base.athinput | a34a960f | no key change (tlim only); run on **athena_v202_f4cefb181** (md5 54c3df35). If v202 FATALs on an unknown key, report the key and stop; do not edit the deck. |
| **BH_bsgnowb** | bsg_diag_nowb.athinput | 4233fa86 | WB off as in the old runs: wellbalance_dynamic = false, force_reference = none, he_gravity = pointmass, force_reference_work removed |
4. For **all seven BSG arms** (BB, BD, BC, BE, BF, BG, BH): commit the hst, raycut / oldvec, log tail, plus the **m1prof (NOTE fy recipe; BSG L_in 5.58e38) of the m1 bins at 0.5 d and at the end**. Refresh the hst about every 20 min. Delete each arm's own rst at the end. Put the deck and binary md5 in info.txt.
5. The v20.4 CM* / V204* arms come last.

-- orion (coordinator, for the USER)
