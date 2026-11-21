# NOTE 2026-11-21ha orion -> caltech: REDO the fresh-production smokes on CUDA 12.9 binaries, incl. rebalance (USER)

USER: "redo those smokes, especially the T19L6.4 ones; can we use rebalance or not". Every earlier smoke and the rebalance
verdicts ran on CUDA 12.2 binaries (B_v203, B_icbal 4a7f72334) and are void.

Builds (CUDA 12.9 + g++ 13, own ccache_cu129, whole node -j32, NOTE gy), smoke-gate each (tiny BSG pass 1 = 4.020e-06):
- B_v203n (have it, md5 8d3f469f).
- B_icbaln = caltech-runs/icbal 4a7f72334 (problem/he_ic_balance_aref = file | dpr), rebuilt with CUDA 12.9.

Arms (armg.sh, freshgo/decks, 1 node each, NR 4; m1prof + m1_vet at the end; hst + raycut; bins deleted after):
T19L6.4 first, 520 cycles each (covers the old ~1.4e4 s collapse window); a only, except b-leg on R0:
  R0_t19   B_v203n  t19l64_wngl_smoke  520  (= production recipe, no rebalance; b = 260 + rst + 260 bitwise)
  R1_t19   B_v203n  t19l64_wngl_smoke  520  problem/he_ic_balance=true                      (analytic, file a_ref)
  R2_t19   B_icbaln t19l64_wngl_smoke  520  problem/he_ic_balance=true problem/he_ic_balance_aref=dpr   (discrete force)
Then, 120 cycles each (a only; b-leg on the WN arm of each star):
  R0_t9    B_v203n  t9l62_wngl_smoke   120   |  R2_t9   B_icbaln t9l62_wngl_smoke  120  he_ic_balance=true aref=dpr
  R0_m35   B_v203n  m35mms_wngl_smoke  120   |  R2_m35  B_icbaln m35mms_wngl_smoke 120  he_ic_balance=true aref=dpr
  R2_bsg   B_icbaln bsg_wngl_smoke     120  he_ic_balance=true aref=dpr   (WN BSG = BV_bsg already running)
  R0_t19c  B_v203n  t19l64_corona_wngl_smoke 120  |  R0_t9c  B_v203n  t9l62_corona_wngl_smoke 120
Keys go IN a deck copy, not the CLI (fresh-start FATAL otherwise). Report per arm: KE_int, L_bot/L_in, L_top/L_in,
dt every 20 cycles, CW/FATAL, clean_init lines, inval fraction above the cut. Results to caltech_results/R/.
Order after the BV arms (gx): R0/R1/R2_t19 first. Test/former-production nodes only; no production launch.

-- Orion coordinator
