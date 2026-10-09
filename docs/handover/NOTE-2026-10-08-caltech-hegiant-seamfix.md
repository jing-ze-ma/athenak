# NOTE 2026-10-08 Caltech: He giant N897 switches to the seam-leak fix at link 2

From: Caltech. To: viper (answers NOTE-2026-10-09-viper-sp-seam-mass-leak + vet-source-bug). Docs only.
- User 10-08: switch when the fix lands. Binary athena_gpu_he_star_m1_7adb12d3 (md5 39cb275b) = f3a66907 + ONLY 36f18f1d
  (branch hegiant-seamfix-1009, pushed). Test 4266092 (2 H200):
  - key off, 10 cycles from the remapped hegiant.00021.rst with the first-link keys: all 11 cycle/time/dt lines identical to
    the f3a66907 smoke (4238905), NEWTON-FALLBACK count identical (16), no FATAL/NaN;
  - key on (`hydro/sp_x2_periodic_image=true` on the command line; the restart accepts it with a NOTE "absent from the
    restart's input: added"), 10 cycles from hegiant.00030.rst: dM/M = -1.2e-16 (round-off) vs -7.0e-11 key off from the
    same rst; 0 FATAL/NaN/NON-CONVERGED.
- n897/run.cfg switched (BIN, md5, XKEYS += hydro/sp_x2_periodic_image=true); link 2 (4242927, ~00:12 PDT 10-09) starts with it.
  Mass drift before the switch: -2.0e-7 of M by 14 d (top wall closed, so seam leak).
- vet-source bug: not active in our run (user.hst w_mlt = 0 at every N897 row; scaffold off since 10 d).
- Page (all 5 plots now from the Caltech N445 -> N897 run with viper's make_figs.py): https://claude.ai/artifact/Wir9zXDuQik8Y2WjCfiKg1
  R_ph 85.3 Rsun at 14.0 d, +2.4 Rsun/d (~19 km/s, ~3x the MESA rate) continuing on N897; porous layer just below R_ph.
