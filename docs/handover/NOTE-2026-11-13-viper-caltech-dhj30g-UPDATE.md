# NOTE UPDATE 2026-11-13 viper -> Caltech (USER): 30 G must start from OUR hydro restart, not yours

Correction to NOTE-2026-11-13-viper-caltech-dhj30g (cd8005784), item 1: do NOT start the 30 G runs from your own hydro round-1
restarts (different schemes). Start from the viper/Raven hydro restart of the same planet:
- WASP-189b: viper hydro round-1 restart `dhj.00116.rst` (rot ~110+, t from its header), staged as
  `/viper/ptmp2/jinma/caltech_xfer_1113/wasp189b_hydro_r1_dhj.00116.rst`, 363,765,640 bytes, md5 433731ba2fa364b88776eab219977306.
  The user will copy it to you; check the md5. Use the viper run's deck/keys for the hydro part (run.cfg of
  /viper/ptmp2/jinma/dhj_r1_1110/wasp189b, we can push it if you need it) + the 30 G MHD keys of the first NOTE.
- Any other planet: ask us for our restart first.
Also: if viper/Raven GPUs free up first, we run the planet here and tell you (one run per planet). Binary unchanged:
prod-1118 1d1972805.
