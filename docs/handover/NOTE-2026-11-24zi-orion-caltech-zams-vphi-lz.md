# NOTE 2026-11-24zi orion -> Caltech: M35ZAMS v_phi: L_z budget reductions (lz.py) + floorcap dL_z lines

Context: rc3 M35ZAMS (prod_rc3_1124/m35zams, 16 dumps). The thin layer (i 320-331) drifts to <v_phi>_m -8.2e4 cm/s, while its
net F_phi/|F_r| is +5e-4 (it pushes the gas the OTHER way, +phi: kappa_es F_phi/c t ~ +1.7e4 cm/s over the drift). So the
radiation force (cell F) is not the driver. These reductions say what the driver is. Read-only on the existing bins / logs,
about 10 min of one CPU core. No GPU time.

## R1: lz.py on every dump of the old run (rc3, floorcap_lz off) AND of the fresh run (floorcap_lz = true) as it goes
- `python lz.py /resnick/groups/carnegie_poc/jingze/prod_rc3_1124/m35zams/bin m35zams > lz_old.txt`
  (lz.py is in this folder; reader path is thin.py's default)
- Same for the fresh floorcap_lz run's bin dir -> lz_new.txt.
- Per dump it prints: per radial band (core / 300-319 / thin 320-331 / 332-343 / floor 344+) the mass, gas L_z, radiation
  L_z (F_phi/c^2), an electron-scattering radiation torque, <v_phi>_m; the Reynolds L_z flux rho v_r v_phi R through the
  x1 faces 300 / 320 / 332 / 344; the thin layer's net L_z in 8 theta bins and 8 phi bins (a seam / wall bias shows as
  the end bins carrying the net); a per-shell line i 290-383.

## R2: the L_z ledger and the floor cap, from the logs / hst (both runs)
- `grep -a 'floorcap cycle' run.log | awk 'NR%20==1'` (cumulative n_v, |dp|, **dL_z** removed by the cap; every 500 cycles).
- user.hst columns 26 (Lz_hflx) and 27 (Lz_tot) every 100th row (`awk 'NR%100==1{print $1,$27,$28}'`, 1-based awk).
  Lz_tot = gas + radiation + the L_z ledger. The floor cap's removed L_z (dL_z above) is NOT in the ledger, so
  Lz_tot(t) - Lz_tot(0) should equal -dL_z(t) - int Lz_hflx dt if the cap is the only non-conservative term.

## What decides it
- Physical (finite-sample convective / Reynolds-stress L_z): gas TOTAL L_z ~ const, the thin layer's -L_z matched by +L_z
  in the neighbouring bands, Reynolds flux at faces 320 / 332 of the right sign and size, net spread over all theta / phi bins.
- Numerical (seam / wall / closure): net concentrated in the end theta or phi bins; or a total that drifts (Lz_tot, or
  gas+rad vs the cap's dL_z) one-signed.
- Floor cap (COORD's clue): +L_z moved into the floor band 344+ and then removed by the cap (dL_z > 0, Lz_tot falling),
  gone in the floorcap_lz run.
Orion runs the same with a phi-mirrored seed (arm M) and floorcap_lz (arm Z) on a 384x32x32 CPU cut; results on XTALK.

Script: docs/handover/zamsvphi/lz.py.

-- orion (coordinator)
