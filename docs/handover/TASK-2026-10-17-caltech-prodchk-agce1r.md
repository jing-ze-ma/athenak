# TASK raven -> Caltech: READ-ONLY surface/background check of the AG Car production agcE1r (analysis only)

From the Raven coordinator (worker S3GRID), USER approved 10-17. Post-processing only: NO GPU job, NO change to agcE1r (do not
cancel, hold, edit run.cfg/deck, or touch its files; read them only). CPU on a login node or a small CPU/debug job is enough.

## Why
Raven diagnosed (rol2u_safety_1015 NOTE, section BGDIAG) that the half-resolution S3 test of AG Car A (268x32x32, same E1 deck
and poly grid as agcE1r, half the radial cells) is out of thermal balance after the ramp (L_top/L_in 0.51 -> 0.74). Its envelope
inflates, and the surface (shell <rho> = 1e-12) moves from 1.12 to 1.56 R_ph, beyond the fine radial zone (fine to ~1.15 R_ph) and
beyond wb_rmax (1.1 R_ph). From ~2.0e8 s the surface is a 1-2 cell density cliff, the hydro leaves eint < 0, and the radiation E
floor (implicit_pos_floor) CREATES energy (up to 1.3e37 erg/s = 6 % L_top). Backward-Euler share rises to ~0.94, gas is ejected
at vceil (1000 km/s), and the sponge stops acting (rho at the top >> he_sponge_dmax). agcE1r has the same coarse zone above 1.15
R_ph (at 2x the S3 radial resolution) and its L_top (1.35e38) is also ~0.48 of L_in (2.84e38). We need to know if and when agcE1r
enters this regime.

## Files (docs/handover/prodchk-1017/, check MD5SUMS)
- prodchk.py: the diagnostic (python3 + numpy only); bin_convert.py: the AthenaK reader (copy of vis/python).
- S3_LB_ref_{hst,dump,sum}.txt: the same script's output on the Raven S3 BASE arm (3e7 -> 4.68e8 s) for comparison.

## Run (in the agcE1r chain dir, e.g. rol2u_1015/e1/chain; all links)
```
cd docs/handover/prodchk-1017
python3 prodchk.py --run <CHAINDIR> --stdout '<CHAINDIR>/<link stdout logs glob, all links incl. agcE1r link 1-3>' \
    --dstride 5e6 --win 1e7 --out caltech_agce1r
```
- <CHAINDIR> must hold agcar3d.user.hst and agcar3d.log (event log); the dumps are read from <CHAINDIR>/bin/*.hydro_w.*.bin
  (else pass --bins '<glob>'). If the hst/log of the earlier chain (agcE1s, t < 2.6e7) are elsewhere, that is fine: start at 2.6e7.
- --stdout may be given several times (one glob each). Gzipped logs are read too.
- It reads one hydro_w dump per 5e6 s (~50-60 dumps of 536x128x128); if that is too slow, use --dstride 1e7.
  For a quick first answer, --nodumps gives the hst/log/stdout tables in a minute.

## Report back (NOTE-2026-10-17-caltech-prodchk-agce1r.md)
- caltech_agce1r_{hst,dump,sum}.txt (contents);
- the current state of agcE1r: running link/job, current t, binary, any change of keys since 10-15 (expect none);
- wall time of the analysis. Nothing else. Please flag at once if r_surf > 1.25, BE share > 0.7, the Bres slope
  > 1e-3 L_top, or the frac(<=2 steps) > 0.5 in any window. Raven decides nothing about production: the USER does.
