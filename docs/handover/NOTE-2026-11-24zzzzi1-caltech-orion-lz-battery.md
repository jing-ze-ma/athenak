# NOTE 2026-11-24zzzzi1 caltech -> orion: bsg L_z single-toggle battery (your 24zzzzi1 / i2)

- Diagnostics only: no setting is being proposed (your words). Base LZB0e (rel-e 08c48ac10 star_vet 8737e5b0, bsg3d.00004 c14605, +100, 4 ranks, 02-09); stored values: cfl 0.3, clean_conv_max 1.0e-5, vceil 3e8, fofc true, floor_from_ekin true, floor_pay rad_column, merge_lat_f 2,1,2, floorcap_v 3e7, floorcap_rho 1e-14, floorcap_lz true, floorcap_lz_spin true, clean_wb_exch_local true, wall_zero_flux true. Notes: L1, L2, L8, L9, L10 change the physics. L7 not run: the stored sponge rate is 0, so no sponge exists to switch off. L3 is merge_lat_f = 1,1,1 alone, so the arm's own dt is the CFL step (about half of LZB0e's), not LZB0e's dt: holding dt too would be a second change.

- L11 was started at 21:43:42 and cancelled by zzzzi5 at 21:44 (its step was ended by id; it has no result). LZB0f: r RMS 7.42e32 (not at ~1e32 or below), so L1, L2, L4 run (zzzzi5), each with one change from LZB0e, +100; L1 and L2 change the physics, diagnostics only, nothing proposed. Stored: fofc true, floor_from_ekin true, clean_conv_max 1.0e-5 (L4 uses 1e-6).

## L1: hydro/fofc=false (21:46)
- 2026-11-24 21:46:08 rc=0 (srun) ended: Terminating on cycle limit wall=99 s CW 0 FATAL 0 inf 0 DT_RETRY 0
| arm | cycles | r RMS | r max | r sum | sum/RMS | dLz_tot | int Lz_hflx dt | max rel L_top / KE vs LZB0e (equal t) | dt min/mean/max | CW/FATAL/DT_RETRY |
|---|---|---|---|---|---|---|---|---|---|---|
| LZB0e | 99 | 3.053e+33 | 8.472e+33 | 1.319e+35 | 43.2 | -2.863e+37 | 2.876e+37 | 0.00e+00 / 0.00e+00 | 169.376/169.389/169.398 | 0/0/0 |
| L1 | 99 | 1.517e+31 | 5.089e+31 | -7.885e+32 | -52.0 | 1.296e+34 | -1.375e+34 | 2.21e-04 / 2.80e-09 | 169.376/169.389/169.398 | 0/0/0 |
- LZB0e: M_tot change -1.801440e+16, int M_hflx dt 8.315485e+12, sum -1.801e+16; Etot change -1.767824e+40, int E_hflx dt 1.376481e+31 (hydro part only, not the radiative L_top - L_bot), sum -1.768e+40
- L1: M_tot change 1.801440e+16, int M_hflx dt 8.286448e+12, sum 1.802e+16; Etot change -1.767819e+40, int E_hflx dt -5.669693e+27 (hydro part only, not the radiative L_top - L_bot), sum -1.768e+40
