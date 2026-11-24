# NOTE 2026-11-24zd orion -> Caltech: RSG R_ph continuity: run rph.py on 08362-08365

Thanks for the c100 hst (dead449ad). Orion's check of it:
- L_top/L_in is continuous across the old -> clean hand-over:
  - new: 0.9256 +- 0.0017 over 155 cycles;
  - old run, last 1e6 s: 0.9239 +- 0.0009;
  - old-binary continuation from the same rst: 0.925.
- Etot follows L_in - L_top, M_tot is flat and Lz drifts 9e-13.
- L_bot steps 1.006 -> 1.049 L_in. That is a diagnostic: the staggered face F includes the wall's 4/3 v E, about 5 % of L_in, measured on the old 08361 dump. A fix is coming in the hst.
- So keep RSG running.

Next, the radius. Please run `rph.py` (attached in this folder) on the restart dumps 08362, 08363, 08364, 08365. Run it on a login node or in a small CPU job: ~30 s and ~6 GB of memory per dump.

```
T=<your RSG <problem>/he_opac_table, rosseland_mesakap_4c0e8b64d26a.txt (md5 fa0f4f420022bd6c0fec4eb1803619b0)>
V=<athenak checkout>/vis/python        # bin_convert.py
D=/resnick/groups/carnegie_poc/jingze/rsg_rc3_1124
python3 rph.py $V $T $D/rsg.user.hst $D/bin/rsg.hydro_w.0836{2,3,4,5}.bin > rph_rc3.txt
```
- The m1 file is found from the hydro_w name.
- L_top for T_eff is read from user.hst, from the row nearest each dump time.

Post `rph_rc3.txt` to caltech_results/RSG_RC3_1124/. It is 4 lines, one per dump, each with:
- median / mean / p10 / p90 of the per-column tau_R = 2/3 radius;
- the shell-averaged-opacity tau = 2/3 and tau = 1 radii;
- R_Teff;
- T_eff.

Orion ran the same script on the old run's dumps 08354-08361 and on the old-binary continuation from the same rst (viper, oldprod1121/rsg, 08362-08366). Compare by time, not by dump number: our continuation's 08362 is at t 7.22477e8, while your 08362 is the restart state at t 7.22398e8.

Reference (old scheme), radii in cm:

| run | dump | t [s] | R23 median | R23 mean | R(tau=1) mean | Rsh23 | R_Teff | T_eff [K] |
|---|---|---|---|---|---|---|---|---|
| old prod | 08359 | 7.222180e+08 | 7.10983e+13 | 7.09487e+13 | 7.07430e+13 | 7.78327e+13 | 7.68555e+13 | 3642.9 |
| old prod | 08360 | 7.223040e+08 | 7.10878e+13 | 7.09385e+13 | 7.07338e+13 | 7.78213e+13 | 7.68672e+13 | 3643.5 |
| old prod | 08361 | 7.223910e+08 | 7.10809e+13 | 7.09267e+13 | 7.07227e+13 | 7.78100e+13 | 7.68780e+13 | 3644.2 |
| old cont. | 08362 | 7.224770e+08 | 7.10749e+13 | 7.09147e+13 | 7.07119e+13 | 7.78003e+13 | 7.68864e+13 | 3644.9 |
| old cont. | 08363 | 7.225640e+08 | 7.10718e+13 | 7.09062e+13 | 7.07041e+13 | 7.77889e+13 | 7.68954e+13 | 3643.4 |
| old cont. | 08364 | 7.226500e+08 | 7.10653e+13 | 7.08943e+13 | 7.06929e+13 | 7.77814e+13 | 7.69016e+13 | 3642.5 |
| old cont. | 08365 | 7.227360e+08 | 7.10616e+13 | 7.08823e+13 | 7.06822e+13 | 7.77723e+13 | 7.69017e+13 | 3644.3 |

- Old trend: R23 mean falls by about 1.2e10 cm per dump (-0.017 %, ~1 d), smoothly.
- A continuous hand-over keeps R23 mean within ~0.05 % of the continuation at the same t.
- For scale, the 11-21 clean switch (bf4f917d4) moved the median by +1.0 % within 1 d.
- If yours is off by more than ~0.1 %, post the result and keep the dumps.

-- orion (coordinator)

Script + old-run reference table: docs/handover/rsg_rph/.

-- orion (coordinator)
