# NOTE 2026-11-24zzzzd8 Caltech -> orion: cs arm L, where is the missing KE (zzzzd8)

Read-only (Bp and L hydro_w c150/300/450/600, hst). Files: caltech_results/CS_KE/ke.txt (all bins, cycles 150-600), ke.py. Full sphere, erg. The dump KE_r and lateral KE reproduce hydro.hst 1-KE and 2+3-KE to all printed digits (lateral with v components contravariant on the unit tangents: v2^2 + v3^2 + 2 v2 v3 e1.e2). **48-image mean: built and used** (the 48 signed permutations map the equiangular N = 32 cell centres onto themselves; bijection checked to 1e-10; vectors rotated), so no 6-panel substitute. v_sym = (1/48) sum R^T v(Rx). Cross term v_sym.(v - v_sym) < 1e-5 of the total everywhere. Because v_sym contains the l = 0 radial flow, "cube" below = v_sym minus the shell-mean v_r.

**1. c600, total KE and its parts**
| r [Rsun] | Bp total | L total | L - Bp | Bp l = 0 | Bp cube (radial) | L cube (radial) | Bp rest (not symmetric) | L rest |
|---|---|---|---|---|---|---|---|---|
| 0-100 | 8.034e39 | 1.393e38 | -98.3 % | 7.3e35 | 7.981e39 (6.237e39) | 1.385e38 (1.241e38) | 5.26e37 | 4.3e34 |
| 100-200 | 9.241e38 | 1.235e37 | -98.7 % | 4.3e34 | 9.188e38 (6.763e38) | 1.234e37 | 5.2e36 | 7e32 |
| 200-400 | 3.529e38 | 1.247e36 | -99.6 % | 2.0e34 | 3.513e38 (2.959e37) | 1.231e36 | 1.5e36 | 8e30 |
| 400-600 | 4.265e37 | 8.28e34 | -99.8 % | 1.6e34 | 4.242e37 | 6.67e34 | 2.1e35 | 2e28 |
| 600-800 | 3.82e36 | 2.73e35 | -92.9 % | 2.68e35 | 3.54e36 | 4.8e33 | 1.9e34 | 1e31 |
| 800-850 | 5.7614e40 | 5.7552e40 | -0.1 % | 5.7612e40 | 8.68e34 | 6.14e34 | 2.33e36 | 2.32e36 |
| 850-900 | 1.1316e38 | 1.1314e38 | -0.0 % | 1.1315e38 | 1.63e33 | 1.61e33 | 1.00e34 | 1.00e34 |
| all | 6.7198e40 | 5.7931e40 | -13.8 % | 5.7839e40 | 9.297e39 (6.943e39) | 1.522e38 | 6.19e37 | 2.38e36 |

- Bp - L total: 9.267e39 erg; l = 0 part 6.2e37 (0.1 % of l = 0), the cube part (l = 0 removed) 9.145e39 = **98.7 %** of it, non-symmetric rest 6.0e37 (0.6 %). Of the cube part, 86 % is in 0-100 Rsun (7.84e39, radial 6.11e39 + lateral 1.73e39), 10 % in 100-200, 4 % in 200-400.
- At c600 the Bp - L KE difference is 98.7 % in the cube-symmetric part, mostly below 100 Rsun. This is a decomposition of the two dumps; it does not say why Bp has it. (hst KE_int, the rule's column, is -15.43 %; its interior region is a subset, 6.010e40 vs 5.083e40; the dump totals are 13.8 %.)
- Earlier cycles, total gap Bp - L: c150 1.38e39 (6.1 % of Bp), c300 2.93e39 (6.3 %), c450 5.46e39 (9.6 %). Cube part 1.38e39, 2.92e39, 5.46e39 in Bp against 1.4e37, 4.7e37, 9.4e37 in L.

**2. KE_int(t), hst (row = cycle, dt 1000 s)**: gap (L - Bp)/Bp: c100 -10.2 %, c200 -5.4, c300 -6.4, c400 -8.6, c500 -11.5, c600 -15.4 %. Absolute gap (Bp - L) 1.08, 1.80, 2.93, 4.49, 6.58, 9.27e39 erg: log-log slope 1.29 fitted over c100-600; local slopes 1.49 (c200-600), 1.79 (c400-600), 1.88 (c500-600), so it is rising towards 2 and is neither constant nor t^4 in this range.

**3. Others at matched cycles (L relative to Bp)**: L_bot +6.9e-4 to +7.5e-4 (c100-600); L_mid -7e-8; L_top -7.9e-5 to -5.9e-6; M_tot +3.4e-8 and Etot +3.3e-8 (constant); Lz_tot 2.523e38 (Bp) vs 2.526e38 (L) at c100, 2.522e38 vs 2.526e38 at c600 (+0.1 to +0.2 %, present from the first row compared); KEr_int -2.8 to -11.9 %. |v| max equal (1.727e6 at 850-900 Rsun); below 100 Rsun 8.1e3 (Bp) and 1.06e3 cm/s (L) at c600.

**4. Photosphere 800-850 Rsun**: total KE -0.1 % (5.7614e40 / 5.7552e40, almost all l = 0), lateral -5.4 % (1.40e35 / 1.32e35), cube part -29 % (8.68e34 / 6.14e34, 1.5e-6 of the KE there). m4 there -24 % is your figure from the table, not remeasured.

-- Caltech
