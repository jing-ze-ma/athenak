# NOTE 2026-11-26j caltech -> orion: axis share against radius, time order, the radial motion; older m35mms runs (26i checks 3 and 2, 26j)

Dumps and hst only, read-only; no run, no twin. Band = k256 16-40 (0.0625-0.156 cycles per cell, 16-6.4 cells per wavelength, the same fraction of the grid in all runs). "x" = axis share (|kx| <= 1 or |ky| <= 1) of the band power divided by its share of the band's modes (0.068 m35mms, bsg; 0.051 t9; 0.139 m35zams, which has fewer modes). Dumps averaged over the last 6 of 110 (m35mms, t to 7.92e5), 48 (t9l62, 2.04e6), 50 (t9mlt, 1.60e6), 55 (m35zams, 7.79e5), 30 (bsg, 4.92e6). T = T_gas / <T_gas> - 1.

## Measured

### 1 (26j): x against radius, from the inner wall up (median of x for v_r / T)
| star | 1.00-1.10 r_in | 1.10-1.20 | 1.20-1.30 | 1.30-1.40 | 1.4-2.0 | convection zone and above | crossing of x = 2 (outward) |
|---|---|---|---|---|---|---|---|
| m35mms | 11.0-12.4 / 11.4-12.3 | 12.5-7.6 / 12.0-6.3 | 1.1-1.2 / 1.3-1.4 | 1.3 / 1.4 | 1.2 | 1.2-1.5 to 1.03 R_ph (convection base N^2 < 0 at 0.924 R_ph = 1.308 r_in) | 1.207 r_in (0.852 R_ph) |
| t9l62 | 8.2-8.4 / 8.4 | 8.6 / 8.5 | 8.6-7.3 / 8.4-5.5 | 1.5 / 1.7 | 1.2 / 1.3 | 1.1-1.4 up to 5.7 R_ph,IC (base at 1.8 r_in) | 1.315 r_in |
| t9mlt | 8.3 / 8.5 | 8.5 / 8.4 | 8.2-5.7 / 8.1-4.3 | 1.8 / 1.7 | 1.1-1.2 / 1.2-1.4 | 1.2-1.5 up to 1.0 R_ph,IC (base at ~1.74 r_in), above the photosphere 1.5-3.3 (T up to 4.5) | 1.324 r_in |
| m35zams | 3.5-5.7 / 0.9-1.0 | 1.4 / 1.7 (1.10-1.15) | - | - | - | 1.4 / 1.8 up to 1.01 R_ph | 1.111 r_in (v only, with 4 earlier short crossings) |
| bsg | 2.8-4.7 / 2.9-3.3 | 5.1-5.2 / 3.5-3.9 | 5.3-4.7 / 3.3-2.3 | 3.6 / 2.2 | 1.9 (1.4-1.6), 1.3 | 1.1-1.7 | 1.481 r_in |
- **Coming from below or above:** x is large from the wall outward in all five stars' stable layers and drops to 1.1-1.5 (the counting value with the noise of the few modes) below the base of the convection zone. It does not fade gradually from the convection zone: the excess ends 0.1-0.5 r_in below that base (m35mms 1.207 against 1.308 r_in; t9l62 1.315 against 1.8; t9mlt 1.324 against 1.74; m35zams 1.111 against 1.138; bsg 1.48 against about 1.9 r_in).
- **Convection zones:** x is 1.1-1.5 in all five, v_r and T: the band is isotropic there at these scales (also m35zams 1.4 and bsg 1.2-1.7). The layers above the photosphere have x up to 3 in t9mlt (T up to 4.5).
- **Amplitude:** the excess is not the same as the ring: m35zams (x 3.5-5.7 of v_r) has it at 1e-8 to 2e-3 km/s in the stable layer, with T isotropic below 1.09 r_in, and bsg at 1e-3 to 4e-3 km/s (x 5, T 3-4): the axis-preferring motion exists there at tiny amplitude; only m35mms, t9l62 and t9mlt have grown it.

### 2 (26j): time order, band rms at the base of the convection zone against the ring (first time the band rms exceeds the level, sampled every 4 dumps m35mms 7.2e3 s, t9l62 4.3e4, t9mlt 3.2e4)
| star | base (r) | 1e-3 km/s | 1e-2 | 1e-1 | 1 | ring (r) | 1e-3 km/s | 1e-2 | 1e-1 | 1 |
|---|---|---|---|---|---|---|---|---|---|---|
| m35mms | 8.83e11 | 2.9e4 s | 3.6e4 | 4.3e4 | 4.3e4 | 7.07e11 | 2.30e5 | 3.02e5 | 4.75e5 | 6.98e5 |
| t9l62 | 3.05e12 | 4.00e5 | 4.0e5 | 4.0e5 | 4.0e5 | 1.97e12 | 4.43e5 | 4.86e5 | 9.18e5 | 1.61e6 |
| t9mlt | 2.95e12 | 1.1e4 | 4.3e4 | 1.7e5 | 2.1e5 | 1.97e12 | 3.67e5 | 6.26e5 | 9.83e5 | - |
- The band at the base of the convection zone is always earlier, by 2.0e5 s (m35mms, 1e-3 km/s), 4e4 s (t9l62, one or two dump spacings) and 3.6e5 s (t9mlt). Detrended correlation of ln base with ln ring: m35mms r = -0.17 at lag 0 (95 % level 0.21), t9mlt -0.17 (0.31; 0.52 at lag -5 dumps, i.e. the base earlier), t9l62 0.84 at lag 0 (both rise together).
- First time above 1e-3 km/s between the ring and the base: m35mms ring 1.047 r_in 2.30e5 s, 1.08 2.45e5, 1.11 3.02e5, 1.145 6.8e5, 1.18 never, 1.21 5.8e4, 1.24 5.0e4, 1.28 4.3e4, 1.31 2.9e4: the ring is not reached by a front moving down from the base (the gap at 1.18 has no 1e-3 km/s at any time, and the ring starts after the base but before its neighbours above the gap). t9mlt: earlier radii first from the base down to 1.16 r_in (1.1e4 to 3.7e5 s, monotonic), t9l62 all within one dump spacing (4.0-4.4e5).
- Band minimum between the ring and the convection (m35mms): 2.4 km/s at 0.74, 9.8e-3 at 0.80, 5.2e-4 at 0.82, 6.7e-5 at 0.84, 1.8e-3 at 0.86, 0.47 at 0.88, 8.6 at 0.92 R_ph; axis share 0.82, 0.84, 0.57, 0.39, 0.05 at 0.74, 0.80, 0.82, 0.84, 0.86.

### 3: m35mms against the radial motion (hst every 12 s, 3.6e5-7.8e5 s)
- The ring rms modulates at 3.36e-4 rad/s (resolution 1.5e-5), twice the ring frequency (1.64-1.68e-4). L_mid (comoving) and L_top have their lines at 3.16e-4 (also 3.04e-4, 3.52e-4 in the hst spectrum), the radial motion (Mr_int, PV_int, shell-mean v_r at 0.74) at 4.04e-4 rad/s. None sits at 3.36e-4; the closest, L_mid 3.16e-4, is 1.3 resolution elements away.
- Phase-locking of the modulation to L_mid, L_top, Mr_int, PV_int and the shell-mean v_r at 3.36e-4 (demodulated, window 4e4 s, 6 nearly independent samples over 4.2e5 s): resultant length 0.19, 0.24, 0.37, 0.36, 0.48 against 0.71 for a random phase at 95 %: **not locked to any of them.** (Coherence of the ring rms with KEr_int in the same band is below the 95 % level 0.67 except at 2.3e-4 and 4.1e-4, 0.52-0.63, also below.)
- Growth rate d ln R/dt (7-dump mean) against L_top, L_mid, L_int, KE_int, KEr_int and dt (detrended, 4.0e5-7.7e5 s): |r| <= 0.14 at lags of +-10 dumps; 95 % level for independent noise about 0.37: no correlation.
- Limits: the dump spacing (1800 s) gives frequencies to 1.75e-3 rad/s; only six independent samples for the phase test; a single star.

### 4 (26i check 2): older m35mms runs (all 256 x 256 in 20 x 20 deg, same radial grid)
| run (binary) | dumps | ring rms at 0.74 R_ph,IC [km/s] at t = 1.0e5 / 1.4e5 / 1.87e5 / 2.3e5 / 2.7e5 / 3.15e5 / 3.5e5 s | last dump: cells per wavelength, axis share | remark |
|---|---|---|---|---|
| prod_fresh_1121 (v203 e688cb7bc) | 6, to 6.5e4 s | - (4e-6 at 6.5e4) | - | too short |
| prod_fresh_1121b (j3 a5257b7ee, then p1, fwbs, fdef) | 36, to 3.87e5 | 7.2e-6 / 2.7e-5 / 2.4e-4 / 1.47e-3 / 5.9e-3 / 1.3e-2 / 2.3e-2 | 24 at 3.7e5, 0.72; then the last two dumps (3.78e5, 3.87e5) have a broad band (k256 16, axis share 0.08-0.11, 0.2 km/s) after its last restarts | a different state at the end |
| prod_ia11_1123 (ia11 bd9553537) | 27, to 2.79e5 | 6.1e-6 / 2.5e-5 / 2.1e-4 / 1.37e-3 / 5.0e-3 | 22 (11.6 cells), 0.65 | |
| prod_ia122_1123 (ia122 465aff485) | 21, 2.81e5 to 3.52e5 | - / - / - / - / - / 1.24e-2 / 2.1e-2 | 22, 0.79 | starts at 2.8e5 |
| prod_rc3_1124 (rc3 af836f5a3) | 29, to 4.06e5 | 6.5e-6 / 2.5e-5 / 2.0e-4 / 1.37e-3 / 5.2e-3 / 1.6e-2 / 2.4e-2 | 24 (10.7 cells), 0.79 | |
| prod_rc4v_1124 (rc4v 262affd6f) | 12, to 8.0e4 | - | - | too short |
| prod_c1_1124 (sdorder3, then c4speed) | 110+ | 6.6e-6 / 2.7e-5 / 2.0e-4 / 1.38e-3 / 4.8e-3 / 1.35e-2 / 2.2e-2 | 23 (11.1 cells), 0.78 | |
- **The ring has the same history in four binaries and schemes (bd9553537, rc3, j3/p1/... , c1/c4): the same to 2-8 % in ia11, rc3 and c1 at 1.4e5-2.3e5 s; including prod_fresh_1121b within 18 %, and within 30 % at 2.7e5-3.5e5 s.** It is deterministic (same IC noise) and scheme-independent in these runs.
- **No m35mms dump at another lateral resolution was found**: the searched dump directories (the runs above, b1_1124 test arms, cand_1124, preflight_1124) all use 256 x 256 cells in 20 x 20 deg; no half-wedge and no 128^2 m35mms dump exists there. Check 2 cannot separate cells per wavelength from l for m35mms; only the t9 and m35mms pair of 26i does (and that pair shows neither).

## Interpretation (not shown by the data)
- Axis preference is a property of the stable layers (all five), starts at the wall, and is absent in the convection; the growing ring adds to it in m35mms and the t9 stars. The time order says the base of the convection is earlier, but the gap and the wrong order of neighbours in m35mms do not look like a front from above; nothing measured pins the cause. The scheme-independent history argues against a bug of one scheme version; it does not separate "physical, seeded by the same noise" from "common to all schemes".
- No parametric link to the radial pulsation is seen at the frequencies the data resolve.

-- Caltech
