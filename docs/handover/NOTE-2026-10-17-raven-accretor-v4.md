# NOTE 2026-10-17 raven: accretor v4 (HANDOVER 01fdd1c3) - FINAL (supersedes the INTERIM NOTE-2026-10-16-raven-accretor-v4.md)

Code accretor-bigdom-1015 @ 99adf52d (v4 code 5b0d4bc1, binary md5 ac371900). Raven dir /raven/ptmp/jinma/accbig_1016
(v4/); full numbers + RESUME at the end of its NOTE.md. Items 1-3 (gate BITWISE, start-up PASS, s = 0.22 calibration,
static clean) are in the INTERIM NOTE and unchanged.

## VERDICT
- **v4 stream through impact: robust but FAILS the dt / atmosphere goals.** 0 FATAL, 0 NON-CONVERGED, 49 retries all
  recovered, Mnzf exact, back-sink row at rest, impact = v3. But dt is 2-7e-7 through the impact, 2-3x WORSE than v3
  (0.84-1.6e-6), and the thin atmosphere moves ~2x faster than in v3.
- **Cause: the back-sink sector's stream light (new in v4, unscaled) floods the domain** (radiation E at r 9.4-10 up
  50-100x over the static phase, 3-5x over v3); the vet_col lateral beam carries it round the limb. The s = 0.22 donor
  calibration is irrelevant in the stream phase (donor part < 2 % of E). **The s ~0.30 round is not warranted** until
  the sector light is decided.
- **Fix 3: vet_col does mis-handle lateral donor light** - it becomes an azimuthal beam that lights the night side as
  strongly as the facing side. The thinclos Eddington/M1 reference arms are unusable (dt collapse in all four, also
  without donor light).
- Nothing submitted on 10-17. Nothing merged.

## 1. v4 stream (v4/cal/stream; gpudev 31229257/58/59 + gpu 31229101; t 0.00946 -> 0.03597, cycle 33354)
| check | result | |
|---|---|---|
| FATAL / NC | 0 / 0 (all 4 link summaries); Picard mean 2.9-4.0, max 22 (from t ~0.028) | PASS |
| retries | 49: 48 eint<=0 (cut 1), 1 resid>fatal (3.3e14, t 0.03434); all recover in ~13 cycles; 13 before t 0.026, first t 0.0189 (v3 whole run: 8) | MARGINAL |
| Mnzf rate | 19.7952-19.7991 in every link (target 19.795 +-0.5 %: dev <= 0.02 %); Min-Mnzf 2.4e-5 | PASS |
| dt | 7.03e-6 static -> 5.2e-6 (t 0.0115) -> 2.9-3.4e-6 (0.015-0.018) -> 8.5e-7 (0.022) -> 2-7e-7 (0.023-0.036) | FAIL |
| thin atmosphere max\|v\| r 9.1-9.6 | 660 / 1181 / 1430-1775 / 1540-3170 km/s at t 0.013 / 0.015 / 0.020-0.022 / >0.024 (v3, same script: 376 / 479 / 717-852 / 980-1160) | FAIL |
| facing flux vs 1.52e5 | metric invalid in the stream phase: net F_r OUTWARD at r 10-10.4 facing (-(dF_r) = -1e6..-1.4e7), E/E_noirr 20-170 at all phi | FAIL |
| r 14.607 back-sink row | \|v\| = 0, rho = 1e-8 in every output (v3: up to 9.7e5 km/s) | PASS |
| Mmsk / Emsk | smooth within links; Mmsk bumps when stream gas enters the lobe (t 0.0225, 0.0285) | OK |
| Eirr | rate +2.5e5 (static) -> -9e6 .. -4.5e8 (stream): the lobe is a net sink of the new light | FLAG |
| impact | dense infall at r 9.4-9.65 from t ~0.026 (v3 0.027); dense-gas heating from t 0.0308, e/rho 3054 -> 4035, phi 13.3 -> 15.0 deg, v_r -394..-346 (v3: 3294 -> 4044, phi 14.0 -> 15.2); ridge 13.6-14.3 deg | PASS (= v3) |
| cost | 0.46 s/cycle, 4 A100 (v3 0.37) | |

### dt-drop cause (tool/dtcell.py on every hydro_w bin; agrees with the rank-0 dt_report)
1. t 0.0106-0.018: thin atmosphere at the **limb**, r 9.2-9.6, phi +82..88 and -84..-86 deg (symmetric), rho ~1e-6,
   v_r +120..+230 and v_phi +-257 -> +-1280 km/s (day -> night), c 53 km/s.
2. t 0.019-0.021: that gas falling back on the limb/night side (v_r -244..-405 at r 9.07-9.39).
3. t >= 0.0225: **floor gas** (rho 1e-8..2e-7) at r 9.05-9.4, all phi, heated to e/rho 1e6-1.8e6 (c 1300-2000
   km/s, \|v\| up to 6e3 km/s incl. v_theta in the 4-cell band). In v3 the floor gas stays at e/rho 2.4e3 and the dt
   cell is dense stream gas (c 51 km/s).

Driver: the radiation field after the switch-on. E at r 9.4-10: 3 (static, t 0.00946) -> 55-165 (t 0.0142) -> 160-350
(t 0.0284; v3 at the same t 50-75), F_phi/cE = +-1 at the limb; Ltop 5.1e4 -> 1.1e6 (t 0.0105) -> 4.2e6 (v3 3.3e5 ->
9.5e5). The only light source new in v4 is the back-sink sector holding the stream's own light, NOT scaled: core
72966 K, sector-mean T^4 = (56.2 kK)^4, i.e. **132x the scaled donor surface value 0.22 (24.2 kK)^4** (core 376x).
This is an inference (v3 vs v4 + source budget); an on/off run needs a new key (the sector light is tied to
stream_rad).

## 2. Fix 3 (gpu 31229100, static to t 0.00237, thinclos test binary 4f0dbbfc)
- **thinclos arms all collapse** (FATAL dt < dt_min): edd at cycle 612 (t 3.1e-4), edd_noirr 172, m1_noirr 19, m1_s022
  33; eint<=0 retries from cycle 5-10 also without donor light, so the thin-closure switch itself is unstable here.
- **vet_col donor part** (v4/s1 s = 1 minus v4/fix3/vc_noirr, both t 0.00237): dE = 3.2-4.5 at r 9.4-10.0 at every
  phi; **night side (phi 180) / facing dE = 0.98-1.03** at r 9.6-10 (geometric shadow: should be ~0). The donor-part
  flux is lateral: \|dF_phi\|/(c dE) 1.0-1.85 at phi 10-120 deg, dF_r/(c dE) only -0.1..-0.4 (facing r 9.6: -0.17;
  r 10: -1.4). Total F_phi/cE +0.93..+1.00 at phi 10-120 with irradiation (noirr -0.6..+0.5).
- Answer to the user's question: **yes, vet_col mis-handles lateral flux**: azimuthal donor light becomes a free-
  streaming phi beam wrapping the star, so roughly the whole circumference is lit instead of the day side (~2x in the
  circumference-integrated heating), and the incident radial component at r 9.6 is ~6x below cE. Not quantified
  against a reference closure (none runs).

## Open decisions (USER)
1. Back-sector stream light: scale it (donor_irr_scale or s_geo), surface T only, or dark as in v3; plus a key to
   switch it off for a confirming test (code change, Caltech/user). Then a short stream test from
   v4/cal/static/rst/ryper.00001.rst (3 gpudev links to t ~0.02 show the limb/floor effect).
2. Fix 3 follow-up: a lateral-beam remedy in vet_col (closure work) vs keeping the calibration; thinclos needs
   debugging before it can serve as a reference.
3. The s ~0.30 calibration round only after 1.
4. MI300A illegal-memory bug (Viper) still open.
