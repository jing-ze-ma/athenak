# T9L6.2 paper-BC run: why dt collapses (11-15, read-only diagnosis)

Run: /viper/ptmp2/jinma/t9l62_1113/prod_viper_paperbc (deck decks/t9l62_final.athinput, binary
athena_he_gpu72_x_0940c7d1_t9fix2 = commit 0940c7d14). Fresh start t = 0 (no restart).
Scripts and per-dump profiles: figs/dtdrop_1115/ (dtdrop.py, summ.py, prof_NNNNN.npz, run63.log, runB.log).

## 1. What sets dt (reproduced from the dumps)

dtdrop.py re-implements Hydro::NewTimeStep (hydro_newdt.cpp, `rad_signal_speed = true`):
c^2 = c_g^2 + (c_eq^2 - c_g^2)(1 - exp(-kappa_T rho dx)), c_eq^2 = Gamma_1(beta) (P_g + P_r)/rho,
P_r = E/3 (Rosseland table for kappa_T). Dump 0 reproduces the hst dt exactly (99.5 s); dump 63 gives
13.6 s vs the logged 13.8 s at that time.

| dump | t [h] | cfl*dt gas only | cfl*dt with rad | limiting cell |
|---|---|---|---|---|
| 0, 56-58 | 0-174 | 123 | 99.2 | base cell i=0 (24.3 Rsun), static |
| 59 | 177 | 107 | 88 | i=66 (38.2 Rsun), v_r = -4.3e7, rho 4.6e-8 |
| 60 | 180 | 94 | 77 | i=56 (35.7 Rsun), v_r = -4.5e7 |
| 61 | 183 | 87 | 61 | i=66 (38.2), rho 3.5e-9, c_eq 5.5e7 (c_g 7e6) |
| 62 | 186 | 87 | 29 | i=60 (36.7) |
| 63 | 189 | 19 | 13.6 | i=56 (35.7), rho 1.9e-10 (IC 1.1e-7), E 4.3e7, c_eq 3.2e8, kT dx 1.2 |

So 99 -> 14 s: the radiation-augmented sound speed in low-density, radiation-dominated voids at
r = 34-38 Rsun (i 48-66, rho_min down 3 decades below the shell mean in 6 h; dt ~ sqrt(rho_min)),
plus |v_r| up to 4.9e7 cm/s in downflows. Gas-only CFL alone would still give 19-87 s.

14 -> 0.2-4 s (cycle 7552 on, t >= 6.807e5 s = 189 h, run.log "dt COLLAPSE" x16): a DIFFERENT
limiter: near-floor cells (rho 1e-18 to 2e-17) at r = 7.6-7.9e12 cm (108-113 Rsun, i 223-227) with
P/rho 1e19-2.5e21 (T 1e11-2e13 K), |v| 1e8-3e8. Coincides with `implicit_wb_vkick_reject ... accepted
(no redo)`: the backward-Euler M1 solve kicks 1.5e3 (cycle 7572) -> 2.9e4 (cycle 7743) cells to
|v| > c (max 3e11 cm/s). First NOT ADMISSIBLE at cycle 7452 (~188.7 h). Last logged dt 0.27 s (cycle 7700).

## 2. Inner boundary and envelope

- Base (i 0-30, 24-30 Rsun): shell-mean rho and T within 1-3.5 % of the IC at all dumps; |v| < 5e5
  through dump 63 (max 3e4-1.6e5 at i 0-16 only in dump 63). No heating, no plume, no mass pile-up.
  The base layer is "still", as J18 say.
- L_in (hst col 7) = lab wall luminosity of the bath solve (he_star_m1 HeStarHist: lbath_f = comoving
  face flux + the row's donor enthalpy flux). Measured: L_in - L_IC/wedge = (4.1-6.5)e35 x sqrt(v1sq_wall)
  over 12-168 h, vs (4/3) E_wall A_wall = 7.6e35: L_in is dominated by (4/3) v_r E advected with the wall
  cells' slow radial drift (+5.8e3 cm/s at 12 h decaying to +2.2e3 at 168 h, then -2.5e3 uniform in
  angle at 183 h). L_in < 0 = the wall cells moving inward, not a hot base draining into the bath.
  L_bot (face is+1) stays +7.9e38 throughout. Hypothesis (unverified): counting advected radiation
  enthalpy through a zero-mass-flux face is a BC bookkeeping/energy artefact; it injected ~6.7e44 erg
  extra (Etot +0.3 % by 168 h, ~1.2 x L_IC t) - small, not the cause of the collapse.
- 32-100 Rsun: supersonic radial flows (v_r -4.9e7 .. +4.5e7; free fall from 80 to 38 Rsun = 5.4e7,
  escape 7.7e7 at 35 Rsun) into the layer where the IC turns radiation dominated (P_r/P_g 1 at 30 Rsun,
  5 at 36, 44 at 60). The IC density inversion (2.5e-9 at 60 -> 9e-9 at 80 Rsun) is being destroyed:
  shell rho at 63.5 Rsun 1.5x IC, at 83.5 0.43x; at 96-103 Rsun 3e3-2e9 x IC (envelope expanding,
  v_r +1.9e7 at 100 Rsun). KE_int 1.7e39 (84 h) -> 4.6e43 (168 h), e-fold ~7 h, 90 % radial.
- Above the IC photosphere (>= 107 Rsun) the IC is a floor atmosphere (rho = dfloor 3.6e-19, T 2.1e4 K,
  no pressure support) which free-falls: v_r = -2.3e7 cm/s at 107-118 Rsun by 174 h. The expanding
  envelope (+1.9e7) meets this infalling floor gas at 108-116 Rsun; that interface is exactly where
  the hot cells (208 with T > 1e7 K in dump 63, 0 in dumps 56-61) and the v > c kicks appear.

## 3. Physical or artefact

- dt 99 -> 14 s (175-189 h): physical-looking. J18 (1809.10187 l.94-100): convection destroys the
  density inversion after ~10 t_dyn (~43 h), clumps rise/expand/cool, the envelope expands with
  ~0.05 Msun/yr until ~400 h. We see the same sequence, ~3x later (seed n2 1e-3). The cost comes from
  our rad_signal_speed CFL in optically-thin-ish voids (kT dx ~1); J18 (Athena++) use the gas CFL,
  so their dt would not see c_eq. Not an inner-BC artefact (base unchanged).
- dt 14 -> 0.3 s (>= 189 h): numerical artefact (hypothesis with strong evidence): infalling floor
  atmosphere x expanding envelope at floor density, where the BE radiation momentum/energy exchange
  per unit mass is huge; the vkick guard accepts > c velocities, vceil 3e8 clips them, T runs to
  1e11-1e13 K. The kicked-cell count grows ~20x in 300 cycles: the run will very likely crash or stall.

## 4. Projection and options (none applied)

- Wall rate (run.log): 1.9-2.3 s/cycle now (3.1 s at dt 80). At dt 14 s: ~4 h wall per model day;
  at dt 0.3-4 s: 12 h to >7 d wall per model day; 611 model hours left. Not viable as is.
- Options: (a) restart from rst 15 (180 h, before hot cells; dt already 77) or rst 14 (168 h) with the
  floor atmosphere fixed: e.g. remove/reset the IC gas above the photosphere (J18 T19L6.4: "after
  removing the initial gas above that region"), a sponge/infall clamp in floor cells, or exclude
  floor-density cells from the M1 momentum kick (vkick_reject redo instead of accept); smoke first.
  (b) cap rad_signal_speed in voids (e.g. rad_signal_speed = false as Athena++) to keep dt ~gas CFL
  (~19-90 s) - changes stability margin, needs a test. (c) switch to the our-wall twin (prod_viper,
  132 h): same IC atmosphere, so it should hit the same floor-interface failure; check before use.
  Keeping the current job going is not useful.
