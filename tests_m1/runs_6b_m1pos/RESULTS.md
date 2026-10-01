# m1-positivity (09-30): positivity of the implicit M1 solve in the He-star wedge plume

Branch `m1-positivity` (off `he-presn-m1` abf547b5), worktree /viper/ptmp2/jinma/wt_m1pos,
run tree /viper/ptmp2/jinma/m1pos_0930 (diag*/ iso*/ smoke/ gate_off/ runs/ hebox/).
Binary (He star, HIP ROCm 7.2): /viper/ptmp2/jinma/m1pos_0930/bin/athena_gpu_he_star_m1_7192a176.

## 1. What fires first (instrumented, read-only: dbg_t2_admiss now also reports the
implicit_vimp fallback and every solve that ends with E <= e_floor or eint <= 0)

ad3d_128 rst 00011 (t = 25850 s, 5.5 turnovers), keys of ad3d_128/r2.sh, 13 cycles (diag0/diag1):
- every hesdirk2 stage (13/13) is redone with backward Euler. The stage is CONVERGED and its
  final E, T are positive (min(E,T) = 0.2177); it is rejected only because the implicit_vimp
  positivity fallback fired on the FIRST Picard pass (vimp_now != impl_vimp).
- that pass-0 linear solve has E < 0 in 281 + 91 cells (min -2.28e8 against E ~ 7e6), in dense
  plume cells (rho/<rho> 1-10, r = 1.46-1.65e11, v_r 1-3e7). Their rows are Z-matrices (0 positive
  off-diagonals, vimp block ~1e-2 of the diagonal); the right-hand side is the problem:
  x1 face term -c dt v_f g0_f = -1.9e8 (base 2.3e7), opacity-Newton rhs +-1e7..9e7.
- g0 = rho (kappa_E E - kappa_P a T^4) of the ITERATE: in these radiation-dominated, stiffly
  coupled cells (c dt rho kappa_P ~ 5e5, G0 ~ 30, chat dt g0 ~ 4e12 >> E) it is the difference
  of two numbers ~1e5 x the converged exchange; at pass 0 it is garbage.
- then (g0 fixed): the opacity-Newton NEIGHBOUR entries keep the right sign but break row
  dominance (sum|off|/diag = 2.4 at r = 1.46e11, 2.0 at the photosphere r = 2.39e11): E < 0.
- then (both fixed, vimp kept): Picard limit cycle at res ~1.5e-7 > tol 1e-8 (iso5/E plog),
  6/8 solves NON-CONVERGED -> stage rejected anyway.
- gas eint <= 0 appears in the OLD vector (after hydro) in 37 near-void cells (rho/<rho> down
  to 2e-5); the radiation write-back makes them positive in all but 1 of 26 solves.
- The transverse operator is NOT the trigger: advection is already upwinded in x2/x3
  (M1SphTransRow / m1_impl_tcell) and every failing row had 0 positive off-diagonals.
  Item 1 of the brief (ap_hll on x2/x3 faces) was therefore not built.

## 2. New keys (all read only when named: absent = bitwise the old code, restart echo too)
- `implicit_g0_exchange = true`: from pass 1 on g0 = -(SRCR - SRCB E^k)/(chat dt), the
  exchange of the previous pass's linearised source row (identical to the pointwise g0 at the
  Picard fixed point).
- `implicit_g0_limit = w` (1.0): pass 0 only, |chat dt g0| <= w (max(E^k, E_old) + (chat/c)
  max(e_gas,0)). Clipping every pass was tried and REJECTED (1.7e7 cell-passes clipped, 5/6 NC).
- `implicit_opac_newton_guard_mode = 6` (bit 4 new): a row's remaining Newton face terms are
  dropped when they raise sum|off|/diag above max(1, value without them) x (1 + guard).
- `implicit_pos_gas = true` (+ `implicit_pos_gas_frac`, 1e-3): written-back gas eint below
  max(e(rho, tfloor), frac x max(e_gas old, 0)) is raised with energy taken from E (down to
  e_floor); e_gas + (c/chat) E per cell exact; counted (cells, energy).
- `implicit_pos_floor = true`: E raised to e_floor is paid by the gas (down to its limit);
  the uncovered part is counted as "energy created".
Counters: one `<rad_m1> m1-positivity:` line in the end-of-run report.

## 3. Short tests from rst 00011 (4 cycles, apudev, binary pos4 = 7192a176)
| arm | keys | Picard mean / NC | stage fallbacks | vimp fallbacks | E<=floor / eint<=0 |
|---|---|---|---|---|---|
| base (diag0, 13 cyc) | none | 19.5 / 2 of 26 | 13/13 | 26 (min E -9.7e8) | 31 / 1 |
| iso4 A | exchange | 16.7 / 0 | 3/3 | 6 (-2.3e8, pass 0) | 0 / 0 |
| iso5 G | all + guard 6, vimp on | 156 / 6 of 8 | 4/4 (NC) | 0 (min E +2086) | 0 / 0 |
| iso6 H | vimp off only | 49 / 2 of 9 | 3 | - | 18 / 4 |
| iso6 J | all + guard 6, vimp off | 28 / 0 of 8 | 1 of 4 | - | 2 / 3 (both limited, 0 created) |

## 4. Gates
### 4a. Off switches bitwise (apudev, HIP ROCm 7.2), keys absent vs he-presn-m1 abf547b5
- He-star wedge 128^2, ad3d_128 rst 00011 -> t 25880, 2 GPUs (job 12046531, gate_off/wb vs wn):
  both hst and rst 00012-00014 IDENTICAL (base binary athena_gpu_he_star_m1_base, new 7192a176).
- He box, hebox_cfl2 H6 rst 00036 (t 59800) -> 59830, 1 GPU (job 12046647, gate_off/bb vs bn):
  both hst and rst 00037-00039 IDENTICAL (athena_box_gpu72_abf547b5 vs athena_box_gpu72_7192a176).
  (A first attempt from rst 00037 = t 61800 ran 0 cycles and is void.)

### 4b. He box cost, GPU (apudev, 2 GPUs, same binary, interleaved C N C N C N, 59800 -> 59880,
job 12046584, hebox/cost): wall (cpu time used, s) C 25.01 / 25.37 / 26.01, N 24.97 / 25.12 / 25.32;
Picard 7.0 per solve and 496 solves in every run.  N clips g0 on pass 0 in 42 cell-passes; the gas
and floor limiters never fire.  Cost: no measurable change (N -1 %, within the C spread).

### 4c. He box accuracy (saturated window 59800 -> 74000, cfl 0.6, C / P (1-ulp cfl noise) / N):
job 12046701 (apu, submitted 23:05, smoke 12046649 clean). PENDING. Analysis:
`/viper/ptmp2/jinma/m1pos_0930/hebox/ana.sh 60000 74000` (copy: hebox_ana.sh) -> ana_refC_*.txt,
conv_*.txt (v1'/v_MLT, r(v1',T'), Fc/F in the FeCZ), kehst_*.txt (KE1, KEh), counters.txt.
Accept if |N - C| <= |P - C| (and the snapshot scatter) for v1'/v_MLT, KE1, KEh, Fc/F, mean T/rho.

### 4d. DECISIVE He-star restarts (apu, 1 node = 2 GPUs each; smokes 12046496, 12046497 clean)
Arms (arms.sh; ON = implicit_g0_exchange + implicit_g0_limit 1 + implicit_pos_gas +
implicit_pos_floor + implicit_opac_newton_guard_mode 6):
| arm | job | state at 23:05 |
|---|---|---|
| off128 (control, r2.sh keys) | 12046524 | FAILED like ad3d_128: last sane cycle 2580 t = 26846 s, then dt 2262 (runaway); E<=floor 1.7e8, eint<=0 1.5e8 cell-solves, vimp fallbacks 440 (min E -inf), stages 12 ok / 234 redone |
| be128 (time_scheme be, vimp off, options off) | 12046527 | FAILED the same way: last sane t = 26847 s, E<=floor 7.6e7 |
| onv128 (ON, vimp on) | 12046525 | running, t = 27970 s (past 26825), dt 3.6, 300 stages non-admissible (Picard NC) |
| onn128 (ON, vimp off) | 12046526 | running, t = 28084 s, dt 3.7, 324 stages non-admissible |
| beon128 (ON, be, vimp off) | 12046528 | running, t = 29752 s, dt 2.4 |
| offh64 (hllc 64^2 control, rst 00012 t = 28200) | 12046529 | pending |
| onnh64 (hllc 64^2, ON, vimp off) | 12046530 | running, t = 32047 s (past the 29200 onset and the 31626 NaN of ad3d_hllc), dt 4.5 |
Targets: 128^2 t >= 32900 s, 64^2 t >= 33000 s. Analysis: `python3 /viper/ptmp2/jinma/m1pos_0930/analyze.py`
(time reached, nan lines, non-admissible stages, end-of-run counters incl. the m1-positivity line).
